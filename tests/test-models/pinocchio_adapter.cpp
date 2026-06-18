/**
 * @file pinocchio_adapter.cpp
 * @brief Pinocchio adapter implementation — bridges test-models to Pinocchio
 * @details This is the ONLY file in tests/test-models/ that includes Pinocchio
 *          headers (following the D-12 compilation firewall pattern established
 *          by sa_adapter.cpp). It implements the PIMPL pattern declared in
 *          pinocchio_adapter.h.
 *
 *          Conversion patterns (PCC-04):
 *          - Matrix4d → pinocchio::SE3 via SE3(rotation, translation)
 *          - mass + com + dense-inertia → pinocchio::InertiaTpl
 *          - JointSpec axis → JointModelRevoluteUnaligned / JointModelPrismaticUnaligned
 *
 *          Model building follows Pinocchio standard practice:
 *          1. model.addJoint(parent_id, joint_model, placement, name)
 *          2. model.appendBodyToJoint(joint_id, inertia, name)
 *
 *          Joint index mapping:
 *          Pinocchio joint 0 = universe (fixed). Our joints 0..dof-1 map to
 *          Pinocchio joints 1..dof. RobotModel parent -1 → Pinocchio parent 0.
 *
 *          Gravity convention:
 *          model.gravity is set to Motion(angular=0, linear=g) before each
 *          call to rnea() or aba(). This matches how SpatialAlgebra's
 *          InverseDynamics/ForwardDynamics accept a gravity vector parameter.
 *
 * @see pinocchio_adapter.h for public API
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include "pinocchio_adapter.h"

// Pinocchio headers — ONLY file that includes these (per compilation firewall)
#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/multibody/joint/joint-revolute-unaligned.hpp>
#include <pinocchio/multibody/joint/joint-prismatic-unaligned.hpp>
#include <pinocchio/multibody/joint/joint-composite.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace test_models {

// =========================================================================
// PIMPL definition — hides all Pinocchio types from the public header
// =========================================================================

struct PinocchioAdapter::Impl {
    pinocchio::Model model;                 ///< Pinocchio kinematic model
    pinocchio::Data data;                   ///< Pinocchio computation cache
    std::vector<pinocchio::JointIndex> jointIds;  ///< Map: our joint i → Pinocchio joint ID
    Eigen::VectorXd q_vec;                  ///< Cached joint positions (size model.nq)
    Eigen::VectorXd v_vec;                  ///< Cached joint velocities (size model.nv)

    /// Default gravity (matching SA convention: -Z by default)
    static constexpr double kDefaultGravity = -9.81;
};

// =========================================================================
// File-scope helpers
// =========================================================================

namespace {

using Vector3d = Eigen::Vector3d;

/**
 * @brief Map a unit axis vector to the nearest Pinocchio cardinal joint type
 * @details For efficiency, Pinocchio provides axis-aligned joint types
 *          (JointModelRX/RY/RZ, JointModelPX/PY/PZ). We fall back to
 *          the unaligned version when the axis is not along a cardinal
 *          direction. For our test models, all axes are cardinal (X, Y, or Z).
 *
 *          Returns the appropriate JointModelVariant for the given axis.
 */
static pinocchio::JointModel buildJointModel(JointType type,
                                              const Vector3d& axis) {
    // Determine axis index: 0=X, 1=Y, 2=Z, -1=custom
    int axisIdx = -1;
    if (axis.isApprox(Vector3d::UnitX(), 1e-12)) axisIdx = 0;
    else if (axis.isApprox(Vector3d::UnitY(), 1e-12)) axisIdx = 1;
    else if (axis.isApprox(Vector3d::UnitZ(), 1e-12)) axisIdx = 2;

    switch (type) {
        case JointType::REVOLUTE:
            if (axisIdx == 0)
                return pinocchio::JointModelRX();
            else if (axisIdx == 1)
                return pinocchio::JointModelRY();
            else if (axisIdx == 2)
                return pinocchio::JointModelRZ();
            else
                return pinocchio::JointModelRevoluteUnaligned(axis);
        case JointType::PRISMATIC:
            if (axisIdx == 0)
                return pinocchio::JointModelPX();
            else if (axisIdx == 1)
                return pinocchio::JointModelPY();
            else if (axisIdx == 2)
                return pinocchio::JointModelPZ();
            else
                return pinocchio::JointModelPrismaticUnaligned(axis);
        case JointType::FIXED:
            // Pinocchio 3.x does not have JointModelFixed in its default
            // joint collection. FIXED joints (0 DOF) are not used in any
            // current test model — throw if encountered.
            throw std::invalid_argument(
                "FIXED joints are not supported by the Pinocchio adapter. "
                "Use REVOLUTE or PRISMATIC joints.");
    }
    // Should never reach here
    throw std::logic_error("Unhandled joint type in buildJointModel");
}

/**
 * @brief Convert Eigen::Matrix4d (homogeneous) to pinocchio::SE3
 * @param T 4×4 homogeneous transform [R | t; 0 | 1]
 * @return Equivalent SE3(rotation, translation)
 */
static pinocchio::SE3 matrix4dToSE3(const Eigen::Matrix4d& T) {
    return pinocchio::SE3(
        T.topLeftCorner<3, 3>(),
        T.topRightCorner<3, 1>()
    );
}

/**
 * @brief Build a pinocchio::Inertia from mass, COM, and dense inertia
 * @param mass Link mass
 * @param com Center of mass in link frame
 * @param rotInertia 3×3 rotational inertia at COM (will be symmetrized)
 * @return Pinocchio Inertia object
 * @note We symmetrize the inertia matrix before passing to Pinocchio,
 *       matching the same pattern used in sa_adapter.cpp (Pitfall 2 guard).
 */
static pinocchio::Inertia buildInertia(double mass,
                                        const Vector3d& com,
                                        const Eigen::Matrix3d& rotInertia) {
    Eigen::Matrix3d symInertia = 0.5 * (rotInertia + rotInertia.transpose());
    return pinocchio::Inertia(mass, com, symInertia);
}

} // anonymous namespace

// =========================================================================
// Destructor — defined here where Impl is complete (PIMPL requirement)
// =========================================================================

PinocchioAdapter::~PinocchioAdapter() = default;

// =========================================================================
// build() factory — per D-05 builder pattern
// =========================================================================

std::unique_ptr<PinocchioAdapter> PinocchioAdapter::build(
    const RobotModel& model) {
    int dof = model.getDOF();

    auto adapter = std::unique_ptr<PinocchioAdapter>(
        new PinocchioAdapter());
    adapter->dof_ = dof;
    adapter->impl_ = std::make_unique<Impl>();

    auto& pinModel = adapter->impl_->model;
    auto& jointIds = adapter->impl_->jointIds;

    // Add joints in topological order
    jointIds.reserve(static_cast<size_t>(dof));

    for (int i = 0; i < dof; ++i) {
        const auto& js = model.joints[static_cast<size_t>(i)];

        // Map parent index: RobotModel parent -1 → Pinocchio universe (0)
        // Pinocchio joints are 1-indexed (0 = universe), so our joint i maps
        // to Pinocchio joint i+1, and parent i-1 maps to Pinocchio parent i.
        pinocchio::JointIndex parentId = (js.parent < 0)
            ? 0
            : static_cast<pinocchio::JointIndex>(js.parent + 1);

        // Build joint model from type and axis
        pinocchio::JointModel jointModel = buildJointModel(js.type, js.jointAxis);

        // Joint placement: parentToJoint transform
        pinocchio::SE3 placement = matrix4dToSE3(js.parentToJoint);

        // Create joint name
        std::string jointName = js.name.empty()
            ? "joint_" + std::to_string(i)
            : js.name;

        // Add joint to model
        pinocchio::JointIndex jid = pinModel.addJoint(
            parentId, jointModel, placement, jointName);
        jointIds.push_back(jid);

        // Add body inertia (identity placement — inertia is in link frame)
        pinocchio::Inertia bodyInertia = buildInertia(
            js.mass, js.com, js.inertia);
        pinModel.appendBodyToJoint(jid, bodyInertia, pinocchio::SE3::Identity());
    }

    // Create computation cache
    adapter->impl_->data = pinocchio::Data(pinModel);

    // Initialize state vectors to zeros
    adapter->impl_->q_vec = Eigen::VectorXd::Zero(pinModel.nq);
    adapter->impl_->v_vec = Eigen::VectorXd::Zero(pinModel.nv);

    // Set default gravity
    adapter->impl_->model.gravity.linear(
        Eigen::Vector3d(0.0, 0.0, Impl::kDefaultGravity));
    adapter->impl_->model.gravity.angular(Eigen::Vector3d::Zero());

    adapter->fkDirty_ = true;

    return adapter;
}

// =========================================================================
// Private helper — autoFK() dirty-flag pattern
// =========================================================================

void PinocchioAdapter::autoFK() {
    if (fkDirty_) {
        forwardKinematics();
        fkDirty_ = false;
    }
}

// =========================================================================
// setState() — with input validation
// =========================================================================

void PinocchioAdapter::setState(const Eigen::VectorXd& q,
                                 const Eigen::VectorXd& qdot) {
    if (q.size() != dof_ || qdot.size() != dof_) {
        throw std::invalid_argument(
            "PinocchioAdapter::setState: vector size mismatch. Expected " +
            std::to_string(dof_) + ", got q=" +
            std::to_string(q.size()) + " qdot=" +
            std::to_string(qdot.size()));
    }

    // Copy to internal state vectors
    impl_->q_vec = q;
    impl_->v_vec = qdot;

    fkDirty_ = true;
}

// =========================================================================
// computeTorques() — RNEA inverse dynamics
// =========================================================================

Eigen::VectorXd PinocchioAdapter::computeTorques(
    const Eigen::VectorXd& qddot,
    const Eigen::Vector3d& gravity) {

    if (qddot.size() != dof_) {
        throw std::invalid_argument(
            "PinocchioAdapter::computeTorques: qddot size mismatch. Expected " +
            std::to_string(dof_) + ", got " +
            std::to_string(qddot.size()));
    }

    autoFK();

    // Set gravity for this computation
    impl_->model.gravity.linear(gravity);
    impl_->model.gravity.angular(Eigen::Vector3d::Zero());

    // Run RNEA
    pinocchio::rnea(impl_->model, impl_->data,
                     impl_->q_vec, impl_->v_vec, qddot);

    // Return joint torques
    return impl_->data.tau;
}

// =========================================================================
// computeAccelerations() — ABA forward dynamics
// =========================================================================

Eigen::VectorXd PinocchioAdapter::computeAccelerations(
    const Eigen::VectorXd& tau,
    const Eigen::Vector3d& gravity) {

    if (tau.size() != dof_) {
        throw std::invalid_argument(
            "PinocchioAdapter::computeAccelerations: tau size mismatch. Expected " +
            std::to_string(dof_) + ", got " +
            std::to_string(tau.size()));
    }

    autoFK();

    // Set gravity for this computation
    impl_->model.gravity.linear(gravity);
    impl_->model.gravity.angular(Eigen::Vector3d::Zero());

    // Run ABA
    pinocchio::aba(impl_->model, impl_->data,
                    impl_->q_vec, impl_->v_vec, tau);

    // Return joint accelerations
    return impl_->data.ddq;
}

// =========================================================================
// forwardKinematics() — recompute world transforms from joint angles
// =========================================================================

void PinocchioAdapter::forwardKinematics() {
    // Compute joint placements and spatial velocities
    pinocchio::forwardKinematics(impl_->model, impl_->data,
                                  impl_->q_vec, impl_->v_vec);
}

// =========================================================================
// getJointTransform(idx) — world-frame joint placement
// =========================================================================

Eigen::Matrix4d PinocchioAdapter::getJointTransform(int idx) const {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "PinocchioAdapter::getJointTransform: index " +
            std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    const_cast<PinocchioAdapter*>(this)->autoFK();

    // oMi[joint_id] gives the world-frame placement of joint id
    const auto& oMi = impl_->data.oMi[impl_->jointIds[static_cast<size_t>(idx)]];
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = oMi.rotation();
    T.topRightCorner<3, 1>() = oMi.translation();
    return T;
}

// =========================================================================
// computeMassMatrix() — CRBA method (O(n²), exact)
// =========================================================================

Eigen::MatrixXd PinocchioAdapter::computeMassMatrix() {
    autoFK();

    // Composite Rigid Body Algorithm
    pinocchio::crba(impl_->model, impl_->data, impl_->q_vec);

    // Symmetrize: CRBA computes only the upper triangular part
    impl_->data.M.triangularView<Eigen::StrictlyLower>() =
        impl_->data.M.transpose().triangularView<Eigen::StrictlyLower>();

    return impl_->data.M;
}

// =========================================================================
// computeGravityTorques() — gravity-only RNEA
// =========================================================================

Eigen::VectorXd PinocchioAdapter::computeGravityTorques(
    const Eigen::Vector3d& gravity) {

    autoFK();

    // Set gravity for this computation
    impl_->model.gravity.linear(gravity);
    impl_->model.gravity.angular(Eigen::Vector3d::Zero());

    // Zero-velocity, zero-acceleration RNEA = gravity torques
    Eigen::VectorXd zeroVec = Eigen::VectorXd::Zero(dof_);
    pinocchio::rnea(impl_->model, impl_->data,
                     impl_->q_vec, zeroVec, zeroVec);

    return impl_->data.tau;
}

// =========================================================================
// computeJointSpaceJacobian(idx) — geometric Jacobian in world frame
// =========================================================================

Eigen::MatrixXd PinocchioAdapter::computeJointSpaceJacobian(int idx) {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "PinocchioAdapter::computeJointSpaceJacobian: index " +
            std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    autoFK();

    // Compute all joint Jacobians
    pinocchio::computeJointJacobians(impl_->model, impl_->data, impl_->q_vec);

    // Get Jacobian for the specified joint in WORLD frame
    // Returns a 6×nv matrix where nv == dof for our 1-DOF joints
    Eigen::MatrixXd J(6, impl_->model.nv);
    J.setZero();

    pinocchio::getJointJacobian(
        impl_->model, impl_->data,
        impl_->jointIds[static_cast<size_t>(idx)],
        pinocchio::WORLD, J);

    // Resize to dof if needed (nv should equal dof for 1-DOF joints)
    if (impl_->model.nv < dof_) {
        J.conservativeResize(6, dof_);
        J.rightCols(dof_ - impl_->model.nv).setZero();
    }

    return J;
}

// =========================================================================
// getLinkCOM(idx) — world-frame center of mass
// =========================================================================

Eigen::Vector3d PinocchioAdapter::getLinkCOM(int idx) const {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "PinocchioAdapter::getLinkCOM: index " +
            std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    const_cast<PinocchioAdapter*>(this)->autoFK();

    // Get the body inertia for this joint and transform its lever (COM)
    // to world frame: worldCOM = oMi * lever (where lever = COM in body frame)
    pinocchio::JointIndex jid = impl_->jointIds[static_cast<size_t>(idx)];
    const auto& oMi = impl_->data.oMi[jid];
    const auto& inertia = impl_->model.inertias[jid];

    // worldCOM = R * localCOM + t
    return oMi.rotation() * inertia.lever() + oMi.translation();
}

// =========================================================================
// getDOF()
// =========================================================================

int PinocchioAdapter::getDOF() const {
    return dof_;
}

} // namespace test_models
