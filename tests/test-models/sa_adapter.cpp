/**
 * @file sa_adapter.cpp
 * @brief SpatialAlgebra adapter implementation — bridges test-models to SA solver
 * @details This is the ONLY file in tests/test-models/ that includes SpatialAlgebra
 *          headers (per D-12 compilation firewall). It implements the PIMPL pattern
 *          declared in sa_adapter.h: struct Impl contains ForwardDynamics and
 *          InverseDynamics solver objects, joint type/axis metadata for forward
 *          kinematics, and a world-frame transform cache.
 *
 *          Conversion patterns follow benchmarks/common/model_factory.cpp:
 *          - Matrix4d → PluckerTransform via Rotation(R) + Vector3d(t)
 *          - mass+com+dense-inertia → RigidBodyInertia with symmetrization (Pitfall 2)
 *          - Vector3d+JointType → MotionVector screw axis (Pitfall 3 convention)
 *
 *          Error handling: input validation at every entry point per threat model
 *          T-20-01 (size checks) and T-20-02 (bounds checks). Exceptions propagated
 *          from SA solver are passed through unchanged.
 *
 * @see sa_adapter.h for public API
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include "sa_adapter.h"

// SpatialAlgebra headers — ONLY file that includes these (per D-12)
#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include "LowerTriangular.h"
#include "Rotation.h"
#include "MotionVector.h"
#include "SpatialVector.h"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>

namespace test_models {

// =========================================================================
// PIMPL definition — hides all SA types from the public header
// =========================================================================

struct SpatialAlgebraAdapter::Impl {
    SpatialAlgebra::ForwardDynamics fd_;       ///< ABA solver
    SpatialAlgebra::InverseDynamics id_;       ///< RNEA solver
    std::vector<JointType> jointTypes_;        ///< Joint type per link
    std::vector<Eigen::Vector3d> jointAxes_;   ///< Joint axis per link
    std::vector<SpatialAlgebra::PluckerTransform> worldX_; ///< World-frame transforms
};

// =========================================================================
// Anonymous namespace — file-scope conversion helpers
// =========================================================================

namespace {

using Vector3d = Eigen::Vector3d;

/**
 * @brief Build a ForwardDynamics solver from a RobotModel
 * @param model Kinematic chain description
 * @param dof Number of degrees of freedom
 * @param[out] jointTypes Filled with JointType per link
 * @param[out] jointAxes Filled with joint axis per link
 * @return Fully constructed ForwardDynamics solver object
 * @details Conversion follows model_factory.cpp pattern:
 *          Matrix4d → PluckerTransform, mass+com+dense-inertia → RigidBodyInertia,
 *          Vector3d+JointType → MotionVector screw axis.
 *          Inertia is symmetrized before conversion per Pitfall 2.
 */
SpatialAlgebra::ForwardDynamics buildFD(const RobotModel& model, int dof,
                                         std::vector<JointType>& jointTypes,
                                         std::vector<Eigen::Vector3d>& jointAxes) {
    SpatialAlgebra::ForwardDynamics fd;
    fd.links.reserve(static_cast<size_t>(dof));

    for (int i = 0; i < dof; ++i) {
        const auto& js = model.joints[static_cast<size_t>(i)];
        SpatialAlgebra::Link link;

        // Parent index
        link.parent = js.parent;

        // Convert Matrix4d → PluckerTransform per Pitfall 3 convention
        link.X = SpatialAlgebra::PluckerTransform(
            SpatialAlgebra::Rotation(js.parentToJoint.topLeftCorner<3, 3>()),
            js.parentToJoint.topRightCorner<3, 1>()
        );

        // Convert mass+com+dense-inertia → RigidBodyInertia per Pitfall 2 (symmetrize)
        Eigen::Matrix3d symInertia = 0.5 * (js.inertia + js.inertia.transpose());
        link.I = SpatialAlgebra::RigidBodyInertia(
            js.mass, js.com,
            SpatialAlgebra::LowerTriangular::fromFullMatrix(symInertia)
        );

        // Convert Vector3d + JointType → MotionVector screw axis
        switch (js.type) {
            case JointType::REVOLUTE:
                link.S = SpatialAlgebra::MotionVector(js.jointAxis, Vector3d::Zero());
                break;
            case JointType::PRISMATIC:
                link.S = SpatialAlgebra::MotionVector(Vector3d::Zero(), js.jointAxis);
                break;
            case JointType::FIXED:
                link.S = SpatialAlgebra::MotionVector(Vector3d::Zero(), Vector3d::Zero());
                break;
        }

        // State initialized to zero
        link.q = 0.0;
        link.qdot = 0.0;
        link.f = SpatialAlgebra::ForceVector(Vector3d::Zero(), Vector3d::Zero());

        fd.links.push_back(std::move(link));

        // Store metadata for forward kinematics
        jointTypes.push_back(js.type);
        jointAxes.push_back(js.jointAxis);
    }

    return fd;
}

/**
 * @brief Build an InverseDynamics solver from a RobotModel
 * @param model Kinematic chain description
 * @param dof Number of degrees of freedom
 * @return Fully constructed InverseDynamics solver object
 * @details Same conversion pattern as buildFD() but using InverseDynamicsLink struct.
 */
SpatialAlgebra::InverseDynamics buildID(const RobotModel& model, int dof) {
    SpatialAlgebra::InverseDynamics id;
    id.links.reserve(static_cast<size_t>(dof));

    for (int i = 0; i < dof; ++i) {
        const auto& js = model.joints[static_cast<size_t>(i)];
        SpatialAlgebra::InverseDynamicsLink link;

        link.parent = js.parent;

        // Matrix4d → PluckerTransform
        link.X = SpatialAlgebra::PluckerTransform(
            SpatialAlgebra::Rotation(js.parentToJoint.topLeftCorner<3, 3>()),
            js.parentToJoint.topRightCorner<3, 1>()
        );

        // mass+com+dense-inertia → RigidBodyInertia with symmetrization
        Eigen::Matrix3d symInertia = 0.5 * (js.inertia + js.inertia.transpose());
        link.I = SpatialAlgebra::RigidBodyInertia(
            js.mass, js.com,
            SpatialAlgebra::LowerTriangular::fromFullMatrix(symInertia)
        );

        // Joint axis mapping
        switch (js.type) {
            case JointType::REVOLUTE:
                link.S = SpatialAlgebra::MotionVector(js.jointAxis, Vector3d::Zero());
                break;
            case JointType::PRISMATIC:
                link.S = SpatialAlgebra::MotionVector(Vector3d::Zero(), js.jointAxis);
                break;
            case JointType::FIXED:
                link.S = SpatialAlgebra::MotionVector(Vector3d::Zero(), Vector3d::Zero());
                break;
        }

        // State initialized to zero
        link.q = 0.0;
        link.qdot = 0.0;
        link.qddot = 0.0;

        id.links.push_back(std::move(link));
    }

    return id;
}

} // anonymous namespace

// =========================================================================
// Destructor — defined here where Impl is complete (PIMPL requirement)
// =========================================================================

SpatialAlgebraAdapter::~SpatialAlgebraAdapter() = default;

// =========================================================================
// build() factory — per D-05 builder pattern
// =========================================================================

std::unique_ptr<SpatialAlgebraAdapter> SpatialAlgebraAdapter::build(
    const RobotModel& model) {
    int dof = model.getDOF();

    auto adapter = std::unique_ptr<SpatialAlgebraAdapter>(
        new SpatialAlgebraAdapter());
    adapter->dof_ = dof;
    adapter->impl_ = std::make_unique<Impl>();

    // Build both solvers from the same RobotModel
    adapter->impl_->fd_ = buildFD(model, dof,
                                   adapter->impl_->jointTypes_,
                                   adapter->impl_->jointAxes_);
    adapter->impl_->id_ = buildID(model, dof);

    // Initialize FK dirty flag
    adapter->fkDirty_ = true;

    return adapter;
}

// =========================================================================
// Private helper — autoFK() dirty-flag pattern
// =========================================================================

void SpatialAlgebraAdapter::autoFK() {
    if (fkDirty_) {
        forwardKinematics();
        fkDirty_ = false;
    }
}

// =========================================================================
// setState() — with input validation (T-20-01)
// =========================================================================

void SpatialAlgebraAdapter::setState(const Eigen::VectorXd& q,
                                      const Eigen::VectorXd& qdot) {
    if (q.size() != dof_ || qdot.size() != dof_) {
        throw std::invalid_argument(
            "setState: vector size mismatch. Expected " +
            std::to_string(dof_) + ", got q=" +
            std::to_string(q.size()) + " qdot=" +
            std::to_string(qdot.size()));
    }

    // Copy to both FD and ID solvers
    for (int i = 0; i < dof_; ++i) {
        impl_->fd_.links[static_cast<size_t>(i)].q = q[i];
        impl_->fd_.links[static_cast<size_t>(i)].qdot = qdot[i];
        impl_->id_.links[static_cast<size_t>(i)].q = q[i];
        impl_->id_.links[static_cast<size_t>(i)].qdot = qdot[i];
    }

    fkDirty_ = true;
}

// =========================================================================
// computeTorques() — RNEA inverse dynamics
// =========================================================================

Eigen::VectorXd SpatialAlgebraAdapter::computeTorques(
    const Eigen::VectorXd& qddot,
    const Eigen::Vector3d& gravity) {

    if (qddot.size() != dof_) {
        throw std::invalid_argument(
            "computeTorques: qddot size mismatch. Expected " +
            std::to_string(dof_) + ", got " +
            std::to_string(qddot.size()));
    }

    // Copy qddot to ID link structs
    for (int i = 0; i < dof_; ++i) {
        impl_->id_.links[static_cast<size_t>(i)].qddot = qddot[i];
    }

    autoFK();

    // SpatialAlgebra::Vector3d == Eigen::Vector3d — no conversion needed
    return impl_->id_.computeTorques(qddot, gravity);
}

// =========================================================================
// computeAccelerations() — ABA forward dynamics
// =========================================================================

Eigen::VectorXd SpatialAlgebraAdapter::computeAccelerations(
    const Eigen::VectorXd& tau,
    const Eigen::Vector3d& gravity) {

    if (tau.size() != dof_) {
        throw std::invalid_argument(
            "computeAccelerations: tau size mismatch. Expected " +
            std::to_string(dof_) + ", got " +
            std::to_string(tau.size()));
    }

    autoFK();

    impl_->fd_.computeAccelerations(tau, gravity);

    // Extract resulting qddot from FD link structs
    Eigen::VectorXd result(dof_);
    for (int i = 0; i < dof_; ++i) {
        result[i] = impl_->fd_.links[static_cast<size_t>(i)].qddot;
    }
    return result;
}

// =========================================================================
// forwardKinematics() — recompute world-frame transforms from joint angles
// =========================================================================

void SpatialAlgebraAdapter::forwardKinematics() {
    impl_->worldX_.clear();
    impl_->worldX_.reserve(static_cast<size_t>(dof_));

    for (int i = 0; i < dof_; ++i) {
        const auto& link = impl_->fd_.links[static_cast<size_t>(i)];
        JointType jt = impl_->jointTypes_[static_cast<size_t>(i)];
        const Vector3d& axis = impl_->jointAxes_[static_cast<size_t>(i)];
        double qi = link.q;

        // Build joint displacement transform based on joint type
        SpatialAlgebra::Rotation jointRot(Eigen::Matrix3d::Identity());
        Vector3d jointTrans = Vector3d::Zero();

        switch (jt) {
            case JointType::REVOLUTE:
                jointRot = SpatialAlgebra::Rotation(
                    Eigen::AngleAxisd(qi, axis));
                break;
            case JointType::PRISMATIC:
                jointTrans = axis * qi;
                break;
            case JointType::FIXED:
                // Identity transform — no change
                break;
        }

        SpatialAlgebra::PluckerTransform X_joint(jointRot, jointTrans);

        // Compose: X_link = link.X * X_joint(q) — joint displacement first
        SpatialAlgebra::PluckerTransform X_link = link.X.multiply(X_joint);

        if (link.parent < 0) {
            // Base link: world transform is just the link transform
            impl_->worldX_.push_back(X_link);
        } else {
            // Compose: worldX[i] = worldX[parent] * X_link
            const auto& worldParent = impl_->worldX_[
                static_cast<size_t>(link.parent)];
            impl_->worldX_.push_back(worldParent.multiply(X_link));
        }
    }
}

// =========================================================================
// getJointTransform(idx) — with bounds check (T-20-02)
// =========================================================================

Eigen::Matrix4d SpatialAlgebraAdapter::getJointTransform(int idx) const {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "getJointTransform: index " + std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    const_cast<SpatialAlgebraAdapter*>(this)->autoFK();

    const auto& worldX = impl_->worldX_[static_cast<size_t>(idx)];
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = worldX.getRotation();
    T.topRightCorner<3, 1>() = worldX.getTranslation();
    return T;
}

// =========================================================================
// computeMassMatrix() — RNEA-column method
// =========================================================================

Eigen::MatrixXd SpatialAlgebraAdapter::computeMassMatrix() {
    autoFK();

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dof_, dof_);

    // Save original qddot values
    std::vector<double> savedQddot(static_cast<size_t>(dof_));
    for (int i = 0; i < dof_; ++i) {
        savedQddot[static_cast<size_t>(i)] =
            impl_->id_.links[static_cast<size_t>(i)].qddot;
    }

    // For each DOF, set unit acceleration and compute torque column
    for (int j = 0; j < dof_; ++j) {
        // Zero all qddot
        for (int k = 0; k < dof_; ++k) {
            impl_->id_.links[static_cast<size_t>(k)].qddot = 0.0;
        }
        // Set unit acceleration for column j
        impl_->id_.links[static_cast<size_t>(j)].qddot = 1.0;

        Eigen::VectorXd e_j = Eigen::VectorXd::Zero(dof_);
        e_j[j] = 1.0;

        Eigen::VectorXd tau_col = impl_->id_.computeTorques(
            e_j, Vector3d::Zero());
        H.col(j) = tau_col;
    }

    // Restore original qddot values
    for (int i = 0; i < dof_; ++i) {
        impl_->id_.links[static_cast<size_t>(i)].qddot =
            savedQddot[static_cast<size_t>(i)];
    }

    return H;
}

// =========================================================================
// computeGravityTorques() — gravity-only inverse dynamics
// =========================================================================

Eigen::VectorXd SpatialAlgebraAdapter::computeGravityTorques(
    const Eigen::Vector3d& gravity) {

    autoFK();

    Eigen::VectorXd qddotZero = Eigen::VectorXd::Zero(dof_);

    // Zero all qddot in ID links
    for (int i = 0; i < dof_; ++i) {
        impl_->id_.links[static_cast<size_t>(i)].qddot = 0.0;
    }

    return impl_->id_.computeTorques(qddotZero, gravity);
}

// =========================================================================
// computeJointSpaceJacobian(idx) — column-by-column geometric Jacobian
// =========================================================================

Eigen::MatrixXd SpatialAlgebraAdapter::computeJointSpaceJacobian(int idx) {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "computeJointSpaceJacobian: index " + std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    autoFK();

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, dof_);

    for (int j = 0; j < dof_; ++j) {
        // Check if link j is an ancestor of (or equal to) link idx
        bool isAncestor = false;
        int current = idx;
        while (current >= 0) {
            if (current == j) {
                isAncestor = true;
                break;
            }
            current = impl_->fd_.links[
                static_cast<size_t>(current)].parent;
        }
        if (!isAncestor) {
            continue; // Column stays zero
        }

        // Get joint screw axis in world frame
        const auto& S_local = impl_->fd_.links[
            static_cast<size_t>(j)].S;
        const auto& worldX_j = impl_->worldX_[
            static_cast<size_t>(j)];

        SpatialAlgebra::MotionVector S_world =
            worldX_j.transformMotion(S_local);

        // Column j = [angular; linear] of the world-frame screw axis
        J.col(j).head<3>() = S_world.getAngular();
        J.col(j).tail<3>() = S_world.getLinear();
    }

    return J;
}

// =========================================================================
// getLinkCOM(idx) — world-frame center of mass (T-20-02)
// =========================================================================

Eigen::Vector3d SpatialAlgebraAdapter::getLinkCOM(int idx) const {
    if (idx < 0 || idx >= dof_) {
        throw std::out_of_range(
            "getLinkCOM: index " + std::to_string(idx) +
            " out of range [0, " + std::to_string(dof_ - 1) + "]");
    }

    const_cast<SpatialAlgebraAdapter*>(this)->autoFK();

    const auto& rbi = impl_->fd_.links[static_cast<size_t>(idx)].I;
    const auto& worldX = impl_->worldX_[static_cast<size_t>(idx)];

    // Transform local COM to world: worldCOM = R * localCOM + t
    const Eigen::Matrix3d& R = worldX.getRotation();
    const Vector3d& t = worldX.getTranslation();

    return R * rbi.getCom() + t;
}

// =========================================================================
// getDOF()
// =========================================================================

int SpatialAlgebraAdapter::getDOF() const {
    return dof_;
}

} // namespace test_models
