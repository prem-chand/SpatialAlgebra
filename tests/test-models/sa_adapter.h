#pragma once

/**
 * @file sa_adapter.h
 * @brief SpatialAlgebra adapter implementing the RobotSolver interface
 * @details This is the sole compilation unit bridging the test-models library
 *          to the SpatialAlgebra solver. All SA types (ForwardDynamics,
 *          InverseDynamics, PluckerTransform, RigidBodyInertia, LowerTriangular)
 *          are hidden behind a PIMPL (pointer-to-implementation) pattern via
 *          std::unique_ptr<Impl>. Consumers of sa_adapter.h see only the
 *          abstract RobotSolver interface.
 *
 *          State Model (per D-04):
 *          RobotModel is the reference kinematic description. The adapter
 *          maintains cached internal state (ForwardDynamics and InverseDynamics
 *          solver objects). setState() copies joint values to internal link
 *          arrays and marks the forward-kinematics cache as dirty via fkDirty_.
 *          The first call to any query method after setState() triggers
 *          autoFK() which runs forwardKinematics() internally.
 *
 *          Builder Pattern (per D-05):
 *          Use the static build(const RobotModel&) factory to construct a
 *          fully initialized adapter from a kinematic chain description.
 *          The same RobotModel can be passed to different solver adapters.
 *
 *          Error Handling:
 *          Computation methods may throw std::invalid_argument (size mismatch)
 *          or std::runtime_error (singular configuration). Exceptions are
 *          propagated directly from the underlying SA solver — the adapter
 *          adds only size-check validation.
 *
 *          This header does NOT include any SpatialAlgebra headers.
 *          Per D-12, sa_adapter.cpp is the ONLY file in tests/test-models/
 *          that includes SpatialAlgebra headers.
 *
 * @see RobotSolver for the abstract interface
 * @see RobotModel for kinematic chain descriptions
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include "robot_solver.h"
#include "robot_model.h"

#include <Eigen/Dense>
#include <memory>

namespace test_models {

/**
 * @brief SpatialAlgebra adapter implementing RobotSolver
 * @details Wraps ForwardDynamics and InverseDynamics solver classes
 *          behind the solver-agnostic RobotSolver interface.
 *          Constructed via build() factory method from a RobotModel
 *          description.
 *
 *          Internal state: maintains ForwardDynamics and InverseDynamics
 *          solver objects as cached internal representations. setState()
 *          updates joint positions and velocities on both solvers.
 *          Computation methods delegate to the appropriate SA solver.
 *
 *          PIMPL pattern: struct Impl is forward-declared here and defined
 *          in sa_adapter.cpp, keeping all SA type dependencies out of the
 *          public header.
 */
class SpatialAlgebraAdapter : public RobotSolver {
public:
    /**
     * @brief Destroy the adapter
     * @details Defined in sa_adapter.cpp where Impl is complete,
     *          required for PIMPL with std::unique_ptr.
     */
    ~SpatialAlgebraAdapter() override;

    /**
     * @brief Factory method per D-05 builder pattern
     * @param model Kinematic chain description
     * @return Fully initialized adapter ready for computation
     * @details Constructs SA solver objects (ForwardDynamics, InverseDynamics)
     *          from the RobotModel description. The adapter is returned as
     *          std::unique_ptr for ownership transfer.
     */
    static std::unique_ptr<SpatialAlgebraAdapter> build(
        const RobotModel& model);

    // -- RobotSolver interface implementation --

    /**
     * @brief Load joint configuration and velocity
     * @details Copies q and qdot to both fd_ and id_ link arrays,
     *          and marks forward kinematics as dirty.
     */
    void setState(const Eigen::VectorXd& q,
                  const Eigen::VectorXd& qdot) override;

    /**
     * @brief Compute joint torques via inverse dynamics (RNEA)
     * @details Delegates to id_.computeTorques() after auto forward
     *          kinematics.
     */
    Eigen::VectorXd computeTorques(
        const Eigen::VectorXd& qddot,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;

    /**
     * @brief Compute joint accelerations via forward dynamics (ABA)
     * @details Delegates to fd_.computeAccelerations() after auto forward
     *          kinematics. Extracts qddot from link structs.
     */
    Eigen::VectorXd computeAccelerations(
        const Eigen::VectorXd& tau,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;

    /**
     * @brief Update internal joint transforms
     * @details Recomputes world-frame transforms from current joint
     *          positions. Called automatically via autoFK() dirty-flag
     *          pattern before any query method.
     */
    void forwardKinematics() override;

    /**
     * @brief Get world-frame transform for a link
     * @details Converts cached PluckerTransform to 4x4 homogeneous matrix.
     */
    Eigen::Matrix4d getJointTransform(int idx) const override;

    /**
     * @brief Compute joint-space inertia matrix H(q)
     * @details Uses RNEA-column method: calls computeTorques() with
     *          unit vectors for each DOF to extract columns of H(q).
     */
    Eigen::MatrixXd computeMassMatrix() override;

    /**
     * @brief Compute gravity compensation torques
     * @details Calls computeTorques with zero acceleration and given gravity.
     */
    Eigen::VectorXd computeGravityTorques(
        const Eigen::Vector3d& gravity) override;

    /**
     * @brief Compute geometric Jacobian for a specific link
     * @details Column-by-column method: each column is the world-frame
     *          joint screw axis for links that affect the queried link.
     */
    Eigen::MatrixXd computeJointSpaceJacobian(int idx) override;

    /**
     * @brief Get link center of mass in world frame
     * @details Transforms local COM (from RigidBodyInertia) to world frame
     *          using the cached world transform.
     */
    Eigen::Vector3d getLinkCOM(int idx) const override;

    /**
     * @brief Number of degrees of freedom
     * @return Cached DOF count from build().
     */
    int getDOF() const override;

private:
    SpatialAlgebraAdapter() = default;  ///< Private; use build() factory

    /**
     * @brief Automatically run forward kinematics if state is dirty
     * @details Called at the start of every query method. If fkDirty_ is
     *          true, calls forwardKinematics() and clears the flag.
     */
    void autoFK();

    struct Impl;                      ///< Forward-declared PIMPL
    std::unique_ptr<Impl> impl_;      ///< Opaque pointer to SA solver objects
    int dof_ = 0;                     ///< Cached DOF count
    bool fkDirty_ = true;             ///< Forward kinematics dirty flag
};

} // namespace test_models
