#pragma once

/**
 * @file robot_solver.h
 * @brief Abstract solver interface for robotics dynamics
 * @details Defines RobotSolver — a pure virtual base class providing a
 *          solver-agnostic API for forward and inverse dynamics, kinematics,
 *          and related computations. This is the abstract contract that all
 *          solver adapters implement.
 *
 *          The interface uses only dynamic Eigen types (VectorXd, MatrixXd)
 *          for joint-space quantities — no templates on scalar type or DOF
 *          count. This keeps the interface simple and avoids template
 *          proliferation across solver implementations.
 *
 *          State management follows a hybrid model:
 *          - setState() loads joint configuration and velocity, marking
 *            internal cached state as dirty.
 *          - Adapter implementations use a dirty-flag pattern: the first
 *            query method after setState() triggers forwardKinematics()
 *            internally.
 *          - Computation methods read from the cached internal state.
 *
 *          Each adapter (SA, Pinocchio, RBDL) provides its own static
 *          build(const RobotModel&) factory method returning a
 *          std::unique_ptr to the concrete adapter type. The RobotSolver
 *          base class does NOT declare build() — the return type varies
 *          per adapter.
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include "robot_model.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace test_models {

/**
 * @brief Abstract solver interface for robotics dynamics
 * @details Pure virtual interface enabling solver-agnostic test models.
 *          Implementations wrap specific solvers (local SA, Pinocchio,
 *          RBDL) behind a uniform API. All joint-space vectors use
 *          dynamic Eigen types (VectorXd, MatrixXd).
 *
 *          State management:
 *          Call setState(q, qdot) before querying dynamics quantities.
 *          Computation methods read from the cached internal state.
 *          Forward kinematics is triggered automatically by adapter
 *          implementations via a dirty-flag pattern when the first
 *          query method is called after a state change.
 *
 *          Error handling:
 *          Methods may throw std::invalid_argument on size mismatch
 *          between input vectors and getDOF(). Methods may throw
 *          std::out_of_range on invalid link indices.
 *
 * @see RobotModel for kinematic chain description
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */
class RobotSolver {
public:
    virtual ~RobotSolver() = default;

    /**
     * @brief Load joint configuration and velocity
     * @param q Joint positions (size must equal getDOF())
     * @param qdot Joint velocities (size must equal getDOF())
     * @details Copies the provided state and marks internal cache as
     *          dirty. Adapter implementations rebuild internal
     *          representations on first subsequent query method call.
     *          This follows the hybrid state model: RobotModel holds
     *          reference configuration; adapter maintains cached
     *          internal state.
     */
    virtual void setState(const Eigen::VectorXd& q,
                          const Eigen::VectorXd& qdot) = 0;

    /**
     * @brief Compute joint torques via inverse dynamics (RNEA)
     * @param qddot Joint accelerations (size must equal getDOF())
     * @param gravity Gravity vector in world frame (default: zero)
     * @return Joint torques tau (size = getDOF())
     * @throws std::invalid_argument if qddot.size() != getDOF()
     */
    virtual Eigen::VectorXd computeTorques(
        const Eigen::VectorXd& qddot,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) = 0;

    /**
     * @brief Compute joint accelerations via forward dynamics (ABA)
     * @param tau Joint torques (size must equal getDOF())
     * @param gravity Gravity vector in world frame (default: zero)
     * @return Joint accelerations qddot (size = getDOF())
     * @throws std::invalid_argument if tau.size() != getDOF()
     * @throws std::runtime_error on singular configuration
     */
    virtual Eigen::VectorXd computeAccelerations(
        const Eigen::VectorXd& tau,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) = 0;

    /**
     * @brief Update internal joint transforms
     * @details Recomputes link-to-world transforms from current joint
     *          positions. Called automatically by adapter
     *          implementations when state is dirty (dirty-flag
     *          pattern). Explicit calls are safe but typically
     *          unnecessary.
     */
    virtual void forwardKinematics() = 0;

    /**
     * @brief Get world-frame transform for a link
     * @param idx Link index (0 to getDOF()-1)
     * @return 4×4 homogeneous transform [R | t; 0 | 1] in world frame
     * @throws std::out_of_range if idx is out of bounds
     */
    virtual Eigen::Matrix4d getJointTransform(int idx) const = 0;

    /**
     * @brief Compute joint-space inertia matrix H(q)
     * @return DOF × DOF symmetric positive-definite matrix
     * @details The mass matrix maps joint accelerations to joint
     *          torques: tau = H(q) * qddot + C(q, qdot) + G(q).
     *          Adapter implementations typically use the RNEA-column
     *          method: call computeTorques() with unit vectors.
     */
    virtual Eigen::MatrixXd computeMassMatrix() = 0;

    /**
     * @brief Compute gravity compensation torques
     * @param gravity Gravity vector in world frame
     * @return Static joint torques (size = getDOF()) for zero velocity
     *         and zero acceleration
     * @details Equivalent to computeTorques(zero_qddot, gravity) but
     *          may be optimized by adapter implementations to skip
     *          velocity-dependent terms.
     */
    virtual Eigen::VectorXd computeGravityTorques(
        const Eigen::Vector3d& gravity) = 0;

    /**
     * @brief Compute geometric Jacobian for a specific link
     * @param idx Link index (0 to getDOF()-1)
     * @return 6 × DOF Jacobian matrix mapping joint velocities to
     *         spatial velocity of link idx in world frame:
     *         v_spatial = J * qdot
     * @throws std::out_of_range if idx is out of bounds
     * @details Column j is the world-frame joint screw axis for
     *          link idx. Columns for joints that do not affect
     *          link idx are zero. Uses cached link transforms
     *          from forwardKinematics().
     */
    virtual Eigen::MatrixXd computeJointSpaceJacobian(int idx) = 0;

    /**
     * @brief Get link center of mass in world frame
     * @param idx Link index (0 to getDOF()-1)
     * @return World-frame COM position (3D vector)
     * @throws std::out_of_range if idx is out of bounds
     */
    virtual Eigen::Vector3d getLinkCOM(int idx) const = 0;

    /**
     * @brief Number of degrees of freedom
     * @return Total joint DOF count (equals number of joints * DOF per joint)
     */
    virtual int getDOF() const = 0;
};

} // namespace test_models
