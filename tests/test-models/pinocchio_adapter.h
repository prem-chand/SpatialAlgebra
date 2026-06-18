#pragma once

/**
 * @file pinocchio_adapter.h
 * @brief Pinocchio adapter implementing the RobotSolver interface
 * @details Implements the RobotSolver abstract interface using the Pinocchio
 *          C++ library (https://stack-of-tasks.github.io/pinocchio/).
 *          All Pinocchio types are hidden behind a PIMPL pattern so that
 *          consumers see only the abstract RobotSolver interface.
 *
 *          Joint type mapping (PCC-04):
 *          - REVOLUTE with axis=X → pinocchio::JointModelRevoluteUnaligned(X)
 *          - PRISMATIC with axis → pinocchio::JointModelPrismaticUnaligned(axis)
 *          - FIXED → pinocchio::JointModelFixed
 *
 *          State Model (per D-04):
 *          RobotModel is the reference kinematic description. The adapter
 *          maintains a pinocchio::Model + pinocchio::Data pair as cached
 *          internal state. setState() copies joint values and marks forward
 *          kinematics as dirty via fkDirty_. The first query method after
 *          setState() triggers autoFK() which runs forward kinematics.
 *
 *          Builder Pattern (per D-05):
 *          Use the static build(const RobotModel&) factory to construct a
 *          fully initialized adapter from a kinematic chain description.
 *
 *          Gravity convention:
 *          The Pinocchio adapter stores model.gravity as a Motion vector
 *          with linear component equal to the gravity parameter passed to
 *          computeTorques/computeAccelerations and zero angular component.
 *          This is set before each dynamics computation call.
 *
 *          This header does NOT include any Pinocchio headers.
 *          All Pinocchio types are hidden behind std::unique_ptr<Impl>.
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
 * @brief Pinocchio adapter implementing RobotSolver
 * @details Wraps pinocchio::Model, pinocchio::Data, and the Pinocchio
 *          dynamics algorithms (rnea, aba, crba, forwardKinematics,
 *          computeJointJacobians) behind the solver-agnostic RobotSolver
 *          interface. Constructed via build() factory method from a
 *          RobotModel description.
 *
 *          Internal state: pinocchio::Model + pinocchio::Data pair
 *          rebuilt from RobotModel on construction. setState() copies
 *          joint configuration and velocity values. Computation methods
 *          call the appropriate Pinocchio algorithm and extract results.
 *
 *          PIMPL pattern: struct Impl is forward-declared here and
 *          defined in pinocchio_adapter.cpp, keeping all Pinocchio type
 *          dependencies out of the public header.
 */
class PinocchioAdapter : public RobotSolver {
public:
    /**
     * @brief Destroy the adapter
     * @details Defined in pinocchio_adapter.cpp where Impl is complete,
     *          required for PIMPL with std::unique_ptr.
     */
    ~PinocchioAdapter() override;

    /**
     * @brief Factory method per D-05 builder pattern
     * @param model Kinematic chain description
     * @return Fully initialized adapter ready for computation
     * @details Constructs a Pinocchio Model+Data pair from the RobotModel
     *          description. The adapter is returned as std::unique_ptr
     *          for ownership transfer.
     */
    static std::unique_ptr<PinocchioAdapter> build(
        const RobotModel& model);

    // -- RobotSolver interface implementation --

    void setState(const Eigen::VectorXd& q,
                  const Eigen::VectorXd& qdot) override;

    Eigen::VectorXd computeTorques(
        const Eigen::VectorXd& qddot,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;

    Eigen::VectorXd computeAccelerations(
        const Eigen::VectorXd& tau,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;

    void forwardKinematics() override;

    Eigen::Matrix4d getJointTransform(int idx) const override;

    Eigen::MatrixXd computeMassMatrix() override;

    Eigen::VectorXd computeGravityTorques(
        const Eigen::Vector3d& gravity) override;

    Eigen::MatrixXd computeJointSpaceJacobian(int idx) override;

    Eigen::Vector3d getLinkCOM(int idx) const override;

    int getDOF() const override;

private:
    PinocchioAdapter() = default;  ///< Private; use build() factory

    /**
     * @brief Automatically run forward kinematics if state is dirty
     * @details Called at the start of every query method. If fkDirty_ is
     *          true, calls forwardKinematics() and clears the flag.
     */
    void autoFK();

    struct Impl;                      ///< Forward-declared PIMPL
    std::unique_ptr<Impl> impl_;      ///< Opaque pointer to Pinocchio objects
    int dof_ = 0;                     ///< Cached DOF count
    bool fkDirty_ = true;             ///< Forward kinematics dirty flag
};

} // namespace test_models
