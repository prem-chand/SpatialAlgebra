#pragma once

/**
 * @file mujoco_adapter.h
 * @brief MuJoCo adapter implementing the RobotSolver interface
 * @details Wraps MuJoCo's C API (mj_rne, mj_forward) behind the abstract
 *          RobotSolver interface for cross-validation of dynamics.
 *
 *          Builds an MJCF model from a RobotModel description at construction
 *          time, then uses MuJoCo's native solvers for inverse and forward
 *          dynamics computations. State is managed via mjData.
 *
 *          This adapter is optional — it requires MuJoCo 3.4+ C library.
 *
 * @see RobotSolver for the abstract interface
 * @see RobotModel for kinematic chain descriptions
 */

#include "robot_solver.h"
#include "robot_model.h"

#include <Eigen/Dense>
#include <memory>

namespace test_models {

/**
 * @brief MuJoCo adapter implementing RobotSolver
 * @details Wraps the MuJoCo physics engine behind a solver-agnostic interface
 *          for cross-validation of dynamics algorithms. Supports serial and
 *          branching kinematic trees with revolute and prismatic joints.
 *
 *          Construction:
 *          Use the static build(const RobotModel&) factory. The adapter
 *          generates MJCF XML from the RobotModel at construction time,
 *          loads it with mj_loadXML, and initializes mjData.
 *
 *          State Management:
 *          setState() copies joint positions/velocities to mjData and marks
 *          the kinematics cache dirty. The first dynamics query triggers
 *          forwardKinematics() via a dirty-flag pattern.
 *
 *          Limitations:
 *          - No free joints. A floating base must use six serial
 *            1-DOF joints, as buildRobotModelFromMJCF() makes them.
 *          - 1-DOF joints only (matching RobotModel convention)
 *          - Branching trees supported via MuJoCo's natural body hierarchy
 */
class MujocoAdapter : public RobotSolver {
public:
    /**
     * @brief Factory: build adapter from RobotModel
     * @param model Kinematic chain description
     * @return Initialized adapter with MuJoCo model loaded
     * @throws std::runtime_error if MuJoCo model compilation fails
     */
    static std::unique_ptr<MujocoAdapter> build(const RobotModel& model);

    ~MujocoAdapter() override;

    // RobotSolver interface
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
    struct Impl;
    std::unique_ptr<Impl> impl_;

    MujocoAdapter() = default;
};

} // namespace test_models
