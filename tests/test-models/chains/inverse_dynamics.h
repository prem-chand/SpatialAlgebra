#pragma once

/**
 * @file inverse_dynamics.h
 * @brief Test model factories for inverse dynamics (RNEA) tests
 * @details Defines factory functions that produce RobotModel instances for
 *          testing the Recursive Newton-Euler Algorithm (RNEA) inverse
 *          dynamics solver. Covers single-link, serial, branching, and
 *          gravity-specific kinematic configurations.
 *
 *          All chains use mass=1.0, inertia=Identity_3×3, and Z-axis
 *          revolute joints (except makeSingleLinkXAxisWithCOMOffset which
 *          uses an X-axis joint for static gravity tests).
 *
 *          Chain topologies match the patterns from
 *          TestInverseDynamics.cpp: single-link, two-link serial, branching
 *          Y-tree, COM offset, and X-axis gravity configurations.
 *
 * @see TestInverseDynamics.cpp (tests/TestInverseDynamics.cpp)
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Single-link Z-revolute chain. mass=1, COM at origin, identity inertia.
 *        Baseline RNEA test model.
 * @details Constructs a 1-DOF robot with a single Z-axis revolute joint and
 *          unit mass, zero COM offset, and identity rotational inertia.
 *          This is the simplest non-trivial RNEA model: for qddot=1.0 and
 *          zero gravity, the expected torque is tau=1.0 (I_z * qddot).
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="base"
 *
 *          Test usage: verify computeTorques() produces tau=I*alpha for
 *          the simplest case. Used by TestInverseDynamics SingleLinkPendulum
 *          and the consistency round-trip tests.
 *
 * @return RobotModel containing a 1-link Z-revolute chain
 */
inline RobotModel makeSingleLinkChain()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "base";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Two-link serial Z-revolute chain with 1m X spacing. mass=1 on each
 *        link, COM at link origins. Matches TestInverseDynamics two-link test
 *        model.
 * @details Constructs a 2-DOF serial chain with 1 meter X-axis spacing
 *          between links. Both links have unit mass, zero COM offset, and
 *          identity rotational inertia.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link1"
 *
 *          Test usage: verify velocity propagation and force propagation
 *          across a multi-link chain. Both joints should produce non-zero
 *          torques for non-zero accelerations.
 *
 * @return RobotModel containing a 2-link serial chain with 1m X spacing
 */
inline RobotModel makeTwoLinkSerialChain()
{
    RobotModel model;

    // Link 0: base link at origin
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d::Zero();
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "link0";
    model.joints.push_back(link0);

    // Link 1: child of link 0, 1m X translation
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitZ();
    link1.type = JointType::REVOLUTE;
    link1.mass = 1.0;
    link1.com = Eigen::Vector3d::Zero();
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "link1";
    model.joints.push_back(link1);

    return model;
}

/**
 * @brief Branching Y configuration: base link with two symmetric children at
 *        +/- 1m X. Matches TestInverseDynamics branching test model. Both
 *        children share parent=0.
 * @details Constructs a 3-DOF branching kinematic tree where a single base
 *          link has two children symmetrically placed at +1m and -1m along
 *          the X axis. Both children have parent=0 (not serial i-1).
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="base"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([+1,0,0]),
 *            name="branch1"
 *          - Link 2: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([-1,0,0]),
 *            name="branch2"
 *
 *          Test usage: verify that the inward pass correctly accumulates
 *          forces from multiple children. Symmetric branches with equal
 *          parameters and accelerations should produce identical joint
 *          torques (tau[1] == tau[2]).
 *
 * @return RobotModel containing a 3-link branching Y configuration
 */
inline RobotModel makeBranchingYConfiguration()
{
    RobotModel model;

    // Link 0: base link at origin
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d::Zero();
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "base";
    model.joints.push_back(link0);

    // Link 1: first branch, child of base, +1m X
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitZ();
    link1.type = JointType::REVOLUTE;
    link1.mass = 1.0;
    link1.com = Eigen::Vector3d::Zero();
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "branch1";
    model.joints.push_back(link1);

    // Link 2: second branch, child of base (same parent as link1), -1m X
    JointSpec link2;
    link2.parent = 0;

    Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
    T2.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T2.topRightCorner<3, 1>() = Eigen::Vector3d(-1.0, 0.0, 0.0);

    link2.parentToJoint = T2;
    link2.jointAxis = Eigen::Vector3d::UnitZ();
    link2.type = JointType::REVOLUTE;
    link2.mass = 1.0;
    link2.com = Eigen::Vector3d::Zero();
    link2.inertia = Eigen::Matrix3d::Identity();
    link2.name = "branch2";
    model.joints.push_back(link2);

    return model;
}

/**
 * @brief Two-link serial chain with COM offset [0,0.1,0] on each link.
 *        Tests RNEA with non-zero COM. Matches TestInverseDynamics velocity
 *        test model.
 * @details Constructs a 2-DOF serial chain where each link has its center
 *          of mass displaced by 0.1m along the Y axis. Non-zero COM produces
 *          coupling between angular and linear components in the spatial
 *          inertia, which affects both the velocity product terms (Coriolis)
 *          in the outward pass and the force propagation in the inward pass.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0.1,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0.1,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link1"
 *
 *          Test usage: compare RNEA torque with non-zero qdot against the
 *          zero-velocity case. Coriolis effect should produce measurably
 *          different torques when qdot != 0 and COM != 0.
 *
 * @return RobotModel containing a 2-link chain with COM offset on each link
 */
inline RobotModel makeTwoLinkWithCOMOffset()
{
    RobotModel model;

    // Link 0: base link with COM offset
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d(0.0, 0.1, 0.0);
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "link0";
    model.joints.push_back(link0);

    // Link 1: child of link 0 with COM offset, 1m X translation
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitZ();
    link1.type = JointType::REVOLUTE;
    link1.mass = 1.0;
    link1.com = Eigen::Vector3d(0.0, 0.1, 0.0);
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "link1";
    model.joints.push_back(link1);

    return model;
}

/**
 * @brief Single-link X-axis revolute chain with COM at [0,0.5,0].
 *        Tests static gravity proportionality. Matches TestInverseDynamics
 *        X-axis gravity model.
 * @details Constructs a 1-DOF robot with an X-axis revolute joint and a
 *          center of mass displaced 0.5m along the Y axis. The X-axis joint
 *          axis combined with the Y-offset COM creates a lever arm for
 *          gravitational torque when gravity acts along -Z.
 *
 *          Physical derivation for static gravity:
 *          - Joint axis: X, S=[1,0,0; 0,0,0]
 *          - COM = (0, 0.5, 0), gravity = (0, 0, -g)
 *          - Torque: tau = COM × F = (0, 0.5, 0) × (0, 0, -m*g)
 *            = (-0.5*m*g, 0, 0)
 *          - Project onto X joint: tau[0] = -0.5 * m * g
 *          - Invariant: |tau[0]|/(g * 0.5) = m = 1.0
 *
 *          Chain layout:
 *          - Link 0: root, X-revolute, mass=1.0, COM=[0,0.5,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="base"
 *
 *          Test usage: verify that static torque (qddot=0) is proportional
 *          to gravity magnitude. At g={0, 5, 10}, torque scales linearly
 *          with g and the ratio tau/g is constant.
 *
 * @return RobotModel containing a 1-link X-axis revolute chain with COM offset
 *
 * @note The jointAxis is Eigen::Vector3d::UnitX() — the only function in
 *       this header that does not use the Z-axis default.
 */
inline RobotModel makeSingleLinkXAxisWithCOMOffset()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitX();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d(0.0, 0.5, 0.0);
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "base";
    model.joints.push_back(link);
    return model;
}

} // namespace test_models
