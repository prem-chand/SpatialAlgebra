#pragma once

/**
 * @file consistency.h
 * @brief Test model factories for round-trip RNEA↔ABA consistency tests
 * @details Defines factory functions that produce RobotModel instances
 *          specifically for testing the consistency between inverse dynamics
 *          (RNEA) and forward dynamics (ABA) solvers. Each consistency chain
 *          is parameter-equivalent to a corresponding chain from
 *          inverse_dynamics.h or forward_dynamics.h — same masses, COMs,
 *          inertias, and transforms — but with distinct link names for
 *          clarity when debugging round-trip failures.
 *
 *          The key design principle: construct ONE RobotModel, pass it to
 *          BOTH the ID and FD adapters. No duplicated chain definitions
 *          needed. The round-trip test then runs RNEA→tau→ABA→qddot and
 *          verifies qddot matches the input.
 *
 *          These models directly support the consistency tests in
 *          TestDynamicsConsistency.cpp: SingleLinkRoundTrip,
 *          ThreeLinkSerialChain, ThreeLinkNonZeroCOM (CR-02),
 *          and BranchingYConfiguration.
 *
 * @see inverse_dynamics.h for ID chain model counterparts
 * @see forward_dynamics.h for FD chain model counterparts
 * @see TestDynamicsConsistency.cpp (tests/TestDynamicsConsistency.cpp)
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Single-link chain for round-trip RNEA↔ABA consistency tests.
 *        Parameter-equivalent to makeSingleLinkChain. Pass to both ID and FD
 *        adapters.
 * @details Constructs a 1-DOF robot with a single Z-axis revolute joint,
 *          unit mass, zero COM offset, and identity rotational inertia.
 *          Structurally identical to makeSingleLinkChain() but uses
 *          "roundtrip_link0" as the link name for clarity in round-trip
 *          test output.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="roundtrip_link0"
 *
 *          Test usage: the RoundTripABARNEA and RoundTripRNEAABA tests in
 *          TestDynamicsConsistency.cpp use this 1-link model. Construct the
 *          ID and FD adapters from the same RobotModel instance, compute
 *          torques via RNEA, then recover accelerations via ABA.
 *
 * @return RobotModel containing a 1-link Z-revolute chain for round-trip tests
 */
inline RobotModel makeSingleLinkConsistencyChain()
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
    link.name = "roundtrip_link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Three-link serial chain for round-trip consistency tests.
 *        Parameter-equivalent to makeThreeLinkSerialChain. Used with both ID
 *        and FD adapters.
 * @details Constructs a 3-DOF serial chain with 1m X-axis spacing between
 *          links. All links have unit mass, zero COM offset, and identity
 *          rotational inertia — identical to makeThreeLinkSerialChain() but
 *          with "rt_" prefixed link names.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="rt_link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="rt_link1"
 *          - Link 2: child of 1, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="rt_link2"
 *
 *          Test usage: the ThreeLinkSerialChain consistency test feeds known
 *          accelerations (qddot=[1, 0.5, 0.25]) through RNEA to get torques,
 *          then passes the torques to ABA to verify the original accelerations
 *          are recovered within tolerance.
 *
 * @return RobotModel containing a 3-link serial chain for round-trip tests
 */
inline RobotModel makeThreeLinkConsistencyChain()
{
    RobotModel model;

    // Link 0: root
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d::Zero();
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "rt_link0";
    model.joints.push_back(link0);

    // Link 1: child of 0, 1m X translation
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
    link1.name = "rt_link1";
    model.joints.push_back(link1);

    // Link 2: child of 1, 1m X translation
    JointSpec link2;
    link2.parent = 1;

    Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
    T2.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T2.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link2.parentToJoint = T2;
    link2.jointAxis = Eigen::Vector3d::UnitZ();
    link2.type = JointType::REVOLUTE;
    link2.mass = 1.0;
    link2.com = Eigen::Vector3d::Zero();
    link2.inertia = Eigen::Matrix3d::Identity();
    link2.name = "rt_link2";
    model.joints.push_back(link2);

    return model;
}

/**
 * @brief Three-link chain with non-zero COM for round-trip consistency.
 *        Matches CR-02 test model. COM [0.1,0,0] on each link.
 * @details Constructs a 3-DOF serial chain where each link has its center
 *          of mass displaced by 0.1m along the X axis — identical to
 *          makeThreeLinkSerialChainNonZeroCOM() but with "nzc_" prefixed
 *          link names for clear identification in consistency test output.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="nzc_link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="nzc_link1"
 *          - Link 2: child of 1, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="nzc_link2"
 *
 *          Test usage: the ThreeLinkSerialChainNonZeroCOM consistency test
 *          (CR-02) exercises the full condensation step of ABA with
 *          non-block-diagonal spatial inertias. A failure here indicates
 *          incomplete articulation in the ABA inward pass.
 *
 * @return RobotModel containing a 3-link chain with non-zero COM for
 *         round-trip tests
 */
inline RobotModel makeThreeLinkNonZeroCOMConsistencyChain()
{
    RobotModel model;

    // Link 0: root with COM offset
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d(0.1, 0.0, 0.0);
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "nzc_link0";
    model.joints.push_back(link0);

    // Link 1: child of 0 with COM offset, 1m X translation
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitZ();
    link1.type = JointType::REVOLUTE;
    link1.mass = 1.0;
    link1.com = Eigen::Vector3d(0.1, 0.0, 0.0);
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "nzc_link1";
    model.joints.push_back(link1);

    // Link 2: child of 1 with COM offset, 1m X translation
    JointSpec link2;
    link2.parent = 1;

    Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
    T2.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T2.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link2.parentToJoint = T2;
    link2.jointAxis = Eigen::Vector3d::UnitZ();
    link2.type = JointType::REVOLUTE;
    link2.mass = 1.0;
    link2.com = Eigen::Vector3d(0.1, 0.0, 0.0);
    link2.inertia = Eigen::Matrix3d::Identity();
    link2.name = "nzc_link2";
    model.joints.push_back(link2);

    return model;
}

/**
 * @brief Branching Y chain for round-trip consistency tests.
 *        Parameter-equivalent to makeBranchingYConfiguration.
 * @details Constructs a 3-DOF branching kinematic tree where a single base
 *          link has two children symmetrically placed at +1m and -1m along
 *          the X axis. Both children share parent=0. Structurally identical
 *          to makeBranchingYConfiguration() but with "y_" prefixed names
 *          for clarity in round-trip test diagnostics.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="y_base"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([+1,0,0]),
 *            name="y_branch1"
 *          - Link 2: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([-1,0,0]),
 *            name="y_branch2"
 *
 *          Test usage: verifies that branching tree RNEA↔ABA round-trip is
 *          correct — the inward pass must correctly accumulate forces from
 *          multiple children. Symmetric branches with equal accelerations
 *          should produce identical torques and accelerations.
 *
 * @return RobotModel containing a 3-link branching Y chain for round-trip
 *         tests
 */
inline RobotModel makeBranchingYConsistencyChain()
{
    RobotModel model;

    // Link 0: base
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d::Zero();
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "y_base";
    model.joints.push_back(link0);

    // Link 1: first branch, +1m X
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
    link1.name = "y_branch1";
    model.joints.push_back(link1);

    // Link 2: second branch (same parent as link1), -1m X
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
    link2.name = "y_branch2";
    model.joints.push_back(link2);

    return model;
}

} // namespace test_models
