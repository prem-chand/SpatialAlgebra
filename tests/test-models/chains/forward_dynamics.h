#pragma once

/**
 * @file forward_dynamics.h
 * @brief Test model factories for forward dynamics (ABA) tests
 * @details Defines factory functions that produce RobotModel instances for
 *          testing the Articulated Body Algorithm (ABA) forward dynamics
 *          solver. Covers serial chains, non-zero COM, rotated transforms,
 *          and gravity-specific configurations.
 *
 *          All chains use mass=1.0, inertia=Identity_3×3, and Z-axis
 *          revolute joints. Link spacing is 1m along the X axis (except
 *          where rotation alters the frame).
 *
 *          These models are designed to be used alongside the inverse
 *          dynamics models from inverse_dynamics.h — the same RobotModel
 *          can be passed to both ID and FD adapters for round-trip
 *          consistency verification.
 *
 * @see TestForwardDynamics.cpp (tests/TestForwardDynamics.cpp)
 * @see inverse_dynamics.h for ID-compatible chain models
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Three-link serial Z-revolute chain with 1m X spacing.
 *        Matches TestForwardDynamics ThreeLinkNumericalValidation model.
 * @details Constructs a 3-DOF serial chain with 1 meter X-axis spacing
 *          between successive links. All links have unit mass, zero COM
 *          offset, and identity rotational inertia.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link1"
 *          - Link 2: child of 1, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link2"
 *
 *          Test usage: verify ABA round-trip with RNEA — feed known
 *          qddot through RNEA to get tau, then feed tau through ABA to
 *          recover qddot. Also tests that base link carries reflected
 *          inertias from all descendants (CondensationReducesInertiaNorm).
 *
 * @return RobotModel containing a 3-link serial chain
 */
inline RobotModel makeThreeLinkSerialChain()
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

    // Link 2: child of link 1, 1m X translation
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
    link2.name = "link2";
    model.joints.push_back(link2);

    return model;
}

/**
 * @brief Three-link serial chain with non-zero COM [0.1,0,0] on each link.
 *        Matches CR-02 non-zero COM test model from TestForwardDynamics.
 * @details Constructs a 3-DOF serial chain where each link has its center
 *          of mass displaced by 0.1m along the X axis. Non-zero COM produces
 *          off-diagonal blocks in the spatial inertia matrix, which exercises
 *          the full condensation step of the ABA inward pass (Phase 14 CR-02
 *          bug fix).
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link1"
 *          - Link 2: child of 1, Z-revolute, mass=1.0, COM=[0.1,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="link2"
 *
 *          Test usage: RNEA→ABA round-trip with non-zero COM. The ThreeLink-
 *          NumericalValidation test in TestForwardDynamics.cpp uses this
 *          exact model (with COM=[0.1,0,0]) to verify round-trip accuracy.
 *
 * @return RobotModel containing a 3-link chain with non-zero COM offsets
 */
inline RobotModel makeThreeLinkSerialChainNonZeroCOM()
{
    RobotModel model;

    // Link 0: base link with COM offset at [0.1, 0, 0]
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d(0.1, 0.0, 0.0);
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
    link1.com = Eigen::Vector3d(0.1, 0.0, 0.0);
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "link1";
    model.joints.push_back(link1);

    // Link 2: child of link 1 with COM offset, 1m X translation
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
    link2.name = "link2";
    model.joints.push_back(link2);

    return model;
}

/**
 * @brief Two-link chain with 90-degree Z rotation on link 1's parent
 *        transform. Tests PluckerTransform usage in forward dynamics.
 *        Matches TestForwardDynamics PluckerTransformUsage test.
 * @details Constructs a 2-DOF serial chain where link 1's parentToJoint
 *          combines a 90° Z-axis rotation with a 1m X-axis translation.
 *          The 90° rotation changes the frame orientation relative to the
 *          parent, which tests that the Plücker transform correctly handles
 *          rotated link frames in the ABA outward and inward passes.
 *
 *          The combined 4×4 homogeneous transform is:
 *          @f$ T = [[R_z(90°) | [1,0,0]^T]; [0 0 0 | 1]] @f$
 *          where R_z(90°) rotates X→Y, Y→-X.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="link0"
 *          - Link 1: child of 0, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3,
 *            parentToJoint=rotateZ(90°)+translate([1,0,0]),
 *            name="link1"
 *
 *          Test usage: verify that ABA produces correct accelerations when
 *          link 1's frame is rotated relative to link 0. For tau[0]=tau[1]=1.0,
 *          the base acceleration qddot[0] should be zero (child reaction
 *          torque exactly cancels applied torque), while qddot[1] is non-zero.
 *
 * @return RobotModel containing a 2-link chain with 90° Z rotation on link 1
 */
inline RobotModel makeTwoLinkWith90DegreeRotation()
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

    // Link 1: child of link 0, 90° Z rotation + 1m X translation
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() =
        Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()).matrix();
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

} // namespace test_models
