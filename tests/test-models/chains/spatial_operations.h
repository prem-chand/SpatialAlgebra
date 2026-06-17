#pragma once

/**
 * @file spatial_operations.h
 * @brief Test model factories for spatial operation tests
 * @details Defines factory functions that produce RobotModel instances for
 *          testing the SpatialOperations static utility class, which provides
 *          crossProductMotion(), crossProductForce(), and transformInertia()
 *          operations on spatial vectors and rigid body inertias.
 *
 *          The spatial operations always preserve the Featherstone convention:
 *          - crossProductMotion (Lie bracket): [ω₁; v₁] × [ω₂; v₂]
 *          - crossProductForce: similar operation on force vectors
 *          - transformInertia: applies a Plücker transform to a
 *            RigidBodyInertia, converting it to a new coordinate frame
 *
 *          Each factory constructs chains with specific joint axes and COM
 *          offsets to produce non-trivial spatial vectors and coupled
 *          inertias for rigorous testing.
 *
 * @see SpatialOperations (include/SpatialOperations.h)
 * @see TestSpatialOperations.cpp (tests/TestSpatialOperations.cpp)
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms, Ch. 2.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Single link with combined 45° X rotation and [1,0,0] translation +
 *        mass=2, COM=[0.5,0,0]. Tests SpatialOperations::transformInertia
 *        behavior.
 * @details Constructs a 1-DOF robot whose parentToJoint encodes a 45° X-axis
 *          rotation combined with a 1m X-axis translation. The link has
 *          mass=2.0 and a COM offset at [0.5, 0, 0], giving non-trivial
 *          rigid body inertia.
 *
 *          The combined 4×4 homogeneous transform is:
 *          @f$ T = [[R_x(45°) | [1,0,0]^T]; [0 0 0 | 1]] @f$
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=2.0, COM=[0.5,0,0],
 *            inertia=I_3×3,
 *            parentToJoint=rotateX(45°)+translate([1,0,0]),
 *            name="transform_rbi"
 *
 *          Test usage: the adapter extracts the Plücker transform from
 *          parentToJoint and the RigidBodyInertia from mass/COM/inertia,
 *          then calls SpatialOperations::transformInertia() to verify
 *          the transformed inertia matches the analytical result.
 *
 * @return RobotModel containing a 1-link chain with non-identity transform
 *         and non-trivial RBI
 */
inline RobotModel makeTransformRBIPair()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    // Build 4×4: 45° X rotation + [1,0,0] translation
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() =
        Eigen::AngleAxisd(M_PI_4, Eigen::Vector3d::UnitX()).matrix();
    T.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link.parentToJoint = T;
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 2.0;
    link.com = Eigen::Vector3d(0.5, 0.0, 0.0);
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "transform_rbi";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Two-link chain with Z and X revolute joints for motion vector cross
 *        product tests. Different joint axes produce non-trivial crossMotion
 *        results.
 * @details Constructs a 2-DOF serial chain where link 0 has a Z-axis revolute
 *          joint and link 1 has an X-axis revolute joint. The different joint
 *          axes produce spatial velocity vectors whose cross product (Lie
 *          bracket) is non-zero in all components, providing comprehensive
 *          coverage for the crossProductMotion() function.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="cross_mv0"
 *          - Link 1: child of 0, X-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="cross_mv1"
 *
 *          Test usage: the adapter constructs MotionVector instances from the
 *          joint screw axes (Z and X). The cross product of these two motion
 *          vectors produces a spatial vector with both angular and linear
 *          components, verifying the full Lie bracket computation.
 *
 * @return RobotModel containing a 2-link chain with Z and X revolute joints
 */
inline RobotModel makeCrossProductMotionFixtures()
{
    RobotModel model;

    // Link 0: Z-axis revolute
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d::Zero();
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "cross_mv0";
    model.joints.push_back(link0);

    // Link 1: X-axis revolute, child of link 0
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitX();
    link1.type = JointType::REVOLUTE;
    link1.mass = 1.0;
    link1.com = Eigen::Vector3d::Zero();
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "cross_mv1";
    model.joints.push_back(link1);

    return model;
}

/**
 * @brief Two-link chain with Z and Y revolute joints and offset COM for force
 *        vector cross product tests. Tests SpatialOperations::crossProductForce
 *        with non-trivial moment arm coupling.
 * @details Constructs a 2-DOF serial chain where link 0 has a Z-axis revolute
 *          joint with COM offset [0.1, 0.2, 0.3] and link 1 has a Y-axis
 *          revolute joint with COM offset [0.3, 0.2, 0.1]. Both links have
 *          mass=2.0. The asymmetric COM offsets create non-trivial moment arm
 *          coupling in the spatial force vectors, producing coupled angular
 *          and linear components in the crossProductForce() result.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=2.0, COM=[0.1,0.2,0.3],
 *            inertia=I_3×3, parentToJoint=I_4×4, name="cross_fv0"
 *          - Link 1: child of 0, Y-revolute, mass=2.0, COM=[0.3,0.2,0.1],
 *            inertia=I_3×3, parentToJoint=translate([1,0,0]),
 *            name="cross_fv1"
 *
 *          Test usage: verify that crossProductForce() correctly computes the
 *          spatial force Lie bracket with asymmetric COM configurations that
 *          couple angular and linear components through the moment arm.
 *
 * @return RobotModel containing a 2-link chain with Z and Y revolute joints
 *         and offset COM
 */
inline RobotModel makeCrossProductForceFixtures()
{
    RobotModel model;

    // Link 0: Z-axis revolute with offset COM
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 2.0;
    link0.com = Eigen::Vector3d(0.1, 0.2, 0.3);
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "cross_fv0";
    model.joints.push_back(link0);

    // Link 1: Y-axis revolute with offset COM, child of link 0
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitY();
    link1.type = JointType::REVOLUTE;
    link1.mass = 2.0;
    link1.com = Eigen::Vector3d(0.3, 0.2, 0.1);
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "cross_fv1";
    model.joints.push_back(link1);

    return model;
}

} // namespace test_models
