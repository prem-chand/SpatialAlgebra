#pragma once

/**
 * @file rotation.h
 * @brief Test model factories for rotation matrix tests
 * @details Defines factory functions that produce RobotModel instances with
 *          specific parentToJoint rotation components for testing Rotation
 *          class operations. Each function constructs a single-link chain
 *          encoding a known 3×3 rotation matrix (identity, 45° X, 90° Z)
 *          for verification against expected rotation properties.
 *
 *          Rotation matrices stored in the 3×3 top-left corner of the 4×4
 *          homogeneous parentToJoint transform. The adapter extracts the
 *          rotation via .topLeftCorner<3,3>() for use with the Rotation
 *          class (which extends Eigen::Matrix3d).
 *
 * @see Rotation (include/Rotation.h)
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Zero rotation transform for rotation comparison baseline
 * @details Constructs a single-link chain with identity rotation and zero
 *          translation. Serves as the baseline for rotation comparison tests
 *          where the identity provides a known reference point.
 *
 *          Chain layout:
 *          - Link 0: root link, identity rotation, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify that the extracted rotation is exactly identity,
 *          and that composed rotations produce expected results when compared
 *          against the identity baseline.
 *
 * @return RobotModel containing a 1-link chain with identity rotation
 */
inline RobotModel makeIdentityRotation()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T.topRightCorner<3, 1>() = Eigen::Vector3d::Zero();

    link.parentToJoint = T;
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief 45-degree X-axis rotation for non-axial rotation tests
 * @details Constructs a single-link chain whose parentToJoint encodes a
 *          45-degree rotation about the X axis with zero translation.
 *          The rotation matrix is:
 *          @f$ R_x(45°) = [[1, 0, 0], [0, cos45, -sin45], [0, sin45, cos45]] @f$
 *
 *          Chain layout:
 *          - Link 0: root link, 45° X rotation, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify rotation matrix properties for a non-axial
 *          rotation: determinant is 1.0, columns are orthonormal,
 *          inverse equals transpose, and composed with itself gives a
 *          90° rotation about X.
 *
 * @return RobotModel containing a 1-link chain with 45° X rotation
 *
 * @note This represents a general non-axis-aligned rotation — the rotation
 *       axis is not one of the coordinate axes (it IS the X axis, but the
 *       rotation angle is not 0, 90, or 180 degrees, making the resulting
 *       matrix non-trivial).
 */
inline RobotModel make45DegreeXRotation()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() =
        Eigen::AngleAxisd(M_PI_4, Eigen::Vector3d::UnitX()).matrix();
    T.topRightCorner<3, 1>() = Eigen::Vector3d::Zero();

    link.parentToJoint = T;
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief 90-degree Z rotation fixture for rotation-only tests
 * @details Constructs a single-link chain whose parentToJoint encodes a
 *          90-degree rotation about the Z axis with zero translation.
 *          Matches Plücker transform test patterns for cross-domain
 *          consistency.
 *
 *          The rotation matrix is:
 *          @f$ R_z(90°) = [[0, -1, 0], [1, 0, 0], [0, 0, 1]] @f$
 *
 *          Chain layout:
 *          - Link 0: root link, 90° Z rotation, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify rotation matrix extraction and properties
 *          for a 90° axis-aligned rotation. The "Fixture" suffix
 *          distinguishes this from make90DegreeZRotation() in
 *          plucker_transforms.h — that function targets Plücker transform
 *          tests, while this one targets the Rotation class directly.
 *
 * @return RobotModel containing a 1-link chain with 90° Z rotation
 *
 * @see make90DegreeZRotation() in plucker_transforms.h for the equivalent
 *      transform-domain fixture
 */
inline RobotModel make90DegreeZRotationFixture()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() =
        Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()).matrix();
    T.topRightCorner<3, 1>() = Eigen::Vector3d::Zero();

    link.parentToJoint = T;
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

} // namespace test_models
