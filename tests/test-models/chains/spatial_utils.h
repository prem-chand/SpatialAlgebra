#pragma once

/**
 * @file spatial_utils.h
 * @brief Test model factories for spatial utility function tests
 * @details Defines factory functions that produce RobotModel instances for
 *          testing spatial utility operations: skew-symmetric matrix
 *          construction, motion vector cross/dot products, and force
 *          vector cross/dot products.
 *
 *          The spatial utilities are free functions (not class methods):
 *          - skew(v): converts 3×1 vector to 3×3 skew-symmetric matrix
 *          - dot(m1, m2): spatial dot product of two MotionVectors
 *          - cross(m1, m2): spatial cross product (Lie bracket)
 *          - These operate on the angular [ω] and linear [v] components
 *            of 6D spatial vectors.
 *
 * @see SpatialUtils (include/SpatialUtils.h)
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Fixture with 90-degree Z rotation for skew-symmetric matrix tests
 * @details Constructs a single-link chain whose parentToJoint encodes a
 *          90° Z-axis rotation. The adapter extracts the rotation matrix
 *          and converts it to skew-symmetric form for testing the skew()
 *          utility function.
 *
 *          Chain layout:
 *          - Link 0: root link, 90° Z rotation, Z-revolute,
 *            mass=1.0, COM=origin, inertia=I_3×3, name="skew_fixture"
 *
 *          Test usage: extract the translation vector from the transform
 *          (or any 3-vector) and verify that skew(v) produces a matrix
 *          with zero diagonal and skew-symmetric off-diagonals, and that
 *          skew(v)*w = v×w (cross product equivalence).
 *
 * @return RobotModel containing a 1-link chain with 90° Z rotation transform
 */
inline RobotModel makeSkewFixture()
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
    link.name = "skew_fixture";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Two-link chain with offset COM for motion vector cross/dot tests
 * @details Constructs a two-link serial chain where link 0 has an offset
 *          COM at [1,0,0], producing non-zero angular-linear coupling in
 *          the spatial velocity. The adapter constructs MotionVector
 *          instances from these links.
 *
 *          Link 0 COM at [1,0,0] means the link's center of mass is offset
 *          1 meter along the X axis from the joint origin. This gives:
 *          @f$ v = ω × c @f$ for pure rotation, producing a non-zero
 *          linear velocity component.
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=1.0, COM=[1,0,0],
 *            inertia=I, name="mv_link0"
 *          - Link 1: child of 0, Z-revolute, 1m X translation,
 *            mass=1.0, COM=origin, inertia=I, name="mv_link1"
 *
 *          Test usage: verify cross(mv1, mv2) anti-commutativity
 *          (cross(a,b) = -cross(b,a)) and dot() distributivity
 *          with non-trivial COM offsets.
 *
 * @return RobotModel containing a 2-link serial chain with offset COM
 */
inline RobotModel makeMotionVectorPair()
{
    RobotModel model;

    // Link 0: root with COM offset at [1,0,0]
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 1.0;
    link0.com = Eigen::Vector3d(1.0, 0.0, 0.0);
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "mv_link0";
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
    link1.name = "mv_link1";
    model.joints.push_back(link1);

    return model;
}

/**
 * @brief Two-link chain for force vector cross/dot product tests
 * @details Constructs a two-link serial chain where both links have mass=2
 *          and COM=[0,1,0], producing non-zero spatial force coupling
 *          for testing force vector operations (cross, dot).
 *
 *          Chain layout:
 *          - Link 0: root, Z-revolute, mass=2.0, COM=[0,1,0],
 *            inertia=I, name="fv_link0"
 *          - Link 1: child of 0, Z-revolute, 1m Y translation,
 *            mass=2.0, COM=[0,1,0], inertia=I, name="fv_link1"
 *
 *          Test usage: verify cross(fv, mv) produces the correct spatial
 *          force transformation (Lie bracket of force and motion), and
 *          that dot(fv, mv) gives the expected power (work rate).
 *
 * @return RobotModel containing a 2-link serial chain for force vector tests
 */
inline RobotModel makeForceVectorPair()
{
    RobotModel model;

    // Link 0: root with mass=2 and COM offset at [0,1,0]
    JointSpec link0;
    link0.parent = -1;
    link0.parentToJoint = Eigen::Matrix4d::Identity();
    link0.jointAxis = Eigen::Vector3d::UnitZ();
    link0.type = JointType::REVOLUTE;
    link0.mass = 2.0;
    link0.com = Eigen::Vector3d(0.0, 1.0, 0.0);
    link0.inertia = Eigen::Matrix3d::Identity();
    link0.name = "fv_link0";
    model.joints.push_back(link0);

    // Link 1: child of link 0, 1m Y translation
    JointSpec link1;
    link1.parent = 0;

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();
    T1.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T1.topRightCorner<3, 1>() = Eigen::Vector3d(0.0, 1.0, 0.0);

    link1.parentToJoint = T1;
    link1.jointAxis = Eigen::Vector3d::UnitZ();
    link1.type = JointType::REVOLUTE;
    link1.mass = 2.0;
    link1.com = Eigen::Vector3d(0.0, 1.0, 0.0);
    link1.inertia = Eigen::Matrix3d::Identity();
    link1.name = "fv_link1";
    model.joints.push_back(link1);

    return model;
}

} // namespace test_models
