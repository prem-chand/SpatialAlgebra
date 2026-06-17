#pragma once

/**
 * @file plucker_transforms.h
 * @brief Test model factories for Plücker coordinate transform tests
 * @details Defines factory functions that produce RobotModel instances with
 *          specific parentToJoint homogeneous transforms for testing
 *          PluckerTransform operations. Each function constructs a single-link
 *          chain whose parentToJoint matrix encodes a known transform
 *          (identity, rotation, translation, or combined) for verification
 *          against expected analytical results.
 *
 *          Plücker transforms use the 6×6 formulation:
 *          @f$ X = [R, 0; -R[t]_\times, R] @f$ for motion and
 *          @f$ X^* = [R^T, [t]_\times R^T; 0, R^T] @f$ for force.
 *          The 4×4 homogeneous matrix in parentToJoint provides the R and t
 *          components that the adapter converts to PluckerTransform.
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms, Ch. 2.
 * @see PluckerTransform (include/PluckerTransform.h)
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Identity transform (R=I, t=0) for PluckerTransform identity tests
 * @details Constructs a single-link chain with an identity parentToJoint
 *          transform. The identity Plücker transform leaves spatial vectors
 *          unchanged — the angular and linear components are preserved.
 *
 *          Chain layout:
 *          - Link 0: root link, parentToJoint=I_4×4, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify that apply() with identity returns the same
 *          MotionVector/ForceVector.
 *
 * @return RobotModel containing a 1-link chain with identity transform
 */
inline RobotModel makeIdentityTransform()
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
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief 90-degree Z-axis rotation transform for Plücker transform tests
 * @details Constructs a single-link chain whose parentToJoint encodes a
 *          90-degree rotation about the Z axis with zero translation.
 *          The rotation matrix is:
 *          @f$ R_z(90°) = [[0, -1, 0], [1, 0, 0], [0, 0, 1]] @f$
 *
 *          Chain layout:
 *          - Link 0: root link, 90° Z rotation, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify that angular components rotate correctly
 *          (e.g., X → Y, Y → -X) and linear components account for
 *          the rotational frame change.
 *
 * @return RobotModel containing a 1-link chain with 90° Z rotation transform
 *
 * @note Uses Eigen::AngleAxisd for rotation construction — matches existing
 *       test fixture patterns in TestPluckerTransform.cpp.
 */
inline RobotModel make90DegreeZRotation()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    // Build 4×4 homogeneous transform: 90° Z rotation, zero translation
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

/**
 * @brief Pure [1,0,0] translation transform with identity rotation
 * @details Constructs a single-link chain whose parentToJoint encodes a
 *          1-meter translation along the X axis with no rotation.
 *          The homogeneous transform is:
 *          @f$ T = [[I | [1,0,0]^T]; [0 0 0 | 1]] @f$
 *
 *          Chain layout:
 *          - Link 0: root link, 1m X translation, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify that the Plücker transform translation
 *          component correctly accounts for the linear offset via the
 *          skew-symmetric cross-product term @f$ -R[t]_\times @f$.
 *
 * @return RobotModel containing a 1-link chain with pure X translation
 */
inline RobotModel makePureTranslation()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    // Build 4×4 homogeneous transform: identity rotation, [1,0,0] translation
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();
    T.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

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
 * @brief Combined 90-degree Z rotation + [1,0,0] translation
 * @details Constructs a single-link chain whose parentToJoint encodes the
 *          composition of a 90° Z rotation followed by a 1-meter X-axis
 *          translation. This is the most general case for testing full
 *          Plücker transform composition.
 *
 *          The homogeneous transform is:
 *          @f$ T = [[R_z(90°) | [1,0,0]^T]; [0 0 0 | 1]] @f$
 *
 *          Chain layout:
 *          - Link 0: root link, 90° Z rot + 1m X trans, Z-revolute,
 *            mass=1.0, COM=origin, identity inertia
 *
 *          Test usage: verify that the Plücker transform correctly composes
 *          rotation and translation — the linear component should be
 *          @f$ -R[t]_\times @f$. Compare individual rotation and translation
 *          contributions against analytical composition.
 *
 * @return RobotModel containing a 1-link chain with combined transform
 *
 * @see make90DegreeZRotation() for rotation-only variant
 * @see makePureTranslation() for translation-only variant
 */
inline RobotModel makeCombinedTransform()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;

    // Build 4×4 homogeneous transform: 90° Z rotation + [1,0,0] translation
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() =
        Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()).matrix();
    T.topRightCorner<3, 1>() = Eigen::Vector3d(1.0, 0.0, 0.0);

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
