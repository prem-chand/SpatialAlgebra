#pragma once

/**
 * @file spatial_vectors.h
 * @brief Test model factories for spatial vector transform tests
 * @details Defines factory functions that produce RobotModel instances for
 *          spatial vector (MotionVector/ForceVector) testing. Each function
 *          constructs a valid kinematic chain description using only Eigen3
 *          types, suitable for cross-solver validation.
 *
 *          All models follow Featherstone's spatial vector algebra conventions
 *          (Featherstone 2008). Models are solver-agnostic — the same
 *          RobotModel can be consumed by any solver adapter.
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Single-link Z-revolute chain for basic spatial vector transform tests
 * @details Constructs a 1-DOF kinematic chain with a single link connected
 *          to the world frame (-1 parent) through a Z-axis revolute joint.
 *          The model uses identity transform, unit mass, and zero COM offset
 *          — the simplest non-trivial spatial vector configuration.
 *
 *          Matches Featherstone Example 2.1 pattern: a single rigid body
 *          with spatial velocity [omega; v] expressed in Plücker coordinates.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=origin,
 *            identity inertia
 *
 * @return RobotModel containing a 1-link chain (1 DOF)
 *
 * @note This is the minimal valid model for testing spatial vector operations
 *       like cross-product, dot-product, and addition/subtraction.
 */
inline RobotModel makeSingleLink()
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

} // namespace test_models
