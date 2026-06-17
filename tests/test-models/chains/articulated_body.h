#pragma once

/**
 * @file articulated_body.h
 * @brief Test model factories for articulated body inertia tests
 * @details Defines factory functions that produce RobotModel instances for
 *          testing ArticulatedBodyInertia operations. An articulated body
 *          inertia generalizes rigid body inertia for multi-body systems:
 *          @f$ I^A = [[I, H]; [H^T, M]] @f$ where:
 *          - I is the 3×3 rotational inertia
 *          - H is the 3×3 coupling matrix
 *          - M is the 3×3 apparent mass matrix
 *
 *          When H=0 and M=I, the ABI reduces to the equivalent RBI —
 *          this is the key reduction property tested by these models.
 *
 *          Each function constructs a single-link chain whose mass, com,
 *          and inertia fields parameterize the articulated body.
 *
 * @see ArticulatedBodyInertia (include/ArticulatedBodyInertia.h)
 * @see RigidBodyInertia (include/RigidBodyInertia.h)
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms, Ch. 7.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Identity-parameter ABI baseline for reduction and behavior tests
 * @details Constructs a single-link chain with unit mass, zero COM offset,
 *          and identity rotational inertia. This configuration produces
 *          an articulated body inertia equivalent to the identity RBI:
 *          I_A = [[I_3, 0]; [0, I_3]].
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, name="abi_identity"
 *
 *          Test usage: verify that apply() on a MotionVector produces
 *          the same result as RBI::apply() for the same parameters,
 *          and that the reduction to RBI (H=0, M=I case) is exact.
 *
 * @return RobotModel containing a 1-link chain with identity ABI
 *
 * @see makeABIReduced() for the H=0 reduction test variant
 */
inline RobotModel makeABIIdentity()
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
    link.name = "abi_identity";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief ABI configuration matching the reduced case (H=0, M=I)
 * @details Constructs a single-link chain with the same physical
 *          parameters as makeABIIdentity() but intended for ABI→RBI
 *          reduction tests. When the coupling matrix H is zero and
 *          the apparent mass M equals the identity, the articulated
 *          body inertia should reduce to the equivalent rigid body
 *          inertia.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, name="abi_reduced"
 *
 *          Test usage: verify that ABI::reduceToRBI() returns the
 *          correct RigidBodyInertia matching the input parameters,
 *          and that the reduced RBI produces identical apply()
 *          results as the full ABI for this H=0 configuration.
 *
 * @return RobotModel containing a 1-link chain with reduced ABI parameters
 *
 * @see makeABIIdentity() for the equivalent baseline model
 */
inline RobotModel makeABIReduced()
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
    link.name = "abi_reduced";
    model.joints.push_back(link);
    return model;
}

} // namespace test_models
