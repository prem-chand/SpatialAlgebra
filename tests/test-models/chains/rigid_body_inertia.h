#pragma once

/**
 * @file rigid_body_inertia.h
 * @brief Test model factories for rigid body inertia tests
 * @details Defines factory functions that produce RobotModel instances with
 *          specific mass, COM, and inertia parameters for testing
 *          RigidBodyInertia operations (apply(), transform(), addition).
 *
 *          Rigid body inertia is the 6×6 spatial inertia matrix:
 *          @f$ I = [[I_c + m[c]_\times[c]_\times^T, m[c]_\times];
 *                   [m[c]_\times^T, mI_3]] @f$
 *          where I_c is the 3×3 rotational inertia at COM and c is the
 *          COM offset vector.
 *
 *          Each function constructs a single-link chain whose mass, com,
 *          and inertia fields define the rigid body. The adapter extracts
 *          these to construct RigidBodyInertia instances.
 *
 * @see RigidBodyInertia (include/RigidBodyInertia.h)
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms, Ch. 2.
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Default rigid body inertia: mass=1, COM at origin, identity inertia
 * @details Constructs a single-link chain with the simplest non-trivial
 *          rigid body inertia: unit mass, zero COM offset, and identity
 *          rotational inertia. This is the baseline for RBI behavior tests.
 *
 *          The resulting spatial inertia matrix is block diagonal:
 *          @f$ I = [[I_3, 0]; [0, I_3]] @f$ (since [c]_\times = 0).
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=[0,0,0],
 *            inertia=I_3×3, name="rbi_default"
 *
 *          Test usage: verify apply() on a MotionVector doubles when mass=1
 *          and inertia is identity, and that transform() preserves the
 *          spatial inertia structure.
 *
 * @return RobotModel containing a 1-link chain with default RBI
 */
inline RobotModel makeRBIDefault()
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
    link.name = "rbi_default";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Configurable COM offset RBI for transform and apply tests
 * @details Constructs a single-link chain with mass=5.0 and a
 *          parameterizable center of mass offset. The COM offset
 *          introduces coupling between angular and linear components
 *          in the spatial inertia matrix via the skew-symmetric
 *          cross-product blocks @f$ m[c]_\times @f$.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=5.0,
 *            COM=[comX,comY,comZ], inertia=I_3×3, name="rbi_offset"
 *
 *          Test usage: verify that the COM offset correctly affects
 *          apply() — a rotation-only motion should produce a non-zero
 *          linear force component proportional to [c]_\times ω.
 *
 * @param comX X-component of center of mass offset (default 1.0)
 * @param comY Y-component of center of mass offset (default 2.0)
 * @param comZ Z-component of center of mass offset (default 3.0)
 *
 * @return RobotModel containing a 1-link chain with offset COM RBI
 */
inline RobotModel makeRBIOffsetCOM(double comX = 1.0,
                                   double comY = 2.0,
                                   double comZ = 3.0)
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 5.0;
    link.com = Eigen::Vector3d(comX, comY, comZ);
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "rbi_offset";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Configurable diagonal inertia RBI for anisotropic tests
 * @details Constructs a single-link chain with configurable mass and
 *          diagonal rotational inertia. The COM is at the origin,
 *          keeping the spatial inertia matrix block-diagonal, but the
 *          rotational component is anisotropic — useful for testing
 *          how different inertia moments affect dynamics.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=mass, COM=[0,0,0],
 *            inertia=diag(ix,iy,iz), name="rbi_diag"
 *
 *          Test usage: verify that different diagonal entries produce
 *          anisotropic responses — a torque about X should produce
 *          different acceleration than the same torque about Z when
 *          ix ≠ iz.
 *
 * @param mass Link mass (default 1.0)
 * @param ix   Rotational inertia about X axis (default 1.0)
 * @param iy   Rotational inertia about Y axis (default 2.0)
 * @param iz   Rotational inertia about Z axis (default 3.0)
 *
 * @return RobotModel containing a 1-link chain with diagonal inertia RBI
 */
inline RobotModel makeRBIDiagonalInertia(double mass = 1.0,
                                          double ix = 1.0,
                                          double iy = 2.0,
                                          double iz = 3.0)
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = mass;
    link.com = Eigen::Vector3d::Zero();

    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    I(0, 0) = ix;
    I(1, 1) = iy;
    I(2, 2) = iz;
    link.inertia = I;
    link.name = "rbi_diag";

    model.joints.push_back(link);
    return model;
}

} // namespace test_models
