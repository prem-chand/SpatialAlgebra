#pragma once

/**
 * @file robot_model.h
 * @brief Pure-data kinematic chain description types using only Eigen3
 * @details Defines JointSpec, RobotModel, and JointType — the zero-dependency
 *          data contracts for the test model library. All types use only Eigen3
 *          and C++ STL types; no SA headers are included.
 *          Designed to be compatible with Featherstone's spatial vector
 *          algebra (Featherstone 2008) while remaining solver-agnostic.
 *
 *          The same RobotModel instance can be passed to multiple solver
 *          adapters (local SA, Pinocchio, RBDL) for cross-validation.
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace test_models {

/**
 * @brief Joint type enumeration
 * @details Defines the three joint types supported by the test model library.
 *          Each type determines how the jointAxis field of JointSpec is
 *          interpreted as a screw axis:
 *          - REVOLUTE: jointAxis is the rotation axis (angular component ω);
 *            the linear component is zero.
 *          - PRISMATIC: jointAxis is the translation direction (linear
 *            component v); the angular component is zero.
 *          - FIXED: rigid 0-DOF connection; jointAxis is ignored.
 */
enum class JointType { REVOLUTE, PRISMATIC, FIXED };

/**
 * @brief Per-link specification for a kinematic chain
 * @details Pure data POD struct describing one rigid body and its connection
 *          to its parent. All fields use only Eigen3 types — no
 *          SA library dependency. Designed to be directly convertible
 *          to solver-specific types (FD/ID link structs,
 *          pinocchio::Model, etc.).
 *
 *          The parentToJoint transform uses a 4×4 homogeneous matrix
 *          [R | t; 0 | 1] convention in the parent frame, where R is the
 *          3×3 rotation and t is the translation vector to the joint origin.
 *          The inertia matrix is a 3×3 dense rotational inertia at the
 *          center of mass (not lower-triangular packed storage).
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */
struct JointSpec {
    int parent = -1; ///< Parent link index (-1 for base link)

    Eigen::Matrix4d parentToJoint =
        Eigen::Matrix4d::Identity(); ///< Homogeneous transform
                                     ///< parent → link, in parent frame
                                     ///< Convention: [R | t; 0 | 1]

    Eigen::Vector3d jointAxis =
        Eigen::Vector3d::UnitZ(); ///< Joint axis direction (unit vector)

    JointType type = JointType::REVOLUTE; ///< Joint type (REVOLUTE,
                                          ///< PRISMATIC, or FIXED)

    double mass = 1.0; ///< Link mass (kg)

    Eigen::Vector3d com =
        Eigen::Vector3d::Zero(); ///< Center of mass in link frame

    Eigen::Matrix3d inertia =
        Eigen::Matrix3d::Identity(); ///< 3×3 dense rotational inertia at
                                     ///< COM (not lower-triangular packed)

    std::string name; ///< Optional link name for debugging and display
};

/**
 * @brief Complete kinematic chain description
 * @details An ordered collection of JointSpec entries describing a robot
 *          model as a kinematic tree. Links must be in topological order:
 *          parent links MUST appear before their children in the joints
 *          vector. For branching trees, multiple children of the same
 *          parent must have consecutive indices.
 *
 *          The degrees of freedom are derived from the number of joints
 *          (1 DOF per joint). This matches all existing SA library
 *          test models which use exclusively 1-DOF joints.
 *
 *          Usage:
 *          @code{.cpp}
 *          RobotModel model;
 *          JointSpec link;
 *          link.parent = -1;
 *          link.parentToJoint = Eigen::Matrix4d::Identity();
 *          link.jointAxis = Eigen::Vector3d::UnitZ();
 *          link.type = JointType::REVOLUTE;
 *          link.mass = 1.0;
 *          link.com = Eigen::Vector3d::Zero();
 *          link.inertia = Eigen::Matrix3d::Identity();
 *          link.name = "base";
 *          model.joints.push_back(link);
 *          int dof = model.getDOF(); // 1
 *          @endcode
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */
struct RobotModel {
    std::vector<JointSpec> joints; ///< Link specifications in topological
                                   ///< order

    /**
     * @brief Number of degrees of freedom
     * @return Total number of joints (1 DOF per joint)
     */
    int getDOF() const { return static_cast<int>(joints.size()); }
};

} // namespace test_models
