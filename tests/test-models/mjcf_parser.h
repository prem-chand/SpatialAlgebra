#pragma once

/**
 * @file mjcf_parser.h
 * @brief Parse MuJoCo MJCF files into RobotModel descriptions
 * @details Loads MJCF XML files through MuJoCo's C API and extracts
 *          kinematic chain data (joints, bodies, inertias) into the
 *          solver-agnostic RobotModel format.
 *
 *          Handles:
 *          - Serial and branching trees
 *          - Fixed, revolute (hinge), and prismatic (slide) joints
 *          - Floating-base robots (a root free joint becomes six virtual joints)
 *          - Body-fixed inertia tensors at COM
 *
 * @see RobotModel for the output format
 * @see https://mujoco.readthedocs.io/en/stable/XMLreference.html
 */

#include "robot_model.h"
#include <mujoco/mujoco.h>
#include <memory>
#include <string>

namespace test_models {

/**
 * @brief Build a RobotModel from an MJCF file loaded through MuJoCo
 * @param xml_path Absolute or relative path to an MJCF XML file
 * @return RobotModel containing 1-DOF joints and explicit zero-DOF fixed bodies
 * @throws std::runtime_error if the file cannot be loaded or parsed
 *
 * The resulting RobotModel uses the MuJoCo body tree directly:
 * - Each body becomes one entry in RobotModel.joints, except the floating base
 * - Additional FREE joints and BALL joints cause an error
 * - Parent indices reference other RobotModel entries (not MuJoCo body IDs)
 * - The chain is in topological order (parents before children)
 * - Armature, damping, friction loss, and joint limits are not kept
 */
RobotModel buildRobotModelFromMJCF(const std::string& xml_path);

} // namespace test_models
