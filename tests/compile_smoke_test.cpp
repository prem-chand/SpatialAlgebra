/**
 * @brief Compile smoke test: verifies test-models library is zero-dependency
 * @details This file tests that robot_model.h and robot_solver.h can be included
 *          and used WITHOUT any SpatialAlgebra headers. If this compiles, the
 *          test-models library satisfies TML-01 (zero SA dependency guarantee).
 *
 *          Also verifies basic type instantiation: RobotModel, JointSpec,
 *          JointType enum, and RobotSolver (abstract class — no construction).
 *          Includes a chains/ header (inverse_dynamics.h) to verify the full
 *          include chain compiles without SpatialAlgebra.
 *
 *          NOTE: SpatialAlgebra.h is NOT included. This is intentional — the
 *          entire point of the smoke test is to prove the compilation firewall.
 */
#include "robot_model.h"
#include "robot_solver.h"
#include "chains/inverse_dynamics.h"
#include <Eigen/Dense>
#include <cassert>
#include <iostream>

int main() {
    using namespace test_models;

    // Verify JointType enum compiles and is usable
    JointType jt = JointType::REVOLUTE;
    assert(jt == JointType::REVOLUTE);
    (void)JointType::PRISMATIC;
    (void)JointType::FIXED;

    // Verify JointSpec struct compiles with all fields
    JointSpec js;
    js.parent = -1;
    js.parentToJoint = Eigen::Matrix4d::Identity();
    js.jointAxis = Eigen::Vector3d::UnitZ();
    js.type = JointType::REVOLUTE;
    js.mass = 1.0;
    js.com = Eigen::Vector3d::Zero();
    js.inertia = Eigen::Matrix3d::Identity();
    js.name = "test_link";

    // Verify RobotModel compiles: push_back and getDOF()
    RobotModel model;
    model.joints.push_back(js);
    assert(model.getDOF() == 1);

    // Verify RobotSolver is an abstract class (cannot instantiate)
    // RobotSolver solver; // would not compile — proves it's abstract

    // Verify factory function from chains/ compiles and produces valid model
    RobotModel chain = makeSingleLinkChain();
    assert(chain.getDOF() == 1);
    assert(chain.joints[0].type == JointType::REVOLUTE);
    assert(chain.joints[0].name == "base");

    // Verify a multi-link chain factory also works
    RobotModel twoLink = makeTwoLinkSerialChain();
    assert(twoLink.getDOF() == 2);
    assert(twoLink.joints[0].parent == -1);  // base link
    assert(twoLink.joints[1].parent == 0);   // child of base

    std::cout << "compile_smoke_test: PASSED — test-models library compiles "
              << "with zero SpatialAlgebra dependency" << std::endl;
    return 0;
}
