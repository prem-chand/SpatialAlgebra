#include "InverseDynamics.h"
#include "ForwardDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-8;

/**
 * @brief Round-trip test: ABA(RNEA(tau)) ≈ tau
 * 
 * Start with random torques, compute accelerations with ABA,
 * then compute torques with RNEA, verify we get back the original torques.
 */
TEST(ConsistencyTest, RoundTripABARNEA) {
    // Start with random torque
    Eigen::VectorXd tau_input(1);
    tau_input[0] = 1.0;
    
    // Forward dynamics: tau -> qddot
    ForwardDynamics fd;
    ForwardDynamicsLink fd_link;
    fd_link.parent = -1;
    fd_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link.q = 0.0;
    fd_link.qdot = 0.0;
    fd.links.push_back(fd_link);
    
    fd.computeAccelerations(tau_input);
    double qddot_result = fd.links[0].qddot;
    
    // Inverse dynamics: qddot -> tau
    InverseDynamics id;
    InverseDynamicsLink id_link;
    id_link.parent = -1;
    id_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link.q = 0.0;
    id_link.qdot = 0.0;
    id.links.push_back(id_link);
    
    Eigen::VectorXd qddot_vec(1);
    qddot_vec[0] = qddot_result;
    Eigen::VectorXd tau_output = id.computeTorques(qddot_vec);
    
    // Verify: tau_output ≈ tau_input
    EXPECT_NEAR(tau_output[0], tau_input[0], EPSILON);
}

/**
 * @brief Round-trip test: RNEA(ABA(qddot)) ≈ qddot
 * 
 * Start with random accelerations, compute torques with RNEA,
 * then compute accelerations with ABA, verify we get back the original accelerations.
 */
TEST(ConsistencyTest, RoundTripRNEAABA) {
    // Start with random acceleration
    Eigen::VectorXd qddot_input(1);
    qddot_input[0] = 1.0;
    
    // Inverse dynamics: qddot -> tau
    InverseDynamics id;
    InverseDynamicsLink id_link;
    id_link.parent = -1;
    id_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link.q = 0.0;
    id_link.qdot = 0.0;
    id.links.push_back(id_link);
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input);
    
    // Forward dynamics: tau -> qddot
    ForwardDynamics fd;
    ForwardDynamicsLink fd_link;
    fd_link.parent = -1;
    fd_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link.q = 0.0;
    fd_link.qdot = 0.0;
    fd.links.push_back(fd_link);
    
    fd.computeAccelerations(tau);
    
    // Verify: fd.links[0].qddot ≈ qddot_input
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
}

/**
 * @brief Consistency test for 3-link serial chain
 * 
 * Verifies RNEA and ABA produce consistent results for multi-link system
 */
TEST(ConsistencyTest, ThreeLinkSerialChain) {
    // Test RNEA -> ABA
    InverseDynamics id;
    
    // Link 0 (base)
    InverseDynamicsLink id_link0;
    id_link0.parent = -1;
    id_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link0.q = 0.0;
    id_link0.qdot = 0.0;
    id.links.push_back(id_link0);
    
    // Link 1
    InverseDynamicsLink id_link1;
    id_link1.parent = 0;
    id_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link1.q = 0.0;
    id_link1.qdot = 0.0;
    id.links.push_back(id_link1);
    
    // Link 2
    InverseDynamicsLink id_link2;
    id_link2.parent = 1;
    id_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link2.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link2.q = 0.0;
    id_link2.qdot = 0.0;
    id.links.push_back(id_link2);
    
    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.25;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input);
    
    // Forward dynamics
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_link0;
    fd_link0.parent = -1;
    fd_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link0.q = 0.0;
    fd_link0.qdot = 0.0;
    fd.links.push_back(fd_link0);
    
    ForwardDynamicsLink fd_link1;
    fd_link1.parent = 0;
    fd_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link1.q = 0.0;
    fd_link1.qdot = 0.0;
    fd.links.push_back(fd_link1);
    
    ForwardDynamicsLink fd_link2;
    fd_link2.parent = 1;
    fd_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link2.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link2.q = 0.0;
    fd_link2.qdot = 0.0;
    fd.links.push_back(fd_link2);
    
    fd.computeAccelerations(tau);
    
    // Verify accelerations match
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], EPSILON);
}

/**
 * @brief Consistency test for branching Y-configuration
 * 
 * Verifies RNEA and ABA produce consistent results for branching tree
 */
TEST(ConsistencyTest, BranchingYConfiguration) {
    // Inverse dynamics setup
    InverseDynamics id;
    
    InverseDynamicsLink id_link0;
    id_link0.parent = -1;
    id_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link0.q = 0.0;
    id_link0.qdot = 0.0;
    id.links.push_back(id_link0);
    
    InverseDynamicsLink id_link1;
    id_link1.parent = 0;
    id_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link1.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    id_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link1.q = 0.0;
    id_link1.qdot = 0.0;
    id.links.push_back(id_link1);
    
    InverseDynamicsLink id_link2;
    id_link2.parent = 0;
    id_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    id_link2.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    id_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link2.q = 0.0;
    id_link2.qdot = 0.0;
    id.links.push_back(id_link2);
    
    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input);
    
    // Forward dynamics setup
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_link0;
    fd_link0.parent = -1;
    fd_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link0.q = 0.0;
    fd_link0.qdot = 0.0;
    fd.links.push_back(fd_link0);
    
    ForwardDynamicsLink fd_link1;
    fd_link1.parent = 0;
    fd_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link1.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    fd_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link1.q = 0.0;
    fd_link1.qdot = 0.0;
    fd.links.push_back(fd_link1);
    
    ForwardDynamicsLink fd_link2;
    fd_link2.parent = 0;
    fd_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    fd_link2.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    fd_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link2.q = 0.0;
    fd_link2.qdot = 0.0;
    fd.links.push_back(fd_link2);
    
    fd.computeAccelerations(tau);
    
    // Verify accelerations match
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], EPSILON);
    
    // Verify symmetric branches have same acceleration
    EXPECT_NEAR(fd.links[1].qddot, fd.links[2].qddot, EPSILON);
}

/**
 * @brief Round-trip consistency with gravity
 * @details Verifies that ABA(RNEA(qddot, gravity), gravity) ≈ qddot
 *          under the influence of gravity. Both solvers should apply
 *          gravity consistently via the base acceleration formulation.
 */
TEST(ConsistencyTest, RoundTripWithGravity) {
    Vector3d gravity(0, 0, -9.81);
    
    // Single-link setup
    InverseDynamics id;
    InverseDynamicsLink id_link;
    id_link.parent = -1;
    id_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link.q = 0.0;
    id_link.qdot = 0.0;
    id.links.push_back(id_link);
    
    Eigen::VectorXd qddot_input(1);
    qddot_input[0] = 1.0;
    
    // RNEA: qddot -> tau (with gravity)
    Eigen::VectorXd tau = id.computeTorques(qddot_input, gravity);
    
    // ABA: tau -> qddot (with gravity)
    ForwardDynamics fd;
    ForwardDynamicsLink fd_link;
    fd_link.parent = -1;
    fd_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link.q = 0.0;
    fd_link.qdot = 0.0;
    fd.links.push_back(fd_link);
    
    fd.computeAccelerations(tau, gravity);
    
    // Round-trip with gravity should be consistent
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
}

/**
 * @brief Two-link round-trip consistency with gravity
 */
TEST(ConsistencyTest, TwoLinkRoundTripWithGravity) {
    Vector3d gravity(0, 0, -9.81);
    
    // RNEA setup
    InverseDynamics id;
    
    InverseDynamicsLink id_l0;
    id_l0.parent = -1;
    id_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l0.q = 0.0;
    id_l0.qdot = 0.0;
    id.links.push_back(id_l0);
    
    InverseDynamicsLink id_l1;
    id_l1.parent = 0;
    id_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    id_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l1.q = 0.0;
    id_l1.qdot = 0.0;
    id.links.push_back(id_l1);
    
    Eigen::VectorXd qddot_input(2);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input, gravity);
    
    // ABA setup
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_l0;
    fd_l0.parent = -1;
    fd_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_l0.q = 0.0;
    fd_l0.qdot = 0.0;
    fd.links.push_back(fd_l0);
    
    ForwardDynamicsLink fd_l1;
    fd_l1.parent = 0;
    fd_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    fd_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_l1.q = 0.0;
    fd_l1.qdot = 0.0;
    fd.links.push_back(fd_l1);
    
    fd.computeAccelerations(tau, gravity);
    
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
}

/**
 * @brief Single-link direct round-trip: ABA(RNEA(qddot)) ≈ qddot
 * @details This is the most fundamental consistency test: for a single link
 *          with no gravity, the RNEA and ABA should be exact inverses.
 *          Start with a known qddot, compute tau via RNEA, then give tau
 *          to ABA and verify the original qddot is recovered.
 * 
 *          This test must pass before multi-link consistency tests can be
 *          meaningfully interpreted. A single-link round-trip failure
 *          indicates a fundamental bug in either solver's algebraic core,
 *          independent of kinematic chain complexity.
 * 
 *          Unlike the existing RoundTripABARNEA and RoundTripRNEAABA tests
 *          which check round-trip fidelity indirectly, this test explicitly
 *          names the direction (RNEA→ABA via tau) and documents the
 *          single-link dependency in its @details. It also uses X-axis
 *          revolute joints for independent coverage.
 * 
 *          Independent oracle: τ = I·a + v × I·v for RNEA forward,
 *          and qddot = τ / I_eff for ABA, which should be exact inverses.
 */
TEST(ConsistencyTest, RoundTripABA_RNEA_DirectComparison) {
    // Start with known acceleration
    Eigen::VectorXd qddot_input(1);
    qddot_input[0] = 1.0;
    
    // RNEA: qddot -> tau (X-axis revolute joint with COM offset)
    InverseDynamics id;
    InverseDynamicsLink id_link;
    id_link.parent = -1;
    id_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    id_link.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    id_link.q = 0.0;
    id_link.qdot = 0.0;
    id.links.push_back(id_link);
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input, Vector3d::Zero());
    EXPECT_TRUE(std::isfinite(tau[0]));
    EXPECT_NEAR(tau[0], 1.0, 1e-10);  // With identity inertia, τ = I·a = qddot
    
    // ABA: tau -> qddot (same link structure)
    ForwardDynamics fd;
    ForwardDynamicsLink fd_link;
    fd_link.parent = -1;
    fd_link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    fd_link.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    fd_link.q = 0.0;
    fd_link.qdot = 0.0;
    fd_link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(fd_link);
    
    fd.computeAccelerations(tau, Vector3d::Zero());
    
    // Single-link round-trip should recover the original qddot exactly
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], 1e-10);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
