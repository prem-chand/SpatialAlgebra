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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
