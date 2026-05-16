#include "InverseDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-10;

/**
 * @brief Test RNEA on single-link pendulum
 * 
 * Verifies: τ = I*α for simple case
 * Single revolute joint around Z-axis with known inertia
 * Expected: tau = I_z * qddot
 */
TEST(InverseDynamicsTest, SingleLinkPendulum) {
    InverseDynamics id;
    
    // Setup: single revolute joint around Z
    InverseDynamicsLink link;
    link.parent = -1;  // Base
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());  // Revolute Z
    link.q = 0.0;
    link.qdot = 0.0;
    
    id.links.push_back(link);
    
    // Apply acceleration: qddot = 1.0 rad/s²
    Eigen::VectorXd qddot(1);
    qddot[0] = 1.0;
    
    // Execute
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    // Verify: tau = I_z * qddot = 1.0 * 1.0 = 1.0 Nm
    EXPECT_NEAR(tau[0], 1.0, EPSILON);
}

/**
 * @brief Test RNEA on two-link serial chain
 * 
 * Verifies velocity propagation and force propagation
 * Two revolute joints around Z-axis with offset
 */
TEST(InverseDynamicsTest, TwoLinkSerialChain) {
    InverseDynamics id;
    
    // Link 0 (base)
    InverseDynamicsLink link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    id.links.push_back(link0);
    
    // Link 1 (child of link 0, offset along X)
    InverseDynamicsLink link1;
    link1.parent = 0;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    id.links.push_back(link1);
    
    // Apply accelerations
    Eigen::VectorXd qddot(2);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    // Verify both joints produce torque
    EXPECT_GT(std::abs(tau[0]), 0.0);
    EXPECT_GT(std::abs(tau[1]), 0.0);
    // Torques should be finite
    EXPECT_TRUE(std::isfinite(tau[0]));
    EXPECT_TRUE(std::isfinite(tau[1]));
}

/**
 * @brief Test RNEA with branching kinematic tree (Y-configuration)
 * 
 * Verifies inward pass accumulates forces from multiple children
 * Three links: base with two children (symmetric branches)
 */
TEST(InverseDynamicsTest, BranchingKinematicTree) {
    InverseDynamics id;
    
    // Base link
    InverseDynamicsLink link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    id.links.push_back(link0);
    
    // Branch 1
    InverseDynamicsLink link1;
    link1.parent = 0;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    id.links.push_back(link1);
    
    // Branch 2 (identical to branch 1)
    InverseDynamicsLink link2;
    link2.parent = 0;
    link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    link2.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link2.q = 0.0;
    link2.qdot = 0.0;
    id.links.push_back(link2);
    
    // Apply equal accelerations to all joints
    Eigen::VectorXd qddot(3);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    qddot[2] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    // Verify: symmetric branches should have same torque
    EXPECT_NEAR(tau[1], tau[2], EPSILON);
    // Base torque should be non-zero (supports both branches)
    EXPECT_GT(std::abs(tau[0]), 0.0);
    EXPECT_TRUE(std::isfinite(tau[0]));
}

/**
 * @brief Test zero acceleration produces zero torque (static equilibrium)
 */
TEST(InverseDynamicsTest, ZeroAccelerationStaticEquilibrium) {
    InverseDynamics id;
    InverseDynamicsLink link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    id.links.push_back(link);
    
    Eigen::VectorXd qddot(1);
    qddot[0] = 0.0;
    
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    EXPECT_NEAR(tau[0], 0.0, EPSILON);
}

/**
 * @brief Test linearity: double acceleration produces double torque
 */
TEST(InverseDynamicsTest, LargeAccelerationProportionalTorque) {
    InverseDynamics id;
    InverseDynamicsLink link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    id.links.push_back(link);
    
    // Test with acceleration = 1.0
    Eigen::VectorXd qddot1(1);
    qddot1[0] = 1.0;
    Eigen::VectorXd tau1 = id.computeTorques(qddot1);
    
    // Test with acceleration = 2.0
    Eigen::VectorXd qddot2(1);
    qddot2[0] = 2.0;
    Eigen::VectorXd tau2 = id.computeTorques(qddot2);
    
    // Linearity: double acceleration = double torque
    EXPECT_NEAR(tau2[0], 2.0 * tau1[0], EPSILON * 2);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
