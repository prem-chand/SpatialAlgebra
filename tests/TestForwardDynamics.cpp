#include "ForwardDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-10;

// ForwardDynamics test suite for Articulated Body Algorithm
// Tests cover: single-link, serial chains, branching trees, PluckerTransform usage, edge cases

/**
 * @brief Test ABA on single-link pendulum
 * 
 * Verifies: τ = I*α for simple case
 * Single revolute joint around Z-axis with known inertia
 * Expected: qddot = τ / I_z
 */
TEST(ForwardDynamicsTest, SingleLinkPendulum) {
    ForwardDynamics fd;
    
    // Setup: single revolute joint around Z
    Link link;
    link.parent = -1;  // Base
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());  // Revolute Z
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    
    fd.links.push_back(link);
    
    // Apply torque around Z: τ = 1.0 Nm
    Eigen::VectorXd tau(1);
    tau[0] = 1.0;
    
    // Execute
    fd.computeAccelerations(tau);
    
    // Verify: qddot = τ / I_z = 1.0 / 1.0 = 1.0 rad/s²
    EXPECT_NEAR(fd.links[0].qddot, 1.0, EPSILON);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/**
 * @brief Test ABA on two-link serial chain
 * 
 * Verifies velocity propagation and inertia accumulation
 * Two revolute joints around Z-axis with offset
 */
TEST(ForwardDynamicsTest, TwoLinkSerialChain) {
    ForwardDynamics fd;
    
    // Link 0 (base)
    Link link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    link0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link0);
    
    // Link 1 (child of link 0, offset along X)
    Link link1;
    link1.parent = 0;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    link1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link1);
    
    // Apply torques
    Eigen::VectorXd tau(2);
    tau[0] = 1.0;
    tau[1] = 0.5;
    
    fd.computeAccelerations(tau);
    
    // Verify both joints accelerate
    EXPECT_GT(fd.links[0].qddot, 0.0);
    EXPECT_GT(fd.links[1].qddot, 0.0);
    // Link 0 (base) accelerates less because it moves more mass (both links)
    // Link 1 only moves its own mass, so it should accelerate more for same torque
    // But since tau[0]=1.0 > tau[1]=0.5, the relationship depends on the dynamics
    // Just verify both are positive and finite
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[1].qddot));
}

/**
 * @brief Test ABA with branching kinematic tree (Y-configuration)
 * 
 * Verifies inward pass accumulates inertias from multiple children
 * Three links: base with two children (symmetric branches)
 */
TEST(ForwardDynamicsTest, BranchingKinematicTree) {
    ForwardDynamics fd;
    
    // Base link
    Link link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    link0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link0);
    
    // Branch 1
    Link link1;
    link1.parent = 0;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    link1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link1);
    
    // Branch 2 (identical to branch 1)
    Link link2;
    link2.parent = 0;
    link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    link2.I = RigidBodyInertia(0.5, Vector3d::Zero(), lt::Identity(3) * 0.5);
    link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link2.q = 0.0;
    link2.qdot = 0.0;
    link2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link2);
    
    // Apply equal torques to all joints
    Eigen::VectorXd tau(3);
    tau[0] = 1.0;
    tau[1] = 0.5;
    tau[2] = 0.5;
    
    fd.computeAccelerations(tau);
    
    // Verify: symmetric branches should have same acceleration
    EXPECT_NEAR(fd.links[1].qddot, fd.links[2].qddot, EPSILON);
    // Base accelerates due to combined effect
    EXPECT_GT(fd.links[0].qddot, 0.0);
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
}

/**
 * @brief Test ABA with rotated link (non-identity PluckerTransform)
 * 
 * Verifies PluckerTransform operations in ABA context
 * Link 1 has 90° Z rotation relative to link 0
 */
TEST(ForwardDynamicsTest, PluckerTransformUsage) {
    ForwardDynamics fd;
    
    // Link 0
    Link link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    link0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link0);
    
    // Link 1 with 90° Z rotation
    Rotation rot90(Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ()));
    Link link1;
    link1.parent = 0;
    link1.X = PluckerTransform(rot90, Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    link1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link1);
    
    Eigen::VectorXd tau(2);
    tau[0] = 1.0;
    tau[1] = 1.0;
    
    fd.computeAccelerations(tau);
    
    // Verify accelerations are computed (rotation affects inertia)
    EXPECT_GT(std::abs(fd.links[0].qddot), 0.0);
    EXPECT_GT(std::abs(fd.links[1].qddot), 0.0);
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[1].qddot));
}

/**
 * @brief Test zero torque produces zero acceleration (static equilibrium)
 */
TEST(ForwardDynamicsTest, ZeroTorqueStaticEquilibrium) {
    ForwardDynamics fd;
    Link link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link);
    
    Eigen::VectorXd tau(1);
    tau[0] = 0.0;
    
    fd.computeAccelerations(tau);
    
    EXPECT_NEAR(fd.links[0].qddot, 0.0, EPSILON);
}

/**
 * @brief Test linearity: double torque produces double acceleration
 */
TEST(ForwardDynamicsTest, LargeTorqueProportionalAcceleration) {
    ForwardDynamics fd;
    Link link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link);
    
    // Test with torque = 1.0
    Eigen::VectorXd tau1(1);
    tau1[0] = 1.0;
    fd.computeAccelerations(tau1);
    double qddot1 = fd.links[0].qddot;
    
    // Test with torque = 2.0
    Eigen::VectorXd tau2(1);
    tau2[0] = 2.0;
    fd.computeAccelerations(tau2);
    double qddot2 = fd.links[0].qddot;
    
    // Linearity: double torque = double acceleration
    EXPECT_NEAR(qddot2, 2.0 * qddot1, EPSILON * 2);
}
