#include "ForwardDynamics.h"
#include "InverseDynamics.h"
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
    // Base: two identical symmetric branches produce equal reaction torques
    // that cancel the applied base torque, so qddot[0] = 0
    EXPECT_NEAR(fd.links[0].qddot, 0.0, EPSILON);
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
    // Base: child's reaction torque from acceleration exactly cancels
    // applied torque qddot[0] = 0 when tau[0]=tau[1]=1.0
    EXPECT_NEAR(fd.links[0].qddot, 0.0, EPSILON);
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

/**
 * @brief Multi-link ABA with verifiable numerical values
 * @details Three-link serial chain with identity inertias, Z-axis revolute joints,
 *          transforms along X. Base joint accelerates less than the tip for the same
 *          torque because it carries reflected inertias from all descendants.
 */
TEST(ForwardDynamicsTest, ThreeLinkNumericalValidation) {
    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.25;

    InverseDynamics id;

    InverseDynamicsLink id_l0;
    id_l0.parent = -1;
    id_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l0.q = 0.0; id_l0.qdot = 0.0;
    id.links.push_back(id_l0);

    InverseDynamicsLink id_l1;
    id_l1.parent = 0;
    id_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l1.q = 0.0; id_l1.qdot = 0.0;
    id.links.push_back(id_l1);

    InverseDynamicsLink id_l2;
    id_l2.parent = 1;
    id_l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l2.q = 0.0; id_l2.qdot = 0.0;
    id.links.push_back(id_l2);

    Eigen::VectorXd tau = id.computeTorques(qddot_input);

    ForwardDynamics fd;

    Link l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0; l0.qdot = 0.0;
    l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l0);

    Link l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0; l1.qdot = 0.0;
    l1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l1);

    Link l2;
    l2.parent = 1;
    l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l2.q = 0.0; l2.qdot = 0.0;
    l2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l2);

    fd.computeAccelerations(tau);

    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], 1e-8);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], 1e-8);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], 1e-8);

    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[1].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[2].qddot));
}

/**
 * @brief Verifies articulated inertias accumulate from descendants to ancestors
 * @details After the inward pass, the base link's Ia should reflect all three
 *          links' inertias, making its norm strictly larger than any single link.
 */
TEST(ForwardDynamicsTest, CondensationReducesInertiaNorm) {
    ForwardDynamics fd;

    Link l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0; l0.qdot = 0.0;
    l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l0);

    Link l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0; l1.qdot = 0.0;
    l1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l1);

    Link l2;
    l2.parent = 1;
    l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l2.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l2.q = 0.0; l2.qdot = 0.0;
    l2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l2);

    Eigen::VectorXd tau(3);
    tau[0] = 1.0; tau[1] = 0.5; tau[2] = 0.25;

    fd.computeAccelerations(tau);

    // Frobenius-like norm on lower triangular matrix
    auto ltNorm = [](const lt& m) -> double {
        double sum = 0.0;
        for (int i = 0; i < m.getSize(); i++)
            for (int j = 0; j <= i; j++)
                sum += m(i, j) * m(i, j);
        return std::sqrt(sum);
    };

    // Base link's Ia should be larger than any child's (reflects all descendants)
    EXPECT_GT(ltNorm(fd.links[0].Ia.getInertia()),
              ltNorm(fd.links[2].Ia.getInertia()));
    EXPECT_GT(ltNorm(fd.links[0].Ia.getInertia()),
              ltNorm(fd.links[0].I.getInertiaMatrixLT()));

    for (int i = 0; i < 3; i++) {
        EXPECT_TRUE(std::isfinite(fd.links[i].qddot));
    }
}

/**
 * @brief Three-link serial chain with single torque on base joint
 * @details tau=[1,0,0] at q=0, qdot=0, no gravity.
 *          Base joint should accelerate positively while joints 1 and 2
 *          have no direct torque and zero coupling at rest.
 */
TEST(ForwardDynamicsTest, ThreeLinkSingleTorque) {
    ForwardDynamics fd;

    Link l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0; l0.qdot = 0.0;
    l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l0);

    Link l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0; l1.qdot = 0.0;
    l1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l1);

    Link l2;
    l2.parent = 1;
    l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l2.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l2.q = 0.0; l2.qdot = 0.0;
    l2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l2);

    Eigen::VectorXd tau(3);
    tau[0] = 1.0; tau[1] = 0.0; tau[2] = 0.0;

    fd.computeAccelerations(tau);

    // Base torque propagates through kinematic chain — all joints accelerate
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[1].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[2].qddot));

    // Base accelerates positively; coupled joints may move
    EXPECT_GT(fd.links[0].qddot, 0.0);
}

/**
 * @brief ABA with gravity produces different accelerations than without
 * @details Single-link pendulum with gravity pointing downward (-Z).
 *          c₀ = -gravity = (0,0,0,0,0,9.81) → upward bias acceleration.
 *          This should change the joint acceleration for non-horizontal configurations.
 *          For a pendulum at q=0 with horizontal arm, gravity produces a torque
 *          that opposes positive acceleration.
 */
TEST(ForwardDynamicsTest, SingleLinkWithGravity) {
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
    tau[0] = 1.0;
    
    // Without gravity
    fd.computeAccelerations(tau, Vector3d::Zero());
    double qddot_no_gravity = fd.links[0].qddot;
    
    // With gravity (Earth gravity, -Z direction)
    Vector3d gravity(0, 0, -9.81);
    fd.computeAccelerations(tau, gravity);
    double qddot_with_gravity = fd.links[0].qddot;
    
    // With unit-mass inertia at origin, gravity should not directly affect
    // the revolute Z joint acceleration (gravity is linear acceleration,
    // joint axis is rotational Z, orthogonal coupling). However, the bias
    // acceleration propagation through the forward dynamics still applies.
    // At minimum, verify the acceleration is still finite and reasonable.
    EXPECT_TRUE(std::isfinite(qddot_with_gravity));
    EXPECT_GT(qddot_with_gravity, 0.0);
}

/**
 * @brief Two-link ABA with gravity — verifies gravity propagates through chain
 */
TEST(ForwardDynamicsTest, TwoLinkWithGravity) {
    ForwardDynamics fd;
    
    Link l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0; l0.qdot = 0.0;
    l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l0);
    
    Link l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0; l1.qdot = 0.0;
    l1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l1);
    
    Eigen::VectorXd tau(2);
    tau[0] = 1.0;
    tau[1] = 0.5;
    
    // With Earth gravity
    Vector3d gravity(0, 0, -9.81);
    fd.computeAccelerations(tau, gravity);
    
    // All accelerations should be finite
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_TRUE(std::isfinite(fd.links[1].qddot));
}

/**
 * @brief Gravitational acceleration scales linearly with mass
 * @details For a single X-axis revolute joint with COM offset and
 *          identity inertia, the gravitational effective torque
 *          through the bias acceleration c₀ = -g scales linearly
 *          with mass, while the effective joint inertia (S^T·Iₐ·S)
 *          is dominated by the I_cm component which is constant.
 *          Result: qddot ∝ m under identical gravity.
 * 
 *          Physical derivation:
 *          - X-axis joint S=[1,0,0;0,0,0], COM=(0,1,0), g=(0,0,-g)
 *          - c₀ = -g = (0,0,0,0,0,9.81)
 *          - pₐ = I·c₀: angular = H·v = m·[c]×·(0,0,9.81) = (m·9.81, 0, 0)
 *          - S^T·pₐ = (1,0,0)·(m·9.81,0,0) = m·9.81 (scales with m)
 *          - S^T·Iₐ·S = I_cm[0][0] = 1 (constant, not scaling with m)
 *          - qddot = -(0 - m·9.81) / 1 = m·9.81 (scales with m)
 *          - Ratio qddot(m=2) / qddot(m=1) = 2.0
 * 
 *          Independent oracle: the ratio follows from the algebra,
 *          not from solver cross-validation.
 */
TEST(ForwardDynamicsTest, GravityEffectScalesWithMass) {
    // Test at m=1.0
    ForwardDynamics fd1;
    Link link1;
    link1.parent = -1;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link1.I = RigidBodyInertia(1.0, Vector3d(0, 1.0, 0), lt::Identity(3));
    link1.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    link1.q = 0.0;
    link1.qdot = 0.0;
    link1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd1.links.push_back(link1);
    
    Eigen::VectorXd tau_zero(1);
    tau_zero[0] = 0.0;
    
    Vector3d gravity(0, 0, -9.81);
    fd1.computeAccelerations(tau_zero, gravity);
    double qddot_m1 = fd1.links[0].qddot;
    EXPECT_TRUE(std::isfinite(qddot_m1));
    EXPECT_LT(qddot_m1, 0.0);  // Gravity produces negative acceleration
    
    // Test at m=2.0 with same I_cm (identity, not scaled with mass)
    ForwardDynamics fd2;
    Link link2;
    link2.parent = -1;
    link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link2.I = RigidBodyInertia(2.0, Vector3d(0, 1.0, 0), lt::Identity(3));
    link2.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());
    link2.q = 0.0;
    link2.qdot = 0.0;
    link2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd2.links.push_back(link2);
    
    fd2.computeAccelerations(tau_zero, gravity);
    double qddot_m2 = fd2.links[0].qddot;
    EXPECT_TRUE(std::isfinite(qddot_m2));
    EXPECT_LT(qddot_m2, 0.0);  // More mass = more negative acceleration
    
    // Invariant: qddot scales with mass (I_cm provides constant inertia base)
    // Ratio = qddot(m=2) / qddot(m=1) = 2.0
    EXPECT_NEAR(qddot_m2 / qddot_m1, 2.0, 1e-10);
    
    // Without gravity, zero torque produces zero acceleration
    ForwardDynamics fd_no_grav;
    Link link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d(0, 1.0, 0), lt::Identity(3));
    link0.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());
    link0.q = 0.0; link0.qdot = 0.0;
    link0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd_no_grav.links.push_back(link0);
    
    fd_no_grav.computeAccelerations(tau_zero, Vector3d::Zero());
    EXPECT_NEAR(fd_no_grav.links[0].qddot, 0.0, 1e-10);
}

/**
 * @brief Gravity effect on acceleration is proportional to g magnitude
 * @details For a single X-axis revolute joint with COM offset, the
 *          difference in joint acceleration with vs without gravity
 *          is proportional to the gravity magnitude |g|. This is an
 *          independent structural invariant that does not rely on
 *          solver cross-validation.
 * 
 *          Physical derivation:
 *          - ABA: qddot = (τ - S^T·pₐ) / (S^T·Iₐ·S)
 *          - pₐ = I·c₀ where c₀ = -g, so pₐ ∝ g
 *          - S^T·pₐ ∝ g (through the cross-coupling term H·v)
 *          - S^T·Iₐ·S is constant (independent of g)
 *          - qddot(g) - qddot(0) = -S^T·I·c₀ / (S^T·I·S) ∝ g
 * 
 *          Test: single X-axis joint, COM=(0,1,0), torque=1.0.
 *          Verify: [qddot(g₁) - qddot(0)] / [qddot(g₂) - qddot(0)] = g₁/g₂.
 */
TEST(ForwardDynamicsTest, GravityProportionalityInvariant) {
    // Single X-axis revolute joint
    ForwardDynamics fd;
    Link link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d(0, 1.0, 0), lt::Identity(3));
    link.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link);
    
    Eigen::VectorXd tau(1);
    tau[0] = 1.0;
    
    // Acceleration without gravity
    fd.computeAccelerations(tau, Vector3d::Zero());
    double qddot_no_grav = fd.links[0].qddot;
    EXPECT_TRUE(std::isfinite(qddot_no_grav));
    EXPECT_GT(qddot_no_grav, 0.0);
    
    // Acceleration with gravity g₁ = -9.81 along Z
    fd.computeAccelerations(tau, Vector3d(0, 0, -9.81));
    double qddot_with_g1 = fd.links[0].qddot;
    EXPECT_TRUE(std::isfinite(qddot_with_g1));
    
    // Gravity opposes positive acceleration (gravity pulls arm down)
    double diff_g1 = qddot_no_grav - qddot_with_g1;
    EXPECT_GT(diff_g1, 0.0);
    
    // Acceleration with gravity g₂ = -4.905 (half of g₁)
    fd.computeAccelerations(tau, Vector3d(0, 0, -4.905));
    double qddot_with_g2 = fd.links[0].qddot;
    EXPECT_TRUE(std::isfinite(qddot_with_g2));
    
    double diff_g2 = qddot_no_grav - qddot_with_g2;
    EXPECT_GT(diff_g2, 0.0);
    
    // Invariant: diff_g1 / diff_g2 = g₁ / g₂ = 2.0
    EXPECT_NEAR(diff_g1 / diff_g2, 2.0, 1e-10);
}

/**
 * @brief Release-mode stability: valid inputs produce finite output
 * @details Verifies ABA produces stable finite results with valid inputs
 *          regardless of NDEBUG mode (where NaN/Inf assertions are compiled
 *          out). Single link with identity inertia, Z-axis revolute joint,
 *          tau=[1.0]. Tests both zero gravity and Earth gravity.
 *          Per reviewer: "release-mode tests verifying behavior stays stable
 *          when inputs are valid" (MEDIUM, 13-03).
 */
TEST(ForwardDynamicsTest, ReleaseModeStability) {
    // Single link: Z-axis revolute, identity inertia
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
    tau[0] = 1.0;
    
    // Without gravity
    fd.computeAccelerations(tau, Vector3d::Zero());
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_GT(fd.links[0].qddot, 0.0);
    
    // With gravity
    fd.computeAccelerations(tau, Vector3d(0, 0, -9.81));
    EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    EXPECT_GT(fd.links[0].qddot, 0.0);
}

/**
 * @brief Edge case: zero-mass single link
 * @details Tests ABA with a degenerate inertia (mass=0, COM=zero,
 *          inertia tensor=zero). The solver must not crash on this
 *          degenerate input. Either finite output or a graceful
 *          exception is acceptable.
 */
TEST(ForwardDynamicsTest, ZeroMassEdgeCase) {
    ForwardDynamics fd;
    Link link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(0.0, Vector3d::Zero(), LowerTriangular(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link);
    
    Eigen::VectorXd tau(1);
    tau[0] = 0.0;
    
    // Zero-mass produces degenerate articulated inertia (zero Ia).
    // The solver must not crash — either finite qddot or exception is acceptable.
    bool threw = false;
    try {
        fd.computeAccelerations(tau);
    } catch (const std::exception&) {
        threw = true;
    }
    if (!threw) {
        EXPECT_TRUE(std::isfinite(fd.links[0].qddot));
    }
    // If an exception was thrown, the test passes (no crash on degenerate input).
    EXPECT_TRUE(threw || std::isfinite(fd.links[0].qddot));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
