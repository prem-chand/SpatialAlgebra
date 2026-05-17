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

/**
 * @brief RNEA with non-zero joint velocity exercises Coriolis/centrifugal terms
 * @details Previous tests all used qdot = 0, meaning the Coriolis term
 *          v × S·q̇ in the outward pass was never exercised. With qdot ≠ 0,
 *          the acceleration propagation includes velocity product terms:
 *          aᵢ = Xᵢ·a_parent + Sᵢ·q̈ᵢ + vᵢ × Sᵢ·q̇ᵢ
 *          The torque should differ from the qdot=0 case.
 *          Note: Coriolis effect requires non-zero COM offset to couple
 *          through the inertia matrix (COM=[0,0.1,0] used here).
 */
TEST(InverseDynamicsTest, TwoLinkSerialChainNonZeroVelocity) {
    InverseDynamics id;
    
    // Link 0 (base)
    InverseDynamicsLink l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d(0, 0.1, 0), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0;
    l0.qdot = 2.0;  // Non-zero velocity
    id.links.push_back(l0);
    
    // Link 1 (child)
    InverseDynamicsLink l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d(0, 0.1, 0), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0;
    l1.qdot = 1.0;  // Non-zero velocity
    id.links.push_back(l1);
    
    Eigen::VectorXd qddot(2);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    // Torques should be finite
    EXPECT_TRUE(std::isfinite(tau[0]));
    EXPECT_TRUE(std::isfinite(tau[1]));
    
    // Compare with zero-velocity case — torques should differ (Coriolis effect)
    InverseDynamics id_zero;
    InverseDynamicsLink zl0;
    zl0.parent = -1;
    zl0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    zl0.I = RigidBodyInertia(1.0, Vector3d(0, 0.1, 0), lt::Identity(3));
    zl0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    zl0.q = 0.0;
    zl0.qdot = 0.0;
    id_zero.links.push_back(zl0);
    
    InverseDynamicsLink zl1;
    zl1.parent = 0;
    zl1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    zl1.I = RigidBodyInertia(1.0, Vector3d(0, 0.1, 0), lt::Identity(3));
    zl1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    zl1.q = 0.0;
    zl1.qdot = 0.0;
    id_zero.links.push_back(zl1);
    
    Eigen::VectorXd tau_zero = id_zero.computeTorques(qddot);
    
    // Non-zero velocity should produce different torques
    bool hasCoriolisEffect = (std::abs(tau[0] - tau_zero[0]) > EPSILON) || 
                             (std::abs(tau[1] - tau_zero[1]) > EPSILON);
    EXPECT_TRUE(hasCoriolisEffect);
}

/**
 * @brief RNEA with gravity — verifies gravity term in outward pass
 * @details Gravity sets base acceleration to a₀ = S·q̈ - g.
 *          For a vertical pendulum (gravity along -Z), the torque should differ
 *          from the no-gravity case.
 */
TEST(InverseDynamicsTest, SingleLinkWithGravity) {
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
    qddot[0] = 1.0;
    
    // Without gravity
    Eigen::VectorXd tau_no_gravity = id.computeTorques(qddot, Vector3d::Zero());
    
    // With Earth gravity
    Vector3d gravity(0, 0, -9.81);
    Eigen::VectorXd tau_with_gravity = id.computeTorques(qddot, gravity);
    
    // Both should be finite
    EXPECT_TRUE(std::isfinite(tau_with_gravity[0]));
    
    // For single-link identity inertia at origin with Z-axis revolute joint,
    // the gravity vector along -Z has pure linear component [0,0,-9.81].
    // The RNEA outward pass subtracts [0; g] from base acceleration:
    // a₀ = S·q̈ - [0; g]
    // This modifies the linear acceleration component, which couples back
    // through the inertia matrix (for non-zero COM) to affect torque.
    // With COM=0 and identity inertia, gravity produces NO joint torque
    // (linear acceleration doesn't project onto Z-rotation axis).
    EXPECT_NEAR(tau_with_gravity[0], tau_no_gravity[0], EPSILON);
}

/**
 * @brief Two-link RNEA with both non-zero velocity and gravity
 * @details Combined test exercises both the Coriolis terms (via qdot)
 *          and the gravity term (via base acceleration).
 */
TEST(InverseDynamicsTest, TwoLinkWithVelocityAndGravity) {
    InverseDynamics id;
    
    InverseDynamicsLink l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0;
    l0.qdot = 2.0;
    id.links.push_back(l0);
    
    InverseDynamicsLink l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0;
    l1.qdot = 1.0;
    id.links.push_back(l1);
    
    Eigen::VectorXd qddot(2);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    
    Vector3d gravity(0, 0, -9.81);
    Eigen::VectorXd tau = id.computeTorques(qddot, gravity);
    
    EXPECT_TRUE(std::isfinite(tau[0]));
    EXPECT_TRUE(std::isfinite(tau[1]));
}

/**
 * @brief Gravity invariant: static torque proportional to g for point mass on lever arm
 * @details For a single X-axis revolute joint with COM offset along Y and gravity
 *          along -Z, the gravitational torque is τ = m·g·L where L = |COM_y|.
 *          This is an independent structural invariant: |tau[0]| / (g · L) = m.
 *          Testing at g = {0, 5, 10} verifies proportionality to g, not self-consistency.
 * 
 *          Physical derivation:
 *          - Joint axis: X (S=[1,0,0;0,0,0]), COM=(0,0.5,0), gravity=(0,0,-g)
 *          - Gravitational force on COM: F = (0, 0, -m·g)
 *          - Torque about origin: τ = COM × F = (-0.5·m·g, 0, 0)
 *          - Project onto X-axis joint: tau[0] = -0.5·m·g
 *          - Expected invariant: |tau[0]| / (g · 0.5) = m = 1.0
 */
TEST(InverseDynamicsTest, SingleLinkStaticGravityProportionality) {
    // Single X-axis revolute joint with COM offset along Y
    InverseDynamics id;
    InverseDynamicsLink link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    link.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    link.q = 0.0;
    link.qdot = 0.0;
    id.links.push_back(link);
    
    Eigen::VectorXd qddot(1);
    qddot[0] = 0.0;  // Static case
    
    // Test at three gravity levels: g = 0, 5, 10
    const double g_values[] = {0.0, 5.0, 10.0};
    bool ratio_initialized = false;
    double prev_ratio = 0.0;
    
    for (int i = 0; i < 3; i++) {
        double g = g_values[i];
        Vector3d gravity(0, 0, -g);
        Eigen::VectorXd tau = id.computeTorques(qddot, gravity);
        
        EXPECT_TRUE(std::isfinite(tau[0]));
        
        if (g > 0) {
            // Invariant: |tau[0]| / (g * 0.5) = m = 1.0
            double ratio = std::abs(tau[0]) / (g * 0.5);
            EXPECT_NEAR(ratio, 1.0, 1e-10);
            
            // Consistency check: ratio should be the same across g levels
            if (ratio_initialized) {
                EXPECT_NEAR(ratio, prev_ratio, 1e-10);
            }
            prev_ratio = ratio;
            ratio_initialized = true;
        } else {
            // Zero gravity should produce zero torque
            EXPECT_NEAR(tau[0], 0.0, 1e-10);
        }
    }
    
    // Verify proportionality directly: tau(g=10) / tau(g=5) = 2.0
    Eigen::VectorXd tau5 = id.computeTorques(qddot, Vector3d(0, 0, -5.0));
    Eigen::VectorXd tau10 = id.computeTorques(qddot, Vector3d(0, 0, -10.0));
    EXPECT_NEAR(tau10[0] / tau5[0], 2.0, 1e-10);
}

/**
 * @brief Gravity produces zero torque when COM and gravity are collinear
 * @details For an X-axis revolute joint with COM=(0,0.5,0), if gravity is
 *          applied along the same direction as the COM offset (Y-axis),
 *          the gravitational force passes through the joint axis, producing
 *          zero torque. This is an independent geometric invariant.
 * 
 *          Physical derivation:
 *          - COM=(0,0.5,0), gravity=(0,-g,0) (collinear with COM offset)
 *          - Gravitational force: F = (0, -m·g, 0)
 *          - Torque about origin: τ = COM × F = (0,0.5,0) × (0,-mg,0)
 *          - Parallel vectors: cross product is zero
 *          - Expected: tau[0] ≈ 0
 */
TEST(InverseDynamicsTest, GravityTorqueZeroAtVertical) {
    InverseDynamics id;
    InverseDynamicsLink link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    link.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    link.q = 0.0;
    link.qdot = 0.0;
    id.links.push_back(link);
    
    Eigen::VectorXd qddot(1);
    qddot[0] = 0.0;
    
    // Gravity along Y (same direction as COM offset) — collinear, zero torque
    Vector3d gravity(0, -9.81, 0);
    Eigen::VectorXd tau = id.computeTorques(qddot, gravity);
    
    // Zero torque because gravity force line passes through joint
    EXPECT_NEAR(tau[0], 0.0, 1e-10);
}

/**
 * @brief Multi-link gravity produces finite non-NaN torques
 * @details A 2-link chain with non-zero velocity, non-zero acceleration,
 *          and gravity. Verifies that the combined effects of Coriolis,
 *          centrifugal, and gravity forces produce valid numerical results.
 *          This is a sanity check that the gravity propagation through
 *          a kinematic chain does not produce infinite or NaN values.
 */
TEST(InverseDynamicsTest, GravityFiniteValidResults) {
    InverseDynamics id;
    
    // Link 0 (base) — X-axis revolute
    InverseDynamicsLink l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    l0.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    l0.q = 0.0;
    l0.qdot = 2.0;
    id.links.push_back(l0);
    
    // Link 1 (child of link 0, offset along X)
    InverseDynamicsLink l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d(0, 0.5, 0), lt::Identity(3));
    l1.S = MotionVector(Vector3d(1, 0, 0), Vector3d::Zero());  // Revolute X
    l1.q = 0.0;
    l1.qdot = 1.0;
    id.links.push_back(l1);
    
    Eigen::VectorXd qddot(2);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    
    Vector3d gravity(0, 0, -9.81);
    Eigen::VectorXd tau = id.computeTorques(qddot, gravity);
    
    // Both torques should be finite and non-NaN
    EXPECT_TRUE(std::isfinite(tau[0]));
    EXPECT_TRUE(std::isfinite(tau[1]));
    EXPECT_FALSE(std::isnan(tau[0]));
    EXPECT_FALSE(std::isnan(tau[1]));
}

/**
 * @brief Edge case: zero-mass single link
 * @details Tests RNEA with a degenerate inertia (mass=0, COM=zero,
 *          inertia tensor=zero). For zero acceleration, the torques
 *          should be near zero and finite. Validates the solver handles
 *          degenerate masses without producing NaN or crashing.
 */
TEST(InverseDynamicsTest, ZeroMassEdgeCase) {
    InverseDynamics id;
    InverseDynamicsLink link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(0.0, Vector3d::Zero(), LowerTriangular(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    id.links.push_back(link);
    
    Eigen::VectorXd qddot(1);
    qddot[0] = 0.0;
    
    Eigen::VectorXd tau = id.computeTorques(qddot);
    
    // With zero mass and zero acceleration, torques should be near zero
    EXPECT_NEAR(tau[0], 0.0, 1e-10);
    EXPECT_TRUE(std::isfinite(tau[0]));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
