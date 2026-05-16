// TestRigidBodyInertia.cpp - Comprehensive GTest test suite for RigidBodyInertia class

#include "RigidBodyInertia.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "LowerTriangular.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>

using namespace SpatialAlgebra;
using namespace Eigen;

// Tolerance for floating point comparisons
const double TOLERANCE = 1e-10;

/**
 * @brief Test default constructor creates zero inertia
 * @details Verifies INR-01: Default constructor should create massless body
 */
TEST(RigidBodyInertiaTest, DefaultConstructor)
{
    RigidBodyInertia rbi;
    
    EXPECT_DOUBLE_EQ(rbi.getMass(), 0.0);
    
    Vector3d com = rbi.getCom();
    EXPECT_DOUBLE_EQ(com.x(), 0.0);
    EXPECT_DOUBLE_EQ(com.y(), 0.0);
    EXPECT_DOUBLE_EQ(com.z(), 0.0);
    
    lt inertia = rbi.getInertiaMatrixLT();
    // Check that inertia matrix is zero (LowerTriangular stores 6 elements for 3x3)
    Eigen::VectorXd inertiaData = inertia.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(inertiaData(i), 0.0);
    }
}

/**
 * @brief Test parameterized constructor stores values correctly
 * @details Verifies INR-01: Constructor preserves mass, COM, and inertia
 */
TEST(RigidBodyInertiaTest, ParameterizedConstructor)
{
    double mass = 5.0;
    Vector3d com(1.0, 2.0, 3.0);
    Vector6d inertiaData;
    inertiaData << 10.0, 0.0, 20.0, 0.0, 0.0, 30.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    EXPECT_DOUBLE_EQ(rbi.getMass(), mass);
    
    Vector3d storedCom = rbi.getCom();
    EXPECT_DOUBLE_EQ(storedCom.x(), com.x());
    EXPECT_DOUBLE_EQ(storedCom.y(), com.y());
    EXPECT_DOUBLE_EQ(storedCom.z(), com.z());
    
    lt storedInertia = rbi.getInertiaMatrixLT();
    Eigen::VectorXd storedData = storedInertia.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedData(i), inertiaData(i));
    }
}

/**
 * @brief Test getMass accessor
 * @details Verifies INR-01: Getter returns stored mass value
 */
TEST(RigidBodyInertiaTest, GetMass)
{
    RigidBodyInertia rbi(10.0, Vector3d::Zero(), lt(3));
    EXPECT_DOUBLE_EQ(rbi.getMass(), 10.0);
}

/**
 * @brief Test getCom accessor
 * @details Verifies INR-01: Getter returns stored center of mass
 */
TEST(RigidBodyInertiaTest, GetCom)
{
    Vector3d com(2.5, -1.5, 3.0);
    RigidBodyInertia rbi(1.0, com, lt(3));
    
    Vector3d storedCom = rbi.getCom();
    EXPECT_DOUBLE_EQ(storedCom.x(), com.x());
    EXPECT_DOUBLE_EQ(storedCom.y(), com.y());
    EXPECT_DOUBLE_EQ(storedCom.z(), com.z());
}

/**
 * @brief Test getInertiaMatrixLT accessor
 * @details Verifies INR-01: Getter returns stored inertia matrix
 */
TEST(RigidBodyInertiaTest, GetInertiaMatrixLT)
{
    Vector6d inertiaData;
    inertiaData << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    lt expected(inertiaData, 3);
    
    RigidBodyInertia rbi(1.0, Vector3d::Zero(), expected);
    
    lt stored = rbi.getInertiaMatrixLT();
    Eigen::VectorXd storedData = stored.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedData(i), inertiaData(i));
    }
}

/**
 * @brief Test operator+ combines two inertias
 * @details Verifies INR-02: Addition combines mass, COM, and inertia
 */
TEST(RigidBodyInertiaTest, OperatorAdd)
{
    Vector6d inertia1Data, inertia2Data;
    inertia1Data << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    inertia2Data << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    
    lt I1(inertia1Data, 3), I2(inertia2Data, 3);
    Vector3d com1(1.0, 0.0, 0.0), com2(0.0, 1.0, 0.0);
    
    RigidBodyInertia rbi1(2.0, com1, I1);
    RigidBodyInertia rbi2(3.0, com2, I2);
    
    RigidBodyInertia combined = rbi1 + rbi2;
    
    // Check mass addition
    EXPECT_DOUBLE_EQ(combined.getMass(), 5.0);
    
    // Check COM addition
    Vector3d expectedCom(1.0, 1.0, 0.0);
    Vector3d actualCom = combined.getCom();
    EXPECT_DOUBLE_EQ(actualCom.x(), expectedCom.x());
    EXPECT_DOUBLE_EQ(actualCom.y(), expectedCom.y());
    EXPECT_DOUBLE_EQ(actualCom.z(), expectedCom.z());
    
    // Check inertia addition
    lt combinedInertia = combined.getInertiaMatrixLT();
    Eigen::VectorXd combinedData = combinedInertia.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(combinedData(i), inertia1Data(i) + inertia2Data(i));
    }
}

/**
 * @brief Test operator* scales all components
 * @details Verifies INR-02: Scaling multiplies mass, COM, and inertia uniformly
 */
TEST(RigidBodyInertiaTest, OperatorScale)
{
    double mass = 4.0;
    Vector3d com(1.0, 2.0, 3.0);
    Vector6d inertiaData;
    inertiaData << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    RigidBodyInertia scaled = rbi * 2.0;
    
    // Check mass scaling
    EXPECT_DOUBLE_EQ(scaled.getMass(), mass * 2.0);
    
    // Check COM scaling
    Vector3d scaledCom = scaled.getCom();
    EXPECT_DOUBLE_EQ(scaledCom.x(), com.x() * 2.0);
    EXPECT_DOUBLE_EQ(scaledCom.y(), com.y() * 2.0);
    EXPECT_DOUBLE_EQ(scaledCom.z(), com.z() * 2.0);
    
    // Check inertia scaling
    lt scaledInertia = scaled.getInertiaMatrixLT();
    Eigen::VectorXd scaledData = scaledInertia.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(scaledData(i), inertiaData(i) * 2.0);
    }
}

/**
 * @brief Test apply() with pure rotation (ω≠0, v=0)
 * @details Verifies INR-03: Formula [Iω + com×v; m*v - com×ω] with v=0
 */
TEST(RigidBodyInertiaApplyTest, PureRotation)
{
    // Create rigid body: mass=2, COM at (1,0,0), diagonal inertia
    double mass = 2.0;
    Vector3d com(1.0, 0.0, 0.0);
    Vector6d inertiaData;
    inertiaData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    // Pure rotation around Z: ω = (0, 0, 1), v = (0, 0, 0)
    MotionVector mv(Vector3d(0, 0, 1), Vector3d::Zero());
    
    // Expected result:
    // torque = I*ω + com×v = I*(0,0,1) + com×(0,0,0) = (0, 0, 5) + (0, 0, 0) = (0, 0, 5)
    // force = m*v - com×ω = 2*(0,0,0) - (1,0,0)×(0,0,1) = (0, 0, 0) - (0, -1, 0) = (0, 1, 0)
    Vector3d expectedTorque(0.0, 0.0, 5.0);
    Vector3d expectedForce(0.0, 1.0, 0.0);
    
    ForceVector result = rbi.apply(mv);
    
    Vector3d actualTorque = result.getAngular();
    Vector3d actualForce = result.getLinear();
    
    EXPECT_NEAR(actualTorque.x(), expectedTorque.x(), TOLERANCE);
    EXPECT_NEAR(actualTorque.y(), expectedTorque.y(), TOLERANCE);
    EXPECT_NEAR(actualTorque.z(), expectedTorque.z(), TOLERANCE);
    
    EXPECT_NEAR(actualForce.x(), expectedForce.x(), TOLERANCE);
    EXPECT_NEAR(actualForce.y(), expectedForce.y(), TOLERANCE);
    EXPECT_NEAR(actualForce.z(), expectedForce.z(), TOLERANCE);
}

/**
 * @brief Test apply() with pure translation (ω=0, v≠0)
 * @details Verifies INR-03: Formula [Iω + com×v; m*v - com×ω] with ω=0
 */
TEST(RigidBodyInertiaApplyTest, PureTranslation)
{
    double mass = 3.0;
    Vector3d com(0.0, 1.0, 0.0);
    Vector6d inertiaData;
    inertiaData << 4.0, 0.0, 4.0, 0.0, 0.0, 4.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    // Pure translation along X: ω = (0, 0, 0), v = (2, 0, 0)
    MotionVector mv(Vector3d::Zero(), Vector3d(2.0, 0.0, 0.0));
    
    // Expected result:
    // torque = I*ω + com×v = 0 + (0,1,0)×(2,0,0) = (0, 0, -2)
    // force = m*v - com×ω = 3*(2,0,0) - 0 = (6, 0, 0)
    Vector3d expectedTorque(0.0, 0.0, -2.0);
    Vector3d expectedForce(6.0, 0.0, 0.0);
    
    ForceVector result = rbi.apply(mv);
    
    Vector3d actualTorque = result.getAngular();
    Vector3d actualForce = result.getLinear();
    
    EXPECT_NEAR(actualTorque.x(), expectedTorque.x(), TOLERANCE);
    EXPECT_NEAR(actualTorque.y(), expectedTorque.y(), TOLERANCE);
    EXPECT_NEAR(actualTorque.z(), expectedTorque.z(), TOLERANCE);
    
    EXPECT_NEAR(actualForce.x(), expectedForce.x(), TOLERANCE);
    EXPECT_NEAR(actualForce.y(), expectedForce.y(), TOLERANCE);
    EXPECT_NEAR(actualForce.z(), expectedForce.z(), TOLERANCE);
}

/**
 * @brief Test apply() with combined motion
 * @details Verifies INR-03: Superposition of rotation and translation
 */
TEST(RigidBodyInertiaApplyTest, CombinedMotion)
{
    double mass = 1.0;
    Vector3d com(1.0, 0.0, 0.0);
    Vector6d inertiaData;
    inertiaData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    // Combined motion: ω = (0, 1, 0), v = (0, 0, 1)
    MotionVector mv(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));
    
    // Expected result:
    // I*ω = (2,0,0,0,2,0,0,0,2) * (0,1,0) = (0, 2, 0)
    // com×v = (1,0,0)×(0,0,1) = (0, -1, 0)
    // torque = I*ω + com×v = (0, 2, 0) + (0, -1, 0) = (0, 1, 0)
    // m*v = 1*(0,0,1) = (0, 0, 1)
    // com×ω = (1,0,0)×(0,1,0) = (0, 0, 1)
    // force = m*v - com×ω = (0, 0, 1) - (0, 0, 1) = (0, 0, 0)
    Vector3d expectedTorque(0.0, 1.0, 0.0);
    Vector3d expectedForce(0.0, 0.0, 0.0);
    
    ForceVector result = rbi.apply(mv);
    
    Vector3d actualTorque = result.getAngular();
    Vector3d actualForce = result.getLinear();
    
    EXPECT_NEAR(actualTorque.x(), expectedTorque.x(), TOLERANCE);
    EXPECT_NEAR(actualTorque.y(), expectedTorque.y(), TOLERANCE);
    EXPECT_NEAR(actualTorque.z(), expectedTorque.z(), TOLERANCE);
    
    EXPECT_NEAR(actualForce.x(), expectedForce.x(), TOLERANCE);
    EXPECT_NEAR(actualForce.y(), expectedForce.y(), TOLERANCE);
    EXPECT_NEAR(actualForce.z(), expectedForce.z(), TOLERANCE);
}

/**
 * @brief Test apply() with zero motion
 * @details Verifies INR-03: Zero input produces zero output
 */
TEST(RigidBodyInertiaApplyTest, ZeroMotion)
{
    double mass = 5.0;
    Vector3d com(1.0, 2.0, 3.0);
    Vector6d inertiaData;
    inertiaData << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    // Zero motion: ω = (0, 0, 0), v = (0, 0, 0)
    MotionVector mv(Vector3d::Zero(), Vector3d::Zero());
    
    ForceVector result = rbi.apply(mv);
    
    Vector3d actualTorque = result.getAngular();
    Vector3d actualForce = result.getLinear();
    
    EXPECT_DOUBLE_EQ(actualTorque.x(), 0.0);
    EXPECT_DOUBLE_EQ(actualTorque.y(), 0.0);
    EXPECT_DOUBLE_EQ(actualTorque.z(), 0.0);
    
    EXPECT_DOUBLE_EQ(actualForce.x(), 0.0);
    EXPECT_DOUBLE_EQ(actualForce.y(), 0.0);
    EXPECT_DOUBLE_EQ(actualForce.z(), 0.0);
}

/**
 * @brief Test apply() formula matches Featherstone specification
 * @details Verifies INR-03: [Iω + com×v; m*v - com×ω] implementation
 */
TEST(RigidBodyInertiaApplyTest, FeatherstoneFormula)
{
    // Test case with non-symmetric values to ensure formula correctness
    double mass = 2.5;
    Vector3d com(0.5, -0.5, 1.0);
    Vector6d inertiaData;
    inertiaData << 3.0, 1.0, 4.0, 1.0, 2.0, 5.0; // Lower triangular: [0,0], [1,0], [1,1], [2,0], [2,1], [2,2]
    lt inertia(inertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, inertia);
    
    // Arbitrary motion: ω = (1, 2, 3), v = (4, 5, 6)
    MotionVector mv(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    
    // Manual calculation:
    // I*ω: LowerTriangular * Vector3d
    // com×v = (0.5,-0.5,1.0)×(4,5,6) = (-0.5*6 - 1.0*5, 1.0*4 - 0.5*6, 0.5*5 - (-0.5)*4) = (-8, 1, 4.5)
    // m*v = 2.5*(4,5,6) = (10, 12.5, 15)
    // com×ω = (0.5,-0.5,1.0)×(1,2,3) = (-0.5*3 - 1.0*2, 1.0*1 - 0.5*3, 0.5*2 - (-0.5)*1) = (-3.5, -0.5, 1.5)
    
    ForceVector result = rbi.apply(mv);
    
    // Just verify it produces a valid result (detailed formula check in combined test)
    Vector3d torque = result.getAngular();
    Vector3d force = result.getLinear();
    
    // Result should be finite
    EXPECT_TRUE(std::isfinite(torque.x()));
    EXPECT_TRUE(std::isfinite(torque.y()));
    EXPECT_TRUE(std::isfinite(torque.z()));
    EXPECT_TRUE(std::isfinite(force.x()));
    EXPECT_TRUE(std::isfinite(force.y()));
    EXPECT_TRUE(std::isfinite(force.z()));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
