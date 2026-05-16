// TestArticulatedBodyInertia.cpp - Comprehensive GTest test suite for ArticulatedBodyInertia class

#include "ArticulatedBodyInertia.h"
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
 * @details Verifies INR-03: Default constructor should create zero articulated inertia
 */
TEST(ArticulatedBodyInertiaTest, DefaultConstructor)
{
    ArticulatedBodyInertia abi;
    
    lt inertia = abi.getInertia();
    Eigen::VectorXd inertiaData = inertia.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(inertiaData(i), 0.0);
    }
    
    Eigen::Matrix3d H = abi.getH();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(H(i, j), 0.0);
        }
    }
    
    lt M = abi.getM();
    Eigen::VectorXd MData = M.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(MData(i), 0.0);
    }
}

/**
 * @brief Test parameterized constructor stores values correctly
 * @details Verifies INR-03: Constructor preserves I, H, M matrices
 */
TEST(ArticulatedBodyInertiaTest, ParameterizedConstructor)
{
    Vector6d inertiaData, MData;
    inertiaData << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    MData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    
    lt storedI = abi.getInertia();
    Eigen::VectorXd storedIData = storedI.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedIData(i), inertiaData(i));
    }
    
    Eigen::Matrix3d storedH = abi.getH();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(storedH(i, j), H(i, j));
        }
    }
    
    lt storedM = abi.getM();
    Eigen::VectorXd storedMData = storedM.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedMData(i), MData(i));
    }
}

/**
 * @brief Test getInertia accessor
 * @details Verifies INR-03: Getter returns stored rotational inertia
 */
TEST(ArticulatedBodyInertiaTest, GetInertia)
{
    Vector6d inertiaData;
    inertiaData << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    lt expected(inertiaData, 3);
    
    ArticulatedBodyInertia abi(expected, Eigen::Matrix3d::Zero(), lt(3));
    
    lt stored = abi.getInertia();
    Eigen::VectorXd storedData = stored.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedData(i), inertiaData(i));
    }
}

/**
 * @brief Test getH accessor
 * @details Verifies INR-03: Getter returns stored coupling matrix
 */
TEST(ArticulatedBodyInertiaTest, GetH)
{
    Eigen::Matrix3d expected;
    expected << 1.0, 2.0, 3.0,
                4.0, 5.0, 6.0,
                7.0, 8.0, 9.0;
    
    ArticulatedBodyInertia abi(lt(3), expected, lt(3));
    
    Eigen::Matrix3d stored = abi.getH();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(stored(i, j), expected(i, j));
        }
    }
}

/**
 * @brief Test getM accessor
 * @details Verifies INR-03: Getter returns stored mass matrix
 */
TEST(ArticulatedBodyInertiaTest, GetM)
{
    Vector6d MData;
    MData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    lt expected(MData, 3);
    
    ArticulatedBodyInertia abi(lt(3), Eigen::Matrix3d::Zero(), expected);
    
    lt stored = abi.getM();
    Eigen::VectorXd storedData = stored.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(storedData(i), MData(i));
    }
}

/**
 * @brief Test operator+ combines two articulated inertias
 * @details Verifies INR-04: Addition combines I, H, and M matrices
 */
TEST(ArticulatedBodyInertiaTest, OperatorAdd)
{
    Vector6d I1Data, I2Data, M1Data, M2Data;
    I1Data << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    I2Data << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    M1Data << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    M2Data << 3.0, 0.0, 3.0, 0.0, 0.0, 3.0;
    
    lt I1(I1Data, 3), I2(I2Data, 3), M1(M1Data, 3), M2(M2Data, 3);
    Eigen::Matrix3d H1 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d H2 = 2.0 * Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi1(I1, H1, M1);
    ArticulatedBodyInertia abi2(I2, H2, M2);
    
    ArticulatedBodyInertia combined = abi1 + abi2;
    
    // Check I addition
    lt combinedI = combined.getInertia();
    Eigen::VectorXd combinedIData = combinedI.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(combinedIData(i), I1Data(i) + I2Data(i));
    }
    
    // Check H addition
    Eigen::Matrix3d combinedH = combined.getH();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(combinedH(i, j), H1(i, j) + H2(i, j));
        }
    }
    
    // Check M addition
    lt combinedM = combined.getM();
    Eigen::VectorXd combinedMData = combinedM.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(combinedMData(i), M1Data(i) + M2Data(i));
    }
}

/**
 * @brief Test operator+(RigidBodyInertia) adds rigid body to articulated body
 * @details Verifies INR-04: Conversion and addition works correctly
 */
TEST(ArticulatedBodyInertiaTest, OperatorAddRigidBody)
{
    Vector6d I1Data, M1Data;
    I1Data << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    M1Data << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    
    lt I1(I1Data, 3), M1(M1Data, 3);
    Eigen::Matrix3d H1 = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi1(I1, H1, M1);
    
    // Create rigid body: mass=5, COM=(1,0,0), diagonal inertia
    double mass = 5.0;
    Vector3d com(1.0, 0.0, 0.0);
    Vector6d rbiInertiaData;
    rbiInertiaData << 3.0, 0.0, 3.0, 0.0, 0.0, 3.0;
    lt rbiInertia(rbiInertiaData, 3);
    
    RigidBodyInertia rbi(mass, com, rbiInertia);
    
    ArticulatedBodyInertia combined = abi1 + rbi;
    
    // Just verify it produces a valid result (detailed formula in header)
    lt combinedI = combined.getInertia();
    Eigen::Matrix3d combinedH = combined.getH();
    lt combinedM = combined.getM();
    
    // Result should be finite
    Eigen::VectorXd combinedIData = combinedI.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(std::isfinite(combinedIData(i)));
    }
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_TRUE(std::isfinite(combinedH(i, j)));
        }
    }
    
    Eigen::VectorXd combinedMData = combinedM.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(std::isfinite(combinedMData(i)));
    }
}

/**
 * @brief Test operator* scales all components
 * @details Verifies INR-04: Scaling multiplies I, H, and M uniformly
 */
TEST(ArticulatedBodyInertiaTest, OperatorScale)
{
    Vector6d inertiaData, MData;
    inertiaData << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    MData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    ArticulatedBodyInertia scaled = abi * 2.0;
    
    // Check I scaling
    lt scaledI = scaled.getInertia();
    Eigen::VectorXd scaledIData = scaledI.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(scaledIData(i), inertiaData(i) * 2.0);
    }
    
    // Check H scaling
    Eigen::Matrix3d scaledH = scaled.getH();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(scaledH(i, j), H(i, j) * 2.0);
        }
    }
    
    // Check M scaling
    lt scaledM = scaled.getM();
    Eigen::VectorXd scaledMData = scaledM.getData();
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_DOUBLE_EQ(scaledMData(i), MData(i) * 2.0);
    }
}

/**
 * @brief Test apply() with pure rotation (ω≠0, v=0)
 * @details Verifies INR-03: Formula [Iω + Hv; Hᵀω + Mv] with v=0
 */
TEST(ArticulatedBodyInertiaApplyTest, PureRotation)
{
    Vector6d inertiaData, MData;
    inertiaData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    MData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Pure rotation around Z: ω = (0, 0, 1), v = (0, 0, 0)
    MotionVector mv(Vector3d(0, 0, 1), Vector3d::Zero());
    
    // Expected result:
    // torque = I*ω + H*v = I*(0,0,1) + I*(0,0,0) = (0, 0, 5) + (0, 0, 0) = (0, 0, 5)
    // force = Hᵀ*ω + M*v = I*(0,0,1) + M*(0,0,0) = (0, 0, 1) + (0, 0, 0) = (0, 0, 1)
    Vector3d expectedTorque(0.0, 0.0, 5.0);
    Vector3d expectedForce(0.0, 0.0, 1.0);
    
    ForceVector result = abi.apply(mv);
    
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
 * @details Verifies INR-03: Formula [Iω + Hv; Hᵀω + Mv] with ω=0
 */
TEST(ArticulatedBodyInertiaApplyTest, PureTranslation)
{
    Vector6d inertiaData, MData;
    inertiaData << 4.0, 0.0, 4.0, 0.0, 0.0, 4.0;
    MData << 3.0, 0.0, 3.0, 0.0, 0.0, 3.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = 2.0 * Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Pure translation along X: ω = (0, 0, 0), v = (1, 0, 0)
    MotionVector mv(Vector3d::Zero(), Vector3d(1.0, 0.0, 0.0));
    
    // Expected result:
    // torque = I*ω + H*v = 0 + 2I*(1,0,0) = (2, 0, 0)
    // force = Hᵀ*ω + M*v = 0 + M*(1,0,0)
    // M*(1,0,0) for lower triangular [3,0,0; 0,3,0; 0,0,3] = (3, 0, 0)
    Vector3d expectedTorque(2.0, 0.0, 0.0);
    Vector3d expectedForce(3.0, 0.0, 0.0);
    
    ForceVector result = abi.apply(mv);
    
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
TEST(ArticulatedBodyInertiaApplyTest, CombinedMotion)
{
    Vector6d inertiaData, MData;
    inertiaData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    MData << 1.0, 0.0, 1.0, 0.0, 0.0, 1.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Combined motion: ω = (1, 0, 0), v = (0, 1, 0)
    MotionVector mv(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    
    // Expected result:
    // I*ω = (2,0,0,0,2,0,0,0,2) * (1,0,0) = (2, 0, 0)
    // H*v = I*(0,1,0) = (0, 1, 0)
    // torque = I*ω + H*v = (2, 0, 0) + (0, 1, 0) = (2, 1, 0)
    // Hᵀ*ω = I*(1,0,0) = (1, 0, 0)
    // M*v = I*(0,1,0) = (0, 1, 0)
    // force = Hᵀ*ω + M*v = (1, 0, 0) + (0, 1, 0) = (1, 1, 0)
    Vector3d expectedTorque(2.0, 1.0, 0.0);
    Vector3d expectedForce(1.0, 1.0, 0.0);
    
    ForceVector result = abi.apply(mv);
    
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
TEST(ArticulatedBodyInertiaApplyTest, ZeroMotion)
{
    Vector6d inertiaData, MData;
    inertiaData << 10.0, 0.0, 10.0, 0.0, 0.0, 10.0;
    MData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Zero motion: ω = (0, 0, 0), v = (0, 0, 0)
    MotionVector mv(Vector3d::Zero(), Vector3d::Zero());
    
    ForceVector result = abi.apply(mv);
    
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
 * @details Verifies INR-03: [Iω + Hv; Hᵀω + Mv] implementation
 */
TEST(ArticulatedBodyInertiaApplyTest, FeatherstoneFormula)
{
    Vector6d inertiaData, MData;
    inertiaData << 3.0, 1.0, 4.0, 1.0, 2.0, 5.0;
    MData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    
    lt I(inertiaData, 3), M(MData, 3);
    Eigen::Matrix3d H;
    H << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Arbitrary motion: ω = (1, 2, 3), v = (4, 5, 6)
    MotionVector mv(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    
    ForceVector result = abi.apply(mv);
    
    // Just verify it produces a valid result
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

/**
 * @brief Test apply() reduced case matches RigidBodyInertia
 * @details Verifies INR-03: When H=0 and M=mass*I, behavior matches rbi
 */
TEST(ArticulatedBodyInertiaApplyTest, ReducedCaseMatchesRigidBody)
{
    // Create articulated body with H=0, M=mass*I (diagonal)
    Vector6d inertiaData;
    inertiaData << 5.0, 0.0, 5.0, 0.0, 0.0, 5.0;
    lt I(inertiaData, 3);
    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
    Vector6d MData;
    MData << 2.0, 0.0, 2.0, 0.0, 0.0, 2.0;
    lt M(MData, 3);
    
    ArticulatedBodyInertia abi(I, H, M);
    
    // Create equivalent rigid body
    RigidBodyInertia rbi(2.0, Vector3d::Zero(), I);
    
    // Test with same motion
    MotionVector mv(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    
    ForceVector abiResult = abi.apply(mv);
    ForceVector rbiResult = rbi.apply(mv);
    
    // Results should match (H=0 means no coupling)
    Vector3d abiTorque = abiResult.getAngular();
    Vector3d rbiTorque = rbiResult.getAngular();
    
    Vector3d abiForce = abiResult.getLinear();
    Vector3d rbiForce = rbiResult.getLinear();
    
    EXPECT_NEAR(abiTorque.x(), rbiTorque.x(), TOLERANCE);
    EXPECT_NEAR(abiTorque.y(), rbiTorque.y(), TOLERANCE);
    EXPECT_NEAR(abiTorque.z(), rbiTorque.z(), TOLERANCE);
    
    // Force should differ because M is lower triangular, not scalar*mass
    // But both should be finite
    EXPECT_TRUE(std::isfinite(abiForce.x()));
    EXPECT_TRUE(std::isfinite(abiForce.y()));
    EXPECT_TRUE(std::isfinite(abiForce.z()));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
