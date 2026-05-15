/**
 * @file TestSpatialUtils.cpp
 * @brief Comprehensive test suite for SpatialUtils and SpatialOperations
 */

#include "SpatialUtils.h"
#include "SpatialOperations.h"
#include <gtest/gtest.h>

using namespace SpatialAlgebra;

/**
 * @brief Test suite for skew-symmetric matrix function
 */
class TestSkew : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Verify skew() creates correct skew-symmetric matrix
 * @details Tests that skew(v) produces matrix S where S*x = v×x
 */
TEST(TestSkew, CreatesSkewSymmetricMatrix) {
    // Arrange
    Vector3d v(1.0, 2.0, 3.0);
    
    // Act
    Eigen::Matrix3d S = skew(v);
    
    // Assert - verify structure: [0, -vz, vy; vz, 0, -vx; -vy, vx, 0]
    EXPECT_DOUBLE_EQ(S(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(S(0, 1), -3.0);
    EXPECT_DOUBLE_EQ(S(0, 2), 2.0);
    EXPECT_DOUBLE_EQ(S(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(S(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(S(1, 2), -1.0);
    EXPECT_DOUBLE_EQ(S(2, 0), -2.0);
    EXPECT_DOUBLE_EQ(S(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(S(2, 2), 0.0);
}

/**
 * @brief Property test: verify skew-symmetry property S + S^T = 0
 */
TEST(TestSkew, Property_SkewSymmetric) {
    // Arrange
    Vector3d v(1.0, 2.0, 3.0);
    Eigen::Matrix3d S = skew(v);
    
    // Act
    Eigen::Matrix3d sum = S + S.transpose();
    
    // Assert - should be zero matrix
    EXPECT_NEAR(sum.norm(), 0.0, 1e-10);
}

/**
 * @brief Test suite for dot product functions
 */
class TestDot : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test SpatialVector dot product
 * @details Verifies: a·b = Σ(angular_i * angular_i) + Σ(linear_i * linear_i)
 */
TEST(TestDot, SpatialVector_SpatialVector) {
    // Arrange
    MotionVector v1(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector v2(Vector3d(2, 3, 4), Vector3d(5, 6, 7));
    
    // Act
    double result = dot(v1, v2);
    
    // Assert: 1*2+2*3+3*4 + 4*5+5*6+6*7 = 20 + 92 = 112
    EXPECT_DOUBLE_EQ(result, 112.0);
}

/**
 * @brief Test MotionVector dot product overload
 */
TEST(TestDot, MotionVector_MotionVector) {
    // Arrange
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    double result = dot(v1, v2);
    
    // Assert: 1*0 + 0*1 + 0*0 + 0*1 + 1*0 + 0*0 = 0
    EXPECT_DOUBLE_EQ(result, 0.0);
}

/**
 * @brief Test ForceVector dot product overload
 */
TEST(TestDot, ForceVector_ForceVector) {
    // Arrange
    ForceVector f1(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    ForceVector f2(Vector3d(2, 3, 4), Vector3d(5, 6, 7));
    
    // Act
    double result = dot(f1, f2);
    
    // Assert: same as SpatialVector = 112
    EXPECT_DOUBLE_EQ(result, 112.0);
}

/**
 * @brief Test MotionVector·ForceVector dot product
 * @details Physical interpretation: power = ω·τ + v·f
 */
TEST(TestDot, MotionVector_ForceVector) {
    // Arrange
    MotionVector mv(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    ForceVector fv(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    double result = dot(mv, fv);
    
    // Assert: 1*0 + 0*1 + 0*0 + 0*1 + 1*0 + 0*0 = 0
    EXPECT_DOUBLE_EQ(result, 0.0);
}

/**
 * @brief Property test: verify dot product commutativity a·b = b·a
 */
TEST(TestDot, Property_Commutativity) {
    // Arrange
    MotionVector v1(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector v2(Vector3d(2, 3, 4), Vector3d(5, 6, 7));
    
    // Act
    double ab = dot(v1, v2);
    double ba = dot(v2, v1);
    
    // Assert
    EXPECT_DOUBLE_EQ(ab, ba);
}

/**
 * @brief Test suite for cross product functions
 */
class TestCross : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test motion cross motion product
 * @details Formula: [ω1×ω2; ω1×v2 + v1×ω2]
 */
TEST(TestCross, Motion_Motion) {
    // Arrange
    MotionVector mv1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    MotionVector mv2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    MotionVector result = cross(mv1, mv2);
    
    // Assert: ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    //         v1×ω2 + ω1×v2 = 0 + 0 = 0
    EXPECT_DOUBLE_EQ(result.getAngular()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(2), 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(2), 0.0);
}

/**
 * @brief Test force cross force product
 * @details Formula: [τ1×τ2 + f1×f2; τ1×f2]
 */
TEST(TestCross, Force_Force) {
    // Arrange
    ForceVector fv1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    ForceVector fv2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    ForceVector result = cross(fv1, fv2);
    
    // Assert: τ1×τ2 = (1,0,0)×(0,1,0) = (0,0,1)
    //         f1×f2 = 0
    //         τ1×f2 = 0
    EXPECT_DOUBLE_EQ(result.getAngular()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(2), 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(2), 0.0);
}

/**
 * @brief Test motion cross force product
 * @details Formula: [ω1×f1; ω1×f2 + v1×f1]
 */
TEST(TestCross, Motion_Force) {
    // Arrange
    MotionVector mv(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    ForceVector fv(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    ForceVector result = cross(mv, fv);
    
    // Assert: ω×f = (1,0,0)×(0,1,0) = (0,0,1)
    //         ω×f_linear + v×f = 0 + 0 = 0
    EXPECT_DOUBLE_EQ(result.getAngular()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()(2), 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(0), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(1), 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()(2), 0.0);
}

/**
 * @brief Property test: verify cross product anti-commutativity a×b = -(b×a)
 */
TEST(TestCross, Property_AntiCommutativity_Motion) {
    // Arrange
    MotionVector mv1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    MotionVector mv2(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    MotionVector ab = cross(mv1, mv2);
    MotionVector ba = cross(mv2, mv1);
    
    // Assert: ab should equal -ba
    EXPECT_NEAR((ab + ba).getAngular().norm(), 0.0, 1e-10);
    EXPECT_NEAR((ab + ba).getLinear().norm(), 0.0, 1e-10);
}

/**
 * @brief Property test: verify force cross product anti-commutativity
 */
TEST(TestCross, Property_AntiCommutativity_Force) {
    // Arrange
    ForceVector fv1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    ForceVector fv2(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    ForceVector ab = cross(fv1, fv2);
    ForceVector ba = cross(fv2, fv1);
    
    // Assert: ab should equal -ba
    EXPECT_NEAR((ab + ba).getAngular().norm(), 0.0, 1e-10);
    EXPECT_NEAR((ab + ba).getLinear().norm(), 0.0, 1e-10);
}

/**
 * @brief Test suite for SpatialOperations static class
 */
class TestSpatialOperations : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test SpatialOperations::crossProductMotion
 * @details Verify it delegates correctly to underlying cross product
 */
TEST(TestSpatialOperations, CrossProductMotion) {
    // Arrange
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    SpatialVector result_static = SpatialOperations::crossProductMotion(v1, v2);
    MotionVector result_free = cross(v1, v2);
    
    // Assert - static method should match free function
    EXPECT_NEAR((result_static.getAngular() - result_free.getAngular()).norm(), 0.0, 1e-10);
    EXPECT_NEAR((result_static.getLinear() - result_free.getLinear()).norm(), 0.0, 1e-10);
}

/**
 * @brief Test SpatialOperations::crossProductForce
 * @details Verify it delegates correctly to underlying cross product
 */
TEST(TestSpatialOperations, CrossProductForce) {
    // Arrange
    MotionVector mv(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    ForceVector fv(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    SpatialVector result_static = SpatialOperations::crossProductForce(mv, fv);
    ForceVector result_free = cross(mv, fv);
    
    // Assert - static method should match free function
    EXPECT_NEAR((result_static.getAngular() - result_free.getAngular()).norm(), 0.0, 1e-10);
    EXPECT_NEAR((result_static.getLinear() - result_free.getLinear()).norm(), 0.0, 1e-10);
}

/**
 * @brief Test SpatialOperations::transformInertia
 * @details Verify inertia transformation using Plucker transform
 */
TEST(TestSpatialOperations, TransformInertia) {
    // Arrange
    // Create a simple rigid body: mass=1, COM at origin, identity inertia
    // LowerTriangular stores packed data - for 3x3 identity, need 6 elements
    // Packed storage: [L00, L10, L11, L20, L21, L22]
    std::vector<double> identityData = {1.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    LowerTriangular I_lt(3);
    // Use operator[] to access packed data
    // For 3x3: indices 0,1,2,3,4,5 correspond to (0,0), (1,0), (1,1), (2,0), (2,1), (2,2)
    
    RigidBodyInertia inertia(1.0, Vector3d::Zero(), I_lt);
    
    // Identity transform (no rotation, no translation)
    Rotation rot = Rotation(Eigen::AngleAxisd(0.0, Vector3d::UnitX()));
    PluckerTransform transform(rot, Vector3d::Zero());
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert - identity transform should preserve inertia
    EXPECT_DOUBLE_EQ(result.getMass(), 1.0);
    EXPECT_NEAR(result.getCom().norm(), 0.0, 1e-10);
}

/**
 * @brief Main entry point for GTest
 */
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
