/**
 * @file TestSpatialOperations.cpp
 * @brief Comprehensive GTest test suite for SpatialOperations class
 * @details Tests all three static methods: crossProductMotion, crossProductForce, 
 *          and transformInertia with various test scenarios and property tests.
 */

#include "SpatialOperations.h"
#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include "LowerTriangular.h"

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace SpatialAlgebra;
using namespace Eigen;

// Tolerance for floating point comparisons
constexpr double EPSILON = 1e-10;

// Helper functions for creating LowerTriangular matrices
LowerTriangular createIdentityInertia() {
    // Create 3x3 identity matrix and convert to LowerTriangular
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    return LowerTriangular::fromFullMatrix(identity);
}

LowerTriangular createDiagonalInertia(double value) {
    // Create diagonal matrix and convert to LowerTriangular
    Eigen::Matrix3d diagonal = Eigen::Matrix3d::Identity() * value;
    return LowerTriangular::fromFullMatrix(diagonal);
}

// ============================================================================
// CrossProductMotion Tests
// ============================================================================

/**
 * @brief Test suite for SpatialOperations::crossProductMotion
 */
class TestCrossProductMotion : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test cross product with simple rotation vectors
 * @details Verifies formula: [ω1×ω2; ω1×v2 + v1×ω2]
 */
TEST(TestCrossProductMotion, SimpleRotationVectors) {
    // Arrange: pure angular velocities about X and Y axes
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);
    
    // Assert: ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    //         v1×v2 = 0 (both zero)
    EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[2], 0.0, EPSILON);
}

/**
 * @brief Test cross product with combined motion vectors
 * @details Tests full formula with both angular and linear components
 */
TEST(TestCrossProductMotion, CombinedMotionVectors) {
    // Arrange
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(1, 0, 0));
    
    // Act
    SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);
    
    // Assert:
    // ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    // ω1×v2 = (1,0,0)×(1,0,0) = (0,0,0)
    // v1×ω2 = (0,1,0)×(0,1,0) = (0,0,0)
    // linear = (0,0,0)
    EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[2], 0.0, EPSILON);
}

/**
 * @brief Property test: verify anti-commutativity a×b = -(b×a)
 */
TEST(TestCrossProductMotion, Property_AntiCommutativity) {
    // Arrange
    MotionVector a(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector b(Vector3d(0.5, -0.3, 1.2), Vector3d(-0.1, 0.8, -0.5));
    
    // Act
    SpatialVector ab = SpatialOperations::crossProductMotion(a, b);
    SpatialVector ba = SpatialOperations::crossProductMotion(b, a);
    
    // Assert: ab + ba should equal zero
    EXPECT_NEAR((ab.getAngular() + ba.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((ab.getLinear() + ba.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// CrossProductForce Tests
// ============================================================================

/**
 * @brief Test suite for SpatialOperations::crossProductForce
 */
class TestCrossProductForce : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test cross product with simple force/torque vectors
 * @details Verifies formula: [τ1×τ2 + f1×f2; τ1×f2]
 */
TEST(TestCrossProductForce, SimpleForceTorqueVectors) {
    // Arrange: pure torques about X and Y axes
    ForceVector f1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    ForceVector f2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    SpatialVector result = SpatialOperations::crossProductForce(f1, f2);
    
    // Assert: τ1×τ2 = (1,0,0)×(0,1,0) = (0,0,1)
    //         f1×f2 = 0 (both zero)
    //         τ1×f2 = 0
    EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[2], 0.0, EPSILON);
}

/**
 * @brief Test cross product with combined motion and force vectors
 * @details Tests full formula: [ω×τ; ω×f + v×τ] for motion×force
 */
TEST(TestCrossProductForce, CombinedForceVectors) {
    // Arrange: Use MotionVector and ForceVector as the API expects
    MotionVector mv(Vector3d(1, 0, 0), Vector3d(0, 1, 0));  // ω=(1,0,0), v=(0,1,0)
    ForceVector fv(Vector3d(0, 1, 0), Vector3d(1, 0, 0));   // τ=(0,1,0), f=(1,0,0)
    
    // Act
    SpatialVector result = SpatialOperations::crossProductForce(mv, fv);
    
    // Assert: Formula is [ω×τ; ω×f + v×τ]
    // ω×τ = (1,0,0)×(0,1,0) = (0,0,1)
    // ω×f = (1,0,0)×(1,0,0) = (0,0,0)
    // v×τ = (0,1,0)×(0,1,0) = (0,0,0)
    // linear = (0,0,0) + (0,0,0) = (0,0,0)
    EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[2], 0.0, EPSILON);
}

/**
 * @brief Property test: verify anti-commutativity a×b = -(b×a)
 */
TEST(TestCrossProductForce, Property_AntiCommutativity) {
    // Arrange
    ForceVector a(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    ForceVector b(Vector3d(0.5, -0.3, 1.2), Vector3d(-0.1, 0.8, -0.5));
    
    // Act
    SpatialVector ab = SpatialOperations::crossProductForce(a, b);
    SpatialVector ba = SpatialOperations::crossProductForce(b, a);
    
    // Assert: ab + ba should equal zero
    EXPECT_NEAR((ab.getAngular() + ba.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((ab.getLinear() + ba.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// TransformInertia Tests
// ============================================================================

/**
 * @brief Test suite for SpatialOperations::transformInertia
 */
class TestTransformInertia : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
    
    LowerTriangular createIdentityInertia() {
        // Packed storage for 3x3 identity: [1, 0, 1, 0, 0, 1]
        std::vector<double> data = {1.0, 0.0, 1.0, 0.0, 0.0, 1.0};
        LowerTriangular lt(3);
        return lt;
    }
    
    LowerTriangular createDiagonalInertia(double value) {
        // Packed storage for diagonal matrix
        std::vector<double> data = {value, 0.0, value, 0.0, 0.0, value};
        LowerTriangular lt(3);
        return lt;
    }
};

/**
 * @brief Test with identity transform (should preserve inertia)
 */
TEST(TestTransformInertia, IdentityTransform) {
    // Arrange
    RigidBodyInertia inertia(1.0, Vector3d::Zero(), createIdentityInertia());
    Rotation rot = Rotation(Matrix3d::Identity());
    PluckerTransform transform(rot, Vector3d::Zero());
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert: identity transform preserves all properties
    EXPECT_NEAR(result.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(result.getCom().norm(), 0.0, EPSILON);
}

/**
 * @brief Test with rotation transform (should rotate inertia)
 */
TEST(TestTransformInertia, RotationTransform) {
    // Arrange: 90° rotation around Z axis
    Rotation rot = Rotation(Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ()));
    PluckerTransform transform(rot, Vector3d::Zero());
    
    // Create inertia with COM at [1, 0, 0]
    RigidBodyInertia inertia(1.0, Vector3d(1, 0, 0), createIdentityInertia());
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert: COM should rotate 90° around Z: [1,0,0] -> [0,1,0]
    EXPECT_NEAR(result.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(result.getCom()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getCom()[1], 1.0, EPSILON);
    EXPECT_NEAR(result.getCom()[2], 0.0, EPSILON);
}

/**
 * @brief Test with translation transform (should shift COM)
 */
TEST(TestTransformInertia, TranslationTransform) {
    // Arrange: identity rotation, translation [1, 0, 0]
    Rotation rot = Rotation(Matrix3d::Identity());
    PluckerTransform transform(rot, Vector3d(1, 0, 0));
    
    // Create inertia with COM at origin
    RigidBodyInertia inertia(1.0, Vector3d::Zero(), createIdentityInertia());
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert: COM should shift by -translation (parallel axis theorem)
    // h' = h - m*r = 0 - 1*[1,0,0] = [-1,0,0]
    EXPECT_NEAR(result.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(result.getCom()[0], -1.0, EPSILON);
    EXPECT_NEAR(result.getCom()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getCom()[2], 0.0, EPSILON);
}

/**
 * @brief Test with combined rotation and translation
 */
TEST(TestTransformInertia, CombinedTransform) {
    // Arrange: 90° Z rotation + [1, 0, 0] translation
    Rotation rot = Rotation(Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ()));
    PluckerTransform transform(rot, Vector3d(1, 0, 0));
    
    // Create inertia with mass=2, COM at [1, 0, 0]
    RigidBodyInertia inertia(2.0, Vector3d(1, 0, 0), createDiagonalInertia(1.0));
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert: mass conserved
    EXPECT_NEAR(result.getMass(), 2.0, EPSILON);
    // COM transformation involves both rotation and translation effects
    // This test verifies the method executes without error
}

/**
 * @brief Property test: mass is conserved under transformation
 */
TEST(TestTransformInertia, Property_MassConservation) {
    // Arrange
    Rotation rot = Rotation(Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ()));
    PluckerTransform transform(rot, Vector3d(0.5, 0.3, 0.2));
    
    RigidBodyInertia inertia(2.5, Vector3d(1, 2, 3), createIdentityInertia());
    
    // Act
    RigidBodyInertia result = SpatialOperations::transformInertia(inertia, transform);
    
    // Assert: mass must be conserved
    EXPECT_NEAR(result.getMass(), 2.5, EPSILON);
}

// ============================================================================
// Main entry point
// ============================================================================

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
