#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include <gtest/gtest.h>

using namespace SpatialAlgebra;

// ============================================================================
// SpatialVector Tests
// ============================================================================

/**
 * @brief Test default constructor creates zero vector
 */
TEST(TestSpatialVector, Constructor)
{
    // Default constructor should create zero vector
    SpatialVector zero;
    EXPECT_DOUBLE_EQ(zero.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(zero.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(zero.getAngular()[2], 0.0);
    EXPECT_DOUBLE_EQ(zero.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(zero.getLinear()[1], 0.0);
    EXPECT_DOUBLE_EQ(zero.getLinear()[2], 0.0);

    // Constructor with components
    Vector3d angular(1.0, 2.0, 3.0);
    Vector3d linear(4.0, 5.0, 6.0);
    SpatialVector v1(angular, linear);

    EXPECT_DOUBLE_EQ(v1.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(v1.getAngular()[1], 2.0);
    EXPECT_DOUBLE_EQ(v1.getAngular()[2], 3.0);
    EXPECT_DOUBLE_EQ(v1.getLinear()[0], 4.0);
    EXPECT_DOUBLE_EQ(v1.getLinear()[1], 5.0);
    EXPECT_DOUBLE_EQ(v1.getLinear()[2], 6.0);
}

/**
 * @brief Test getters return correct components
 */
TEST(TestSpatialVector, Getters)
{
    Vector3d angular(1.0, 2.0, 3.0);
    Vector3d linear(4.0, 5.0, 6.0);
    SpatialVector v1(angular, linear);

    Vector3d retrievedAngular = v1.getAngular();
    Vector3d retrievedLinear = v1.getLinear();

    EXPECT_DOUBLE_EQ(retrievedAngular[0], 1.0);
    EXPECT_DOUBLE_EQ(retrievedAngular[1], 2.0);
    EXPECT_DOUBLE_EQ(retrievedAngular[2], 3.0);
    EXPECT_DOUBLE_EQ(retrievedLinear[0], 4.0);
    EXPECT_DOUBLE_EQ(retrievedLinear[1], 5.0);
    EXPECT_DOUBLE_EQ(retrievedLinear[2], 6.0);
}

/**
 * @brief Test component-wise addition
 */
TEST(TestSpatialVector, Addition)
{
    SpatialVector v1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    SpatialVector v2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    SpatialVector sum = v1 + v2;

    EXPECT_DOUBLE_EQ(sum.getAngular()[0], 3.0);
    EXPECT_DOUBLE_EQ(sum.getAngular()[1], 5.0);
    EXPECT_DOUBLE_EQ(sum.getAngular()[2], 7.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[0], 9.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[1], 11.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[2], 13.0);
}

/**
 * @brief Test component-wise subtraction
 */
TEST(TestSpatialVector, Subtraction)
{
    SpatialVector v1(Vector3d(5.0, 7.0, 9.0), Vector3d(10.0, 12.0, 14.0));
    SpatialVector v2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    SpatialVector diff = v1 - v2;

    EXPECT_DOUBLE_EQ(diff.getAngular()[0], 3.0);
    EXPECT_DOUBLE_EQ(diff.getAngular()[1], 4.0);
    EXPECT_DOUBLE_EQ(diff.getAngular()[2], 5.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[0], 5.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[1], 6.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[2], 7.0);
}

/**
 * @brief Test scalar multiplication
 */
TEST(TestSpatialVector, ScalarMultiplication)
{
    SpatialVector v1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));

    SpatialVector scaled = v1 * 2.5;

    EXPECT_DOUBLE_EQ(scaled.getAngular()[0], 2.5);
    EXPECT_DOUBLE_EQ(scaled.getAngular()[1], 5.0);
    EXPECT_DOUBLE_EQ(scaled.getAngular()[2], 7.5);
    EXPECT_DOUBLE_EQ(scaled.getLinear()[0], 10.0);
    EXPECT_DOUBLE_EQ(scaled.getLinear()[1], 12.5);
    EXPECT_DOUBLE_EQ(scaled.getLinear()[2], 15.0);
}

/**
 * @brief Test dot product: ω1·ω2 + v1·v2
 */
TEST(TestSpatialVector, DotProduct)
{
    SpatialVector v1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    SpatialVector v2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    // ω1·ω2 = 1*2 + 2*3 + 3*4 = 2 + 6 + 12 = 20
    // v1·v2 = 4*5 + 5*6 + 6*7 = 20 + 30 + 42 = 92
    // Total = 20 + 92 = 112
    double dotProduct = v1.dot(v2);

    EXPECT_DOUBLE_EQ(dotProduct, 112.0);
}

/**
 * @brief Property test: Cross product anti-commutativity a×b = -(b×a)
 * @details Verifies the fundamental mathematical property that cross product
 *          is anti-commutative. This is a critical invariant for spatial algebra.
 * @see Featherstone 2008, Chapter 2
 */
TEST(TestSpatialVector, CrossProductAntiCommutativity)
{
    SpatialVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    SpatialVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    SpatialVector a_cross_b = a.crossMotion(b);
    SpatialVector b_cross_a = b.crossMotion(a);
    SpatialVector neg_b_cross_a = b_cross_a * -1.0;

    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[0], neg_b_cross_a.getAngular()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[1], neg_b_cross_a.getAngular()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[2], neg_b_cross_a.getAngular()[2]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[0], neg_b_cross_a.getLinear()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[1], neg_b_cross_a.getLinear()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[2], neg_b_cross_a.getLinear()[2]);
}

/**
 * @brief Property test: Dot product commutativity a·b = b·a
 * @details Verifies that the spatial dot product is commutative,
 *          a fundamental property of inner product spaces.
 */
TEST(TestSpatialVector, DotProductCommutativity)
{
    SpatialVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    SpatialVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    double a_dot_b = a.dot(b);
    double b_dot_a = b.dot(a);

    EXPECT_DOUBLE_EQ(a_dot_b, b_dot_a);
}

/**
 * @brief Textbook example: Featherstone 2008, Chapter 2, Example 2.1
 * @details Representative spatial vector example showing rotation about X
 *          combined with translation along Y, and rotation about Y with
 *          translation along Z. Verifies cross product calculation.
 */
TEST(TestSpatialVector, FeatherstoneExample2_1)
{
    // v1: rotation about X (1,0,0), translation along Y (0,1,0)
    SpatialVector v1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    // v2: rotation about Y (0,1,0), translation along Z (0,0,1)
    SpatialVector v2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));

    // Cross product: [ω1×ω2; ω1×v2 + v1×ω2]
    // ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    // ω1×v2 = (1,0,0)×(0,0,1) = (0,-1,0)
    // v1×ω2 = (0,1,0)×(0,1,0) = (0,0,0)
    // Linear = (0,-1,0) + (0,0,0) = (0,-1,0)
    // Result: [(0,0,1); (0,-1,0)]
    SpatialVector result = v1.crossMotion(v2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], -1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}

// ============================================================================
// MotionVector Tests
// ============================================================================

/**
 * @brief Test MotionVector constructors
 */
TEST(TestMotionVector, Constructor)
{
    // Default constructor
    MotionVector zero;
    EXPECT_DOUBLE_EQ(zero.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(zero.getLinear()[0], 0.0);

    // Constructor with components
    MotionVector mv1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    EXPECT_DOUBLE_EQ(mv1.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(mv1.getLinear()[1], 5.0);

    // Constructor from SpatialVector
    SpatialVector sv(Vector3d(2.0, 4.0, 6.0), Vector3d(8.0, 10.0, 12.0));
    MotionVector mv2(sv);
    EXPECT_DOUBLE_EQ(mv2.getAngular()[0], 2.0);
    EXPECT_DOUBLE_EQ(mv2.getLinear()[1], 10.0);
}

/**
 * @brief Test MotionVector getters
 */
TEST(TestMotionVector, Getters)
{
    MotionVector mv(Vector3d(3.0, 6.0, 9.0), Vector3d(12.0, 15.0, 18.0));

    Vector3d angular = mv.getAngular();
    Vector3d linear = mv.getLinear();

    EXPECT_DOUBLE_EQ(angular[0], 3.0);
    EXPECT_DOUBLE_EQ(angular[1], 6.0);
    EXPECT_DOUBLE_EQ(angular[2], 9.0);
    EXPECT_DOUBLE_EQ(linear[0], 12.0);
    EXPECT_DOUBLE_EQ(linear[1], 15.0);
    EXPECT_DOUBLE_EQ(linear[2], 18.0);
}

/**
 * @brief Test MotionVector arithmetic operations
 */
TEST(TestMotionVector, Operations)
{
    MotionVector mv1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    MotionVector mv2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));

    // Addition
    MotionVector sum = mv1 + mv2;
    EXPECT_DOUBLE_EQ(sum.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(sum.getAngular()[1], 1.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[1], 1.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[2], 1.0);

    // Subtraction
    MotionVector diff = mv1 - mv2;
    EXPECT_DOUBLE_EQ(diff.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(diff.getAngular()[1], -1.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[1], 1.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[2], -1.0);

    // Scalar multiplication
    MotionVector scaled = mv1 * 3.0;
    EXPECT_DOUBLE_EQ(scaled.getAngular()[0], 3.0);
    EXPECT_DOUBLE_EQ(scaled.getLinear()[1], 3.0);
}

/**
 * @brief Test MotionVector crossMotion operation
 * Formula: [ω1×ω2; ω1×v2 + v1×ω2]
 */
TEST(TestMotionVector, CrossMotion)
{
    MotionVector mv1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0));
    MotionVector mv2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 0.0));

    // ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    // ω1×v2 = (1,0,0)×(0,0,0) = (0,0,0)
    // v1×ω2 = (0,0,0)×(0,1,0) = (0,0,0)
    // Result: [(0,0,1); (0,0,0)]
    MotionVector result = mv1.crossMotion(mv2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}

/**
 * @brief Test MotionVector crossMotion anti-commutativity: a×b = -(b×a)
 */
TEST(TestMotionVector, CrossMotionAntiCommutativity)
{
    MotionVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    MotionVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    MotionVector a_cross_b = a.crossMotion(b);
    MotionVector b_cross_a = b.crossMotion(a);
    MotionVector neg_b_cross_a = b_cross_a * -1.0;

    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[0], neg_b_cross_a.getAngular()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[1], neg_b_cross_a.getAngular()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[2], neg_b_cross_a.getAngular()[2]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[0], neg_b_cross_a.getLinear()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[1], neg_b_cross_a.getLinear()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[2], neg_b_cross_a.getLinear()[2]);
}

/**
 * @brief Test MotionVector crossMotion with non-zero linear components (bug detection)
 * Formula: [ω1×ω2; ω1×v2 + v1×ω2]
 * This test exposes the bug where linear.cross(other.linear) is used incorrectly
 */
TEST(TestMotionVector, CrossMotionWithLinearComponents)
{
    // ω1 = (1,0,0), v1 = (0,1,0)
    // ω2 = (0,1,0), v2 = (0,0,1)
    MotionVector mv1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    MotionVector mv2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));

    // Expected result using correct formula:
    // ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    // ω1×v2 = (1,0,0)×(0,0,1) = (0,-1,0)
    // v1×ω2 = (0,1,0)×(0,1,0) = (0,0,0)
    // Linear result = (0,-1,0) + (0,0,0) = (0,-1,0)
    MotionVector result = mv1.crossMotion(mv2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], -1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}

/**
 * @brief Test MotionVector dot product
 */
TEST(TestMotionVector, DotProduct)
{
    MotionVector mv1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    MotionVector mv2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    // Same as SpatialVector: ω1·ω2 + v1·v2 = 112
    double dotProduct = mv1.dot(mv2);

    EXPECT_DOUBLE_EQ(dotProduct, 112.0);
}

/**
 * @brief Property test: Cross product distributivity a×(b+c) = a×b + a×c
 * @details Verifies that cross product distributes over vector addition,
 *          a fundamental property of spatial algebra.
 * @see Featherstone 2008, Chapter 2
 */
TEST(TestMotionVector, CrossProductDistributivity)
{
    MotionVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    MotionVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));
    MotionVector c(Vector3d(3.0, 4.0, 5.0), Vector3d(6.0, 7.0, 8.0));

    // Left side: a × (b + c)
    MotionVector b_plus_c = b + c;
    MotionVector left = a.crossMotion(b_plus_c);

    // Right side: a × b + a × c
    MotionVector a_cross_b = a.crossMotion(b);
    MotionVector a_cross_c = a.crossMotion(c);
    MotionVector right = a_cross_b + a_cross_c;

    EXPECT_DOUBLE_EQ(left.getAngular()[0], right.getAngular()[0]);
    EXPECT_DOUBLE_EQ(left.getAngular()[1], right.getAngular()[1]);
    EXPECT_DOUBLE_EQ(left.getAngular()[2], right.getAngular()[2]);
    EXPECT_DOUBLE_EQ(left.getLinear()[0], right.getLinear()[0]);
    EXPECT_DOUBLE_EQ(left.getLinear()[1], right.getLinear()[1]);
    EXPECT_DOUBLE_EQ(left.getLinear()[2], right.getLinear()[2]);
}

/**
 * @brief Property test: Scalar multiplication property (k*a)×b = k*(a×b)
 * @details Verifies that scalar multiplication commutes with cross product,
 *          demonstrating linearity of the cross product operation.
 */
TEST(TestMotionVector, ScalarMultiplicationProperty)
{
    double k = 2.5;
    MotionVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    MotionVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    // Left side: (k*a) × b
    MotionVector k_a = a * k;
    MotionVector left = k_a.crossMotion(b);

    // Right side: k * (a × b)
    MotionVector a_cross_b = a.crossMotion(b);
    MotionVector right = a_cross_b * k;

    EXPECT_DOUBLE_EQ(left.getAngular()[0], right.getAngular()[0]);
    EXPECT_DOUBLE_EQ(left.getAngular()[1], right.getAngular()[1]);
    EXPECT_DOUBLE_EQ(left.getAngular()[2], right.getAngular()[2]);
    EXPECT_DOUBLE_EQ(left.getLinear()[0], right.getLinear()[0]);
    EXPECT_DOUBLE_EQ(left.getLinear()[1], right.getLinear()[1]);
    EXPECT_DOUBLE_EQ(left.getLinear()[2], right.getLinear()[2]);
}

/**
 * @brief Textbook example: Twist representing screw motion
 * @details Creates a twist with rotation about Z-axis and translation along X-axis,
 *          representing a screw motion. Verifies cross product with another twist.
 * @see Featherstone 2008, Chapter 2 - Screw Theory
 */
TEST(TestMotionVector, TwistExample)
{
    // Twist: ω = (0, 0, 1) rad/s (rotation about Z)
    //        v = (1, 0, 0) m/s (translation along X)
    MotionVector twist1(Vector3d(0.0, 0.0, 1.0), Vector3d(1.0, 0.0, 0.0));
    
    // Second twist: rotation about X, translation along Y
    MotionVector twist2(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));

    // Cross product: [ω1×ω2; ω1×v2 + v1×ω2]
    // ω1×ω2 = (0,0,1)×(1,0,0) = (0,1,0)
    // ω1×v2 = (0,0,1)×(0,1,0) = (-1,0,0)
    // v1×ω2 = (1,0,0)×(1,0,0) = (0,0,0)
    // Linear = (-1,0,0) + (0,0,0) = (-1,0,0)
    // Result: [(0,1,0); (-1,0,0)]
    MotionVector result = twist1.crossMotion(twist2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 1.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], -1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}

// ============================================================================
// ForceVector Tests
// ============================================================================

/**
 * @brief Test ForceVector constructors
 */
TEST(TestForceVector, Constructor)
{
    // Default constructor
    ForceVector zero;
    EXPECT_DOUBLE_EQ(zero.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(zero.getLinear()[0], 0.0);

    // Constructor with components
    ForceVector fv1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    EXPECT_DOUBLE_EQ(fv1.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(fv1.getLinear()[1], 5.0);

    // Constructor from SpatialVector
    SpatialVector sv(Vector3d(2.0, 4.0, 6.0), Vector3d(8.0, 10.0, 12.0));
    ForceVector fv2(sv);
    EXPECT_DOUBLE_EQ(fv2.getAngular()[0], 2.0);
    EXPECT_DOUBLE_EQ(fv2.getLinear()[1], 10.0);
}

/**
 * @brief Test ForceVector getters
 */
TEST(TestForceVector, Getters)
{
    ForceVector fv(Vector3d(3.0, 6.0, 9.0), Vector3d(12.0, 15.0, 18.0));

    Vector3d torque = fv.getAngular();
    Vector3d force = fv.getLinear();

    EXPECT_DOUBLE_EQ(torque[0], 3.0);
    EXPECT_DOUBLE_EQ(torque[1], 6.0);
    EXPECT_DOUBLE_EQ(torque[2], 9.0);
    EXPECT_DOUBLE_EQ(force[0], 12.0);
    EXPECT_DOUBLE_EQ(force[1], 15.0);
    EXPECT_DOUBLE_EQ(force[2], 18.0);
}

/**
 * @brief Test ForceVector arithmetic operations
 */
TEST(TestForceVector, Operations)
{
    ForceVector fv1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    ForceVector fv2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));

    // Addition
    ForceVector sum = fv1 + fv2;
    EXPECT_DOUBLE_EQ(sum.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(sum.getAngular()[1], 1.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[1], 1.0);
    EXPECT_DOUBLE_EQ(sum.getLinear()[2], 1.0);

    // Subtraction
    ForceVector diff = fv1 - fv2;
    EXPECT_DOUBLE_EQ(diff.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(diff.getAngular()[1], -1.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[1], 1.0);
    EXPECT_DOUBLE_EQ(diff.getLinear()[2], -1.0);

    // Scalar multiplication
    ForceVector scaled = fv1 * 3.0;
    EXPECT_DOUBLE_EQ(scaled.getAngular()[0], 3.0);
    EXPECT_DOUBLE_EQ(scaled.getLinear()[1], 3.0);
}

/**
 * @brief Test ForceVector crossForce operation
 * @details Formula: [τ1×τ2; τ1×f2 - τ2×f1]
 *          This test verifies the basic crossForce operation with pure torques.
 */
TEST(TestForceVector, CrossForce)
{
    ForceVector fv1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0));
    ForceVector fv2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 0.0));

    // τ1×τ2 = (1,0,0)×(0,1,0) = (0,0,1)
    // τ1×f2 - τ2×f1 = (1,0,0)×(0,0,0) - (0,1,0)×(0,0,0) = (0,0,0)
    // Result: [(0,0,1); (0,0,0)]
    ForceVector result = fv1.crossForce(fv2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}

/**
 * @brief Test ForceVector dot product
 */
TEST(TestForceVector, DotProduct)
{
    ForceVector fv1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    ForceVector fv2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    // Same as SpatialVector: τ1·τ2 + f1·f2 = 112
    double dotProduct = fv1.dot(fv2);

    EXPECT_DOUBLE_EQ(dotProduct, 112.0);
}

/**
 * @brief Property test: Cross force anti-commutativity a×b = -(b×a)
 * @details Verifies that the force cross product is anti-commutative,
 *          a fundamental property of spatial force algebra.
 * @see Featherstone 2008, Chapter 2
 */
TEST(TestForceVector, CrossForceAntiCommutativity)
{
    ForceVector a(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    ForceVector b(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));

    ForceVector a_cross_b = a.crossForce(b);
    ForceVector b_cross_a = b.crossForce(a);
    ForceVector neg_b_cross_a = b_cross_a * -1.0;

    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[0], neg_b_cross_a.getAngular()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[1], neg_b_cross_a.getAngular()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[2], neg_b_cross_a.getAngular()[2]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[0], neg_b_cross_a.getLinear()[0]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[1], neg_b_cross_a.getLinear()[1]);
    EXPECT_DOUBLE_EQ(a_cross_b.getLinear()[2], neg_b_cross_a.getLinear()[2]);
}

/**
 * @brief Textbook example: Wrench representing force and torque
 * @details Creates a wrench with torque about X-axis and force along Y-axis,
 *          representing a typical wrench in robotic manipulation.
 * @see Featherstone 2008, Chapter 2 - Wrenches
 */
TEST(TestForceVector, WrenchExample)
{
    // Wrench: τ = (1, 0, 0) N⋅m (torque about X)
    //         f = (0, 1, 0) N (force along Y)
    ForceVector wrench1(Vector3d(1.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0));
    
    // Second wrench: torque about Y, force along Z
    ForceVector wrench2(Vector3d(0.0, 1.0, 0.0), Vector3d(0.0, 0.0, 1.0));

    // Cross product: [τ1×τ2; τ1×f2 - τ2×f1]
    // τ1×τ2 = (1,0,0)×(0,1,0) = (0,0,1)
    // τ1×f2 = (1,0,0)×(0,0,1) = (0,-1,0)
    // τ2×f1 = (0,1,0)×(0,1,0) = (0,0,0)
    // Linear = (0,-1,0) - (0,0,0) = (0,-1,0)
    // Result: [(0,0,1); (0,-1,0)]
    ForceVector result = wrench1.crossForce(wrench2);

    EXPECT_DOUBLE_EQ(result.getAngular()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[1], 0.0);
    EXPECT_DOUBLE_EQ(result.getAngular()[2], 1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[0], 0.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[1], -1.0);
    EXPECT_DOUBLE_EQ(result.getLinear()[2], 0.0);
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
