// TestRotation.cpp - Comprehensive GTest test suite for Rotation class

#include "Rotation.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>

using namespace SpatialAlgebra;
using namespace Eigen;

// Tolerance for floating point comparisons
const double TOLERANCE = 1e-10;

/**
 * @brief Test default constructor creates identity matrix
 * @details Verifies ROT-01: Default rotation should be identity (no rotation)
 */
TEST(RotationTest, DefaultConstructor)
{
    Rotation rot;
    
    // Check that it equals identity matrix
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == j)
            {
                EXPECT_DOUBLE_EQ(rot(i, j), 1.0);
            }
            else
            {
                EXPECT_DOUBLE_EQ(rot(i, j), 0.0);
            }
        }
    }
}

/**
 * @brief Test construction from 3x3 matrix
 * @details Verifies ROT-01: Matrix constructor preserves values
 */
TEST(RotationTest, FromMatrix)
{
    Matrix3d matrix;
    matrix << 0.0, -1.0, 0.0,
              1.0,  0.0, 0.0,
              0.0,  0.0, 1.0;
    
    Rotation rot(matrix);
    
    // Verify all elements are copied correctly
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(rot(i, j), matrix(i, j));
        }
    }
}

/**
 * @brief Test construction from angle-axis representation
 * @details Verifies ROT-02: 90-degree rotation around Z axis
 */
TEST(RotationTest, FromAngleAxis)
{
    // 90-degree rotation around Z axis
    AngleAxisd angleAxis(M_PI / 2.0, Vector3d::UnitZ());
    Rotation rot(angleAxis);
    
    // Expected rotation matrix for 90° around Z
    // [ cos(90°)  -sin(90°)  0 ]   [ 0  -1  0 ]
    // [ sin(90°)   cos(90°)  0 ] = [ 1   0  0 ]
    // [    0         0       1 ]   [ 0   0  1 ]
    EXPECT_NEAR(rot(0, 0), 0.0, TOLERANCE);
    EXPECT_NEAR(rot(0, 1), -1.0, TOLERANCE);
    EXPECT_NEAR(rot(0, 2), 0.0, TOLERANCE);
    EXPECT_NEAR(rot(1, 0), 1.0, TOLERANCE);
    EXPECT_NEAR(rot(1, 1), 0.0, TOLERANCE);
    EXPECT_NEAR(rot(1, 2), 0.0, TOLERANCE);
    EXPECT_DOUBLE_EQ(rot(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(rot(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(rot(2, 2), 1.0);
}

/**
 * @brief Test construction from quaternion representation
 * @details Verifies ROT-02: Quaternion to matrix conversion
 */
TEST(RotationTest, FromQuaternion)
{
    // 45-degree rotation around X axis
    AngleAxisd angleAxis(M_PI / 4.0, Vector3d::UnitX());
    Quaterniond quat(angleAxis);
    Rotation rot(quat);
    
    // Expected rotation matrix for 45° around X
    double cos45 = std::cos(M_PI / 4.0);
    double sin45 = std::sin(M_PI / 4.0);
    
    EXPECT_DOUBLE_EQ(rot(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(rot(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(rot(0, 2), 0.0);
    EXPECT_NEAR(rot(1, 1), cos45, TOLERANCE);
    EXPECT_NEAR(rot(1, 2), -sin45, TOLERANCE);
    EXPECT_NEAR(rot(2, 1), sin45, TOLERANCE);
    EXPECT_NEAR(rot(2, 2), cos45, TOLERANCE);
}

/**
 * @brief Test setFromAngleAxis setter
 * @details Verifies ROT-02: Setter updates rotation correctly
 */
TEST(RotationTest, SetFromAngleAxis)
{
    Rotation rot; // Start with identity
    
    // Set to 180-degree rotation around Y axis
    AngleAxisd angleAxis(M_PI, Vector3d::UnitY());
    rot.setFromAngleAxis(angleAxis);
    
    // Expected rotation matrix for 180° around Y
    // [ cos(180°)   0   sin(180°) ]   [ -1   0   0 ]
    // [     0       1      0      ] = [  0   1   0 ]
    // [ -sin(180°)  0   cos(180°) ]   [  0   0  -1 ]
    EXPECT_NEAR(rot(0, 0), -1.0, TOLERANCE);
    EXPECT_DOUBLE_EQ(rot(0, 1), 0.0);
    EXPECT_NEAR(rot(0, 2), 0.0, TOLERANCE);
    EXPECT_DOUBLE_EQ(rot(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(rot(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(rot(1, 2), 0.0);
    EXPECT_NEAR(rot(2, 0), 0.0, TOLERANCE);
    EXPECT_DOUBLE_EQ(rot(2, 1), 0.0);
    EXPECT_NEAR(rot(2, 2), -1.0, TOLERANCE);
}

/**
 * @brief Test setFromQuaternion setter
 * @details Verifies ROT-02: Quaternion setter updates rotation
 */
TEST(RotationTest, SetFromQuaternion)
{
    Rotation rot; // Start with identity
    
    // 90-degree rotation around X axis
    AngleAxisd angleAxis(M_PI / 2.0, Vector3d::UnitX());
    Quaterniond quat(angleAxis);
    rot.setFromQuaternion(quat);
    
    // Expected rotation matrix for 90° around X
    EXPECT_DOUBLE_EQ(rot(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(rot(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(rot(0, 2), 0.0);
    EXPECT_NEAR(rot(1, 1), 0.0, TOLERANCE);
    EXPECT_NEAR(rot(1, 2), -1.0, TOLERANCE);
    EXPECT_NEAR(rot(2, 1), 1.0, TOLERANCE);
    EXPECT_NEAR(rot(2, 2), 0.0, TOLERANCE);
}

/**
 * @brief Test conversion to angle-axis representation
 * @details Verifies ROT-03: Matrix to angle-axis conversion
 */
TEST(RotationTest, ToAngleAxis)
{
    // Create known rotation: 60° around Z
    double angle = M_PI / 3.0; // 60 degrees
    AngleAxisd original(angle, Vector3d::UnitZ());
    Rotation rot(original);
    
    // Convert back to angle-axis
    AngleAxisd recovered = rot.toAngleAxis();
    
    // Verify axis is Z (within tolerance for sign ambiguity)
    EXPECT_NEAR(std::abs(recovered.axis().dot(Vector3d::UnitZ())), 1.0, TOLERANCE);
    
    // Verify angle (may differ by 2π or sign)
    double recoveredAngle = std::abs(recovered.angle());
    double expectedAngle = angle;
    if (recoveredAngle > M_PI)
    {
        recoveredAngle = 2 * M_PI - recoveredAngle;
    }
    EXPECT_NEAR(recoveredAngle, expectedAngle, TOLERANCE);
}

/**
 * @brief Test conversion to quaternion representation
 * @details Verifies ROT-03: Matrix to quaternion conversion
 */
TEST(RotationTest, ToQuaternion)
{
    // Create known rotation: 120° around axis (1,1,1)/sqrt(3)
    Vector3d axis = Vector3d(1, 1, 1).normalized();
    double angle = 2.0 * M_PI / 3.0; // 120 degrees
    AngleAxisd original(angle, axis);
    Rotation rot(original);
    
    // Convert to quaternion
    Quaterniond quat = rot.toQuaternion();
    
    // Verify quaternion is normalized
    EXPECT_NEAR(quat.norm(), 1.0, TOLERANCE);
    
    // Convert back to matrix and verify
    Matrix3d recovered = quat.toRotationMatrix();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(rot(i, j), recovered(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test inverse operation
 * @details Verifies ROT-04: Inverse equals transpose for rotation matrices
 */
TEST(RotationTest, Inverse)
{
    // Create arbitrary rotation
    AngleAxisd angleAxis(M_PI / 6.0, Vector3d(1, 2, 3).normalized());
    Rotation rot(angleAxis);
    
    // Compute inverse
    Rotation inv = rot.inverse();
    
    // Verify R * R^(-1) = I
    Matrix3d product = static_cast<Matrix3d>(rot) * static_cast<Matrix3d>(inv);
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == j)
            {
                EXPECT_NEAR(product(i, j), 1.0, TOLERANCE);
            }
            else
            {
                EXPECT_NEAR(product(i, j), 0.0, TOLERANCE);
            }
        }
    }
}

/**
 * @brief Test transpose operation
 * @details Verifies ROT-04: Transpose equals inverse for rotation matrices
 */
TEST(RotationTest, Transpose)
{
    // Create arbitrary rotation
    AngleAxisd angleAxis(M_PI / 4.0, Vector3d::UnitY());
    Rotation rot(angleAxis);
    
    // Compute transpose
    Rotation trans = rot.transpose();
    
    // Verify transpose equals inverse
    Rotation inv = rot.inverse();
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(trans(i, j), inv(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test rotation composition (operator*)
 * @details Verifies ROT-04: Rotation multiplication composes correctly
 */
TEST(RotationTest, OperatorMultiplyRotation)
{
    // Two successive 90° rotations around Z should give 180°
    AngleAxisd rot90(M_PI / 2.0, Vector3d::UnitZ());
    Rotation r1(rot90);
    Rotation r2(rot90);
    
    Rotation combined = r1 * r2;
    
    // Expected: 180° rotation around Z
    AngleAxisd rot180(M_PI, Vector3d::UnitZ());
    Rotation expected(rot180);
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(combined(i, j), expected(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test vector transformation (operator*)
 * @details Verifies ROT-04: Vector rotation works correctly
 */
TEST(RotationTest, OperatorMultiplyVector)
{
    // 90° rotation around Z
    AngleAxisd angleAxis(M_PI / 2.0, Vector3d::UnitZ());
    Rotation rot(angleAxis);
    
    // Transform vector (1, 0, 0) - should become (0, 1, 0)
    Vector3d vec(1.0, 0.0, 0.0);
    Vector3d rotated = rot * vec;
    
    EXPECT_NEAR(rotated.x(), 0.0, TOLERANCE);
    EXPECT_NEAR(rotated.y(), 1.0, TOLERANCE);
    EXPECT_DOUBLE_EQ(rotated.z(), 0.0);
}

/**
 * @brief Test matrix multiplication (operator*)
 * @details Verifies matrix multiplication with general 3x3 matrices
 */
TEST(RotationTest, OperatorMultiplyMatrix)
{
    // 90° rotation around Z
    AngleAxisd angleAxis(M_PI / 2.0, Vector3d::UnitZ());
    Rotation rot(angleAxis);
    
    // Multiply with arbitrary matrix
    Matrix3d mat;
    mat << 1.0, 2.0, 3.0,
           4.0, 5.0, 6.0,
           7.0, 8.0, 9.0;
    
    Matrix3d result = rot * mat;
    
    // Expected result: R * M
    Matrix3d expected = static_cast<Matrix3d>(rot) * mat;
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(result(i, j), expected(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test orthogonality preservation
 * @details Verifies ROT-04: R * R^T = I for all rotation matrices
 */
TEST(RotationTest, Orthogonality)
{
    // Test with various rotations
    std::vector<AngleAxisd> testRotations = {
        AngleAxisd(M_PI / 6.0, Vector3d::UnitX()),
        AngleAxisd(M_PI / 4.0, Vector3d::UnitY()),
        AngleAxisd(M_PI / 3.0, Vector3d::UnitZ()),
        AngleAxisd(M_PI / 2.0, Vector3d(1, 1, 1).normalized())
    };
    
    for (const auto &aa : testRotations)
    {
        Rotation rot(aa);
        Rotation trans = rot.transpose();
        
        // Verify R * R^T = I
        Matrix3d product = static_cast<Matrix3d>(rot) * static_cast<Matrix3d>(trans);
        
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (i == j)
                {
                    EXPECT_NEAR(product(i, j), 1.0, TOLERANCE) << "Failed orthogonality check for rotation";
                }
                else
                {
                    EXPECT_NEAR(product(i, j), 0.0, TOLERANCE) << "Failed orthogonality check for rotation";
                }
            }
        }
    }
}

/**
 * @brief Test determinant equals +1 for valid rotations
 * @details Verifies ROT-04: Proper rotation matrices have determinant +1
 */
TEST(RotationTest, Determinant)
{
    // Test with various rotations
    std::vector<AngleAxisd> testRotations = {
        AngleAxisd(0.0, Vector3d::UnitX()),           // Identity
        AngleAxisd(M_PI / 4.0, Vector3d::UnitY()),    // 45° around Y
        AngleAxisd(M_PI / 2.0, Vector3d::UnitZ()),    // 90° around Z
        AngleAxisd(M_PI, Vector3d(1, 2, 3).normalized()) // 180° around arbitrary axis
    };
    
    for (const auto &aa : testRotations)
    {
        Rotation rot(aa);
        double det = rot.determinant();
        EXPECT_NEAR(det, 1.0, TOLERANCE) << "Determinant should be +1 for valid rotation matrix";
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
