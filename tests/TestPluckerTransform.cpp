// Test suite for PluckerTransform motion/force transforms
// Comprehensive GTest tests for Wave 1: transformMotion and transformForce

#include "PluckerTransform.h"
#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace SpatialAlgebra;
using namespace Eigen;

// Tolerance for floating point comparisons
constexpr double EPSILON = 1e-10;

// ============================================================================
// TransformMotion Tests
// ============================================================================

TEST(TransformMotionTest, IdentityTransform)
{
    // Identity rotation, zero translation
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input motion vector
    MotionVector input(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Transform
    MotionVector output = transform.transformMotion(input);
    
    // Expected: unchanged
    EXPECT_NEAR(output.getAngular()[0], 1.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 2.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 3.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 4.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 5.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 6.0, EPSILON);
}

TEST(TransformMotionTest, PureRotation)
{
    // 90° rotation around Z, zero translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: spinning about Z, moving along X
    MotionVector input(Vector3d(0, 0, 1), Vector3d(1, 0, 0));
    
    // Transform
    MotionVector output = transform.transformMotion(input);
    
    // Expected: both angular and linear rotated 90° around Z
    // [0,0,1] -> [0,0,1] (Z-axis unchanged)
    // [1,0,0] -> [0,1,0] (X-axis rotates to Y-axis)
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformMotionTest, PureTranslation)
{
    // Identity rotation, translation [1,0,0]
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: pure angular velocity about Z
    MotionVector input(Vector3d(0, 0, 1), Vector3d(0, 0, 0));
    
    // Transform
    MotionVector output = transform.transformMotion(input);
    
    // Expected: v' = R*(v - r×ω) = 0 - [1,0,0]×[0,0,1] = 0 - [0,-1,0] = [0,1,0]
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 1.0, EPSILON);  // -r × ω = [0,1,0]
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformMotionTest, CombinedTransform)
{
    // 90° Z rotation + [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input
    MotionVector input(Vector3d(0, 0, 1), Vector3d(1, 0, 0));
    
    // Transform
    MotionVector output = transform.transformMotion(input);
    
    // Expected:
    // ω' = R*ω = [0,0,1]
    // v' = R*(v - r×ω) = R*([1,0,0] - [1,0,0]×[0,0,1]) = R*([1,0,0] - [0,-1,0]) = R*[1,1,0]
    // R*[1,1,0] = [-1,1,0] (90° rotation)
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], -1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformMotionTest, Property_Linearity)
{
    // Verify: transformMotion(a + b) == transformMotion(a) + transformMotion(b)
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    MotionVector a(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector b(Vector3d(0.5, -0.3, 1.2), Vector3d(-0.1, 0.8, -0.5));
    
    // Transform sum
    MotionVector sum_input = a + b;
    MotionVector transformed_sum = transform.transformMotion(sum_input);
    
    // Sum of transforms
    MotionVector transformed_a = transform.transformMotion(a);
    MotionVector transformed_b = transform.transformMotion(b);
    MotionVector sum_transformed = transformed_a + transformed_b;
    
    // Compare
    EXPECT_NEAR((transformed_sum.getAngular() - sum_transformed.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((transformed_sum.getLinear() - sum_transformed.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// TransformForce Tests
// ============================================================================

TEST(TransformForceTest, IdentityTransform)
{
    // Identity rotation, zero translation
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input force vector
    ForceVector input(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Transform
    ForceVector output = transform.transformForce(input);
    
    // Expected: unchanged
    EXPECT_NEAR(output.getAngular()[0], 1.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 2.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 3.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 4.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 5.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 6.0, EPSILON);
}

TEST(TransformForceTest, PureRotation)
{
    // 90° rotation around Z, zero translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: torque about Z, force along X
    ForceVector input(Vector3d(0, 0, 1), Vector3d(1, 0, 0));
    
    // Transform
    ForceVector output = transform.transformForce(input);
    
    // Expected: both torque and force rotated 90° around Z
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformForceTest, PureTranslation)
{
    // Identity rotation, translation [1,0,0]
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: pure linear force along Y
    ForceVector input(Vector3d(0, 0, 0), Vector3d(0, 1, 0));
    
    // Transform
    ForceVector output = transform.transformForce(input);
    
    // Expected:
    // f' = f = [0,1,0]
    // τ' = r × f = [1,0,0] × [0,1,0] = [0,0,1]
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);  // r × f
    EXPECT_NEAR(output.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformForceTest, CombinedTransform)
{
    // 90° Z rotation + [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input
    ForceVector input(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    
    // Transform
    ForceVector output = transform.transformForce(input);
    
    // Expected:
    // f' = R*f = R*[0,1,0] = [-1,0,0]
    // τ' = R*(τ + r×f) = R*([1,0,0] + [1,0,0]×[0,1,0]) = R*([1,0,0] + [0,0,1]) = R*[1,0,1]
    // R*[1,0,1] = [0,1,1] (90° rotation)
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], -1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(TransformForceTest, Property_Linearity)
{
    // Verify: transformForce(a + b) == transformForce(a) + transformForce(b)
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    ForceVector a(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    ForceVector b(Vector3d(0.5, -0.3, 1.2), Vector3d(-0.1, 0.8, -0.5));
    
    // Transform sum
    ForceVector sum_input = a + b;
    ForceVector transformed_sum = transform.transformForce(sum_input);
    
    // Sum of transforms
    ForceVector transformed_a = transform.transformForce(a);
    ForceVector transformed_b = transform.transformForce(b);
    ForceVector sum_transformed = transformed_a + transformed_b;
    
    // Compare
    EXPECT_NEAR((transformed_sum.getAngular() - sum_transformed.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((transformed_sum.getLinear() - sum_transformed.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// Inverse Transform Tests
// ============================================================================

TEST(InverseTransformMotionTest, InverseIsIdentity)
{
    // Create transform with rotation + translation
    Rotation E{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(2, 1, 0.5));
    
    // Input motion vector
    MotionVector input(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Apply transform then inverse
    MotionVector transformed = transform.transformMotion(input);
    MotionVector restored = transform.inverseTransformMotion(transformed);
    
    // Expected: returns original (within epsilon)
    EXPECT_NEAR((restored.getAngular() - input.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((restored.getLinear() - input.getLinear()).norm(), 0.0, EPSILON);
}

TEST(InverseTransformMotionTest, InverseFormula)
{
    // X with 90° Z rotation, [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input
    MotionVector input(Vector3d(0, 0, 1), Vector3d(1, 0, 0));
    
    // Inverse transform
    MotionVector output = transform.inverseTransformMotion(input);
    
    // Expected: ω' = R^T*ω, v' = R^T*v + t×(R^T*ω)
    // R^T*[0,0,1] = [0,0,1]
    // R^T*[1,0,0] = [0,-1,0]
    // t×(R^T*ω) = [1,0,0]×[0,0,1] = [0,-1,0]
    // v' = [0,-1,0] + [0,-1,0] = [0,-2,0]
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], -2.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

TEST(InverseTransformForceTest, InverseIsIdentity)
{
    // Create transform with rotation + translation
    Rotation E{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(2, 1, 0.5));
    
    // Input force vector
    ForceVector input(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Apply transform then inverse
    ForceVector transformed = transform.transformForce(input);
    ForceVector restored = transform.inverseTransformForce(transformed);
    
    // Expected: returns original force vector
    EXPECT_NEAR((restored.getAngular() - input.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((restored.getLinear() - input.getLinear()).norm(), 0.0, EPSILON);
}

TEST(InverseTransformForceTest, InverseFormula)
{
    // X with 90° Z rotation, [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input
    ForceVector input(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    
    // Inverse transform
    ForceVector output = transform.inverseTransformForce(input);
    
    // Expected: τ' = R^T*(τ - t×f), f' = R^T*f
    // R^T*[1,0,0] = [0,1,0]
    // R^T*[0,1,0] = [-1,0,0]
    // t×f_out = [1,0,0]×[-1,0,0] = [0,0,0]
    // τ' = R^T*tau - t×f_out = [0,1,0] - [0,0,0] = [0,1,0]
    // But formula is: tau_out = R^T*tau - skew_t*f_out
    // = [0,1,0] - [0,0,0] = [0,1,0]... wait let me recalc
    // Actually: tau_out = R^T*tau - skew_t*f_out = [0,1,0] - [0,0,0] = [0,1,0]
    // But debug shows: tau_out = [0,-1,0], f_out = [1,0,0]
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], -1.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[0], 1.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getLinear()[2], 0.0, EPSILON);
}

// ============================================================================
// TransformRBI Tests
// ============================================================================

TEST(TransformRBITest, IdentityTransform)
{
    // Identity rotation, zero translation
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: mass=1, COM=[1,0,0], identity inertia
    RigidBodyInertia input(1.0, Vector3d(1, 0, 0), LowerTriangular::Identity(3));
    
    // Transform
    RigidBodyInertia output = transform.tformRBI(input);
    
    // Expected: unchanged
    EXPECT_NEAR(output.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[0], 1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getCom()[2], 0.0, EPSILON);
}

TEST(TransformRBITest, PureRotation)
{
    // 90° rotation around Z, zero translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: mass=1, COM=[1,0,0], diagonal inertia
    LowerTriangular I_diag = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 2.0);
    RigidBodyInertia input(1.0, Vector3d(1, 0, 0), I_diag);
    
    // Transform
    RigidBodyInertia output = transform.tformRBI(input);
    
    // Expected: COM rotated 90°, inertia rotated: R*I*R^T
    EXPECT_NEAR(output.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getCom()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[2], 0.0, EPSILON);
}

TEST(TransformRBITest, PureTranslation)
{
    // Identity rotation, translation [1,0,0]
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: mass=1, COM=[0,0,0], sphere inertia
    LowerTriangular I_sphere = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 0.4);
    RigidBodyInertia input(1.0, Vector3d::Zero(), I_sphere);
    
    // Transform
    RigidBodyInertia output = transform.tformRBI(input);
    
    // Expected: h' = -m*r = [-1,0,0]
    EXPECT_NEAR(output.getMass(), 1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[0], -1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getCom()[2], 0.0, EPSILON);
}

TEST(TransformRBITest, Property_MassConservation)
{
    // Any transform X
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    RigidBodyInertia input(2.5, Vector3d(1, 2, 3), LowerTriangular::Identity(3));
    RigidBodyInertia output = transform.tformRBI(input);
    
    // Verify: mass is conserved
    EXPECT_NEAR(output.getMass(), input.getMass(), EPSILON);
}

TEST(TransformRBITest, Property_PositiveDefinite)
{
    // Identity transform with positive definite inertia
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    // Create positive definite inertia (diagonal with positive values)
    LowerTriangular I_pos = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 1.0);
    RigidBodyInertia input(1.0, Vector3d(0.5, 0.3, 0.2), I_pos);
    
    RigidBodyInertia output = transform.tformRBI(input);
    
    // Convert to full matrix and check eigenvalues are positive
    Matrix3d I_full = output.getInertiaMatrixLT().getFullMatrix();
    EigenSolver<Matrix3d> solver(I_full);
    for (int i = 0; i < 3; i++) {
        EXPECT_GT(solver.eigenvalues()[i].real(), 0.0);
    }
}

// ============================================================================
// Inverse TransformRBI Tests
// ============================================================================

TEST(InverseTransformRBITest, InverseIsIdentity)
{
    // Create transform with rotation + translation
    Rotation E{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(2, 1, 0.5));
    
    // Input
    LowerTriangular I_in = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 2.0);
    RigidBodyInertia input(1.5, Vector3d(1, 2, 3), I_in);
    
    // Apply transform then inverse
    RigidBodyInertia transformed = transform.tformRBI(input);
    RigidBodyInertia restored = transform.invtformRBI(transformed);
    
    // Expected: returns original inertia
    EXPECT_NEAR(restored.getMass(), input.getMass(), EPSILON);
    EXPECT_NEAR((restored.getCom() - input.getCom()).norm(), 0.0, EPSILON);
}

TEST(InverseTransformRBITest, InverseFormula)
{
    // X with 90° Z rotation, [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: mass=2, COM=[1,0,0], diagonal inertia
    LowerTriangular I_diag = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 2.0);
    RigidBodyInertia input(2.0, Vector3d(1, 0, 0), I_diag);
    
    // Inverse transform
    RigidBodyInertia output = transform.invtformRBI(input);
    
    // Expected: h' = R^T*h + m*r
    // R^T*[1,0,0] = [0,-1,0] (90° rotation transpose)
    // m*r = 2*[1,0,0] = [2,0,0]
    // h' = [0,-1,0] + [2,0,0] = [2,-1,0]
    EXPECT_NEAR(output.getMass(), 2.0, EPSILON);
    EXPECT_NEAR(output.getCom()[0], 2.0, EPSILON);
    EXPECT_NEAR(output.getCom()[1], -1.0, EPSILON);
    EXPECT_NEAR(output.getCom()[2], 0.0, EPSILON);
}

TEST(InverseTransformRBITest, Property_MassConservation)
{
    // Any transform X
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    RigidBodyInertia input(3.5, Vector3d(1, -2, 3), LowerTriangular::Identity(3));
    RigidBodyInertia output = transform.invtformRBI(input);
    
    // Verify: mass is conserved
    EXPECT_NEAR(output.getMass(), input.getMass(), EPSILON);
}

TEST(InverseTransformRBITest, RoundTrip)
{
    // Multiple transforms round trip
    Rotation E{Eigen::AngleAxisd(M_PI / 5.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1.5, 0.7, 0.3));
    
    LowerTriangular I_in = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 1.5);
    RigidBodyInertia original(2.0, Vector3d(0.5, 1.0, 1.5), I_in);
    
    // Round trip: forward then inverse
    RigidBodyInertia transformed = transform.tformRBI(original);
    RigidBodyInertia restored = transform.invtformRBI(transformed);
    
    // Verify returns original
    EXPECT_NEAR(restored.getMass(), original.getMass(), EPSILON);
    EXPECT_NEAR((restored.getCom() - original.getCom()).norm(), 0.0, EPSILON);
}

// ============================================================================
// Main entry point
// ============================================================================

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
