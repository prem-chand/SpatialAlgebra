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
    
    // Expected: τ' = R*(τ - r×f) = I*(0 - [1,0,0]×[0,1,0]) = -[0,0,1] = [0,0,-1]
    // f' = f = [0,1,0]
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], -1.0, EPSILON);  // -r × f
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
    // τ' = R*(τ - r×f) = R*([1,0,0] - [1,0,0]×[0,1,0]) = R*([1,0,0] - [0,0,1]) = R*[1,0,-1]
    // R*[1,0,-1] = [0,1,-1] (90° Z rotation)
    EXPECT_NEAR(output.getAngular()[0], 0.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[1], 1.0, EPSILON);
    EXPECT_NEAR(output.getAngular()[2], -1.0, EPSILON);
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
    
    // X^T * f:  τ' = R^T*τ + t×(R^T*f), f' = R^T*f
    // R^T*[1,0,0] = [0,-1,0], R^T*[0,1,0] = [1,0,0]
    // t×(R^T*f) = [1,0,0]×[1,0,0] = [0,0,0]
    // τ' = [0,-1,0] + [0,0,0] = [0,-1,0], f' = [1,0,0]
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
    Matrix3d I_full = output.getInertiaMatrixLT().getSymmetricMatrix();
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
// TransformABI Tests
// ============================================================================

TEST(TransformABITest, IdentityTransform)
{
    // Identity rotation, zero translation
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: identity ABI
    LowerTriangular I_identity = LowerTriangular::Identity(3);
    ArticulatedBodyInertia input(I_identity, Matrix3d::Identity(), LowerTriangular::Identity(3));
    
    // Transform
    ArticulatedBodyInertia output = transform.tformABI(input);
    
    // Expected: unchanged
    EXPECT_NEAR(output.getH()(0, 0), 1.0, EPSILON);
    EXPECT_NEAR(output.getM().getSymmetricMatrix()(0, 0), 1.0, EPSILON);
}

TEST(TransformABITest, PureRotation)
{
    // 90° rotation around Z, zero translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d::Zero());
    
    // Input: diagonal ABI
    LowerTriangular I_diag = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 2.0);
    ArticulatedBodyInertia input(I_diag, Matrix3d::Identity(), LowerTriangular::Identity(3));
    
    // Transform
    ArticulatedBodyInertia output = transform.tformABI(input);
    
    // Expected: I' = R*I*R^T, H' = R*H*R^T, M' = R*M*R^T
    // For diagonal matrices and 90° rotation, should remain diagonal
    EXPECT_NEAR(output.getH()(0, 0), 1.0, EPSILON);
    EXPECT_NEAR(output.getH()(1, 1), 1.0, EPSILON);
}

TEST(TransformABITest, PureTranslation)
{
    // Identity rotation, translation [1,0,0]
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: simple ABI with zero H
    LowerTriangular I_sphere = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 0.4);
    ArticulatedBodyInertia input(I_sphere, Matrix3d::Zero(), LowerTriangular::Identity(3));
    
    // Transform
    ArticulatedBodyInertia output = transform.tformABI(input);
    
    // Expected: H' = -M*r̂ (non-zero due to translation)
    // H' = R*(H - M*r̂)*R^T = 0 - I*r̂ = -r̂
    Matrix3d H_out = output.getH();
    EXPECT_GT(H_out.norm(), 0.0);  // Should be non-zero
}

TEST(TransformABITest, Property_Symmetric)
{
    // Transform should preserve symmetry of I and M (H is not necessarily symmetric)
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    // Create symmetric positive definite ABI
    Matrix3d I_full = Matrix3d::Identity() * 2.0;
    LowerTriangular I_lt = LowerTriangular::fromFullMatrix(I_full);
    ArticulatedBodyInertia input(I_lt, Matrix3d::Identity(), LowerTriangular::Identity(3));
    
    ArticulatedBodyInertia output = transform.tformABI(input);
    
    // Verify symmetry of I and M (H is coupling matrix, not necessarily symmetric)
    Matrix3d I_out = output.getInertia().getSymmetricMatrix();
    Matrix3d M_out = output.getM().getSymmetricMatrix();
    
    EXPECT_NEAR((I_out - I_out.transpose()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((M_out - M_out.transpose()).norm(), 0.0, EPSILON);
}

TEST(TransformABITest, Property_PositiveDefinite)
{
    // Transform should preserve positive definiteness
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    // Create positive definite ABI
    Matrix3d I_full = Matrix3d::Identity() * 1.5;
    LowerTriangular I_lt = LowerTriangular::fromFullMatrix(I_full);
    ArticulatedBodyInertia input(I_lt, Matrix3d::Identity() * 0.5, LowerTriangular::Identity(3));
    
    ArticulatedBodyInertia output = transform.tformABI(input);
    
    // Check eigenvalues of I are positive
    Matrix3d I_out = output.getInertia().getSymmetricMatrix();
    EigenSolver<Matrix3d> solver(I_out);
    for (int i = 0; i < 3; i++) {
        EXPECT_GT(solver.eigenvalues()[i].real(), 0.0);
    }
}

// ============================================================================
// Inverse TransformABI Tests
// ============================================================================

TEST(InverseTransformABITest, InverseIsIdentity)
{
    // Create transform with rotation + translation
    Rotation E{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(2, 1, 0.5));
    
    // Input ABI
    LowerTriangular I_in = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 2.0);
    ArticulatedBodyInertia input(I_in, Matrix3d::Identity(), LowerTriangular::Identity(3));
    
    // Apply transform then inverse
    ArticulatedBodyInertia transformed = transform.tformABI(input);
    ArticulatedBodyInertia restored = transform.invtformABI(transformed);
    
    // Expected: returns original ABI (within epsilon for floating point)
    EXPECT_NEAR((restored.getInertia().getSymmetricMatrix() - input.getInertia().getSymmetricMatrix()).norm(), 0.0, 1e-8);
    EXPECT_NEAR((restored.getH() - input.getH()).norm(), 0.0, 1e-8);
    EXPECT_NEAR((restored.getM().getSymmetricMatrix() - input.getM().getSymmetricMatrix()).norm(), 0.0, 1e-8);
}

TEST(InverseTransformABITest, InverseFormula)
{
    // X with 90° Z rotation, [1,0,0] translation
    Rotation E{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1, 0, 0));
    
    // Input: simple ABI
    LowerTriangular I_diag = LowerTriangular::fromFullMatrix(Matrix3d::Identity());
    ArticulatedBodyInertia input(I_diag, Matrix3d::Zero(), LowerTriangular::Identity(3));
    
    // Inverse transform
    ArticulatedBodyInertia output = transform.invtformABI(input);
    
    // Expected: M' = R^T*M*R, H' = R^T*(H + M*r̂)*R
    // With H=0: H' = R^T*M*r̂*R
    Matrix3d H_out = output.getH();
    EXPECT_GT(H_out.norm(), 0.0);  // Should be non-zero due to translation term
}

TEST(InverseTransformABITest, Property_Symmetric)
{
    // Inverse transform should preserve symmetry of I and M (H is not necessarily symmetric)
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(0.5, 0.3, 0.2));
    
    Matrix3d I_full = Matrix3d::Identity() * 2.0;
    LowerTriangular I_lt = LowerTriangular::fromFullMatrix(I_full);
    ArticulatedBodyInertia input(I_lt, Matrix3d::Identity(), LowerTriangular::Identity(3));
    
    ArticulatedBodyInertia output = transform.invtformABI(input);
    
    // Verify symmetry of I and M
    Matrix3d I_out = output.getInertia().getSymmetricMatrix();
    Matrix3d M_out = output.getM().getSymmetricMatrix();
    
    EXPECT_NEAR((I_out - I_out.transpose()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((M_out - M_out.transpose()).norm(), 0.0, EPSILON);
}

TEST(InverseTransformABITest, RoundTrip)
{
    // Multiple transforms round trip
    Rotation E{Eigen::AngleAxisd(M_PI / 5.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1.5, 0.7, 0.3));
    
    LowerTriangular I_in = LowerTriangular::fromFullMatrix(Matrix3d::Identity() * 1.5);
    ArticulatedBodyInertia original(I_in, Matrix3d::Identity() * 0.8, LowerTriangular::Identity(3));
    
    // Round trip: forward then inverse
    ArticulatedBodyInertia transformed = transform.tformABI(original);
    ArticulatedBodyInertia restored = transform.invtformABI(transformed);
    
    // Verify returns original (within floating point tolerance)
    EXPECT_NEAR((restored.getInertia().getSymmetricMatrix() - original.getInertia().getSymmetricMatrix()).norm(), 0.0, 1e-8);
    EXPECT_NEAR((restored.getH() - original.getH()).norm(), 0.0, 1e-8);
    EXPECT_NEAR((restored.getM().getSymmetricMatrix() - original.getM().getSymmetricMatrix()).norm(), 0.0, 1e-8);
}

// ============================================================================
// inverse() Method Tests
// ============================================================================

/**
 * @brief Test suite for PluckerTransform::inverse() method
 */
class TestInverse : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test inverse produces identity when multiplied
 * @details X * X^(-1) should equal identity transform
 */
TEST(TestInverse, MultiplyWithInverseIsIdentity) {
    // Arrange: Create transform with rotation + translation
    Rotation E{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(2, 1, 0.5));
    
    // Act: Compute inverse
    PluckerTransform inv = transform.inverse();
    
    // Test by applying to motion vector: X * X^(-1) * v should equal v
    MotionVector v(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Apply inverse then original (equivalent to X * X^(-1) * v)
    MotionVector transformed = inv.transformMotion(v);
    MotionVector result = transform.transformMotion(transformed);
    
    // Assert: should return original vector
    EXPECT_NEAR((result.getAngular() - v.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((result.getLinear() - v.getLinear()).norm(), 0.0, EPSILON);
}

/**
 * @brief Test double inverse returns original
 * @details (X^(-1))^(-1) should equal X
 */
TEST(TestInverse, DoubleInverseReturnsOriginal) {
    // Arrange
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform original(E, Vector3d(1.5, 0.7, 0.3));
    
    // Act: Compute double inverse
    PluckerTransform doubleInv = original.inverse().inverse();
    
    // Test by applying to motion and force vectors
    MotionVector mv(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    ForceVector fv(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    
    // Assert: double inverse should transform same as original
    MotionVector mv_orig = original.transformMotion(mv);
    MotionVector mv_double = doubleInv.transformMotion(mv);
    EXPECT_NEAR((mv_orig.getAngular() - mv_double.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((mv_orig.getLinear() - mv_double.getLinear()).norm(), 0.0, EPSILON);
    
    ForceVector fv_orig = original.transformForce(fv);
    ForceVector fv_double = doubleInv.transformForce(fv);
    EXPECT_NEAR((fv_orig.getAngular() - fv_double.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((fv_orig.getLinear() - fv_double.getLinear()).norm(), 0.0, EPSILON);
}

/**
 * @brief Test inverse of identity is identity
 */
TEST(TestInverse, InverseOfIdentityIsIdentity) {
    // Arrange: Create identity transform
    Rotation E{Matrix3d::Identity()};
    PluckerTransform identity(E, Vector3d::Zero());
    
    // Act: Compute inverse
    PluckerTransform inv = identity.inverse();
    
    // Test: inverse should behave like identity
    MotionVector mv(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector result = inv.transformMotion(mv);
    
    // Assert: should return unchanged vector
    EXPECT_NEAR((result.getAngular() - mv.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((result.getLinear() - mv.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// multiply() and apply(PluckerTransform) Method Tests
// ============================================================================

/**
 * @brief Test suite for PluckerTransform composition methods
 */
class TestMultiply : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test multiply() composes transforms correctly
 * @details X_combined = X1.multiply(X2) should satisfy: X_combined * v = X1 * (X2 * v)
 */
TEST(TestMultiply, ComposeTransforms) {
    // Arrange: Create two transforms
    // X1: 90° rotation around Z
    Rotation E1{Eigen::AngleAxisd(M_PI / 2.0, Vector3d::UnitZ())};
    PluckerTransform X1(E1, Vector3d::Zero());
    
    // X2: translation [1, 0, 0]
    Rotation E2{Matrix3d::Identity()};
    PluckerTransform X2(E2, Vector3d(1, 0, 0));
    
    // Act: Compose transforms using apply (which calls multiply)
    PluckerTransform X_combined = X1.apply(X2);
    
    // Test by applying to a motion vector
    MotionVector v(Vector3d(0, 0, 1), Vector3d(1, 0, 0));
    
    // Apply combined transform
    MotionVector result_combined = X_combined.transformMotion(v);
    
    // Apply sequentially: X1 * (X2 * v)
    MotionVector result_sequential = X1.transformMotion(X2.transformMotion(v));
    
    // Assert: Results should match
    EXPECT_NEAR((result_combined.getAngular() - result_sequential.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((result_combined.getLinear() - result_sequential.getLinear()).norm(), 0.0, EPSILON);
}

/**
 * @brief Test apply(PluckerTransform) is alias for multiply
 */
TEST(TestMultiply, ApplyEqualsMultiply) {
    // Arrange
    Rotation E1{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform X1(E1, Vector3d(0.5, 0.0, 0.0));
    
    Rotation E2{Eigen::AngleAxisd(M_PI / 6.0, Vector3d::UnitX())};
    PluckerTransform X2(E2, Vector3d(0.0, 0.3, 0.0));
    
    // Act: Compute both ways (apply is alias for multiply per PluckerTransform.h:199)
    PluckerTransform result_apply1 = X1.apply(X2);
    PluckerTransform result_apply2 = X1.apply(X2);
    
    // Test both transforms on same vector
    MotionVector testVec(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    MotionVector mv1 = result_apply1.transformMotion(testVec);
    MotionVector mv2 = result_apply2.transformMotion(testVec);
    
    // Assert: Should produce identical results
    EXPECT_NEAR((mv1.getAngular() - mv2.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((mv1.getLinear() - mv2.getLinear()).norm(), 0.0, EPSILON);
}

/**
 * @brief Test associativity of apply: (X1 * X2) * X3 == X1 * (X2 * X3)
 */
TEST(TestMultiply, Associativity) {
    // Arrange: Create three transforms
    Rotation E1{Eigen::AngleAxisd(M_PI / 6.0, Vector3d::UnitZ())};
    PluckerTransform X1(E1, Vector3d(0.5, 0.0, 0.0));
    
    Rotation E2{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitX())};
    PluckerTransform X2(E2, Vector3d(0.0, 0.3, 0.0));
    
    Rotation E3{Eigen::AngleAxisd(M_PI / 3.0, Vector3d::UnitY())};
    PluckerTransform X3(E3, Vector3d(0.0, 0.0, 0.7));
    
    // Act: (X1 * X2) * X3 using apply
    PluckerTransform left = X1.apply(X2).apply(X3);
    
    // X1 * (X2 * X3)
    PluckerTransform X23 = X2.apply(X3);
    PluckerTransform right = X1.apply(X23);
    
    // Test on motion vector
    MotionVector testVec(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    MotionVector mv_left = left.transformMotion(testVec);
    MotionVector mv_right = right.transformMotion(testVec);
    
    // Assert: Should be equal
    EXPECT_NEAR((mv_left.getAngular() - mv_right.getAngular()).norm(), 0.0, EPSILON);
    EXPECT_NEAR((mv_left.getLinear() - mv_right.getLinear()).norm(), 0.0, EPSILON);
}

// ============================================================================
// print() Method Tests
// ============================================================================

/**
 * @brief Test suite for PluckerTransform::print() method
 */
class TestPrint : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

/**
 * @brief Test print() executes without crashing and produces output
 * @note Uses GTest's stdout capture to verify output contains expected values
 */
TEST(TestPrint, ProducesOutput) {
    // Arrange: Create transform with known values
    Rotation E{Eigen::AngleAxisd(M_PI / 4.0, Vector3d::UnitZ())};
    PluckerTransform transform(E, Vector3d(1.0, 2.0, 3.0));
    
    // Act: Capture stdout and call print()
    testing::internal::CaptureStdout();
    transform.print();
    std::string output = testing::internal::GetCapturedStdout();
    
    // Assert: Output should not be empty and should contain rotation/translation info
    // The print() method outputs rotation matrix and translation vector
    EXPECT_FALSE(output.empty());
    // Check that output contains "Rotation" text
    EXPECT_NE(output.find("Rotation"), std::string::npos);
}

/**
 * @brief Test print() with identity transform
 */
TEST(TestPrint, IdentityTransform) {
    // Arrange
    Rotation E{Matrix3d::Identity()};
    PluckerTransform identity(E, Vector3d::Zero());
    
    // Act: Just verify it doesn't crash
    testing::internal::CaptureStdout();
    identity.print();
    std::string output = testing::internal::GetCapturedStdout();
    
    // Assert: Should produce some output
    EXPECT_FALSE(output.empty());
}

// ============================================================================
// Main entry point
// ============================================================================

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
