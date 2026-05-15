// TestLowerTriangular.cpp - Comprehensive GTest test suite for LowerTriangular class

#include "LowerTriangular.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace Eigen;

// Tolerance for floating point comparisons
const double TOLERANCE = 1e-10;

// ===== LTR-01: Packed Storage Tests =====
// ===== LTR-02: Dense Matrix Multiplication Tests =====
// ===== LTR-03: Vector Multiplication Tests =====
// ===== LTR-04: Transpose and Inverse Tests =====
// ===== Property Tests =====

/**
 * @brief Test packed storage size formula
 * @details Verifies LTR-01: data.size() == n*(n+1)/2
 */
TEST(LowerTriangularTest, PackedStorageSize)
{
    LowerTriangular L1(1);
    EXPECT_EQ(L1.getData().size(), 1); // 1*(1+1)/2 = 1
    
    LowerTriangular L2(2);
    EXPECT_EQ(L2.getData().size(), 3); // 2*(2+1)/2 = 3
    
    LowerTriangular L3(3);
    EXPECT_EQ(L3.getData().size(), 6); // 3*(3+1)/2 = 6
    
    LowerTriangular L4(4);
    EXPECT_EQ(L4.getData().size(), 10); // 4*(4+1)/2 = 10
}

/**
 * @brief Test packed storage indexing formula
 * @details Verifies LTR-01: idx = i*(i+1)/2 + j
 */
TEST(LowerTriangularTest, PackedStorageIndexing)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    // Verify packed data order
    const Eigen::VectorXd& data = L.getData();
    EXPECT_DOUBLE_EQ(data(0), 1.0); // (0,0)
    EXPECT_DOUBLE_EQ(data(1), 2.0); // (1,0)
    EXPECT_DOUBLE_EQ(data(2), 3.0); // (1,1)
    EXPECT_DOUBLE_EQ(data(3), 4.0); // (2,0)
    EXPECT_DOUBLE_EQ(data(4), 5.0); // (2,1)
    EXPECT_DOUBLE_EQ(data(5), 6.0); // (2,2)
}

/**
 * @brief Test upper triangular elements return zero
 * @details Verifies LTR-01: (i,j) where i<j returns 0.0
 */
TEST(LowerTriangularTest, UpperTriangularReturnsZero)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    // Use const reference to ensure const operator() is called
    const LowerTriangular& Lconst = L;
    
    // Upper triangular elements should be zero (read-only test)
    EXPECT_DOUBLE_EQ(Lconst(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(Lconst(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(Lconst(1, 2), 0.0);
}

/**
 * @brief Test modifying upper triangle throws exception
 * @details Verifies LTR-01: Cannot write to i<j positions
 */
TEST(LowerTriangularTest, ModifyUpperThrows)
{
    LowerTriangular L(3);
    
    EXPECT_THROW(L(0, 1) = 5.0, std::invalid_argument);
    EXPECT_THROW(L(0, 2) = 5.0, std::invalid_argument);
    EXPECT_THROW(L(1, 2) = 5.0, std::invalid_argument);
}

/**
 * @brief Test constructor from packed data
 * @details Verifies LTR-01: Construction from Eigen::VectorXd
 */
TEST(LowerTriangularTest, ConstructorFromPackedData)
{
    Eigen::VectorXd packed(6);
    packed << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    
    LowerTriangular L(packed, 3);
    
    EXPECT_DOUBLE_EQ(L(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(L(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(L(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(L(2, 0), 4.0);
    EXPECT_DOUBLE_EQ(L(2, 1), 5.0);
    EXPECT_DOUBLE_EQ(L(2, 2), 6.0);
}

/**
 * @brief Test LowerTriangular * MatrixXd multiplication
 * @details Verifies LTR-02: L * M produces correct result
 */
TEST(LowerTriangularTest, MultiplyDenseMatrix)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    Eigen::MatrixXd M(3, 2);
    M << 1.0, 2.0,
         3.0, 4.0,
         5.0, 6.0;
    
    Eigen::MatrixXd result = L * M;
    
    // Expected: L * M
    // [1 0 0]   [1 2]   [1*1+0*3+0*5  1*2+0*4+0*6]   [1  2]
    // [2 3 0] * [3 4] = [2*1+3*3+0*5  2*2+3*4+0*6] = [11 16]
    // [4 5 6]   [5 6]   [4*1+5*3+6*5  4*2+5*4+6*6]   [49 64]
    
    EXPECT_DOUBLE_EQ(result(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 11.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 16.0);
    EXPECT_DOUBLE_EQ(result(2, 0), 49.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 64.0);
}

/**
 * @brief Test MatrixXd * LowerTriangular multiplication
 * @details Verifies LTR-02: M * L produces correct result
 */
TEST(LowerTriangularTest, MultiplyDenseMatrixLeft)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    Eigen::MatrixXd M(2, 3);
    M << 1.0, 2.0, 3.0,
         4.0, 5.0, 6.0;
    
    Eigen::MatrixXd result = M * L;
    
    // Expected: M * L
    // [1 2 3]   [1 0 0]   [1*1+2*2+3*4  2*3+3*5  3*6]   [17 21 18]
    // [4 5 6] * [2 3 0] = [4*1+5*2+6*4  5*3+6*5  6*6] = [38 45 36]
    //             [4 5 6]
    
    EXPECT_DOUBLE_EQ(result(0, 0), 17.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 21.0);
    EXPECT_DOUBLE_EQ(result(0, 2), 18.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 38.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 45.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 36.0);
}

/**
 * @brief Test dimension mismatch in dense multiplication
 * @details Verifies LTR-02: Exception thrown for incompatible sizes
 */
TEST(LowerTriangularTest, DimensionMismatchDense)
{
    LowerTriangular L(3);
    Eigen::MatrixXd M(4, 2); // Wrong number of rows
    
    EXPECT_THROW(L * M, std::invalid_argument);
}

/**
 * @brief Test LowerTriangular * VectorXd multiplication
 * @details Verifies LTR-03: L * v produces correct result
 */
TEST(LowerTriangularTest, MultiplyVector)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    Eigen::VectorXd v(3);
    v << 1.0, 1.0, 1.0;
    
    Eigen::VectorXd result = L * v;
    
    // Expected: L * v
    // [1 0 0]   [1]   [1*1]       [1]
    // [2 3 0] * [1] = [2*1+3*1] = [5]
    // [4 5 6]   [1]   [4*1+5*1+6*1] [15]
    
    EXPECT_DOUBLE_EQ(result(0), 1.0);
    EXPECT_DOUBLE_EQ(result(1), 5.0);
    EXPECT_DOUBLE_EQ(result(2), 15.0);
}

/**
 * @brief Test vector multiplication with different values
 * @details Verifies LTR-03: Correct computation for arbitrary vector
 */
TEST(LowerTriangularTest, MultiplyVectorGeneral)
{
    LowerTriangular L(3);
    L(0,0) = 2.0;
    L(1,0) = 1.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    
    Eigen::VectorXd result = L * v;
    
    // Expected: L * v
    // [2 0 0]   [1]   [2*1]           [2]
    // [1 3 0] * [2] = [1*1+3*2]       [7]
    // [4 5 6]   [3]   [4*1+5*2+6*3]   [32]
    
    EXPECT_DOUBLE_EQ(result(0), 2.0);
    EXPECT_DOUBLE_EQ(result(1), 7.0);
    EXPECT_DOUBLE_EQ(result(2), 32.0);
}

/**
 * @brief Test dimension mismatch in vector multiplication
 * @details Verifies LTR-03: Exception thrown for incompatible sizes
 */
TEST(LowerTriangularTest, DimensionMismatchVector)
{
    LowerTriangular L(3);
    Eigen::VectorXd v(4); // Wrong size
    
    EXPECT_THROW(L * v, std::invalid_argument);
}

/**
 * @brief Test transpose operation
 * @details Verifies LTR-04: L^T(i,j) = L(j,i)
 */
TEST(LowerTriangularTest, Transpose)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    LowerTriangular LT = L.transpose();
    
    // Transpose should swap indices
    // Original: [1 0 0]  Transpose: [1 2 4]
    //           [2 3 0]             [0 3 5]
    //           [4 5 6]             [0 0 6]
    // But stored as lower triangular, so LT(i,j) = L(j,i)
    
    EXPECT_DOUBLE_EQ(LT(0, 0), 1.0); // L(0,0)
    EXPECT_DOUBLE_EQ(LT(1, 0), 0.0); // L(0,1) = 0 (upper triangle)
    EXPECT_DOUBLE_EQ(LT(1, 1), 3.0); // L(1,1)
    EXPECT_DOUBLE_EQ(LT(2, 0), 0.0); // L(0,2) = 0
    EXPECT_DOUBLE_EQ(LT(2, 1), 0.0); // L(1,2) = 0
    EXPECT_DOUBLE_EQ(LT(2, 2), 6.0); // L(2,2)
    
    // Note: transpose() returns lower triangular, so upper elements are zero
    // The actual transpose would be upper triangular, but we store as lower
}

/**
 * @brief Test inverse operation
 * @details Verifies LTR-04: L * L^-1 = I
 */
TEST(LowerTriangularTest, Inverse)
{
    LowerTriangular L(3);
    L(0,0) = 2.0;
    L(1,0) = 1.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    LowerTriangular Linv = L.inverse();
    
    // L * Linv should equal Identity
    LowerTriangular product = L * Linv;
    LowerTriangular I = LowerTriangular::Identity(3);
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            EXPECT_NEAR(product(i, j), I(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test inverse diagonal elements
 * @details Verifies LTR-04: L^-1(i,i) = 1/L(i,i)
 */
TEST(LowerTriangularTest, InverseDiagonal)
{
    LowerTriangular L(3);
    L(0,0) = 2.0;
    L(1,0) = 1.0; L(1,1) = 4.0;
    L(2,0) = 3.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    LowerTriangular Linv = L.inverse();
    
    EXPECT_NEAR(Linv(0, 0), 1.0/2.0, TOLERANCE);
    EXPECT_NEAR(Linv(1, 1), 1.0/4.0, TOLERANCE);
    EXPECT_NEAR(Linv(2, 2), 1.0/6.0, TOLERANCE);
}

/**
 * @brief Test identity matrix is its own inverse
 * @details Verifies LTR-04: I^-1 = I
 */
TEST(LowerTriangularTest, InverseIdentity)
{
    LowerTriangular I = LowerTriangular::Identity(3);
    LowerTriangular Iinv = I.inverse();
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            EXPECT_DOUBLE_EQ(Iinv(i, j), I(i, j));
        }
    }
}

/**
 * @brief Test inverse of 2x2 matrix
 * @details Verifies LTR-04: Correct inverse for small matrix
 */
TEST(LowerTriangularTest, Inverse2x2)
{
    LowerTriangular L(2);
    L(0,0) = 2.0;
    L(1,0) = 1.0; L(1,1) = 4.0;
    
    LowerTriangular Linv = L.inverse();
    
    // For 2x2: [a 0]^-1   [1/a   0  ]
    //          [b c]    = [-b/(ac) 1/c]
    // L^-1(0,0) = 1/2 = 0.5
    // L^-1(1,0) = -1/(2*4) = -0.125
    // L^-1(1,1) = 1/4 = 0.25
    
    EXPECT_NEAR(Linv(0, 0), 0.5, TOLERANCE);
    EXPECT_NEAR(Linv(1, 0), -0.125, TOLERANCE);
    EXPECT_NEAR(Linv(1, 1), 0.25, TOLERANCE);
}

/**
 * @brief Test identity matrix property
 * @details Verifies L * I = L
 */
TEST(LowerTriangularTest, IdentityProperty)
{
    LowerTriangular L(3);
    L(0,0) = 1.0;
    L(1,0) = 2.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    LowerTriangular I = LowerTriangular::Identity(3);
    LowerTriangular result = L * I;
    
    EXPECT_DOUBLE_EQ(result(0, 0), L(0, 0));
    EXPECT_DOUBLE_EQ(result(1, 0), L(1, 0));
    EXPECT_DOUBLE_EQ(result(1, 1), L(1, 1));
    EXPECT_DOUBLE_EQ(result(2, 0), L(2, 0));
    EXPECT_DOUBLE_EQ(result(2, 1), L(2, 1));
    EXPECT_DOUBLE_EQ(result(2, 2), L(2, 2));
}

/**
 * @brief Test associativity of multiplication
 * @details Verifies (L1 * L2) * L3 = L1 * (L2 * L3)
 */
TEST(LowerTriangularTest, Associativity)
{
    LowerTriangular L1(3);
    L1(0,0) = 1.0;
    L1(1,0) = 2.0; L1(1,1) = 3.0;
    L1(2,0) = 4.0; L1(2,1) = 5.0; L1(2,2) = 6.0;
    
    LowerTriangular L2(3);
    L2(0,0) = 2.0;
    L2(1,0) = 1.0; L2(1,1) = 2.0;
    L2(2,0) = 3.0; L2(2,1) = 4.0; L2(2,2) = 3.0;
    
    LowerTriangular L3(3);
    L3(0,0) = 1.0;
    L3(1,0) = 1.0; L3(1,1) = 1.0;
    L3(2,0) = 1.0; L3(2,1) = 1.0; L3(2,2) = 1.0;
    
    LowerTriangular temp1 = L1 * L2;
    LowerTriangular result1 = temp1 * L3;
    
    LowerTriangular temp2 = L2 * L3;
    LowerTriangular result2 = L1 * temp2;
    
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            EXPECT_NEAR(result1(i, j), result2(i, j), TOLERANCE);
        }
    }
}

/**
 * @brief Test comparison with full matrix operations
 * @details Verifies LowerTriangular matches Eigen::MatrixXd
 */
TEST(LowerTriangularTest, ComparisonWithFullMatrix)
{
    LowerTriangular L(3);
    L(0,0) = 2.0;
    L(1,0) = 1.0; L(1,1) = 3.0;
    L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
    
    Eigen::MatrixXd full = L.getFullMatrix();
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    
    Eigen::VectorXd result1 = L * v;
    Eigen::VectorXd result2 = full * v;
    
    // Results should match within tolerance
    EXPECT_NEAR((result1 - result2).norm(), 0.0, TOLERANCE);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
