#pragma once

/**
 * @file LowerTriangular.h
 * @brief Memory-efficient implementation of lower triangular matrices
 * @details This file provides a specialized implementation of lower triangular matrices
 *          optimized for both memory usage and computational efficiency. It uses a
 *          packed storage format that only stores the non-zero elements, making it
 *          particularly suitable for applications in spatial algebra and robotics.
 *
 *          Storage Format:
 *          For an n×n matrix, elements are stored in a 1D array in row-major order:
 *          [ a11,
 *            a21, a22,
 *            a31, a32, a33,
 *            ...
 *          ]
 *          Total storage: n(n+1)/2 elements
 *
 *          Performance Features:
 *          - O(n²/2) storage instead of O(n²)
 *          - Implicit zeros in upper triangle
 *          - OpenMP parallelization for matrix operations
 *          - Efficient matrix-vector and matrix-matrix products
 *
 * Example usage:
 * @code{.cpp}
 *     // Create a 3x3 lower triangular matrix
 *     LowerTriangular L(3);
 *
 *     // Set elements (only lower triangle)
 *     L(0,0) = 1.0;
 *     L(1,0) = 2.0; L(1,1) = 3.0;
 *     L(2,0) = 4.0; L(2,1) = 5.0; L(2,2) = 6.0;
 *
 *     // Matrix operations
 *     LowerTriangular L2(3);
 *     auto product = L * L2;        // Matrix-matrix multiplication
 *     auto scaled = L * 2.0;        // Scalar multiplication
 *     auto result = L * someVector; // Matrix-vector multiplication
 * @endcode
 *
 * @note This implementation is particularly optimized for the needs of spatial
 *       algebra computations, where lower triangular matrices frequently appear
 *       in inertia tensors and articulated body calculations.
 */

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

/**
 * @brief Memory-efficient lower triangular matrix class
 * @details Implements a lower triangular matrix using packed storage format,
 *          providing both memory efficiency and fast operations. The class
 *          ensures that only the lower triangular elements are stored and
 *          modified, while upper triangular elements are implicitly zero.
 *
 *          Mathematical Properties:
 *          - All elements above main diagonal are zero
 *          - Matrix-matrix product of lower triangular matrices is lower triangular
 *          - Transpose converts lower to upper triangular
 *
 *          Implementation Features:
 *          - Packed storage format for memory efficiency
 *          - Bounds checking in debug mode
 *          - Exception safety for invalid operations
 *          - OpenMP parallelization for performance
 */
class LowerTriangular
{
private:
    Eigen::VectorXd data; ///< Packed storage of lower triangular elements
    int n;                ///< Matrix dimension

    /**
     * @brief Maps 2D matrix indices to 1D packed storage index
     * @param i Row index (0-based)
     * @param j Column index (0-based)
     * @return Linear index in packed storage
     * @details Computes the index in the packed storage array for element (i,j).
     *          The mapping follows the formula: idx = i*(i+1)/2 + j
     *          This assumes row-major storage of lower triangular elements.
     */
    inline int getIndex(int i, int j) const noexcept
    {
        return (i * (i + 1)) / 2 + j;
    }

    /**
     * @brief Internal helper for scalar multiplication
     * @param scalar Scaling factor
     * @return New scaled matrix
     * @details Efficiently implements scaling by operating directly on packed data
     */
    LowerTriangular scale(double scalar) const
    {
        LowerTriangular result(n);
        result.data = data * scalar;
        return result;
    }

public:
    /**
     * @brief Constructs a zero-initialized lower triangular matrix
     * @param size Matrix dimension (n×n)
     * @details Allocates memory for n(n+1)/2 elements and initializes them to zero.
     *          The size parameter determines both rows and columns of the matrix.
     */
    explicit LowerTriangular(int size) : n(size)
    {
        data.resize(size * (size + 1) / 2);
        data.setZero();
    }

    /**
     * @brief Constructs from packed data array
     * @param packedData Vector containing packed lower triangular elements
     * @param size Matrix dimension
     * @throws std::invalid_argument if data size doesn't match n(n+1)/2
     * @details Creates a matrix using pre-packed data. The data must be in
     *          row-major order and contain exactly n(n+1)/2 elements.
     */
    LowerTriangular(const Eigen::VectorXd &packedData, int size)
        : data(packedData), n(size)
    {
        if (packedData.size() != size * (size + 1) / 2)
        {
            throw std::invalid_argument("Invalid data size for lower triangular matrix");
        }
    }

    /**
     * @brief Matrix element access (const)
     * @param i Row index (0-based)
     * @param j Column index (0-based)
     * @return Matrix element value
     * @throws std::out_of_range if indices are out of bounds (in debug mode)
     * @details Returns 0 for upper triangular elements (i < j)
     */
    inline double operator()(int i, int j) const
    {
#ifndef NDEBUG
        if (i >= n || j >= n || i < 0 || j < 0)
            throw std::out_of_range("Index out of bounds");
#endif
        return (i < j) ? 0.0 : data[getIndex(i, j)];
    }

    /**
     * @brief Matrix element access (non-const)
     * @param i Row index (0-based)
     * @param j Column index (0-based)
     * @return Reference to matrix element
     * @throws std::out_of_range if indices are out of bounds (in debug mode)
     * @throws std::invalid_argument if attempting to modify upper triangle
     * @details Provides write access to lower triangular elements only
     */
    inline double &operator()(int i, int j)
    {
#ifndef NDEBUG
        if (i >= n || j >= n || i < 0 || j < 0)
            throw std::out_of_range("Index out of bounds");
#endif
        if (i < j)
            throw std::invalid_argument("Cannot modify upper triangular elements");
        return data[getIndex(i, j)];
    }

    /**
     * @brief Access to packed storage data
     * @return Const reference to internal packed data vector
     */
    inline const Eigen::VectorXd &getData() const noexcept { return data; }

    /**
     * @brief Get matrix dimension
     * @return Size of the square matrix
     */
    inline int getSize() const noexcept { return n; }

    /**
     * @brief Matrix-matrix multiplication (this * other)
     * @param other Right-hand side lower triangular matrix
     * @return Result of multiplication
     * @throws std::invalid_argument if matrix dimensions don't match
     * @details Implements optimized multiplication for lower triangular matrices.
     *          The result is guaranteed to be lower triangular.
     *          Uses OpenMP for parallel computation when available.
     */
    LowerTriangular operator*(const LowerTriangular &other) const
    {
        if (n != other.n)
            throw std::invalid_argument("Matrix size mismatch");

        LowerTriangular result(n);
#pragma omp parallel for collapse(2)
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                double sum = 0.0;
                for (int k = j; k <= i; ++k)
                {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    /**
     * @brief Matrix-matrix multiplication with dense matrix (L * M)
     * @param other Dense matrix to multiply with
     * @return Result as dense matrix
     * @throws std::invalid_argument if dimensions don't match
     * @details Efficiently implements multiplication with a dense matrix,
     *          taking advantage of the lower triangular structure.
     */
    Eigen::MatrixXd operator*(const Eigen::MatrixXd &other) const
    {
        if (n != other.rows())
            throw std::invalid_argument("Matrix dimension mismatch");

        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, other.cols());
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                result.row(i) += (*this)(i, j) * other.row(j);
            }
        }
        return result;
    }

    /**
     * @brief Dense matrix multiplication with lower triangular (M * L)
     * @param lhs Dense matrix on the left
     * @param rhs Lower triangular matrix on the right
     * @return Result as dense matrix
     * @throws std::invalid_argument if dimensions don't match
     */
    friend Eigen::MatrixXd operator*(const Eigen::MatrixXd &lhs, const LowerTriangular &rhs)
    {
        if (lhs.cols() != rhs.n)
            throw std::invalid_argument("Matrix dimension mismatch");

        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(lhs.rows(), rhs.n);
        for (int j = 0; j < rhs.n; ++j)
        {
            for (int k = 0; k <= j; ++k)
            {
                result.col(j) += rhs(j, k) * lhs.col(k);
            }
        }
        return result;
    }

    /**
     * @brief Scalar multiplication (matrix * scalar)
     * @param scalar Scaling factor
     * @return Scaled matrix
     */
    LowerTriangular operator*(double scalar) const
    {
        return scale(scalar);
    }

    /**
     * @brief Scalar multiplication (scalar * matrix)
     * @param scalar Scaling factor
     * @param matrix Matrix to scale
     * @return Scaled matrix
     */
    friend LowerTriangular operator*(double scalar, const LowerTriangular &matrix)
    {
        return matrix.scale(scalar);
    }

    /**
     * @brief Matrix addition
     * @param other Matrix to add
     * @return Sum of matrices
     * @throws std::invalid_argument if dimensions don't match
     */
    LowerTriangular operator+(const LowerTriangular &other) const
    {
        if (n != other.n)
            throw std::invalid_argument("Matrix size mismatch");

        LowerTriangular result(n);
        result.data = data + other.data;
        return result;
    }

    /**
     * @brief Matrix transpose
     * @return Transposed matrix
     * @details Creates a new lower triangular matrix representing
     *          the transpose of this matrix.
     */
    LowerTriangular transpose() const
    {
        LowerTriangular result(n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                result(i, j) = (*this)(j, i);
            }
        }
        return result;
    }

    /**
     * @brief Create identity matrix
     * @param size Matrix dimension
     * @return Identity matrix of given size
     * @details Static factory method to create an identity matrix
     */
    static LowerTriangular Identity(int size)
    {
        LowerTriangular result(size);
        for (int i = 0; i < size; ++i)
        {
            result(i, i) = 1.0;
        }
        return result;
    }

    /**
     * @brief Stream output operator
     * @param os Output stream
     * @param matrix Matrix to output
     * @return Modified output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const LowerTriangular &matrix)
    {
        for (int i = 0; i < matrix.n; ++i)
        {
            for (int j = 0; j < matrix.n; ++j)
            {
                os << std::setw(8) << matrix(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    }

    /**
     * @brief Convert to full matrix representation
     * @return Dense matrix containing all elements
     * @details Creates a full n×n Eigen::MatrixXd from the packed storage format.
     *          The resulting matrix contains zeros in the upper triangular part
     *          and the stored values in the lower triangular part.
     *
     * Example usage:
     * @code{.cpp}
     *     LowerTriangular L(3);
     *     L(1,0) = 2.0; L(2,0) = 3.0; L(2,1) = 4.0;
     *     Eigen::MatrixXd full = L.getFullMatrix();
     *     // full = [0 0 0]
     *     //        [2 0 0]
     *     //        [3 4 0]
     * @endcode
     */
    Eigen::MatrixXd getFullMatrix() const
    {
        Eigen::MatrixXd fullMatrix = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                fullMatrix(i, j) = (*this)(i, j);
            }
        }
        return fullMatrix;
    }

    /**
     * @brief Create lower triangular matrix from full matrix
     * @param fullMatrix Dense matrix to convert
     * @return Lower triangular matrix containing lower elements
     * @throws std::invalid_argument if input matrix is not square
     * @details Creates a LowerTriangular from a dense matrix by extracting
     *          the lower triangular elements. The upper triangular part
     *          of the input matrix is ignored.
     *
     * Example usage:
     * @code{.cpp}
     *     Eigen::MatrixXd M(2,2);
     *     M << 1, 2,
     *          3, 4;
     *     LowerTriangular L = LowerTriangular::fromFullMatrix(M);
     *     // L contains: [1 0]
     *     //            [3 4]
     * @endcode
     *
     * @note This method only considers the lower triangular part of the input
     *       matrix. Any values in the upper triangular part are ignored.
     */
    static LowerTriangular fromFullMatrix(const Eigen::MatrixXd &fullMatrix)
    {
        if (fullMatrix.rows() != fullMatrix.cols())
        {
            throw std::invalid_argument("Input matrix must be square");
        }

        int size = fullMatrix.rows();
        LowerTriangular result(size);

        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                result(i, j) = fullMatrix(i, j);
            }
        }
        return result;
    }
};