#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

/**
 * @brief A memory-efficient lower triangular matrix implementation
 * @details Stores only the lower triangular part of the matrix in a packed 1D vector format.
 *          For an n×n matrix, storage requires n(n+1)/2 elements.
 */
class LowerTriangular
{
private:
    Eigen::VectorXd data; ///< Packed storage of lower triangular elements
    int n;                ///< Matrix dimension

    /**
     * @brief Converts 2D coordinates to packed 1D index
     * @param i Row index
     * @param j Column index
     * @return Linear index in packed storage
     */
    inline int getIndex(int i, int j) const noexcept
    {
        return (i * (i + 1)) / 2 + j;
    }

    /**
     * @brief Helper method to scale the matrix
     * @param scalar The scaling factor
     * @return New scaled matrix
     */
    LowerTriangular scale(double scalar) const
    {
        LowerTriangular result(n);
        result.data = data * scalar;
        return result;
    }

public:
    /**
     * @brief Constructor initializing the lower triangular matrix with given values.
     * @param values Array of 6 values representing the lower triangular matrix.
     */
    // LowerTriangular(std::initializer_list<double> values)
    // {
    //     std::copy(values.begin(), values.end(), data);
    // }

    /**
     * @brief Constructs a zero-initialized lower triangular matrix
     * @param size Matrix dimension
     */
    explicit LowerTriangular(int size) : n(size)
    {
        data.resize(size * (size + 1) / 2);
        data.setZero();
    }

    /**
     * @brief Constructs from packed data
     * @param packedData Vector containing packed lower triangular elements
     * @param size Matrix dimension
     * @throws std::invalid_argument if data size doesn't match n(n+1)/2
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
     * @param i Row index
     * @param j Column index
     * @return Matrix element value
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
     * @param i Row index
     * @param j Column index
     * @throws std::invalid_argument if accessing upper triangular part
     * @return Reference to matrix element
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

    /// @brief Returns the packed data vector
    inline const Eigen::VectorXd &getData() const noexcept { return data; }

    /// @brief Returns the matrix dimension
    inline int getSize() const noexcept { return n; }

    /**
     * @brief Computes matrix-vector product
     * @param vec Input vector
     * @return Result vector
     * @throws std::invalid_argument if vector size doesn't match matrix dimension
     */
    Eigen::VectorXd operator*(const Eigen::VectorXd &vec) const
    {
        if (vec.size() != n)
            throw std::invalid_argument("Vector size mismatch");

        Eigen::VectorXd result = Eigen::VectorXd::Zero(n);
#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            result(i) = data.segment(getIndex(i, 0), i + 1).dot(vec.head(i + 1));
        }
        return result;
    }

    /**
     * @brief Computes matrix-matrix product
     * @param other Another lower triangular matrix
     * @return Result matrix
     * @throws std::invalid_argument if matrix dimensions don't match
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
                result(i, j) = data.segment(getIndex(i, 0), i + 1).dot(other.data.segment(getIndex(j, 0), j + 1));
            }
        }
        return result;
    }

    /**
     * @brief Multiplies matrix by a scalar (matrix * scalar)
     * @param scalar The scaling factor
     * @return Scaled matrix
     */
    LowerTriangular operator*(double scalar) const
    {
        return scale(scalar);
    }

    /**
     * @brief Multiplies matrix by a scalar (scalar * matrix)
     * @param scalar The scaling factor
     * @param matrix The matrix to scale
     * @return Scaled matrix
     */
    friend LowerTriangular operator*(double scalar, const LowerTriangular &matrix)
    {
        return matrix.scale(scalar);
    }

    /**
     * @brief Computes matrix transpose
     * @return Transposed matrix
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
     * @brief Adds two lower triangular matrices
     * @param other The matrix to add
     * @return Result of addition
     * @throws std::invalid_argument if matrix dimensions don't match
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
     * @brief Subtracts two lower triangular matrices
     * @param other The matrix to subtract
     * @return Result of subtraction
     * @throws std::invalid_argument if matrix dimensions don't match
     */
    LowerTriangular operator-(const LowerTriangular &other) const
    {
        if (n != other.n)
            throw std::invalid_argument("Matrix size mismatch");

        LowerTriangular result(n);
        result.data = data - other.data;
        return result;
    }

    /**
     * @brief Negates the matrix (unary minus)
     * @return Negated matrix
     */
    LowerTriangular operator-() const
    {
        LowerTriangular result(n);
        result.data = -data;
        return result;
    }

    /**
     * @brief Creates an identity matrix of specified size
     * @param size Dimension of the matrix
     * @return Identity matrix in lower triangular format
     * @throws std::invalid_argument if size is negative
     *
     * Example:
     *   auto I = LowerTriangular::Identity(3); // Creates 3x3 identity matrix
     */
    static LowerTriangular Identity(int size)
    {
        if (size < 0)
            throw std::invalid_argument("Matrix size cannot be negative");

        LowerTriangular result(size);
        for (int i = 0; i < size; ++i)
        {
            result(i, i) = 1.0; // Set diagonal elements to 1
        }
        return result;
    }

    // get full symmetric matrix
    Eigen::MatrixXd getFullMatrix() const
    {
        Eigen::MatrixXd fullMatrix = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                fullMatrix(i, j) = (*this)(i, j);
                fullMatrix(j, i) = (*this)(i, j);
            }
        }
        return fullMatrix;
    }

    /**
     * @brief Stream insertion operator for formatted output
     * @param os Output stream
     * @param matrix Matrix to output
     * @return Reference to output stream
     *
     * Example:
     *   LowerTriangular L(2);
     *   std::cout << L; // Prints full matrix format
     */
    friend std::ostream &operator<<(std::ostream &os, const LowerTriangular &matrix)
    {
        os << "LowerTriangular " << matrix.n << "x" << matrix.n << " matrix:\n";
        for (int i = 0; i < matrix.n; ++i)
        {
            for (int j = 0; j < matrix.n; ++j)
            {
                os << std::setw(10) << std::setprecision(4) << matrix(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    }

    /**
     * @brief Stream extraction operator for reading matrix data
     * @param is Input stream
     * @param matrix Matrix to populate
     * @return Reference to input stream
     * @throws std::runtime_error if input format is invalid
     *
     * Expected input format:
     * size
     * element_1_1
     * element_2_1 element_2_2
     * ...
     */
    friend std::istream &operator>>(std::istream &is, LowerTriangular &matrix)
    {
        int size;
        if (!(is >> size))
        {
            throw std::runtime_error("Failed to read matrix size");
        }

        matrix = LowerTriangular(size);
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                if (!(is >> matrix(i, j)))
                {
                    throw std::runtime_error("Failed to read matrix element at (" +
                                             std::to_string(i) + "," +
                                             std::to_string(j) + ")");
                }
            }
        }
        return is;
    }
};