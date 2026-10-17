#ifndef LOWER_TRIANGULAR_H
#define LOWER_TRIANGULAR_H

/**
 * @file LowerTriangular.h
 * @brief Fixed-size 3×3 lower-triangle storage for symmetric matrices
 * @details Stores the six canonical lower-triangle values in a fixed-size Eigen
 *          matrix. Symmetric operations reconstruct the matching upper triangle.
 */

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include "SpatialVector.h"

namespace SpatialAlgebra {

class LowerTriangular
{
private:
    Eigen::Matrix3d data; ///< Lower triangle is canonical; upper triangle is unused
    mutable Eigen::Matrix3d symmetricData;

public:
    /** @brief Zero-initialized matrix (always 3×3, size ignored for backward compat) */
    LowerTriangular() : data(Eigen::Matrix3d::Zero()), symmetricData(Eigen::Matrix3d::Zero()) {}

    /** @brief Legacy constructor: size is ignored (always 3×3) */
    explicit LowerTriangular(int /*size*/) : data(Eigen::Matrix3d::Zero()), symmetricData(Eigen::Matrix3d::Zero()) {}

    /** @brief Legacy packed-data constructor (size must be 3, data must have 6 elements) */
    LowerTriangular(const Eigen::VectorXd& packed, int size) : data(Eigen::Matrix3d::Zero()), symmetricData(Eigen::Matrix3d::Zero()) {
        if (packed.size() >= 6) {
            data(0,0) = packed(0);
            data(1,0) = packed(1); data(1,1) = packed(2);
            data(2,0) = packed(3); data(2,1) = packed(4); data(2,2) = packed(5);
        }
    }

    /** @brief Identity matrix */
    static LowerTriangular Identity(int /*size*/) {
        LowerTriangular result;
        result.data = Eigen::Matrix3d::Identity();
        return result;
    }

    /** @brief From full symmetric matrix (symmetrized to handle numerical noise) */
    static LowerTriangular fromFullMatrix(const Eigen::Matrix3d& full) {
        LowerTriangular result;
        Eigen::Matrix3d symmetric = 0.5 * (full + full.transpose());
        result.data.template triangularView<Eigen::Lower>() =
            symmetric.template triangularView<Eigen::Lower>();
        return result;
    }

    /** @brief Element access — (i,j) and (j,i) read the same lower-triangle value */
    inline double operator()(int i, int j) const {
        return i >= j ? data(i, j) : data(j, i);
    }

    /** @brief Element write access — (i,j) and (j,i) write the same lower-triangle value */
    inline double& operator()(int i, int j) {
        return i >= j ? data(i, j) : data(j, i);
    }

    /** @brief Access internal Matrix3d (symmetric) */
    inline const Eigen::Matrix3d& getSymmetricMatrix() const noexcept {
        symmetricData = data.selfadjointView<Eigen::Lower>();
        return symmetricData;
    }

    /** @brief Return packed lower-triangular storage (backward compat, 6 elements) */
    inline Eigen::VectorXd getData() const noexcept {
        Eigen::VectorXd packed(6);
        packed(0) = data(0,0);
        packed(1) = data(1,0); packed(2) = data(1,1);
        packed(3) = data(2,0); packed(4) = data(2,1); packed(5) = data(2,2);
        return packed;
    }

    /** @brief Return full 3×3 matrix */
    inline Eigen::MatrixXd getFullMatrix() const {
        return getSymmetricMatrix();
    }

    /** @brief Get matrix dimension (always 3) */
    inline int getSize() const noexcept { return 3; }

    /** @brief Symmetric matrix-vector multiply */
    inline Eigen::VectorXd multiplySymmetric(const Eigen::VectorXd& v) const {
        return getSymmetricMatrix() * v;
    }

    /** @brief 3-vector specialized symmetric multiply */
    inline Vector3d multiplySymmetric(const Vector3d& v) const {
        return getSymmetricMatrix() * v;
    }

    /** @brief Lower-triangular-only multiply (for completeness, not commonly used) */
    Vector3d operator*(const Vector3d& v) const {
        Vector3d result;
        result(0) = data(0, 0) * v(0);
        result(1) = data(1, 0) * v(0) + data(1, 1) * v(1);
        result(2) = data(2, 0) * v(0) + data(2, 1) * v(1) + data(2, 2) * v(2);
        return result;
    }

    Eigen::VectorXd operator*(const Eigen::VectorXd& v) const;

    /** @brief Matrix-matrix multiply (this * other) */
    LowerTriangular operator*(const LowerTriangular& other) const;

    /** @brief Multiplication with dense matrix */
    Eigen::MatrixXd operator*(const Eigen::MatrixXd& other) const;

    friend Eigen::MatrixXd operator*(const Eigen::MatrixXd& lhs, const LowerTriangular& rhs);

    /** @brief Scalar multiply */
    LowerTriangular operator*(double scalar) const {
        LowerTriangular result;
        result.data = data * scalar;
        return result;
    }

    friend LowerTriangular operator*(double scalar, const LowerTriangular& matrix) {
        return matrix * scalar;
    }

    /** @brief Matrix addition */
    LowerTriangular operator+(const LowerTriangular& other) const {
        LowerTriangular result;
        result.data = data + other.data;
        return result;
    }

    /** @brief Matrix subtraction */
    LowerTriangular operator-(const LowerTriangular& other) const {
        LowerTriangular result;
        result.data = data - other.data;
        return result;
    }

    /** @brief Matrix transpose — returns LowerTriangular with swapped indices */
    LowerTriangular transpose() const {
        LowerTriangular result;
        result.data(0,0) = data(0,0);
        result.data(1,0) = data(0,1); result.data(1,1) = data(1,1);
        result.data(2,0) = data(0,2); result.data(2,1) = data(1,2); result.data(2,2) = data(2,2);
        // Mirror
        result.data(0,1) = result.data(1,0);
        result.data(0,2) = result.data(2,0);
        result.data(1,2) = result.data(2,1);
        return result;
    }

    /** @brief Compute inverse */
    LowerTriangular inverse() const;

    /** @brief Output */
    friend std::ostream& operator<<(std::ostream& os, const LowerTriangular& matrix) {
        os << matrix.data;
        return os;
    }
};

using lt = LowerTriangular;

} // namespace SpatialAlgebra

#endif // LOWER_TRIANGULAR_H
