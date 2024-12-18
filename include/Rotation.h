#ifndef ROTATION_H
#define ROTATION_H

/**
 * @file Rotation.h
 * @brief Class representing 3D rotations with various representations and operations
 * @details This file defines the Rotation class which provides a comprehensive interface
 *          for handling 3D rotations. It supports multiple rotation representations including:
 *          - 3x3 rotation matrices (inherits from Eigen::Matrix3d)
 *          - Angle-axis representation (via Eigen::AngleAxisd)
 *          - Quaternions (via Eigen::Quaterniond)
 * 
 *          The class ensures proper conversion between these representations and provides
 *          operations like composition, inversion, and vector transformation.
 * 
 * @note All rotation matrices are guaranteed to be orthogonal with determinant +1
 * 
 * Example usage:
 * @code{.cpp}
 *     // Create a rotation from angle-axis representation
 *     Eigen::AngleAxisd angleAxis(M_PI/4, Eigen::Vector3d::UnitZ()); // 45° around Z
 *     Rotation rot1(angleAxis);
 * 
 *     // Create another rotation from a quaternion
 *     Eigen::Quaterniond quat(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
 *     Rotation rot2(quat);
 * 
 *     // Compose rotations
 *     Rotation combined = rot1 * rot2;
 * 
 *     // Transform a vector
 *     Eigen::Vector3d vec(1, 0, 0);
 *     Eigen::Vector3d rotated = combined * vec;
 * @endcode
 */

#include <Eigen/Dense>

/**
 * @brief Class for handling 3D rotations with multiple representations
 * @details The Rotation class extends Eigen::Matrix3d to provide a robust interface
 *          for working with 3D rotations. It maintains the underlying rotation matrix
 *          representation while offering convenient conversions to and from other
 *          common rotation representations.
 * 
 *          Key features:
 *          - Automatic handling of rotation composition
 *          - Conversion between different rotation representations
 *          - Vector transformation operations
 *          - Proper handling of rotation inverse and transpose
 * 
 * @warning When constructing from matrices, the input should be a valid rotation matrix.
 *          The class does not currently enforce orthogonality or proper determinant.
 */
class Rotation : public Eigen::Matrix3d
{
public:
    /** 
     * @brief Default constructor initializing to identity rotation
     * @details Creates a rotation matrix representing no rotation (identity matrix)
     */
    Rotation();

    /**
     * @brief Construct from 3x3 rotation matrix
     * @param matrix Input rotation matrix
     * @warning The input matrix should be a valid rotation matrix (orthogonal with det = +1)
     */
    Rotation(const Eigen::Matrix3d &matrix);

    /**
     * @brief Construct from angle-axis representation
     * @param angleAxis Angle-axis rotation (angle in radians)
     * @details Converts the angle-axis representation to the internal rotation matrix format.
     *          The angle is specified in radians, and the axis should be a unit vector.
     */
    Rotation(const Eigen::AngleAxisd &angleAxis);

    /**
     * @brief Construct from quaternion representation
     * @param quaternion Unit quaternion representing rotation
     * @details Converts the quaternion representation to the internal rotation matrix format.
     *          The quaternion should be normalized (unit quaternion).
     */
    Rotation(const Eigen::Quaterniond &quaternion);

    /**
     * @brief Set rotation from angle-axis representation
     * @param angleAxis Angle-axis rotation
     * @details Updates the internal rotation matrix based on the provided angle-axis representation.
     *          The angle should be in radians and the axis should be a unit vector.
     */
    void setFromAngleAxis(const Eigen::AngleAxisd &angleAxis);

    /**
     * @brief Set rotation from quaternion representation
     * @param quaternion Unit quaternion
     * @details Updates the internal rotation matrix based on the provided quaternion.
     *          The quaternion should be normalized.
     */
    void setFromQuaternion(const Eigen::Quaterniond &quaternion);

    /**
     * @brief Convert to angle-axis representation
     * @return Angle-axis representation of the rotation
     * @details The returned angle is in radians and the axis is a unit vector
     */
    Eigen::AngleAxisd toAngleAxis() const;

    /**
     * @brief Convert to quaternion representation
     * @return Unit quaternion representation of the rotation
     * @details The returned quaternion is guaranteed to be normalized
     */
    Eigen::Quaterniond toQuaternion() const;

    /**
     * @brief Compute inverse rotation
     * @return Inverse of the current rotation
     * @details For a rotation matrix R, its inverse equals its transpose: R^(-1) = R^T
     */
    Rotation inverse() const;

    /**
     * @brief Compute matrix transpose
     * @return Transpose of the current rotation matrix
     * @details For a valid rotation matrix, transpose equals inverse
     */
    Rotation transpose() const;

    /**
     * @brief Multiply with another rotation
     * @param other Right-hand side rotation
     * @return Combined rotation (this * other)
     * @details Composes two rotations: result = this * other
     *          The resulting rotation represents applying 'other' first, then 'this'
     */
    Rotation operator*(const Rotation &other) const;

    /**
     * @brief Apply rotation to a vector
     * @param vector 3D vector to rotate
     * @return Rotated vector
     * @details Transforms the input vector by this rotation: result = R * v
     */
    Eigen::Vector3d operator*(const Eigen::Vector3d &vector) const
    {
        return static_cast<const Eigen::Matrix3d&>(*this) * vector;
    }

    /**
     * @brief Multiply rotation matrix with a general 3x3 matrix from the right
     * @param matrix Right-hand side 3x3 matrix
     * @return Result of matrix multiplication
     * @note Result may not be a valid rotation matrix
     * @details Performs standard matrix multiplication: result = R * M
     */
    Eigen::Matrix3d operator*(const Eigen::Matrix3d &matrix) const;

    /**
     * @brief Friend operator for multiplication with a general 3x3 matrix from the left
     * @param matrix Left-hand side 3x3 matrix
     * @param rotation Right-hand side rotation
     * @return Result of matrix multiplication
     * @note Result may not be a valid rotation matrix
     * @details Performs standard matrix multiplication: result = M * R
     */
    friend Eigen::Matrix3d operator*(const Eigen::Matrix3d &matrix, const Rotation &rotation);
};

#endif // ROTATION_H
