#ifndef ROTATION_H
#define ROTATION_H

/**
 * @file Rotation.h
 * @brief Class representing 3D rotations with various representations
 */

#include <Eigen/Dense>

/**
 * @brief Class for handling 3D rotations with multiple representations
 * 
 * Extends Eigen::Matrix3d to provide additional functionality for rotation operations
 * including conversions between different rotation representations.
 */
class Rotation : public Eigen::Matrix3d
{
public:
    /** @brief Default constructor initializing to identity rotation */
    Rotation();

    /**
     * @brief Construct from 3x3 rotation matrix
     * @param matrix Input rotation matrix
     */
    Rotation(const Eigen::Matrix3d &matrix);

    /**
     * @brief Construct from angle-axis representation
     * @param angleAxis Angle-axis rotation
     */
    Rotation(const Eigen::AngleAxisd &angleAxis);

    /**
     * @brief Construct from quaternion representation
     * @param quaternion Unit quaternion representing rotation
     */
    Rotation(const Eigen::Quaterniond &quaternion);

    /**
     * @brief Set rotation from angle-axis representation
     * @param angleAxis Angle-axis rotation
     */
    void setFromAngleAxis(const Eigen::AngleAxisd &angleAxis);

    /**
     * @brief Set rotation from quaternion representation
     * @param quaternion Unit quaternion
     */
    void setFromQuaternion(const Eigen::Quaterniond &quaternion);

    /** @brief Convert to angle-axis representation */
    Eigen::AngleAxisd toAngleAxis() const;

    /** @brief Convert to quaternion representation */
    Eigen::Quaterniond toQuaternion() const;

    /** @brief Compute inverse rotation */
    Rotation inverse() const;

    /** @brief Compute matrix transpose */
    Rotation transpose() const;

    /**
     * @brief Multiply with another rotation
     * @param other Right-hand side rotation
     * @return Combined rotation
     */
    Rotation operator*(const Rotation &other) const;

    /**
     * @brief Apply rotation to a vector
     * @param vector 3D vector to rotate
     * @return Rotated vector
     */
    Eigen::Vector3d operator*(const Eigen::Vector3d &vector) const;
};

#endif // ROTATION_H
