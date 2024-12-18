#include "Rotation.h"

// Default constructor
Rotation::Rotation() : Eigen::Matrix3d(Eigen::Matrix3d::Identity()) {}

// Constructor from Eigen::Matrix3d
Rotation::Rotation(const Eigen::Matrix3d &matrix) : Eigen::Matrix3d(matrix) {}

// Constructor from angle-axis representation
Rotation::Rotation(const Eigen::AngleAxisd &angleAxis)
{
    *this = angleAxis.toRotationMatrix();
}

// Constructor from quaternion representation
Rotation::Rotation(const Eigen::Quaterniond &quaternion)
{
    *this = quaternion.toRotationMatrix();
}

// Function to set rotation from angle-axis
void Rotation::setFromAngleAxis(const Eigen::AngleAxisd &angleAxis)
{
    *this = angleAxis.toRotationMatrix();
}

// Function to set rotation from quaternion
void Rotation::setFromQuaternion(const Eigen::Quaterniond &quaternion)
{
    *this = quaternion.toRotationMatrix();
}

// Function to get the angle-axis representation
Eigen::AngleAxisd Rotation::toAngleAxis() const
{
    return Eigen::AngleAxisd(*this);
}

// Function to get the quaternion representation
Eigen::Quaterniond Rotation::toQuaternion() const
{
    return Eigen::Quaterniond(*this);
}

// inverse
Rotation Rotation::inverse() const
{
    return Rotation(this->transpose());
}

// transpose
Rotation Rotation::transpose() const
{
    return Rotation(Eigen::Matrix3d::transpose());
}

// operator* for Rotation-Rotation multiplication
Rotation Rotation::operator*(const Rotation &other) const
{
    return Rotation(static_cast<Eigen::Matrix3d>(*this) * static_cast<Eigen::Matrix3d>(other));
}
