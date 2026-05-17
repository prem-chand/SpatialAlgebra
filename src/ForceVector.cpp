#include "ForceVector.h"
#include "SpatialUtils.h"

namespace SpatialAlgebra
{
    // Constructors
    ForceVector::ForceVector() : SpatialVector() {}

    ForceVector::ForceVector(const Vector3d &angular, const Vector3d &linear)
        : SpatialVector(angular, linear) {}

    ForceVector::ForceVector(const SpatialVector &other)
        : SpatialVector(other) {}

    // Basic operations
    ForceVector ForceVector::operator+(const ForceVector &other) const
    {
        return ForceVector(this->angular + other.angular, this->linear + other.linear);
    }

    ForceVector ForceVector::operator-(const ForceVector &other) const
    {
        return ForceVector(this->angular - other.angular, this->linear - other.linear);
    }

    ForceVector ForceVector::operator*(double scalar) const
    {
        return ForceVector(this->angular * scalar, this->linear * scalar);
    }

    // Cross-product operations
    ForceVector ForceVector::crossForce(const ForceVector &other) const
    {
        return cross(*this, other);
    }

    // dot product
    double ForceVector::dot(const ForceVector &other) const
    {
        return this->angular.dot(other.angular) + this->linear.dot(other.linear);
    }

    // Printing
    void ForceVector::print() const
    {
        std::cout << "Angular: " << this->angular.transpose() << ", Linear: " << this->linear.transpose() << std::endl;
    }

} // namespace SpatialAlgebra
