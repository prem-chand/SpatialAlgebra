#include "MotionVector.h"

namespace SpatialAlgebra
{

    // Constructors
    MotionVector::MotionVector() : SpatialVector() {}

    MotionVector::MotionVector(const Vector3d &angular, const Vector3d &linear)
        : SpatialVector(angular, linear) {}

    MotionVector::MotionVector(const SpatialVector &other)
        : SpatialVector(other) {}

    // Accessors
    Vector3d MotionVector::getAngular() const
    {
        return angular;
    }

    Vector3d MotionVector::getLinear() const
    {
        return linear;
    }

    // Basic operations
    MotionVector MotionVector::operator+(const MotionVector &other) const
    {
        return MotionVector(angular + other.angular, linear + other.linear);
    }

    MotionVector MotionVector::operator-(const MotionVector &other) const
    {
        return MotionVector(angular - other.angular, linear - other.linear);
    }

    MotionVector MotionVector::operator*(double scalar) const
    {
        return MotionVector(angular * scalar, linear * scalar);
    }

    // Cross-product operations
    MotionVector MotionVector::crossMotion(const MotionVector &other) const
    {
        return MotionVector(angular.cross(other.angular), linear.cross(other.linear));
    }

    MotionVector MotionVector::crossForce(const MotionVector &other) const
    {
        return MotionVector(angular.cross(other.linear), linear.cross(other.angular));
    }

    // dot product
    double MotionVector::dot(const MotionVector &other) const
    {
        return angular.dot(other.angular) + linear.dot(other.linear);
    }

    // Printing
    void MotionVector::print() const
    {
        std::cout << "Angular: " << angular.transpose() << ", Linear: " << linear.transpose() << std::endl;
    }

} // namespace SpatialAlgebra
