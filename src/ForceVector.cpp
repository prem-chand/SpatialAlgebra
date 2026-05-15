#include "ForceVector.h"

namespace SpatialAlgebra
{
    // Constructors
    ForceVector::ForceVector() : SpatialVector() {}

    ForceVector::ForceVector(const Vector3d &angular, const Vector3d &linear)
        : SpatialVector(angular, linear) {}

    ForceVector::ForceVector(const SpatialVector &other)
        : SpatialVector(other) {}

    // Accessors
    Vector3d ForceVector::getAngular() const
    {
        return this->angular;
    }

    Vector3d ForceVector::getLinear() const
    {
        return this->linear;
    }

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
    ForceVector ForceVector::crossMotion(const ForceVector &other) const
    {
        // Motion cross product: [ω1×ω2; ω1×v2 + v1×ω2]
        return ForceVector(
            this->angular.cross(other.angular),
            this->linear.cross(other.angular) + this->angular.cross(other.linear)
        );
    }

    ForceVector ForceVector::crossForce(const ForceVector &other) const
    {
        // Force cross product: [τ1×τ2; τ1×f2 - τ2×f1]
        // This ensures anti-commutativity: a×b = -(b×a)
        return ForceVector(
            this->angular.cross(other.angular),
            this->angular.cross(other.linear) - other.angular.cross(this->linear)
        );
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
