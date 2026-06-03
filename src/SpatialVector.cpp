#include "SpatialVector.h"
#include "SpatialUtils.h"
#include <Eigen/Dense>
#include <iostream>

using namespace SpatialAlgebra;
using namespace Eigen;

// constructor
SpatialVector::SpatialVector() : angular(Vector3d::Zero()), linear(Vector3d::Zero()) {}

SpatialVector::SpatialVector(const Vector3d &a, const Vector3d &l)
    : angular(a), linear(l)
{
#ifndef NDEBUG
    if (a.hasNaN() || l.hasNaN() ||
        a.array().isInf().any() || l.array().isInf().any()) {
        std::cerr << "WARNING: NaN or Inf detected in SpatialVector constructor\n";
    }
#endif
}

SpatialVector::SpatialVector(const SpatialVector &other)
    : angular(other.angular), linear(other.linear) {}

const Vector3d& SpatialVector::getAngular() const
{
    return angular;
}

const Vector3d& SpatialVector::getLinear() const
{
    return linear;
}

SpatialVector SpatialVector::operator+(const SpatialVector &other) const
{
    return SpatialVector((angular + other.angular).eval(), (linear + other.linear).eval());
}

SpatialVector SpatialVector::operator-(const SpatialVector &other) const
{
    return SpatialVector((angular - other.angular).eval(), (linear - other.linear).eval());
}

SpatialVector SpatialVector::operator*(double scalar) const
{
    return SpatialVector((angular * scalar).eval(), (linear * scalar).eval());
}

SpatialVector SpatialAlgebra::SpatialVector::crossMotion(const SpatialVector &other) const
{
    // V1 x V2 = [w1 x w2, w1 x v2 + v1 x w2]
    return SpatialVector(angular.cross(other.angular), linear.cross(other.angular) + angular.cross(other.linear));
}

SpatialVector SpatialAlgebra::SpatialVector::crossForce(const SpatialVector &other) const
{
    // Force×Force cross product: [τ1×τ2 + f1×f2; τ1×f2 - τ2×f1]
    return SpatialVector(
        angular.cross(other.angular) + linear.cross(other.linear),
        angular.cross(other.linear) - other.angular.cross(linear)
    );
}

double SpatialVector::dot(const SpatialVector &other) const
{
    return angular.dot(other.angular) + linear.dot(other.linear);
}

void SpatialVector::print() const
{
    std::cout << "Angular: " << angular.transpose() << std::endl;
    std::cout << "Linear: " << linear.transpose() << std::endl;
}