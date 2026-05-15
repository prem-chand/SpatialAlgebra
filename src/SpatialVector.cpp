#include "SpatialVector.h"
#include <Eigen/Dense>
#include <iostream>

using namespace SpatialAlgebra;
using namespace Eigen;

// constructor
SpatialVector::SpatialVector() : angular(Vector3d::Zero()), linear(Vector3d::Zero()) {}

SpatialVector::SpatialVector(const Vector3d &a, const Vector3d &l)
    : angular(Vector3d(a.data())), linear(Vector3d(l.data())) {}

SpatialVector::SpatialVector(const SpatialVector &other)
    : angular(other.angular), linear(other.linear) {}

Vector3d SpatialVector::getAngular() const
{
    return {angular[0], angular[1], angular[2]};
}

Vector3d SpatialVector::getLinear() const
{
    return {linear[0], linear[1], linear[2]};
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
    // Force cross product: [τ1×τ2 + f1×f2; τ1×f2] per Featherstone
    return SpatialVector(
        angular.cross(other.angular) + linear.cross(other.linear),
        angular.cross(other.linear)
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