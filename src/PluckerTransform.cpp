#include "SpatialVector.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include "SpatialUtils.h"
#include <Eigen/Dense>
#include <iostream>

using namespace SpatialAlgebra;
using namespace Eigen;

// constructor
PluckerTransform::PluckerTransform(const Rotation &r, const Vector3d &t)
    : rotation(r), translation(t) {}

SpatialVector PluckerTransform::transformMotion(const SpatialVector &vec) const
{
    // vec = [ω, v]
    // X = [R, 0; -R[t]x, R]

    // return X * vec = [Rω, R(v - [t]xω)]

    // mv(E^(T)omega,E^(T)v+r xxE^(T)omega)
    // $\operatorname{mv}\left(\boldsymbol{E}^{\mathrm{T}} \boldsymbol{\omega}, \boldsymbol{E}^{\mathrm{T}} \boldsymbol{v}+\boldsymbol{r} \times \boldsymbol{E}^{\mathrm{T}} \boldsymbol{\omega}\right)$

    Vector3d transformedAngular = rotation * vec.getAngular();
    Vector3d transformedLinear = rotation * (vec.getLinear() + skew(translation) * vec.getAngular());

    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector PluckerTransform::transformForce(const SpatialVector &vec) const
{
    // vec = [τ, f]
    // X = [R, 0; -R[t]x, R]
    // return X^{-T} * vec = [R(τ - [t]xf), Rf]

    Vector3d transformedAngular = rotation * (vec.getAngular() - skew(translation) * vec.getLinear());
    Vector3d transformedLinear = rotation * vec.getLinear();
    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector SpatialAlgebra::PluckerTransform::inverseTransformMotion(const SpatialVector &vec) const
{
    // vec = [ω, v]
    // X = [R, 0; -R[t]x, R]

    // return X^{-1} * vec = [R^Tω, R^Tv + [t]xR^Tω)]
    Vector3d transformedAngular = rotation.transpose() * vec.getAngular();
    Vector3d transformedLinear = rotation.transpose() * vec.getLinear() + skew(translation) * transformedAngular;
    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector SpatialAlgebra::PluckerTransform::inverseTransformForce(const SpatialVector &vec) const
{
    // vec = [τ, f]
    // X = [R, 0; -R[t]x, R]

    // return X^{T} * vec = [R^Tτ + [t]xR^Tf, R^Tf ]
    Vector3d transformedLinear = rotation.transpose() * vec.getLinear();
    Vector3d transformedAngular = rotation.transpose() * vec.getAngular() + skew(translation) * transformedLinear;

    return SpatialVector(transformedAngular, transformedLinear);
}

RigidBodyInertia SpatialAlgebra::PluckerTransform::tformRBI(const RigidBodyInertia &Ihat) const
{
    return RigidBodyInertia();
}

RigidBodyInertia SpatialAlgebra::PluckerTransform::invtformRBI(const RigidBodyInertia &Ihat) const
{
    return RigidBodyInertia();
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::tformABI(const ArticulatedBodyInertia &Ia) const
{
    return ArticulatedBodyInertia();
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::invtformABI(const ArticulatedBodyInertia &Ia) const
{
    return ArticulatedBodyInertia();
}

PluckerTransform PluckerTransform::inverse() const
{
    // X = [R, 0; -R[t]x, R]
    // return X^-1 = plux([R^T, -Rx])

    Rotation invRotation = rotation.transpose();
    Vector3d invTranslation = -invRotation * translation;
    return PluckerTransform(invRotation, invTranslation);
}

auto PluckerTransform::multiply(const PluckerTransform &X) const
{
    // X1 = [R1, 0; -R1[t1]x, R1]
    // X2 = [R2, 0; -R2[t2]x, R2]
    // X1 * X2 = [R1R2, 0; -R1R2[t2]x + R1[t1]x, R1R2]

    // TODO: how to ensure product to 2 rotation matrices is still a rotation matrix upto finite precision?
    Rotation newRotation = rotation * X.rotation;
    Vector3d newTranslation = X.translation + X.rotation.transpose() * translation;
    return PluckerTransform(newRotation, newTranslation);
}

PluckerTransform PluckerTransform::apply(const PluckerTransform &X) const
{
    return multiply(X);
}

void PluckerTransform::print() const
{
    std::cout << "Rotation matrix:" << std::endl;
    std::cout << rotation << std::endl;
    std::cout << "Translation vector: " << translation.transpose() << std::endl;
}