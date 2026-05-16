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
    // Featherstone (2008) Eq 2.43: ω' = R*ω, v' = R*(v - r×ω)

    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation) * vec.getAngular();
    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation) * (vec.getLinear() - skew(translation) * vec.getAngular());

    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector PluckerTransform::transformForce(const SpatialVector &vec) const
{
    // vec = [τ, f]
    // X = [R, 0; -R[t]x, R]
    // return X^{-T} * vec = [R(τ + [t]xf), Rf]
    // Featherstone (2008) Eq 2.44: τ' = R*(τ + r×f), f' = R*f

    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation) * (vec.getAngular() + skew(translation) * vec.getLinear());
    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation) * vec.getLinear();
    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector SpatialAlgebra::PluckerTransform::inverseTransformMotion(const SpatialVector &vec) const
{
    // vec = [ω, v]
    // X^{-1} = [R^T, 0; skew(t)*R^T, R^T]
    // return X^{-1} * vec = [R^T*ω, R^T*v + t×(R^T*ω)]
    // Featherstone (2008) Eq 2.46
    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getAngular();
    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getLinear() + skew(translation) * transformedAngular;
    return SpatialVector(transformedAngular, transformedLinear);
}

SpatialVector SpatialAlgebra::PluckerTransform::inverseTransformForce(const SpatialVector &vec) const
{
    // vec = [τ, f]
    // X^{-T} = [R^T, -skew(t)*R^T; 0, R^T]
    // return X^{T} * vec = [R^T*τ - skew(t)*R^T*f, R^T*f]
    //                    = [R^T*(τ - t×f), R^T*f]
    // Featherstone (2008) Eq 2.46
    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getLinear();
    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getAngular() - skew(translation) * transformedLinear;

    return SpatialVector(transformedAngular, transformedLinear);
}

RigidBodyInertia SpatialAlgebra::PluckerTransform::tformRBI(const RigidBodyInertia &Ihat) const
{
    // Featherstone (2008) Eq 2.52: I' = X*I*X^T
    // m' = m
    // h' = R*(h - m*r)
    // I' = R*(I + r̂*ĥ + ĥ*r̂)*R^T where ĥ = h - m*r
    auto m = Ihat.getMass();
    auto h = Ihat.getCom();
    auto I = Ihat.getInertiaMatrixLT();

    auto y = h - m * translation;  // y = h - m*r = ĥ
    auto h_new = static_cast<const Eigen::Matrix3d &>(rotation) * y;

    auto Z = I + LowerTriangular::fromFullMatrix(skew(translation) * skew(y) + skew(y) * skew(translation));
    lt I_new = LowerTriangular::fromFullMatrix(static_cast<const Eigen::Matrix3d &>(rotation) * Z * static_cast<const Eigen::Matrix3d &>(rotation).transpose());

    return RigidBodyInertia(m, h_new, I_new);
}

RigidBodyInertia SpatialAlgebra::PluckerTransform::invtformRBI(const RigidBodyInertia &Ihat) const
{
    // Inverse transform: I' = X^{-1}*I*X^{-T}
    // Using inverse transform parameters: R' = R^T, t' = -R^T*t
    // m' = m
    // h' = R^T*h + m*t
    // I' = R^T*I*R - (r̂*ĥ + ĥ*r̂) where ĥ = h_new
    auto m = Ihat.getMass();
    auto h = Ihat.getCom();
    auto I = Ihat.getInertiaMatrixLT();

    auto h_new = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * h + m * translation;
    auto I1 = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * I * static_cast<const Eigen::Matrix3d &>(rotation);
    auto I2 = skew(translation) * skew(h_new);
    auto I3 = skew(h_new) * skew(translation);

    lt I_new = LowerTriangular::fromFullMatrix(I1 - I2 - I3);

    return RigidBodyInertia(m, h_new, I_new);
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::tformABI(const ArticulatedBodyInertia &Ia) const
{
    // Featherstone (2008) Eq 7.16: I'a = X * Ia * X^T
    // 
    // 6x6 block structure: Ia = [I,   H;
    //                           H^T, M]
    // where I and M are symmetric (3x3), H is general (3x3)
    //
    // Spatial motion transform: X = [R,     0;
    //                                -R*r̂,  R]
    // where r̂ = skew(translation)
    //
    // We compute I'a = X * Ia * X^T directly using 6x6 matrices.
    
    auto M = Ia.getM();
    auto H = Ia.getH();
    auto Inertia = Ia.getInertia();

    const Matrix3d& R = static_cast<const Eigen::Matrix3d &>(rotation);
    Matrix3d r_skew = skew(translation);  // r̂ = skew(t)
    
    // Build 6x6 articulated body inertia matrix
    // I and M are symmetric, so we need to reconstruct the full symmetric matrix
    // from LowerTriangular storage (copy lower triangle to upper)
    Matrix3d I_full = Inertia.getFullMatrix();
    I_full = 0.5 * (I_full + I_full.transpose());  // Symmetrize
    
    Matrix3d M_full = M.getFullMatrix();
    M_full = 0.5 * (M_full + M_full.transpose());  // Symmetrize
    
    // Order: [angular; linear] so Ia = [I, H; H^T, M]
    MatrixXd Ia_6x6 = MatrixXd::Zero(6, 6);
    Ia_6x6.block<3,3>(0,0) = I_full;
    Ia_6x6.block<3,3>(0,3) = H;
    Ia_6x6.block<3,3>(3,0) = H.transpose();
    Ia_6x6.block<3,3>(3,3) = M_full;
    
    // Build 6x6 spatial transform X = [R, 0; -R*r̂, R]
    MatrixXd X_6x6 = MatrixXd::Zero(6, 6);
    X_6x6.block<3,3>(0,0) = R;
    X_6x6.block<3,3>(0,3) = Matrix3d::Zero();
    X_6x6.block<3,3>(3,0) = -R * r_skew;
    X_6x6.block<3,3>(3,3) = R;
    
    // Compute I'a = X * Ia * X^T
    MatrixXd Ia_prime = X_6x6 * Ia_6x6 * X_6x6.transpose();
    
    // Extract blocks - I' and M' should be symmetric
    Matrix3d I_new_full = Ia_prime.block<3,3>(0,0);
    Matrix3d M_new_full = Ia_prime.block<3,3>(3,3);
    // Symmetrize to handle numerical errors
    I_new_full = 0.5 * (I_new_full + I_new_full.transpose());
    M_new_full = 0.5 * (M_new_full + M_new_full.transpose());
    
    lt I_new = LowerTriangular::fromFullMatrix(I_new_full);
    Matrix3d H_new = Ia_prime.block<3,3>(0,3);
    lt M_new = LowerTriangular::fromFullMatrix(M_new_full);
    
    return ArticulatedBodyInertia(I_new, H_new, M_new);
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::invtformABI(const ArticulatedBodyInertia &Ia) const
{
    // Inverse transform: I' = X^{-1} * Ia * X^{-T}
    // 
    // For X = [R, 0; -R*r̂, R], the inverse is:
    // X^{-1} = [R^T, 0; r̂*R^T, R^T]
    // and X^{-T} = (X^{-1})^T = [R, -R*r̂; 0, R^T]
    //
    // We compute I' = X^{-1} * Ia * X^{-T} directly using 6x6 matrices.
    
    auto M = Ia.getM();
    auto H = Ia.getH();
    auto Inertia = Ia.getInertia();
    
    const Matrix3d& R = static_cast<const Eigen::Matrix3d &>(rotation);
    Matrix3d r_skew = skew(translation);  // r̂ = skew(t)
    
    // Symmetrize I and M from LowerTriangular storage
    Matrix3d I_full = Inertia.getFullMatrix();
    I_full = 0.5 * (I_full + I_full.transpose());
    
    Matrix3d M_full = M.getFullMatrix();
    M_full = 0.5 * (M_full + M_full.transpose());
    
    // Build 6x6 articulated body inertia matrix
    MatrixXd Ia_6x6 = MatrixXd::Zero(6, 6);
    Ia_6x6.block<3,3>(0,0) = I_full;
    Ia_6x6.block<3,3>(0,3) = H;
    Ia_6x6.block<3,3>(3,0) = H.transpose();
    Ia_6x6.block<3,3>(3,3) = M_full;
    
    // Build 6x6 inverse spatial transform X^{-1} = [R^T, 0; r̂*R^T, R^T]
    MatrixXd X_inv_6x6 = MatrixXd::Zero(6, 6);
    X_inv_6x6.block<3,3>(0,0) = R.transpose();
    X_inv_6x6.block<3,3>(0,3) = Matrix3d::Zero();
    X_inv_6x6.block<3,3>(3,0) = r_skew * R.transpose();
    X_inv_6x6.block<3,3>(3,3) = R.transpose();
    
    // Compute I' = X^{-1} * Ia * X^{-T}
    MatrixXd X_inv_T = X_inv_6x6.transpose();
    MatrixXd Ia_prime = X_inv_6x6 * Ia_6x6 * X_inv_T;
    
    // Extract blocks and symmetrize
    Matrix3d I_new_full = Ia_prime.block<3,3>(0,0);
    Matrix3d M_new_full = Ia_prime.block<3,3>(3,3);
    I_new_full = 0.5 * (I_new_full + I_new_full.transpose());
    M_new_full = 0.5 * (M_new_full + M_new_full.transpose());
    
    lt I_new = LowerTriangular::fromFullMatrix(I_new_full);
    Matrix3d H_new = Ia_prime.block<3,3>(0,3);
    lt M_new = LowerTriangular::fromFullMatrix(M_new_full);
    
    return ArticulatedBodyInertia(I_new, H_new, M_new);
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
    Vector3d newTranslation = X.translation + static_cast<const Eigen::Matrix3d &>(X.rotation.transpose()) * translation;
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
