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
    // X^{-T} = [R, -R*[t]x; 0, R]
    // return X^{-T} * vec = [R*(τ - r×f), R*f]
    // Featherstone (2008) Eq 2.44

    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation) * vec.getLinear();
    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation) *
        (vec.getAngular() - skew(translation) * vec.getLinear());
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
    // X = [R, 0; -R*r̂, R]
    // X^T = [R^T, r̂*R^T; 0, R^T] where r̂ = skew(translation)
    // return X^T * vec = [R^T*τ + r̂*R^T*f, R^T*f]
    //           = [R^T*τ + t×(R^T*f), R^T*f]
    Vector3d transformedLinear = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getLinear();
    Vector3d transformedAngular = static_cast<const Eigen::Matrix3d &>(rotation.transpose()) * vec.getAngular() + skew(translation) * transformedLinear;

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

    Eigen::Matrix3d Z = I.getSymmetricMatrix()
        + skew(translation) * skew(y) + skew(y) * skew(translation);
    Eigen::Matrix3d I_rotated = static_cast<const Eigen::Matrix3d &>(rotation)
        * Z * static_cast<const Eigen::Matrix3d &>(rotation).transpose();
    lt I_new = LowerTriangular::fromFullMatrix(I_rotated);

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
    auto I1 = static_cast<const Eigen::Matrix3d &>(rotation.transpose())
        * I.getSymmetricMatrix() * static_cast<const Eigen::Matrix3d &>(rotation);
    auto I2 = skew(translation) * skew(h_new);
    auto I3 = skew(h_new) * skew(translation);

    lt I_new = LowerTriangular::fromFullMatrix(I1 - I2 - I3);

    return RigidBodyInertia(m, h_new, I_new);
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::tformABI(const ArticulatedBodyInertia &Ia) const
{
    // I' = X * Ia * X^T using block 3×3 formulas (avoids 6×6 temporaries)
    //
    // X = [R, 0; -R·skew(t), R]
    // Ia = [I, H; H^T, M]
    //
    // Block result:
    //   I'_00 = R·I·R^T
    //   I'_01 = R·(I·skew(t) + H)·R^T
    //   I'_11 = R·(-skew(t)·I·skew(t) - skew(t)·H + H^T·skew(t) + M)·R^T

    auto M_lt = Ia.getM();
    auto H = Ia.getH();
    auto I_lt = Ia.getInertia();

    const Matrix3d& R = static_cast<const Eigen::Matrix3d &>(rotation);
    Matrix3d r_skew = skew(translation);  // skew(t)
    Matrix3d Rt = R.transpose();

    // Get symmetric matrices
    Matrix3d I_sym = I_lt.getSymmetricMatrix();
    Matrix3d M_sym = M_lt.getSymmetricMatrix();

    // I'_00 = R·I·R^T
    Matrix3d I_new = R * I_sym * Rt;

    // I'_01 = R·(I·skew(t) + H)·R^T
    Matrix3d I_skew = I_sym * r_skew;      // I·skew(t)
    Matrix3d temp01 = I_skew + H;           // I·skew(t) + H
    Matrix3d H_new = R * temp01 * Rt;       // R·(I·skew(t)+H)·R^T

    // I'_11 = R·(-skew(t)·I·skew(t) - skew(t)·H + H^T·skew(t) + M)·R^T
    Matrix3d skew_I = r_skew * I_sym;       // skew(t)·I
    Matrix3d skew_I_skew = skew_I * r_skew; // skew(t)·I·skew(t)
    Matrix3d skew_H = r_skew * H;           // skew(t)·H
    Matrix3d Ht_skew = H.transpose() * r_skew; // H^T·skew(t)
    Matrix3d temp11 = -skew_I_skew - skew_H + Ht_skew + M_sym;
    Matrix3d M_new = R * temp11 * Rt;

    // Symmetrize rotational and mass blocks
    I_new = 0.5 * (I_new + I_new.transpose());
    M_new = 0.5 * (M_new + M_new.transpose());

    return ArticulatedBodyInertia(
        LowerTriangular::fromFullMatrix(I_new),
        H_new,
        LowerTriangular::fromFullMatrix(M_new));
}

ArticulatedBodyInertia SpatialAlgebra::PluckerTransform::invtformABI(const ArticulatedBodyInertia &Ia) const
{
    // I' = X^(-1) * Ia * X^(-T) using block 3×3 formulas
    //
    // X^(-1) = [R^T, 0; skew(t)·R^T, R^T]
    // X^(-T) = [R, -R·skew(t); 0, R]
    //
    // Block result:
    //   I'_00 = R^T·I·R
    //   I'_01 = R^T·H·R - R^T·I·R·skew(t)
    //   I'_11 = -skew(t)·R^T·I·R·skew(t) + skew(t)·R^T·H·R
    //           - R^T·H^T·R·skew(t) + R^T·M·R

    auto M_lt = Ia.getM();
    auto H = Ia.getH();
    auto I_lt = Ia.getInertia();

    const Matrix3d& R = static_cast<const Eigen::Matrix3d &>(rotation);
    Matrix3d r_skew = skew(translation);
    Matrix3d Rt = R.transpose();

    Matrix3d I_sym = I_lt.getSymmetricMatrix();
    Matrix3d M_sym = M_lt.getSymmetricMatrix();

    // I'_00 = R^T·I·R
    Matrix3d I_new = Rt * I_sym * R;

    // Shared intermediate: R^T·I·R (already in I_new)
    Matrix3d Rt_I_R = I_new;

    // I'_01 = R^T·H·R - R^T·I·R·skew(t)
    Matrix3d Rt_I_R_skew = Rt_I_R * r_skew;
    Matrix3d Rt_H_R = Rt * H * R;
    Matrix3d H_new = Rt_H_R - Rt_I_R_skew;

    // I'_11 = -skew(t)·R^T·I·R·skew(t) + skew(t)·R^T·H·R
    //         - R^T·H^T·R·skew(t) + R^T·M·R
    Matrix3d skew_Rt_I_R = r_skew * Rt_I_R;
    Matrix3d skew_Rt_I_R_skew = skew_Rt_I_R * r_skew;
    Matrix3d skew_Rt_H_R = r_skew * Rt_H_R;
    Matrix3d Rt_Ht_R = Rt * H.transpose() * R;
    Matrix3d Rt_Ht_R_skew = Rt_Ht_R * r_skew;
    Matrix3d Rt_M_R = Rt * M_sym * R;

    Matrix3d temp11 = -skew_Rt_I_R_skew + skew_Rt_H_R - Rt_Ht_R_skew + Rt_M_R;
    Matrix3d M_new = temp11;

    // Symmetrize
    I_new = 0.5 * (I_new + I_new.transpose());
    M_new = 0.5 * (M_new + M_new.transpose());

    return ArticulatedBodyInertia(
        LowerTriangular::fromFullMatrix(I_new),
        H_new,
        LowerTriangular::fromFullMatrix(M_new));
}

PluckerTransform PluckerTransform::inverse() const
{
    // X = [R, 0; -R[t]x, R]
    // X^{-1} = [R^T, 0; [t]x*R^T, R^T] = [R^T, 0; -R^T*[-R*t]x, R^T]
    // So the inverse has: R_inv = R^T, t_inv = -R*t (not -R^T*t!)

    const Matrix3d& R = static_cast<const Matrix3d &>(rotation);
    Rotation invRotation = rotation.transpose();
    Vector3d invTranslation = -R * translation;  // -R*t, not -R^T*t
    return PluckerTransform(invRotation, invTranslation);
}

PluckerTransform PluckerTransform::multiply(const PluckerTransform &X) const
{
    // X1 = [R1, 0; -R1[t1]x, R1]
    // X2 = [R2, 0; -R2[t2]x, R2]
    // X1 * X2 (apply X2 then X1): R = R1*R2, t = t2 + R2^T*t1
    // This ensures: (X1*X2)*v = X1*(X2*v)

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
