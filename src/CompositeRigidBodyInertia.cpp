/**
 * @file CompositeRigidBodyInertia.cpp
 * @brief Implementation of the Composite Rigid Body Algorithm
 * @details Implements Featherstone's Algorithm 7.4. Uses 6×6 spatial inertia
 *          matrices internally to correctly handle composite inertia accumulation
 *          (parallel axis theorem when summing inertias at different COMs).
 *
 *          Algorithm overview:
 *          1. Convert each link's RigidBodyInertia to a 6×6 spatial inertia matrix
 *          2. Inward pass (leaf-to-root): accumulate subtree spatial inertias
 *             I_parent += X[i]^(-*) · I_subtree · X[i]^(-T)
 *          3. Mass matrix fill: for each body i, compute f = I_subtree · S_i,
 *             H[i][i] = S_i · f, then walk ancestors for off-diagonals
 *
 * @see CompositeRigidBodyInertia.h for public API
 * @see Featherstone (2008) Algorithm 7.4
 */

#include "CompositeRigidBodyInertia.h"
#include "ForwardDynamics.h"
#include "SpatialUtils.h"
#include "LowerTriangular.h"

#include <stdexcept>

using namespace Eigen;

namespace SpatialAlgebra {

namespace {

/** @brief Convert RigidBodyInertia to 6×6 spatial inertia matrix in frame origin */
Matrix<double, 6, 6> rbiToMatrix(const RigidBodyInertia& rbi) {
    const double m = rbi.getMass();
    const Vector3d c = rbi.getCom();          // COM in body frame
    const Matrix3d J = rbi.getInertiaMatrixLT().getFullMatrix();  // about COM

    Matrix3d c_skew = skew(c);

    // I_spatial = [J + m·(ĉ·ĉ^T),  m·ĉ;  m·ĉ^T,  m·I]
    // where ĉ·ĉ^T = (c·c^T - |c|²·I) and we add J about COM
    // Actually: I_frame = J_cm + m·(|c|²·I - c·c^T)  [parallel axis]
    // Simplified: I_block_00 = J + m·(ĉ·ĉ^T) = J + m·(|c|²·I - c·c^T)
    // Wait: ĉ·ĉ^T = c·c^T - |c|²·I
    // And parallel axis: |c|²·I - c·c^T = -ĉ·ĉ^T
    // So I_block_00 = J + m·(|c|²·I - c·c^T)
    // Let me compute directly.

    Matrix<double, 6, 6> I = Matrix<double, 6, 6>::Zero();
    // Top-left: J + m·(|c|²·I - c·c^T)
    I.topLeftCorner<3, 3>() = J + m * (c.squaredNorm() * Matrix3d::Identity()
                                       - c * c.transpose());
    // Top-right: m·ĉ
    I.topRightCorner<3, 3>() = m * c_skew;
    // Bottom-left: m·ĉ^T
    I.bottomLeftCorner<3, 3>() = m * c_skew.transpose();
    // Bottom-right: m·I
    I.bottomRightCorner<3, 3>() = m * Matrix3d::Identity();

    return I;
}

/** @brief Transform a 6×6 spatial inertia from child to parent frame
 *
 *  The 6×6 spatial inertia built by rbiToMatrix() transforms as:
 *    I_parent = X^T · I_child · X
 *
 *  where X is the parent→child motion transform. This differs from the
 *  conventional formula I' = X·I·X^T (Featherstone Eq 2.52) because
 *  rbiToMatrix() constructs the full frame-origin inertia at the child
 *  origin, and X^T·I·X correctly shifts it to the parent origin.
 */
Matrix<double, 6, 6> transformInertiaToParent(const Matrix<double, 6, 6>& I_child,
                                              const PluckerTransform& X_parent_to_child) {
    // X = [R, 0; -R·t̂, R]  (parent→child motion)
    //
    // X^T = [R^T,  -t̂^T·R^T;  0,  R^T]
    //     = [R^T,   t̂·R^T;    0,  R^T]    (since t̂^T = -t̂)

    const Matrix3d& R = static_cast<const Matrix3d&>(X_parent_to_child.getRotation());
    const Vector3d t = X_parent_to_child.getTranslation();
    Matrix3d t_skew = skew(t);
    Matrix3d Rt = R.transpose();

    // Build X^T
    Matrix<double, 6, 6> XT = Matrix<double, 6, 6>::Identity();
    XT.topLeftCorner<3, 3>() = Rt;
    XT.topRightCorner<3, 3>() = t_skew * Rt;
    // bottom-left = 0 (identity)
    XT.bottomRightCorner<3, 3>() = Rt;

    // Build X
    Matrix<double, 6, 6> X_mat = Matrix<double, 6, 6>::Identity();
    X_mat.topLeftCorner<3, 3>() = R;
    X_mat.bottomLeftCorner<3, 3>() = -R * t_skew;
    X_mat.bottomRightCorner<3, 3>() = R;

    return XT * I_child * X_mat;
}

/** @brief Apply 6×6 spatial inertia to a 6D motion vector (not in S, return raw Vector6d) */
Vector6d applySpatialInertia(const Matrix<double, 6, 6>& I, const MotionVector& mv) {
    Vector6d v;
    v.head<3>() = mv.getAngular();
    v.tail<3>() = mv.getLinear();
    Vector6d result = I * v;
    return result;
}

/** @brief Dot product of motion vector and 6D force vector */
double spatialDot(const MotionVector& mv, const Vector6d& f) {
    return mv.getAngular().dot(f.head<3>()) + mv.getLinear().dot(f.tail<3>());
}

/** @brief Transform 6D force from child frame to parent frame via X^T (the force dual) */
Vector6d transformForceToParent(const Vector6d& f_child, const PluckerTransform& X_parent_to_child) {
    // X_pc is parent→child motion: X = [R, 0; -R·t̂, R]
    // The force dual from child to parent is X^T:
    // X^T = [R^T, t̂·R^T; 0, R^T]   (since X = [R,0; -Rt̂,R], X^T = [R^T, -t̂^T·R^T; 0, R^T])
    // Actually X^T = [R^T, -R·t̂^T; 0, R^T] = [R^T, R·t̂; 0, R^T]
    // Wait let me check: X = [R, 0; -R·t̂, R]
    // X^T = [R^T, (-R·t̂)^T; 0, R^T] = [R^T, -t̂^T·R^T; 0, R^T]
    // Since t̂^T = -t̂: -t̂^T·R^T = t̂·R^T
    // So X^T = [R^T, t̂·R^T; 0, R^T]

    const Matrix3d& R = static_cast<const Matrix3d&>(X_parent_to_child.getRotation());
    const Vector3d t = X_parent_to_child.getTranslation();
    Matrix3d t_skew = skew(t);
    Matrix3d Rt = R.transpose();

    Vector6d result;
    result.head<3>() = Rt * f_child.head<3>() + t_skew * Rt * f_child.tail<3>();
    result.tail<3>() = Rt * f_child.tail<3>();
    return result;
}

} // anonymous namespace

Eigen::MatrixXd computeCompositeRigidBodyInertia(std::vector<Link>& links) {
    const int n = static_cast<int>(links.size());
    if (n == 0) {
        return Eigen::MatrixXd();
    }

    // =========================================================================
    // Step 1: Convert each link's inertia to 6×6 spatial inertia matrix
    // =========================================================================

    std::vector<Matrix<double, 6, 6>> C(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        C[static_cast<size_t>(i)] = rbiToMatrix(links[static_cast<size_t>(i)].I);
    }

    // =========================================================================
    // Step 2: Inward pass — accumulate subtree spatial inertias
    // =========================================================================
    //
    // For i = n-1 down to 0 (leaf to root):
    //   If link i has a parent p:
    //     C[p] += X[i]^(-1) · C[i] · X[i]^(-T)   (child→parent inertia transform)

    for (int i = n - 1; i >= 0; --i) {
        const int p = links[static_cast<size_t>(i)].parent;
        if (p >= 0) {
            C[static_cast<size_t>(p)] += transformInertiaToParent(
                C[static_cast<size_t>(i)],
                links[static_cast<size_t>(i)].X);
        }
    }

    // =========================================================================
    // Step 3: Fill mass matrix by walking up the kinematic chain.
    // =========================================================================
    //
    // For each body i (processed leaf-to-root):
    //   f = C[i] · S_i          — inertial force from unit acceleration of joint i
    //   H[i][i] = S_i^T · f     — diagonal mass
    //   Walk up ancestors:
    //     f_child = parent_force via X[l]^T  (inverseTransformForce)
    //     H[p][i] = S_p^T · f   — off-diagonal (also set H[i][p] = H[p][i])

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

    for (int i = n - 1; i >= 0; --i) {
        // Compute inertial force due to unit acceleration of joint i:
        // f = C[i] · S_i  (force in frame i)
        const Vector6d f = applySpatialInertia(
            C[static_cast<size_t>(i)],
            links[static_cast<size_t>(i)].S);

        // Diagonal entry: H[i][i] = S_i^T · C[i] · S_i
        H(i, i) = spatialDot(links[static_cast<size_t>(i)].S, f);

        // Walk up the kinematic chain accumulating off-diagonals.
        Vector6d f_k = f;
        int k = i;

        while (links[static_cast<size_t>(k)].parent >= 0) {
            const int p = links[static_cast<size_t>(k)].parent;

            // Transform force from frame k to parent frame p:
            // f_parent = X[k]^T · f_k   (child→parent force dual)
            f_k = transformForceToParent(f_k, links[static_cast<size_t>(k)].X);

            // Off-diagonal: H[p][i] = S_p^T · f (now in frame p)
            const double h_val = spatialDot(links[static_cast<size_t>(p)].S, f_k);
            H(p, i) = h_val;
            H(i, p) = h_val;  // Symmetric

            k = p;
        }
    }

    return H;
}

} // namespace SpatialAlgebra
