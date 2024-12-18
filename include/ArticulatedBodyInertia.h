#ifndef ARTIC_BODY_INERTIA_H
#define ARTIC_BODY_INERTIA_H

/**
 * @file ArticulatedBodyInertia.h
 * @brief Class representing articulated body inertia properties in spatial algebra
 * @details This file implements the generalized inertia of articulated bodies,
 *          which extends the concept of rigid body inertia to handle systems with
 *          internal degrees of freedom. This representation is crucial for efficient
 *          articulated body algorithms (ABA) in robotics and multibody dynamics.
 * 
 *          Mathematical Representation:
 *          An articulated body inertia Ia is represented by three components in
 *          6x6 block matrix form:
 *          Ia = [ I    H;
 *                Hᵀ    M ]
 *          where:
 *          - I: 3x3 rotational inertia matrix (kg⋅m²)
 *          - H: 3x3 coupling matrix (kg⋅m)
 *          - M: 3x3 mass matrix (kg)
 * 
 *          Key Properties:
 *          - Symmetric positive definite
 *          - Includes coupling effects between joints
 *          - Reduces to rigid body inertia for single bodies
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create components for articulated body inertia
 *     lt inertia(3);  // Rotational inertia
 *     Eigen::Matrix3d H = Eigen::Matrix3d::Identity();  // Coupling
 *     lt M(3);  // Mass matrix
 *     ArticulatedBodyInertia abi(inertia, H, M);
 *     
 *     // Add a rigid body inertia to it
 *     RigidBodyInertia rbi(1.0, Vector3d(1,0,0), lt(3));
 *     ArticulatedBodyInertia combined = abi + rbi;
 *     
 *     // Apply to a motion vector
 *     MotionVector vel(Vector3d(0,0,1), Vector3d(1,0,0));
 *     ForceVector f = combined.apply(vel);
 * @endcode
 * 
 * @see RigidBodyInertia for single body inertia representation
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7
 */

#include <array>
#include <iostream>
#include "RigidBodyInertia.h"
#include "SpatialUtils.h"
#include "LowerTriangular.h"
#include <Eigen/Geometry>

using Vector6d = Eigen::Matrix<double, 6, 1>;
using lt = LowerTriangular;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing inertial properties of an articulated body
     * @details The ArticulatedBodyInertia class provides a specialized representation
     *          for the composite inertia of articulated systems, including the
     *          effects of internal degrees of freedom and kinematic constraints.
     * 
     *          Physical Properties:
     *          - Symmetric positive definite (when physically valid)
     *          - Includes inertial coupling between joints
     *          - Transforms covariantly under coordinate changes
     * 
     *          Implementation Features:
     *          - Efficient block matrix structure
     *          - Supports composition with rigid body inertias
     *          - Optimized for articulated body algorithms
     * 
     * @note This representation is particularly efficient for the articulated body
     *       algorithm (ABA) and related computations in multibody dynamics.
     */
    class ArticulatedBodyInertia
    {
    private:
        lt Inertia;        ///< Rotational inertia matrix in lower triangular form (kg⋅m²)
        Eigen::Matrix3d H;  ///< Coupling matrix between rotational and linear motion (kg⋅m)
        lt M;              ///< Mass matrix in lower triangular form (kg)

    public:
        /**
         * @brief Construct articulated body inertia
         * @param inertia Rotational inertia matrix (kg⋅m²)
         * @param h Coupling matrix (kg⋅m)
         * @param M Mass matrix (kg)
         * @details Creates an articulated body inertia with specified components.
         *          The matrices must satisfy physical constraints to represent
         *          a valid inertia.
         * @warning The resulting matrix must be symmetric positive definite
         */
        ArticulatedBodyInertia(const lt &inertia, const Eigen::Matrix3d &h, const lt &M) 
            : Inertia(inertia), H(h), M(M) {}

        /**
         * @brief Default constructor creating zero inertia
         * @details Initializes all components to zero. This represents a
         *          massless articulated system and is useful as a starting
         *          point for accumulating inertias.
         */
        ArticulatedBodyInertia() 
            : Inertia(lt(Vector6d::Zero(),3)), 
              H(Eigen::Matrix3d::Zero()), 
              M(lt(Vector6d::Zero(),3)) {}

        /**
         * @brief Get the rotational inertia matrix
         * @return Rotational inertia in lower triangular form (kg⋅m²)
         */
        inline const lt &getInertia() const { return Inertia; }

        /**
         * @brief Get the coupling matrix
         * @return Coupling matrix between rotational and linear motion (kg⋅m)
         */
        inline const Eigen::Matrix3d &getH() const { return H; }

        /**
         * @brief Get the mass matrix
         * @return Mass matrix in lower triangular form (kg)
         */
        inline const lt &getM() const { return M; }

        /**
         * @brief Add two articulated body inertias
         * @param other ArticulatedBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         * @details Implements parallel combination of articulated body inertias.
         *          This operation is fundamental in the articulated body algorithm
         *          for combining inertias of connected bodies.
         */
        inline ArticulatedBodyInertia operator+(const ArticulatedBodyInertia &other) const
        {
            return ArticulatedBodyInertia(Inertia + other.Inertia, 
                                        H + other.H, 
                                        M + other.M);
        }

        /**
         * @brief Add rigid body inertia to articulated body inertia
         * @param other RigidBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         * @details Converts the rigid body inertia to articulated body form
         *          and combines it with this inertia. This is useful when
         *          building composite articulated bodies from rigid components.
         */
        inline ArticulatedBodyInertia operator+(const RigidBodyInertia &other) const
        {
            return ArticulatedBodyInertia(M + lt::Identity(3) * other.getMass(),
                                        H + skew(other.getCom()),
                                        Inertia + other.getInertiaMatrixLT());
        }

        /**
         * @brief Scale articulated body inertia
         * @param scalar Scaling factor (dimensionless)
         * @return Scaled ArticulatedBodyInertia
         * @details Uniformly scales all components of the inertia. This is
         *          useful for implementing symmetric scaling of articulated
         *          systems or applying distribution factors.
         */
        inline ArticulatedBodyInertia operator*(double scalar) const
        {
            return ArticulatedBodyInertia(Inertia * scalar,
                                        H * scalar,
                                        M * scalar);
        }

        /**
         * @brief Apply articulated body inertia to motion vector
         * @param mv Motion vector to apply
         * @return Resulting force vector
         * @details Computes the force and torque resulting from the motion of
         *          the articulated body. The operation implements:
         *          f = Ia * v = [Iω + Hv; Hᵀω + Mv]
         */
        inline fv apply(const mv &mv) const;

        /**
         * @brief Print inertia properties
         * @details Outputs the rotational inertia, coupling matrix, and mass
         *          matrix components in a human-readable format with units.
         */
        inline void print() const;
    };

    /** @brief Type alias for more concise notation */
    using abi = ArticulatedBodyInertia;

} // namespace SpatialAlgebra

#endif // ARTIC_BODY_INERTIA_H
