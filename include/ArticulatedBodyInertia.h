#ifndef ARTIC_BODY_INERTIA_H
#define ARTIC_BODY_INERTIA_H

/**
 * @file ArticulatedBodyInertia.h
 * @brief Class representing articulated body inertia properties
 */

#include <array>
#include <iostream>
#include "RigidBodyInertia.h"
#include "SpatialUtils.h"
#include "LowerTriangular.h"
#include <Eigen/Geometry> // Include for RotationMatrix

using Vector6d = Eigen::Matrix<double, 6, 1>;
using lt = LowerTriangular;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing inertial properties of an articulated body
     *
     * Stores and manages generalized inertia matrix components for
     * articulated body dynamics calculations.
     */
    class ArticulatedBodyInertia
    {
    private:
        lt Inertia;      ///< Generalized inertia vector
        Eigen::Matrix3d H; ///< Coupling matrix
        lt M;            ///< Mass matrix components

    public:
        /**
         * @brief Construct articulated body inertia
         * @param inertia Generalized inertia vector
         * @param h Coupling matrix
         * @param M Mass matrix components
         */
        ArticulatedBodyInertia(const lt &inertia, const Eigen::Matrix3d &h, const lt &M) : Inertia(inertia), H(h), M(M) {}

        /** @brief Default constructor creating zero inertia */
        ArticulatedBodyInertia() : Inertia(lt(Vector6d::Zero(),3)), H(Eigen::Matrix3d::Zero()), M(lt(Vector6d::Zero(),3)) {}

        // getters
        inline const lt &getInertia() const { return Inertia; }
        inline const Eigen::Matrix3d &getH() const { return H; }
        inline const lt &getM() const { return M; }

        /**
         * @brief Add two articulated body inertias
         * @param other ArticulatedBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         */
        inline ArticulatedBodyInertia operator+(const ArticulatedBodyInertia &other) const
        {
            return ArticulatedBodyInertia(Inertia + other.Inertia, H + other.H, M + other.M);
        }

        /**
         * @brief Add rigid body inertia to articulated body inertia
         * @param other RigidBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         */
        inline ArticulatedBodyInertia operator+(const RigidBodyInertia &other) const
        {
            return ArticulatedBodyInertia(M + lt::Identity(3) * other.getMass(), H + skew(other.getCom()), Inertia + other.getInertiaMatrixLT());
        }

        inline ArticulatedBodyInertia operator*(double scalar) const
        {
            return ArticulatedBodyInertia(Inertia * scalar, H * scalar, M * scalar);
        }

        /**
         * @brief Apply articulated body inertia to motion vector
         * @param mv Motion vector to apply
         * @return Resulting force vector
         */
        inline fv apply(const mv &mv) const;

        /** @brief Print inertia properties */
        inline void print() const;
    };

    using abi = ArticulatedBodyInertia;

} // namespace SpatialAlgebra

#endif
