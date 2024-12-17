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
#include <Eigen/Geometry> // Include for RotationMatrix

using Vector6d = Eigen::Matrix<double, 6, 1>;

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
        Vector6d Inertia;     ///< Generalized inertia vector
        Eigen::Matrix3d H;    ///< Coupling matrix
        Vector6d M;           ///< Mass matrix components

    public:
        /**
         * @brief Construct articulated body inertia
         * @param inertia Generalized inertia vector
         * @param h Coupling matrix
         * @param M Mass matrix components
         */
        ArticulatedBodyInertia(const Vector6d &inertia, const Eigen::Matrix3d &h, const Vector6d &M);

        /** @brief Default constructor creating zero inertia */
        ArticulatedBodyInertia();

        // getters
        inline const Vector6d &getInertia() const { return Inertia; }
        inline const Eigen::Matrix3d &getH() const { return H; }
        inline const Vector6d &getM() const { return M; }

        /**
         * @brief Add two articulated body inertias
         * @param other ArticulatedBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         */
        inline ArticulatedBodyInertia operator+(const ArticulatedBodyInertia &other) const;

        /**
         * @brief Add rigid body inertia to articulated body inertia
         * @param other RigidBodyInertia to add
         * @return Combined ArticulatedBodyInertia
         */
        inline ArticulatedBodyInertia operator+(const RigidBodyInertia &other) const;

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
