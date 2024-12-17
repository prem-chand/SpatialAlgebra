#ifndef RIGID_BODY_INERTIA_H
#define RIGID_BODY_INERTIA_H

/**
 * @file RigidBodyInertia.h
 * @brief Class representing rigid body inertia properties
 */

#include <array>
#include <iostream>
#include "SpatialVector.h"
#include "ForceVector.h"
#include "MotionVector.h"
#include <Eigen/Geometry> // Include for RotationMatrix

using Vector3d = Eigen::Matrix<double, 3, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing the inertial properties of a rigid body
     * 
     * Stores and manages mass, center of mass, and rotational inertia matrix
     * for rigid body dynamics calculations.
     */
    class RigidBodyInertia
    {
    private:
        double mass;              ///< Mass of the rigid body
        Vector3d com;            ///< Center of mass position vector
        Vector6d inertiaMatrixLT; ///< Lower triangular part of 3x3 rotational inertia matrix

    public:
        /**
         * @brief Construct a rigid body inertia
         * @param mass Mass of the body
         * @param com Center of mass position
         * @param inertiaMatrixLT Lower triangular part of inertia matrix
         */
        inline RigidBodyInertia(double mass, const Vector3d &com, const Vector6d &inertiaMatrixLT)
            : mass(mass), com(com), inertiaMatrixLT(inertiaMatrixLT) {}

        /** @brief Default constructor creating zero inertia */
        inline RigidBodyInertia()
        {
            RigidBodyInertia(0.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
        }

        // Getters
        inline double getMass() const { return mass; }
        inline const Vector3d &getCom() const { return com; }
        inline const Vector6d &getInertiaMatrixLT() const { return inertiaMatrixLT; }

        /**
         * @brief Add two rigid body inertias
         * @param other RigidBodyInertia to add
         * @return Combined RigidBodyInertia
         */
        inline RigidBodyInertia operator+(const RigidBodyInertia &other) const
        {
            // const double sumMass = mass + other.mass;
            // const Vector3d sumCom = com + other.com;
            return RigidBodyInertia(mass + other.mass, com + other.com, inertiaMatrixLT + other.inertiaMatrixLT);
        }

        /**
         * @brief Scale rigid body inertia
         * @param scalar Scaling factor
         * @return Scaled RigidBodyInertia 
         */
        inline RigidBodyInertia operator*(double scalar) const // Changed to pass by value
        {
            return RigidBodyInertia(mass * scalar, com * scalar, inertiaMatrixLT * scalar);
        }

        /**
         * @brief Apply rigid body inertia to motion vector
         * @param mv Motion vector to apply
         * @return Resulting force vector
         * @details Implements: [Iω + com × v; m*v - com × ω]
         */
        inline ForceVector apply(const MotionVector &mv) const
        {
            // rbi = [m, com, I_LT], mv = [ω, v]
            // rbi.apply(mv) = [Iω + com x v; m*v - com x ω]
            Vector3d omega = mv.getAngular();
            Vector3d v = mv.getLinear();
            Vector3d comCrossV = com.cross(v);
            Vector3d comCrossOmega = com.cross(omega);
            Vector3d Iomega = inertiaMatrixLT.head<3>().cwiseProduct(omega) + inertiaMatrixLT.tail<3>().cwiseProduct(omega);
            Vector3d force = mass * v - comCrossOmega;
            Vector3d torque = Iomega + comCrossV;
            return ForceVector(torque, force);
        }

        /** @brief Print inertia properties */
        inline void print() const
        {
            std::cout << "Mass: " << mass << '\n'
                      << "Center of mass: " << com.transpose() << '\n'
                      << "Inertia matrix (LT): " << inertiaMatrixLT.transpose() << std::endl;
        }
    };

    using rbi = RigidBodyInertia;

} // namespace SpatialAlgebra

#endif
