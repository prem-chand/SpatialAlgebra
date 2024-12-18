#ifndef RIGID_BODY_INERTIA_H
#define RIGID_BODY_INERTIA_H

/**
 * @file RigidBodyInertia.h
 * @brief Class representing spatial rigid body inertia properties
 * @details This file implements the spatial inertia of rigid bodies, which combines
 *          mass, center of mass, and rotational inertia in a unified framework.
 *          The implementation uses a compact representation that allows efficient
 *          computation of dynamics quantities.
 * 
 *          Mathematical Representation:
 *          A spatial inertia I is represented by three components:
 *          - m: Mass scalar (kg)
 *          - h: First moment of mass / mass times COM (kg⋅m)
 *          - I: Rotational inertia matrix about COM (kg⋅m²)
 * 
 *          The 6x6 spatial inertia matrix form is:
 *          I = [ I    h×;
 *               h×ᵀ   mE ]
 *          where:
 *          - h× is the skew-symmetric matrix of h
 *          - E is the 3x3 identity matrix
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create inertia for a 1kg point mass at (1,0,0)
 *     double mass = 1.0;
 *     Vector3d com(1, 0, 0);
 *     lt inertia(3);  // Zero rotational inertia
 *     RigidBodyInertia rbi(mass, com, inertia);
 *     
 *     // Apply to a motion vector to get forces/torques
 *     MotionVector vel(Vector3d(0,0,1), Vector3d(1,0,0));
 *     ForceVector f = rbi.apply(vel);
 *     
 *     // Combine two inertias
 *     RigidBodyInertia total = rbi + rbi;
 * @endcode
 * 
 * @see MotionVector for the velocity representation
 * @see ForceVector for the force/torque representation
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2.4
 */

#include <array>
#include <iostream>
#include "SpatialVector.h"
#include "ForceVector.h"
#include "MotionVector.h"
#include "LowerTriangular.h"
#include <Eigen/Geometry>

using Vector3d = Eigen::Matrix<double, 3, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using lt = LowerTriangular;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing the spatial inertia of a rigid body
     * @details The RigidBodyInertia class provides a compact and efficient representation
     *          of rigid body inertial properties for spatial algebra computations.
     * 
     *          Physical Properties:
     *          - Symmetric positive definite (when physically valid)
     *          - Transforms covariantly under coordinate changes
     *          - Maps motion vectors to force vectors
     * 
     *          Implementation Features:
     *          - Stores rotational inertia in lower triangular form for efficiency
     *          - Provides direct access to mass, COM, and inertia components
     *          - Implements efficient spatial operations
     * 
     * @note The rotational inertia is stored about the center of mass frame
     *       to minimize numerical errors and simplify transformations.
     */
    class RigidBodyInertia
    {
    private:
        double mass;        ///< Mass of the rigid body (kg)
        Vector3d com;       ///< Center of mass position vector (m)
        lt inertiaMatrixLT; ///< Lower triangular part of rotational inertia matrix (kg⋅m²)

    public:
        /**
         * @brief Construct a rigid body inertia
         * @param mass Mass of the body (kg)
         * @param com Center of mass position (m)
         * @param inertiaMatrixLT Lower triangular part of rotational inertia matrix (kg⋅m²)
         * @details Creates a spatial inertia with specified mass properties. The rotational
         *          inertia should be specified about the center of mass.
         * @warning The inertia matrix must be positive definite for physical validity
         */
        inline RigidBodyInertia(double mass, const Vector3d &com, const lt &inertiaMatrixLT)
            : mass(mass), com(com), inertiaMatrixLT(inertiaMatrixLT) {}

        /**
         * @brief Default constructor creating zero inertia
         * @details Initializes all components to zero, representing a massless body.
         *          This is useful as a starting point for accumulating inertias.
         */
        RigidBodyInertia()
            : mass(0.0), com(Vector3d::Zero()), inertiaMatrixLT(lt(3)) {}

        /**
         * @brief Get the mass of the rigid body
         * @return Mass value in kg
         */
        inline double getMass() const { return mass; }

        /**
         * @brief Get the center of mass position
         * @return Center of mass vector in m
         */
        inline const Vector3d &getCom() const { return com; }

        /**
         * @brief Get the rotational inertia matrix
         * @return Lower triangular part of inertia matrix in kg⋅m²
         */
        inline const lt &getInertiaMatrixLT() const { return inertiaMatrixLT; }

        /**
         * @brief Add two rigid body inertias
         * @param other RigidBodyInertia to add
         * @return Combined RigidBodyInertia
         * @details Implements parallel combination of inertias, useful for
         *          computing composite body inertias. The result represents
         *          the inertia of a single rigid body equivalent to the two
         *          original bodies moving together.
         */
        inline RigidBodyInertia operator+(const RigidBodyInertia &other) const
        {
            return RigidBodyInertia(mass + other.mass, 
                                  com + other.com, 
                                  inertiaMatrixLT + other.inertiaMatrixLT);
        }

        /**
         * @brief Scale rigid body inertia
         * @param scalar Scaling factor (dimensionless)
         * @return Scaled RigidBodyInertia
         * @details Uniformly scales all inertial properties. This is useful for
         *          implementing symmetric scaling of bodies or applying mass
         *          distribution factors.
         */
        inline RigidBodyInertia operator*(double scalar) const
        {
            return RigidBodyInertia(mass * scalar, 
                                  com * scalar, 
                                  inertiaMatrixLT * scalar);
        }

        /**
         * @brief Apply rigid body inertia to motion vector
         * @param mv Motion vector to apply
         * @return Resulting force vector
         * @details Computes the force and torque resulting from the motion of
         *          the rigid body. The operation implements:
         *          f = I * v = [Iω + h×v; mv - h×ω]
         *          where:
         *          - ω is angular velocity
         *          - v is linear velocity
         *          - h = m*com is the first moment of mass
         * @note This operation is the fundamental mapping between motion and force
         *       vectors that defines the spatial inertia.
         */
        inline ForceVector apply(const MotionVector &mv) const
        {
            Vector3d omega = mv.getAngular();
            Vector3d v = mv.getLinear();
            Vector3d comCrossV = com.cross(v);
            Vector3d comCrossOmega = com.cross(omega);
            Vector3d Iomega = inertiaMatrixLT * omega;
            Vector3d force = mass * v - comCrossOmega;
            Vector3d torque = Iomega + comCrossV;
            return ForceVector(torque, force);
        }

        /**
         * @brief Print inertia properties
         * @details Outputs the mass, center of mass position, and rotational
         *          inertia matrix in a human-readable format.
         */
        inline void print() const
        {
            std::cout << "Mass: " << mass << " kg\n"
                      << "Center of mass: " << com.transpose() << " m\n"
                      << "Inertia matrix (kg⋅m²):\n" << inertiaMatrixLT << std::endl;
        }
    };

    /** @brief Type alias for more concise notation */
    using rbi = RigidBodyInertia;

} // namespace SpatialAlgebra

#endif // RIGID_BODY_INERTIA_H
