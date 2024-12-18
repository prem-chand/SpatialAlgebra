#ifndef FORCE_VECTOR_H
#define FORCE_VECTOR_H

/**
 * @file ForceVector.h
 * @brief Class representing spatial force vectors (wrenches) in rigid body dynamics
 * @details This file implements force vectors (also known as wrenches) in spatial algebra,
 *          which unify the representation of torques and forces in 3D space.
 *          Force vectors are essential for describing the complete force state
 *          acting on rigid bodies in a compact and geometrically meaningful way.
 * 
 *          Mathematical Representation:
 *          A force vector f = [τ; f] consists of:
 *          - τ: Torque/moment vector (N⋅m)
 *          - f: Linear force vector (N)
 *          
 *          Physical Interpretation:
 *          - τ represents the moment or torque around the reference point
 *          - f represents the linear force acting at the reference point
 *          - Together they form a "wrench" that completely describes force action
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create a pure torque around Z-axis (1 N⋅m)
 *     ForceVector torque(Vector3d(0,0,1), Vector3d(0,0,0));
 *     
 *     // Create a pure force along X-axis (1 N)
 *     ForceVector force(Vector3d(0,0,0), Vector3d(1,0,0));
 *     
 *     // Combine forces
 *     ForceVector combined = torque + force;
 *     
 *     // Compute force cross product (for force transformation)
 *     ForceVector crossed = torque.crossForce(force);
 * @endcode
 * 
 * @see MotionVector for the dual space representation of velocities
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2.
 */

#include <iostream>

#include "SpatialVector.h"

// Force vector class inherits from SpatialVector
namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial force vectors (wrenches)
     * @details The ForceVector class specializes SpatialVector for representing
     *          force quantities in spatial algebra. It provides type-safe operations
     *          specific to force vectors and ensures correct geometric interpretation.
     * 
     *          Key Concepts:
     *          - Force vectors form a vector space
     *          - They transform contravariantly under coordinate changes
     *          - The cross product represents force transformation in rigid bodies
     *          - Dot product with motion vectors yields power
     * 
     *          Common Applications:
     *          - Representing joint forces in articulated systems
     *          - Describing external loads on rigid bodies
     *          - Force propagation through kinematic chains
     * 
     * @note Force vectors and motion vectors are dual to each other, meaning they
     *       transform differently under coordinate changes and their dot product
     *       yields physically meaningful power.
     */
    class ForceVector : public SpatialVector
    {
    public:
        /** 
         * @brief Default constructor creating zero force
         * @details Initializes a force vector with zero torque and linear force
         */
        ForceVector();

        /**
         * @brief Construct from angular and linear components
         * @param angular Torque/moment component (N⋅m)
         * @param linear Force component (N)
         * @details Creates a force vector representing a combined torque and
         *          force action. The angular component represents torque around
         *          an axis, while the linear component represents force along an axis.
         */
        ForceVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Construct from spatial vector
         * @param other Spatial vector to convert
         * @details Converts a generic spatial vector to a force vector. This is useful
         *          when working with algorithms that operate on generic spatial vectors.
         */
        ForceVector(const SpatialVector &other);

        /**
         * @brief Get the torque/moment component
         * @return 3D torque vector (N⋅m)
         */
        Vector3d getAngular() const;

        /**
         * @brief Get the force component
         * @return 3D force vector (N)
         */
        Vector3d getLinear() const;

        /**
         * @brief Add two force vectors
         * @param other Force vector to add
         * @return Combined force vector
         * @details Represents the superposition of two force actions. The result describes
         *          the combined effect of both forces applied simultaneously.
         */
        ForceVector operator+(const ForceVector &other) const;

        /**
         * @brief Subtract two force vectors
         * @param other Force vector to subtract
         * @return Difference force vector
         * @details Represents the relative force between two force states
         */
        ForceVector operator-(const ForceVector &other) const;

        /**
         * @brief Scale a force vector
         * @param scalar Scaling factor
         * @return Scaled force vector
         * @details Useful for scaling forces or representing fractional loads
         */
        ForceVector operator*(double scalar) const;

        /**
         * @brief Compute the motion cross product (crm operation)
         * @param other Force vector to cross with
         * @return Resulting force vector
         * @details This operation is used in the transformation of forces when
         *          changing reference points or coordinate frames.
         */
        ForceVector crossMotion(const ForceVector &other) const;

        /**
         * @brief Compute the force cross product (crf operation)
         * @param other Force vector to cross with
         * @return Resulting force vector
         * @details This operation represents how one force affects another force
         *          through the rigid body kinematics.
         *          The result follows the formula: [τ1×τ2 + f1×f2; τ1×f2]
         */
        ForceVector crossForce(const ForceVector &other) const;

        /**
         * @brief Compute dot product with another force vector
         * @param other Force vector to dot with
         * @return Scalar result
         * @warning The dot product between two force vectors does not have a direct
         *          physical interpretation. Use dot product with motion vectors for
         *          computing power.
         */
        double dot(const ForceVector &other) const;

        /**
         * @brief Print the force vector components
         * @details Outputs the torque (τ) and force (f) components in a readable format
         */
        void print() const;
    };

    /** @brief Type alias for more concise notation */
    using fv = ForceVector;

} // namespace SpatialAlgebra

#endif