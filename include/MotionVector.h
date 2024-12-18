#ifndef MOTION_VECTOR_H
#define MOTION_VECTOR_H

/**
 * @file MotionVector.h
 * @brief Class representing spatial motion vectors (twists) in rigid body dynamics
 * @details This file implements motion vectors (also known as twists) in spatial algebra,
 *          which unify the representation of angular and linear velocities in 3D space.
 *          Motion vectors are essential for describing the instantaneous motion state
 *          of rigid bodies in a compact and geometrically meaningful way.
 * 
 *          Mathematical Representation:
 *          A motion vector m = [ω; v] consists of:
 *          - ω: Angular velocity vector (rad/s)
 *          - v: Linear velocity vector (m/s)
 *          
 *          Physical Interpretation:
 *          - ω represents the instantaneous angular velocity around an axis
 *          - v represents the linear velocity of the body's origin
 *          - Together they form a "twist" that completely describes rigid body motion
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create a pure rotation around Z-axis (1 rad/s)
 *     MotionVector angular_motion(Vector3d(0,0,1), Vector3d(0,0,0));
 *     
 *     // Create a pure translation along X-axis (1 m/s)
 *     MotionVector linear_motion(Vector3d(0,0,0), Vector3d(1,0,0));
 *     
 *     // Combine motions
 *     MotionVector combined = angular_motion + linear_motion;
 *     
 *     // Compute motion cross product (for velocity propagation)
 *     MotionVector crossed = angular_motion.crossMotion(linear_motion);
 * @endcode
 * 
 * @see ForceVector for the dual space representation of forces and torques
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2.
 */

#include <array>
#include <iostream>
#include "SpatialVector.h"

namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial motion vectors (twists)
     * @details The MotionVector class specializes SpatialVector for representing
     *          motion quantities in spatial algebra. It provides type-safe operations
     *          specific to motion vectors and ensures correct geometric interpretation.
     * 
     *          Key Concepts:
     *          - Motion vectors form a vector space
     *          - They transform covariantly under coordinate changes
     *          - The cross product represents velocity propagation in rigid bodies
     *          - Dot product with force vectors yields power
     * 
     *          Common Applications:
     *          - Representing joint velocities in articulated systems
     *          - Describing rigid body motion in space
     *          - Velocity propagation through kinematic chains
     * 
     * @note Motion vectors and force vectors are dual to each other, meaning they
     *       transform differently under coordinate changes and their dot product
     *       yields physically meaningful power.
     */
    class MotionVector : public SpatialVector
    {
    public:
        /** 
         * @brief Default constructor creating zero motion
         * @details Initializes a motion vector with zero angular and linear velocity
         */
        MotionVector();

        /**
         * @brief Construct from angular and linear components
         * @param angular Angular velocity component (rad/s)
         * @param linear Linear velocity component (m/s)
         * @details Creates a motion vector representing a combined rotational and
         *          translational motion. The angular component represents rotation
         *          around an axis, while the linear component represents translation.
         */
        MotionVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Construct from spatial vector
         * @param other Spatial vector to convert
         * @details Converts a generic spatial vector to a motion vector. This is useful
         *          when working with algorithms that operate on generic spatial vectors.
         */
        MotionVector(const SpatialVector &other);

        /**
         * @brief Get the angular velocity component
         * @return 3D angular velocity vector (rad/s)
         */
        Vector3d getAngular() const;

        /**
         * @brief Get the linear velocity component
         * @return 3D linear velocity vector (m/s)
         */
        Vector3d getLinear() const;

        /**
         * @brief Add two motion vectors
         * @param other Motion vector to add
         * @return Combined motion vector
         * @details Represents the superposition of two motions. The result describes
         *          the combined effect of both motions applied simultaneously.
         */
        MotionVector operator+(const MotionVector &other) const;

        /**
         * @brief Subtract two motion vectors
         * @param other Motion vector to subtract
         * @return Difference motion vector
         * @details Represents the relative motion between two motion states
         */
        MotionVector operator-(const MotionVector &other) const;

        /**
         * @brief Scale a motion vector
         * @param scalar Scaling factor
         * @return Scaled motion vector
         * @details Useful for scaling velocities or representing fractional motions
         */
        MotionVector operator*(double scalar) const;

        /**
         * @brief Compute the motion cross product (crm operation)
         * @param other Motion vector to cross with
         * @return Resulting motion vector
         * @details This operation is crucial for velocity propagation in rigid body
         *          systems. It represents how one motion affects another motion.
         *          The result follows the formula: [ω1×ω2; ω1×v2 + v1×ω2]
         */
        MotionVector crossMotion(const MotionVector &other) const;

        /**
         * @brief Compute the force cross product (crf operation)
         * @param other Motion vector to cross with
         * @return Resulting motion vector
         * @details This operation is the transpose of the adjoint map and is used
         *          in force transformation calculations.
         */
        MotionVector crossForce(const MotionVector &other) const;

        /**
         * @brief Compute dot product with another motion vector
         * @param other Motion vector to dot with
         * @return Scalar result
         * @warning The dot product between two motion vectors does not have a direct
         *          physical interpretation. Use dot product with force vectors for
         *          computing power.
         */
        double dot(const MotionVector &other) const;

        /**
         * @brief Print the motion vector components
         * @details Outputs the angular (ω) and linear (v) components in a readable format
         */
        void print() const;
    };

    /** @brief Type alias for more concise notation */
    using mv = MotionVector;

} // namespace SpatialAlgebra

#endif // MOTION_VECTOR_H