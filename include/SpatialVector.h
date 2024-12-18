#ifndef SPATIAL_VECTOR_H
#define SPATIAL_VECTOR_H

/**
 * @file SpatialVector.h
 * @brief Implementation of 6D spatial vectors for rigid body dynamics
 * @details This file provides the foundation for spatial algebra operations used in rigid body
 *          dynamics. Spatial vectors are 6D vectors that combine angular and linear components
 *          to represent motion (twists) and forces (wrenches) in a unified framework.
 * 
 *          Mathematical Background:
 *          A spatial vector ŝ consists of two 3D vectors: [ω; v]
 *          - For motion vectors (twists): ω represents angular velocity, v represents linear velocity
 *          - For force vectors (wrenches): ω represents torque/moment, v represents linear force
 * 
 *          Key Operations:
 *          - Addition/Subtraction: Component-wise operations
 *          - Scalar Multiplication: Scales both components
 *          - Cross Products: Special 6D operations preserving motion/force duality
 *          - Dot Product: Inner product between motion and force vectors
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create a spatial vector
 *     Vector3d angular(0, 0, 1);  // Z-axis rotation
 *     Vector3d linear(1, 0, 0);   // X-axis translation
 *     SpatialVector vec1(angular, linear);
 *     
 *     // Perform operations
 *     SpatialVector vec2(Vector3d(1,0,0), Vector3d(0,1,0));
 *     SpatialVector sum = vec1 + vec2;
 *     SpatialVector crossed = vec1.crossMotion(vec2);
 *     double dotProduct = vec1.dot(vec2);
 * @endcode
 * 
 * @note This implementation follows the spatial algebra formulation by Roy Featherstone
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Springer.
 */

#include <Eigen/Dense>

using Vector3d = Eigen::Matrix<double, 3, 1>;

namespace SpatialAlgebra
{
    /**
     * @brief Base class for spatial vectors in 6D space
     * @details The SpatialVector class implements the fundamental operations for 6D spatial
     *          vectors used in rigid body dynamics. It serves as the base class for both
     *          MotionVector (twists) and ForceVector (wrenches) classes.
     * 
     *          Mathematical Properties:
     *          - Dimension: 6 (3 for angular + 3 for linear)
     *          - Addition forms a vector space
     *          - Cross products follow special rules for motion/force duality
     * 
     *          Implementation Notes:
     *          - Uses Eigen::Vector3d for efficient 3D vector operations
     *          - Stores angular and linear components separately for clarity
     *          - Provides both motion and force cross products
     * 
     * @warning The distinction between motion and force vectors is crucial for correct
     *          physical interpretation and operations. Use the appropriate derived class
     *          (MotionVector or ForceVector) for type safety.
     */
    class SpatialVector
    {
    protected:
        Vector3d angular; ///< Angular component (ω for motion vectors, τ for force vectors)
        Vector3d linear;  ///< Linear component (v for motion vectors, f for force vectors)

    public:
        /** 
         * @brief Default constructor creating a zero spatial vector
         * @details Initializes both angular and linear components to zero vectors
         */
        SpatialVector();

        /**
         * @brief Construct a spatial vector from angular and linear components
         * @param angular The angular component (3D vector)
         * @param linear The linear component (3D vector)
         * @details Creates a spatial vector with specified components. The interpretation
         *          depends on whether it's used as a motion or force vector:
         *          - Motion vector: angular = ω (angular velocity), linear = v (linear velocity)
         *          - Force vector: angular = τ (torque), linear = f (force)
         */
        SpatialVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Copy constructor
         * @param other The spatial vector to copy from
         * @details Performs a deep copy of both angular and linear components
         */
        SpatialVector(const SpatialVector &other);

        /**
         * @brief Get the angular component
         * @return The angular 3D vector
         * @details Returns the top 3D vector component (ω for motion, τ for force)
         */
        Vector3d getAngular() const;

        /**
         * @brief Get the linear component
         * @return The linear 3D vector
         * @details Returns the bottom 3D vector component (v for motion, f for force)
         */
        Vector3d getLinear() const;

        /**
         * @brief Add two spatial vectors
         * @param other The spatial vector to add
         * @return The sum as a new spatial vector
         * @details Performs component-wise addition: [ω1+ω2; v1+v2]
         */
        SpatialVector operator+(const SpatialVector &other) const;

        /**
         * @brief Subtract two spatial vectors
         * @param other The spatial vector to subtract
         * @return The difference as a new spatial vector
         * @details Performs component-wise subtraction: [ω1-ω2; v1-v2]
         */
        SpatialVector operator-(const SpatialVector &other) const;

        /**
         * @brief Scalar multiplication
         * @param scalar The scalar value to multiply with
         * @return The scaled spatial vector
         * @details Scales both components: [k*ω; k*v]
         */
        SpatialVector operator*(double scalar) const;

        /**
         * @brief Compute the motion cross product
         * @param other The spatial vector to cross with
         * @return The resulting spatial vector
         * @details Implements the motion cross product operation:
         *          [ω1×ω2; ω1×v2 + v1×ω2]
         * @note This operation is specific to motion vectors (twists)
         */
        SpatialVector crossMotion(const SpatialVector &other) const;

        /**
         * @brief Compute the force cross product
         * @param other The spatial vector to cross with
         * @return The resulting spatial vector
         * @details Implements the force cross product operation:
         *          [ω1×ω2 + v1×v2; ω1×v2]
         * @note This operation is specific to force vectors (wrenches)
         */
        SpatialVector crossForce(const SpatialVector &other) const;

        /**
         * @brief Compute the dot product with another spatial vector
         * @param other The spatial vector to compute dot product with
         * @return The scalar dot product result
         * @details Computes the spatial dot product:
         *          ω1·ω2 + v1·v2
         * @note The dot product between a motion and force vector has physical meaning
         *       (power), while dot products between the same type do not
         */
        double dot(const SpatialVector &other) const;

        /**
         * @brief Print the spatial vector components
         * @details Outputs the angular and linear components in a readable format
         */
        void print() const;
    };

} // namespace SpatialAlgebra

#endif // SPATIAL_VECTOR_H
