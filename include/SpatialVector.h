#ifndef SPATIAL_VECTOR_H
#define SPATIAL_VECTOR_H

#include <Eigen/Dense>

using Vector3d = Eigen::Matrix<double, 3, 1>;

namespace SpatialAlgebra
{
    /**
     * @brief Base class for spatial vectors in 6D space
     *
     * This class implements the basic operations for 6D spatial vectors used in
     * rigid body dynamics. It serves as a base class for both MotionVector and
     * ForceVector classes.
     */
    class SpatialVector
    {
    protected:
        Vector3d angular; ///< Angular component (ω for motion, τ for force)
        Vector3d linear;  ///< Linear component (v for motion, f for force)

    public:
        /** @brief Default constructor creating a zero spatial vector */
        SpatialVector();

        /**
         * @brief Construct a spatial vector from angular and linear components
         * @param angular The angular component (3D vector)
         * @param linear The linear component (3D vector)
         */
        SpatialVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Copy constructor
         * @param other The spatial vector to copy from
         */
        SpatialVector(const SpatialVector &other);

        /**
         * @brief Get the angular component
         * @return The angular 3D vector
         */
        Vector3d getAngular() const;

        /**
         * @brief Get the linear component
         * @return The linear 3D vector
         */
        Vector3d getLinear() const;

        /**
         * @brief Add two spatial vectors
         * @param other The spatial vector to add
         * @return The sum as a new spatial vector
         */
        SpatialVector operator+(const SpatialVector &other) const;

        /**
         * @brief Subtract two spatial vectors
         * @param other The spatial vector to subtract
         * @return The difference as a new spatial vector
         */
        SpatialVector operator-(const SpatialVector &other) const;

        /**
         * @brief Scalar multiplication
         * @param scalar The scalar value to multiply with
         * @return The scaled spatial vector
         */
        SpatialVector operator*(double scalar) const;

        /**
         * @brief Compute the motion cross product
         * @param other The spatial vector to cross with
         * @return The resulting spatial vector
         */
        SpatialVector crossMotion(const SpatialVector &other) const;

        /**
         * @brief Compute the force cross product
         * @param other The spatial vector to cross with
         * @return The resulting spatial vector
         */
        SpatialVector crossForce(const SpatialVector &other) const;

        /**
         * @brief Compute the dot product with another spatial vector
         * @param other The spatial vector to compute dot product with
         * @return The scalar dot product result
         */
        double dot(const SpatialVector &other) const;

        /**
         * @brief Print the spatial vector components
         */
        void print() const;
    };

} // namespace SpatialAlgebra

#endif
