#ifndef MOTION_VECTOR_H
#define MOTION_VECTOR_H

/**
 * @file MotionVector.h
 * @brief Class representing spatial motion vectors
 */

#include <array>
#include <iostream>

#include "SpatialVector.h"

// Motion vector class inherits from SpatialVector

namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial motion vectors
     * 
     * Specialized spatial vector for representing motion quantities
     * like angular and linear velocities.
     */
    class MotionVector : public SpatialVector
    {
    public:
        /** @brief Default constructor creating zero motion */
        MotionVector();

        /**
         * @brief Construct from angular and linear components
         * @param angular Angular velocity component
         * @param linear Linear velocity component
         */
        MotionVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Construct from spatial vector
         * @param other Spatial vector to convert
         */
        MotionVector(const SpatialVector &other);

        // Accessors
        Vector3d getAngular() const;
        Vector3d getLinear() const;

        // Basic operations
        MotionVector operator+(const MotionVector &other) const;
        MotionVector operator-(const MotionVector &other) const;
        MotionVector operator*(double scalar) const;

        // Cross-product operations
        MotionVector crossMotion(const MotionVector &other) const; // crm equivalent
        MotionVector crossForce(const MotionVector &other) const;  // crf equivalent

        // dot product
        double dot(const MotionVector &other) const;

        // Printing
        void print() const;
    };

    using mv = MotionVector;

} // namespace SpatialAlgebra

#endif