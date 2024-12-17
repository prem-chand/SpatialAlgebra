#ifndef FORCE_VECTOR_H
#define FORCE_VECTOR_H

/**
 * @file ForceVector.h
 * @brief Class representing spatial force vectors
 */

#include <iostream>

#include "SpatialVector.h"

// Force vector class inherits from SpatialVector
namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial force vectors
     * 
     * Specialized spatial vector for representing force quantities
     * like torques and linear forces.
     */
    class ForceVector : public SpatialVector
    {
    public:
        /** @brief Default constructor creating zero force */
        ForceVector();

        /**
         * @brief Construct from angular and linear components
         * @param angular Torque component
         * @param linear Force component
         */
        ForceVector(const Vector3d &angular, const Vector3d &linear);

        /**
         * @brief Construct from spatial vector
         * @param other Spatial vector to convert
         */
        ForceVector(const SpatialVector &other);

        // Accessors
        Vector3d getAngular() const;
        Vector3d getLinear() const;

        // Basic operations
        ForceVector operator+(const ForceVector &other) const;
        ForceVector operator-(const ForceVector &other) const;
        ForceVector operator*(double scalar) const;

        // Cross-product operations
        ForceVector crossMotion(const ForceVector &other) const; // crm equivalent
        ForceVector crossForce(const ForceVector &other) const;  // crf equivalent

        // dot product
        double dot(const ForceVector &other) const;

        // Printing
        void print() const;
    };

    using fv = ForceVector;

} // namespace SpatialAlgebra

#endif