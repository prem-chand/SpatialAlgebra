#ifndef SPATIAL_OPERATIONS_H
#define SPATIAL_OPERATIONS_H

/**
 * @file SpatialOperations.h
 * @brief Static utility class for spatial algebra operations
 */

#include "SpatialVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include <Eigen/Geometry> // Include for RotationMatrix

namespace SpatialAlgebra {

/**
 * @brief Static class providing utility operations for spatial vectors and transforms
 */
class SpatialOperations {
public:
    /**
     * @brief Computes the motion cross product between two spatial vectors
     * @param v1 First spatial vector
     * @param v2 Second spatial vector
     * @return Result of motion cross product
     */
    static SpatialVector crossProductMotion(const SpatialVector& v1, const SpatialVector& v2);

    /**
     * @brief Computes the force cross product between two spatial vectors
     * @param v Force spatial vector
     * @param f Motion spatial vector
     * @return Result of force cross product
     */
    static SpatialVector crossProductForce(const SpatialVector& v, const SpatialVector& f);

    /**
     * @brief Transforms rigid body inertia using a Plucker transform
     * @param inertia Input rigid body inertia
     * @param transform Plucker transform to apply
     * @return Transformed rigid body inertia
     */
    static RigidBodyInertia transformInertia(const RigidBodyInertia& inertia, 
                                            const PluckerTransform& transform);
};

}  // namespace SpatialAlgebra

#endif
