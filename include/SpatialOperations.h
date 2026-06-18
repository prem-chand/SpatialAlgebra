#ifndef SPATIAL_OPERATIONS_H
#define SPATIAL_OPERATIONS_H

/**
 * @file SpatialOperations.h
 * @brief Static utility class for spatial algebra operations
 */

#include "MotionVector.h"
#include "ForceVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include <Eigen/Geometry>

namespace SpatialAlgebra {

class SpatialOperations {
public:
    static SpatialVector crossProductMotion(const MotionVector& v1, const MotionVector& v2);

    static SpatialVector crossProductForce(const MotionVector& v, const ForceVector& f);

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
