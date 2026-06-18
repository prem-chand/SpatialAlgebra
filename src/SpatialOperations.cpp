/**
 * @file SpatialOperations.cpp
 * @brief Implementation of static utility class for spatial algebra operations
 */

#include "SpatialOperations.h"
#include "SpatialUtils.h"

namespace SpatialAlgebra {

SpatialVector SpatialOperations::crossProductMotion(const MotionVector& v1, const MotionVector& v2) {
    return cross(v1, v2);
}

SpatialVector SpatialOperations::crossProductForce(const MotionVector& v, const ForceVector& f) {
    return cross(v, f);
}

RigidBodyInertia SpatialOperations::transformInertia(const RigidBodyInertia& inertia, 
                                                     const PluckerTransform& transform) {
    // Delegate to PluckerTransform::tformRBI method
    // Formula: I' = X * I * X^T (transformed to new coordinate frame)
    return transform.tformRBI(inertia);
}

} // namespace SpatialAlgebra
