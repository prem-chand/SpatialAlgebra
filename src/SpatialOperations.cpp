/**
 * @file SpatialOperations.cpp
 * @brief Implementation of static utility class for spatial algebra operations
 */

#include "SpatialOperations.h"
#include "SpatialUtils.h"

namespace SpatialAlgebra {

SpatialVector SpatialOperations::crossProductMotion(const SpatialVector& v1, const SpatialVector& v2) {
    // Delegate to MotionVector cross product implementation
    // Cast to MotionVector since cross motion product is defined for motion vectors
    const MotionVector& mv1 = static_cast<const MotionVector&>(v1);
    const MotionVector& mv2 = static_cast<const MotionVector&>(v2);
    return cross(mv1, mv2);
}

SpatialVector SpatialOperations::crossProductForce(const SpatialVector& v, const SpatialVector& f) {
    // Delegate to ForceVector cross product implementation
    // Cast to appropriate types: motion vector cross force vector
    const MotionVector& mv = static_cast<const MotionVector&>(v);
    const ForceVector& fv = static_cast<const ForceVector&>(f);
    return cross(mv, fv);
}

RigidBodyInertia SpatialOperations::transformInertia(const RigidBodyInertia& inertia, 
                                                     const PluckerTransform& transform) {
    // Delegate to PluckerTransform::tformRBI method
    // Formula: I' = X * I * X^T (transformed to new coordinate frame)
    return transform.tformRBI(inertia);
}

} // namespace SpatialAlgebra
