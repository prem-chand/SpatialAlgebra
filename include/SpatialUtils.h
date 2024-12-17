#ifndef SPATIAL_UTILS_H
#define SPATIAL_UTILS_H

/**
 * @file SpatialUtils.h
 * @brief Utility functions for spatial algebra operations
 * 
 * This file contains various utility functions used in spatial algebra computations,
 * including cross products, dot products, and matrix operations.
 */

#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"

namespace SpatialAlgebra {

/**
 * @brief Creates a 3x3 skew-symmetric matrix from a 3D vector
 * @param v Input 3D vector
 * @return Skew-symmetric matrix representation
 * @details The resulting matrix S satisfies S*x = v×x for any vector x
 */
inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) noexcept {
    Eigen::Matrix3d m;
    m << 0, -v(2), v(1),
         v(2), 0, -v(0),
         -v(1), v(0), 0;
    return m;
}

/**
 * @brief Computes the spatial dot product between two spatial vectors
 * @param v1 First spatial vector
 * @param v2 Second spatial vector
 * @return Scalar dot product result
 */
inline double dot(const SpatialVector& v1, const SpatialVector& v2) noexcept {
    return v1.getAngular().dot(v2.getAngular()) + 
           v1.getLinear().dot(v2.getLinear());
}

/**
 * @brief Computes the spatial dot product between two motion vectors
 * @param v1 First motion vector
 * @param v2 Second motion vector
 * @return Scalar dot product result
 */
inline double dot(const MotionVector& v1, const MotionVector& v2) noexcept {
    return dot(static_cast<const SpatialVector&>(v1), 
              static_cast<const SpatialVector&>(v2));
}

/**
 * @brief Computes the spatial dot product between two force vectors
 * @param v1 First force vector
 * @param v2 Second force vector
 * @return Scalar dot product result
 */
inline double dot(const ForceVector& v1, const ForceVector& v2) noexcept {
    return dot(static_cast<const SpatialVector&>(v1), 
              static_cast<const SpatialVector&>(v2));
}

inline double dot(const MotionVector& v1, const ForceVector& v2) noexcept {
    return dot(static_cast<const SpatialVector&>(v1), 
              static_cast<const SpatialVector&>(v2));
}

/**
 * @brief Computes the spatial cross product between a motion vector and a force vector
 * @param v1 Motion vector
 * @param v2 Force vector
 * @return Resulting force vector
 * @details Implements: crf(v1, v2) = [w1×f1 + v1×f2, w1×f1]
 */
inline ForceVector cross(const MotionVector& v1, const ForceVector& v2) noexcept {
    const Vector3d& w1 = v1.getAngular();
    const Vector3d& v1_lin = v1.getLinear();
    const Vector3d& f1 = v2.getAngular();
    const Vector3d& f2 = v2.getLinear();
    
    return ForceVector(
        w1.cross(f1), 
        w1.cross(f2) + v1_lin.cross(f1)
    );
};

} // namespace SpatialAlgebra
#endif // SPATIAL_UTILS_H