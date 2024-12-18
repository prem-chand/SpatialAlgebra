#ifndef PLUCKER_TRANSFORM_H
#define PLUCKER_TRANSFORM_H

/**
 * @file PluckerTransform.h
 * @brief Class representing coordinate transformations in Plücker coordinates
 * @details This file implements rigid body coordinate transformations using Plücker
 *          coordinates, which are essential for spatial algebra operations in robotics
 *          and multibody dynamics. These transforms combine rotation and translation
 *          in a way that preserves the dual relationship between motion and force vectors.
 * 
 *          Mathematical Representation:
 *          A Plücker transform X = [E  0; -E*r̂  E] where:
 *          - E: 3x3 rotation matrix
 *          - r̂: 3x3 skew-symmetric matrix from translation vector r
 *          - [A; B] denotes a 6x6 block matrix
 * 
 *          Key Properties:
 *          - Preserves the duality between motion and force vectors
 *          - Forms a Lie group (SE(3)) with composition as the group operation
 *          - Transforms spatial quantities consistently with rigid body motion
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     // Create a transform: 90° rotation around Z + translation along X
 *     Rotation rot = Rotation(Eigen::AngleAxisd(M_PI/2, Vector3d::UnitZ()));
 *     Vector3d trans(1, 0, 0);
 *     PluckerTransform X(rot, trans);
 *     
 *     // Transform a motion vector
 *     MotionVector vel(Vector3d(0,0,1), Vector3d(1,0,0));
 *     MotionVector transformed = X.transformMotion(vel);
 *     
 *     // Transform a force vector
 *     ForceVector force(Vector3d(0,0,0), Vector3d(1,0,0));
 *     ForceVector transformed_force = X.transformForce(force);
 * @endcode
 * 
 * @see MotionVector for the representation of spatial velocities
 * @see ForceVector for the representation of spatial forces
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2.3
 */

#include "SpatialVector.h"
#include "Rotation.h"
#include "ForceVector.h"
#include "MotionVector.h"
#include "RigidBodyInertia.h"   
#include "ArticulatedBodyInertia.h"
#include <Eigen/Geometry>

using Vector3d = Eigen::Matrix<double, 3, 1>;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial coordinate transformations
     * @details The PluckerTransform class implements rigid body transformations in
     *          Plücker coordinates, providing a consistent framework for transforming
     *          spatial vectors and inertias between different coordinate frames.
     * 
     *          Mathematical Properties:
     *          - Composition rule: X1 * X2 represents sequential transforms
     *          - Inverse transform: X^(-1) exists and X * X^(-1) = identity
     *          - Motion transform: v' = Xm * v
     *          - Force transform: f' = Xf * f = X^(-T) * f
     * 
     *          Implementation Notes:
     *          - Stores rotation and translation separately for efficiency
     *          - Provides both motion and force transformation methods
     *          - Handles both direct and inverse transformations
     *          - Supports transformation of inertia matrices
     * 
     * @warning The translation vector represents the origin of the new frame
     *          expressed in the old frame coordinates.
     */
    class PluckerTransform
    {
    private:
        Rotation rotation;    ///< 3x3 Rotation matrix component of the transform
        Vector3d translation; ///< 3D Translation vector from old to new origin

    public:
        /**
         * @brief Construct a Plücker transform
         * @param rotation Rotation component (SO(3) matrix)
         * @param translation Translation vector (3D)
         * @details Creates a spatial transform that rotates and translates spatial vectors.
         *          The translation is applied after the rotation in the rotated frame.
         */
        PluckerTransform(const Rotation &rotation, const Vector3d &translation);

        /**
         * @brief Transform a motion vector
         * @param vec Input motion vector
         * @return Transformed motion vector
         * @details Applies the spatial transformation to a motion vector:
         *          v' = [E    0; -E*r̂  E] * [ω; v]
         */
        SpatialVector transformMotion(const SpatialVector &vec) const;

        /**
         * @brief Transform a force vector
         * @param vec Input force vector
         * @return Transformed force vector
         * @details Applies the spatial transformation to a force vector:
         *          f' = [E    E*r̂; 0    E] * [n; f]
         */
        SpatialVector transformForce(const SpatialVector &vec) const;

        /**
         * @brief Inverse transform a motion vector
         * @param vec Input motion vector
         * @return Inverse transformed motion vector
         * @details Applies the inverse spatial transformation: X^(-1) * v
         */
        SpatialVector inverseTransformMotion(const SpatialVector &vec) const;

        /**
         * @brief Inverse transform a force vector
         * @param vec Input force vector
         * @return Inverse transformed force vector
         * @details Applies the inverse spatial transformation: X^(-T) * f
         */
        SpatialVector inverseTransformForce(const SpatialVector &vec) const;

        /**
         * @brief Transform rigid body inertia
         * @param Ihat Input rigid body inertia
         * @return Transformed inertia
         * @details Transforms a rigid body inertia to a new coordinate frame:
         *          I' = X * I * X^T
         */
        RigidBodyInertia tformRBI(const RigidBodyInertia &Ihat) const;

        /**
         * @brief Inverse transform rigid body inertia
         * @param Ihat Input rigid body inertia
         * @return Inverse transformed inertia
         * @details Applies the inverse transformation: X^(-1) * I * X^(-T)
         */
        RigidBodyInertia invtformRBI(const RigidBodyInertia &Ihat) const;

        /**
         * @brief Transform articulated body inertia
         * @param Ia Input articulated body inertia
         * @return Transformed articulated body inertia
         * @details Transforms an articulated body inertia to a new coordinate frame
         */
        ArticulatedBodyInertia tformABI(const ArticulatedBodyInertia &Ia) const;

        /**
         * @brief Inverse transform articulated body inertia
         * @param Ia Input articulated body inertia
         * @return Inverse transformed articulated body inertia
         */
        ArticulatedBodyInertia invtformABI(const ArticulatedBodyInertia &Ia) const;
        
        /**
         * @brief Compute inverse transform
         * @return Inverse Plücker transform
         * @details Returns X^(-1) such that X * X^(-1) = identity
         */
        PluckerTransform inverse() const;

        /**
         * @brief Compose with another transform
         * @param X Right-hand side transform
         * @return Combined transform
         * @details Computes the sequential application of transforms: this * X
         */
        auto multiply(const PluckerTransform &X) const;

        /**
         * @brief Apply transform to a force vector
         * @param v Input force vector
         * @return Transformed force vector
         */
        auto apply(const fv& v) const {
            return transformForce(v);
        }

        /**
         * @brief Apply transform to a motion vector
         * @param v Input motion vector
         * @return Transformed motion vector
         */
        auto apply(const mv& v) const {
            return transformMotion(v);
        }

        /**
         * @brief Compose with another transform
         * @param X Right-hand side transform
         * @return Combined transform
         */
        PluckerTransform apply(const PluckerTransform &X) const;

        /**
         * @brief Print transform components
         * @details Outputs the rotation matrix and translation vector
         */
        void print() const;
    };

    /** @brief Type alias for more concise notation */
    using plux = PluckerTransform;

} // namespace SpatialAlgebra

#endif // PLUCKER_TRANSFORM_H
