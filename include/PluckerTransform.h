#ifndef PLUCKER_TRANSFORM_H
#define PLUCKER_TRANSFORM_H

/**
 * @file PluckerTransform.h
 * @brief Class representing coordinate transformations in Plucker coordinates
 */

#include "SpatialVector.h"
#include "Rotation.h"
#include "ForceVector.h"
#include "MotionVector.h"
#include "RigidBodyInertia.h"   
#include "ArticulatedBodyInertia.h"
#include <Eigen/Geometry> // Include for RotationMatrix

using Vector3d = Eigen::Matrix<double, 3, 1>;

namespace SpatialAlgebra
{
    /**
     * @brief Class representing spatial coordinate transformations
     * 
     * Implements rigid body transformations in Plucker coordinates,
     * consisting of rotation and translation components.
     */
    class PluckerTransform
    {
    private:
        Rotation rotation;    ///< 3x3 Rotation matrix
        Vector3d translation; ///< 3D Translation vector

    public:
        /**
         * @brief Construct a Plucker transform
         * @param rotation Rotation component
         * @param translation Translation component
         */
        PluckerTransform(const Rotation &rotation, const Vector3d &translation);

        /**
         * @brief Transform a motion vector
         * @param vec Input motion vector
         * @return Transformed motion vector
         */
        SpatialVector transformMotion(const SpatialVector &vec) const;

        /**
         * @brief Transform a force vector
         * @param vec Input force vector
         * @return Transformed force vector
         */
        SpatialVector transformForce(const SpatialVector &vec) const;

        /**
         * @brief Inverse transform a motion vector
         * @param vec Input motion vector
         * @return Inverse transformed motion vector
         */
        SpatialVector inverseTransformMotion(const SpatialVector &vec) const;

        /**
         * @brief Inverse transform a force vector
         * @param vec Input force vector
         * @return Inverse transformed force vector
         */
        SpatialVector inverseTransformForce(const SpatialVector &vec) const;

        /**
         * @brief Transform rigid body inertia
         * @param Ihat Input rigid body inertia
         * @return Transformed inertia
         */
        RigidBodyInertia tformRBI(const RigidBodyInertia &Ihat) const;

        /**
         * @brief Inverse transform rigid body inertia
         * @param Ihat Input rigid body inertia
         * @return Inverse transformed inertia
         */
        RigidBodyInertia invtformRBI(const RigidBodyInertia &Ihat) const;

        // transform abi
        ArticulatedBodyInertia tformABI(const ArticulatedBodyInertia &Ia) const;
        ArticulatedBodyInertia invtformABI(const ArticulatedBodyInertia &Ia) const;
        
        /** @brief Compute inverse transform */
        PluckerTransform inverse() const;

        // multiply by X
        // PluckerTransform multiply(const PluckerTransform &X) const;
        auto multiply(const PluckerTransform &X) const;

        auto apply(const fv& v) const{
            return transformMotion(v);
        }
        auto apply(const mv& v) const{
            return transformForce(v);
        }

        PluckerTransform apply(const PluckerTransform &X) const;

        // template <typename T>
        // T apply(const T&) const{

        // }

        /** @brief Print transform components */
        void print() const;

        // testing
        auto get(int &a) { return a; }
        auto get(double &b)
        {
            char c{'c'};
            return c;
        }
    };

    using plux = PluckerTransform;

} // namespace SpatialAlgebra

#endif
