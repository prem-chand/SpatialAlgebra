#ifndef INVERSE_DYNAMICS_H
#define INVERSE_DYNAMICS_H

/**
 * @file InverseDynamics.h
 * @brief Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics computation
 * @details This file implements Featherstone's Recursive Newton-Euler Algorithm (RNEA)
 *          for computing joint torques from joint motion. The RNEA is an O(n) recursive
 *          algorithm that efficiently solves the inverse dynamics problem for serial
 *          and branching kinematic chains.
 * 
 *          Algorithm Overview:
 *          1. Outward pass (base to tip): Propagate velocities and accelerations
 *          2. Inward pass (tip to base): Propagate forces, compute joint torques
 * 
 *          Mathematical Foundation:
 *          The RNEA solves the equation: τ = H(q)q̈ + C(q,q̇) + G(q)
 *          where:
 *          - H(q): Joint space inertia matrix
 *          - C(q,q̇): Coriolis and centrifugal forces
 *          - G(q): Gravity forces (not implemented, add external forces)
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     InverseDynamics id;
 *     
 *     // Setup kinematic chain
 *     Link link;
 *     link.parent = -1;  // Base link
 *     link.X = PluckerTransform(Rotation::Identity(), Vector3d::Zero());
 *     link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
 *     link.S = MotionVector(Vector3d(0,0,1), Vector3d::Zero());  // Revolute Z
 *     link.q = 0.0;
 *     link.qdot = 0.0;
 *     link.qddot = 1.0;
 *     id.links.push_back(link);
 *     
 *     // Compute torques
 *     Eigen::VectorXd qddot(1);
 *     qddot[0] = 1.0;
 *     
 *     Eigen::VectorXd tau = id.computeTorques(qddot);
 * @endcode
 * 
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7
 * @see PluckerTransform for coordinate transformations
 */

#include <vector>
#include <Eigen/Dense>
#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"

namespace SpatialAlgebra
{
    /**
     * @brief Represents a single link in a kinematic chain for inverse dynamics
     * @details The InverseDynamicsLink structure contains all necessary information for RNEA:
     *          - Kinematic properties (parent, transform, joint axis)
     *          - Inertial properties (rigid body inertia)
     *          - State variables (position, velocity, acceleration)
     *          - Intermediate quantities for RNEA (spatial velocity, spatial acceleration)
     * 
     *          Memory Layout:
     *          All spatial vectors use 6D representation: [angular; linear]
     *          Transforms are stored as PluckerTransform (rotation + translation)
     * 
     * @note Parent index -1 indicates the base link (no parent)
     */
    struct InverseDynamicsLink
    {
        int parent;                     ///< Parent link index (-1 for base)
        PluckerTransform X;             ///< Transform from parent to this link
        RigidBodyInertia I;             ///< Rigid body inertia
        MotionVector S;                 ///< Joint motion axis (screw axis)
        
        double q;                       ///< Joint position
        double qdot;                    ///< Joint velocity
        double qddot;                   ///< Joint acceleration (input)
        
        MotionVector v;                 ///< Spatial velocity (computed in outward pass)
        MotionVector a;                 ///< Spatial acceleration (computed in outward pass)
        
        /**
         * @brief Default constructor
         * @details Initializes all fields to zero/identity values
         */
        InverseDynamicsLink() : parent(-1), 
                 X(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero()),
                 I(1.0, Vector3d::Zero(), lt::Identity(3)),
                 S(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
                 q(0.0), qdot(0.0), qddot(0.0),
                 v(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
                 a(MotionVector(Vector3d::Zero(), Vector3d::Zero())) {}
    };

    /**
     * @brief Inverse dynamics solver using Recursive Newton-Euler Algorithm
     * @details The InverseDynamics class implements Featherstone's RNEA (Algorithm 7.1)
     *          for efficient O(n) computation of joint torques. The algorithm
     *          handles both serial chains and branching kinematic trees.
     * 
     *          Algorithm Complexity:
     *          - Time: O(n) where n is the number of links
     *          - Space: O(n) for storing link states
     * 
     *          Usage Pattern:
     *          1. Setup kinematic tree by populating links vector
     *          2. Set joint states (q, qdot, qddot)
     *          3. Call computeTorques(qddot) to solve for joint torques
     * 
     * @note Links must be ordered such that parents appear before children
     * @note For branching trees, all children of a link must have consecutive indices
     */
    class InverseDynamics
    {
    public:
        std::vector<InverseDynamicsLink> links;  ///< Kinematic tree (ordered: parents before children)

        /**
         * @brief Compute joint torques from joint accelerations
         * @param qddot Vector of joint accelerations (must match links.size())
         * @return Vector of joint torques
         * @details Main entry point for inverse dynamics computation.
         *          Executes the two-phase RNEA:
         *          1. Outward pass: propagate velocities and accelerations
         *          2. Inward pass: propagate forces, compute joint torques
         * 
         *          Mathematical formulation:
         *          Outward: vᵢ = Xᵢ⁻¹·v_parent + Sᵢ·q̇ᵢ
         *                   aᵢ = Xᵢ⁻¹·a_parent + Sᵢ·q̈ᵢ + vᵢ × Sᵢ·q̇ᵢ
         *          Inward:  fᵢ = Iᵢ·aᵢ + vᵢ × Iᵢ·vᵢ
         *                   τᵢ = fᵢ·Sᵢ
         * 
         * @param gravity Gravity vector (default zero) for gravity-aware dynamics
         * @throws std::invalid_argument if qddot.size() != links.size()
         * @throws std::invalid_argument if qddot contains NaN or Inf values
         */
        Eigen::VectorXd computeTorques(const Eigen::VectorXd& qddot, const Vector3d& gravity = Vector3d::Zero());

    private:
        Vector3d gravity;  ///< Gravity vector (default zero, set by computeTorques)

        /**
         * @brief Outward pass: propagate velocities and accelerations
         * @details Iterates from base (index 0) to tip (index n-1).
         *          For each link:
         *          - Transform parent velocity to current frame
         *          - Add joint velocity contribution: v = X⁻¹·v_parent + S·q̇
         *          - Compute acceleration: a = X⁻¹·a_parent + S·q̈ + v × S·q̇
         * 
         *          Base link special case:
         *          - v₀ = S₀·q̇₀ (no parent velocity to transform)
         *          - a₀ = S₀·q̈₀ (no parent acceleration or velocity products)
         */
        void outwardPass();

        /**
         * @brief Inward pass: propagate forces and compute joint torques
         * @return Vector of joint torques
         * @details Iterates from tip (index n-1) to base (index 0).
         *          For each link:
         *          - Compute spatial force: f = I·a + v × I·v
         *          - Add transformed child forces if present
         *          - Project force onto joint axis: τ = f·S
         * 
         * @note Children must be processed before their parent
         */
        Eigen::VectorXd inwardPass();
    };
}

#endif // INVERSE_DYNAMICS_H
