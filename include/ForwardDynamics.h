#ifndef FORWARD_DYNAMICS_H
#define FORWARD_DYNAMICS_H

/**
 * @file ForwardDynamics.h
 * @brief Articulated Body Algorithm (ABA) for forward dynamics computation
 * @details This file implements Featherstone's Articulated Body Algorithm (ABA)
 *          for computing joint accelerations from applied torques. The ABA is
 *          an O(n) recursive algorithm that efficiently solves the forward
 *          dynamics problem for serial and branching kinematic chains.
 * 
 *          Algorithm Overview:
 *          1. Outward pass (base to tip): Propagate velocities, compute bias accelerations
 *          2. Inward pass (tip to base): Accumulate articulated inertias, solve for accelerations
 * 
 *          Mathematical Foundation:
 *          The ABA solves the equation: τ = H(q)q̈ + C(q,q̇) + G(q)
 *          where:
 *          - H(q): Joint space inertia matrix
 *          - C(q,q̇): Coriolis and centrifugal forces
 *          - G(q): Gravity forces (not implemented, add external forces)
 * 
 * Example usage:
 * @code{.cpp}
 *     using namespace SpatialAlgebra;
 *     
 *     ForwardDynamics fd;
 *     
 *     // Setup kinematic chain
 *     Link link;
 *     link.parent = -1;  // Base link
 *     link.X = PluckerTransform(Rotation::Identity(), Vector3d::Zero());
 *     link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
 *     link.S = MotionVector(Vector3d(0,0,1), Vector3d::Zero());  // Revolute Z
 *     link.q = 0.0;
 *     link.qdot = 0.0;
 *     fd.links.push_back(link);
 *     
 *     // Apply torques
 *     Eigen::VectorXd tau(1);
 *     tau[0] = 1.0;
 *     
 *     // Compute accelerations
 *     fd.computeAccelerations(tau);
 *     
 *     // Access results
 *     double qddot = fd.links[0].qddot;
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
#include "ArticulatedBodyInertia.h"

using Vector3d = Eigen::Matrix<double, 3, 1>;
using Matrix3d = Eigen::Matrix<double, 3, 3>;
using lt = LowerTriangular;

namespace SpatialAlgebra
{
    /**
     * @brief Represents a single link in a kinematic chain
     * @details The Link structure contains all necessary information for ABA:
     *          - Kinematic properties (parent, transform, joint axis)
     *          - Inertial properties (rigid body inertia)
     *          - State variables (position, velocity, acceleration)
     *          - Intermediate quantities for ABA (velocities, forces, articulated inertia)
     * 
     *          Memory Layout:
     *          All spatial vectors use 6D representation: [angular; linear]
     *          Transforms are stored as PluckerTransform (rotation + translation)
     * 
     * @note Parent index -1 indicates the base link (no parent)
     */
    struct Link
    {
        int parent;                     ///< Parent link index (-1 for base)
        PluckerTransform X;             ///< Transform from parent to this link
        RigidBodyInertia I;             ///< Rigid body inertia
        MotionVector S;                 ///< Joint motion axis (screw axis)
        
        double q;                       ///< Joint position
        double qdot;                    ///< Joint velocity
        double qddot;                   ///< Joint acceleration (output)
        
        MotionVector v;                 ///< Spatial velocity (computed in outward pass)
        MotionVector c;                 ///< Bias acceleration (Coriolis/centrifugal)
        ForceVector f;                  ///< Spatial force (external forces)
        
        ArticulatedBodyInertia Ia;      ///< Articulated body inertia (computed in inward pass)
        ForceVector pa;                 ///< Bias force (computed in inward pass)
        
        /**
         * @brief Default constructor
         * @details Initializes all fields to zero/identity values
         */
        Link() : parent(-1), 
                 X(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero()),
                 I(1.0, Vector3d::Zero(), lt::Identity(3)),
                 S(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
                 q(0.0), qdot(0.0), qddot(0.0),
                 v(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
                 c(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
                 f(ForceVector(Vector3d::Zero(), Vector3d::Zero())),
                 Ia(lt::Identity(3), Eigen::Matrix3d::Zero(), lt::Identity(3)),
                 pa(ForceVector(Vector3d::Zero(), Vector3d::Zero())) {}
    };

    /**
     * @brief Forward dynamics solver using Articulated Body Algorithm
     * @details The ForwardDynamics class implements Featherstone's ABA (Algorithm 7.3)
     *          for efficient O(n) computation of joint accelerations. The algorithm
     *          handles both serial chains and branching kinematic trees.
     * 
     *          Algorithm Complexity:
     *          - Time: O(n) where n is the number of links
     *          - Space: O(n) for storing link states
     * 
     *          Usage Pattern:
     *          1. Setup kinematic tree by populating links vector
     *          2. Set joint states (q, qdot) and external forces (f)
     *          3. Call computeAccelerations(tau) to solve for qddot
     * 
     * @note Links must be ordered such that parents appear before children
     * @note For branching trees, all children of a link must have consecutive indices
     */
    class ForwardDynamics
    {
    public:
        std::vector<Link> links;  ///< Kinematic tree (ordered: parents before children)

        /**
         * @brief Compute joint accelerations from applied torques
         * @param tau Vector of joint torques (must match links.size())
         * @details Main entry point for forward dynamics computation.
         *          Executes the three-phase ABA:
         *          1. Outward pass: propagate velocities, compute bias accelerations
         *          2. Inward pass: accumulate articulated inertias and bias forces
         *          3. Solve: compute joint accelerations from torques
         * 
         *          Mathematical formulation:
         *          q̈ = (τ - Sᵀ·pₐ) / (Sᵀ·Iₐ·S)
         *          where Iₐ is articulated inertia and pₐ is bias force
         * 
         * @throws std::invalid_argument if tau.size() != links.size()
         * @throws std::runtime_error if denominator is near zero (singular configuration)
         */
        void computeAccelerations(const Eigen::VectorXd& tau);

    private:
        /**
         * @brief Outward pass: propagate velocities and compute bias accelerations
         * @details Iterates from base (index 0) to tip (index n-1).
         *          For each link:
         *          - Transform parent velocity to current frame
         *          - Add joint velocity contribution: v = X⁻¹·v_parent + S·q̇
         *          - Compute bias acceleration: c = v × S · q̇
         * 
         *          Base link special case:
         *          - v₀ = S₀·q̇₀ (no parent velocity to transform)
         *          - c₀ = 0 (no velocity product terms)
         */
        void outwardPass();

        /**
         * @brief Inward pass: accumulate articulated inertias and bias forces
         * @param tau Vector of joint torques
         * @details Iterates from tip (index n-1) to base (index 0).
         *          For each link:
         *          - Initialize Ia with rigid body inertia I
         *          - Add transformed child articulated inertias
         *          - Compute bias force: pₐ = Iₐ·c + f
         *          - Transform Ia and pₐ to parent frame
         * 
         *          After inward pass, solves for joint accelerations:
         *          q̈ = (τ - Sᵀ·pₐ) / (Sᵀ·Iₐ·S)
         * 
         * @note Children must be processed before their parent
         */
        void inwardPass(const Eigen::VectorXd& tau);
    };
}

#endif // FORWARD_DYNAMICS_H
