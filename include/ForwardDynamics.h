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
 *          2. Inward pass (tip to base): Single tip-to-base sweep with condensation per
 *             Featherstone Algorithm 7.3 — for each link, accumulate children contributions,
 *             solve qddot = (tau - S^T*pa)/(S^T*Ia*S), condense Ia/pa, pass to parent
 *          3. Forward pass (base to tip): Compute spatial accelerations a = X*a_parent + c + S*qddot
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

    // Type alias to distinguish Link type for ForwardDynamics
    using ForwardDynamicsLink = Link;

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
         *          Executes the ABA:
         *          1. Outward pass: propagate velocities, compute bias accelerations
         *          2. Inward pass (single sweep with condensation): accumulate articulated
         *             inertias, solve qddot with condensation, propagate to parent
         *          3. Forward pass: compute spatial accelerations from qddot
         * 
         *          Mathematical formulation:
         *          q̈ = (τ - Sᵀ·pₐ) / (Sᵀ·Iₐ·S)
         *          where Iₐ is articulated inertia and pₐ is bias force
         * 
         * @param gravity Gravity vector (default zero) for gravity-aware dynamics
         * @throws std::invalid_argument if tau.size() != links.size()
         * @throws std::runtime_error if denominator is near zero (singular configuration)
         */
        void computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity = Vector3d::Zero());

    private:
        Vector3d gravity;  ///< Gravity vector (default zero, set by computeAccelerations)

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
         * @brief Inward pass: single tip-to-base sweep with condensation per Featherstone Algorithm 7.3
         * @param tau Vector of joint torques
         * @details Single tip-to-base sweep computing articulated inertias and solving joint
         *          accelerations. For each link (tip to base):
         *          - Children's condensed Ia/pa have been accumulated from prior iterations
         *          - Solve qddot = (tau - S^T*pa) / (S^T*Ia*S) as the FINAL joint acceleration
         *          - Condense Ia: I_A' = I_A - (I_A*S) * D^{-1} * (I_A*S)^T if parent exists
         *          - Condense pa: p_A' = p_A + I_A*S * qddot if parent exists
         *          - Transform condensed Ia/pa to parent frame
         * 
         *          After the sweep, a forward pass computes spatial accelerations:
         *          a_i = X_i * a_parent + c_i + S_i * qddot_i
         * 
         *          The condensation step (Featherstone Algorithm 7.3) removes the joint motion
         *          subspace before propagating to the parent, preventing overestimation of the
         *          reflected subtree inertia.
         * 
         * @note Base link (parent == -1) skips condensation — no parent to pass to
         * @note Children must be processed before their parent (tip-to-base order)
         * @throws std::runtime_error if denominator D = S^T*Ia*S is near zero
         * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7, Algorithm 7.3
         */
        void inwardPass(const Eigen::VectorXd& tau);
    };
}

#endif // FORWARD_DYNAMICS_H
