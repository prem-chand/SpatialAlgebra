#ifndef COMPOSITE_RIGID_BODY_INERTIA_H
#define COMPOSITE_RIGID_BODY_INERTIA_H

/**
 * @file CompositeRigidBodyInertia.h
 * @brief Composite Rigid Body Algorithm (CRBA) for joint-space inertia matrix
 * @details Implements Featherstone's Composite Rigid Body Algorithm (Algorithm 7.4)
 *          for O(n²) computation of the joint-space inertia matrix H(q). The CRBA
 *          accumulates subtree rigid body inertias from leaves to root via an inward
 *          pass, then fills the mass matrix by walking up the kinematic chain for
 *          each body.
 *
 *          Complexity: O(n²) operations
 *          - Inward pass (composite inertias): O(n) — single leaf-to-root sweep
 *          - Mass matrix fill: O(n²/2) — for each body, walk to root
 *
 *          Reference:
 *          Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 *          Springer, Boston, MA. Algorithm 7.4 (Composite Rigid-Body Inertia).
 *
 * @see ForwardDynamics.h for the Link struct definition
 * @see RigidBodyInertia.h for the inertia representation
 */

#include <Eigen/Dense>
#include <vector>

namespace SpatialAlgebra {

struct Link;

/**
 * @brief Compute H(q) using the Composite Rigid Body Algorithm
 * @param links Kinematic chain in topological order (parents before children).
 *              Each link must have valid X (parent→child Plücker transform),
 *              S (joint motion subspace), and I (rigid body inertia).
 *              These are typically populated by forwardKinematics().
 * @return n×n symmetric positive-definite mass matrix
 * @throws std::invalid_argument if links is empty (returns empty matrix)
 *
 * Algorithm:
 *   1. Initialize composite inertias: C[i] = I[i]
 *   2. Inward pass (leaf to root): C[parent] += X[i]^(-*) · C[i] (child→parent)
 *   3. Mass matrix fill: for each body, F = C[i] · S[i], then walk ancestors
 *      accumulating H[j][i] = S[j]^T · F_j where F_j is F transformed to frame j
 *
 * @note Requires forwardKinematics() to have been called before to ensure
 *       X and S are valid for the current joint configuration q.
 * @note Only supports serial chains (branching not tested). The algorithm
 *       extends naturally to branching trees by summing contributions from
 *       all children in the inward pass.
 */
Eigen::MatrixXd computeCompositeRigidBodyInertia(std::vector<Link>& links);

} // namespace SpatialAlgebra

#endif // COMPOSITE_RIGID_BODY_INERTIA_H
