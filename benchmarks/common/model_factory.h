#pragma once

/**
 * @file model_factory.h
 * @brief Unified factory for creating arbitrary n-DOF ForwardDynamics and InverseDynamics solver objects
 * @details Provides a ModelFactory class that constructs serial and branching kinematic chains
 *          with configurable joint axes, inertia, and transforms. Shared utility consumed by
 *          all benchmark domains (ABA, RNEA, core microbenchmarks).
 *
 * Design decisions (from Phase 16 CONTEXT.md):
 * - D-01: Unified API returning both ForwardDynamics and InverseDynamics solver objects
 * - D-02: Configurable joint axis per joint (Z, X, Y, or arbitrary screw axis per link)
 * - D-03: Support both serial chains and branching (Y-shaped) configurations
 * - D-04: Link inertia configurable via parameters (mass, COM, inertia tensor)
 */

#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "MotionVector.h"

#include <Eigen/Dense>
#include <vector>

namespace SpatialAlgebra::Bench {

/**
 * @brief Configuration for a single joint/link in a kinematic chain
 * @details Provides default values that create a standard Z-revolute link with
 *          1 kg mass, COM at origin, and half-meter translation along X.
 *          All fields are individually configurable per D-02 and D-04.
 */
struct JointConfig {
    MotionVector axis;                    ///< Screw axis (default: Z revolute)
    Eigen::Vector3d translation;          ///< Transform from parent to this link
    double mass;                          ///< Per-link mass
    Eigen::Vector3d com;                  ///< Center of mass in link frame

    /**
     * @brief Default constructor with standard Z-revolute link defaults
     * @details Initializes:
     *          - axis: Z revolute MotionVector(Vector3d(0,0,1), Vector3d::Zero())
     *          - translation: Vector3d::UnitX() * 0.5 (half-meter along X)
     *          - mass: 1.0 kg
     *          - com: Vector3d::Zero() (at origin)
     */
    JointConfig()
        : axis(MotionVector(Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero())),
          translation(Eigen::Vector3d::UnitX() * 0.5),
          mass(1.0),
          com(Eigen::Vector3d::Zero()) {}
};

/**
 * @brief Unified factory for constructing ForwardDynamics and InverseDynamics solver objects
 * @details Supports serial chains (uniform or per-link config), branching (Y-shaped) chains,
 *          and both solver types via separate createFD/createID methods. All methods return
 *          fully-populated solver objects with links initialized to zero state.
 *
 *          Usage:
 *          @code
 *          ModelFactory factory;
 *          ForwardDynamics fd = factory.createFD(6);  // 6-DOF Z-revolute serial chain
 *          ForwardDynamics fd2 = factory.createFD(6, JointConfig());  // explicit config
 *          ForwardDynamics fd3 = factory.createFDBranching(10, 4);    // Y-shaped at link 4
 *          ForwardDynamics fd4 = factory.createFD(configs);           // per-link configs
 *          @endcode
 */
class ModelFactory {
public:
    /**
     * @brief Create a serial-chain ForwardDynamics solver
     * @param nDOF Number of degrees of freedom (links) in the chain
     * @param cfg Joint configuration applied uniformly to all links
     * @return ForwardDynamics solver with nDOF links initialized to zero state
     */
    ForwardDynamics createFD(int nDOF, const JointConfig& cfg = JointConfig());

    /**
     * @brief Create a serial-chain InverseDynamics solver
     * @param nDOF Number of degrees of freedom (links) in the chain
     * @param cfg Joint configuration applied uniformly to all links
     * @return InverseDynamics solver with nDOF links initialized to zero state
     */
    InverseDynamics createID(int nDOF, const JointConfig& cfg = JointConfig());

    /**
     * @brief Create a branching (Y-shaped) ForwardDynamics solver
     * @param nDOF Total number of links (must be >= branchPoint + 3)
     * @param branchPoint Index where the chain splits into two branches
     * @param cfg Joint configuration applied uniformly to all links
     * @return ForwardDynamics solver with a Y-shaped kinematic tree
     * @note branchPoint must be in [0, nDOF-3] to leave room for at least one link per branch
     */
    ForwardDynamics createFDBranching(int nDOF, int branchPoint,
                                      const JointConfig& cfg = JointConfig());

    /**
     * @brief Create a branching (Y-shaped) InverseDynamics solver
     * @param nDOF Total number of links (must be >= branchPoint + 3)
     * @param branchPoint Index where the chain splits into two branches
     * @param cfg Joint configuration applied uniformly to all links
     * @return InverseDynamics solver with a Y-shaped kinematic tree
     * @note branchPoint must be in [0, nDOF-3] to leave room for at least one link per branch
     */
    InverseDynamics createIDBranching(int nDOF, int branchPoint,
                                      const JointConfig& cfg = JointConfig());

    /**
     * @brief Create a ForwardDynamics solver with per-link configurations
     * @param configs Vector of JointConfig, one per link (size determines nDOF)
     * @return ForwardDynamics solver with each link configured individually
     * @details Enables explicit per-link specification of joint axes, inertia, and transforms
     *          per D-04. Each link i uses configs[i] for its parameters.
     */
    ForwardDynamics createFD(const std::vector<JointConfig>& configs);

    /**
     * @brief Create an InverseDynamics solver with per-link configurations
     * @param configs Vector of JointConfig, one per link (size determines nDOF)
     * @return InverseDynamics solver with each link configured individually
     * @details Enables explicit per-link specification of joint axes, inertia, and transforms
     *          per D-04. Each link i uses configs[i] for its parameters.
     */
    InverseDynamics createID(const std::vector<JointConfig>& configs);
};

}  // namespace SpatialAlgebra::Bench
