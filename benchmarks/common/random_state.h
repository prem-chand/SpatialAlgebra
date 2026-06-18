#pragma once

/**
 * @file random_state.h
 * @brief Deterministic random joint state generator for benchmarks
 * @details Provides a RandomState class that generates random joint positions
 *          and velocities using a fixed-seed Mersenne Twister (seed 42) for
 *          fully reproducible benchmark results across runs.
 *
 * Design decisions (from Phase 16 CONTEXT.md):
 * - D-09: Joint positions uniform random in [-pi/2, pi/2]
 * - D-10: Joint velocities uniform random in [-5, 5] rad/s
 * - D-11: Fixed random seed 42 for reproducibility
 */

#include <random>
#include <vector>

namespace SpatialAlgebra::Bench {

/**
 * @brief Deterministic random joint state generator
 * @details Generates random joint positions and velocities for use in benchmark
 *          setup. Uses std::mt19937 with a fixed seed (42) to ensure identical
 *          random sequences across all benchmark runs, enabling reproducible
 *          performance measurements.
 *
 *          Usage:
 *          @code
 *          RandomState rng;
 *          std::vector<double> q = rng.randomPositions(6);     // 6-DOF positions
 *          std::vector<double> qdot = rng.randomVelocities(6); // 6-DOF velocities
 *          @endcode
 */
class RandomState {
public:
    /**
     * @brief Construct with fixed seed 42 for reproducibility (D-11)
     */
    RandomState() : rng_(42) {}

    /**
     * @brief Generate nDOF random joint positions
     * @param nDOF Number of joints
     * @return Vector of nDOF positions, each uniform random in [-pi/2, pi/2] (D-09)
     */
    std::vector<double> randomPositions(int nDOF);

    /**
     * @brief Generate nDOF random joint velocities
     * @param nDOF Number of joints
     * @return Vector of nDOF velocities, each uniform random in [-5, 5] rad/s (D-10)
     */
    std::vector<double> randomVelocities(int nDOF);

    // Phase 17: Torque/acceleration generation (D-16, D-17, D-19)

    /**
     * @brief Generate nDOF random joint torques
     * @param nDOF Number of joints
     * @return Vector of nDOF torques, each uniform random in [-10, 10] Nm (D-16)
     */
    std::vector<double> randomTorques(int nDOF);

    /**
     * @brief Generate nDOF random joint accelerations
     * @param nDOF Number of joints
     * @return Vector of nDOF accelerations, each uniform random in [-5, 5] rad/s² (D-17)
     */
    std::vector<double> randomAccelerations(int nDOF);

    /**
     * @brief Zero-alloc fill n random positions into pre-allocated buffer (D-19)
     * @param buf Pre-allocated buffer of size n
     * @param n Number of positions to fill
     * @details Distribution matches randomPositions for identical sequences.
     */
    void fillPositions(double* buf, int n);

    /**
     * @brief Zero-alloc fill n random velocities into pre-allocated buffer (D-19)
     * @param buf Pre-allocated buffer of size n
     * @param n Number of velocities to fill
     * @details Distribution matches randomVelocities for identical sequences.
     */
    void fillVelocities(double* buf, int n);

    /**
     * @brief Zero-alloc fill n random torques into pre-allocated buffer (D-19)
     * @param buf Pre-allocated buffer of size n
     * @param n Number of torques to fill
     * @details Distribution matches randomTorques for identical sequences.
     */
    void fillTorques(double* buf, int n);

    /**
     * @brief Zero-alloc fill n random accelerations into pre-allocated buffer (D-19)
     * @param buf Pre-allocated buffer of size n
     * @param n Number of accelerations to fill
     * @details Distribution matches randomAccelerations for identical sequences.
     */
    void fillAccelerations(double* buf, int n);

private:
    std::mt19937 rng_;  ///< Mersenne Twister PRNG with fixed seed
};

}  // namespace SpatialAlgebra::Bench
