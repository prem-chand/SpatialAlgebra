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

private:
    std::mt19937 rng_;  ///< Mersenne Twister PRNG with fixed seed
};

}  // namespace SpatialAlgebra::Bench
