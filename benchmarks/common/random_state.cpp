/**
 * @file random_state.cpp
 * @brief Implementation of RandomState class
 * @details Implements random joint position and velocity generation with
 *          uniform distributions in the ranges specified by D-09 and D-10.
 *          All sequences are deterministic with fixed seed 42 (D-11).
 */

#include "random_state.h"
#include <cmath>

namespace SpatialAlgebra::Bench {

std::vector<double> RandomState::randomPositions(int nDOF) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    std::vector<double> positions(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        positions[i] = dist(rng_);
    }
    return positions;
}

std::vector<double> RandomState::randomVelocities(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> velocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        velocities[i] = dist(rng_);
    }
    return velocities;
}

}  // namespace SpatialAlgebra::Bench
