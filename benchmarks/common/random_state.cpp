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

/**
 * @brief Generate nDOF random joint torques
 * @param nDOF Number of joints
 * @return Vector of nDOF torques in [-10, 10] Nm (D-16)
 */
std::vector<double> RandomState::randomTorques(int nDOF) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    std::vector<double> torques(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        torques[i] = dist(rng_);
    }
    return torques;
}

/**
 * @brief Generate nDOF random joint accelerations
 * @param nDOF Number of joints
 * @return Vector of nDOF accelerations in [-5, 5] rad/s² (D-17)
 */
std::vector<double> RandomState::randomAccelerations(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> accelerations(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        accelerations[i] = dist(rng_);
    }
    return accelerations;
}

/**
 * @brief Zero-alloc fill n random positions into pre-allocated buffer (D-19)
 * @param buf Pre-allocated buffer of size n
 * @param n Number of positions to fill
 * @details Distribution matches randomPositions for identical sequences.
 */
void RandomState::fillPositions(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    for (int i = 0; i < n; ++i) {
        buf[i] = dist(rng_);
    }
}

/**
 * @brief Zero-alloc fill n random velocities into pre-allocated buffer (D-19)
 * @param buf Pre-allocated buffer of size n
 * @param n Number of velocities to fill
 * @details Distribution matches randomVelocities for identical sequences.
 */
void RandomState::fillVelocities(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i) {
        buf[i] = dist(rng_);
    }
}

/**
 * @brief Zero-alloc fill n random torques into pre-allocated buffer (D-19)
 * @param buf Pre-allocated buffer of size n
 * @param n Number of torques to fill
 * @details Distribution matches randomTorques for identical sequences.
 */
void RandomState::fillTorques(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < n; ++i) {
        buf[i] = dist(rng_);
    }
}

/**
 * @brief Zero-alloc fill n random accelerations into pre-allocated buffer (D-19)
 * @param buf Pre-allocated buffer of size n
 * @param n Number of accelerations to fill
 * @details Distribution matches randomAccelerations for identical sequences.
 */
void RandomState::fillAccelerations(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i) {
        buf[i] = dist(rng_);
    }
}

}  // namespace SpatialAlgebra::Bench
