/**
 * @file bench_rnea_stub.cpp
 * @brief RNEA inverse dynamics benchmark
 * @details Implements BM_RNEA_InverseDynamics with DOF sweep n=1..20.
 *          Each iteration generates fresh random accelerations via
 *          fillAccelerations (zero-alloc, D-19) in range [-5, 5] rad/s² (D-17).
 *          Joint state (q, qdot) is set once outside the timed loop.
 *
 *          Implementation:
 *          1. Create model via ModelFactory::createID(nDOF)
 *          2. Generate random state via RandomState
 *          3. Apply random state to solver links
 *          4. For (auto _ : state): fillAccelerations + computeTorques
 *          5. benchmark::DoNotOptimize(result)
 *
 *          Note: computeTorques allocates a new Eigen::VectorXd each call
 *          (InverseDynamics.cpp:57), which is inside the timed loop. This is
 *          inherent to the current API and represents the baseline behavior.
 */
#include <benchmark/benchmark.h>
#include "InverseDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF) {
    ModelFactory factory;
    RandomState rng;
    InverseDynamics id = factory.createID(nDOF);

    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        id.links[i].q = q[i];
        id.links[i].qdot = qdot[i];
    }

    Eigen::VectorXd qddot(nDOF);

    for (auto _ : state) {
        rng.fillAccelerations(qddot.data(), nDOF);
        Eigen::VectorXd tau = id.computeTorques(qddot);
        benchmark::DoNotOptimize(tau[0]);
    }
}
