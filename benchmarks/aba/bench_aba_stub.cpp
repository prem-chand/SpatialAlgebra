/**
 * @file bench_aba_stub.cpp
 * @brief ABA forward dynamics benchmark
 * @details Implements BM_ABA_ForwardDynamics with DOF sweep n=1..20.
 *          Each iteration generates fresh random torques via fillTorques
 *          (zero-alloc, D-19) in range [-10, 10] Nm (D-16).
 *          Joint state (q, qdot) is set once outside the timed loop.
 *
 *          Implementation:
 *          1. Create model via ModelFactory::createFD(nDOF)
 *          2. Generate random state via RandomState
 *          3. Apply random state to solver links
 *          4. For (auto _ : state): fillTorques + computeAccelerations
 *          5. benchmark::DoNotOptimize(result)
 */
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    ModelFactory factory;
    RandomState rng;
    ForwardDynamics fd = factory.createFD(nDOF);

    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        fd.links[i].q = q[i];
        fd.links[i].qdot = qdot[i];
    }

    Eigen::VectorXd tau(nDOF);

    for (auto _ : state) {
        rng.fillTorques(tau.data(), nDOF);
        fd.computeAccelerations(tau);
        benchmark::DoNotOptimize(fd.links[0].qddot);
    }
}
