/**
 * @file bench_aba_stub.cpp
 * @brief ABA forward dynamics benchmark stub
 * @details STUB — Phase 17 fills implementation with actual benchmark logic.
 *          Current stub suppresses unused parameters to satisfy linker.
 *          Forward declarations:
 *          - BM_ABA_ForwardDynamics(benchmark::State&, int) registered in bench_all.cpp
 *
 *          Implementation pattern (Phase 17):
 *          1. Create model via ModelFactory::createFD(nDOF)
 *          2. Generate random state via RandomState
 *          3. Apply random state to solver
 *          4. For (auto _ : state): computeAccelerations(tau)
 *          5. benchmark::DoNotOptimize(result)
 */
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement ABA benchmark with DOF sweep
    // Pattern: create model via ModelFactory::createFD(nDOF),
    //          generate random state via RandomState,
    //          apply state, time computeAccelerations in loop,
    //          DoNotOptimize(result)
    (void)state;
    (void)nDOF;
}
