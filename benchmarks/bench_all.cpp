/**
 * @file bench_all.cpp
 * @brief Benchmark entry point with programmatic registration
 * @details Provides a custom main() function that registers all benchmark
 *          functions programmatically via benchmark::RegisterBenchmark with
 *          a DOF sweep (n=1..20) for each benchmark domain.
 *
 *          Registration pattern (19-RESEARCH.md §Pattern 2):
 *          - Each benchmark is registered with a name like "BM_{Domain}/{N}DOF"
 *          - Runtime filtering: bench_all --benchmark_filter="ABA"
 *          - All 4 domains (ABA, RNEA, PluckerTransform, CrossProduct) are
 *            registered in a single DOF loop for a total of 80 benchmarks.
 *
 *          @see Phase 17: fills stub implementations with actual benchmark logic
 */
#include <benchmark/benchmark.h>
#include "common/model_factory.h"
#include "common/random_state.h"

// Forward declarations for benchmark functions (defined in domain stubs)
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF);
void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF);
void BM_PluckerTransform(benchmark::State& state, int nDOF);
void BM_CrossProduct(benchmark::State& state, int nDOF);

int main(int argc, char** argv) {
    // Register all benchmarks with DOF sweep n=1..20 per D-05
    // D-05 specifies single unified executable with --benchmark_filter support
    for (int n = 1; n <= 20; ++n) {
        benchmark::RegisterBenchmark(
            ("BM_ABA_ForwardDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_ABA_ForwardDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_RNEA_InverseDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_RNEA_InverseDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_PluckerTransform/" + std::to_string(n) + "DOF").c_str(),
            BM_PluckerTransform, n
        );
        benchmark::RegisterBenchmark(
            ("BM_CrossProduct/" + std::to_string(n) + "DOF").c_str(),
            BM_CrossProduct, n
        );
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
