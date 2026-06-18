/**
 * @file bench_all.cpp
 * @brief Benchmark entry point with programmatic registration
 * @details Provides a custom main() function that registers all benchmark
 *          functions programmatically via benchmark::RegisterBenchmark with
 *          a DOF sweep (n=1..20) for solver benchmarks and individual
 *          registrations for core microbenchmarks.
 */
#include <benchmark/benchmark.h>

// Forward declarations for solver benchmarks (DOF sweep)
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF);
void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF);

// Forward declarations for core microbenchmarks (18 individual, no DOF param)
void BM_TransformMotion(benchmark::State& state);
void BM_TransformForce(benchmark::State& state);
void BM_InverseTransformMotion(benchmark::State& state);
void BM_InverseTransformForce(benchmark::State& state);
void BM_TformRBI(benchmark::State& state);
void BM_InvtformRBI(benchmark::State& state);
void BM_TformABI(benchmark::State& state);
void BM_InvtformABI(benchmark::State& state);
void BM_CrossMvMv(benchmark::State& state);
void BM_CrossMvFv(benchmark::State& state);
void BM_CrossFvFv(benchmark::State& state);
void BM_RBIApply(benchmark::State& state);
void BM_RBIAdd(benchmark::State& state);
void BM_RBIScale(benchmark::State& state);
void BM_ABIApply(benchmark::State& state);
void BM_ABIAdd(benchmark::State& state);
void BM_ABIAddRBI(benchmark::State& state);
void BM_ABIScale(benchmark::State& state);

int main(int argc, char** argv) {
    // Register solver benchmarks with DOF sweep n=1..20 per D-05
    for (int n = 1; n <= 20; ++n) {
        benchmark::RegisterBenchmark(
            ("BM_ABA_ForwardDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_ABA_ForwardDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_RNEA_InverseDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_RNEA_InverseDynamics, n
        );
    }

    // Register core microbenchmarks (18 individual, no DOF sweep — O(1) operations)
    benchmark::RegisterBenchmark("BM_TransformMotion", BM_TransformMotion);
    benchmark::RegisterBenchmark("BM_TransformForce", BM_TransformForce);
    benchmark::RegisterBenchmark("BM_InverseTransformMotion", BM_InverseTransformMotion);
    benchmark::RegisterBenchmark("BM_InverseTransformForce", BM_InverseTransformForce);
    benchmark::RegisterBenchmark("BM_TformRBI", BM_TformRBI);
    benchmark::RegisterBenchmark("BM_InvtformRBI", BM_InvtformRBI);
    benchmark::RegisterBenchmark("BM_TformABI", BM_TformABI);
    benchmark::RegisterBenchmark("BM_InvtformABI", BM_InvtformABI);
    benchmark::RegisterBenchmark("BM_CrossMvMv", BM_CrossMvMv);
    benchmark::RegisterBenchmark("BM_CrossMvFv", BM_CrossMvFv);
    benchmark::RegisterBenchmark("BM_CrossFvFv", BM_CrossFvFv);
    benchmark::RegisterBenchmark("BM_RBIApply", BM_RBIApply);
    benchmark::RegisterBenchmark("BM_RBIAdd", BM_RBIAdd);
    benchmark::RegisterBenchmark("BM_RBIScale", BM_RBIScale);
    benchmark::RegisterBenchmark("BM_ABIApply", BM_ABIApply);
    benchmark::RegisterBenchmark("BM_ABIAdd", BM_ABIAdd);
    benchmark::RegisterBenchmark("BM_ABIAddRBI", BM_ABIAddRBI);
    benchmark::RegisterBenchmark("BM_ABIScale", BM_ABIScale);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
