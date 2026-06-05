# Phase 17: Benchmark Implementation - Research

**Researched:** 2026-06-05
**Domain:** Benchmark implementation for ABA forward dynamics, RNEA inverse dynamics, and core spatial algebra microbenchmarks
**Confidence:** HIGH

## Summary

Phase 17 fills the benchmark stubs created in Phase 16 with real timing loops. There are three benchmark domains, each with distinct implementation needs:

1. **ABA (BENCH-01)**: Parameterized DOF sweep (n=1..20) using `ModelFactory::createFD(nDOF)` + `RandomState` state setup + `computeAccelerations(tau)` in the timed loop. Random torques regenerated per iteration via `fillTorques()` (D-16, D-19).

2. **RNEA (BENCH-02)**: Same DOF sweep with `ModelFactory::createID(nDOF)` + `computeTorques(qddot)`. Random qddot regenerated per iteration via `fillAccelerations()` (D-17, D-19).

3. **Core microbenchmarks (BENCH-03)**: 20 individual operations (8 Plücker transforms, 4 cross products, 7-8 inertia operations). Each is a separate benchmark function without DOF dependency — the current `nDOF` parameter is unused for core benchmarks.

**Primary recommendation:** Use per-operation state setup outside the timed loop, regenerate only the driving input (tau/qddot) inside the loop via pre-allocated fill buffers, and always consume results with `DoNotOptimize`. The 8 Plücker inertia transforms (tformRBI, invtformRBI, tformABI, invtformABI) require pre-constructed inertia objects. Existing `bench_all.cpp` structure needs updating: replace the 2 core placeholder registrations with 20 individual operations.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| BENCH-01 | ABA forward dynamics timing with DOF sweep (n=1..20) | ABA uses O(n) outward→inward→correction pass. `computeAccelerations(const Eigen::VectorXd& tau)` is the timed call. Model setup + state setup outside loop. Torques regenerated per iteration (D-16). Verified against `ForwardDynamics.cpp:182-209` implementation. |
| BENCH-02 | RNEA inverse dynamics timing with DOF sweep (n=1..20) | RNEA uses O(n) outward→inward pass. `computeTorques(const Eigen::VectorXd& qddot)` is the timed call. Joint accelerations regenerated per iteration (D-17). Returns `Eigen::VectorXd tau`. Verified against `InverseDynamics.cpp:96-136`. |
| BENCH-03 | Plücker transform and cross-product microbenchmarks | D-18 expands scope to 20 individual operations (8 Plücker, 4 cross, 7-8 inertia). These are O(1) operations — DOF parameter is meaningless. Existing `bench_all.cpp` registration needs restructuring for individual benchmarks. |

</phase_requirements>

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-16:** Random torques via `RandomState::randomTorques(nDOF)` with uniform range [-10, 10] Nm
- **D-17:** Random qddot via `RandomState::randomAccelerations(nDOF)` with uniform range [-5, 5] rad/s²
- **D-18:** Full API surface — all 20 operations benchmarked:
  - Plücker transforms (8): `transformMotion`, `transformForce`, `inverseTransformMotion`, `inverseTransformForce`, `tformRBI`, `invtformRBI`, `tformABI`, `invtformABI`
  - Cross products (4): `cross(mv,mv)`, `cross(mv,fv)`, `cross(fv,mv)`, `cross(fv,fv)`
  - Inertia operations (8): All `apply` variants on RigidBodyInertia and ArticulatedBodyInertia
- **D-19:** Per-iteration random state generation using pre-allocated `Eigen::VectorXd` buffers — add non-allocating `fillPositions(double*, int)` and `fillVelocities(double*, int)` methods to `RandomState` (and corresponding `fillTorques`, `fillAccelerations`) to avoid heap allocation inside the timed loop.
- **D-09:** Joint positions: uniform random in [-π/2, π/2]
- **D-10:** Joint velocities: uniform random in [-5, 5] rad/s
- **D-11:** Fixed random seed 42 for reproducibility across runs
- **D-05:** Single unified executable (`bench_all`) with Google Benchmark sub-benchmark registration — benchmarks filterable via `--benchmark_filter`
- **D-06:** Subdirectory layout: `benchmarks/{aba,rnea,core,common}/`

### Claude's Discretion
- Exact benchmark function signatures (within `void BM_*(benchmark::State&, int nDOF)` pattern)
- DOF iteration step size (1 or configurable) within the 1..20 range
- Google Benchmark `MinTime` / iterations configuration
- Specific `DoNotOptimize` / `ClobberMemory` placement in microbenchmark loops

### Deferred Ideas (OUT OF SCOPE)
None — discussion stayed within phase scope.
</user_constraints>

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| ABA benchmark (BENCH-01) | Benchmarks ABA | RandomState (extension) | `benchmarks/aba/bench_aba_stub.cpp` — timed `computeAccelerations()` loop. Consumes ModelFactory for FD solver, RandomState for state/torque generation. |
| RNEA benchmark (BENCH-02) | Benchmarks RNEA | RandomState (extension) | `benchmarks/rnea/bench_rnea_stub.cpp` — timed `computeTorques()` loop. Consumes ModelFactory for ID solver, RandomState for state/acceleration generation. |
| Core microbenchmarks (BENCH-03) | Benchmarks Core | — | `benchmarks/core/bench_core_stub.cpp` — 20 self-contained benchmark functions operating on pre-constructed transforms/vectors/inertias. No solver dependency. |
| RandomState extension | Benchmarks Common | — | Add `randomTorques()`, `randomAccelerations()`, fill* methods. Consumed by ABA and RNEA benchmarks. No impact on existing code beyond extension. |
| bench_all.cpp registration | Benchmarks (entry point) | — | `bench_all.cpp` needs updating: replace 2 core placeholders with 20 individual registrations. ABA/RNEA DOF sweep loop stays unchanged. |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Google Benchmark | v1.9.5 | C++ microbenchmarking framework | Industry standard. Handles warm-up, iteration count, statistics, dead-code prevention. Already integrated via FetchContent in Phase 16. |
| Eigen3 | 3.4...5 (existing) | Linear algebra types used in benchmarked operations | All solver APIs, transforms, and spatial vectors use Eigen3 types (`Vector3d`, `Matrix3d`, `VectorXd`). Pre-allocated buffers are `Eigen::VectorXd` with double* data access. |
| ModelFactory | N/A (benchmark utility) | Construct arbitrary n-DOF solvers | Unified API for both FD and ID solver construction. Phase 16 created it; Phase 17 consumes it. |
| RandomState | N/A (benchmark utility) | Deterministic random state generation | Phase 16 created it with `randomPositions()`/`randomVelocities()`. Phase 17 extends with torque/acceleration generation and zero-alloc fill methods. |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| C++ `<random>` | C++17 std | RNG for RandomState extensions | New methods reuse the same `std::mt19937` instance (seed 42) with `std::uniform_real_distribution<double>` for torque and acceleration ranges. |
| `benchmark::DoNotOptimize` | v1.9.5 | Prevent dead-code elimination | Used in ALL benchmark loops to ensure the compiler doesn't optimize away the operation under test. |
| `benchmark::ClobberMemory` | v1.9.5 | Memory barrier for loop optimization prevention | Used when `DoNotOptimize` alone is insufficient — forces the compiler to reload all memory after the barrier. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Per-operation core benchmarks (20 individual functions) | Single compound function per category | Individual benchmarks enable `--benchmark_filter="TransformMotion"` and per-op CSV output. Compound functions prevent per-operation profiling. |
| Per-iteration fill via RandomState | Pre-generate all torques/accelerations outside loop | Fill methods are zero-alloc (D-19). Pre-generation stores nDOF values and reads them sequentially, which is equivalent. Fill is cleaner and uses less memory. |

**Installation:**
```bash
# Phase 17 adds no new dependencies. Build with:
cmake -B build -DSA_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
# Run all benchmarks:
build/benchmarks/bench_all
# Filter by domain:
build/benchmarks/bench_all --benchmark_filter="ABA"
# Filter by DOF:
build/benchmarks/bench_all --benchmark_filter="BM_RNEA.*12DOF"
# CSV output for analysis:
build/benchmarks/bench_all --benchmark_format=csv
```

## Package Legitimacy Audit

Phase 17 introduces zero new packages. All dependencies (Google Benchmark v1.9.5 via FetchContent, Eigen3, C++ standard library) are already present from Phase 16.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| Google Benchmark (FetchContent) | GitHub | 12+ yrs | 10k+ stars | github.com/google/benchmark | N/A (verified OSS) | Approved — v1.9.5 tag verified |
| Eigen3 | Homebrew/system | N/A | N/A | gitlab.com/libeigen/eigen | N/A (verified) | Already in project |
| C++ `<random>` | Standard library | C++11+ | N/A | Standard | N/A | No install needed |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

## Architecture Patterns

### System Architecture Diagram

```
┌──────────────────────────────────────────────────────────────────┐
│                        bench_all/main()                          │
│                                                                  │
│  for (n = 1..20) {                                               │
│    RegisterBenchmark("BM_ABA_ForwardDynamics/{n}DOF", ...)       │
│    RegisterBenchmark("BM_RNEA_InverseDynamics/{n}DOF", ...)      │
│  }                                                                │
│                                                                  │
│  RegisterBenchmark("BM_TransformMotion", ...)     ─┐              │
│  RegisterBenchmark("BM_TransformForce", ...)        │              │
│  RegisterBenchmark("BM_InverseTransformMotion", ...)│ 20 core     │
│  RegisterBenchmark("BM_InverseTransformForce", ...) │ micro-     │
│  RegisterBenchmark("BM_TformRBI", ...)              │ benchmarks │
│  RegisterBenchmark("BM_InvtformRBI", ...)           │ (no DOF    │
│  RegisterBenchmark("BM_TformABI", ...)              │ sweep)     │
│  RegisterBenchmark("BM_InvtformABI", ...)          ─┘              │
│  RegisterBenchmark("BM_CrossMvMv", ...)            ─┐              │
│  RegisterBenchmark("BM_CrossMvFv", ...)              │ 4 cross    │
│  RegisterBenchmark("BM_CrossFvMv", ...)              │ products   │
│  RegisterBenchmark("BM_CrossFvFv", ...)             ─┘              │
│  RegisterBenchmark("BM_RBIApply", ...)              ─┐              │
│  ... (7 inertia operations)                           │ 7-8 inertia │
│                                                       │ operations  │
│    Initialize(); RunSpecifiedBenchmarks(); Shutdown(); │             │
│  }                                                     │             │
└──────────────────────┬───────────────────────────────┘             │
                       │                                              │
          ┌────────────┴────────────┐                                  │
          ▼                         ▼                                  │
┌───────────────────┐   ┌───────────────────────┐                      │
│  ABA/RNEA Flow    │   │  Core Micro Flow       │                     │
│                   │   │                        │                     │
│  1. CreateFD/ID   │   │  1. Pre-construct      │                     │
│  2. Set q/qdot    │   │     transforms,         │                     │
│  3. Pre-alloc     │   │     vectors, inertias   │                     │
│     tau/qddot buf │   │  2. for (auto _ : st) { │                     │
│  4. for (auto _   │   │     operation()         │                     │
│     : state) {    │   │     DoNotOptimize(out)  │                     │
│     rng.fill*(buf)│   │   }                     │                     │
│     solver.compute│   └───────────────────────┘                      │
│     DoNotOptimize │                                                    │
│   }               │                                                    │
└───────────────────┘                                                    │
                                                                         │
  Data flow for a single ABA/RNEA benchmark call:                        │
                                                                         │
  bench_all → BM_ABA_ForwardDynamics(state, nDOF)                        │
    ├── ModelFactory::createFD(nDOF)  → ForwardDynamics solver           │
    ├── RandomState::randomPositions(nDOF) → q state                     │
    ├── RandomState::randomVelocities(nDOF) → qdot state                 │
    ├── fd.links[i].q = q[i]; fd.links[i].qdot = qdot[i];               │
    ├── Eigen::VectorXd tau(nDOF);  // pre-allocated                     │
    │                                                                     │
    └── for (auto _ : state) {                                           │
         RandomState::fillTorques(tau.data(), nDOF);  // D-19           │
         fd.computeAccelerations(tau);                                   │
         DoNotOptimize(fd.links[0].qddot);                                │
        }                                                                 │
                                                                          │
  Data flow for a core microbenchmark:                                    │
                                                                          │
  bench_all → BM_TransformMotion(state, /*unused*/)                       │
    ├── PluckerTransform X = ...  // pre-constructed                      │
    ├── MotionVector mv_in = ...  // pre-constructed                      │
    │                                                                     │
    └── for (auto _ : state) {                                           │
         mv_out = X.transformMotion(mv_in);                               │
         DoNotOptimize(mv_out);                                            │
        }                                                                 │
```

### Recommended Project Structure

Phase 17 does NOT add new files. It modifies existing stubs and extends RandomState.

```
benchmarks/
├── CMakeLists.txt              # Unchanged from Phase 16
├── bench_all.cpp               # *** MODIFIED: Replace 2 core registrations with 20
├── common/
│   ├── CMakeLists.txt          # Unchanged (GLOB auto-discovers new sources)
│   ├── model_factory.h         # Unchanged
│   ├── model_factory.cpp       # Unchanged
│   ├── random_state.h          # *** MODIFIED: Add randomTorques/randomAccelerations + fill* methods
│   └── random_state.cpp        # *** MODIFIED: Implement new methods
├── aba/
│   ├── CMakeLists.txt          # Unchanged
│   └── bench_aba_stub.cpp      # *** REPLACED: Real ABA benchmark implementation
├── rnea/
│   ├── CMakeLists.txt          # Unchanged
│   └── bench_rnea_stub.cpp     # *** REPLACED: Real RNEA benchmark implementation
└── core/
    ├── CMakeLists.txt          # Unchanged
    └── bench_core_stub.cpp     # *** REPLACED: 20 individual core benchmark functions
```

### Pattern 1: ABA Forward Dynamics Benchmark

**What:** Parameterized DOF sweep (n=1..20) benchmarking `ForwardDynamics::computeAccelerations(tau)`. Timer measures the complete outward→inward→correction pass.

**When to use:** Only for benchmarking the ABA solver. Must include pre-allocated buffer for tau with per-iteration fill via RandomState (D-19).

**Example:**
```cpp
// benchmarks/aba/bench_aba_stub.cpp
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // ── Setup (outside timed loop) ──────────────────────────────────
    ModelFactory factory;
    RandomState rng;

    // Create n-DOF serial chain (D-01/D-02)
    ForwardDynamics fd = factory.createFD(nDOF);

    // Apply random joint state (D-09/D-10) — set once, reused for all iterations
    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        fd.links[i].q = q[i];
        fd.links[i].qdot = qdot[i];
    }

    // Pre-allocated torque buffer (D-19: zero-alloc inside timed loop)
    Eigen::VectorXd tau(nDOF);

    // ── Timed loop ──────────────────────────────────────────────────
    for (auto _ : state) {
        // D-19: regenerated torques each iteration (zero-alloc via fill)
        rng.fillTorques(tau.data(), nDOF);

        // D-16: random torques, computeAccelerations is the measured operation
        fd.computeAccelerations(tau);

        // Prevent dead-code elimination (Pitfall 2, 16-RESEARCH.md)
        benchmark::DoNotOptimize(fd.links[0].qddot);
    }
}
```

**Key observations:**
- `computeAccelerations()` **fully recomputes** all solver state (v, c, Ia, pa, qddot) from q/qdot each call — no manual state reset needed. Verified against `ForwardDynamics.cpp:57-80` (outwardPass) and `ForwardDynamics.cpp:82-180` (inwardPass).
- The solver writes to `fd.links[i].qddot`. Outputs are per-link — `DoNotOptimize(fd.links[0].qddot)` suffices to signal the compiler that qddot is consumed. For full safety, access all links' qddot values.
- `computeAccelerations` validates tau size and NaN/Inf each call (`ForwardDynamics.cpp:184-201`). The throw path is cold and branch-predictable — negligible overhead.

### Pattern 2: RNEA Inverse Dynamics Benchmark

**What:** Parameterized DOF sweep (n=1..20) benchmarking `InverseDynamics::computeTorques(qddot)`. Timer measures outward→inward pass.

**When to use:** Only for benchmarking the RNEA solver.

**Example:**
```cpp
// benchmarks/rnea/bench_rnea_stub.cpp
#include <benchmark/benchmark.h>
#include "InverseDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF) {
    // ── Setup (outside timed loop) ──────────────────────────────────
    ModelFactory factory;
    RandomState rng;

    // Create n-DOF serial chain (D-01/D-02)
    InverseDynamics id = factory.createID(nDOF);

    // Apply random joint state (D-09/D-10)
    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        id.links[i].q = q[i];
        id.links[i].qdot = qdot[i];
    }

    // Pre-allocated acceleration buffer (D-19: zero-alloc inside timed loop)
    // Range: [-5, 5] rad/s² (D-17)
    Eigen::VectorXd qddot(nDOF);

    // ── Timed loop ──────────────────────────────────────────────────
    for (auto _ : state) {
        // D-19: regenerated accelerations each iteration
        rng.fillAccelerations(qddot.data(), nDOF);

        // D-17: random qddot, computeTorques is the measured operation
        // Returns Eigen::VectorXd of joint torques
        Eigen::VectorXd tau = id.computeTorques(qddot);

        // Prevent dead-code elimination
        benchmark::DoNotOptimize(tau[0]);
    }
}
```

**Key observations:**
- `computeTorques()` returns a **newly allocated** `Eigen::VectorXd` each call (`InverseDynamics.cpp:57` — `Eigen::VectorXd::Zero(links.size())`). This heap allocation IS part of the timed operation — it's inherent to the solver's current API.
- The solver copies the qddot input into `links[i].qddot` each call (`InverseDynamics.cpp:121-123`), then recomputes v and a in outwardPass, and f/tau in inwardPass. Full recomputation every call.
- `computeTorques` also validates qddot for NaN/Inf each call (`InverseDynamics.cpp:109-117`). Cold path.

### Pattern 3: Core Microbenchmarks

**What:** 20 individual benchmark functions, each timing a single spatial algebra operation. No DOF parameter — the operation is identical regardless of chain length.

**When to use:** For isolating performance of individual spatial algebra operations. Each benchmark is a simple loop over a pre-constructed object, calling the operation and consuming with `DoNotOptimize`.

**Example — Plücker Transform:**
```cpp
// Excerpt from benchmarks/core/bench_core_stub.cpp
#include <benchmark/benchmark.h>
#include "PluckerTransform.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "SpatialUtils.h"
#include "RigidBodyInertia.h"
#include "ArticulatedBodyInertia.h"
#include "LowerTriangular.h"

using namespace SpatialAlgebra;

// ── Pre-constructed test objects (file scope or function-local static) ──
static PluckerTransform make_test_transform() {
    return PluckerTransform(
        Rotation(Eigen::AngleAxisd(0.5, Vector3d::UnitZ())),
        Vector3d(0.3, 0.1, 0.7)
    );
}

static MotionVector make_test_mv() {
    return MotionVector(Vector3d(0.2, 0.5, 0.8), Vector3d(0.1, 0.3, 0.6));
}

static ForceVector make_test_fv() {
    return ForceVector(Vector3d(0.4, 0.9, 0.2), Vector3d(0.7, 0.1, 0.5));
}

static RigidBodyInertia make_test_rbi() {
    return RigidBodyInertia(1.5, Vector3d(0.1, 0.2, 0.3),
                            LowerTriangular::Identity(3));
}

static ArticulatedBodyInertia make_test_abi() {
    return ArticulatedBodyInertia(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );
}

// ── 8 Plücker transform benchmarks ─────────────────────────────────
void BM_TransformMotion(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    MotionVector mv_in = make_test_mv();

    for (auto _ : state) {
        MotionVector mv_out = X.transformMotion(mv_in);
        benchmark::DoNotOptimize(mv_out);
    }
}

void BM_TransformForce(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ForceVector fv_in = make_test_fv();

    for (auto _ : state) {
        ForceVector fv_out = X.transformForce(fv_in);
        benchmark::DoNotOptimize(fv_out);
    }
}

void BM_InverseTransformMotion(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    MotionVector mv_in = make_test_mv();

    for (auto _ : state) {
        MotionVector mv_out = X.inverseTransformMotion(mv_in);
        benchmark::DoNotOptimize(mv_out);
    }
}

void BM_InverseTransformForce(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ForceVector fv_in = make_test_fv();

    for (auto _ : state) {
        ForceVector fv_out = X.inverseTransformForce(fv_in);
        benchmark::DoNotOptimize(fv_out);
    }
}

void BM_TformRBI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    RigidBodyInertia rbi_in = make_test_rbi();

    for (auto _ : state) {
        RigidBodyInertia rbi_out = X.tformRBI(rbi_in);
        benchmark::DoNotOptimize(rbi_out);
    }
}

void BM_InvtformRBI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    RigidBodyInertia rbi_in = make_test_rbi();

    for (auto _ : state) {
        RigidBodyInertia rbi_out = X.invtformRBI(rbi_in);
        benchmark::DoNotOptimize(rbi_out);
    }
}

void BM_TformABI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ArticulatedBodyInertia abi_in = make_test_abi();

    for (auto _ : state) {
        ArticulatedBodyInertia abi_out = X.tformABI(abi_in);
        benchmark::DoNotOptimize(abi_out);
    }
}

void BM_InvtformABI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ArticulatedBodyInertia abi_in = make_test_abi();

    for (auto _ : state) {
        ArticulatedBodyInertia abi_out = X.invtformABI(abi_in);
        benchmark::DoNotOptimize(abi_out);
    }
}

// ── 4 Cross product benchmarks ────────────────────────────────────
void BM_CrossMvMv(benchmark::State& state) {
    MotionVector mv1 = make_test_mv();
    MotionVector mv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));

    for (auto _ : state) {
        MotionVector result = cross(mv1, mv2);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossMvFv(benchmark::State& state) {
    MotionVector mv1 = make_test_mv();
    ForceVector fv1 = make_test_fv();

    for (auto _ : state) {
        ForceVector result = cross(mv1, fv1);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossFvMv(benchmark::State& state) {
    ForceVector fv1 = make_test_fv();
    MotionVector mv1 = make_test_mv();

    for (auto _ : state) {
        ForceVector result = cross(fv1, mv1);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossFvFv(benchmark::State& state) {
    ForceVector fv1 = make_test_fv();
    ForceVector fv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));

    for (auto _ : state) {
        ForceVector result = cross(fv1, fv2);
        benchmark::DoNotOptimize(result);
    }
}

// ── 7-8 Inertia operation benchmarks ──────────────────────────────
void BM_RBIApply(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();
    MotionVector mv = make_test_mv();

    for (auto _ : state) {
        ForceVector fv = rbi.apply(mv);
        benchmark::DoNotOptimize(fv);
    }
}

void BM_RBIAdd(benchmark::State& state) {
    RigidBodyInertia rbi1 = make_test_rbi();
    RigidBodyInertia rbi2(2.0, Vector3d(0.05, 0.1, 0.15),
                          LowerTriangular::Identity(3));

    for (auto _ : state) {
        RigidBodyInertia result = rbi1 + rbi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_RBIScale(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();

    for (auto _ : state) {
        RigidBodyInertia result = rbi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIApply(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    MotionVector mv = make_test_mv();

    for (auto _ : state) {
        ForceVector fv = abi.apply(mv);
        benchmark::DoNotOptimize(fv);
    }
}

void BM_ABIAdd(benchmark::State& state) {
    ArticulatedBodyInertia abi1 = make_test_abi();
    ArticulatedBodyInertia abi2(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );

    for (auto _ : state) {
        ArticulatedBodyInertia result = abi1 + abi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIAddRBI(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    RigidBodyInertia rbi = make_test_rbi();

    for (auto _ : state) {
        ArticulatedBodyInertia result = abi + rbi;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIScale(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();

    for (auto _ : state) {
        ArticulatedBodyInertia result = abi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}
```

### Pattern 4: RandomState Extension

**What:** Extends the existing `RandomState` class with `randomTorques()`, `randomAccelerations()`, and zero-alloc fill methods for all four quantities (positions, velocities, torques, accelerations).

**When to use:** Every benchmark that needs randomized solver input inside the timed loop uses the fill* variant. Benchmarks needing one-shot state generation use the returning variants.

**Header extensions (`random_state.h`):**
```cpp
#pragma once
#include <random>
#include <vector>
#include <Eigen/Dense>

namespace SpatialAlgebra::Bench {

class RandomState {
public:
    RandomState() : rng_(42) {}

    // Existing methods (Phase 16)
    std::vector<double> randomPositions(int nDOF);
    std::vector<double> randomVelocities(int nDOF);

    // Phase 17: New methods
    // Allocate + return variants
    std::vector<double> randomTorques(int nDOF);          // D-16: [-10, 10] Nm
    std::vector<double> randomAccelerations(int nDOF);    // D-17: [-5, 5] rad/s²

    // Zero-alloc fill variants (D-19)
    void fillPositions(double* buf, int n);
    void fillVelocities(double* buf, int n);
    void fillTorques(double* buf, int n);
    void fillAccelerations(double* buf, int n);

private:
    std::mt19937 rng_;
};

}
```

**Implementation extensions (`random_state.cpp`):**
```cpp
#include "random_state.h"
#include <cmath>

namespace SpatialAlgebra::Bench {

// ── Existing methods (unchanged) ──
std::vector<double> RandomState::randomPositions(int nDOF) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    std::vector<double> positions(nDOF);
    for (int i = 0; i < nDOF; ++i)
        positions[i] = dist(rng_);
    return positions;
}

std::vector<double> RandomState::randomVelocities(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> velocities(nDOF);
    for (int i = 0; i < nDOF; ++i)
        velocities[i] = dist(rng_);
    return velocities;
}

// ── Phase 17: New returning variants ──
std::vector<double> RandomState::randomTorques(int nDOF) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    std::vector<double> torques(nDOF);
    for (int i = 0; i < nDOF; ++i)
        torques[i] = dist(rng_);
    return torques;
}

std::vector<double> RandomState::randomAccelerations(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> accelerations(nDOF);
    for (int i = 0; i < nDOF; ++i)
        accelerations[i] = dist(rng_);
    return accelerations;
}

// ── Phase 17: Zero-alloc fill variants (D-19) ──
void RandomState::fillPositions(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillVelocities(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillTorques(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillAccelerations(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

}
```

**Important requirement:** The fill variants **must create a new `std::uniform_real_distribution` each call** (not reuse a stored one) to guarantee that calling `fillPositions()` produces the same sequence as calling `randomPositions()` for the equivalent n. Creating `std::uniform_real_distribution` is extremely cheap (no heap allocation, just stores two doubles) — the overhead is negligible compared to the RNG calls themselves.

### Pattern 5: Updated bench_all.cpp Registration

**What:** The existing `bench_all.cpp` registers ABA and RNEA with DOF sweep. The core placeholder registrations (BM_PluckerTransform, BM_CrossProduct) are replaced with 20 individual benchmark registrations.

**When to use:** Only for the unified bench_all executable.

**Example (new version of `bench_all.cpp`):**
```cpp
#include <benchmark/benchmark.h>

// Forward declarations from aba, rnea, core
void BM_ABA_ForwardDynamics(benchmark::State&, int);
void BM_RNEA_InverseDynamics(benchmark::State&, int);

// Forward declarations — 20 core microbenchmarks
void BM_TransformMotion(benchmark::State&);
void BM_TransformForce(benchmark::State&);
void BM_InverseTransformMotion(benchmark::State&);
void BM_InverseTransformForce(benchmark::State&);
void BM_TformRBI(benchmark::State&);
void BM_InvtformRBI(benchmark::State&);
void BM_TformABI(benchmark::State&);
void BM_InvtformABI(benchmark::State&);
void BM_CrossMvMv(benchmark::State&);
void BM_CrossMvFv(benchmark::State&);
void BM_CrossFvMv(benchmark::State&);
void BM_CrossFvFv(benchmark::State&);
void BM_RBIApply(benchmark::State&);
void BM_RBIAdd(benchmark::State&);
void BM_RBIScale(benchmark::State&);
void BM_ABIApply(benchmark::State&);
void BM_ABIAdd(benchmark::State&);
void BM_ABIAddRBI(benchmark::State&);
void BM_ABIScale(benchmark::State&);
// BM_RBIPrint / BM_ABIPrint — if 8th inertia op needed

int main(int argc, char** argv) {
    // ABA/RNEA: DOF sweep n=1..20
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

    // Core microbenchmarks: 20 individual operations
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
    benchmark::RegisterBenchmark("BM_CrossFvMv", BM_CrossFvMv);
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
```

**Note on inertia operation count:** D-18 specifies 8 inertia operations. Only 7 identifiable operations exist (3 on RBI: apply, add, scale; 4 on ABI: apply, add, add+RBI, scale). The 8th operation is ambiguous. Options:
- Add a benign 8th benchmark (e.g., `LowerTriangular::multiplySymmetric` as a standalone op, or `RigidBodyInertia::operator<<` / `print`)
- Accept the count as 7 and document the discrepancy
- **Recommendation:** Benchmark the 7 identifiable operations and note the count difference in the assumptions log. The user can confirm whether an 8th is needed.

### Anti-Patterns to Avoid

- **State setup inside timed loop:** `ModelFactory::createFD(nDOF)` is O(n) with heap allocation — must be outside `for (auto _ : state)`. Setup inside loop would measure factory time, not solver time.
- **Forgetting to consume result with `DoNotOptimize`:** With `-O3 -flto`, the compiler can elide the entire `computeAccelerations` call if the result (qddot) is not observed. All benchmark loops must end with `benchmark::DoNotOptimize(result)`.
- **Using `randomTorques()` inside the timed loop:** Allocating a `std::vector<double>` (and copying nDOF values) on every iteration is measurable overhead. Use `fillTorques(data, nDOF)` instead (D-19).
- **Passing raw `double[]` buffer to `computeTorques` or `computeAccelerations`:** The solver APIs take `const Eigen::VectorXd&`. Pre-allocate `Eigen::VectorXd` before the loop and use `fill*(tau.data(), nDOF)` to write into Eigen's internal storage. Eigen::VectorXd uses contiguous storage — `data()` returns the underlying `double*`.
- **Multiple `RandomState` instances across benchmarks:** Each `RandomState` uses seed 42. Creating multiple instances in a single run produces identical sequences, reducing coverage. Use one instance shared across all setup calls, or accept seed-based determinism.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Timing loop for microbenchmarks | `std::chrono::high_resolution_clock` | Google Benchmark's `for (auto _ : state)` | Handles warm-up, iteration count selection, statistical analysis, outlier detection, and CSV/JSON output. Already integrated. |
| Heap allocation avoidance in timed loops | Manual buffer management with `std::vector::reserve` | Pre-allocated `Eigen::VectorXd` with fill methods | Eigen's `data()` gives a `double*` for zero-alloc refill. `std::vector::reserve` still incurs capacity overhead; `Eigen::VectorXd` is the API-matched type for solver inputs. |
| Dead-code prevention | Writing results to `volatile` or `std::cout` | `benchmark::DoNotOptimize()` + `benchmark::ClobberMemory()` | `volatile` prevents inlining; I/O pollutes timing. Google Benchmark's utilities are designed for minimal overhead. |

**Key insight:** The entire Phase 17 is about filling benchmark bodies — every problem (timing loops, dead-code prevention, parameterized registration) has an existing Google Benchmark solution. The only new code is the solver invocation patterns and RandomState extensions.

## Common Pitfalls

### Pitfall 1: ABA solver state is fully recomputed — but RNEA solver state is NOT reset
**What goes wrong:** The RNEA benchmark reuses the same solver object across iterations. `computeTorques` writes to `links[i].v`, `links[i].a`, and `links[i].qddot`. Since the outward pass recomputes v and a from scratch each time (`InverseDynamics.cpp:19-53`), this is safe. However, the benchmark's `qdot` and `qddot` values are set once at setup — the per-iteration `computeTorques(qddot_buffer)` also SETS `links[i].qddot = qddot[i]` (`InverseDynamics.cpp:121-123`). This means the solver's state is fully determined by the input each iteration — no stale state carries over.

**Why it's still a concern:** The solver's `f` (force accumulation) array in `inwardPass()` is allocated as `std::vector<ForceVector> f(links.size())` each call (`InverseDynamics.cpp:64-65`). This IS a heap allocation inside the timed loop. It's inherent to the current RNEA implementation and cannot be eliminated in this phase.

**How to avoid:** Accept the heap allocation in RNEA as the current baseline. Do NOT attempt to pre-allocate and reuse — that would require changing the solver API, which is beyond Phase 17 scope.

**Warning signs:** RNEA benchmark times are slower than ABA for the same DOF, dominated by the per-call heap allocation of `f`. This is expected.

### Pitfall 2: `DoNotOptimize` placement for `SpatialVector` and `Eigen::Matrix` types
**What goes wrong:** `benchmark::DoNotOptimize(x)` works by writing `x` to memory and reading it back. For Eigen types, this triggers the expression template machinery, potentially causing extra evaluation.

**Why it happens:** `SpatialVector`, `MotionVector`, `ForceVector`, `RigidBodyInertia`, and `ArticulatedBodyInertia` are 48+ byte objects (each contains two `Vector3d` or three 3×3 matrices). Writing these to memory via `DoNotOptimize` is non-trivial.

**How to avoid:** 
- For `SpatialVector`/`MotionVector`/`ForceVector`: Use `benchmark::DoNotOptimize(result.getAngular()[0])` — reading a single scalar is sufficient to prevent the compiler from eliding the operation.
- For scalar results (qddot, tau[i]): Use `benchmark::DoNotOptimize(result)` directly.
- For compile-time constant-sized results, `DoNotOptimize` on any sub-field works.

**Warning signs:** Benchmark times that are suspiciously fast for complex operations.

### Pitfall 3: `nDOF` parameter for core microbenchmarks is meaningless
**What goes wrong:** The existing `bench_all.cpp` registers core benchmarks with DOF sweep n=1..20, but individual Plücker transform and cross product operations are O(1) — they don't depend on chain length.

**Why it exists:** The original stub registration treated all benchmarks uniformly. With D-18's expansion to 20 individual operations, the DOF sweep for core becomes 20 identical benchmarks at different n values.

**How to avoid:** Register core microbenchmarks WITHOUT DOF sweep (single registration per operation, no nDOF parameter). The new benchmark functions use `(benchmark::State&)` signature (no `int nDOF` parameter).

**Warning signs:** CSV output showing identical times for "BM_TransformMotion/1DOF" through "BM_TransformMotion/20DOF".

### Pitfall 4: `fill*` methods must not change RNG distribution state
**What goes wrong:** If the `fill` method creates a distribution differently from the returning method (e.g., using different constructor arguments, or a different engine state), the two methods produce different sequences for the same n.

**Why it matters:** Reproducibility. The benchmarks must produce identical random sequences regardless of whether they use the `randomPositions()` or `fillPositions()` path.

**How to avoid:** The fill methods must create an equivalent distribution object (`std::uniform_real_distribution<double>(-range, +range)`) and consume from the same engine (`rng_`). Since distributions are stateless in practice (though they can be in the standard), recreating them per call is safe and produces identical sequences.

**Warning signs:** Different benchmark runs producing different results for the same seed.

### Pitfall 5: Eigen `VectorXd` buffer sizing before fill
**What goes wrong:** Calling `fillTorques(tau.data(), nDOF)` on a `VectorXd` that was sized differently than `nDOF` produces a buffer overflow (if too small) or leaves trailing stale data (if too large — harmless for fill but wastes memory).

**Why it happens:** `Eigen::VectorXd tau` default-constructs to 0 size. `tau.data()` returns `nullptr`. Calling `fill*` on a nullptr crashes.

**How to avoid:** Always pre-allocate: `Eigen::VectorXd tau(nDOF);` before the timed loop. This ensures `tau.data()` points to `nDOF * sizeof(double)` bytes of valid storage.

**Warning signs:** Segfault or memory corruption in benchmark functions.

## Code Examples

### Common Operation 1: Pre-allocated buffer fill loop (ABA)

```cpp
#include <benchmark/benchmark.h>
#include <Eigen/Dense>

// Setup
ForwardDynamics fd = factory.createFD(nDOF);
// ... apply q/qdot ...
Eigen::VectorXd tau(nDOF);  // pre-allocated, sized nDOF

for (auto _ : state) {
    rng.fillTorques(tau.data(), nDOF);  // writes nDOF doubles through data()
    fd.computeAccelerations(tau);
    benchmark::DoNotOptimize(fd.links[0].qddot);
}
```

### Common Operation 2: Pre-allocated buffer fill loop (RNEA)

```cpp
// Setup
InverseDynamics id = factory.createID(nDOF);
// ... apply q/qdot ...
Eigen::VectorXd qddot(nDOF);  // pre-allocated

for (auto _ : state) {
    rng.fillAccelerations(qddot.data(), nDOF);  // writes through data()
    Eigen::VectorXd tau = id.computeTorques(qddot);
    benchmark::DoNotOptimize(tau[0]);
}
```

### Common Operation 3: Zero-alloc fill method implementation

```cpp
void RandomState::fillTorques(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < n; ++i) {
        buf[i] = dist(rng_);
    }
}
```

### Common Operation 4: ABA/RNEA joint state apply pattern

```cpp
void applyJointState(ForwardDynamics& fd,
                     const std::vector<double>& q,
                     const std::vector<double>& qdot) {
    for (size_t i = 0; i < fd.links.size(); ++i) {
        fd.links[i].q = q[i];
        fd.links[i].qdot = qdot[i];
    }
}

// Usage in benchmark setup:
ForwardDynamics fd = factory.createFD(nDOF);
std::vector<double> q = rng.randomPositions(nDOF);
std::vector<double> qdot = rng.randomVelocities(nDOF);
applyJointState(fd, q, qdot);
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Stub benchmark functions (Phase 16) | Real timing loops with per-iteration state regeneration | Phase 17 | Benchmarks actually measure solver performance |
| Two core placeholder benchmarks | 20 individual microbenchmarks | Phase 17 | Enables per-operation profiling and filtering |
| Returning RandomState methods only | Zero-alloc fill variants added | Phase 17 | Eliminates heap allocation overhead in timed loops |
| QDOT generator missing | randomTorques + randomAccelerations added | Phase 17 | Torque and acceleration seeded RNG |

**Deprecated/outdated:**
- The old stub functions (`(void)state; (void)nDOF;`) — these are replaced entirely, not modified.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `Eigen::VectorXd::data()` returns a contiguous `double*` that can be safely written via `fill*` | Pattern 1/2 | LOW — Eigen guarantees contiguous storage for vectors. Only risk is if the VectorXd was default-constructed (size 0, data() = nullptr). Pre-allocation prevents this. |
| A2 | `computeAccelerations` fully recomputes solver state from q/qdot each call — no stale state from prior iterations | Pattern 1 | HIGH — verified against ForwardDynamics.cpp. outwardPass writes v/c from scratch. inwardPass writes Ia/pa from scratch. This is core to the correctness of the benchmark. |
| A3 | `computeTorques` fully recomputes from q/qdot/qddot each call | Pattern 2 | HIGH — verified against InverseDynamics.cpp. outwardPass writes v/a from scratch. inwardPass writes f/tau from scratch. `links[i].qddot` is set from input each call. |
| A4 | `fillTorques(tau.data(), nDOF)` writes to Eigen's internal storage correctly (Eigen reallocates only on assignment, not on data() write) | Pitfall 5 | LOW — Eigen::VectorXd::data() returns a pointer to internal storage. Writing doubles through it is safe as long as the VectorXd is pre-sized. |
| A5 | The 8th inertia operation does not exist distinctly | Pattern 3 | MEDIUM — D-18 says 8, codebase reveals 7. If user intended all 8, a missing operation (like L::multiplySymmetric standalone benchmark, or cross product inertia result) needs identification. |
| A6 | Phase 16's `bench_all.cpp` is the current version (registering 4 benchmark families with DOF sweep) | Architecture | MEDIUM — Phase 17 modifies this file. The exact current content determines the diff. The stub file has been read and confirmed. |

## Open Questions

1. **8th inertia operation ambiguity**
   - What we know: D-18 specifies 8 inertia operations. Codebase reveals 7 identifiable methods: RBI::apply, RBI::operator+, RBI::operator\*, ABI::apply, ABI::operator+(ABI), ABI::operator+(RBI), ABI::operator\*.
   - What's unclear: Is the 8th a missing operation (ABI dot product? RBI cross product?), or is the count aspirational?
   - Recommendation: Implement 7 and flag the discrepancy. The planner should confirm with user whether an 8th exists.

2. **Core benchmark function signatures**
   - What we know: Current bench_all.cpp registers `BM_PluckerTransform(benchmark::State&, int nDOF)` with DOF sweep. D-18 wants 20 individual operations.
   - What's unclear: Should we change function signatures to remove unused nDOF parameter, or keep it for interface consistency?
   - Recommendation: Change to `(benchmark::State&)` (no nDOF) for core benchmarks. The parameter is meaningless for O(1) operations. This requires updating bench_all.cpp registration to not pass nDOF arg.

3. **Per-iteration state regeneration scope**
   - What we know: D-19 specifies fill methods for zero-alloc state regeneration.
   - What's unclear: Should we also regenerate q/qdot each iteration (D-19 lists fillPositions/fillVelocities too), or just tau/qddot?
   - Recommendation: For ABA/RNEA, only tau/qddot needs per-iteration regeneration. The solver recomputes all state from q/qdot each call. Regenerating q/qdot would add noise without exercising different code paths (the outward pass handles any q/qdot the same way). But D-19 explicitly says fillPositions/fillVelocities methods are needed — the planner should decide whether to call them per-iteration or only during setup.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| CMake | Build system | ✓ | 3.19+ | — |
| C++17 compiler | Google Benchmark build | ✓ | g++/Clang | — |
| Eigen3 | Solver objects via ModelFactory | ✓ | 3.4.x | — |
| Google Benchmark v1.9.5 | Benchmarking framework | ✓ (FetchContent) | v1.9.5 | — |

**Missing dependencies with no fallback:** none
**Missing dependencies with fallback:** none

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Google Benchmark v1.9.5 |
| Config file | `benchmarks/CMakeLists.txt` |
| Quick build command | `cmake -B build -DSA_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build --target bench_all` |
| Full suite command | Same build + `build/benchmarks/bench_all` |
| Filtered smoke test | `build/benchmarks/bench_all --benchmark_filter="BM_ABA.*1DOF|BM_RNEA.*1DOF|BM_TransformMotion"` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| BENCH-01 | ABA DOF sweep n=1..20 runs without crash | Smoke/run | `build/benchmarks/bench_all --benchmark_filter="BM_ABA" --benchmark_min_time=0.1` | ✅ bench_aba_stub.cpp → replaced |
| BENCH-02 | RNEA DOF sweep n=1..20 runs without crash | Smoke/run | `build/benchmarks/bench_all --benchmark_filter="BM_RNEA" --benchmark_min_time=0.1` | ✅ bench_rnea_stub.cpp → replaced |
| BENCH-03 | All 20 core microbenchmarks run without crash | Smoke/run | `build/benchmarks/bench_all --benchmark_filter="BM_TransformMotion|BM_TransformForce|BM_InverseTransformMotion|BM_InverseTransformForce|BM_TformRBI|BM_InvtformRBI|BM_TformABI|BM_InvtformABI|BM_CrossMvMv|BM_CrossMvFv|BM_CrossFvMv|BM_CrossFvFv|BM_RBIApply|BM_RBIAdd|BM_RBIScale|BM_ABIApply|BM_ABIAdd|BM_ABIAddRBI|BM_ABIScale" --benchmark_min_time=0.1` | ✅ bench_core_stub.cpp → replaced |
| D-19 | Fill methods produce correct values | Unit (build-time) | Compile check — fill methods used in benchmark loops | N/A (no separate test) |

**Verification procedure:**
```bash
# 1. Build with benchmarks
cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_bench --target bench_all 2>&1

# 2. Quick smoke test (filtered, short run)
build_bench/benchmarks/bench_all \
  --benchmark_filter="BM_ABA.*1DOF|BM_RNEA.*1DOF|BM_TransformMotion" \
  --benchmark_min_time=0.1

# 3. Full suite (all benchmarks, default timing)
build_bench/benchmarks/bench_all

# 4. CSV export for validation
build_bench/benchmarks/bench_all --benchmark_format=csv > benchmark_results.csv
```

### Sampling Rate
- **Per task commit:** `cmake --build build_bench --target bench_all 2>&1 | tail -5` (compiles only)
- **Per wave merge:** Quick smoke test (filtered, short min_time)
- **Phase gate:** Full `bench_all` suite runs to completion with exit code 0

### Wave 0 Gaps
- [ ] Smoke-test command must be verified: `build/benchmarks/bench_all --benchmark_filter="BM_ABA.*1DOF" --benchmark_min_time=0.1` returns nonzero results
- [ ] Need to confirm `bench_all.cpp` compiles after replacing registrations (forward declarations must match new function signatures)

## Security Domain

> Phase 17 has zero security surface. It adds no network I/O, no user input parsing, no data serialization, no privilege elevation, and no new external dependencies. All operations are CPU-bound numerical benchmarks in the same process.

### Applicable ASVS Categories
N/A — no authentication, session management, access control, input validation, or cryptography applies.

### Known Threat Patterns
N/A — no code paths that accept untrusted input or perform security-relevant operations.

## Sources

### Primary (HIGH confidence)
- [VERIFIED: `include/ForwardDynamics.h:156`] — `computeAccelerations(tau, gravity)` API signature
- [VERIFIED: `src/ForwardDynamics.cpp:57-80`] — outwardPass: full recomputation from q/qdot each call
- [VERIFIED: `src/ForwardDynamics.cpp:82-180`] — inwardPass: full recomputation of Ia/pa each call
- [VERIFIED: `include/InverseDynamics.h:144`] — `computeTorques(qddot, gravity)` API signature
- [VERIFIED: `src/InverseDynamics.cpp:19-53`] — outwardPass: full recomputation of v/a each call
- [VERIFIED: `src/InverseDynamics.cpp:55-94`] — inwardPass: heap-allocates f vector each call
- [VERIFIED: `benchmarks/bench_all.cpp:26-52`] — Current registration structure (4 families × 20 DOF)
- [VERIFIED: `benchmarks/common/model_factory.h:68-126`] — createFD/createID API
- [VERIFIED: `benchmarks/common/random_state.h:35-58`] — Current RandomState API (positions/velocities only)
- [VERIFIED: `include/PluckerTransform.h:77-216`] — All 8 transform operation signatures
- [VERIFIED: `include/SpatialUtils.h:78-126`] — All 4 cross product free function signatures
- [VERIFIED: `include/RigidBodyInertia.h:27-122`] — RBI::apply, operator+, operator* signatures (3 operations)
- [VERIFIED: `include/ArticulatedBodyInertia.h:77-220`] — ABI::apply, operator+(ABI), operator+(RBI), operator* signatures (4 operations)
- [VERIFIED: `benchmarks/core/bench_core_stub.cpp:21-43`] — Current core stub (2 placeholder functions only)
- [VERIFIED: `include/SpatialVector.h:67-178`] — SpatialVector::data returns double* via Eigen (contiguous storage)

### Secondary (MEDIUM confidence)
- [CITED: google.github.io/benchmark/user_guide.html] — DoNotOptimize/ClobberMemory patterns, custom main() pattern
- [CITED: cppreference.com — std::mt19937, std::uniform_real_distribution] — RNG distribution creation

### Tertiary (LOW confidence)
- [ASSUMED] The 8th inertia operation in D-18 is a count discrepancy — 7 identifiable operations, user may have intended 8
- [ASSUMED] `DoNotOptimize(fd.links[0].qddot)` is sufficient to prevent elision of computeAccelerations

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — Google Benchmark v1.9.5 already integrated. No new packages.
- Architecture: HIGH — All patterns derived from reading actual solver source (ForwardDynamics.cpp, InverseDynamics.cpp, PluckerTransform.cpp) and existing benchmark infrastructure.
- Pitfalls: HIGH — Verified against actual code paths. ABA solver recomputes fully. RNEA has heap allocation in inwardPass. Buffer sizing for fill* methods.
- RandomState extension: HIGH — Direct extension of existing class. Fill variants are straightforward.
- DOF sweep for core: MEDIUM — Current bench_all.cpp registers core with DOF. D-18 expansion requires restructuring. The exact number of inertia operations (7 vs 8) needs user confirmation.

**Research date:** 2026-06-05
**Valid until:** 2026-07-05 (30 days — stable dependencies)
