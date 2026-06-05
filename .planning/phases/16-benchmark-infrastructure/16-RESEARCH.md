# Phase 16: Benchmark Infrastructure - Research

**Researched:** 2026-06-05
**Domain:** C++ benchmark build system integration and shared utilities
**Confidence:** HIGH

## Summary

This phase integrates Google Benchmark v1.9.5 into the SpatialAlgebra build system via FetchContent (reusing the established GTest pattern), creates the `benchmarks/` directory structure, and builds shared utilities (model factory, random state generator) that Phases 17-18 will consume. Everything is guarded behind `SA_BUILD_BENCHMARKS` (default OFF) so benchmarks have zero impact on normal builds.

The key technical decisions from context: single unified `bench_all` executable with programmatic `RegisterBenchmark` calls (not `BENCHMARK_MAIN()` macro), FetchContent scoped to `benchmarks/CMakeLists.txt`, and a unified model factory API returning both `ForwardDynamics` and `InverseDynamics` solver objects for arbitrary n-DOF chains. All benchmark source files are created in this phase as stubs — Phases 17-18 fill in the actual benchmark logic.

**Primary recommendation:** Follow the exact FetchContent pattern from `CMakeLists.txt:17-23` (GTest declarations) for Google Benchmark, scoped to `benchmarks/CMakeLists.txt`. Use `RegisterBenchmark` in a custom `main()` for the unified executable. Place the model factory in `benchmarks/common/` since it is shared across all benchmark domains.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| BINF-01 | Google Benchmark v1.9.5 integrated via FetchContent | Verified v1.9.5 release (2026-01-21) uses Git tag `v1.9.5`. FetchContent pattern established at `CMakeLists.txt:17-23`. Requires `BENCHMARK_ENABLE_TESTING OFF` to avoid GTest dependency cascade. Targets: `benchmark::benchmark` and `benchmark::benchmark_main`. |
| BINF-02 | `benchmarks/` directory with `SA_BUILD_BENCHMARKS` guard (default OFF) | Follow `ENABLE_COVERAGE` pattern at `CMakeLists.txt:42`. Guard wraps `add_subdirectory(benchmarks)` in root CMakeLists.txt. `option(SA_BUILD_BENCHMARKS OFF ...)`. |
| BINF-03 | Shared benchmark utilities (model factory, random state generator) | `benchmarks/common/model_factory.h/.cpp` for unified factory (per D-01/D-02/D-03). `benchmarks/common/random_state.h/.cpp` for deterministic random state (per D-09/D-10/D-11). |

</phase_requirements>

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Unified API returning both ForwardDynamics and InverseDynamics solver objects — caller specifies which solver type to construct
- **D-02:** Configurable joint axis per joint (Z, X, Y, or arbitrary screw axis per link)
- **D-03:** Support both serial chains and branching (Y-shaped) configurations
- **D-04:** Link inertia configurable via parameters (mass range, COM range, inertia tensor range) — allows uniform, random, or explicit per-link specification
- **D-05:** Single unified executable (`bench_all`) with Google Benchmark sub-benchmark registration — benchmarks filterable at runtime via `--benchmark_filter`
- **D-06:** Subdirectories per domain: `benchmarks/aba/`, `benchmarks/rnea/`, `benchmarks/core/` (Plücker, cross-product), `benchmarks/common/` (shared utilities)
- **D-07:** Google Benchmark FetchContent lives in `benchmarks/CMakeLists.txt` — main CMakeLists.txt only adds `add_subdirectory(benchmarks)` guarded by `SA_BUILD_BENCHMARKS`
- **D-08:** Phase 16 creates all benchmark source files (stubs for Phases 17-18 to fill with benchmark logic)
- **D-09:** Joint positions: uniform random in [-π/2, π/2]
- **D-10:** Joint velocities: uniform random in [-5, 5] rad/s
- **D-11:** Fixed random seed 42 for reproducibility across runs
- **D-12:** Both Debug and Release build profiles supported
- **D-13:** Default build type inherits from parent CMake project (no forced Release)
- **D-14:** LTO (`-flto`) enabled for benchmark targets
- **D-15:** Guard option `SA_BUILD_BENCHMARKS` (default OFF) in top-level CMakeLists.txt

### the agent's Discretion
- Specific benchmark executable name (within `bench_all` convention)
- Exact subdirectory structure within each domain directory
- CMake minimum version for `benchmarks/CMakeLists.txt` (inherit from parent)
- Google Benchmark version pin (v1.9.5 per requirements) and FetchContent URL details
- Shared utility function signatures and header file organization

### Deferred Ideas (OUT OF SCOPE)
None — discussion stayed within phase scope.

</user_constraints>

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Google Benchmark FetchContent integration | Build System (CMake) | — | Dependency management at build configure time; no runtime component |
| `SA_BUILD_BENCHMARKS` guard option | Build System (CMake) | — | Standard CMake `option()` pattern, scoped to root CMakeLists.txt |
| LTO optimization for benchmark targets | Build System (CMake) | — | `target_compile_options` + `target_link_options` applied at target level |
| Model factory (unified API) | Benchmarks Common (utility) | — | Pure C++ production code lives in `benchmarks/common/`; consumed by all benchmark domains |
| Random state generator | Benchmarks Common (utility) | — | Pure C++ utility in `benchmarks/common/`; consumed by all benchmark domains |
| `bench_all` executable registration | Build System (CMake) | Benchmarks Common (main) | CMake target links all domain source files; custom `main()` uses `RegisterBenchmark` |
| ABA benchmark stubs | Benchmarks ABA | — | Source files in `benchmarks/aba/`, filled in Phase 17 |
| RNEA benchmark stubs | Benchmarks RNEA | — | Source files in `benchmarks/rnea/`, filled in Phase 17 |
| Core microbenchmark stubs | Benchmarks Core | — | Source files in `benchmarks/core/`, filled in Phase 17 |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Google Benchmark | v1.9.5 | C++ microbenchmarking framework | Industry standard for C++ benchmarks; statistical analysis, warm-up cycles, iteration counting, parameterized registration. FetchContent integration matches existing GTest pattern. [VERIFIED: github.com/google/benchmark tag v1.9.5] |
| Eigen3 | 3.4...5 (existing) | Linear algebra backend used in factory | Factory constructs `PluckerTransform`, `Rotation`, `MotionVector`, `ForceVector` which all depend on Eigen3 types. Already in project. |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| C++ `<random>` | C++17 std | Deterministic random state generation | Used for generating random joint positions/velocities with fixed seed 42 (D-09 through D-11). `std::mt19937` + `std::uniform_real_distribution`. The `<random>` header is standard library — no external dependency. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `RegisterBenchmark` (programmatic) | `BENCHMARK()` macro (static) | Macro registration requires file-scope static initializers; programmatic `RegisterBenchmark` enables DOF-swept registration (n=1..20) and return-value chaining for args. D-05 mandates `--benchmark_filter` support, which both approaches provide, but programmatic registration is the only way to parameterize the DOF sweep. |
| FetchContent in `benchmarks/CMakeLists.txt` | FetchContent in root `CMakeLists.txt` | D-07 locks this decision. Scoping to benchmarks/CMakeLists.txt keeps the root CMakeLists.txt clean and means benchmarks deps are only fetched when `SA_BUILD_BENCHMARKS=ON`. |

**Installation:**
```bash
# No manual installation needed — FetchContent handles everything.
# Google Benchmark v1.9.5 is auto-downloaded and built when:
cmake -B build -DSA_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
```

**Version verification:**
```bash
# Google Benchmark v1.9.5 release verified: 2026-01-21
# GitHub tag: v1.9.5
# Repo: https://github.com/google/benchmark
# Requires C++17 to build (matches project's C++17 requirement)
# Targets: benchmark::benchmark, benchmark::benchmark_main
```

## Package Legitimacy Audit

> Phase 16 has no external packages beyond those introduced via FetchContent. Google Benchmark v1.9.5 is a well-established project (10k+ stars, 340+ contributors, first release 2013). No package manager installs are needed.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| Google Benchmark (FetchContent) | GitHub | 12+ yrs | 10k+ stars | github.com/google/benchmark | N/A (verified OSS) | Approved — v1.9.5 tag verified on GitHub |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

## Architecture Patterns

### System Architecture Diagram

```
cmake -B build -DSA_BUILD_BENCHMARKS=ON
         │
         ▼
┌──────────────────────┐
│   CMakeLists.txt      │
│   (root)              │
│                       │
│  option(SA_BUILD_     │──── D-15: Guard (default OFF)
│   BENCHMARKS OFF)     │
│                       │
│  if(SA_BUILD_         │
│   BENCHMARKS)         │──── D-07: Delegates to benchmarks/
│    add_subdirectory(  │
│     benchmarks)       │
│  endif()              │
└──────────────────────┘
         │
         ▼
┌──────────────────────────────────────────────┐
│  benchmarks/CMakeLists.txt                    │
│                                               │
│  include(FetchContent)                        │
│  set(BENCHMARK_ENABLE_TESTING OFF)            │
│  set(BENCHMARK_ENABLE_INSTALL OFF)            │
│  FetchContent_Declare(googlebenchmark         │──── D-07: FetchContent lives here
│    GIT_REPOSITORY ...                         │
│    GIT_TAG v1.9.5)                            │─── BINF-01: v1.9.5 pinned
│  FetchContent_MakeAvailable(googlebenchmark)  │
│                                               │
│  add_subdirectory(common)                     │─── Shared utilities
│  add_subdirectory(aba)                        │─── Stubs for Phase 17
│  add_subdirectory(rnea)                       │─── Stubs for Phase 17
│  add_subdirectory(core)                       │─── Stubs for Phase 17
│                                               │
│  add_executable(bench_all ...)                │─── D-05: Single unified executable
│  target_link_libraries(bench_all              │
│    SpatialAlgebra                             │─── Links against libSpatialAlgebra.a
│    benchmark::benchmark                       │
│    Eigen3::Eigen)                             │
│  target_compile_options(bench_all             │─── D-14: LTO enabled
│    PRIVATE -O3 -DNDEBUG -flto)               │
└──────────────────────────────────────────────┘
         │
         ▼
┌──────────────────────────────────────────────────────┐
│  bench_all (single executable)                        │
│                                                       │
│  main() {                                             │
│    RegisterBenchmark("BM_ABA_ForwardDynamics",        │── D-05: Programmatic registration
│                      BM_ABA_ForwardDynamics, nDOF)    │    for DOF sweep (n=1..20)
│    RegisterBenchmark("BM_RNEA_InverseDynamics",       │
│                      BM_RNEA_InverseDynamics, nDOF)   │
│    RegisterBenchmark("BM_PluckerTransform", ...)      │
│    RegisterBenchmark("BM_CrossProduct", ...)          │
│                                                       │
│    benchmark::Initialize(&argc, argv);                │
│    benchmark::RunSpecifiedBenchmarks();               │
│    benchmark::Shutdown();                             │
│  }                                                     │
│                                                       │
│  # Runtime filtering:                                  │
│  $ ./bench_all --benchmark_filter="ABA"               │── Filters by name substring
│  $ ./bench_all --benchmark_filter="BM_RNEA.*6"       │── RNEA with 6 DOF
│  $ ./bench_all --benchmark_format=csv                 │── CSV output
└──────────────────────────────────────────────────────┘

Data Flow for a single benchmark run:

  bench_all/main()
       │
       ▼
  ┌──────────────┐     ┌──────────────────────┐
  │ RandomState  │────▶│ ModelFactory          │
  │  (seed 42)   │     │  .createFD(nDOF, ...) │
  │              │     │  .createID(nDOF, ...) │
  │ q  ∈ [-π/2,  │     │                      │
  │     π/2]     │     │ Returns populated     │
  │ qdot ∈ [-5,  │     │ ForwardDynamics or    │
  │     5]       │     │ InverseDynamics obj   │
  └──────────────┘     └──────────────────────┘
                               │
                               ▼
                      ┌──────────────────┐
                      │ Benchmark Loop    │
                      │ for (auto _ :     │
                      │      state) {     │
                      │   solver->compute │
                      │   ...             │
                      │ }                 │
                      └──────────────────┘
```

### Recommended Project Structure
```
benchmarks/
├── CMakeLists.txt              # FetchContent + top-level build rules + bench_all
├── bench_all.cpp               # main() with RegisterBenchmark registrations
├── common/
│   ├── CMakeLists.txt          # Build shared_bench_utils static library
│   ├── model_factory.h         # Unified factory API declaration
│   ├── model_factory.cpp       # Factory implementation
│   ├── random_state.h          # Deterministic random state generator
│   └── random_state.cpp        # Random state implementation
├── aba/
│   ├── CMakeLists.txt          # Add sources to bench_all (or library)
│   └── bench_aba.cpp           # ABA benchmark (stub for Phase 17)
├── rnea/
│   ├── CMakeLists.txt          # Add sources to bench_all (or library)
│   └── bench_rnea.cpp          # RNEA benchmark (stub for Phase 17)
└── core/
    ├── CMakeLists.txt          # Add sources to bench_all (or library)
    └── bench_core.cpp          # Core microbenchmarks (stub for Phase 17)
```

### Pattern 1: FetchContent for Google Benchmark
**What:** Reuse the established GTest FetchContent pattern from `CMakeLists.txt:17-23` for Google Benchmark, scoped to `benchmarks/CMakeLists.txt`. Must disable benchmark's own test/install targets to avoid redundant GTest dependency.

**When to use:** This is the only dependency integration pattern — Google Benchmark always via FetchContent in `benchmarks/CMakeLists.txt` (D-07).

**Example:**
```cmake
# benchmarks/CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(SpatialAlgebraBenchmarks)

# Suppress benchmark internal test/install targets
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
)
FetchContent_MakeAvailable(googlebenchmark)

# ... add_subdirectory calls, bench_all target ...
```
Source: [VERIFIED: google/benchmark README.md] — pattern confirmed in official docs.

### Pattern 2: Programmatic Benchmark Registration with DOF Sweep
**What:** Use `benchmark::RegisterBenchmark(name, fn, args...)` to register benchmarks programmatically in `main()`, enabling DOF-parameterized registration (n=1..20). This replaces the `BENCHMARK_MAIN()` macro pattern.

**When to use:** Required whenever the benchmark argument (e.g., DOF count) must vary programmatically. D-05 mandates this approach for the unified `bench_all` executable.

**Example:**
```cpp
// bench_all.cpp
#include <benchmark/benchmark.h>
#include "common/model_factory.h"
#include "common/random_state.h"

// Forward declaration of benchmark function (defined in aba/bench_aba.cpp)
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF);

int main(int argc, char** argv) {
    // Register ABA benchmarks for DOF sweep n=1..20
    for (int n = 1; n <= 20; ++n) {
        benchmark::RegisterBenchmark(
            ("BM_ABA_ForwardDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_ABA_ForwardDynamics, n
        );
    }
    // Same pattern for RNEA, core microbenchmarks...

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
```
Source: [CITED: google.github.io/benchmark/user_guide.html §Using RegisterBenchmark]

### Pattern 3: Unified Model Factory
**What:** A factory class in `benchmarks/common/` that constructs arbitrary n-DOF kinematic chains. Returns populated `ForwardDynamics` or `InverseDynamics` solver objects based on a template parameter or enum.

**When to use:** Every benchmark that needs a kinematic chain must go through the factory — never construct solver objects manually in benchmark code.

**Example:**
```cpp
// benchmarks/common/model_factory.h
#pragma once
#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include <Eigen/Dense>

namespace SpatialAlgebra::Bench {

// Joint configuration for a single link
struct JointConfig {
    MotionVector axis;     // Screw axis (default: Z revolute)
    Vector3d translation;  // Transform from parent (default: along X)
    double mass = 1.0;
    Vector3d com = Vector3d::Zero();
};

class ModelFactory {
public:
    // Create serial chain of nDOF links, all with same JointConfig
    ForwardDynamics createFD(int nDOF, const JointConfig& cfg = JointConfig());
    InverseDynamics createID(int nDOF, const JointConfig& cfg = JointConfig());
    
    // For branching Y-shaped: nDOF total, split at branchPoint
    ForwardDynamics createFDBranching(int nDOF, int branchPoint,
                                       const JointConfig& cfg = JointConfig());
    InverseDynamics createIDBranching(int nDOF, int branchPoint,
                                       const JointConfig& cfg = JointConfig());
    
    // Configure inertia per-link for random/uniform/explicit settings
    // ... mass range, COM range, inertia tensor range
};

}
```

### Pattern 4: Deterministic Random State Generator
**What:** A utility that generates random joint positions and velocities using a fixed-seed Mersenne Twister.

**When to use:** Every benchmark uses this to populate joint states before timing the solver.

**Example:**
```cpp
// benchmarks/common/random_state.h
#pragma once
#include <random>
#include <vector>

namespace SpatialAlgebra::Bench {

class RandomState {
public:
    RandomState() : rng_(42) {}  // D-11: Fixed seed
    
    // Generate nDOF random positions uniform in [-π/2, π/2]
    std::vector<double> randomPositions(int nDOF);
    
    // Generate nDOF random velocities uniform in [-5, 5]
    std::vector<double> randomVelocities(int nDOF);
    
private:
    std::mt19937 rng_;
};

}
```
Source: [CITED: cppreference.com — std::mt19937, std::uniform_real_distribution]

### Anti-Patterns to Avoid
- **Using `BENCHMARK_MAIN()` macro in `bench_all`:** The `BENCHMARK_MAIN()` macro generates a `main()` that only runs benchmarks registered via `BENCHMARK()` macros. Since benchmarks are in separate translation units and registered programmatically, a custom `main()` with `benchmark::Initialize`/`RunSpecifiedBenchmarks`/`Shutdown` is required.
- **Forcing Release build type:** D-13 explicitly says inherit build type from parent. `-O3 -DNDEBUG` is applied via `target_compile_options` on the `bench_all` target, not by overriding `CMAKE_BUILD_TYPE`.
- **Passing compile flags via `CMAKE_CXX_FLAGS`:** Add options via `target_compile_options(bench_all PRIVATE ...)` to limit scope to benchmark targets only.
- **Placing FetchContent in root CMakeLists.txt:** D-07 scopes it to `benchmarks/CMakeLists.txt`, keeping root clean.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Microbenchmark timing, statistics, warm-up | Custom timer with `std::chrono` | Google Benchmark | Google Benchmark handles statistical analysis (mean, median, stddev), automatic iteration count, warm-up cycles, CPU vs real-time measurement, and CSV/JSON output. Custom timer logic introduces measurement errors from warm-up transients, clock resolution issues, and lacks statistical rigor. |
| Parameterized benchmark registration | Manually written main() with repeated calls | `RegisterBenchmark(name, fn, args...)` | Google Benchmark's programmatic registration handles name generation, argument passing, and `--benchmark_filter` integration. |
| Deterministic RNG | `rand()` or `std::default_random_engine` | `std::mt19937` with explicit seed | `rand()` has poor distribution quality and is not reproducible across platforms. `std::default_random_engine` is implementation-defined. `std::mt19937` with seed 42 gives identical sequences across all standard libraries. |

**Key insight:** Google Benchmark eliminates an entire class of measurement pitfalls (compiler optimizations of unused results, warm-up transients, clock resolution aliasing, statistical noise). The library has 10+ years of battle testing. Custom timing frameworks are universally worse.

## Common Pitfalls

### Pitfall 1: Google Benchmark pulls in GTest dependency via FetchContent
**What goes wrong:** Google Benchmark's own tests depend on GTest. When fetched via FetchContent with default settings, it tries to also fetch GTest, potentially conflicting with the project's GTest setup at `CMakeLists.txt:15-24`.
**Why it happens:** Benchmark's CMakeLists.txt conditionally includes GTest for its own test suite.
**How to avoid:** Always set `BENCHMARK_ENABLE_TESTING OFF` before `FetchContent_MakeAvailable(googlebenchmark)`. Also set `BENCHMARK_ENABLE_INSTALL OFF` to suppress install targets.
**Warning signs:** Build errors about missing GTest, or duplicate GTest target definitions.

### Pitfall 2: Dead-code elimination in Release mode
**What goes wrong:** The compiler optimizes away the entire benchmark loop because the result is unused, giving artificially fast (zero-cost) timings.
**Why it happens:** `-O3 -DNDEBUG` with LTO is aggressive about eliminating dead code.
**How to avoid:** Use Google Benchmark's `state.PauseTiming()`/`state.ResumeTiming()` for setup/teardown, and always consume the result (e.g., `benchmark::DoNotOptimize(result)`, `benchmark::ClobberMemory()`). The standard pattern `for (auto _ : state) { output = fn(input); benchmark::DoNotOptimize(output); }` prevents this.
**Warning signs:** Benchmark time of ~0 ns or times that don't scale with DOF when they should.

### Pitfall 3: LTO flags conflict during Debug builds
**What goes wrong:** Adding `-flto` unconditionally causes linker errors or warnings on Debug builds where LTO is unexpected.
**Why it happens:** LTO is typically Release-only. Some compilers warn or error when mixing LTO objects with non-LTO objects.
**How to avoid:** Conditionally apply `-flto` based on build type, or apply it only to benchmark targets (not the library). Since benchmarks are always compiled at `-O3`, LTO is safe — just document that benchmarks are built at `-O3` regardless of parent build type.
**Warning signs:** Linker warnings about LTO type mismatches, or link failures mentioning `lto`.

### Pitfall 4: Model factory creating expensive objects in the timed loop
**What goes wrong:** Factory construction time leaks into benchmark measurements if the factory is called inside the timed loop.
**Why it happens:** Setup code placed inside `for (auto _ : state)` is measured as part of the benchmark.
**How to avoid:** Always construct solver objects and populate link states OUTSIDE the timing loop. Use `state.PauseTiming()` during setup if setup must happen in the benchmark function scope. The standard pattern: setup before loop, `for (auto _ : state) { timed_operation(); }`, teardown after.
**Warning signs:** Benchmark times that include setup overhead, scaling disproportionately with DOF vs the actual solver.

### Pitfall 5: Google Benchmark v1.9.5 requires C++17 to build
**What goes wrong:** The build fails with C++11/14 errors from benchmark's own source.
**Why it happens:** From the official README: "The library can be used with C++11. However, it requires C++17 to build."
**How to avoid:** The project already requires C++17 (`CMakeLists.txt:4: set(CMAKE_CXX_STANDARD 17)`). No special handling needed — the parent project's standard propagates to the fetched content.
**Warning signs:** Template errors in benchmark headers during build.

### Pitfall 6: `SpatialVector` and type aliases are not available by default in benchmark code
**What goes wrong:** Benchmark code in `benchmarks/` fails to compile because it doesn't include the right SpatialAlgebra headers.
**Why it happens:** The `benchmarks/` directory is separate from `include/`. Benchmarks link against the compiled library but also need header includes for type definitions.
**How to avoid:** The `bench_all` target links against `SpatialAlgebra` (the static library). Include paths are inherited from the parent project via `include_directories(include)` at `CMakeLists.txt:9`. The `benchmarks/common/model_factory.h` should `#include "ForwardDynamics.h"` and `#include "InverseDynamics.h"` directly — the project's top-level include directory is available.
**Warning signs:** "No such file" errors for `SpatialVector.h`, `ForwardDynamics.h`, etc.

## Code Examples

### Common Operation 1: Factory creating a serial chain (ForwardDynamics)
```cpp
// benchmarks/common/model_factory.cpp
#include "model_factory.h"

namespace SpatialAlgebra::Bench {

ForwardDynamics ModelFactory::createFD(int nDOF, const JointConfig& cfg) {
    ForwardDynamics fd;
    fd.links.reserve(nDOF);
    
    for (int i = 0; i < nDOF; ++i) {
        ForwardDynamics::Link link;
        link.parent = i - 1;  // -1 for base, i-1 for children
        
        // Transform: translate along X by cfg.translation
        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );
        
        // Inertia: uniform mass, COM at origin
        link.I = RigidBodyInertia(
            cfg.mass,
            cfg.com,
            lt::Identity(3)  // Identity inertia tensor
        );
        
        // Joint axis: configurable screw axis
        link.S = cfg.axis;
        
        link.q = 0.0;
        link.qdot = 0.0;
        link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
        
        fd.links.push_back(link);
    }
    
    return fd;
}

}  // namespace SpatialAlgebra::Bench
```
### Common Operation 2: Random state generation
```cpp
// benchmarks/common/random_state.cpp
#include "random_state.h"

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
```
### Common Operation 3: Applying random state to solver
```cpp
// Called during benchmark setup (outside timed loop)
void applyRandomState(ForwardDynamics& fd, const std::vector<double>& q,
                      const std::vector<double>& qdot) {
    for (size_t i = 0; i < fd.links.size(); ++i) {
        fd.links[i].q = q[i];
        fd.links[i].qdot = qdot[i];
    }
}
```

### Common Operation 4: FetchContent in benchmarks/CMakeLists.txt
```cmake
# benchmarks/CMakeLists.txt
cmake_minimum_required(VERSION 3.19)
project(SpatialAlgebraBenchmarks LANGUAGES CXX)

# Suppress Google Benchmark's internal test/install targets
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
)
FetchContent_MakeAvailable(googlebenchmark)

# Set C++17 (inherit from parent, but benchmark needs it)
set(CMAKE_CXX_STANDARD 17)

# Shared utilities library
add_subdirectory(common)

# Domain benchmark directories (each has CMakeLists.txt)
add_subdirectory(aba)
add_subdirectory(rnea)
add_subdirectory(core)

# Collect all benchmark sources
set(BENCH_SOURCES
    bench_all.cpp
    ${ABA_SOURCES}
    ${RNEA_SOURCES}
    ${CORE_SOURCES}
)

# Unified benchmark executable
add_executable(bench_all ${BENCH_SOURCES})

target_link_libraries(bench_all
    PRIVATE
    SpatialAlgebra
    benchmark::benchmark
    Eigen3::Eigen
    bench_common  # shared utilities library
)

# Release-level optimization regardless of parent build type (D-12/D-13)
target_compile_options(bench_all PRIVATE -O3 -DNDEBUG)
target_compile_options(bench_all PRIVATE -flto)       # D-14
target_link_options(bench_all PRIVATE -flto)
```
### Common Operation 5: Guard in root CMakeLists.txt
```cmake
# In root CMakeLists.txt, after existing test/example sections:

# Build benchmarks (opt-in, must be explicitly enabled)
option(SA_BUILD_BENCHMARKS "Build performance benchmarks" OFF)
if(SA_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
```
Source: Pattern from `CMakeLists.txt:42` (ENABLE_COVERAGE option pattern).

### Common Operation 6: Stub benchmark source file
```cpp
// benchmarks/aba/bench_aba.cpp (STUB — Phase 17 fills implementation)
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

// Forward declaration for bench_all.cpp registration
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement ABA benchmark with DOF sweep
    // 1. Create model via ModelFactory::createFD(nDOF)
    // 2. Generate random state via RandomState
    // 3. Apply random state to solver
    // 4. For (auto _ : state): computeAccelerations(random_tau)
    // 5. benchmark::DoNotOptimize(result)
    (void)state;
    (void)nDOF;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `BENCHMARK_MAIN()` macro | Custom `main()` with `RegisterBenchmark` | v1.5+ (Google Benchmark evolved over time) | Manual registration enables DOF-swept parameterization — critical for this project |
| Custom timing with `std::chrono` | Google Benchmark FetchContent integration | This phase | Eliminates measurement boilerplate and common timing errors |
| Inline RNG with `rand()` | `std::mt19937` with fixed seed | C++11 (2011) | Deterministic, portable, high-quality random sequences |

**Deprecated/outdated:**
- `BENCHMARK_MAIN()`: Not deprecated, but inappropriate for programmatic DOF-swept registration. Use custom `main()`.
- Manual GTest linking for benchmarks: Google Benchmark v1.9+ provides `benchmark::benchmark_main` for standalone main, but we use custom main.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Google Benchmark v1.9.5 FetchContent URL pattern uses `GIT_REPOSITORY` + `GIT_TAG` | Standard Stack | Low — could alternatively use URL download; both work with FetchContent |
| A2 | `BENCHMARK_ENABLE_TESTING OFF` + `BENCHMARK_ENABLE_INSTALL OFF` suppress all unwanted targets | Pattern 1 | Low — verified in multiple community examples; if wrong, plan adds extra CMake toggle |
| A3 | `-flto` flag works on macOS with Apple Clang | Build Config | Low — Apple Clang supports `-flto`; if not, CMake's `cmake_policy` or `check_cxx_compiler_flag` can gate it |
| A4 | Project's `include_directories(include)` at root propagates to `benchmarks/` subdirectory | Pitfall 6 | Medium — `add_subdirectory` inherits parent scope's `include_directories`; if not, `target_include_directories(bench_all PRIVATE ${CMAKE_SOURCE_DIR}/include)` would be needed |

## Open Questions

1. **Exact subdirectory CMakeLists.txt structure within each domain**
   - What we know: D-06 mandates `benchmarks/aba/`, `benchmarks/rnea/`, `benchmarks/core/`, `benchmarks/common/`
   - What's unclear: Whether each domain subdirectory has its own internal CMakeLists.txt (adding sources to a parent variable) or sources are listed centrally in `benchmarks/CMakeLists.txt`
   - Recommendation: Each subdirectory should have its own CMakeLists.txt that appends sources to a parent-scope list variable (e.g., `set(ABA_SOURCES bench_aba.cpp PARENT_SCOPE)`), matching the `examples/CMakeLists.txt` pattern.

2. **LTO flag compatibility across compilers**
   - What we know: `-flto` works with GCC and Clang; D-14 mandates LTO
   - What's unclear: Whether MSVC or other compilers need `/GL` / `/LTCG` instead
   - Recommendation: Use `-flto` (sufficient for macOS with Apple Clang and Linux with GCC). If cross-platform support is added later, use `check_cxx_compiler_flag` to probe.

3. **bench_common library vs direct source inclusion**
   - What we know: Shared utilities must be available to all benchmark domains
   - What's unclear: Whether `bench_common` should be a STATIC library (requiring separate CMake target) or sources included directly in `bench_all`
   - Recommendation: Use a STATIC library `bench_common` (via `add_library(bench_common STATIC ...)`) in `benchmarks/common/CMakeLists.txt`. This matches the library approach used for `SpatialAlgebra` itself and avoids ODR issues if shared utilities grow.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| CMake | FetchContent for Google Benchmark | ✓ | 3.19+ (project minimum) | — |
| C++17 compiler | Google Benchmark build requirement | ✓ | g++/Clang via project build | — |
| Git | FetchContent GIT_REPOSITORY fetch | ✓ | — | — |
| Eigen3 | Model factory (via SpatialAlgebra) | ✓ | 3.4.x (installed) | — |

**Missing dependencies with no fallback:** none
**Missing dependencies with fallback:** none

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Google Test v1.12.1 (via FetchContent), Google Benchmark v1.9.5 |
| Config file | CMakeLists.txt (root) + CMakeLists.txt (benchmarks/) |
| Quick run command | `cmake --build build && ctest --output-on-failure` |
| Full suite command | same |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| BINF-01 | Google Benchmark FetchContent integration | Build/smoke | `cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON && cmake --build build_bench --target bench_all` | ❌ Wave 0 |
| BINF-02 | benchmarks/ exists with SA_BUILD_BENCHMARKS guard | Build | `cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON && test -f build_bench/benchmarks/bench_all` | ❌ Wave 0 |
| BINF-03 | Factory creates solver objects correctly | Unit | Compile and run `bench_all` (stub functions link and run) | ❌ Wave 0 |

**Note:** Phase 16's benchmarks are stubs — they compile and link but do nothing at runtime. Verification of BINF-03 is a compile-and-link check, not a behavioral test. Actual benchmark execution is verified in Phase 17.

### Sampling Rate
- **Per task commit:** `cmake --build build --target bench_all 2>&1 | tail -5`
- **Per wave merge:** `cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON && cmake --build build_bench --target bench_all`
- **Phase gate:** `build_bench/benchmarks/bench_all` exists and runs with exit code 0

### Wave 0 Gaps
- [ ] Need build verification script that compiles with `-DSA_BUILD_BENCHMARKS=ON` and confirms `bench_all` executable exists
- [ ] Need to verify `SA_BUILD_BENCHMARKS=OFF` (default) does NOT build benchmarks

## Security Domain

> `security_enforcement` is not explicitly set in config.json (absent = enabled by default). However, Phase 16 involves no network I/O, no user input, no data serialization, and no privilege elevation. The security surface is effectively zero — this is a pure infrastructure/build-system phase with no runtime attack surface.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | — |
| V3 Session Management | no | — |
| V4 Access Control | no | — |
| V5 Input Validation | no | — |
| V6 Cryptography | no | — |

### Known Threat Patterns

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| FetchContent MITM (dependency substitution) | Tampering | Google Benchmark v1.9.5 via GitHub HTTPS GIT_REPOSITORY with pinned GIT_TAG. Auto-negotiates TLS on first clone; subsequent builds use local cache. |
| Code execution via FetchContent build | Elevation of Privilege | Google Benchmark is a well-established project (10k+ stars, 340+ contributors, Apache 2.0). Build is sandboxed by CMake's sub-build mechanism. No network access after initial clone. |

**Risk acceptance:** Google Benchmark's build system executes arbitrary CMake code from the fetched repository. This is standard industry practice for all FetchContent dependencies and accepted for well-established OSS projects. The v1.9.5 git tag is pinned for immutability.

## Sources

### Primary (HIGH confidence)
- [VERIFIED: github.com/google/benchmark tag v1.9.5] — Release 2026-01-21, C++17 build requirement, FetchContent pattern, RegisterBenchmark API
- [VERIFIED: CMakeLists.txt:17-23] — Existing FetchContent pattern for GTest (template for Google Benchmark)
- [VERIFIED: CMakeLists.txt:42-46] — Existing `option()` pattern for ENABLE_COVERAGE (template for SA_BUILD_BENCHMARKS)
- [VERIFIED: examples/CMakeLists.txt] — Subdirectory CMakeLists.txt pattern with `find_package(Eigen3)`, executable setup, link targets
- [VERIFIED: ForwardDynamics.h:79-111] — ForwardDynamics::Link struct API for factory output
- [VERIFIED: InverseDynamics.h:75-100] — InverseDynamicsLink struct API for factory output

### Secondary (MEDIUM confidence)
- [CITED: google.github.io/benchmark/user_guide.html] — RegisterBenchmark usage, custom main() pattern, DoNotOptimize/ClobberMemory patterns
- [CITED: cppreference.com — std::mt19937, std::uniform_real_distribution] — Standard C++ random utilities

### Tertiary (LOW confidence)
- [ASSUMED] `-flto` compatibility with Apple Clang on macOS — standard flag, but not verified in this specific build environment

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Google Benchmark v1.9.5 confirmed on GitHub, FetchContent pattern identical to existing GTest usage
- Architecture: HIGH - All patterns derived from existing codebase (CMakeLists.txt options, examples/CMakeLists.txt structure, ForwardDynamics/InverseDynamics APIs)
- Pitfalls: HIGH - Well-documented Google Benchmark integration issues (TESTING/INSTALL suppression, dead-code elimination, LTO)
- Build system integration: HIGH - Directly reuses established CMake patterns

**Research date:** 2026-06-05
**Valid until:** 2026-07-05 (30 days — Google Benchmark releases are stable but new versions may appear)
