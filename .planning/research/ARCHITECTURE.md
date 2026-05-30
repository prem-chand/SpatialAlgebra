# Architecture Research

**Domain:** C++ Robotics Dynamics — Benchmarks & Examples  
**Researched:** 2026-05-30  
**Confidence:** HIGH

---

## Standard Architecture

### Benchmark Executable Structure

The industry-standard pattern for C++ performance benchmarks is a dedicated `benchmarks/` directory with its own `CMakeLists.txt`, kept entirely separate from tests. This follows the same separation-of-concerns logic as tests — benchmarks are not tests, they are measurement harnesses.

**Google Benchmark** is the de facto standard (developed by Google, widely adopted in C++ ecosystem, header-only style with compiled library). It handles:
- Statistical sampling (warm-up, iteration count selection)
- Output formatting (console tables, JSON export)
- Parameterized benchmarks (`Arg()`, `Range()`, `Ranges()` for multi-dimensional sweeps)
- Fixture support via `benchmark::Fixture` for setup/teardown

Reference: [Google Benchmark GitHub](https://github.com/google/benchmark/) (HIGH confidence), [Integrating with CMake](https://www.studyplan.dev/google-benchmark/integrating-google-benchmark) (MEDIUM confidence — tutorial source, verified pattern matches multiple real projects).

### Comparison Benchmark Architecture

Comparison benchmarks (SpatialAlgebra vs RBDL vs Pinocchio) follow a **parameterized harness** pattern:

1. **Model factory** — Construct equivalent kinematic chains in each library (same DOF count, link masses, transforms)
2. **Input generator** — Produce random joint states (q, qdot, qddot) and torques
3. **Measurement loop** — Call equivalent algorithm (ABA/RNEA) in each library, measuring wall time
4. **DOF sweep** — Repeat for n = 2, 4, 6, 8, 12, 18 links (parameterized via `Range()`)

This matches the methodology from the IROS 2019 benchmarking paper (Neuman et al., "Benchmarking and Workload Analysis of Robot Dynamics Algorithms"). RBDL and Pinocchio are both body-coordinate libraries using Eigen, making them directly comparable to SpatialAlgebra.

### Robot Example Structure

Robot examples follow the existing pattern in `examples/`:
- Single-file `.cpp` with a `main()` that constructs a kinematic chain
- Uses `ForwardDynamics` or `InverseDynamics` classes
- Outputs computed results to stdout

New robot examples should demonstrate real-world applicability:
- **2-link planar arm**: Z-axis revolute joints in XY plane (evolution of existing `dynamics.cpp`)
- **3-link spatial arm**: Three revolute joints about different axes (Z, Y, Y) — creates a spatial reach

### CMake Integration for External Libraries

Two distinct patterns exist:

**Pattern A: With-installed-library (RBDL/Pinocchio)**
- Optional `find_package()` guarded by CMake options (e.g., `-DSA_BUILD_BENCHMARKS_COMPARISON=ON`)
- RBDL uses `FindRBDL.cmake` module (no Config mode); install via `brew install rbdl` or custom path
- Pinocchio uses `find_package(pinocchio REQUIRED)` + `target_link_libraries(... pinocchio::pinocchio)` (modern CMake)

**Pattern B: FetchContent (Google Benchmark)**
- Download and build from source at configure time
- No system install required
- Well-suited for CI environments

---

## Recommended Project Structure

```
SpatialAlgebra/
├── CMakeLists.txt                  # Root — add_subdirectory for benchmarks, examples
├── include/                        # (unchanged)
├── src/                            # (unchanged)
├── tests/                          # (unchanged)
├── examples/                       # (existing, add new examples)
│   ├── CMakeLists.txt              # (modified)
│   ├── basic_vectors.cpp
│   ├── transforms.cpp
│   ├── inertia.cpp
│   ├── dynamics.cpp                # Existing 2-link example
│   ├── robot_2link_planar.cpp      # NEW: Realistic 2-link planar arm
│   └── robot_3link_spatial.cpp     # NEW: 3-link spatial arm
├── benchmarks/                     # NEW
│   ├── CMakeLists.txt              # Google Benchmark fetch + executables
│   ├── CMakePresets.json           # Release-mode preset for benchmarks
│   ├── BenchmarkUtils.h            # Shared helpers (model factories, random input gen)
│   ├── bench_aba.cpp               # ABA forward dynamics benchmarks
│   ├── bench_rnea.cpp              # RNEA inverse dynamics benchmarks
│   ├── bench_plucker.cpp           # Plücker transform micro-benchmarks
│   ├── bench_cross_product.cpp     # Cross product micro-benchmarks
│   └── bench_comparison.cpp        # SpatialAlgebra vs RBDL vs Pinocchio
└── cmake/                          # NEW (optional dependency modules)
    ├── FindRBDL.cmake              # Copy from RBDL examples
    └── FindPinocchio.cmake         # (only if system doesn't provide Config mode)
```

### Structure Rationale

**Why separate `benchmarks/` from `tests/`:**
- Different dependency: benchmarks need `benchmark::benchmark`; tests need `GTest::GTest`
- Different build profile: benchmarks must be built in Release/RelWithDebInfo for meaningful results; tests can run in Debug
- Different invocation: `ctest` runs tests; `./benchmarks/bench_aba --benchmark_format=json` runs benchmarks
- Prevents benchmark dependencies leaking into test infrastructure

**Why optional dependency for RBDL/Pinocchio:**
- RBDL and Pinocchio are heavy dependencies (~5MB+ each, with Boost/URDF/other transitive deps)
- Not all users need comparison benchmarks; they are CI-only or developer-only
- Guarded by `SA_BUILD_COMPARISON_BENCHMARKS` (default OFF)

**Why `cmake/` modules directory:**
- RBDL does NOT provide CMake Config mode — it requires a custom `FindRBDL.cmake` module
- Keeping modules in `cmake/` is standard CMake practice and keeps root `CMakeLists.txt` clean

---

## Architectural Patterns

### Pattern 1: Parameterized Benchmark (ABA Scale Test)

```cpp
// benchmarks/bench_aba.cpp
#include <benchmark/benchmark.h>
#include "SpatialAlgebra.h"
using namespace SpatialAlgebra;

// Factory: creates an n-link serial chain
std::vector<Link> CreateSerialChain(int n_links) {
    std::vector<Link> links;
    for (int i = 0; i < n_links; ++i) {
        Link link;
        link.parent = i - 1;
        // ... set inertia, transform, joint axis
        links.push_back(link);
    }
    return links;
}

static void BM_ABA_ForwardDynamics(benchmark::State& state) {
    const int n = state.range(0);
    ForwardDynamics fd;
    fd.links = CreateSerialChain(n);
    Eigen::VectorXd tau = Eigen::VectorXd::Random(n);

    for (auto _ : state) {
        fd.computeAccelerations(tau);
        benchmark::DoNotOptimize(fd.links[0].qddot);
    }
}
BENCHMARK(BM_ABA_ForwardDynamics)->RangeMultiplier(2)->Range(2, 18);
```

Uses `Range()` to sweep DOF count (2, 4, 8, 18) — Google Benchmark handles multiplication automatically.

### Pattern 2: Comparison Benchmark (SpatialAlgebra vs RBDL)

```cpp
// benchmarks/bench_comparison.cpp (guarded by SA_BUILD_COMPARISON_BENCHMARKS)
#include <benchmark/benchmark.h>

// Forward declare comparison harness (implemented in separate TU per library)
void BM_ABA_SpatialAlgebra(benchmark::State& state, int n_links);
void BM_ABA_RBDL(benchmark::State& state, int n_links);

static void BM_ABA_Comparison(benchmark::State& state) {
    const int n = state.range(0);
    // Dispatch to both implementations in the same benchmark
    // (or use separate benchmarks and compare via JSON output)
}
```

**Better approach:** Separate benchmarks per library, compare via `--benchmark_format=json` output. Each benchmark TU links only against its own library — no cross-linkage needed.

### Pattern 3: Micro-benchmark (Plücker Transform)

```cpp
static void BM_PluckerTransformMotion(benchmark::State& state) {
    PluckerTransform X(Rotation::Random(), Vector3d::Random());
    MotionVector mv(Vector3d::Random(), Vector3d::Random());

    for (auto _ : state) {
        MotionVector result = X.transformMotion(mv);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(BM_PluckerTransformMotion);
```

### Pattern 4: Robot Example Structure

```cpp
// examples/robot_3link_spatial.cpp
#include "ForwardDynamics.h"
#include <iostream>

using namespace SpatialAlgebra;

int main() {
    ForwardDynamics robot;
    
    // Link 1: base — Z-axis rotation
    Link l1;
    l1.parent = -1;
    l1.X = PluckerTransform(Rotation::Identity(), Vector3d::Zero());
    l1.I = RigidBodyInertia(1.0, Vector3d::Zero(), /* ... */);
    l1.S = MotionVector(Vector3d::UnitZ(), Vector3d::Zero());
    robot.links.push_back(l1);
    
    // Link 2: Y-axis rotation, offset along Z
    // Link 3: Y-axis rotation, offset along X
    
    Eigen::VectorXd tau(3);
    tau << 1.0, 0.5, 0.3;
    robot.computeAccelerations(tau, Vector3d(0, 0, -9.81));
    
    for (size_t i = 0; i < robot.links.size(); ++i)
        std::cout << "Joint " << i << ": " << robot.links[i].qddot << " rad/s^2\n";
    
    return 0;
}
```

---

## Data Flow

### Benchmark Data Flow

```
Input Generator                    Benchmark Loop                     Output
(random q, qdot, qddot, tau)  →   ForwardDynamics::computeAccelerations()  →  qddot (double)
(parameter: n_links)              InverseDynamics::computeTorques()         →  tau (VectorXd)
                                  PluckerTransform::transformMotion()       →  SpatialVector
                                                                                + wall time
```

Key invariants:
- Random seeds are fixed (`std::mt19937` with constant seed) for reproducibility
- Input generation happens OUTSIDE the measured loop (setup excluded from timing)
- `benchmark::DoNotOptimize()` prevents dead-code elimination on outputs
- Physical correctness NOT checked in benchmarks (that's what tests are for)

### CI Integration Flow

```
PR push → CI workflow
  ├── Eigen 3.4 + g++ → build + test + coverage
  ├── Eigen 3.4 + clang++ → build + test
  ├── Eigen 5.x + g++ → build + test (NEW matrix entry)
  └── Eigen 5.x + clang++ → build + test (NEW matrix entry)

Benchmarks are NOT run in CI (variance from shared CPU). Instead:
  - Benchmarks are compiled in CI (build check only)
  - Full benchmark runs are manual or scheduled on dedicated hardware
```

---

## Integration Points

### CMakeLists.txt Changes

**Root CMakeLists.txt** — four changes:

1. **Eigen 5.x version detection** — Replace the hard-coded `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` with a version-range aware pattern that works with both Eigen 3.3+ and 5.x:

```cmake
# Eigen version detection: supports Eigen 3.3+ through 5.x
# Eigen 5.0+ changed CMake version compatibility, requiring range syntax
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.19")
  find_package(Eigen3 3.3...5 REQUIRED NO_MODULE)
else()
  find_package(Eigen3 REQUIRED NO_MODULE)
  # Verify minimum version manually (CMake range syntax unavailable)
  if(Eigen3_VERSION VERSION_LESS 3.3)
    message(FATAL_ERROR "Eigen3 version >= 3.3 required (found ${Eigen3_VERSION})")
  endif()
endif()
```

This follows the pattern adopted by PCL (PR #6354), VowpalWabbit (PR #4728), and others in 2025-2026. The range syntax `3.3...5` means ">= 3.3 and < 6.0.0" which covers both Eigen 3.4.x and Eigen 5.x.

2. **Add `benchmarks/` subdirectory** (guarded by option):

```cmake
option(SA_BUILD_BENCHMARKS "Build performance benchmarks" OFF)
option(SA_BUILD_COMPARISON_BENCHMARKS "Build comparison benchmarks (requires RBDL/Pinocchio)" OFF)

# ... existing library targets ...

if(SA_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
```

3. **Add `cmake/` to module path** (for FindRBDL.cmake):

```cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
```

4. **Raise minimum CMake version to 3.19** (to use range syntax natively) OR keep the version-guard pattern above. Decision: keep `cmake_minimum_required(VERSION 3.10)` and use the version-guard pattern, because bumping to 3.19 excludes Ubuntu 18.04 users.

**benchmarks/CMakeLists.txt** — new file:

```cmake
cmake_minimum_required(VERSION 3.10)

# Google Benchmark — try system install first, fallback to FetchContent
find_package(benchmark QUIET)
if(NOT benchmark_FOUND)
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.1
  )
  FetchContent_MakeAvailable(googlebenchmark)
endif()

# Common benchmark utilities (header-only)
# (no target needed — header-only utility)

# ABA benchmark
add_executable(bench_aba bench_aba.cpp)
target_link_libraries(bench_aba PRIVATE SpatialAlgebra benchmark::benchmark benchmark::benchmark_main)

# RNEA benchmark
add_executable(bench_rnea bench_rnea.cpp)
target_link_libraries(bench_rnea PRIVATE SpatialAlgebra benchmark::benchmark benchmark::benchmark_main)

# Micro-benchmarks
add_executable(bench_plucker bench_plucker.cpp)
target_link_libraries(bench_plucker PRIVATE SpatialAlgebra benchmark::benchmark benchmark::benchmark_main)

add_executable(bench_cross_product bench_cross_product.cpp)
target_link_libraries(bench_cross_product PRIVATE SpatialAlgebra benchmark::benchmark benchmark::benchmark_main)

# Comparison benchmarks (optional, requires RBDL and/or Pinocchio)
if(SA_BUILD_COMPARISON_BENCHMARKS)
  find_package(RBDL QUIET)
  find_package(pinocchio QUIET)

  if(RBDL_FOUND)
    add_executable(bench_comparison_rbdl bench_comparison_rbdl.cpp)
    target_link_libraries(bench_comparison_rbdl PRIVATE SpatialAlgebra RBDL::rbdl benchmark::benchmark benchmark::benchmark_main)
  endif()

  if(pinocchio_FOUND)
    add_executable(bench_comparison_pinocchio bench_comparison_pinocchio.cpp)
    target_link_libraries(bench_comparison_pinocchio PRIVATE SpatialAlgebra pinocchio::pinocchio benchmark::benchmark benchmark::benchmark_main)
  endif()
endif()
```

**examples/CMakeLists.txt** — add two new example targets:

```cmake
# New robot examples
add_executable(example_robot_2link robot_2link_planar.cpp)
target_link_libraries(example_robot_2link SpatialAlgebra Eigen3::Eigen)

add_executable(example_robot_3link robot_3link_spatial.cpp)
target_link_libraries(example_robot_3link SpatialAlgebra Eigen3::Eigen)
```

**CI workflow (`.github/workflows/ci.yml`)** — add Eigen 5.x matrix entries:

```yaml
strategy:
  matrix:
    eigen: [3.4, 5.0]           # NEW: test both Eigen versions
    os: [ubuntu-latest, macos-latest]
    compiler: [g++, clang++]
    exclude:                     # Optional: reduce matrix size
      - eigen: 5.0
        os: macos-latest
        compiler: clang++
```

Add `DEPS="eigen3"` for Eigen 3.4 and appropriate package for Eigen 5.x. On Ubuntu, Eigen 5.x is available via `libeigen3-dev` (5.0+) or manual install. On macOS, `brew install eigen` picks up Homebrew's version (currently 5.0.1).

### File Summary: Modified vs New

| Action | File | Reason |
|--------|------|--------|
| **MODIFY** | `CMakeLists.txt` | Eigen version detection, benchmarks subdirectory, cmake module path |
| **MODIFY** | `.github/workflows/ci.yml` | Add Eigen 5.x build matrix entries |
| **MODIFY** | `examples/CMakeLists.txt` | Add new robot example targets |
| **NEW** | `benchmarks/CMakeLists.txt` | Google Benchmark integration |
| **NEW** | `benchmarks/BenchmarkUtils.h` | Shared benchmark helpers (model factories) |
| **NEW** | `benchmarks/bench_aba.cpp` | ABA forward dynamics benchmarks |
| **NEW** | `benchmarks/bench_rnea.cpp` | RNEA inverse dynamics benchmarks |
| **NEW** | `benchmarks/bench_plucker.cpp` | Plücker transform micro-benchmarks |
| **NEW** | `benchmarks/bench_cross_product.cpp` | Cross product micro-benchmarks |
| **NEW** | `benchmarks/bench_comparison.cpp` | SpatialAlgebra vs RBDL vs Pinocchio (guarded) |
| **NEW** | `examples/robot_2link_planar.cpp` | Realistic 2-link planar arm example |
| **NEW** | `examples/robot_3link_spatial.cpp` | 3-link spatial arm example |
| **NEW** | `cmake/FindRBDL.cmake` | Module for finding RBDL (required because RBDL has no Config mode) |
| **UNCHANGED** | `include/*.h`, `src/*.cpp`, `tests/*.cpp` | Core library untouched |

### RBDL/Pinocchio Integration Details

**RBDL CMake:**
- RBDL does NOT install CMake config files (no `RBDLConfig.cmake`)
- Must use custom `FindRBDL.cmake` module (available from [RBDL examples](https://github.com/rbdl/rbdl/tree/master/examples/simple))
- Alternatively, use a thin find module that checks `rbdl/rbdl.h` and `-lrbdl`
- On macOS: `brew install rbdl` installs headers to `/usr/local/include/rbdl/` and library to `/usr/local/lib/librbdl.dylib`
- CMake variables set by FindRBDL.cmake: `RBDL_INCLUDE_DIR`, `RBDL_LIBRARY`
- Important: RBDL depends on Eigen (same as SpatialAlgebra), so Eigen must be found first

**Pinocchio CMake:**
- Modern CMake: `find_package(pinocchio REQUIRED)` + `target_link_libraries(... pinocchio::pinocchio)`
- Pinocchio 3+ splits into sublibraries; `pinocchio::pinocchio` is a meta-target linking all of them
- Available via `brew install pinocchio` or `conda install pinocchio`
- Transitive dependencies: Boost (filesystem, serialization), urdfdom (if URDF support enabled)
- For comparison benchmarks, only core dynamics algorithms needed (no URDF parsing)
- Pinocchio requires CMake 3.10+ (increased from the Eigen 3.0.5 dependency)
- Pinocchio API: `pinocchio::Model`, `pinocchio::Data`, `pinocchio::aba()`, `pinocchio::rnea()`

---

## Scalability Considerations

| Concern | At n=2 links | At n=10 links | At n=100 links |
|---------|--------------|---------------|----------------|
| **ABA time** | O(n) ~1μs | O(n) ~5μs | O(n) ~50μs |
| **Test granularity** | Single benchmark | Range sweep | Range sweep or discrete |
| **RBDL comparison** | 1:1 match feasible | 1:1 match feasible | Setup complexity grows linearly |
| **Pinocchio comparison** | 1:1 match feasible | 1:1 match feasible | Model creation scripted |
| **CI execution** | Not in CI (variance) | Not in CI (variance) | Not in CI (variance) |
| **Memory** | <1KB links | ~5KB links | ~50KB links |

Benchmarks are NOT run in CI due to CPU frequency scaling, hyperthread, and shared-host variance. CI compiles benchmarks (build check only). Dedicated benchmark runs use `--benchmark_format=json` for archival.

---

## Key Decisions

| Decision | Rationale |
|----------|-----------|
| Separate `benchmarks/` directory | Different dep, different build profile, different invocation |
| FetchContent for Google Benchmark (with system fallback) | No system install required; works on all CI platforms |
| `SA_BUILD_COMPARISON_BENCHMARKS` default OFF | RBDL/Pinocchio are heavy deps; not needed for library users |
| Eigen 5.x range syntax with CMake version guard | Supports both CMake <3.19 and >=3.19; covers Eigen 3.3 through 5.x |
| Micro-benchmarks + macro-benchmarks | Plücker/cross-product measure primitive ops; ABA/RNEA measure algorithm-level perf |
| Comparison via separate executables (not in-process) | Each TU links only its library; no symbol conflicts; easier to debug |
| Keep `cmake_minimum_required(VERSION 3.10)` | Don't force CMake upgrade on downstream users; use version guard for range syntax |

---

## Sources

- **Google Benchmark CMake integration**: [studyplan.dev](https://www.studyplan.dev/google-benchmark/integrating-google-benchmark), [GitHub](https://github.com/google/benchmark/) — HIGH confidence
- **Eigen 5.x CMake changes**: [Eigen 5.0 Changelog](https://gitlab.com/libeigen/eigen/-/commit/549bf8c75b6aae071cde2f28aa48f16ee3ae60b0), [PCL PR #6354](https://github.com/PointCloudLibrary/pcl/pull/6354), [VowpalWabbit PR #4728](https://github.com/VowpalWabbit/vowpal_wabbit/pull/4728) — HIGH confidence
- **RBDL CMake**: [RBDL examples](https://github.com/rbdl/rbdl/tree/master/examples/simple), [RBDL GitHub](https://github.com/rbdl/rbdl) — HIGH confidence
- **Pinocchio CMake**: [Pinocchio CMakeLists.txt](https://github.com/stack-of-tasks/pinocchio/blob/master/CMakeLists.txt), [Pinocchio minimal example](https://github.com/stack-of-tasks/pinocchio-minimal) — HIGH confidence
- **Benchmarking methodology for robotics dynamics**: Neuman et al., IROS 2019 — MEDIUM confidence (academic paper, not project-specific)
- **Eigen 5.x version detection pattern**: PCL issue #6351, [Eigen CMake Guide](https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html) — HIGH confidence
