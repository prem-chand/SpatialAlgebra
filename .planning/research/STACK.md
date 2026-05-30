# Stack Research: Performance Benchmarks & Eigen 5.x

**Domain:** C++ Robotics Dynamics — Performance Benchmarks  
**Researched:** 2026-05-30  
**Overall confidence:** HIGH

---

## Executive Summary

Adding performance benchmarking capability and Eigen 5.x CI compatibility to SpatialAlgebra requires three additions to the stack: **Google Benchmark** (v1.9.5) for microbenchmarking dynamics algorithms, a new `benchmarks/` directory with its own `CMakeLists.txt` using FetchContent, and changes to the Eigen version specification in `CMakeLists.txt` to use the range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)`.

RBDL and Pinocchio both use Eigen 3.x as their sole math backend. SpatialAlgebra benchmarks should compare ABA/RNEA performance against RBDL on identical kinematic chain models (joint types, masses, inertias) rather than against Pinocchio, because Pinocchio's template-heavy design and optional analytical derivatives make apples-to-apples comparisons harder. RBDL's `addons/benchmark/` provides the reference pattern: a custom timer class over many samples, parametric model depth, and algorithm-level timing. We should match this pattern using Google Benchmark's `->Range()` parameterization.

**Key change to existing CMakeLists.txt:** Remove the bare `find_package(Eigen3 REQUIRED NO_MODULE)` and replace with the version-range-aware form introduced in Eigen 3.4.1.

---

## Recommended Stack

### Core Benchmarking Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Google Benchmark | v1.9.5 | C++ microbenchmark library | Industry standard; CMake FetchContent integration; parameterized benchmarks via `->Range()`; statistical rigor (warmup, iteration tuning, outlier rejection). Released Jan 2026. |

### Alternative (Not Recommended for This Use Case)

| Technology | Version | Purpose | Why Not |
|------------|---------|---------|----------|
| ankerl::nanobench | v4.3.11 | Single-header microbenchmark | Excellent for quick ad-hoc perf measurements, but lacks the CI-friendly structured output and parameterized multi-DOF range testing needed for comparing against RBDL. The token-thin integration is outweighed by Google Benchmark's ecosystem. |
| Custom timing loops (RBDL's approach) | N/A | Hand-written timer | RBDL's `addons/benchmark/Timer.h` uses `std::chrono::high_resolution_clock`. This works but provides no statistical rigor, no warmup, no outlier detection. Google Benchmark's `State` abstraction handles all of this. |

### Comparison Libraries (External Dependencies)

| Library | Version | Purpose | How We Use It |
|---------|---------|---------|---------------|
| RBDL | v3.3.1 | Reference dynamics library | Build identical kinematic chain models, time identical ABA/RNEA calls. RBDL's benchmark model generator (`addons/benchmark/model_generator.cc`) creates planar trees with configurable DOF — we replicate those models in SpatialAlgebra. |
| Pinocchio | 3.x (devel) | Advanced dynamics library | Compare against only if a credible common model (e.g. URDF) can be loaded. Pinocchio's template metaprogramming gives it a structural performance advantage that is not directly comparable to SpatialAlgebra's hand-written classes. **Recommend comparing against RBDL only for v1.2.** |

### What NOT to Add

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| Pinocchio as direct benchmark target in v1.2 | Template-heavy design makes comparison unfair and hard to interpret; Pinocchio also uses code-gen (CppADCodeGen) for joint models, which SpatialAlgebra doesn't do. | RBDL only — same Featherstone algorithms, same Eigen backend, comparable design philosophy. |
| Google Test for benchmarks | GTest has `benchmark::State`-like features but they're not meant for performance measurement. Mixing test and benchmark concerns is a well-known anti-pattern. | Google Benchmark in a separate `benchmarks/` directory. |
| OpenMP for benchmark parallelization | The OpenMP-removal decision from LowerTriangular (PROJECT.md) means we must not introduce threading in benchmarks. | Single-threaded benchmarks only; report wall-clock time per call. |
| Python benchmarking | The Python RNEA is a standalone educational reference, not part of the performance-critical path. | C++ benchmarks only. |
| CI coverage on benchmark runs | Benchmarks are stochastic and machine-dependent; failing CI on a 5% regression is noisy. | Manual `cmake --build build/benchmarks && ./build/benchmarks/...` workflow; tracked in a RELEASES.md or BENCHMARKS.md. |

---

## Installation

### Homebrew (macOS Development)

```bash
# Already installed
brew install eigen googletest

# New for benchmarks
brew install google-benchmark
```

### CI Dependencies (GitHub Actions)

```yaml
# Add to existing ci.yml install steps
- name: Install benchmark dependencies (macOS)
  if: runner.os == 'macOS'
  run: brew install google-benchmark

- name: Install benchmark dependencies (Ubuntu)
  if: runner.os == 'Linux'
  run: sudo apt-get install -y libbenchmark-dev
```

### CMake FetchContent (Fallback)

For users who don't want a system install, add this to `benchmarks/CMakeLists.txt`:

```cmake
include(FetchContent)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(googlebenchmark)
```

---

## Alternatives Considered

| Aspect | Google Benchmark | nanobench | Custom Timer |
|--------|-----------------|-----------|--------------|
| **Integration** | FetchContent or `find_package` | Single-header `#include` | Zero deps |
| **Statistical rigor** | Iteration tuning, warmup, min time, outlier rejection | Error % reporting, minimum samples | None |
| **Parameterized ranges** | `->Range(8, 64)` for DOF sweep | Manual loop | Manual loop |
| **Output formats** | Console, CSV, JSON, bencher-compatible | Console table, CSV | Custom printf |
| **CI readiness** | `--benchmark_format=json`, exit code on regression | CSV output | None |
| **Build time impact** | Compiles ~30s first time | ~2.4s (header) | None |
| **Ecosystem** | Bencher.dev, historical tracking, community extensions | Minimal | None |
| **RBDL pattern match** | Does not match RBDL's approach exactly | Closer (single file timer) | Exact match |

**Recommendation:** Use Google Benchmark. The build time is paid once per rebuild of the benchmark target (not the main library). The statistical rigor and output formats are essential for credible vs-RBDL comparisons. RBDL's own benchmark approach (custom timer, many samples) is naive — we should not replicate it, we should supersede it.

---

## Eigen 5.x Compatibility

### The Problem

Eigen 5.0.0 (September 2025) and 5.0.1 (November 2025) changed their CMake version-compatibility range. The current `CMakeLists.txt` has:

```cmake
find_package(Eigen3 REQUIRED NO_MODULE)
```

This works with both Eigen 3.x and 5.x, but:
1. It does not specify a minimum version
2. Eigen 5.x's `Eigen3Config.cmake` reports its version as 5.0.x, so `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` may **fail** with Eigen 5.x if the package's `version.cmake` uses stricter comparison

### The Solution

Eigen 3.4.1+ introduced a version range syntax:

```cmake
find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)
```

This means: "Accept any version >= 3.4.1 but < 6.0.0." It will match both Eigen 3.4.x and Eigen 5.0.x.

**Change required in `CMakeLists.txt` line 12:**
- Current: `find_package(Eigen3 REQUIRED NO_MODULE)`
- New: `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)`

### Breaking Changes in Eigen 5.x to Watch For

| Change | Impact on SpatialAlgebra | Mitigation |
|--------|--------------------------|------------|
| CMake build system modernized — older properties removed | Low, because we use `NO_MODULE` (config-mode) | Already works; the range syntax handles this |
| Requires C++14 (not C++03) | None — we target C++17 | No action needed |
| All LGPL licensed code removed (Constrained Conjugate Gradient) | None — we don't use unsupported modules | No action needed |
| `EIGEN_HAS_CXX11` macros removed | None — we use `__cplusplus >= 201703L` | Verify with a CI build on Eigen 5.x |

### CI Matrix Addition

Add an `eigen-version` dimension to the CI matrix (or a separate job):

```yaml
matrix:
  os: [ubuntu-latest, macos-latest]
  compiler: [g++, clang++]
  eigen-version: [system, 5.0]  # 'system' = apt/brew default (3.4), '5.0' = explicit 5.x
```

For the Eigen 5.x jobs:
```yaml
- name: Install Eigen 5.x (macOS)
  if: matrix.eigen-version == '5.0' && runner.os == 'macOS'
  run: |
    brew tap libeigen/eigen
    brew install eigen@5.0

- name: Install Eigen 5.x (Ubuntu)
  if: matrix.eigen-version == '5.0' && runner.os == 'Linux'
  run: |
    git clone --depth 1 --branch 5.0.1 https://gitlab.com/libeigen/eigen.git /tmp/eigen5
    cmake -B /tmp/eigen5/build -S /tmp/eigen5
    sudo cmake --install /tmp/eigen5/build
```

---

## Benchmark Directory Structure

```
benchmarks/
├── CMakeLists.txt          # FetchContent + benchmark targets
├── BenchmarkUtils.h        # Model creation helpers (shared by benchmarks)
├── BenchmarkModels.h       # Canonical test chains (2-link, 3-link, 6-link spatial)
├── RBDLCompat.h            # Thin wrappers to build equivalent RBDL models
├── ABA_Benchmark.cpp       # Forward dynamics benchmarks
├── RNEA_Benchmark.cpp      # Inverse dynamics benchmarks
├── Plucker_Benchmark.cpp   # Transform microbenchmarks
└── SpatialVector_Benchmark.cpp  # Vector operation microbenchmarks
```

### `benchmarks/CMakeLists.txt` Pattern

```cmake
set(CMAKE_CXX_STANDARD 17)

# Find Eigen3 with version range (3.4.x through 5.x)
find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)

# Google Benchmark via FetchContent (fallback)
include(FetchContent)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googlebenchmark)

# RBDL: try system find, skip if not found (informational only)
find_package(RBDL QUIET)

# Main benchmark executable
add_executable(BenchmarkABA ABA_Benchmark.cpp)
target_link_libraries(BenchmarkABA
    SpatialAlgebra
    benchmark::benchmark
    benchmark::benchmark_main
)

add_executable(BenchmarkRNEA RNEA_Benchmark.cpp)
target_link_libraries(BenchmarkRNEA
    SpatialAlgebra
    benchmark::benchmark
    benchmark::benchmark_main
)
```

---

## Integration Points with Existing Stack

| Existing Component | Integration | Changes Required |
|-------------------|-------------|------------------|
| **CMakeLists.txt** root | Add `add_subdirectory(benchmarks)` behind an option like `BUILD_BENCHMARKS` | +3 lines |
| **CMakeLists.txt** Eigen find_package | Change `REQUIRED NO_MODULE` → `3.4...5 REQUIRED NO_MODULE` | Line 12 change |
| **GitHub Actions CI** | No automatic benchmark runs; separate manual workflow | None |
| **Include headers** | `BenchmarkUtils.h` reuses `SpatialAlgebra.h` umbrella header | None |
| **GTest tests** | Separate directory; Google Benchmark is in `benchmarks/` | No GTest dependency |
| **RBDL** | Optional system dependency; `find_package(RBDL QUIET)` — benchmarks that compare against RBDL are guarded by `#ifdef RBDL_FOUND` | New CMake module or QUIET find |

---

## Design Decisions

### Why Google Benchmark over nanobench or custom timers

1. **Parameterized ranges** — `->Range(2, 64)` on DOF count is the primary benchmark axis. Google Benchmark handles this natively; nanobench requires manual loops.
2. **Statistical rigor** — automatic iteration tuning and warmup eliminate cold-cache bias, which matters for small dynamics computations (~1-20 µs per call).
3. **CI-friendly output** — `--benchmark_format=json` produces parseable output for comparison across runs. nanobench's CSV output works but is less structured.
4. **Ecosystem compatibility** — Bencher.dev and historical tracking tools integrate with Google Benchmark's output format natively.

### Why RBDL-only comparison (not Pinocchio)

Pinocchio uses Eigen extensively but employs C++ template metaprogramming, expression template optimizations, and CppADCodeGen for joint models. Direct performance comparison with SpatialAlgebra (hand-written classes, no code generation) would either be misleading or require disabling Pinocchio's optimizations. RBDL shares the same Featherstone algorithms, same Eigen backend, and similar design philosophy (hand-written C++ classes, no heavy metaprogramming).

### Why benchmarks are not in CI

Benchmark results depend on CPU model, clock speed, thermal throttling, and system load. A CI job running on shared GitHub Actions runners produces non-reproducible results. Instead, benchmarks should be:
- Documented in a `BENCHMARKS.md` with the machine spec used
- Run manually before releases
- Timed on a reference machine with CPU frequency pinned

---

## Sources

| Source | URL | Confidence |
|--------|-----|------------|
| Google Benchmark v1.9.5 release | https://github.com/google/benchmark/releases | HIGH |
| Google Benchmark CMake integration guide | https://chromium.googlesource.com/external/github.com/google/benchmark/ | HIGH |
| Homebrew google-benchmark formula | https://formulae.brew.sh/formula/google-benchmark | HIGH |
| nanobench v4.3.11 docs | https://nanobench.ankerl.com/tutorial.html | HIGH |
| Eigen 5.0 release notes | https://libeigen.gitlab.io/releases/5.0/ | HIGH |
| Eigen CMake version range syntax | https://libeigen.gitlab.io/eigen/docs-5.0/TopicCMakeGuide.html | HIGH |
| Eigen 5.0.1 tags | https://gitlab.com/libeigen/eigen/-/tags/5.0.1 | HIGH |
| RBDL benchmarking architecture | https://deepwiki.com/rbdl/rbdl/7.2-benchmarking | HIGH |
| RBDL v3.3.1 release | https://github.com/rbdl/rbdl/releases | HIGH |
| Pinocchio Eigen dependency | https://github.com/stack-of-tasks/pinocchio/blob/devel/CMakeLists.txt | HIGH |
| Pinocchio vs RBDL performance paper | https://hal-laas.archives-ouvertes.fr/hal-01866228 | MEDIUM |
