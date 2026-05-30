# Feature Research

**Domain:** C++ Robotics Dynamics — Benchmarks & Examples  
**Researched:** 2026-05-30  
**Confidence:** HIGH

## Feature Landscape

### Table Stakes (Users Expect These)

These are features that any credible robotics dynamics benchmark MUST include. Missing any would make the results unconvincing.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Forward dynamics (ABA) timing** | Primary use case — torque→acceleration. Pinocchio and RBDL both benchmark this as the core metric. | Low | Uses existing `ForwardDynamics::computeAccelerations()` |
| **Inverse dynamics (RNEA) timing** | Second primary use case — acceleration→torque. Both reference libraries include it. | Low | Uses existing `InverseDynamics::computeTorques()` |
| **Average time per call (microseconds)** | Standard reporting unit. Pinocchio benchmarks report 1e-3 × (ns total / NBT) = µs per call. | Low | Simple division |
| **100,000+ iterations per benchmark** | Industry standard (NBT=100k in Pinocchio benchmarks). Enough iterations to drown out noise. | Low | Loop count configurable |
| **Random joint state per iteration** | Avoids warm-state bias. Pinocchio regenerates random `q`, `qdot`, `qddot` per iteration, pausing the timer for randomization. | Low | Use `benchmark::DoNotOptimize` + `state.PauseTiming()` |
| **Release-mode compilation (`-O3 -DNDEBUG`)** | Debug builds are ~10× slower and produce meaningless timing. | Low | CMake `CMAKE_BUILD_TYPE=Release` |
| **Report system/HW info (CPU, cache)** | Google Benchmark auto-reports CPU MHz and cache sizes. Critical for reproducibility. | Low | Built into Google Benchmark library |
| **Comparison vs RBDL and/or Pinocchio** | Users need to know "Is SpatialAlgebra competitive?" — the entire point of benchmarking. | Medium | Requires installing reference libraries OR publishing comparable standalone numbers |
| **Multiple chain sizes (n=2, n=3, n=10+)** | Timing should scale O(n) as predicted. Single-size benchmarks hide algorithmic issues. | Medium | Need programmatic chain generation for arbitrary n |
| **Gravity and no-gravity modes** | Gravity changes the pass structure (affects bias forces). Both modes should be benchmarked. | Low | Existing API supports both via `gravity` parameter |

**Source:** Pinocchio benchmark suite (`pinocchio-benchmarks/src/benchmarks-pinocchio.cpp`, `pinocchio-benchmarks/src/pinocchio-benchmark.cpp`), uses 100k iterations, random states, reports µs per call. RBDL's announcement (2012) reported 0.79×10⁻⁵s for ID, 1.85×10⁻⁵s for FD on a 31-DOF model — same methodology.

### Differentiators (Competitive Advantage)

Features that set SpatialAlgebra's benchmarks apart and demonstrate unique value.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Granular per-pass timing** | Report outward pass vs inward pass timing separately, not just total. Identifies where optimizations matter. | Medium | Requires instrumenting `outwardPass()` and `inwardPass()` |
| **Per-link scaling plot** | Show timing grows linearly with chain length (n=1 through n=20). Proves O(n) complexity. | Medium | Repeated runs with different chain sizes |
| **Round-trip consistency verification in benchmark** | Verify RNEA↔ABA round-trip passes (abs(tau - RNEA(ABA(tau))) < ε) during benchmark. Ensures performance gains don't break correctness. | Low | Reuse existing `TestDynamicsConsistency` logic |
| **CI-based performance regression detection** | Catch slowdowns before they ship. Use Google Benchmark's comparison mode to compare against previous run. | High | Requires storing baseline results, `--benchmark_out` with JSON |
| **OpenMP/no-OpenMP comparison** | Demonstrate that (unlike RBDL's OpenMP variants) the library doesn't need threading for good single-core perf. | Low | Already removed OpenMP — just document |
| **Benchmark with real URDF-derived mass properties** | Use realistic inertia values from example robots (Franka Emika, KUKA LWR) for credibility. | Medium | URDF parsing itself is out of scope, but hand-coded equivalents are feasible |
| **Memory bandwidth measurement** | Show that the library's compact data layout (LowerTriangular packed storage) reduces cache pressure. | Medium | Requires `perf stat -e cache-misses` alongside timing |

### Anti-Features

What to explicitly NOT include in the benchmark suite.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Micro-benchmark of individual math ops** (cross product, spatial add) | These are dominated by Eigen's performance, not SpatialAlgebra's. Micro-benchmarks mask higher-level algorithmic wins. | Benchmark at solver level (ABA, RNEA) — that's what users care about |
| **Single-iteration timing** | Too much noise from CPU frequency scaling, cache state, context switches. Useless numbers. | Use 100k+ iterations, Google Benchmark's adaptive iteration count |
| **Debug-mode benchmarks** | Debug builds can be 10-100× slower and give false impression of performance issues. | Only benchmark Release builds; guard with `#ifndef NDEBUG` |
| **Over-optimized warm-state only** | Always timing the same q/qdot values — cache prediction makes results unrealistic. | Random state per iteration (Pinocchio pattern) |
| **Comparing with different CPU governors** | `powersave` vs `performance` governor changes results by 30-50%. | Document governor used; prefer `performance` |
| **CRBA (Composite Rigid Body Algorithm) benchmarking** | SpatialAlgebra doesn't have CRBA implemented. Attempting to benchmark it would be misleading and out of scope. | Focus on ABA and RNEA — the two algorithms the library implements |
| **Python RNEA (rnea.py) performance** | Python overhead dominates; the standalone RNEA is a reference/validation tool, not a performance target. | Benchmark only the C++ library |

### Benchmark Architecture

```
benchmarks/
  CMakeLists.txt              # Google Benchmark FetchContent
  benchmark_rnea.cpp          # RNEA timing (inverse dynamics)
  benchmark_aba.cpp           # ABA timing (forward dynamics)
  benchmark_rnea_gravity.cpp  # RNEA with gravity
  benchmark_aba_gravity.cpp   # ABA with gravity
  benchmark_scaling.cpp       # O(n) scaling: n=1..20
  benchmark_roundtrip.cpp     # Verify round-trip during benchmark
  benchmark_all.cpp           # Optional: run all benchmarks in one binary
  include/
    benchmark_models.h        # Shared chain definitions for all benchmarks
```

### Real-World Example Architecture

```
examples/
  dynamics.cpp               # Existing: 2-link planar arm
  examples/
    robot_2link_planar.cpp    # NEW: 2-link planar with gravity + trajectory
    robot_3link_spatial.cpp   # NEW: 3-link spatial (RRR) arm demo  
```

**2-link planar arm** — Revolute Z + Revolute Z, planar motion, XY plane. Clear example of coupling effects.
**3-link spatial arm** — Revolute Z + Revolute Y + Revolute Z (RPY-like wrist). Shows 3D rotation effects.

### Benchmark Execution Flow

```
main()
  └─ Report CPU info (Google Benchmark auto-report)
  └─ For each chain size [n=1, 2, 3, 5, 10, 20]:
  │    └─ Build n-link serial chain
  │    └─ BM_RNEA(n):
  │    │    for _ in state:
  │    │      state.PauseTiming()
  │    │      q, qdot, qddot = random_state(n)
  │    │      state.ResumeTiming()
  │    │      computeTorques(q, qdot, qddot)
  │    └─ BM_ABA(n):
  │         for _ in state:
  │           state.PauseTiming()
  │           q, qdot, tau = random_state(n)
  │           state.ResumeTiming()
  │           computeAccelerations(tau)
  └─ BM_RoundTrip(n):
       Verify abs(tau - RNEA(ABA(tau))) < EPSILON
```

## Feature Dependencies

### Dependency Graph

```
CR-02 Bug Fix
  └─ Required by: All benchmarks (benchmarking known-broken code is useless)
  └─ Required by: Robot examples (examples must produce correct dynamics)

Benchmarks (Google Benchmark library)
  └─ Depends on: CR-02 Fix
  └─ Depends on: CMake FetchContent for Google Benchmark
  └─ Depends on: Existing ForwardDynamics, InverseDynamics APIs
  └─ Depends on: benchmark_models.h (shared chain builders)
  └─ Produces: benchmark binaries in build/benchmarks/

Real-World Examples
  └─ Depends on: CR-02 Fix (otherwise examples produce wrong accelerations)
  └─ Depends on: Existing ForwardDynamics, InverseDynamics APIs
  └─ Depends on: Existing example infrastructure (CMakeLists.txt in examples/)
  └─ Adds: examples/robot_2link_planar.cpp, examples/robot_3link_spatial.cpp

Benchmark Models
  └─ Depends on: Existing ForwardDynamicsLink, InverseDynamicsLink structs
  └─ Same structure as tests/TestDynamicsConsistency.cpp (the 3-link chains)
  └─ Generalizes: buildSerialChain(n) → creates n-link chain with configurable params

Eigen 5.x CI Compatibility
  └─ Depends on: Existing CMakeLists.txt find_package(Eigen3) adjustment
  └─ Independent of: Benchmarks and examples (can be done in parallel)
  └─ Verified by: CI matrix addition (ubuntu × g++ × Eigen 5.x)
```

### Key Dependency Rule

> **CR-02 must be fixed BEFORE any benchmark or example work.**  
> Without correct multi-link dynamics, benchmarks measure wrong results, and examples demonstrate incorrect behavior. The `TestDynamicsConsistency.cpp` tests (3-link, branching, gravity) must all pass before any performance or example work is credible.

## MVP Definition

### v1.2 Launch (Minimum Viable)

| Feature | Priority | Rationale |
|---------|----------|-----------|
| CR-02 bug fix in ABA inward pass | **BLOCKING** | All other work depends on correct multi-link dynamics |
| `benchmark_aba.cpp` (ABA timing, n=1..5) | P0 | Core forward dynamics perf — the library's marquee algorithm |
| `benchmark_rnea.cpp` (RNEA timing, n=1..5) | P0 | Core inverse dynamics perf |
| `benchmark_models.h` (shared chain builders) | P0 | Required by both benchmarks above |
| Google Benchmark CMake integration | P0 | Required by all benchmarks |
| Eigen 5.x CI matrix entry | P0 | Maintain compatibility guarantee |
| `robot_2link_planar.cpp` example | P1 | Simple real-world example, reinforces v1.1 gravity work |
| `robot_3link_spatial.cpp` example | P1 | Spatial RRR arm — demonstrates 3D dynamics |
| Gravity benchmarks (RNEA+ABA with gravity) | P1 | Gravity changes computation; users need both numbers |
| Round-trip consistency check in benchmarks | P1 | Catch regressions early |
| Comparison numbers vs RBDL (published) | P2 | Valuable but can be added post-launch |
| Scaling benchmarks (n=1..20 O(n) plot) | P2 | Nice-to-have, not blocking release |
| CI regression tracking | P2 | Infrastructure investment, not feature |

### Deferred to v1.3+

| Feature | Why Deferred |
|---------|-------------|
| URDF model loading | Requires external dependency; hand-coded models sufficient for v1.2 |
| Memory bandwidth profiling | `perf` tooling too platform-specific for CI |
| RBDL integration benchmarks | Installing RBDL in CI requires addon complexity |
| Real-time performance tests (jitter, worst-case) | Beyond scope of basic library benchmarking |
| Warm/cold cache benchmarking | Too subtle for initial release; useful for optimization-focused v1.3 |

### Model Parameters for Benchmarks

All benchmark chains use physically realistic parameters:

```
For each link in chain:
  mass  = 1.0 kg
  COM   = (0, 0.5, 0)  // offset from joint axis
  I_xx = 0.1, I_yy = 0.1, I_zz = 0.1  // diagonal inertia
  joint_type = revolute (axis varies per model)
  transform = translate along X by 1.0m × parent
```

**2-link planar arm example:**
- Link 0: parent=-1, axis=Z, X=identity
- Link 1: parent=0, axis=Z, X=[I, (1,0,0)]

**3-link spatial arm example:**
- Link 0: parent=-1, axis=Z
- Link 1: parent=0, axis=Y
- Link 2: parent=1, axis=Z

### Expected Performance Ranges

Based on RBDL published numbers (31-DOF, i7 920, 2012):
- ID (RNEA): ~7.9 µs for 31 DOF → ~0.25 µs/DOF
- FD (ABA): ~18.5 µs for 31 DOF → ~0.60 µs/DOF

For SpatialAlgebra on modern hardware (Apple M-series, 2026):
- Expected: ~0.1-0.3 µs/DOF for RNEA, ~0.2-0.5 µs/DOF for ABA
- Degradation threshold: >1 µs/DOF for either solver should flag investigation

**Confidence:** MEDIUM — these are inferred from RBDL's published numbers (i7 920, 2012) scaled by ~10× for modern hardware, before actual measurement.

## Sources

- Pinocchio benchmark suite source code: `https://github.com/stack-of-tasks/pinocchio-benchmarks` (C++17, Google Benchmark, 100k iterations, random states) — HIGH confidence
- RBDL original announcement (2012): i7 920, ID=7.9µs, FD=18.5µs for 31-DOF model — MEDIUM confidence (older hardware)
- RBDL documentation: `https://rbdl.github.io/` — HIGH confidence
- Google Benchmark user guide: `https://github.com/google/benchmark/blob/main/docs/user_guide.md` — HIGH confidence
- Featherstone "Rigid Body Dynamics Algorithms" — HIGH confidence
- `test/TestDynamicsConsistency.cpp` — existing 3-link and branching chain test patterns — HIGH confidence
- `examples/dynamics.cpp` — existing 2-link planar example structure — HIGH confidence
