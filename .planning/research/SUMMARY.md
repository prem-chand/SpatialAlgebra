# Project Research Summary

**Project:** SpatialAlgebra — v1.2 (Benchmarks & Examples)
**Domain:** C++ Robotics Dynamics — Performance Benchmarking, Stability Fixes & CI Compatibility
**Researched:** 2026-05-30
**Confidence:** HIGH

## Executive Summary

SpatialAlgebra v1.2 adds performance benchmarking capability, Eigen 5.x CI compatibility, and real-world robot examples to a C++17 library implementing Featherstone's spatial vector algebra for rigid body dynamics. The research confirms three parallel workstreams: (1) fixing the CR-02 multi-link ABA inward-pass bug which corrupts all forward dynamics for chains with 3+ joints, (2) adding Eigen 5.x CMake version-range compatibility to prevent build breakage on macOS where Homebrew ships Eigen 5.x, and (3) building a `benchmarks/` directory with Google Benchmark v1.9.5 for ABA/RNEA microbenchmarks with parameterized DOF sweep.

The critical finding is a strict ordering dependency: **CR-02 must be fixed before any benchmark or example work**. Benchmarking a buggy solver produces meaningless numbers, and robot examples that demonstrate incorrect physics damage library credibility. The Eigen 5.x CI work is the only stream that can proceed in parallel with the bug fix.

Three key risks emerge: (1) the CR-02 root cause is a subtle initialization-order bug where bias forces (pa) are computed using rigid-body inertia I_i instead of the final articulated inertia I_A that includes accumulated children per Featherstone Algorithm 7.3, (2) Eigen 5.x's CMake version pin breakage will silently fail on `brew install eigen` updates unless addressed with the range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)`, and (3) cross-library comparison benchmarks against RBDL v3.3.1 carry a frame-convention mismatch risk that requires exact numerical equivalence verification before any performance comparison is credible. The recommended approach is to fix CR-02 first, add Eigen 5.x CI in parallel, then build benchmarks, examples, and finally optional RBDL comparison.

## Key Findings

### Recommended Stack

Adding performance benchmarking and Eigen 5.x compatibility requires minimal new dependencies — three additions to the existing stack, zero changes to the core library.

**Core technologies:**
- **C++17**: Existing language standard — unchanged, no migration needed
- **Eigen 3.4...5**: Linear algebra backend — switch from version-pinned `find_package(Eigen3 REQUIRED NO_MODULE)` to range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` to support both Eigen 3.4.x and 5.x
- **Google Benchmark v1.9.5**: C++ microbenchmark library — FetchContent integration with parameterized `->Range()` for DOF sweep, statistical rigor (warmup, iteration tuning, outlier rejection), JSON output for archival
- **RBDL v3.3.1**: Reference dynamics library for comparison benchmarks — optional dependency guarded by `SA_BUILD_COMPARISON_BENCHMARKS` (default OFF); RBDL shares Featherstone algorithms and Eigen backend for the fairest comparison
- **Pinocchio**: Explicitly NOT recommended as a v1.2 benchmark target — template-heavy design makes apples-to-apples comparisons misleading; deferred to v1.3+

**What NOT to add:**
- nanobench or custom timers (Google Benchmark's ecosystem wins on CI-readiness and parameterized ranges)
- Python benchmarking (Python RNEA is a standalone reference, not on the performance-critical path)
- CI-based benchmark runs (benchmarks are stochastic and machine-dependent — manual runs only)
- OpenMP threading (already removed from LowerTriangular; single-threaded benchmarks only)

### Expected Features

**Must have (table stakes) — P0:**
- Forward dynamics (ABA) timing with parameterized chain lengths (n=1..5+) — core metric
- Inverse dynamics (RNEA) timing — second primary use case
- Google Benchmark CMake integration via FetchContent — all benchmarks depend on this
- Random joint state per iteration (avoids warm-state bias)
- Release-mode compilation (`-O3 -DNDEBUG`) — debug builds produce meaningless numbers
- Multiple chain sizes (n=2, 3, 6, 10, 20) — establish O(n) scaling curve
- Eigen 5.x CI matrix entry — maintain compatibility guarantee
- CR-02 bug fix in ABA inward pass (BLOCKING all other features)

**Should have (competitive) — P1/P2:**
- Per-pass timing (outward vs inward) to identify optimization targets
- Round-trip consistency verification during benchmarks (RNEA∘ABA ≈ identity)
- Gravity and no-gravity benchmark modes
- Robot examples: `robot_2link_planar.cpp` and `robot_3link_spatial.cpp`
- Comparison numbers vs RBDL (published with machine spec)
- Scaling benchmarks (n=1..20 O(n) plot)

**Defer (v1.3+):**
- URDF model loading (requires external dependency; hand-coded models sufficient for v1.2)
- Memory bandwidth profiling (`perf` platform-specific)
- CI regression tracking (Google Benchmark comparison mode; infrastructure investment)
- Real-time performance tests (jitter, worst-case timing)
- Warm/cold cache benchmarking

### Architecture Approach

The core library (`include/`, `src/`, `tests/`) is **untouched** — no changes to existing classes, interfaces, or data structures. All new work lives in new directories that add functionality without modifying the library's API or ABI.

**New/modified components:**
1. **`CMakeLists.txt` (root, MODIFIED)** — Change Eigen `find_package` to range syntax `3.4...5`; add `SA_BUILD_BENCHMARKS` option (default OFF); add `cmake/` module path via `list(APPEND CMAKE_MODULE_PATH)`
2. **`benchmarks/` (NEW directory)** — Separate CMakeLists.txt with Google Benchmark FetchContent; contains `BenchmarkUtils.h` (shared model factories), `bench_aba.cpp`, `bench_rnea.cpp`, `bench_plucker.cpp`, `bench_cross_product.cpp`, and `bench_comparison.cpp` (guarded by `SA_BUILD_COMPARISON_BENCHMARKS`)
3. **`cmake/` (NEW directory)** — `FindRBDL.cmake` module (RBDL has no CMake Config mode); keeps root CMakeLists.txt clean
4. **`examples/` (MODIFIED)** — New targets `example_robot_2link` and `example_robot_3link` with physically realistic link parameters
5. **`.github/workflows/ci.yml` (MODIFIED)** — Add Eigen 5.x matrix dimension (`eigen: [3.4, 5.0]`) with manual Eigen 5.0.1 fetch/install for Ubuntu jobs

**Key architectural decisions:**
- Separate `benchmarks/` from `tests/` — different dependencies (benchmark::benchmark vs GTest), different build profile (Release vs Debug), different invocation pattern
- Google Benchmark via FetchContent with system-install fallback — no mandatory system dependency
- Comparison benchmarks as separate executables (not in-process) — each links only against its library; no symbol conflicts
- Keep `cmake_minimum_required(VERSION 3.10)` with CMake version guard for range syntax — don't force CMake upgrade on downstream users

### Critical Pitfalls

1. **CR-02 Inward Pass — Bias force (pa) initialized before child inertia accumulation (CRITICAL):** The ABA inward pass in `src/ForwardDynamics.cpp` uses a two-phase approach (init all links, then accumulate children) that computes `pa[i] = I_i × c_i` using rigid-body inertia I_i instead of the final articulated inertia I_A. Per Featherstone Algorithm 7.3, each link's bias force must be computed with the FULL I_A that includes all children. Single and 2-link chains pass (one accumulation hop), but 3+ link chains produce incorrect accelerations. **Fix: Restructure as a single pass from tip to base, computing pa_i AFTER children are accumulated.**

2. **Eigen 5.x CMake `find_package` version pin breakage (HIGH):** The current `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` rejects Eigen 5.x because CMake's default version logic interprets `3.3` as requiring `< 4.0.0`. Eigen 5.x is `5.0.0 ≥ 3.3` but fails the implicit upper bound. **Fix: Use `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` with a CMake 3.19+ guard, or drop the version pin entirely.**

3. **Benchmarking an unfixed bug (CRITICAL ordering):** Adding performance benchmarks while CR-02 is active measures "how fast the wrong answer is computed." All multi-link results are invalidated when CR-02 is fixed. **Hard rule: Benchmarks must not be added until all 4 consistency tests pass.** The acceptance gate is `TestDynamicsConsistency.cpp` showing 0 failures.

4. **Frame convention mismatch in RBDL/Pinocchio comparison (HIGH):** Each library has different transform direction APIs (parent→child vs child→parent), different floating-base conventions, and different joint screw definitions. Naive cross-validation measures "difference due to convention" rather than performance. **Prevention: Establish exact numerical equivalence (to 1e-12) on single-link before any multi-link comparison. Document the convention mapping explicitly. Only compare fixed-base models.**

5. **Small-N benchmark misleading extrapolation (MEDIUM):** Testing only 2-3 link chains makes O(n) scaling look flat. **Prevention: Always benchmark across N = {2, 3, 6, 10, 20, 50} to establish the scaling curve.**

## Implications for Roadmap

### Phase 1: Fix CR-02 ABA Inward Pass
**Rationale:** The single blocking dependency for ALL other work. Without correct multi-link dynamics, benchmarks measure wrong numbers, examples demonstrate incorrect physics, and comparison against RBDL is meaningless. Per Featherstone Algorithm 7.3, the inward pass must be restructured as a single tip→base pass where bias forces are computed using the FINAL articulated inertia I_A that includes all children.
**Delivers:** Correct forward dynamics for chains with 3+ links; all 4 consistency tests passing (ThreeLinkSerialChain, BranchingYConfiguration, TwoLinkRoundTrip, ThreeLinkNumericalValidation); existing 156 tests preserved.
**Addresses:** BLOCKING dependency for all P0 features.
**Avoids:** Pitfall 1 (CR-02 root cause), Pitfall 3 (benchmarking unfixed bug).
**Research flag:** Well-documented — Featherstone Algorithm 7.3 provides exact pseudocode. No deeper research needed.
**Confidence:** HIGH — root cause confirmed via codebase analysis, fix strategy verified against textbook.

### Phase 2: Eigen 5.x CI Compatibility
**Rationale:** Can proceed in parallel with Phase 1 (no code dependency). Prevents build breakage when Homebrew updates Eigen from 3.4.x to 5.x. The fix is a single-line CMake change plus CI matrix expansion.
**Delivers:** `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` in root CMakeLists.txt; Eigen 5.x matrix entries in CI (ubuntu × g++ × Eigen 5.0.1, macos × clang++ × Eigen 5.0.1); compile-time verification that existing Eigen usage compiles cleanly under 5.x.
**Uses:** Eigen version range syntax (CMake 3.19+); Eigen 5.0.1 source from GitLab for Linux CI; `brew install eigen@5` for macOS CI.
**Avoids:** Pitfall 2 (Eigen version pin breakage).
**Research flag:** Standard pattern — dozens of projects (PCL, VowpalWabbit) have applied the same fix. Skip research-phase.
**Confidence:** HIGH — exact fix documented by Eigen team and adopted by PCL (#6354).

### Phase 3: Benchmark Infrastructure Setup
**Rationale:** Requires Phases 1-2 complete (correct solver + working build with Eigen 5.x). Creates the `benchmarks/` directory, CMakeLists.txt with Google Benchmark FetchContent, and shared utility headers.
**Delivers:** `benchmarks/CMakeLists.txt` with Google Benchmark v1.9.5 via FetchContent; `benchmarks/BenchmarkUtils.h` (model factory, random state generator); `benchmarks/BenchmarkModels.h` (canonical test chains); root CMakeLists.txt option `SA_BUILD_BENCHMARKS` (default OFF); build-verifies in CI.
**Addresses:** P0 feature "Google Benchmark CMake integration."
**Uses:** Google Benchmark v1.9.5, Eigen 3.4...5 range syntax, existing SpatialAlgebra library.
**Implements:** `benchmarks/` directory architecture.
**Avoids:** Pitfall 4 (mixing benchmark/test targets), Pitfall 5 (single-size benchmarks — model factory supports arbitrary n).
**Research flag:** Standard CMake FetchContent + Google Benchmark pattern. Well-documented. Skip research-phase.

### Phase 4: Benchmark Implementation (ABA, RNEA, Microbenchmarks)
**Rationale:** Depends on Phase 3 infrastructure being in place. Implements the actual benchmark executables with parameterized DOF ranges (n=1..20), random joint state per iteration, gravity/no-gravity modes, and round-trip verification.
**Delivers:** `bench_aba.cpp` (forward dynamics timing, n=2,4,8,16 via `->Range()`); `bench_rnea.cpp` (inverse dynamics timing, same range); `bench_plucker.cpp` and `bench_cross_product.cpp` (microbenchmarks); JSON-formatted benchmark output; round-trip consistency check integrated.
**Addresses:** P0 features "ABA timing," "RNEA timing"; P1 features "gravity modes," "round-trip check," "multiple chain sizes."
**Uses:** Google Benchmark `->Range()`, `DoNotOptimize()`, `PauseTiming()`/`ResumeTiming()`.
**Implements:** Parameterized benchmark pattern from ARCHITECTURE.md.
**Avoids:** Pitfall 3 (benchmarking unfixed bug — Phase 1 gates this), Pitfall 5 (small-N misleading extrapolation — tests n=1..20), performance traps (debug mode, dead-code elimination, setup-in-timing).
**Research flag:** Well-documented — Google Benchmark user guide provides all patterns. Skip research-phase.

### Phase 5: Robot Examples
**Rationale:** Depends on Phase 1 (correct solver). Adds real-world applicability demonstrations that produce physically correct output. No dependency on benchmark infrastructure (can proceed in parallel with Phases 3-4).
**Delivers:** `examples/robot_2link_planar.cpp` (Z-Z revolute, XY planar motion, gravity + trajectory); `examples/robot_3link_spatial.cpp` (Z-Y-Z revolute, 3D RRR arm); updated `examples/CMakeLists.txt` with new targets; validation that examples produce physically correct output.
**Addresses:** P1 features "real-world examples."
**Implements:** Robot example architecture pattern from ARCHITECTURE.md.
**Avoids:** Pitfall 1 (CR-02 buggy solver in examples — gated by Phase 1).
**Research flag:** Standard pattern — evolves existing `examples/dynamics.cpp`. Skip research-phase.

### Phase 6: RBDL Comparison Benchmarks (Optional)
**Rationale:** Depends on all prior phases. Requires installing RBDL v3.3.1, building equivalent chain models, validating convention equivalence to 1e-12, then running performance comparison. Guarded by `SA_BUILD_COMPARISON_BENCHMARKS` (default OFF).
**Delivers:** `bench_comparison.cpp` (guarded); `cmake/FindRBDL.cmake` module; convention mapping documentation; published comparison numbers (machine spec + benchmark methodology).
**Addresses:** P2 feature "comparison vs RBDL."
**Uses:** RBDL v3.3.1, Google Benchmark JSON comparison.
**Avoids:** Pitfall 4 (frame convention mismatch — requires exact numerical validation before comparison).
**Research flag:** HIGH — needs convention mapping research, RBDL API verification, and equivalence test design. Recommend `/gsd-plan-phase --research-phase` during planning.

### Phase Ordering Rationale

```
Phase 1 (CR-02 Fix) ──────────────────────┐
                                          │
Phase 2 (Eigen 5.x CI) ── parallel ───────┤
                                          │
         ┌─────────────────────────────────┘
         ▼
Phase 3 (Benchmark Infra) ───┐
                             │
Phase 4 (Benchmarks) ────────┤
                             │
Phase 5 (Examples) ──────────┤  parallel
                             │
Phase 6 (RBDL Comparison) ───┘
```

- **Phases 1-2 are independent** — CR-02 fix and Eigen 5.x CI have no code dependency on each other and can execute in parallel
- **Phases 3-6 all depend on Phase 1** — without correct multi-link dynamics, none of these produce meaningful output
- **Phases 3 vs 5 are independent** — benchmarks and examples have no dependency on each other
- **Phase 6 depends on everything** — comparison benchmarks need correct solver (Phase 1), benchmark infrastructure (Phase 3), and RBDL integration (new dependency)

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 6 (RBDL Comparison):** Frame convention mapping between SpatialAlgebra, RBDL, and Pinocchio requires careful API analysis. Different transform directions, joint screw conventions, and floating-base handling need exact documentation. Recommend `/gsd-plan-phase --research-phase`.

Phases with standard patterns (skip research-phase):
- **Phase 1 (CR-02 Fix):** Featherstone Algorithm 7.3 provides exact pseudocode. Root cause confirmed. Straightforward implementation.
- **Phase 2 (Eigen 5.x CI):** Industry-standard fix; dozens of projects have applied it. Single CMake line + CI matrix expansion.
- **Phase 3 (Benchmark Infra):** Google Benchmark + FetchContent is a well-documented pattern.
- **Phase 4 (Benchmarks):** Google Benchmark user guide provides all needed patterns.
- **Phase 5 (Examples):** Evolves existing `examples/dynamics.cpp` pattern.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Google Benchmark v1.9.5 is released and documented; Eigen 5.x range syntax is official; RBDL v3.3.1 is stable. All sources are primary (official repos, release notes, CMake docs). |
| Features | HIGH | Pinocchio and RBDL benchmark methodologies are well-documented and directly applicable. Feature priorities validated against both reference libraries' approaches. |
| Architecture | HIGH | Separate `benchmarks/` directory follows industry standard. Eigen 5.x CMake guard pattern verified against PCL and VowpalWabbit PRs. No core library modification needed. |
| Pitfalls | HIGH | CR-02 root cause confirmed via codebase analysis of `src/ForwardDynamics.cpp`. Eigen version pin breakage documented by PCL #6351. Frame convention mismatch documented in Pinocchio issues. |

**Overall confidence:** HIGH — all research areas have primary sources (official docs, codebase analysis, published benchmarks, existing PRs). The two areas with MEDIUM confidence are: (1) expected performance ranges (inferred from RBDL 2012 numbers scaled to modern hardware — needs actual measurement), and (2) RBDL/Pinocchio convention details (need API-level verification during Phase 6).

### Gaps to Address

- **Expected performance ranges:** Estimated from RBDL's 2012 i7-920 numbers scaled ~10× for modern Apple M-series. These are rough estimates and must be replaced with actual measurements during Phase 4. The degradation threshold (>1 µs/DOF) is a reasonable starting point but may need adjustment.
- **RBDL `FindRBDL.cmake` module:** RBDL does NOT provide CMake Config mode. The research identified that a custom module is needed but did not write one. Phase 6 planning must include writing `cmake/FindRBDL.cmake` based on RBDL's example `CMakeLists.txt`.
- **Pinocchio API specifics:** Research recommends against using Pinocchio for v1.2, but if v1.3 reconsiders, the exact API mapping (JointModel variants, SE3 act vs SpatialAlgebra transform) needs deeper investigation.
- **Floating base decision:** SpatialAlgebra is fixed-base only. All comparison benchmarks must explicitly document this. Phase 6 planning must ensure RBDL/Pinocchio benchmarks use fixed-base models too.
- **Compiler autovectorization differences:** Clang vs GCC can produce 2× differences on Eigen expressions. The benchmark methodology must report per-compiler results separately. This should be documented in `BENCHMARKS.md` during Phase 4.

## Sources

### Primary (HIGH confidence)
- Google Benchmark v1.9.5 — GitHub releases, official CMake integration guide
- Eigen 5.0 release notes — libeigen.gitlab.io, official CMake range syntax documentation
- Eigen 5.0.1 tags — GitLab repository tags
- RBDL v3.3.1 — GitHub releases, DeepWiki mathematical foundation
- CR-02 root cause — Codebase analysis: `src/ForwardDynamics.cpp:56-93`
- Existing CI matrix — `.github/workflows/ci.yml` (current: 4-matrix, no Eigen 5.x)
- Existing test gaps — `TestDynamicsConsistency.cpp` (156/158 passing, 3 CR-02 failures)
- Featherstone Algorithm 7.3 — Featherstone (2008) Rigid Body Dynamics Algorithms, §7.2.1
- Eigen 5.x version detection pattern — PCL PR #6354, VowpalWabbit PR #4728
- Google Benchmark optimization barriers — Official docs

### Secondary (MEDIUM confidence)
- RBDL benchmarking architecture — DeepWiki documentation
- Pinocchio vs RBDL performance paper — HAL archives (2018 IROS paper)
- Pinocchio Eigen dependency tree — GitHub source inspection
- RBDL/Pinocchio convention mismatch — Pinocchio Issue #1721
- Expected performance ranges — Inferred from RBDL 2012 i7-920 numbers scaled for modern hardware
- Frame convention mapping — Cross-referenced between Featherstone, RBDL docs, and Pinocchio docs

---

*Research completed: 2026-05-30*
*Ready for roadmap: yes*
