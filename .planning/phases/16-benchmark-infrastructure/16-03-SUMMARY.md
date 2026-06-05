---
phase: 16-benchmark-infrastructure
plan: 03
subsystem: testing
tags: google-benchmark, benchmark-infrastructure, stubs
requires:
  - phase: 16-benchmark-infrastructure
    provides: benchmarks/CMakeLists.txt, bench_common library, domain directories
provides:
  - bench_all.cpp with programmatic RegisterBenchmark DOF sweep
  - Domain stub files for ABA, RNEA, and core microbenchmarks
  - Full linkable bench_all executable for Phase 17 benchmark implementation
affects:
  - 17-benchmark-aba
  - 18-benchmark-rnea
  - 19-benchmark-comparison
tech-stack:
  added: []
  patterns:
    - Programmatic RegisterBenchmark with DOF sweep (n=1..20)
    - Stub functions with (void)state;(void)nDOF; suppress pattern
    - target_sources() for post-add_executable source appending
key-files:
  created:
    - benchmarks/aba/bench_aba_stub.cpp
    - benchmarks/rnea/bench_rnea_stub.cpp
    - benchmarks/core/bench_core_stub.cpp
    - benchmarks/bench_all.cpp
  modified:
    - benchmarks/CMakeLists.txt
    - benchmarks/common/model_factory.cpp
key-decisions:
  - "Stub body uses (void)state;(void)nDOF; suppress pattern — empty stubs would otherwise trigger Google Benchmark's 'did not run' diagnostic"
  - "bench_all.cpp added via target_sources() after add_executable() because file didn't exist when BENCH_SOURCES was assembled in Plan 01 (CMake 4.x validates sources at configure time)"
requirements-completed:
  - BINF-03
duration: 28 min
completed: 2026-06-05
---

# Phase 16 Plan 03: Benchmark Stub Files + Entry Point Summary

**Programmatic RegisterBenchmark DOF sweep stubs for ABA, RNEA, and core microbenchmarks, producing a linkable bench_all executable with 80 registered benchmarks**

## Performance

- **Duration:** 28 min
- **Started:** 2026-06-05T13:19:54Z
- **Completed:** 2026-06-05T13:48:04Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments

- Created 3 domain stub files (ABA, RNEA, core microbenchmarks) with correct function signatures, includes, and Doxygen comments
- Created bench_all.cpp with custom main(), 4 forward declarations, and DOF sweep n=1..20 for all benchmark types
- Modified benchmarks/CMakeLists.txt to append bench_all.cpp via target_sources() (CMake 4.x validation constraint)
- Verified full build produces a linkable bench_all executable
- Verified 80 benchmarks registered and filterable via --benchmark_filter
- Verified SA_BUILD_BENCHMARKS=OFF excludes bench_all from normal builds (no regression)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create domain benchmark stub files** - `184ddc8` (feat)
2. **Task 2: Create bench_all.cpp with RegisterBenchmark** - `7b0c17d` (feat)
3. **Deviation: Fix ForwardDynamics::Link in model_factory** - `4571581` (fix)

**Plan metadata:** (committed below as SUMMARY)

## Files Created/Modified

- `benchmarks/aba/bench_aba_stub.cpp` - ABA forward dynamics stub (BM_ABA_ForwardDynamics)
- `benchmarks/rnea/bench_rnea_stub.cpp` - RNEA inverse dynamics stub (BM_RNEA_InverseDynamics)
- `benchmarks/core/bench_core_stub.cpp` - Core microbenchmarks stub (BM_PluckerTransform, BM_CrossProduct)
- `benchmarks/bench_all.cpp` - Entry point with DOF sweep registration for all 4 benchmark types
- `benchmarks/CMakeLists.txt` - Added target_sources(bench_all PRIVATE bench_all.cpp)
- `benchmarks/common/model_factory.cpp` - Fixed `ForwardDynamics::Link` → `Link` (3 occurrences, blocking build)

## Decisions Made

- **Stub pattern**: `(void)state; (void)nDOF;` suppresses unused parameter warnings but leaves empty benchmark bodies. Phase 17 will fill with actual `for (auto _ : state)` loops.
- **target_sources pattern**: bench_all.cpp is appended via `target_sources()` (not added to the initial BENCH_SOURCES) because CMake 4.x validates source files at configure time and the file didn't exist when `add_executable(bench_all)` ran in Plan 01.
- **Forward declarations in bench_all.cpp**: All 4 benchmark function signatures exposed as `void BM_*(benchmark::State&, int)` — these must remain in sync with stub definitions.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed ForwardDynamics::Link qualification in model_factory.cpp**
- **Found during:** Task 2 (build verification)
- **Issue:** `model_factory.cpp` used `ForwardDynamics::Link` which fails to compile because `Link` is a namespace-scope struct in `SpatialAlgebra`, not a nested type inside `ForwardDynamics`. This pre-existing bug from Plan 16-02 blocked the entire bench_all build.
- **Fix:** Changed `ForwardDynamics::Link` → `Link` in all 3 FD factory methods.
- **Files modified:** benchmarks/common/model_factory.cpp
- **Verification:** Full bench_all build succeeds with 0 errors.
- **Committed in:** `4571581` (fix commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Fix was essential for build to succeed. No scope creep — corrected a type qualification error in pre-existing code.

## Issues Encountered

- None — the plan executed cleanly after the auto-fixed compilation error in model_factory.cpp.

## Next Phase Readiness

- bench_all executable ready for Phase 17 benchmark implementation
- 80 registered benchmarks (4 domains × 20 DOF values) ready to receive real logic
- Stubs provide correct function signatures — Phase 17 fills `for (auto _ : state)` loops with actual solver calls, model factory usage, and random state generation
- bench_all can be tested now with `--benchmark_filter` to verify Phase 17 changes incrementally

---

*Phase: 16-benchmark-infrastructure*
*Completed: 2026-06-05*
