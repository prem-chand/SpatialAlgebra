---
phase: 16-benchmark-infrastructure
plan: 01
subsystem: build-system
tags: cmake, fetchcontent, google-benchmark, benchmarks, lto

# Dependency graph
requires:
  - phase: 15-eigen5-ci
    provides: CMake 3.19+ compatible build system with FetchContent pattern
provides:
  - SA_BUILD_BENCHMARKS guard (default OFF) in root CMakeLists.txt
  - benchmarks/ directory with per-domain subdirectories (aba/, rnea/, core/, common/)
  - Google Benchmark v1.9.5 integrated via FetchContent in benchmarks/CMakeLists.txt
  - bench_all executable target with bench_common static utility library
  - LTO optimization (-O3 -DNDEBUG -flto) on all benchmark targets
affects: [17-aba-benchmarks, 18-rnea-benchmarks, 19-core-benchmarks]

# Tech tracking
tech-stack:
  added: [Google Benchmark v1.9.5 via FetchContent]
  patterns:
    - Per-domain subdirectory CMakeLists.txt with PARENT_SCOPE source variables
    - FetchContent scoped to benchmarks/CMakeLists.txt (not root)
    - -O3 -DNDEBUG per-target compilation (not CMAKE_BUILD_TYPE override)
    - LTO enabled via target_compile_options / target_link_options

key-files:
  created:
    - benchmarks/CMakeLists.txt — FetchContent + bench_all target
    - benchmarks/common/CMakeLists.txt — bench_common STATIC library
    - benchmarks/aba/CMakeLists.txt — ABA_SOURCES PARENT_SCOPE
    - benchmarks/rnea/CMakeLists.txt — RNEA_SOURCES PARENT_SCOPE
    - benchmarks/core/CMakeLists.txt — CORE_SOURCES PARENT_SCOPE
    - benchmarks/_init.cpp — CMake 4.x placeholder (removed by Plan 03)
    - benchmarks/common/_init.cpp — CMake 4.x placeholder (removed by Plan 02)
  modified:
    - CMakeLists.txt — added SA_BUILD_BENCHMARKS guard block

key-decisions:
  - "FetchContent for Google Benchmark in benchmarks/CMakeLists.txt (scoped), not root CMakeLists.txt"
  - "Per-domain subdirectories with PARENT_SCOPE source collection variables"
  - "init.cpp placeholders required because CMake 4.2.0 does not allow add_executable/add_library with zero sources"
  - "target_compile_options(-O3 -DNDEBUG -flto) applied per-target, not as CMAKE_BUILD_TYPE override"
  - "Single bench_all executable with bench_common static library for shared utilities"

patterns-established:
  - "Benchmark subdirectory pattern: file(GLOB) + set(VAR PARENT_SCOPE) for source collection"
  - "Common utility library as separate STATIC target (bench_common) to avoid ODR issues"
  - "LTO on benchmark targets only, not project-wide"

requirements-completed: [BINF-01, BINF-02]

# Metrics
duration: 149 min
completed: 2026-06-05
---

# Phase 16: Benchmark Infrastructure — Plan 01 Summary

**Google Benchmark v1.9.5 FetchContent integration, SA_BUILD_BENCHMARKS guard (default OFF), and full benchmarks/ directory structure with per-domain CMakeLists.txt files**

## Performance

- **Duration:** 149 min (includes Google Benchmark library download and build)
- **Started:** 2026-06-05T10:45:21Z
- **Completed:** 2026-06-05T13:14:27Z
- **Tasks:** 3
- **Files modified:** 9

## Accomplishments

- Added `option(SA_BUILD_BENCHMARKS OFF)` guard to root CMakeLists.txt, wrapping `add_subdirectory(benchmarks)` — benchmarks are opt-in with zero impact on default builds
- Integrated Google Benchmark v1.9.5 via FetchContent in `benchmarks/CMakeLists.txt` with `BENCHMARK_ENABLE_TESTING OFF` and `BENCHMARK_ENABLE_INSTALL OFF` to prevent transitive GTest dependency conflict
- Created four per-domain subdirectories (`common/`, `aba/`, `rnea/`, `core/`) each with its own CMakeLists.txt using `file(GLOB)` and `PARENT_SCOPE` source collection
- Created `bench_all` executable target and `bench_common` static library, both with `-O3 -DNDEBUG -flto` optimization profile and LTO linking
- CMake configure with `-DSA_BUILD_BENCHMARKS=ON` succeeds — Google Benchmark downloads, all targets created, and `bench_all` compiles with placeholder sources

## Task Commits

Each task was committed atomically:

1. **Task 1: Add SA_BUILD_BENCHMARKS guard to root CMakeLists.txt** — `2e18201` (feat)
2. **Task 2: Create benchmarks/ directory and benchmarks/CMakeLists.txt** — `18ddc0a` (feat)
3. **Task 3: Create subdirectory CMakeLists.txt files (common, aba, rnea, core)** — `5f33ebd` (feat)

**Plan metadata:** Pending SUMMARY.md commit

## Files Created/Modified

- `CMakeLists.txt` — Inserted `option(SA_BUILD_BENCHMARKS OFF)` + guarded `add_subdirectory(benchmarks)` block
- `benchmarks/CMakeLists.txt` — FetchContent for Google Benchmark v1.9.5, subdirectory includes, bench_all target, LTO flags
- `benchmarks/common/CMakeLists.txt` — `bench_common` STATIC library via `file(GLOB)` with `SpatialAlgebra` and `Eigen3::Eigen` linkage
- `benchmarks/aba/CMakeLists.txt` — `file(GLOB)` + `ABA_SOURCES PARENT_SCOPE`
- `benchmarks/rnea/CMakeLists.txt` — `file(GLOB)` + `RNEA_SOURCES PARENT_SCOPE`
- `benchmarks/core/CMakeLists.txt` — `file(GLOB)` + `CORE_SOURCES PARENT_SCOPE`
- `benchmarks/_init.cpp` — Placeholder source for `bench_all` target (CMake 4.x requires ≥1 source; removed by Plan 03)
- `benchmarks/common/_init.cpp` — Placeholder source for `bench_common` target (removed by Plan 02)

## Decisions Made

- **FetchContent scoping:** Google Benchmark FetchContent lives in `benchmarks/CMakeLists.txt`, keeping root CMakeLists.txt clean and ensuring benchmark dependencies are only fetched when `SA_BUILD_BENCHMARKS=ON`
- **Source collection pattern:** Each subdirectory CMakeLists.txt uses `file(GLOB)` + `PARENT_SCOPE` variable (e.g., `ABA_SOURCES`), collected in parent as `BENCH_SOURCES`
- **CMake 4.x compatibility:** Requires `_init.cpp` placeholder files because CMake 4.2.0 does not allow `add_executable()` or `add_library()` with zero sources (plan incorrectly assumed CMake 3.19+ would allow it)
- **Optimization approach:** `-O3 -DNDEBUG -flto` per-target via `target_compile_options()` and `target_link_options()` rather than `CMAKE_BUILD_TYPE` override, preserving parent build type inheritance

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] CMake 4.x requires at least one source per target**
- **Found during:** Task 2 verification (cmake configuration)
- **Issue:** CMake 4.2.0 errors with "No SOURCES given to target" when `add_executable()` or `add_library()` receives zero source files. The plan's assertion that "CMake allows add_executable with empty source lists in 3.19+" is incorrect for CMake 4.2.0. Both `bench_all` and `bench_common` had no `.cpp` sources during Plan 01 (sources are added by Plans 02 and 03).
- **Fix:** Created minimal `_init.cpp` placeholder files in `benchmarks/` (for `bench_all`) and `benchmarks/common/` (for `bench_common`), and added `if(NOT BENCH_SOURCES)` fallback in `benchmarks/CMakeLists.txt` to use the placeholder when subdirectories produce no sources.
- **Files modified:** `benchmarks/CMakeLists.txt`, `benchmarks/common/CMakeLists.txt`, `benchmarks/_init.cpp` (created), `benchmarks/common/_init.cpp` (created)
- **Verification:** `cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON && cmake --build build_bench --target bench_all` succeeds — Google Benchmark library compiles, placeholder compiles, bench_all links
- **Committed in:** `5f33ebd` (Task 3 commit, placeholder files added alongside subdirectory CMakeLists.txt)

---

**Total deviations:** 1 auto-fixed (1 blocking fix for CMake 4.x compatibility)
**Impact on plan:** Minor — two tiny placeholder files added to the benchmarks/ tree. They will be removed automatically by CMake GLOB once Plans 02 and 03 add real source files. No scope creep — the plan's CMake architecture is preserved unchanged.

## Issues Encountered

- **CMake 4.2.0 empty source restriction:** The plan assumed `add_executable(name ${empty_list})` works in CMake 3.19+. Testing confirmed CMake 4.2.0 errors on this pattern even with `cmake_policy(SET CMP0112 NEW)`. Fixed with placeholder source files as documented above.
- **Google Benchmark FetchContent duration:** First `cmake -B build_bench` took ~25 seconds to clone and build Google Benchmark's library. Subsequent runs reuse local CMake cache and are nearly instant.

## Next Phase Readiness

- Build system skeleton complete — Plan 02 can add `model_factory.cpp` and `random_state.cpp` to `benchmarks/common/` without touching any CMake files
- Plan 03 can add `bench_all.cpp` and domain stubs to `benchmarks/aba/`, `benchmarks/rnea/`, `benchmarks/core/` — the CMake `file(GLOB)` will auto-discover new `.cpp` files on re-configure
- Default build (`SA_BUILD_BENCHMARKS=OFF` by default) is unaffected — no benchmark artifacts leak into normal builds
- **Note:** The `_init.cpp` placeholders in `benchmarks/` and `benchmarks/common/` should be removed by Plans 02 and 03 respectively when real source files are added

---
*Phase: 16-benchmark-infrastructure*
*Plan: 01*
*Completed: 2026-06-05*
