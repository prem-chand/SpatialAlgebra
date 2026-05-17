---
phase: 13-production-readiness
plan: 05
subsystem: build
tags: cmake, gtest, fetchcontent, github-actions, ci, namespace, openmp, umbrella-header

requires:
  - phase: 13-01
    provides: Research on namespace cleanup, CMake changes, CI patterns

provides:
  - Namespace pollution fix: `using Vector3d` moved from global scope into `namespace SpatialAlgebra`
  - OpenMP dependency removal from LowerTriangular.h
  - Umbrella header `include/SpatialAlgebra.h` for single-include access
  - CMake FetchContent fallback for GTest (no system-installed GTest required)
  - CMake ENABLE_COVERAGE option for CI coverage builds
  - GitHub Actions CI workflow with 4-matrix build (ubuntu/macos × g++/clang++)
  - CodeCov upload on ubuntu+g++ builds
  - Empty stub files deleted from compilation

affects:
  - downstream code must use `SpatialAlgebra::Vector3d` instead of global `::Vector3d`
  - users of `using namespace SpatialAlgebra;` are unaffected

tech-stack:
  added:
    - GitHub Actions CI (`.github/workflows/ci.yml`)
    - codecov/codecov-action@v4 for coverage upload
  patterns:
    - Use `#include "SpatialAlgebra.h"` as a single-include entry point
    - `Vector3d` canonical declaration in `SpatialVector.h` inside `namespace SpatialAlgebra`
    - CMake FetchContent for optional test dependency fallback
    - CI workflow for automated build verification on push/PR

key-files:
  created:
    - include/SpatialAlgebra.h (umbrella header, 12 includes in dependency order)
    - .github/workflows/ci.yml (CI workflow, 4-matrix, coverage upload)
    - .planning/phases/13-production-readiness/13-05-deferred-items.md
  modified:
    - include/SpatialVector.h (Vector3d moved into namespace)
    - include/LowerTriangular.h (removed duplicate Vector3d, removed OpenMP pragma, added SpatialVector.h include)
    - CMakeLists.txt (FetchContent fallback, stub removal, coverage option)
  deleted:
    - src/RigidBodyInertia.cpp (empty stub)
    - src/ArticulatedBodyInertia.cpp (empty stub)

key-decisions:
  - "LowerTriangular.h must include SpatialVector.h to access Vector3d typedef (build fix — not in plan)"
  - ".github/ is in .gitignore — CI workflow must be force-added on future changes"
  - "Eigen3 already had version-pin removed — no change needed (D-18 already applied)"

patterns-established:
  - "Canonical Vector3d declared once in SpatialVector.h inside namespace SpatialAlgebra"
  - "Umbrella header SpatialAlgebra.h for single-include convenience"
  - "CMake FetchContent fallback when system GTest is not installed"
  - "ENABLE_COVERAGE option for CI-only coverage flags"
  - "CI workflow with push/PR triggers, 4-matrix, coverage on ubuntu+g++"

requirements-completed:
  - VEC-01

duration: ~45min
completed: 2026-05-17
---

# Phase 13 Plan 05: Code Quality, Build Improvements & CI Setup Summary

**Namespace cleanup (Vector3d into SpatialAlgebra), OpenMP removal, umbrella header, CMake GTest FetchContent fallback, stub file removal, and GitHub Actions CI workflow with 4-matrix build and code coverage upload**

## Performance

- **Duration:** ~45 min
- **Started:** 2026-05-17
- **Completed:** 2026-05-17
- **Tasks:** 3 (all auto)
- **Files modified:** 8 (3 created, 3 modified, 2 deleted)

## Accomplishments

- **Namespace cleanup (D-23):** `using Vector3d` moved from global scope into `namespace SpatialAlgebra` in `SpatialVector.h`. Duplicate removed from `LowerTriangular.h`. Eliminates ODR hazard and global namespace pollution.
- **OpenMP removal (D-26):** `#pragma omp parallel for collapse(2)` removed from `LowerTriangular::operator*`. Removes hidden linkage dependency on OpenMP runtime.
- **Umbrella header (D-28):** `include/SpatialAlgebra.h` created with 12 public headers in dependency order using `#ifndef`/`#define` include guard.
- **Empty stub removal (D-27):** `src/RigidBodyInertia.cpp` and `src/ArticulatedBodyInertia.cpp` deleted. CMake `list(REMOVE_ITEM)` excludes them from the source glob.
- **GTest FetchContent fallback (D-25):** `find_package(GTest QUIET)` with FetchContent fallback to `release-1.12.1.zip`. Library builds without system-installed GTest.
- **Coverage option (D-11):** `ENABLE_COVERAGE` option (default OFF) adds `--coverage` flags when enabled.
- **GitHub Actions CI (D-09 through D-12):** 4-matrix workflow (ubuntu/macos × g++/clang++) with cmake configure/build/test. Coverage upload on ubuntu+g++ via `codecov/codecov-action@v4`.

## Task Commits

Each task was committed atomically:

1. **Task 1: Namespace cleanup, OpenMP removal, umbrella header, empty stub removal** — `38e6915` (fix)
2. **Task 2: CMake changes — GTest FetchContent fallback, stub removal, coverage flags** — `4316cd8` (fix)
3. **Task 3: Create GitHub Actions CI workflow** — `cc8c0a5` (feat)

## Files Created/Modified

- `include/SpatialVector.h` — `using Vector3d` moved inside `namespace SpatialAlgebra`
- `include/LowerTriangular.h` — Removed duplicate `Vector3d`, removed `#pragma omp`, added `#include "SpatialVector.h"`
- `include/SpatialAlgebra.h` — **NEW** Umbrella header with 12 includes
- `CMakeLists.txt` — GTest FetchContent fallback, stub exclusion, coverage option
- `.github/workflows/ci.yml` — **NEW** CI workflow with 4-matrix build and coverage upload
- `src/RigidBodyInertia.cpp` — **DELETED** (empty stub)
- `src/ArticulatedBodyInertia.cpp` — **DELETED** (empty stub)
- `.planning/phases/13-production-readiness/13-05-deferred-items.md` — **NEW** Pre-existing test failure log

## Decisions Made

- **LowerTriangular.h must include SpatialVector.h:** The plan assumed `LowerTriangular.h` always received `Vector3d` transitively, but `LowerTriangular.cpp` includes only `LowerTriangular.h`, which doesn't (and shouldn't) depend on `SpatialVector.h`. Added explicit `#include "SpatialVector.h"` to resolve the dependency — a deviation from the plan (Rule 3 — blocking).
- **Force-add CI workflow:** `.github/` is in `.gitignore`, so the CI file must be force-added (`git add -f`).
- **Eigen3 version pin already removed:** Line 12 already says `find_package(Eigen3 REQUIRED NO_MODULE)` without a version — D-18 was already applied before this phase.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 — Blocking] Add missing include to LowerTriangular.h**
- **Found during:** Task 1 (post-edit build verification)
- **Issue:** Removing the `using Vector3d` from `LowerTriangular.h` (inside `namespace SpatialAlgebra`) broke the build because `LowerTriangular.h` uses `Vector3d` in its public API (e.g., `operator*(const Vector3d &v)`), and `LowerTriangular.cpp` includes only `LowerTriangular.h` — not `SpatialVector.h`.
- **Fix:** Added `#include "SpatialVector.h"` to `LowerTriangular.h` so it has access to the canonical `SpatialAlgebra::Vector3d` declaration. No circular dependency — `SpatialVector.h` does not include `LowerTriangular.h`.
- **Files modified:** `include/LowerTriangular.h`
- **Verification:** Full clean build succeeds (0 errors), all tests link and run.
- **Committed in:** `4316cd8` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking fix)
**Impact on plan:** The fix is necessary for the build to compile. No scope creep — `LowerTriangular.h` genuinely depends on `SpatialVector.h` for the `Vector3d` type used in its interface, and the plan's assumption of transitive inclusion was incorrect for direct compilation units.

## Issues Encountered

- **Pre-existing test failures in TestDynamicsConsistency:** 3 sub-tests fail (`ThreeLinkSerialChain`, `BranchingYConfiguration`, `TwoLinkRoundTripWithGravity`) due to ABA/RNEA roundtrip accuracy for multi-link chains with gravity. These are pre-existing dynamics algorithm issues not caused by phase 13 changes. Logged in `13-05-deferred-items.md`.
- **.github/ in .gitignore:** The `.github/` directory is gitignored, so the CI workflow file must be force-added. This is a project-level config that should be updated if CI is now an expected part of the repo.

## Deferred Items

- **TestDynamicsConsistency failures:** 3 pre-existing dynamics algorithm accuracy issues. See `13-05-deferred-items.md` for details. These should be addressed in a dedicated algorithm bug-fix phase.

## Next Phase Readiness

- Code quality fixes complete (namespace, OpenMP, umbrella header, stub cleanup)
- Build system improvements complete (FetchContent fallback, coverage option)
- CI workflow ready for first run on next push/PR to main
- Ready for next plan in phase 13

## Self-Check: PASSED

**File existence:**
- `include/SpatialAlgebra.h` — ✅ FOUND
- `.github/workflows/ci.yml` — ✅ FOUND
- `include/SpatialVector.h` — ✅ FOUND (modified)
- `include/LowerTriangular.h` — ✅ FOUND (modified)
- `CMakeLists.txt` — ✅ FOUND (modified)
- `.planning/phases/13-production-readiness/13-05-SUMMARY.md` — ✅ FOUND

**Stub deletions:**
- `src/RigidBodyInertia.cpp` — ✅ DELETED
- `src/ArticulatedBodyInertia.cpp` — ✅ DELETED

**Commits verified:**
- `38e6915` — Task 1: namespace cleanup, OpenMP removal, umbrella header, stub removal
- `4316cd8` — Task 2: CMake changes (FetchContent, stub exclusion, coverage)
- `cc8c0a5` — Task 3: CI workflow creation
- `293d288` — SUMMARY.md metadata commit

**All 4 commits present and accounted for.**

---
*Phase: 13-production-readiness*
*Completed: 2026-05-17*
