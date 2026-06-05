---
phase: 15-eigen-5-x-ci
plan: 01
subsystem: infra
tags: [eigen, cmake, ci, github-actions, build-system]

requires:
  - phase: 13-production-readiness
    provides: Existing 4-matrix CI workflow and CMakeLists.txt foundation

provides:
  - Eigen version range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` accepting both 3.4.x and 5.x
  - 8-job CI matrix (2 OS × 2 compilers × 2 Eigen versions) with fail-fast disabled
  - Local Eigen 5.0.1 compilation and test verification — zero errors, zero warnings, all 11 tests pass
  - Updated AGENTS.md with version range syntax guidance

affects: [build-system, ci-configuration]

tech-stack:
  added: []
  patterns:
    - "Version range syntax: find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)"
    - "CI matrix expansion with eigen-version dimension"

key-files:
  created: []
  modified:
    - CMakeLists.txt
    - .github/workflows/ci.yml
    - AGENTS.md

key-decisions:
  - "Used version range syntax 3.4...5 instead of bare version pin (D-01)"
  - "Built Eigen from source on macOS for both versions; apt for Ubuntu 3.4, source for 5.0 (D-03)"
  - "Coverage gated on eigen-version == '3.4' (D-05)"

patterns-established:
  - "CMake version range syntax for Eigen dependency declaration"
  - "CI matrix with OS × compiler × eigen-version for cross-product testing"

requirements-completed: [CI-01, CI-02]

duration: 16 min
completed: 2026-06-05
---

# Phase 15: Eigen 5.x CI — Summary

**Expanded CI matrix from 4 to 8 jobs with Eigen version range syntax, verified zero-error compilation and all 11 tests passing under Eigen 5.0.1**

## Performance

- **Duration:** 16 min
- **Started:** 2026-06-05T05:00:00Z
- **Completed:** 2026-06-05T05:16:00Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Updated CMakeLists.txt with `cmake_minimum_required(VERSION 3.19)` and `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` — accepts both Eigen 3.4.x and 5.x
- Updated AGENTS.md to reflect the version range syntax solution, replacing the old version mismatch workaround
- Built Eigen 5.0.1 from source locally and verified the entire project compiles with zero errors and zero warnings under Eigen 5.x
- All 11 CTest tests pass under Eigen 5.0.1 — no source code changes needed (Eigen 5.x is fully backward-compatible for this codebase)
- Expanded CI matrix to 8 jobs with `eigen-version: [3.4, 5.0]`, `fail-fast: false`, and per-version install steps
- Coverage upload gated on `eigen-version == '3.4'`

## Task Commits

1. **Task 1: Update CMakeLists.txt with version range syntax and bump minimum CMake version** - `faadfd8` (build)
2. **Task 2: Verify compilation under Eigen 5.x and expand CI matrix to 8 jobs** - `928aae1` (ci)

**Plan metadata:** (commits below)

## Files Created/Modified

- `CMakeLists.txt` — Bumped CMake minimum to 3.19, added version range syntax to `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)`
- `.github/workflows/ci.yml` — Expanded to 8 jobs with eigen-version dimension, conditional install steps, fail-fast disabled
- `AGENTS.md` — Replaced Eigen version mismatch workaround with version range syntax guidance

## Decisions Made

- Followed all D-01 through D-05 from the phase context as locked decisions
- Eigen 5.x installed from source on both macOS and Ubuntu for determinism (per D-03)
- Eigen 3.4 from apt on Ubuntu, from source on macOS (Homebrew no longer has a 3.4 tap)
- Coverage gated on eigen-version == '3.4' per D-05

## Deviations from Plan

None — plan executed exactly as written. Eigen 5.0.1 compiled the entire codebase with zero errors and zero warnings, confirming no preprocessor guards were needed.

## Issues Encountered

- `cmake --install` for Eigen failed when building BLAS static library. Workaround: installed headers and cmake config files manually (Eigen is header-only, so only headers and CMake config are needed).

## User Setup Required

None — no external service configuration required. CI runners will download and build Eigen versions from source.

## Next Phase Readiness

- CI-01 and CI-02 complete
- Eigen 5.x compatibility verified locally — zero code changes needed
- Ready for CI push to confirm GitHub Actions runs all 8 jobs green
- Phase complete, ready for next step in v1.2 milestone

---
*Phase: 15-eigen-5-x-ci*
*Completed: 2026-06-05*
