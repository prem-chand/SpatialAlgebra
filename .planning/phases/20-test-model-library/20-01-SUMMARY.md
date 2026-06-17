---
phase: 20-test-model-library
plan: 01
subsystem: testing
tags: [eigen3, cmake, interface-target, robotics, spatial-algebra, header-only]

# Dependency graph
requires: []
provides:
  - Zero-dependency test model library (INTERFACE target test_models, Eigen3 only)
  - RobotModel, JointSpec, JointType data types
  - RobotSolver abstract base class with 10 pure virtual methods
  - CMake infrastructure for test_models, sa_test_adapter, compile_smoke_test
affects: [20-test-model-library plans 02-04, 21-test-refinement, 22-pinocchio-cpp]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - INTERFACE CMake target for zero-dependency header-only library
    - PIMPL pattern for adapter header (scheduled for Plan 04)
    - Builder pattern via static build() factory (adapter-specific, per D-05)
    - Pure-data POD structs for kinematic chain description (per D-07/D-08)

key-files:
  created:
    - tests/test-models/robot_model.h
    - tests/test-models/robot_solver.h
    - tests/test-models/CMakeLists.txt
    - tests/test-models/sa_adapter.cpp
  modified:
    - CMakeLists.txt

key-decisions:
  - "INTERFACE target test_models links only Eigen3::Eigen — build-system-enforced compilation firewall"
  - "Derived getDOF() from joints.size() per RESEARCH.md Option B (1 DOF per joint)"
  - "Dynamic Eigen types (VectorXd/MatrixXd) for all joint-space quantities — no templates"
  - "build() factory NOT declared in RobotSolver base — each adapter provides its own with concrete return type"

patterns-established:
  - "Compilation firewall: INTERFACE CMake target prevents transitive SA inclusion accidentally"
  - "Adapter pattern: abstract RobotSolver + concrete adapter (sa_adapter separated from test_models)"

requirements-completed: [TML-01, TML-03, TML-05]

# Metrics
duration: 4min 23s
completed: 2026-06-17
---

# Phase 20 Plan 01: Test Model Library Data Types + CMake Infrastructure

**Zero-dependency robot_model.h (JointSpec, RobotModel, JointType), abstract RobotSolver interface (10 pure virtual methods), and CMake INTERFACE target with Eigen3-only compilation firewall**

## Performance

- **Duration:** 4 min 23 sec
- **Started:** 2026-06-17T03:20:55Z
- **Completed:** 2026-06-17T03:25:18Z
- **Tasks:** 3
- **Files modified:** 5 (4 created, 1 modified)

## Accomplishments

- Created `robot_model.h` — JointType enum (REVOLUTE, PRISMATIC, FIXED), JointSpec POD struct (8 fields), RobotModel struct with `getDOF()` — all in namespace `test_models` with zero SpatialAlgebra includes
- Created `robot_solver.h` — RobotSolver abstract base class with 10 pure virtual methods (setState, computeTorques, computeAccelerations, forwardKinematics, getJointTransform, computeMassMatrix, computeGravityTorques, computeJointSpaceJacobian, getLinkCOM, getDOF) — dynamic Eigen types, no templates
- Created `tests/test-models/CMakeLists.txt` — INTERFACE target `test_models` (Eigen3::Eigen only), `sa_test_adapter` library (test_models + SpatialAlgebra), `compile_smoke_test` executable (test_models only), install rules (TML-05)
- Updated root `CMakeLists.txt` — `add_subdirectory(tests/test-models)` placed after `enable_testing()`, before existing test targets; all 11 existing test executables unchanged
- Doxygen `@brief`/`@details` on all public declarations per project conventions
- `cmake -B build` configures successfully — compilation firewall verified at build-system level

## Task Commits

Each task was committed atomically:

1. **Task 1: Create robot_model.h** — `3554ac0` (feat)
2. **Task 2: Create robot_solver.h** — `8356599` (feat)
3. **Task 3: Create test-models CMakeLists.txt + update root CMakeLists.txt** — `68bd915` (feat)

## Files Created/Modified

- `tests/test-models/robot_model.h` — JointType enum, JointSpec struct (8 fields), RobotModel struct with getDOF()
- `tests/test-models/robot_solver.h` — RobotSolver abstract base class, 10 pure virtual methods
- `tests/test-models/CMakeLists.txt` — INTERFACE test_models, sa_test_adapter, compile_smoke_test, install rules
- `tests/test-models/sa_adapter.cpp` — Placeholder source for sa_test_adapter (implementation deferred to Plan 04)
- `CMakeLists.txt` — Added `add_subdirectory(tests/test-models)` after `enable_testing()`

## Decisions Made

None — followed plan as specified. All design decisions (D-01 through D-15) were locked in CONTEXT.md and RESEARCH.md.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Created placeholder sa_adapter.cpp for CMake configure-time source check**
- **Found during:** Task 3 (CMakeLists.txt creation)
- **Issue:** Plan assumed "CMake only evaluates the source list at build time, not configure time" — but `add_library(sa_test_adapter sa_adapter.cpp)` requires `sa_adapter.cpp` to exist at configure time. CMake configure (`cmake -B build`) failed with "Cannot find source file: sa_adapter.cpp".
- **Fix:** Created a minimal placeholder `tests/test-models/sa_adapter.cpp` with a comment-only content so CMake configure succeeds. The file serves as a compilation unit for the `sa_test_adapter` target; its implementation is deferred to Plan 04. The empty `.o` file produces an empty `.a` archive with no symbols (harmless — full implementation replaces it in Plan 04).
- **Files modified:** `tests/test-models/sa_adapter.cpp` (new)
- **Verification:** `cmake -B build` succeeds; `cmake --build build --target sa_test_adapter` produces `libsa_test_adapter.a` (empty archive, as expected)
- **Committed in:** `68bd915` (part of Task 3 commit)

---

**Total deviations:** 1 auto-fixed (Rule 3 blocking)
**Impact on plan:** Minimal — the placeholder file is intentionally empty and its existence is documented for Plan 04 to replace. No scope creep or architectural change.

## Issues Encountered

None beyond the CMake source-check issue documented above.

## Known Stubs

- `tests/test-models/sa_adapter.cpp` — Placeholder source file for `sa_test_adapter` CMake target. Contains only a comment. Full adapter implementation in Plan 04.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- `test_models` INTERFACE target is ready for Plans 02-03 to add `chains/` model factory headers
- `sa_test_adapter` target is declared and configures — ready for Plan 04 to add adapter implementation
- `compile_smoke_test` target is declared — ready for Plan 04 to update `tests/compile_smoke_test.cpp` with test_models-specific includes
- Existing test targets are unchanged and continue to build/pass normally

---

*Phase: 20-test-model-library*
*Completed: 2026-06-17*
