---
phase: 20-test-model-library
plan: 04
subsystem: tests/test-models
tags: [spatial-algebra, adapter, pimpl, compilation-firewall, rnea, aba, jacobian]

# Dependency graph
requires:
  - 20-01 (robot_model.h, robot_solver.h, CMakeLists.txt infrastructure)
  - 20-02 (simple chains: spatial_vectors, plucker_transforms, rotation, etc.)
  - 20-03 (dynamics chains: inverse_dynamics, forward_dynamics, consistency)
provides:
  - sa_adapter.h: SpatialAlgebraAdapter PIMPL declaration (zero SA types visible)
  - sa_adapter.cpp: Full 10-method adapter implementation bridging test-models ↔ SpatialAlgebra
  - compile_smoke_test.cpp: Updated zero-dependency verification smoke test
affects:
  - 21-test-refinement (existing tests can now use test_models + sa_test_adapter)
  - 22-pinocchio-cpp (Pinocchio adapter follows same pattern)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - PIMPL (pointer-to-implementation) — struct Impl hides all SA types
    - Dirty-flag pattern for lazy forward kinematics (fkDirty_ + autoFK)
    - RNEA-column method for computeMassMatrix (unit vector per DOF)
    - Column-by-column geometric Jacobian with ancestor checking
    - Symmetrized inertia conversion (Pitfall 2 guard)
    - Joint axis → MotionVector screw axis mapping (REVOLUTE/PRISMATIC/FIXED)
key-files:
  created:
    - tests/test-models/sa_adapter.h
  modified:
    - tests/test-models/sa_adapter.cpp
    - tests/compile_smoke_test.cpp
decisions:
  - "ForwardKinematics auto-triggers via dirty flag: setState marks dirty, first query calls autoFK"
  - "Mass matrix uses RNEA-column method (O(n²)) — acceptable for n ≤ 10 test models"
  - "Jacobian uses column-by-column geometric method with world-frame screw axis transformation"
  - "Joint metadata (types, axes) stored in Impl alongside solvers for FK reconstruction"
  - "build() factory constructs BOTH FD and ID solvers from single RobotModel"

patterns-established:
  - "PIMPL enforcement: unique_ptr<Impl> with destructor in .cpp — no SA types in header"
  - "Input validation at all entry points: size checks (T-20-01), bounds checks (T-20-02)"
  - "Compilation firewall: compile_smoke_test links ONLY test_models (INTERFACE), zero SA symbols"

requirements-completed: [TML-04]

# Metrics
duration: 7 min 23 sec
completed: 2026-06-17
---

# Phase 20 Plan 04: SpatialAlgebra Adapter Implementation + Compilation Firewall Verification

**PIMPL SpatialAlgebraAdapter bridging test-models to SA solvers, with all 10 RobotSolver methods, RNEA-column mass matrix, geometric Jacobian, and zero-dependency compilation firewall verified via smoke test**

## Performance

- **Duration:** 7 min 23 sec
- **Started:** 2026-06-17
- **Completed:** 2026-06-17
- **Tasks:** 3
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments

- Created `sa_adapter.h` — SpatialAlgebraAdapter class declaration with PIMPL pattern (`struct Impl` + `std::unique_ptr<Impl>`), extending RobotSolver with all 10 `override` methods, `build()` static factory, `fkDirty_` dirty flag, and `autoFK()` private helper. Zero SpatialAlgebra includes — only `robot_solver.h`, `robot_model.h`, `<memory>`, and `<Eigen/Dense>`.

- Implemented `sa_adapter.cpp` — Full 512-line adapter implementation:
  - **Impl struct:** `ForwardDynamics fd_`, `InverseDynamics id_`, `jointTypes_`, `jointAxes_`, `worldX_` transform cache
  - **`build()` factory:** Constructs both FD and ID solvers from the same RobotModel using file-scope conversion helpers
  - **`setState()`:** Input validation with `std::invalid_argument` on size mismatch (T-20-01), copies q/qdot to both solvers
  - **`computeTorques()`:** RNEA delegation with size check and auto forward kinematics
  - **`computeAccelerations()`:** ABA delegation, extracts qddot from FD link structs
  - **`forwardKinematics()`:** Recomputes world-frame transforms from joint angles — revolute creates `Rotation(AngleAxisd(q, axis))`, prismatic applies `axis * q` translation, fixed is identity
  - **`getJointTransform()`:** Bounds check (T-20-02), converts cached `PluckerTransform` to 4×4 homogeneous matrix
  - **`computeMassMatrix()`:** RNEA-column method — saves/restores qddot, sets unit vectors per DOF, extracts torque columns
  - **`computeGravityTorques()`:** Zero-acceleration RNEA with given gravity vector
  - **`computeJointSpaceJacobian()`:** Column-by-column geometric Jacobian — transforms joint screw axis to world frame, ancestor checking via parent walk
  - **`getLinkCOM()`:** Transforms local COM to world frame via `R * com + t`
  - **`autoFK()`:** Dirty-flag pattern: if `fkDirty_`, calls `forwardKinematics()` and clears flag
  - JointSpec→Link conversion: Matrix4d→PluckerTransform, **symmetrized inertia** (`0.5 * (M + M^T)`) before `LowerTriangular::fromFullMatrix()` (Pitfall 2), joint axis → MotionVector screw axis (Pitfall 3 convention)

- Updated `compile_smoke_test.cpp` — Replaced all SA-dependent code with test-models-only includes: `robot_model.h`, `robot_solver.h`, `chains/inverse_dynamics.h`. Instantiates `JointType`, `JointSpec`, `RobotModel`, calls `makeSingleLinkChain()` and `makeTwoLinkSerialChain()` factory functions. Links ONLY against `test_models` (INTERFACE, Eigen3 only — zero SpatialAlgebra symbols).

## Task Commits

Each task was committed atomically:

1. **Task 1: Create sa_adapter.h** — `bf7e4af` (feat)
2. **Task 2: Create sa_adapter.cpp** — `71254e2` (feat)
3. **Task 3: Update compile_smoke_test.cpp** — `a7f018d` (feat)

## Files Created/Modified

- `tests/test-models/sa_adapter.h` — PIMPL SpatialAlgebraAdapter declaration, 10 override methods, build() factory, fkDirty_ + autoFK()
- `tests/test-models/sa_adapter.cpp` — Full adapter implementation: Impl struct, conversion helpers, all 10 method bodies, input validation, symmetrized inertia, RNEA-column mass matrix, geometric Jacobian
- `tests/compile_smoke_test.cpp` — Updated to include robot_model.h + robot_solver.h + chains/ only, zero SpatialAlgebra symbols

## Deviations from Plan

None — plan executed exactly as written. All conversion patterns, method implementations, and threat mitigations followed the plan specifications precisely.

## Threat Mitigation Status

All 5 mitigated threats implemented:

| Threat | Mitigation | Status |
|--------|-----------|--------|
| T-20-01 (DoS: size mismatch) | `std::invalid_argument` in setState, computeTorques, computeAccelerations | ✅ Implemented |
| T-20-02 (Tampering: bounds) | `std::out_of_range` in getJointTransform, getLinkCOM, computeJointSpaceJacobian | ✅ Implemented |
| T-20-03 (Info Disclosure: NaN/Inf) | Accept — no sanitization | ✅ Accept |
| T-20-04 (DoS: large model) | Accept — std::vector::reserve may throw bad_alloc | ✅ Accept |
| T-20-05 (Tampering: PIMPL lifecycle) | `std::unique_ptr<Impl>` prevents double-free; build() returns unique_ptr | ✅ Implemented |
| T-20-06 (Tampering: compilation firewall) | compile_smoke_test links ONLY test_models (INTERFACE, zero SA) | ✅ Verified |

## Verification Results

### Build Verification
- ✅ `cmake --build build --target sa_test_adapter` — compiles without errors
- ✅ `cmake --build build --target compile_smoke_test` — compiles without errors
- ✅ `cmake --build build` (full build) — all targets succeed

### Runtime Verification
- ✅ `./build/tests/test-models/compile_smoke_test` — exits 0, prints "PASSED"
- ✅ `cd build && ctest --output-on-failure` — **all 12 tests pass, 0 failures**
  - TestSpatialVector, TestPluckerTransform, TestRotation, TestLowerTriangular
  - TestSpatialUtils, TestRigidBodyInertia, TestArticulatedBodyInertia
  - TestForwardDynamics, TestSpatialOperations, TestInverseDynamics
  - TestDynamicsConsistency, compile_smoke_test

### Compilation Firewall Verification
- ✅ `grep "#include.*SpatialAlgebra.h" tests/compile_smoke_test.cpp` — returns 0 (only in comments)
- ✅ `grep "^#include" tests/test-models/sa_adapter.h` — only robot_solver.h, robot_model.h, Eigen/Dense, memory
- ✅ `compile_smoke_test` links ONLY `test_models` INTERFACE target (Eigen3::Eigen only, zero SA symbols)

## Issues Encountered

None — the adapter compiled and all tests passed on the first build attempt.

## Known Stubs

None — all methods are fully implemented with complete input validation, exception handling, and computation logic.

## User Setup Required

None — no external service configuration required.

---

*Phase: 20-test-model-library*
*Completed: 2026-06-17*
