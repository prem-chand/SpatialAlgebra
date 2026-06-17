---
phase: 20-test-model-library
plan: 02
subsystem: test-models/chains
tags: [test-models, factory-functions, header-only, simple-domains]
dependency-graph:
  requires:
    - robot_model.h (JointSpec, RobotModel)
    - RobotModel, JointSpec types
  provides:
    - chains/spatial_vectors.h (1 factory)
    - chains/plucker_transforms.h (4 factories)
    - chains/rotation.h (3 factories)
    - chains/lower_triangular.h (3 factories)
    - chains/rigid_body_inertia.h (3 factories)
    - chains/articulated_body.h (2 factories)
    - chains/spatial_utils.h (3 factories)
  affects:
    - Plan 03 (inverse_dynamics, forward_dynamics, consistency, spatial_operations chain models)
    - Plan 04 (SA adapter)
tech-stack:
  added: []
  patterns:
    - header-only inline factory functions
    - Eigen3 + robot_model.h only (zero SpatialAlgebra)
    - Doxygen @brief/@details/@return on all functions
    - namespace test_models
    - RobotModel return-by-value with field-initializer pattern
key-files:
  created:
    - tests/test-models/chains/spatial_vectors.h
    - tests/test-models/chains/plucker_transforms.h
    - tests/test-models/chains/rotation.h
    - tests/test-models/chains/lower_triangular.h
    - tests/test-models/chains/rigid_body_inertia.h
    - tests/test-models/chains/articulated_body.h
    - tests/test-models/chains/spatial_utils.h
  modified: []
decisions:
  - D-11 applied: one header per domain under chains/
  - D-13 applied: header-only (all functions inline), no .cpp files
  - Factory naming follows "make" prefix convention with descriptive suffixes
  - multi-link chains use serial parent indexing (i-1 for link i)
metrics:
  duration: "\u22485 min"
  completed-date: 2026-06-17
---

# Phase 20 Plan 02: Simple-Domain Test Model Factory Functions Summary

7 header-only factory headers providing 19 `RobotModel`-returning functions for the simple test domains (1–7), with zero SpatialAlgebra dependency.

## Task-by-Task Results

### Task 1: spatial_vectors, plucker_transforms, rotation chain models

**Commit:** `0eae5a3`
**Files:** 3 files created, 403 lines

| Header | Functions | Count |
|--------|-----------|-------|
| `chains/spatial_vectors.h` | `makeSingleLink()` — 1-link Z-revolute chain (Featherstone Ex 2.1) | 1 |
| `chains/plucker_transforms.h` | `makeIdentityTransform()`, `make90DegreeZRotation()`, `makePureTranslation()`, `makeCombinedTransform()` — 4 transform variants | 4 |
| `chains/rotation.h` | `makeIdentityRotation()`, `make45DegreeXRotation()`, `make90DegreeZRotationFixture()` — 3 rotation variants | 3 |

All transforms use `Eigen::AngleAxisd` for rotation construction and `Matrix4d` for homogeneous assembly. The "Fixture" suffix on `make90DegreeZRotationFixture()` distinguishes rotation-domain tests from the transform-domain `make90DegreeZRotation()` in `plucker_transforms.h`.

### Task 2: lower_triangular, rigid_body_inertia, articulated_body, spatial_utils chain models

**Commit:** `58bd708`
**Files:** 4 files created, 572 lines

| Header | Functions | Count |
|--------|-----------|-------|
| `chains/lower_triangular.h` | `makeIdentityLT()`, `makeDiagonalLT(value=2.0)`, `makeArbitraryLT()` — matrix [[1,0,0],[2,3,0],[4,5,6]] | 3 |
| `chains/rigid_body_inertia.h` | `makeRBIDefault()`, `makeRBIOffsetCOM(comX,comY,comZ)`, `makeRBIDiagonalInertia(mass,ix,iy,iz)` | 3 |
| `chains/articulated_body.h` | `makeABIIdentity()`, `makeABIReduced()` — ABI→RBI reduction test models | 2 |
| `chains/spatial_utils.h` | `makeSkewFixture()`, `makeMotionVectorPair()` (2-link, COM=[1,0,0]), `makeForceVectorPair()` (2-link, mass=2, COM=[0,1,0]) | 3 |

Multi-link chains (`makeMotionVectorPair`, `makeForceVectorPair`) use correct serial parent indexing (i-1) and inter-link translation via `topRightCorner` of the 4×4 homogeneous transform.

## Verification Results

| Check | Result |
|-------|--------|
| 7 headers in `chains/` | ✅ PASS |
| Zero "SpatialAlgebra" references | ✅ PASS (0 matches) |
| All use `#pragma once` | ✅ PASS (7/7) |
| All in `namespace test_models` | ✅ PASS (7/7) |
| Doxygen `@brief` on every function | ✅ PASS (19 functions + 7 file headers) |
| All functions return `RobotModel` by value | ✅ PASS (19 `inline RobotModel` declarations) |
| Includes: only `robot_model.h` and `Eigen/Dense` | ✅ PASS |

**Total factory functions:** 1 + 4 + 3 + 3 + 3 + 2 + 3 = **19**

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Removed "SpatialAlgebra" mention from spatial_vectors.h Doxygen comment**
- **Found during:** Task 1 verification
- **Issue:** The file-level Doxygen comment in `spatial_vectors.h` mentioned "SpatialAlgebra, Pinocchio, or RBDL adapters", violating the zero-SA-references constraint.
- **Fix:** Replaced with "any solver adapter" — retains the cross-solver intent without naming specific implementations.
- **Files modified:** `tests/test-models/chains/spatial_vectors.h`
- **Commit:** Fixed before initial commit (`0eae5a3`)

## Known Stubs

None. All 19 factory functions return fully-initialized `RobotModel` instances with concrete Eigen3 values — no placeholder data, no TODO/FIXME markers, no empty/default stubs.

## Threat Flags

None. All threats in the plan's `<threat_model>` are disposition `accept` — no new security surface introduced beyond what was documented. The 7 headers are pure compile-time artifact generators with no runtime I/O, network access, or external dependency beyond Eigen3.

## Self-Check

- [x] `tests/test-models/chains/spatial_vectors.h` exists
- [x] `tests/test-models/chains/plucker_transforms.h` exists
- [x] `tests/test-models/chains/rotation.h` exists
- [x] `tests/test-models/chains/lower_triangular.h` exists
- [x] `tests/test-models/chains/rigid_body_inertia.h` exists
- [x] `tests/test-models/chains/articulated_body.h` exists
- [x] `tests/test-models/chains/spatial_utils.h` exists
- [x] Commit `0eae5a3` exists (Task 1)
- [x] Commit `58bd708` exists (Task 2)
- [x] Zero "SpatialAlgebra" in any chain header
- [x] 19 `inline RobotModel` function declarations total

**Self-Check: PASSED**
