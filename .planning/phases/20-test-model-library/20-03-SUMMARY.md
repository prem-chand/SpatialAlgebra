---
phase: 20-test-model-library
plan: 03
subsystem: test-models/chains
tags: [test-models, factory-functions, inverse-dynamics, forward-dynamics, consistency, spatial-operations]
requires:
  - 20-02 (simple-domain chains)
provides:
  - tests/test-models/chains/inverse_dynamics.h
  - tests/test-models/chains/forward_dynamics.h
  - tests/test-models/chains/consistency.h
  - tests/test-models/chains/spatial_operations.h
affects:
  - all 11 test domains now have chain model representation
tech-stack:
  added: []
  patterns: [inline-factory-functions, Eigen3-only, pure-data-pod]
key-files:
  created:
    - tests/test-models/chains/inverse_dynamics.h
    - tests/test-models/chains/forward_dynamics.h
    - tests/test-models/chains/consistency.h
    - tests/test-models/chains/spatial_operations.h
  modified: []
decisions:
  - Branching Y uses parent=0 for both children (not serial i-1) — matches existing test models
  - Consistency chains use distinct name prefixes (rt_, nzc_, y_) but identical mass/COM/inertia values
  - makeSingleLinkXAxisWithCOMOffset uses jointAxis=UnitX() — only non-Z joint in inverse_dynamics.h
  - makeTwoLinkWith90DegreeRotation uses AngleAxisd(M_PI_2, UnitZ()) for rotation assembly
metrics:
  duration: ""
  completed_date: 2026-06-17
---

# Phase 20 Plan 03: Dynamics-Domain Chain Model Factories Summary

**One-liner:** Created 15 dynamics-domain factory functions across 4 headers (inverse_dynamics, forward_dynamics, consistency, spatial_operations) using only Eigen3 and robot_model.h, completing TML-02 coverage for all 11 test domains.

## Tasks Completed

| Task | Name | Commits | Files |
|------|------|---------|-------|
| 1 | Create inverse_dynamics.h and forward_dynamics.h | 19f84bc | `tests/test-models/chains/inverse_dynamics.h`, `tests/test-models/chains/forward_dynamics.h` |
| 2 | Create consistency.h and spatial_operations.h | e6135e4 | `tests/test-models/chains/consistency.h`, `tests/test-models/chains/spatial_operations.h` |

## Factories Created

### inverse_dynamics.h (5 functions)
- `makeSingleLinkChain()` — 1-link Z-revolute, baseline RNEA model
- `makeTwoLinkSerialChain()` — 2-link serial, 1m X spacing
- `makeBranchingYConfiguration()` — 3-link Y-tree, both children parent=0
- `makeTwoLinkWithCOMOffset()` — 2-link with COM=[0,0.1,0] on each
- `makeSingleLinkXAxisWithCOMOffset()` — X-axis revolute, COM=[0,0.5,0]

### forward_dynamics.h (3 functions)
- `makeThreeLinkSerialChain()` — 3-link serial, 1m X spacing
- `makeThreeLinkSerialChainNonZeroCOM()` — 3-link with COM=[0.1,0,0] on all links
- `makeTwoLinkWith90DegreeRotation()` — 2-link with 90° Z rotation on link 1

### consistency.h (4 functions)
- `makeSingleLinkConsistencyChain()` — same as makeSingleLinkChain, prefix "roundtrip_"
- `makeThreeLinkConsistencyChain()` — same as makeThreeLinkSerialChain, prefix "rt_"
- `makeThreeLinkNonZeroCOMConsistencyChain()` — same as makeThreeLinkSerialChainNonZeroCOM, prefix "nzc_"
- `makeBranchingYConsistencyChain()` — same as makeBranchingYConfiguration, prefix "y_"

### spatial_operations.h (3 functions)
- `makeTransformRBIPair()` — 45° X rotation+translation, mass=2, COM=[0.5,0,0]
- `makeCrossProductMotionFixtures()` — Z+X revolute joints for motion cross product
- `makeCrossProductForceFixtures()` — Z+Y revolute joints with offset COM for force cross product

## Verification Results

| Check | Result |
|-------|--------|
| 11 headers in chains/ | ✅ PASS (7 simple + 4 dynamics) |
| Zero SpatialAlgebra includes | ✅ PASS across all 11 headers |
| namespace test_models | ✅ PASS in all 11 headers |
| #pragma once | ✅ PASS in all 4 new headers |
| Consistency parameter equivalence | ✅ PASS (same mass, COM, inertia values) |
| Branching Y parent=0 | ✅ PASS (both inverse_dynamics.h and consistency.h) |
| X-axis joint in makeSingleLinkXAxisWithCOMOffset | ✅ PASS (jointAxis=UnitX) |
| makeTwoLinkWith90DegreeRotation assembly | ✅ PASS (AngleAxisd+M_PI_2+UnitZ) |
| Doxygen on all functions | ✅ PASS (all 15 functions) |

## Deviations from Plan

None — plan executed exactly as written. All 15 factory functions match the specifications in RESEARCH.md Domains 8-11.

## Known Stubs

None. All factory functions return fully-specified RobotModel instances with concrete mass, COM, inertia, and transform values — no placeholder data or empty fields.

## Threat Flags

None. All files are pure data factory functions that construct RobotModel instances from compile-time constants. No network endpoints, auth paths, file access patterns, or schema changes at trust boundaries beyond what is documented in the plan's `<threat_model>`.

## TML-02 Requirement Status

TML-02 requires "kinematic chain definitions covering all 11 test domains." With this plan:
- **7 simple domains** (Plan 02): spatial_vectors, plucker_transforms, rotation, lower_triangular, rigid_body_inertia, articulated_body, spatial_utils
- **4 dynamics domains** (Plan 03): inverse_dynamics, forward_dynamics, consistency, spatial_operations
- **Total: 11/11 domains covered** ✅ TML-02 satisfied

## Self-Check: PASSED

- [x] tests/test-models/chains/inverse_dynamics.h exists
- [x] tests/test-models/chains/forward_dynamics.h exists
- [x] tests/test-models/chains/consistency.h exists
- [x] tests/test-models/chains/spatial_operations.h exists
- [x] Commit 19f84bc exists (Task 1)
- [x] Commit e6135e4 exists (Task 2)
