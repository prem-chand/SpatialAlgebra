---
phase: 14-cr-02-bug-fix
plan: 01
status: complete
execution_date: 2026-06-04
test_status: all_pass
---

# Phase 14: CR-02 Bug Fix — Summary

**Plan:** 01 — ABA cross product fixes (root cause)

## One-Liner

Fixed cross product formulas in SpatialUtils.h (`cross(MotionVector, ForceVector)` and `cross(ForceVector, ForceVector)`) that corrupted bias forces for 3+ link chains, resolving all multi-link consistency test failures without restructuring the ABA inward pass.

## What Actually Happened

The phase plan proposed restructuring the ABA inward pass to a single tip-to-base sweep. However, an LLM agent independently diagnosed and fixed the actual root cause:

1. **BUG-01 — `cross(MotionVector, ForceVector)`: crf formula wrong** (SpatialUtils.h)
   - The `v₁×f₂` term was placed in the linear component instead of angular
   - Per Featherstone §2.9: `crf = [ω×τ + v×f; ω×f]`
   - This corrupted the bias force `pa` for intermediate links, causing incorrect qddot corrections in Phase 3 of `inwardPass()`

2. **BUG-02 — `cross(ForceVector, ForceVector)`: missing anti-commutativity term** (SpatialUtils.h)
   - Missing `- τ₂×f₁` in the linear component
   - This broke the anti-commutativity property `cross(a,b) = -cross(b,a)` for force×force products

3. **ISSUE-05 — InverseDynamics::inwardPass O(n²) → O(n)** (src/InverseDynamics.cpp)
   - Tip-to-base sweep replaced child-scanning loop, exploiting topological ordering invariant

4. **BUG-03 — Virtual destructor on SpatialVector** (include/SpatialVector.h)

5. **Namespace cleanup** — Removed global `using` declarations from headers

The ABA `inwardPass()` in `ForwardDynamics.cpp` was **not modified** — the 3-phase approach with partial qddot + correction is correct when the cross product formulas are right.

## Test Results

| Test Executable | Tests | Status |
|----------------|-------|--------|
| TestSpatialVector | all | ✅ PASS |
| TestPluckerTransform | all | ✅ PASS |
| TestRotation | all | ✅ PASS |
| TestLowerTriangular | all | ✅ PASS |
| TestSpatialUtils | all | ✅ PASS |
| TestRigidBodyInertia | all | ✅ PASS |
| TestArticulatedBodyInertia | all | ✅ PASS |
| TestForwardDynamics | 15 | ✅ PASS |
| TestSpatialOperations | all | ✅ PASS |
| TestInverseDynamics | all | ✅ PASS |
| TestDynamicsConsistency | 7 | ✅ PASS |
| **ctest (total)** | **11/11** | **✅ 100% PASS** |

## Key Tests Now Passing

- `ThreeLinkSerialChain` — ABA(RNEA(qddot)) ≈ qddot for 3-link serial chain
- `BranchingYConfiguration` — ABA(RNEA(qddot)) ≈ qddot for branching Y-configuration
- `TwoLinkRoundTripWithGravity` — 2-link chain with gravity round-trip
- `ThreeLinkNumericalValidation` — Numerical validation of 3-link qddot ordering
- `CondensationReducesInertiaNorm` — Condensed Ia norm < uncondensed Ia norm
- `ThreeLinkSingleTorque` — tau=[1,0,0] base joint accelerates positively

## Revision to Phase Plan

The original plan's core premise (ABA inward pass needs restructuring) was incorrect. The plan should not be executed — it would undo valid working code. The actual fix was entirely in cross product formulas.

## Files Changed

- `include/SpatialUtils.h` — crf formula + force×force anti-commutativity fix
- `include/SpatialVector.h` — virtual destructor
- `src/SpatialVector.cpp` — crossForce fix
- `src/InverseDynamics.cpp` — O(n²)→O(n) tip-to-base sweep
- `include/LowerTriangular.h` — remove global `using` 
- `include/Rotation.h` — remove global `using`
- `tests/TestLowerTriangular.cpp` — add `using namespace SpatialAlgebra`
- `tests/TestRotation.cpp` — add `using namespace SpatialAlgebra`
- `tests/TestSpatialOperations.cpp` — update expected crf values

## Phase Artefact Status

- [x] CONTEXT.md — gathered (2026-05-30)
- [x] RESEARCH.md — researched (2026-05-30)
- [x] PLAN.md — planned (2026-05-30, superseded by external fixes)
- [x] VALIDATION.md — created (2026-06-04)
- [x] VERIFICATION.md — plan checker analysis (2026-06-04)
- [x] SUMMARY.md — this file

**Phase goal met:** All multi-link consistency tests pass. All 11 test executables pass.
