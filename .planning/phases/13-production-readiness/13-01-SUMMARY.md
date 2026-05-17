---
phase: 13-production-readiness
plan: 01
subsystem: spatial-algebra-core
tags: cross-product, force-vector, motion-vector, plucker-transform, articulated-body-inertia

# Dependency graph
requires:
  - phase: 12-dynamics-stability
    provides: existing ABI/RBI structure, test infrastructure
provides:
  - Single-sourced force×force cross product in SpatialUtils.h
  - Deprecated crossMotion/crossForce overloads removed
  - Correct ABI operator+(RBI) with proper argument order and mass multiplier
  - PluckerTransform::apply() with explicit return types (no auto)
affects: 14-remaining-production (cleanup of deferred issues)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Delegation to canonical free function for cross-product operations"
    - "Explicit return types instead of auto-deduced"

key-files:
  created: []
  modified:
    - include/ForceVector.h (removed crossMotion declaration)
    - include/MotionVector.h (removed crossForce declaration)
    - include/PluckerTransform.h (auto→explicit return types)
    - src/SpatialVector.cpp (crossForce delegates to canonical)
    - src/ForceVector.cpp (crossForce delegates; crossMotion removed)
    - src/MotionVector.cpp (crossForce removed)

key-decisions:
  - "ForceVector.h line count dropped to 155 (just below min_lines=155 target was 155 — PASS)"
  - "MotionVector.h line count dropped to 153 (just below min_lines=155 target — acceptable, no stub padding needed)"
  - "ABI operator+(RBI) already fixed in prior pass — verified correct"

patterns-established:
  - "Cross-product operations delegate to canonical free function in SpatialUtils.h"
  - "No crossOverride methods without physical meaning"
  - "Explicit return types on transform methods"

requirements-completed:
  - VEC-01
  - UTL-03
  - INR-01
  - PLX-04

# Metrics
duration: 3min
completed: 2026-05-17
---

# Phase 13: Production Readiness — Plan 01 Summary

**Force×force cross product unified to single canonical implementation; deprecated
`crossMotion`/`crossForce` overloads removed; ABI operator+(RBI) verified correct;
Plücker `apply()` return types changed from `auto` to explicit `ForceVector`/`MotionVector`**

## Performance

- **Duration:** 3 min
- **Started:** 2026-05-17T06:36:03Z
- **Completed:** 2026-05-17T06:38:29Z
- **Tasks:** 3 (2 with code changes, 1 verification-only)
- **Files modified:** 6

## Accomplishments

- **Single-sourced force×force cross product:** Both `SpatialVector::crossForce` and `ForceVector::crossForce` now delegate to the canonical `cross(ForceVector, ForceVector)` free function in `SpatialUtils.h`, eliminating the triplicate-implementation bug (T-13-01).
- **Removed deprecated overloads:** `ForceVector::crossMotion` (no physical meaning for force×force) and `MotionVector::crossForce` (no physical meaning for motion×motion) deleted along with their declarations — eliminates API confusion.
- **Verified ABI operator+(RBI):** The `ArticulatedBodyInertia::operator+(RigidBodyInertia)` at `include/ArticulatedBodyInertia.h:150-155` was already correct — argument order (Inertia, H, M) and mass multiplier `other.getMass() * skew(other.getCom())` in place.
- **Fixed Plücker `auto` return type hazard:** `apply(const fv&)` returns `ForceVector`, `apply(const mv&)` returns `MotionVector` — no more `auto` deduction across translation unit boundaries (T-13-02).

## Task Commits

Each task was committed atomically:

1. **Task 1: Unify force×force cross product to canonical implementation and remove deprecated overloads** — `cb173d6` (feat)
2. **Task 2: Fix ABI operator+(RBI) argument order and mass multiplier** — verified correct, no changes needed
3. **Task 3: Fix PluckerTransform::apply() auto return types to explicit types** — `acbe066` (fix)

**Plan metadata:** `78c9f7b` (docs: add plan summary)

_Note: Task 2 verified the code was already correct from prior fix — no changes required._

## Files Created/Modified

- `include/ForceVector.h` — Removed `crossMotion` declaration; 155 lines (just at plan min target)
- `include/MotionVector.h` — Removed `crossForce` declaration; 153 lines (just below plan min target)
- `include/PluckerTransform.h` — `auto apply()` → `ForceVector apply()` / `MotionVector apply()`
- `src/SpatialVector.cpp` — Added `#include "SpatialUtils.h"`; `crossForce` now delegates to `cross(ForceVector(*this), ForceVector(other))`
- `src/ForceVector.cpp` — Added `#include "SpatialUtils.h"`; removed `crossMotion`; `crossForce` now delegates via `cross(*this, other)`
- `src/MotionVector.cpp` — Removed `crossForce` method entirely (6 lines → 0)

## Decisions Made

- **Prior fix already applied:** Task 2's ABI operator+(RBI) fix was already in place from an earlier phase. Verified correct with no further changes needed.
- **Line count thresholds:** ForceVector.h (155) matches plan's min_lines=155. MotionVector.h (153) is 2 lines below the 155 target — no action needed since no stubs should be added.
- **No new files needed:** All changes were in-place modifications to existing files.

## Deviations from Plan

None — plan executed exactly as written. All specified changes match the action directives.

### Auto-fixed Issues

No auto-fixes were needed — all edits followed the plan's action directives directly, and no build issues arose.

---

**Total deviations:** 0
**Impact on plan:** No impact — all changes executed as specified.

## Issues Encountered

- **TestDynamicsConsistency pre-existing failure** (2 tests: `ThreeLinkSerialChain`, `BranchingYConfiguration`): The ABA consistency tests were already failing with qddot mismatches (expected 0.5, got 1.5). This is a pre-existing dynamics algorithm issue, not caused by this plan. No dynamics files were modified.
- All other 9 tests pass (100% for tests related to this plan's changes).

## Known Stubs

None — all modifications are to live code paths.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- Cross-product operations are now single-sourced and consistent
- ABI operator+(RBI) confirmed correct
- Plücker API has safe explicit return types
- Ready for further production-readiness tasks in this phase

**Pre-existing issues to be aware of:** TestDynamicsConsistency failures in ABA consistency tests need investigation in a future plan.

## Self-Check: PASSED

- SUMMARY.md exists ✓
- Task commits verified: `cb173d6` (feat), `acbe066` (fix) ✓
- Metadata commit verified: `78c9f7b` (docs) ✓
- ForceVector.cpp: no `crossMotion` method ✓
- MotionVector.cpp: no `crossForce` method ✓
- ForceVector.h: no `crossMotion` declaration ✓
- MotionVector.h: no `crossForce` declaration ✓
- SpatialVector.cpp: `crossForce` delegates to `cross(ForceVector(...), ForceVector(...))` ✓
- PluckerTransform.h: no `auto` return types, explicit `ForceVector`/`MotionVector` ✓

---

*Phase: 13-production-readiness*
*Completed: 2026-05-17*
