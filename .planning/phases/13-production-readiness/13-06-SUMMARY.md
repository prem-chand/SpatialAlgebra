---
phase: 13-production-readiness
plan: 06
subsystem: testing
tags: gravity, dynamics, conventions, rnea, aba, invariants

# Dependency graph
requires:
  - phase: 13-production-readiness
    plan: 02
    provides: gravity API in ForwardDynamics/InverseDynamics
  - phase: 13-production-readiness
    plan: 04
    provides: multi-link dynamics tests
provides:
  - MATHEMATICAL_CONVENTIONS.md formal specification
  - Gravity invariant test oracles for RNEA and ABA
  - Single-link direct round-trip consistency test
affects: [13-07, verification, future phases relying on gravity dynamics]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Independent structural invariants for gravity verification (not solver cross-validation)
    - Proportionality ratios as verification oracles

key-files:
  created:
    - MATHEMATICAL_CONVENTIONS.md
  modified:
    - tests/TestInverseDynamics.cpp
    - tests/TestForwardDynamics.cpp
    - tests/TestDynamicsConsistency.cpp

key-decisions:
  - "Gravity test oracles use independently verifiable invariants (proportionality, collinearity, mass-scaling) instead of solver round-trip self-consistency"
  - "Tolerance strategy documented: 1e-10 for unit tests, 1e-8 for multi-step dynamics"

patterns-established: []

requirements-completed: [UTL-03, ABA-01, ABA-02, TST-07]

# Metrics
duration: 8min
completed: 2026-05-17
---

# Phase 13 Plan 06: Mathematical Conventions & Gravity Oracles Summary

**Formal conventions document, 6 new independently verifiable gravity invariant tests across RNEA, ABA, and round-trip consistency**

## Performance

- **Duration:** 8 min
- **Started:** 2026-05-17T07:53:57Z
- **Completed:** 2026-05-17T08:02:41Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments
- Created `MATHEMATICAL_CONVENTIONS.md` (211 lines) formalizing cross-product formulas, gravity conventions (c₀ = -g for ABA, a₀ = S·q̈ − [0; g] for RNEA), Plücker transform conventions, tolerance strategy, and API change log
- Added 3 gravity invariant tests to TestInverseDynamics covering static equilibrium proportionality, zero-torque collinear gravity, and multi-chain validity
- Added 2 ABA gravity tests to TestForwardDynamics covering mass-scaling and gravity proportionality invariants
- Added single-link direct round-trip test to TestDynamicsConsistency establishing RNEA→ABA consistency gate

## Task Commits

Each task was committed atomically:

1. **Task 1: Create MATHEMATICAL_CONVENTIONS.md** - `db8b68a` (docs)
2. **Task 2: Add gravity invariant tests to TestInverseDynamics** - `0885baf` (test)
3. **Task 3: Add ABA gravity and round-trip tests** - `6799735` (test)

## Files Created/Modified
- `MATHEMATICAL_CONVENTIONS.md` - Formal conventions: cross products, gravity, Plücker transforms, tolerance, change log (211 lines)
- `tests/TestInverseDynamics.cpp` - Added 3 gravity invariant tests (8→11 tests, all pass)
- `tests/TestForwardDynamics.cpp` - Added 2 ABA gravity tests (9→11 tests, all pass)
- `tests/TestDynamicsConsistency.cpp` - Added 1 single-link round-trip test (6→7 tests, 4 pass + 3 pre-existing CR-02 failures)

## Decisions Made
- Gravity tests use independent structural invariants (proportionality, collinearity, mass-scaling) rather than solver cross-validation, directly addressing the HIGH review concern about "consistency-based validations"
- Mass-scaling test in ABA uses identity I_cm providing a constant inertia denominator, so qddot ratio cleanly reflects mass ratio (2:1)
- Gravity proportionality test verifies (qddot(g₁)−qddot(0))/(qddot(g₂)−qddot(0)) = g₁/g₂, invariant to mass/com

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed ratio comparison logic in SingleLinkStaticGravityProportionality**
- **Found during:** Task 2 (first test run)
- **Issue:** Initial loop attempted to compare first non-zero-g ratio against uninitialized sentinel value -1
- **Fix:** Added `ratio_initialized` boolean flag to gate consistency comparisons
- **Files modified:** tests/TestInverseDynamics.cpp
- **Verification:** Test passes at all 3 gravity levels
- **Committed in:** 0885baf (part of Task 2 commit)

**2. [Rule 2 - Missing Info] Corrected test formula for ABA gravity mass-scaling**
- **Found during:** Task 3 (first test run)
- **Issue:** Initial test used `lt::Identity(3) * 0.0` (zero I_cm) which produced zero denominator in ABI::apply, causing "Near-zero inertia" exception. ABI stores H and M separately without the Steiner term m·[c]×·[c]×ᵀ, requiring non-zero I_cm for joint axes with angular components.
- **Fix:** Used `lt::Identity(3)` as specified in plan. I_cm provides constant denominator so qddot ∝ m.
- **Files modified:** tests/TestForwardDynamics.cpp
- **Verification:** GravityEffectScalesWithMass test passes
- **Committed in:** 6799735 (part of Task 3 commit)

**3. [Rule 3 - Blocking] Changed joint axis from Z to X in GravityProportionalityInvariant**
- **Found during:** Task 3 (first test run)
- **Issue:** Z-axis joint with gravity parallel to joint axis produces zero gravitational torque. The gravity vector (0,0,-9.81) along Z cannot create torque about a Z-axis revolute joint at q=0 regardless of COM offset.
- **Fix:** Changed both tests to use X-axis revolute joints where gravity creates torque about the joint axis through COM=(0,1,0) offset
- **Files modified:** tests/TestForwardDynamics.cpp
- **Verification:** Both GravityEffectScalesWithMass and GravityProportionalityInvariant pass
- **Committed in:** 6799735 (part of Task 3 commit)

---

**Total deviations:** 3 auto-fixed (1 bug, 1 missing info, 1 blocking)
**Impact on plan:** All fixes correct physics errors in test design. No scope creep.

## Issues Encountered
- Plan verification section undercounted new tests (said 5, actually 6: 3 in TestInverseDynamics + 2 in TestForwardDynamics + 1 in TestDynamicsConsistency). Task spec correctly listed 2 for ForwardDynamics.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Mathematical conventions formally documented for all future development
- Gravity invariants verified independently for both RNEA and ABA
- Single-link round-trip gate established before multi-link consistency tests
- Ready for 13-07 (remaining production readiness items)

---

*Phase: 13-production-readiness*
*Completed: 2026-05-17*
