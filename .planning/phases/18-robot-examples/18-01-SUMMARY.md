---
phase: 18-robot-examples
plan: 01
subsystem: core
tags: [dynamics, forward-dynamics, inverse-dynamics, rnea, aba, ur5, examples]

requires:
  - phase: 09-integration-tests
    provides: dynamics consistency testing patterns
  - phase: 12-dynamics-consistency
    provides: forward/inverse dynamics solvers with cross-validation
  - phase: 13-production-readiness
    provides: NaN/Inf validation, solver hardening

provides:
  - Two working robot dynamics example executables (2-link Z-Z, 3-link Z-Y-Z)
  - UR5-derived link parameter model for serial chain dynamics
  - Documentation of two pre-existing solver limitations:
    * ABA bug: incorrect force propagation with non-zero COM
    * RNEA limitation: fixed transforms (q-independent)

affects: [library refactoring phases, ABA bug fix]

tech-stack:
  added: []
  patterns: [UR5-derived parameters, cross-validation pattern, solver finding documentation]

key-files:
  created:
    - examples/example_robot_2link.cpp
    - examples/example_robot_3link.cpp
  modified:
    - examples/CMakeLists.txt

key-decisions:
  - "3-link Z-Y-Z cross-validation FAILS due to pre-existing ABA bug — example documents finding instead of passing"
  - "RNEA solver uses fixed transforms (q-independent) — results valid only at home configuration"
  - "2-link Z-Z cross-validation PASSES but only because tau_g=0 (gravity parallel to both axes)"
  - "UR5-derived parameters (mass 8.393 kg, COM 0.2125 m, inertia 0.3564 kg·m^2) for realistic physics"

patterns-established:
  - "Robot example structure: model setup → FD tests → ID gravity tests → cross-validation → summary"
  - "UR5-derived anthropomorphic arm parameters for Z-Y-Z serial chain"
  - "Cross-validation used as solver consistency check"

requirements-completed: [EX-01, EX-02]

duration: 75min
completed: 2026-06-06
---

# Phase 18: Robot Examples Summary

**Two robot dynamics executables (2-link Z-Z and 3-link Z-Y-Z arms) with UR5-derived parameters, cross-validation revealing two pre-existing solver limitations: ABA multi-link COM propagation bug and RNEA fixed-transform limitation**

## Performance

- **Duration:** 75 min
- **Started:** 2026-06-06T11:30:00Z
- **Completed:** 2026-06-06T12:45:00Z
- **Tasks:** 3 (split into 4 commits)
- **Files modified:** 3

## Accomplishments
- Created 2-link Z-Z planar arm example with FD, ID gravity, and cross-validation (all passing)
- Created 3-link Z-Y-Z anthropomorphic arm example with UR5-derived parameters for realistic dynamics
- Registered both executables in CMakeLists.txt with library and include dependencies
- Discovered and documented two pre-existing solver limitations:
  1. **ABA bug (ForwardDynamics.cpp):** Incorrect force propagation for multi-link chains with non-zero COM offsets — ID→FD round-trip fails when both COM offset and multiple links are present
  2. **RNEA limitation (InverseDynamics.cpp):** Fixed transform X does not update with joint position q — results valid only at home configuration

## Task Commits

Each task was committed atomically:

1. **Task 1: Create 2-link Z-Z arm example** — `da21d47` (feat)
2. **Task 2: Create 3-link Z-Y-Z arm example** — `735ca30` (feat)
3. **Task 3a: Register examples in CMakeLists.txt** — `a9186c5` (chore)
4. **Task 3b: Update 3-link with accurate solver findings** — `59e0c92` (docs)

**Plan metadata:** (final summary commit to follow — not yet part of plan)

## Files Created/Modified
- `examples/example_robot_2link.cpp` — 2-link Z-Z planar arm dynamics example with UR5 shoulder parameters
- `examples/example_robot_3link.cpp` — 3-link Z-Y-Z spatial arm dynamics example with full UR5 arm parameters
- `examples/CMakeLists.txt` — Added `example_robot_2link` and `example_robot_3link` executables

## Decisions Made

1. **Document solver findings rather than hiding them** — The 3-link cross-validation FAILS due to a pre-existing ABA bug. Rather than modifying library code (prohibited by scope), the example honestly documents the bug with numerical evidence. This is more valuable than a synthetic pass.
2. **UR5-derived parameters for physical realism** — Using actual UR5 masses (3.7 kg, 8.393 kg, 2.33 kg) and COM positions makes the examples useful for practitioners comparing with real robot data.
3. **Separate FD/ID model setup** — ForwardDynamics uses `Link` structs, InverseDynamics uses `InverseDynamicsLink` structs. Both are populated from the same physical parameters but registered in parallel vectors.
4. **Cross-validation as primary consistency check** — FD(ID(0,g),g) ≈ 0 is the gold standard for solver correctness. The 2-link example passes (degenerate case, tau_g=0). The 3-link fails, revealing the ABA bug.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Documented pre-existing library bugs found during cross-validation**
- **Found during:** Task 2 / Task 3b (3-link example cross-validation)
- **Issue:** Plan claimed cross-validation would pass (|qddot| < 1e-10) for the 3-link Z-Y-Z arm. It fails with qddot_J2 = 69.63 rad/s^2 instead of 0, revealing two pre-existing library limitations: ABA bug with non-zero COM multi-link propagation, and RNEA fixed-transform limitation.
- **Fix:** Updated the 3-link example to honestly document the findings — added a LIBRARY FINDINGS section explaining both issues, used accurate torque values, and removed incorrect physical interpretation commentary.
- **Files modified:** examples/example_robot_3link.cpp
- **Verification:** Build succeeds, all 11 tests pass, output correctly documents both issues
- **Committed in:** 59e0c92 (docs commit)

**2. [Rule 3 - Blocking] Fixed gravity torque interpretation — RNEA uses fixed transforms**
- **Found during:** Task 3b (writing gravity test commentary)
- **Issue:** Original 45° test claimed torque was "reduced from horizontal case" but the RNEA solver produces identical results at all q positions because transforms X are fixed at model setup time.
- **Fix:** Added explicit note about the RNEA limitation. Gravity test 2 now shows the same torque and explains why.
- **Files modified:** examples/example_robot_3link.cpp
- **Verification:** Compiled output confirms both tests show -13.27 N·m and clearly documents the limitation.
- **Committed in:** 59e0c92 (docs commit)

---

**Total deviations:** 2 auto-fixed (1 missing critical, 1 blocking)
**Impact on plan:** Both deviations improve example quality and honesty. The ABA bug finding is a significant discovery that should inform future refactoring phases.

## Issues Encountered

1. **ABA solver round-trip inconsistency:** ID→FD round-trip fails for multi-link chains with non-zero COM. Extensive diagnostic testing confirmed:
   - Works for single joints with zero COM ✓
   - Fails for single Y-axis joint with COM offset under gravity ✗
   - Fails for 2-link Z-Z with non-zero COM even without gravity ✗
   - Root cause: ABA pass has incorrect force/acceleration propagation
2. **RNEA fixed-transform limitation:** Transforms X are never recomputed as a function of joint position q. All results are valid only at home configuration (q=0 for all joints).
3. **Physical interpretation accuracy:** Initial hand-calculated gravity torque (~30.7 N·m) differed from solver output (-13.27 N·m) due to COM z-offsets and sign convention. Corrected with actual solver values.

## Known Stubs

None — both examples are complete and produce meaningful output.

## Threat Flags

None — no new security-relevant surface introduced.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Two robot example executables ready for demonstration and testing
- Significant solver bugs discovered:
  - ABA needs force propagation fix for non-zero COM multi-link chains
  - RNEA needs transform update mechanism for non-zero q positions
- Both issues should be addressed in a library refactoring phase before more complex robot examples (e.g., 6-DOF UR5, floating base)
- The 2-link Z-Z example works fully and is suitable for CI demonstration

---

*Phase: 18-robot-examples*
*Completed: 2026-06-06*

## Self-Check: PASSED

All created/modified files found:
- `examples/example_robot_2link.cpp` ✓
- `examples/example_robot_3link.cpp` ✓
- `examples/CMakeLists.txt` ✓
- `.planning/phases/18-robot-examples/18-01-SUMMARY.md` ✓

All commits found:
- `da21d47` ✓ (2-link example)
- `735ca30` ✓ (3-link example)
- `a9186c5` ✓ (CMakeLists.txt)
- `59e0c92` ✓ (3-link update)

Built executables found:
- `build/examples/example_robot_2link` ✓
- `build/examples/example_robot_3link` ✓

Tests (11/11): PASS ✓
