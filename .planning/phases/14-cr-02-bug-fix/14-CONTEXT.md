# Phase 14: CR-02 Bug Fix - Context

**Gathered:** 2026-06-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Correct forward dynamics for chains with 3+ joints by implementing the missing articulation (condensation) step from Featherstone Algorithm 7.3 in the ABA inward pass. All multi-link consistency tests must pass with non-zero COM, and all 11 test executables must pass with no regressions.

**Confirmed by v1.2 (Phase 18):** Robot examples with UR5-derived parameters demonstrated that the ABA bug manifests as ID→FD round-trip failure specifically for multi-link chains with non-zero COM. The 2-link Z-Z arm passes cross-validation, the 3-link Z-Y-Z arm fails — isolating the issue to chains with 3+ joints.

**Critical finding (2026-06-17):** All 22 existing ABA/FD/consistency tests currently pass because every multi-link test uses `Vector3d::Zero()` for COM. With zero COM, the `skew(com)*mass` coupling terms are zero, the inertia matrix is block-diagonal, and the Phase 3 correction bug is not triggered. The fix must be developed against tests that actually exercise the bug — requiring non-zero COM test models.

</domain>

<decisions>
## Implementation Decisions

### Test Strategy — TDD (Red → Green → Verify)
- **Red (add failing tests first):** Add non-zero COM variants to `TestDynamicsConsistency.cpp` that fail with current code:
  - Three-link serial chain with COM=[0.1, 0, 0] per link: ID→FD round-trip should recover original qddot to within EPSILON
  - Branching Y-configuration with non-zero COM: same round-trip check
  - Two-link chain with gravity and non-zero COM: ID→FD round-trip
  - These tests must FAIL before the code fix (they'll pass with zero COM, fail with non-zero COM)
- **Green (fix inwardPass):** Restructure to single sweep — once tests fail, the fix is validated when they turn green
- **Verify:** All existing tests continue to pass (no regressions), new non-zero COM tests pass, robot example cross-validation passes
- **Existing tests:** CondensationReducesInertiaNorm and ThreeLinkSingleTorque already exist and will continue to pass — they're inequality/property checks that are still useful
- **Test file scope:** Only `tests/TestDynamicsConsistency.cpp` and `tests/TestForwardDynamics.cpp` modified — no new test files

### Condensation Implementation Approach
- Condensation math implemented inline in `inwardPass()` — single-use, tightly coupled to ABA
- Use existing `Ia.apply(S)` to compute `IaS` (no new method needed on ArticulatedBodyInertia)
- Compute `D = S^T * IaS` as scalar via free function `dot(MotionVector, ForceVector)` (1-DOF joints, revolute), use `1/D` for inversion

### Inward Pass Restructuring
- Single tip-to-base sweep replacing the current 3-phase design. For each link (tip to base):
  1. **Accumulate** — children's condensed Ia/pa already added in prior iterations
  2. **Solve qddot** — `qddot = (tau - S^T*pa) / (S^T*Ia*S)` — this IS the final answer, no correction step
  3. **Condense** — remove joint DOF: `Ia -= (Ia*S)*inv(D)*(Ia*S)^T`, `pa += Ia*S*qddot`
  4. **Pass to parent** — transform condensed Ia/pa to parent frame via `transformInertiaToParent` (kept, computes `X^T*Ia*X`) and `inverseTransformForce`
- Skip condensation for base link (parent == -1) — no parent to pass to
- Keep gravity as `c[0] = -g` on base link bias acceleration (already correct per Featherstone D-07)
- **Forward spatial acceleration pass** lives inside `inwardPass()` after the tip-to-base sweep: `a[i] = X*a_parent + c[i] + S[i]*qddot[i]` — no correction term
- **Critical ordering:** Accumulate BEFORE solving qddot. Solve BEFORE condensing. Condense BEFORE passing to parent.

### Include File Boundary
- **Header changes permitted:** Only `include/ForwardDynamics.h` — Doxygen comment updates on `inwardPass()`, `computeAccelerations()`, and class-level docs to reflect the single sweep design
- Zero changes to `include/ArticulatedBodyInertia.h`, `include/SpatialUtils.h`, `include/LowerTriangular.h`, `include/PluckerTransform.h`, or any other header
- `transformInertiaToParent` helper in anonymous namespace (src/ForwardDynamics.cpp:22-53) is kept — it computes `X^T * Ia * X`, the correct child-to-parent inertia transform for ABA. Do NOT replace with `invtformABI()`

### the agent's Discretion
- Exact non-zero COM values for new test models (e.g., COM=[0.1,0,0] or use Phase 18 UR5-style parameters)
- Exact qddot expected values in ThreeLinkNumericalValidation (use round-trip ID→FD check rather than hardcoded expected values — more robust)
- Doxygen wording specifics (within "single tip-to-base sweep with condensation" description)
- Loop variable naming and intermediate variable naming within the restructured inwardPass()

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Mathematical Reference
- Featherstone, R. (2008). Rigid Body Dynamics Algorithms. **Chapter 7, Algorithm 7.3** — Articulated Body Algorithm with condensation step.
- `.planning/phases/14-cr-02-bug-fix/14-RESEARCH.md` — Full research: condensation formula derivation, Featherstone algorithm trace, reference implementation pseudocode, pitfalls.
- `.planning/phases/12-dynamics-consistency/12-01-SUMMARY.md` — Prior fix context: ABA gravity support, cross-product unification.

### Source Files — Primary Modification Target
- `src/ForwardDynamics.cpp` lines 82-179 — Current inwardPass() with 3-phase structure (initialize → partial qdd with condensation → correct qdd). Replace with single sweep.
- `src/ForwardDynamics.cpp` lines 22-53 — `transformInertiaToParent()` helper (KEEP — computes `X^T * Ia * X`, correct for child-to-parent ABA propagation)

### Source Files — Read Only (API Reference)
- `include/ForwardDynamics.h` — Link struct fields (Ia, pa, S, parent, X, c, v, f, qddot), class Doxygen to update
- `include/ArticulatedBodyInertia.h` — `apply(const MotionVector&)` → ForceVector, `operator+`, `operator*`
- `include/SpatialUtils.h` — `dot(MotionVector, ForceVector)` → double (free function)
- `include/LowerTriangular.h` — `fromFullMatrix(const Eigen::MatrixXd&)` → LowerTriangular
- `include/PluckerTransform.h` — `transformMotion()`, `inverseTransformForce()`, `invtformABI()`

### Test Files
- `tests/TestForwardDynamics.cpp` — 15 tests, modify ThreeLinkNumericalValidation (add non-zero COM, use round-trip check)
- `tests/TestDynamicsConsistency.cpp` — 7 tests, add non-zero COM variants of ThreeLinkSerialChain and BranchingYConfiguration

### Cross-Validation Reference
- `.planning/phases/18-robot-examples/18-01-SUMMARY.md` — v1.2 confirmation: 2-link passes cross-validation, 3-link fails — isolates bug to chains with 3+ joints and non-zero COM

### Existing Test Baseline (2026-06-17)
- ALL existing tests pass with current buggy code (22/22: 15 FD + 7 consistency)
- Reason: every multi-link test uses `RigidBodyInertia(mass, Vector3d::Zero(), ...)` — zero COM avoids triggering the correction bug
- New tests with non-zero COM must FAIL before the code fix is valid

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `ArticulatedBodyInertia::apply(const MotionVector&)` returns `ForceVector` — used for `IaS = Ia.apply(S)`
- `cross(MotionVector, ForceVector)` free function in SpatialUtils.h — used for velocity-product bias forces
- `PluckerTransform::invtformABI(const ArticulatedBodyInertia&)` — transforms child Ia to parent frame
- `PluckerTransform::inverseTransformForce(const ForceVector&)` — transforms child pa to parent frame
- Existing `dot(MotionVector, ForceVector)` free function — used for `S^T * IaS` and `S^T * pa`
- `LowerTriangular::fromFullMatrix()` — wraps `t * t.transpose()` outer products into LT storage

### Established Patterns
- ABA algorithm in `src/ForwardDynamics.cpp` — `outwardPass()` (lines 57-80) and `inwardPass(const VectorXd& tau)` (lines 82-179)
- Current inward pass: Phase 1 initializes Ia/pa (base→tip), Phase 2 accumulates+condenses+propagates (tip→base), Phase 3 corrects qddot (base→tip)
- GTest for all tests, `EXPECT_DOUBLE_EQ` for exact values, `EXPECT_NEAR` with EPSILON=1e-8 for computed values
- Link struct holds both rigid body inertia (I) and articulated body inertia (Ia) simultaneously
- Test pattern for round-trip: ID(qddot_input) → tau → FD(tau) → qddot_output ≈ qddot_input
- Existing ThreeLinkSerialChain test at TestDynamicsConsistency.cpp:100-177 is the template for new non-zero COM variants

### Integration Points
- `inwardPass()` in `src/ForwardDynamics.cpp` lines 82-179 — primary modification target, replace with single sweep
- `outwardPass()` in `src/ForwardDynamics.cpp` lines 57-80 — NO changes, gravity handling is correct
- `computeAccelerations()` in `src/ForwardDynamics.cpp` lines 182-209 — NO changes, calls restructured inwardPass
- `include/ForwardDynamics.h` — Doxygen updates only, NO struct or API changes
- `tests/TestDynamicsConsistency.cpp` — add non-zero COM variants for 3-link and branching tests
- `tests/TestForwardDynamics.cpp` — update ThreeLinkNumericalValidation to use non-zero COM and round-trip check

</code_context>

<specifics>
## Specific Ideas

- The condensation formula: `I_A -= (I_A*S) * inv(S^T*I_A*S) * (I_A*S)^T` then `p_A += I_A*S * qddot + (I_A*S) * inv(S^T*I_A*S) * (tau - S^T*p_A)`
- For 1-DOF revolute joints, `inv(S^T*I_A*S)` is simply `1/D` where `D = S^T * I_A * S`
- The current Phase 3 correction `S^T * Ia_unc * (a_parent + c)` algebraically double-counts: the term `S^T * Ia_unc * c` is already included in `pa` (initialized as `Ia*c + cross(v, Ia*v) + f`), so it gets subtracted twice
- The single sweep eliminates this by computing qddot from fully accumulated Ia/pa only — no correction needed per Featherstone Algorithm 7.3
- Remove storage for `Ia_unc` (vector<ArticulatedBodyInertia>) and `D_store` (vector<double>) — no longer needed without correction step
- For non-zero COM tests: use `RigidBodyInertia(mass, Vector3d(com_x, 0, 0), lt::Identity(3))` — COM offset along link X axis matches the existing 1m link separation pattern
- Round-trip test pattern: ID(qddot_input) → computeTorques → tau, FD(tau) → computeAccelerations → qddot_output, check |qddot_output[i] - qddot_input[i]| < EPSILON for each joint

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

### Related Bugs (separate phases)
- **RNEA fixed-transform limitation** (Phase 19 or future): InverseDynamics uses fixed X transforms that don't update with joint position q — valid only at home configuration. Documented but not fixed.
- **RNEA non-zero qdot tests** (Phase 13, done): RNEA tests with known non-zero qdot values already exist, validating Coriolis/centrifugal computation.

</deferred>

---

*Phase: 14-cr-02-bug-fix*
*Context gathered: 2026-06-17*
