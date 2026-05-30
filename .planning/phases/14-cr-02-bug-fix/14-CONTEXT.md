# Phase 14: CR-02 Bug Fix - Context

**Gathered:** 2026-05-30
**Status:** Ready for planning

<domain>
## Phase Boundary

Correct forward dynamics for chains with 3+ joints by implementing the missing articulation (condensation) step from Featherstone Algorithm 7.3 in the ABA inward pass. All 4 multi-link consistency tests must pass, and all 156 existing tests must continue to pass with no regressions.

</domain>

<decisions>
## Implementation Decisions

### Condensation Implementation Approach
- Condensation math implemented inline in `inwardPass()` — single-use, tightly coupled to ABA
- Use existing `Ia.apply(S)` to compute `IaS` (no new method needed on ArticulatedBodyInertia)
- Compute `D = S^T * IaS` as scalar (1-DOF joints), use `1/D` for inversion

### Inward Pass Restructuring
- Single tip-to-base sweep: initialize Ia/pa, accumulate children contributions, compute qddot, condense Ia/pa, pass to parent
- Keep gravity as `c = -g` on base link bias acceleration (already correct per Featherstone D-07)
- Skip condensation for base link (parent == -1) — no parent to pass to

### Test Strategy
- Update `ThreeLinkNumericalValidation` with exact expected qddot values after fix
- Add condensation unit test: 3-link chain, verify condensed Ia norm < uncondensed Ia norm
- Add test with specific tau=[1,0,0] on 3-link chain, verify qddot[0] < qddot[2]

### the agent's Discretion
All implementation details not covered above are at the agent's discretion.

</decisions>

<code_context>
## Existing Code Insights

### Reusable Assets
- `ArticulatedBodyInertia::apply(const MotionVector&)` returns `ForceVector` — used for `IaS = Ia.apply(S)`
- `cross(v, IaV)` free function in SpatialUtils.h — used for velocity-product bias forces
- `PluckerTransform::invtformABI(const ArticulatedBodyInertia&)` — transforms child Ia to parent frame
- `PluckerTransform::inverseTransformForce(const ForceVector&)` — transforms child pa to parent frame
- Existing `dot(MotionVector, ForceVector)` free function — used for S^T * pa and S^T * IaS
- `LowerTriangular` has `add(const LowerTriangular&)`, `scale(double)`, `subtract(const LowerTriangular&)`

### Established Patterns
- ABA algorithm in `src/ForwardDynamics.cpp` — `outwardPass()` and `inwardPass(const VectorXd& tau)`
- Three-phase inward pass: initialize (base-to-tip), accumulate (tip-to-base), solve (tip-to-base)
- GTest for all tests, `EXPECT_DOUBLE_EQ` for floats, `EXPECT_NE` for inequality
- Link struct holds both rigid body inertia (I) and articulated body inertia (Ia) simultaneously

### Integration Points
- `inwardPass()` in `ForwardDynamics.cpp` lines 53-126 — primary modification target
- `ThreeLinkNumericalValidation` test at TestForwardDynamics.cpp lines 261-312 — update expected values
- `TestDynamicsConsistency.cpp` — multi-link round-trip tests will auto-verify after fix (no changes needed)
- `ArticulatedBodyInertia.h` — no changes needed, existing API sufficient for condensation math

</code_context>

<specifics>
## Specific Ideas

- The condensation formula: `I_A -= (I_A*S) * inv(S^T*I_A*S) * (I_A*S)^T` then `p_A += I_A * S * qddot + (I_A*S) * inv(S^T*I_A*S) * (tau - S^T*p_A)`
- For 1-DOF joints, `inv(S^T*I_A*S)` is simply `1/D` where `D = S^T * I_A * S`
- The condensation reduces reflected inertia at the parent — without it, the parent overestimates mass
- Phase restructuring: merge the 3-phase inward pass into a single tip-to-base sweep

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>
