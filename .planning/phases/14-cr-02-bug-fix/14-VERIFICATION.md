# Phase 14: CR-02 Bug Fix — Plan Verification

**Verified:** 2026-06-04
**Plans checked:** 1 (14-01-PLAN.md)
**Status:** ISSUES FOUND — 1 BLOCKER, 1 WARNING, 1 INFO

---

## Goal-Backward Analysis

### Phase Goal (from ROADMAP.md)

> Correct forward dynamics for chains with 3+ joints — all 4 multi-link consistency tests pass (ThreeLinkSerialChain, BranchingYConfiguration, TwoLinkRoundTripWithGravity, ThreeLinkNumericalValidation). All 156 existing tests continue to pass with no regressions.

### What MUST be true for the goal to be achieved:

1. **Condensation step is correctly implemented**: Each link's articulated inertia Ia must have its joint DOF removed (condensed) before being passed to the parent. Formula: `Ia' = Ia - (Ia*S)*(1/D)*(Ia*S)^T`, `pa' = pa + Ia*S*qdd`.

2. **Single sweep eliminates Phase 3 correction**: The 3-phase approach (initialization → partial qdd with condensation → correction) must become a single tip-to-base sweep where qddot is computed as the FINAL value from fully-accumulated Ia/pa.

3. **Child→parent inertia transform is correct**: Condensed Ia must be correctly transformed from child frame to parent frame before accumulation.

4. **No regressions**: All existing tests must continue passing.

---

## Dimension 1: Requirement Coverage

| Requirement | Plans | Status |
|-------------|-------|--------|
| BFIX-01 | 14-01 | Covered |

**Result: ✅ PASS** — Requirement BFIX-01 is listed in plan 14-01 frontmatter and addressed by both tasks.

---

## Dimension 2: Task Completeness

| Task | Type | Files | Action | Verify | Done |
|------|------|-------|--------|--------|------|
| 1 | auto | ✅ | ✅ | ✅ | ✅ |
| 2 | auto | ✅ | ✅ | ✅ | ✅ |

**Result: ✅ PASS** — Both tasks have all required elements. Actions are specific and detailed. Verify commands are runnable. Done criteria are measurable.

---

## Dimension 3: Dependency Correctness

Single plan `depends_on: []`, Wave 1. No dependency graph issues.

**Result: ✅ PASS**

---

## Dimension 4: Key Links Planned

**Result: ✅ PASS** — All key links are documented:
- `Ia.apply(S)` → `ForceVector IaS` ✓
- `dot(S, IaS)` / `dot(S, pa)` via SpatialUtils.h ✓
- `fromFullMatrix` for outer product corrections ✓
- `invtformABI` / `inverseTransformForce` for parent propagation ✓

---

## Dimension 5: Scope Sanity

1 plan, 2 tasks, ~2 files modified. Well within thresholds.

**Result: ✅ PASS**

---

## Dimension 6: Verification Derivation

Truths are user-observable (test outcomes like "ThreeLinkSerialChain passes"), not implementation details. Artifacts map to truths. Key links connect artifacts.

**Result: ✅ PASS**

---

## Dimension 7: Context Compliance

### Locked Decisions (from 14-CONTEXT.md)

| Decision | Plan Implementation | Status |
|----------|-------------------|--------|
| D-01: Condensation math inline | ✓ Implemented in Task 1 action | ✅ |
| D-02: Use existing API | ✓ Ia.apply(S), dot(S, IaS) | ✅ |
| D-03: Single sweep | ✓ Task 1 restructures to single sweep | ✅ |
| D-04: No Phase 3 correction | ✓ Explicitly removed | ✅ |
| D-05: Gravity as c=-g | ✓ Unchanged | ✅ |
| D-06: Skip condensation for base | ✓ parent != -1 guard | ✅ |
| D-07: No include/ changes | ✓ Only src/ + header Doxygen | ✅ |

**Result: ✅ PASS** — All locked decisions are addressed.

---

## Dimension 7b: Scope Reduction Detection

No scope reduction language found. The plan delivers the full condensation step per Featherstone Algorithm 7.3.

**Result: ✅ PASS**

---

## Dimension 7c: Architectural Tier Compliance

All changes are in the API/Backend tier (C++ library), matching the Architectural Responsibility Map.

**Result: ✅ PASS**

---

## Dimension 8: Nyquist Compliance

VALIDATION.md exists at `.planning/phases/14-cr-02-bug-fix/14-VALIDATION.md`.

### Check 8a — Automated Verify Presence

| Task | Plan | Wave | Automated Command | Status |
|------|------|------|-------------------|--------|
| 1 | 14-01 | 1 | `cmake -B build 2>&1 \| tail -3 && cmake --build build 2>&1 \| tail -5` | ✅ |
| 2 | 14-01 | 1 | `cmake --build build 2>&1 \| tail -3 && cd build && ctest --output-on-failure 2>&1 \| tail -30` | ✅ |

### Check 8b — Feedback Latency

Build commands are fast (~30s). No watch-mode flags.

### Check 8c — Sampling Continuity

Only 2 tasks total, no window of 3 consecutive without verify. ✅

### Check 8d — Wave 0 Completeness

VALIDATION.md lists two Wave 0 test gaps:
- `*Condensation*` filter — condensation unit test
- `*TauOrdering*` filter — tau=[1,0,0] ordering test

However, **these tests already exist** in `tests/TestForwardDynamics.cpp`:
- `CondensationReducesInertiaNorm` (line 310)
- `ThreeLinkSingleTorque` (line 371)

These were written during a prior phase and are already compiled and passing. They do NOT need Wave 0 creation. The VALIDATION.md markdown has stale "❌ W0" markers that should be updated to "✅", but this is a documentation issue, not a plan issue.

**Result: ✅ PASS** — Nyquist compliance is met. VALIDATION.md frontmatter should be updated (`nyquist_compliant: true`, `wave_0_complete: true`) but this does not block execution.

---

## Dimension 9: Cross-Plan Data Contracts

Single plan — no cross-plan data contracts to verify.

**Result: ✅ SKIPPED** (single plan)

---

## Dimension 10: AGENTS.md Compliance

From `AGENTS.md`:
- Build: `cmake -B build && cmake --build build` — plan uses this ✓
- Test: `cmake --build build && cd build && ctest --output-on-failure` — plan uses this ✓
- Code conventions: Doxygen comments, `#ifndef` guards, type aliases — plan respects these ✓
- No CI workflows, linter, formatter — plan doesn't introduce any ✓

**Result: ✅ PASS**

---

## Dimension 11: Research Resolution

RESEARCH.md exists. The `## Open Questions` section (lines 311-321) is present but does NOT have the `(RESOLVED)` suffix. While each question has a recommendation, they are not explicitly marked as "RESOLVED":

```markdown
## Open Questions
1. **Exact expected qddot values for ThreeLinkNumericalValidation after fix**
   - Recommendation: Use round-trip consistency...
2. **Contribution ordering for branching trees...**
   - Recommendation: The order is deterministic...
```

The missing `(RESOLVED)` marker is a documentation issue. The questions themselves are addressed (with recommendations and mitigations) in the research.

**Result: ⚠️ WARNING** — Open Questions section lacks `(RESOLVED)` suffix and individual inline `RESOLVED` markers. Fix before research is considered complete, but does not block the plan.

---

## Dimension 12: Pattern Compliance

No `PATTERNS.md` exists for Phase 14.

**Result: ✅ SKIPPED** (no PATTERNS.md)

---

## ⛔ BLOCKER: `invtformABI` vs `transformInertiaToParent` — Incorrect Transform Replacement

### Issue

The plan replaces the anonymous `transformInertiaToParent` function with `PluckerTransform::invtformABI`, assuming these are equivalent transforms. **They are not.**

### Evidence

**`transformInertiaToParent(X, Ia)`** (ForwardDynamics.cpp lines 22-53):
- Builds `X_6x6 = [R, 0; -R*r̂, R]` (motion transform parent→child)
- Computes `X^T * Ia * X`

**`X.invtformABI(Ia)`** (PluckerTransform.cpp lines 165-215):
- Builds `X_inv_6x6 = [R^T, 0; r̂*R^T, R^T]` (inverse motion transform)
- Computes `X^{-1} * Ia * X^{-T}`

These formulas are mathematically different. Verified with identity Ia and R=I, r=[1,0,0]:

```
transformInertiaToParent:  [I - r̂²,  r̂;  -r̂,  I]
invtformABI:               [I,      -r̂;   r̂,   I - r̂²]
```

The off-diagonal blocks have opposite signs; diagonal blocks differ between I and I-r̂². These are **not** interchangeable.

### Why This Matters

The current ABA code uses `transformInertiaToParent` for child→parent inertia propagation and passes **all 2-link round-trip tests** (TwoLinkRoundTripWithGravity at EPSILON=1e-8). The custom function was written specifically for this purpose rather than using `invtformABI`, indicating the original developer recognized they are not equivalent.

Replacing it with `invtformABI` will:
1. **Produce different numerical results** for all chain lengths, including 2-link
2. **Break the currently-passing 2-link round-trip test** unless the two formulas happen to coincide for the specific test configuration (unlikely given the mathematical difference)
3. If it happens to produce "correct" results for multi-link, it will do so for the wrong reasons, making future debugging impossible

### Mathematical Derivation of Correct Transform

From Featherstone, the correct child→parent inertia transform is:

**Ia_parent = X^T * Ia_child * X**

where `X` is the motion transform from parent to child (`links[i].X`).

Proof: 
- `a_child = X * a_parent` (motion transform)
- `f_child = X^{-T} * f_parent` (force transform)  
- `f_child = Ia_child * a_child = Ia_child * X * a_parent`
- `f_parent = X^T * f_child = X^T * Ia_child * X * a_parent`
- Therefore `Ia_parent = X^T * Ia_child * X`

### Fix

The plan should **keep using `transformInertiaToParent`** (or an equivalent inline computation using `X^T * Ia * X`) for the Ia parent propagation. The anonymous namespace function should be preserved, not removed.

The `invtformABI` call on line 196 of the plan action (step 4 of the single sweep) should be replaced with:
```cpp
ArticulatedBodyInertia IaTransformed =
    transformInertiaToParent(links[i].X, links[i].Ia);
```

**Do NOT remove the anonymous namespace** in Task 2 since `transformInertiaToParent` is still needed.

Alternatively, if the planner believes `invtformABI` is correct, it must provide mathematical proof and verify that the 2-link round-trip test (TwoLinkRoundTripWithGravity) continues to pass within EPSILON=1e-8 after the switch.

---

## Summary

| Dimension | Result |
|-----------|--------|
| 1. Requirement Coverage | ✅ PASS |
| 2. Task Completeness | ✅ PASS |
| 3. Dependency Correctness | ✅ PASS |
| 4. Key Links Planned | ✅ PASS |
| 5. Scope Sanity | ✅ PASS |
| 6. Verification Derivation | ✅ PASS |
| 7. Context Compliance | ✅ PASS |
| 7b. Scope Reduction | ✅ PASS |
| 7c. Architectural Tier Compliance | ✅ PASS |
| 8. Nyquist Compliance | ✅ PASS (with VALIDATION.md doc update needed) |
| 9. Cross-Plan Data Contracts | ⏭️ SKIPPED |
| 10. AGENTS.md Compliance | ✅ PASS |
| 11. Research Resolution | ⚠️ WARNING |
| 12. Pattern Compliance | ⏭️ SKIPPED |

---

## Issues

### Blockers (must fix)

**1. [architectural_assumption] `invtformABI` used in place of `transformInertiaToParent`**
- **Plan:** 14-01
- **Task:** 1 (child→parent Ia propagation, step 4 of single sweep)
- **Severity:** BLOCKER
- **Description:** Plan replaces `transformInertiaToParent(X, Ia)` (computes `X^T * Ia * X`) with `X.invtformABI(Ia)` (computes `X^{-1} * Ia * X^{-T}`). These are mathematically different for Plücker transforms with non-zero translation. The custom function was specifically written for the ABA context and produces correct results for 2-link chains. Replacing it will change numerical results for ALL chain lengths.
- **Fix:** Revert to using `transformInertiaToParent` for the Ia child→parent transform. Do not remove the anonymous namespace function in Task 2.
- **Location in plan:** Task 1 action, step 4 ("Pass to parent"): `links[i].X.invtformABI(links[i].Ia)`

### Warnings (should fix)

**1. [research_resolution] Open Questions section not marked RESOLVED**
- **Plan:** Phase-level (RESEARCH.md)
- **Severity:** WARNING
- **Description:** `## Open Questions` section at RESEARCH.md lines 311-321 lacks the `(RESOLVED)` suffix and inline resolution markers. While recommendations are provided, the section should be formally marked resolved.
- **Fix:** Add `(RESOLVED)` to the section heading and inline `RESOLVED` markers to individual questions.

### Info (suggestions)

**1. [validation_doc] VALIDATION.md has stale Wave 0 gaps**
- **Plan:** Phase-level (VALIDATION.md)
- **Severity:** INFO
- **Description:** VALIDATION.md lists two tests as "❌ W0" (Wave 0 missing), but both tests (`CondensationReducesInertiaNorm`, `ThreeLinkSingleTorque`) already exist in `tests/TestForwardDynamics.cpp`. Update frontmatter: `nyquist_compliant: true`, `wave_0_complete: true`.
- **Fix:** Update VALIDATION.md frontmatter and per-task verification map.

---

## Recommendation

**Return to planner for revision.** The BLOCKER issue must be resolved before execution. The fix is straightforward: keep `transformInertiaToParent` instead of switching to `invtformABI`. Do not remove the anonymous namespace.

### If the planner disagrees and believes `invtformABI` is correct:

The planner must provide:
1. A mathematical proof that `X^{-1} * Ia * X^{-T} = X^T * Ia * X` for Plücker transforms (unlikely to exist for non-zero translation)
2. Evidence that TwoLinkRoundTripWithGravity continues to pass at EPSILON=1e-8 after the switch
3. An explanation for why the original ABA code used a custom function instead of `invtformABI`
