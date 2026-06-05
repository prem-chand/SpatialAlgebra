# Phase 12: Dynamics Consistency — Execution Summary

## Status: ⚠️ PARTIALLY COMPLETE (BF-02 deferred to v1.2)

**Duration:** ~2 hours  
**Commits:** 1

---

## Summary

Fixed ABA bias acceleration propagation, but RNEA↔ABA consistency for multi-link chains remains unresolved.

### Fix Applied

**ForwardDynamics::outwardPass()** — Added bias acceleration propagation from parent to child:
```cpp
// BEFORE (buggy):
links[i].c = cross(links[i].v, links[i].S) * links[i].qdot;

// AFTER (fixed):
MotionVector cParent = links[parent].c;
links[i].c = links[i].X.inverseTransformMotion(cParent) + 
             cross(links[i].v, links[i].S) * links[i].qdot;
```

This matches Featherstone Algorithm 7.3: `cᵢ = Xᵢ⁻¹·c_parent + vᵢ × Sᵢ·q̇ᵢ`

### Test Results

| Test | Before | After |
|------|--------|-------|
| ConsistencyTest.RoundTripABARNEA | PASS | PASS ✅ |
| ConsistencyTest.RoundTripRNEAABA | PASS | PASS ✅ |
| ConsistencyTest.ThreeLinkSerialChain | FAIL | FAIL ⚠️ |
| ConsistencyTest.BranchingYConfiguration | FAIL | FAIL ⚠️ |

**Overall:** 2/4 tests passing (50%)

---

## Root Cause Analysis

**Issue:** RNEA and ABA produce inconsistent results for 3+ link chains.

**Investigation:**
- Single-link consistency works (both produce correct results)
- Multi-link fails: RNEA computes incorrect torques, ABA computes incorrect accelerations
- Example: For 3-link chain with qddot=[1, 0.5, 0.25]:
  - Expected tau ≈ [4.25, 2.25, 1.75]
  - RNEA gives tau = [-1.75, 0.75, 1.75]
  - ABA recovers qddot = [-1.75, 0.75, 1.75] (matches tau, not original qddot)

**Suspected Causes:**
1. **Force propagation direction** — RNEA may be transforming forces in wrong direction
2. **Frame conventions** — RNEA and ABA may use different spatial vector frame conventions
3. **Articulated inertia accumulation** — ABA tformABI may have subtle issues

**Impact:** Forward and inverse dynamics disagree for multi-link systems.

---

## v1.1 Final Status

| Requirement | Status | Tests |
|-------------|--------|-------|
| BF-01: ABI Transforms | ✅ Complete | 5/5 fixed |
| BF-02: Multi-link Consistency | ⚠️ Deferred | 2/4 passing |
| ST-01: 100% Test Pass Rate | ⚠️ Blocked | 156/158 (98.7%) |

---

## Recommendation for v1.2

**Priority:** High — blocks production use for multi-link robots

**Approach:**
1. Derive RNEA and ABA equations side-by-side from Featherstone
2. Verify frame conventions match (body vs spatial)
3. Add intermediate debug output to trace force/acceleration propagation
4. Consider adding analytical test cases with known solutions

---

*Completed: 2026-05-16*
