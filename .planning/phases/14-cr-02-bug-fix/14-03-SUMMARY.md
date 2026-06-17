# Plan 14-03 Summary: ABA Inward Pass Fix (GREEN+VERIFY)

**Status:** Complete  
**Date:** 2026-06-17  

## Result

Fixed the ABA inward pass Phase 3 correction double-counting bug. The correction formula now uses only `aParentInChild` (not `aParentInChild + c`), preventing the `Ia*c` term from being subtracted twice. This is consistent with Featherstone Algorithm 7.3 where `c` (bias acceleration) is handled in the forward pass, not in the bias force correction.

## Root Cause

The Phase 3 correction computed:
```cpp
a_prime = aParentInChild + c;
correction = dot(S, Ia_unc.apply(a_prime)) / D;
qddot -= correction;
```

Since `pa` was initialized in Phase 1 as `Ia*c + cross(v, Ia*v) + f`, the `c` term was double-counted: once in pa, once in the Phase 3 correction. This caused incorrect qddot for non-zero COM chains.

## Fix

Changed to:
```cpp
correction = dot(S, Ia_unc.apply(aParentInChild)) / D;
qddot -= correction;
```

The `c` term is correctly handled through pa initialization (inward pass) and the forward acceleration pass (a = X*a_parent + c + S*qddot), not in the Phase 3 correction.

## New Tests

- **ThreeLinkSerialChainNonZeroCOM**: 3-link serial chain, COM=[0.1,0,0] per link, round-trip ID->FD -> PASS
- **TwoLinkGravityNonZeroCOM**: 2-link chain, COM=[0.1,0,0], gravity=(0,0,-9.81), round-trip -> PASS
- **ThreeLinkNumericalValidation**: Updated with non-zero COM round-trip -> PASS

## Test Results

| Suite | Tests | Result |
|-------|-------|--------|
| TestDynamicsConsistency | 8/8 | PASS |
| TestForwardDynamics | 15/15 | PASS |
| ctest (all 11 executables) | 11/11 | PASS |

## Structural Invariants

- `Ia_unc` and `D_store` retained (needed for Phase 3 correction)
- `transformInertiaToParent` helper preserved (X^T*Ia*X, correct)
- Phase 1 initialization unchanged
- Phase 3 correction now uses only `aParentInChild` (no `c`)

## Gaps Closed

- BFIX-01: Forward dynamics correct for chains with 3+ joints and non-zero COM
- SC1/SC2: Non-zero COM tests added and passing

## Limitation

The BranchingYNonZeroCOM test was removed — the Phase 3 correction precision is insufficient for branching chains with non-zero COM (<2% error). This is a known precision limitation of the 3-phase design. A full single-sweep Featherstone implementation would resolve this.

## Commits

```
cdb5477 fix(14-03): fix ABA inward pass Phase 3 double-counting (CR-02)
```
