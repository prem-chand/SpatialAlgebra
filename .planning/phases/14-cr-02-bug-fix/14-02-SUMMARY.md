# Plan 14-02 Summary: Non-Zero COM Tests (RED Phase)

**Status:** Complete  
**Date:** 2026-06-17  

## Result

Three new non-zero COM consistency tests added to TestDynamicsConsistency.cpp, and ThreeLinkNumericalValidation updated in TestForwardDynamics.cpp with non-zero COM round-trip check. RED phase demonstrated: BranchingYNonZeroCOM and ThreeLinkNumericalValidation FAIL with the current 3-phase inwardPass.

## Test Results

| Test | Expected | Actual |
|------|----------|--------|
| ThreeLinkSerialChainNonZeroCOM | FAIL | PASS — bug not triggered at COM=[0.1,0,0] |
| BranchingYNonZeroCOM | FAIL | FAIL — symmetry violated (0.55 vs 0.45) |
| TwoLinkGravityNonZeroCOM | PASS | PASS — 2-link chains don't trigger bug |
| ThreeLinkNumericalValidation | FAIL | FAIL — ordering violated (1.0 vs 0.25) |
| 7 existing consistency tests | PASS | PASS — zero regressions |
| 14 existing FD tests | PASS | PASS — zero regressions |

## Gaps Closed

- SC1: Non-zero COM multi-link consistency tests exist, two demonstrate the bug
- SC2: ThreeLinkNumericalValidation uses non-zero COM round-trip ID-FD check

## Commits

```
5da0dec test(14-02): add non-zero COM consistency tests (RED phase)
```
