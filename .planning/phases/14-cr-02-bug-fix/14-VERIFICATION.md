---
phase: 14-cr-02-bug-fix
verified: 2026-06-17T00:00:00Z
status: gaps_found
score: 2/10 must-haves verified
overrides_applied: 0
overrides: []
re_verification:
  previous_status: ISSUES_FOUND
  previous_blocker: "invtformABI vs transformInertiaToParent (PLAN checker — moot after implementation divergence)"
  gaps_remaining:
    - "ROADMAP SC1: No non-zero COM multi-link consistency tests exist"
    - "ROADMAP SC2: ThreeLinkNumericalValidation not updated for non-zero COM round-trip"
    - "ROADMAP SC3: inwardPass() still 3-phase with Ia_unc, D_store, Phase 3 correction"
    - "ROADMAP SC5: pa initialized before child accumulation (Featherstone ordering violation)"
  regressions: []
gaps:
  - truth: "All multi-link consistency tests with non-zero COM pass (ThreeLinkSerialChainNonZeroCOM, BranchingYNonZeroCOM, TwoLinkGravityNonZeroCOM)"
    status: failed
    reason: "No non-zero COM multi-link tests exist in TestDynamicsConsistency.cpp or TestForwardDynamics.cpp. All 7 consistency tests + 15 FD tests use Vector3d::Zero() for COM. The bug path (Phase 3 correction double-counting) is only triggered with non-zero COM per CONTEXT.md lines 12-13."
    artifacts:
      - path: tests/TestDynamicsConsistency.cpp
        issue: "8 tests total; ZERO use non-zero COM in multi-link chains"
      - path: tests/TestForwardDynamics.cpp
        issue: "16 tests total; 3 single-link tests use non-zero COM (GravityEffectScalesWithMass, GravityProportionalityInvariant, SingleLinkGravityNoZeroCOM); ZERO multi-link tests use non-zero COM"
    missing:
      - "Add ThreeLinkSerialChainNonZeroCOM test (3-link serial chain, COM=[0.1,0,0] per link, round-trip ID→FD check)"
      - "Add BranchingYNonZeroCOM test (branching Y, non-zero COM, round-trip + symmetric branch check)"
      - "Add TwoLinkGravityNonZeroCOM test (2-link with gravity, non-zero COM, round-trip check)"

  - truth: "ThreeLinkNumericalValidation updated with non-zero COM round-trip check — passes within EPSILON=1e-8"
    status: failed
    reason: "ThreeLinkNumericalValidation (TestForwardDynamics.cpp line 265) still uses RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3)) for all 3 links. Test asserts only finiteness and qddot ordering — no round-trip ID→FD check."
    artifacts:
      - path: tests/TestForwardDynamics.cpp
        issue: "ThreeLinkNumericalValidation at line 265: Vector3d::Zero() COM, no round-trip consistency check per PLAN Task 1 requirements"
    missing:
      - "Replace Vector3d::Zero() with Vector3d(0.1, 0, 0) on all 3 links"
      - "Implement round-trip: ID(qddot_input) → tau → FD(tau) → qddot_output → EXPECT_NEAR(qddot_output, qddot_input, 1e-8)"
      - "Verify qddot[0] < qddot[2] (base carries more reflected inertia)"

  - truth: "ABA inward pass restructured to single tip-to-base sweep per Featherstone Algorithm 7.3 (no Phase 3 correction)"
    status: failed
    reason: "inwardPass() in src/ForwardDynamics.cpp retains the 3-phase structure. Ia_unc vector exists at line 107, D_store at line 106, Phase 3 correction loop at lines 160-179. The comment at line 128 reads 'This is a PARTIAL acceleration (assumes a_parent == 0)', confirming the correction-based approach was not replaced."
    artifacts:
      - path: src/ForwardDynamics.cpp
        issue: "Lines 106-107: Ia_unc + D_store storage. Lines 160-179: Phase 3 correction loop (a_prime, correction/D, qddot adjustment). Line 128: 'PARTIAL acceleration' comment."
    missing:
      - "Remove Ia_unc vector (line 107) and D_store vector (line 106)"
      - "Remove Phase 3 correction loop (lines 160-179)"
      - "Compute FINAL qddot during the tip-to-base sweep (no correction step)"
      - "Add forward spatial acceleration pass after the sweep: a[i] = X*a_parent + c[i] + S[i]*qddot[i]"

  - truth: "Bias forces (pa) computed after child inertia accumulation per Featherstone Algorithm 7.3"
    status: failed
    reason: "Current inwardPass() Phase 1 (lines 84-102) initializes pa for ALL links base-to-tip BEFORE any child accumulation. Featherstone Algorithm 7.3 requires pa computation AFTER children's condensed contributions are accumulated. The single-sweep design (compute pa in the tip-to-base loop after accumulation) resolves this."
    artifacts:
      - path: src/ForwardDynamics.cpp
        issue: "Lines 84-102: Phase 1 initializes Ia and pa for all links before Phase 2 tip-to-base sweep. Per Featherstone, pa should incorporate children's condensed contributions first."
    missing:
      - "Move pa initialization into the tip-to-base sweep, after children accumulation"
      - "Or restructure to single sweep where accumulation, qddot solve, condensation, and pa update happen per link in order"

  - truth: "Cross product formulas in SpatialUtils.h are correct per Featherstone §2.9"
    status: verified
    reason: "cross(MotionVector, ForceVector) at SpatialUtils.h:78 computes [ω×τ + v×f; ω×f] — correct per Featherstone §2.9 crf formula. cross(ForceVector, ForceVector) at line 116 computes [τ₁×τ₂ + f₁×f₂; τ₁×f₂ - τ₂×f₁] — includes anti-commutativity term. These were the Phase 14-01 cross product fixes."
    artifacts:
      - path: include/SpatialUtils.h
        issue: ""

  - truth: "All 11 test executables pass with zero regressions (existing tests)"
    status: verified
    reason: "ctest --output-on-failure reports 11/11 tests pass, 0 fail, 0 skipped. All existing tests (15 FD + 7 consistency) pass. However, these are zero-COM tests that do not exercise the Phase 3 correction bug path."
    artifacts:
      - path: build/ (all test executables)
        issue: ""

deferred: []
---

# Phase 14: CR-02 Bug Fix — Re-Verification Report

**Phase Goal (ROADMAP.md):** Correct forward dynamics for chains with 3+ joints and non-zero COM — TDD with Red→Green→Verify cycle
**Verified:** 2026-06-17
**Status:** gaps_found
**Re-verification:** Yes — after divergence from original plan

## Prior Verification Context

The previous VERIFICATION.md (2026-06-04) was a **PLAN CHECKER** — it audited 14-01-PLAN.md before execution. Its only BLOCKER found: the plan proposed replacing `transformInertiaToParent` with `invtformABI`, which are mathematically different transforms.

**What changed:** The 14-01-SUMMARY.md reports that execution DIVERGED from the plan — the ABA `inwardPass()` was NOT restructured. Instead, an independent diagnosis found cross product formula bugs in `SpatialUtils.h` as the root cause. The `invtformABI` replacement was never implemented, making the previous BLOCKER moot (it applied to a plan that was not followed).

**This re-verification checks the ACTUAL codebase against ROADMAP success criteria, not against the superseded plan.**

## Goal Achievement

### Roadmap Success Criteria

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| SC1 | All multi-link consistency tests with non-zero COM pass | ✗ FAILED | No non-zero COM multi-link tests exist anywhere in codebase |
| SC2 | ThreeLinkNumericalValidation updated with non-zero COM round-trip | ✗ FAILED | Still uses `Vector3d::Zero()` for all COM, no round-trip check |
| SC3 | ABA inward pass restructured to single sweep (no Phase 3 correction) | ✗ FAILED | 3-phase structure intact: Ia_unc (line 107), D_store (line 106), Phase 3 correction (lines 160-179) |
| SC4 | All 11 test executables pass with zero regressions | ✓ PASSED | ctest: 11/11 passed, 0 failed, 0 skipped |
| SC5 | Bias forces computed after child accumulation per Featherstone 7.3 | ✗ FAILED | pa initialized Phase 1 (lines 84-102) before any child accumulation in Phase 2 |

**Score:** 1/5 roadmap success criteria met

### Plan Must-Have Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Non-zero COM 3-link serial chain round-trip ≈ qddot | ✗ FAILED | No test exists |
| 2 | Non-zero COM branching Y round-trip ≈ qddot | ✗ FAILED | No test exists |
| 3 | Non-zero COM 2-link with gravity round-trip ≈ qddot | ✗ FAILED | No test exists |
| 4 | Non-zero COM 3-link: base qddot < tip qddot | ✗ FAILED | No test exists |
| 5 | All 11 test executables pass with zero failures | ✓ VERIFIED | 11/11 ctest pass |
| 6 | inwardPass() is single sweep — no Ia_unc, no D_store, no Phase 3 | ✗ FAILED | All three artifacts present |

**Score:** 1/6 plan must-haves verified

### Overall Score: 2/10 must-haves verified

## What WAS Implemented (Cross Product Fixes)

The following were fixed in Phase 14-01 (cross product root cause fix) and are correctly implemented:

| Fix | File | Line | Status |
|-----|------|------|--------|
| `cross(MotionVector, ForceVector)` — crf formula | `include/SpatialUtils.h` | 78 | ✓ Correct: `[ω×τ + v×f; ω×f]` per Featherstone §2.9 |
| `cross(ForceVector, ForceVector)` — anti-commutativity | `include/SpatialUtils.h` | 116 | ✓ Correct: `[τ₁×τ₂ + f₁×f₂; τ₁×f₂ - τ₂×f₁]` |
| Virtual destructor on SpatialVector | `include/SpatialVector.h` | — | ✓ Added |
| Namespace cleanup | `include/LowerTriangular.h`, `include/Rotation.h` | — | ✓ Global `using` removed |
| InverseDynamics O(n²)→O(n) | `src/InverseDynamics.cpp` | — | ✓ Tip-to-base sweep |

These fixes are necessary but insufficient to meet the ROADMAP goal. They fixed bias force corruption that affected ALL chains, but the Phase 3 correction bug (specific to 3+ joint chains with non-zero COM per CONTEXT.md) was not addressed.

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `tests/TestDynamicsConsistency.cpp` | 9 tests (7 existing + 2 new non-zero COM) | ⚠️ PARTIAL | 8 tests, all zero-COM multi-link; missing ThreeLinkSerialChainNonZeroCOM, BranchingYNonZeroCOM |
| `tests/TestForwardDynamics.cpp` | 18 tests (16 existing + 2 updated/new non-zero COM) | ⚠️ PARTIAL | 16 tests; ThreeLinkNumericalValidation not updated for non-zero COM; no round-trip check |
| `src/ForwardDynamics.cpp` | Single sweep inwardPass() (180+ lines) | ✗ NOT DONE | 3-phase structure intact (210 lines); Ia_unc/D_store/Phase 3 present |
| `include/ForwardDynamics.h` | Updated Doxygen for single sweep | ✗ NOT DONE | Still describes 3-phase ABA |

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `inwardPass()` | `Ia.apply(S)` | `ForceVector IaS` | ✓ WIRED | `src/ForwardDynamics.cpp:114` |
| `inwardPass()` | `dot(S, IaS)` | `SpatialUtils.h:65` | ✓ WIRED | `src/ForwardDynamics.cpp:115` |
| `inwardPass()` | `cross(MotionVector, ForceVector)` | `SpatialUtils.h:78` | ✓ WIRED | `src/ForwardDynamics.cpp:96` |
| `inwardPass()` | `transformInertiaToParent` | Anonymous namespace helper | ✓ WIRED | Defined line 22, called line 149 |
| `inwardPass()` | `lt::fromFullMatrix` | Outer product condensation | ✓ WIRED | `src/ForwardDynamics.cpp:139-141` |
| `inwardPass()` | `inverseTransformForce` | Pa parent propagation | ✓ WIRED | `src/ForwardDynamics.cpp:152` |
| New non-zero COM tests → `computeTorques`/`computeAccelerations` | Round-trip ID↔FD | (tests don't exist) | ✗ NOT_WIRED | No non-zero COM multi-link tests to connect |

## Data-Flow Trace (Level 4)

Only existing (zero-COM) data paths are verified. The non-zero COM paths cannot be traced because the tests don't exist.

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|--------------|--------|-------------------|--------|
| `inwardPass()` qddot | `links[i].qddot` | Phase 3 correction: `u/D - correction/D_store` | ✓ (zero-COM chains) | ⚠️ FLOWING but unvalidated for non-zero COM |
| `inwardPass()` pa | `links[i].pa` | Phase 1 init + Phase 2 condensation | ✓ (zero-COM chains) | ⚠️ FLOWING but initialized before accumulation (Featherstone ordering violation) |
| `cross(MotionVector, ForceVector)` | bias force crf term | `pa = Ia*c + cross(v, Ia*v) + f` | ✓ Fixed formula | ✓ FLOWING (Phase 14-01 fix) |

## Anti-Patterns Found

No debt markers (TBD, FIXME, XXX, TODO, HACK) found in any files modified by this phase. No placeholder or stub patterns detected.

## Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Build succeeds | `cmake -B build && cmake --build build` | Exit 0, all targets built | ✓ PASS |
| All 11 test executables pass | `cd build && ctest --output-on-failure` | 11/11 passed, 0 failed | ✓ PASS |
| ThreeLinkSerialChain (zero-COM) passes | `./TestDynamicsConsistency --gtest_filter=*ThreeLinkSerialChain*` | PASSED | ✓ PASS |
| BranchingYConfiguration (zero-COM) passes | `./TestDynamicsConsistency --gtest_filter=*Branching*` | PASSED | ✓ PASS |
| TwoLinkRoundTripWithGravity (zero-COM) passes | `./TestDynamicsConsistency --gtest_filter=*Gravity*` | PASSED | ✓ PASS |
| Non-zero COM 3-link round-trip | (test does not exist) | N/A | ? SKIP |

## Requirements Coverage

| Requirement | Source | Description | Status | Evidence |
|-------------|--------|-------------|--------|----------|
| BFIX-01 | `.planning/milestones/v1.2-REQUIREMENTS.md` | CR-02 ABA inward pass restructured per Featherstone Algorithm 7.3 — all 4 multi-link consistency tests pass | ✗ BLOCKED | inwardPass() not restructured; 4 multi-link tests pass only with zero COM; no non-zero COM tests exist. REQUIREMENTS.md v1.2 still shows BFIX-01 as "Pending" |

## Human Verification Required

No human verification items needed — all gaps are objectively detectable in codebase.

## Gaps Summary

The phase has a **fundamental execution gap**: the ROADMAP.md and PLAN.md both specify a TDD-driven restructuring of the ABA `inwardPass()` with non-zero COM tests, but the implementation treated cross product formula fixes as sufficient. Critical evidence:

1. **No non-zero COM multi-link tests** — The CONTEXT.md (lines 12-13) explicitly states: *"All 22 existing ABA/FD/consistency tests currently pass because every multi-link test uses Vector3d::Zero() for COM. With zero COM, the skew(com)*mass coupling terms are zero, the inertia matrix is block-diagonal, and the Phase 3 correction bug is not triggered."* Without non-zero COM tests, there is no evidence that forward dynamics works for chains with 3+ joints and non-zero COM.

2. **Phase 3 correction intact** — The `inwardPass()` still has the 3-phase structure with `Ia_unc`, `D_store`, and Phase 3 correction. The comment at line 128 reads *"This is a PARTIAL acceleration (assumes a_parent == 0)"* — confirming the correction-based approach was not replaced.

3. **Featherstone ordering violation** — `pa` is initialized for ALL links in Phase 1 (base-to-tip) BEFORE children's contributions are accumulated in Phase 2 (tip-to-base). Featherstone Algorithm 7.3 requires computing pa AFTER accumulation.

**The cross product fixes (SpatialUtils.h) are correct and valuable** — they fixed bias force corruption that would have affected any chain. But they are a prerequisite fix, not the complete solution. The ROADMAP goal requires non-zero COM multi-link tests that pass, which don't exist.

### Root Cause

The 14-01-SUMMARY.md states: *"The original plan's core premise (ABA inward pass needs restructuring) was incorrect. The plan should not be executed — it would undo valid working code."* This conclusion was reached without non-zero COM test evidence. The cross product fix may be sufficient for correctness, or it may not — we cannot tell without the non-zero COM tests that the PLAN specified.

---

_Verified: 2026-06-17_
_Verifier: the agent (gsd-verifier)_
