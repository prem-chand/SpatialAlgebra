---
phase: 13
reviewers: [codex]
reviewed_at: 2026-05-17T07:23:00Z
plans_reviewed: [13-01-PLAN.md, 13-02-PLAN.md, 13-03-PLAN.md, 13-04-PLAN.md, 13-05-PLAN.md]
---

# Cross-AI Plan Review — Phase 13

## Codex Review (gpt-5.4)

### 13-01-PLAN.md

**Summary:** This is the most important plan in the phase because it targets correctness bugs in core algebra and a likely ABI transformation defect. The direction is good: centralizing the force cross-product logic reduces drift, and fixing ArticulatedBodyInertia composition before expanding tests is the right order.

**Strengths:**
- Consolidates duplicated cross-product logic into one implementation point
- Prioritizes correctness fixes before adding dependent tests
- Identifies a realistic type/ODR hazard in PluckerTransform::apply()
- Keeps scope contained to algebraic core and transform cleanup

**Concerns:**
- HIGH: Removing MotionVector::crossForce and ForceVector::crossMotion may be an API break if downstream callers use them
- HIGH: The plan does not explicitly require derivation-backed tests for the unified cross-product formula
- MEDIUM: ABI operator+ review is underspecified; "fix if needed" is too weak for a known correctness-sensitive area
- MEDIUM: "ODR/linkage hazard" is not tied to a concrete failure mode

**Suggestions:**
- Require explicit algebraic tests for force×force, including anti-symmetry and known-hand-computed cases
- Deprecate overloads first unless the phase explicitly allows breaking API changes
- Expand Task 2 into concrete checks: argument order, COM shift term, mass scaling, symmetry preservation

**Risk:** MEDIUM-HIGH

---

### 13-02-PLAN.md

**Summary:** Addresses a real production-readiness gap. Featherstone framing is correct, backward-compatible defaults are reasonable. Risk is narrow file scope and ambiguous gravity conventions.

**Strengths:**
- Targets an essential missing capability for real robot dynamics
- Uses Featherstone-consistent base acceleration formulation
- Preserves backward compatibility with zero-gravity defaults
- Couples ABA and RNEA changes so dynamics stay aligned

**Concerns:**
- HIGH: Adding optional parameters only in ForwardDynamics/InverseDynamics may miss callers, wrappers, tests that assume old signature
- HIGH: Gravity frame convention is unspecified — world-frame vs base-frame gravity direction must be nailed down
- MEDIUM: Defaulting to Vector3d::Zero() can silently mask missing gravity in production
- MEDIUM: No mention of validating single-link static equilibrium

**Suggestions:**
- Define gravity semantics explicitly: coordinate frame, sign convention
- Add tests for zero torque under free fall, static holding torque under gravity
- Consider overloads or a config object if API clarity suffers

**Risk:** MEDIUM

---

### 13-03-PLAN.md

**Summary:** Useful but modest in direct production impact. Debug NaN/Inf guards and test helper fixes improve diagnosis and test quality.

**Strengths:**
- Adds cheap developer-facing correctness guards in high-value math paths
- Fixes suspicious test helpers instead of building more tests on bad fixtures
- Correctly depends on 13-01 before extending cross-force coverage
- Keeps runtime cost limited to debug builds

**Concerns:**
- MEDIUM: #ifndef NDEBUG assertions do nothing in release builds
- MEDIUM: Guard coverage is selective; only constructor and apply()
- MEDIUM: Fixing helper generators could invalidate existing expectations

**Suggestions:**
- Reframe as developer diagnostics and test correctness, not production hardening
- Audit all helper usages and update affected expected values
- Add release-mode tests verifying behavior stays stable when inputs are valid

**Risk:** LOW-MEDIUM

---

### 13-04-PLAN.md

**Summary:** This is the plan that proves the phase actually achieved its goal. Multi-link ABA, non-zero velocity RNEA, and gravity round-trip tests are exactly the kind of system-level coverage missing from many math libraries.

**Strengths:**
- Directly tests multi-body behavior rather than only unit-level algebra
- Covers Coriolis/centrifugal terms
- Adds cross-consistency checks between ABA and RNEA
- Well-ordered after algebra and gravity fixes

**Concerns:**
- HIGH: Reference values source is underspecified — tests may be self-confirming
- HIGH: Round-trip consistency alone is insufficient; ABA and RNEA can agree while both being wrong
- MEDIUM: Only "2+ link chains" may still be too narrow
- MEDIUM: No mention of edge cases like zero mass/inertia

**Suggestions:**
- Use independently derived expected results (textbook examples, symbolic derivation)
- Include static gravity equilibrium case, dynamic non-zero qdot case, asymmetric chain with COM offset
- Add tolerance strategy explicitly with justified epsilons
- Include legacy zero-gravity regression cases

**Risk:** MEDIUM-HIGH

---

### 13-05-PLAN.md

**Summary:** Mixes useful productionization work with some likely scope creep. CI, dependency cleanup, build-system hardening clearly belong. Namespace/umbrella/guard cleanup may be lower-value style items.

**Strengths:**
- Adds CI, required for real production-readiness
- Addresses Eigen version friction called out in project instructions
- Adds GTest fallback reducing contributor setup pain
- Coverage option plus CI matrix improves cross-platform confidence

**Concerns:**
- HIGH: Mixed bag of build fixes, API surface changes, style cleanup — broadest scope
- MEDIUM: Moving Vector3d into SpatialAlgebra may be source-breaking for consumers
- MEDIUM: Deleting stub .cpp files needs CMake/test discovery verification
- MEDIUM: CodeCov upload introduces third-party integration dependency

**Suggestions:**
- Split into build/CI hardening vs namespace/API cleanup tracks
- Prioritize only items tied directly to production-readiness
- Make CI stages incremental: build + test first, coverage second
- Add consumer-style build test in CI that installs/includes public headers

**Risk:** MEDIUM (HIGH scope-creep risk if not trimmed)

---

### Overall Assessment

The phase is directionally strong. Plans 13-01, 13-02, and 13-04 are the core and map to the stated goal. The main quality gap across the set is insufficient explicitness around conventions and validation sources: cross-product signs, gravity frame semantics, and independent dynamics references all need to be nailed down.

### Cross-Plan Strengths
- Dependency ordering is mostly sensible: algebra fixes before dependent tests
- Covers both implementation and verification
- Backward compatibility considered in gravity API design
- CI/build reliability recognized as part of production readiness

### Cross-Plan Concerns
- HIGH: No plan explicitly defines mathematical conventions that tests and APIs must share
- HIGH: Several validations rely on "consistency" rather than independent correctness oracles
- MEDIUM: Some cleanup items look opportunistic rather than phase-critical
- MEDIUM: Runtime safety still limited (debug assertions not enough for production)
- LOW: No explicit documentation updates for changed APIs

### Overall Risk Assessment

MEDIUM-HIGH — The phase can achieve its goal, but only with discipline about correctness oracles and avoiding letting cleanup work displace math validation.

---

*Review by Codex CLI (gpt-5.4) — 2026-05-17*
