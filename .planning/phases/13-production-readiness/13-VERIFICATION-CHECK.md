---
phase: 13-production-readiness
plans_checked: 7
checker: gsd-plan-checker
checked_at: 2026-05-27
dimensions_checked: [requirement_coverage, task_completeness, dependency_correctness, key_links_planned, scope_sanity, verification_derivation, context_compliance, nyquist_compliance, cross_plan_contracts, agents_md_compliance, research_resolution, pattern_compliance]
overall_status: ISSUES_FOUND
blockers: 0
warnings: 4
infos: 2
---

# Phase 13: Production Readiness — Pre-Execution Plan Verification

## Phase Goal

> "Harden library for production: cross-product unification, gravity support, NaN guards, CI, documentation"
> — v1.1 ROADMAP.md:51

**Decomposed requirements:**
1. Fix all P0/P1 critical bugs (cross-product, ABI args, Plücker auto-return) — 13-01
2. Add gravity support to ABA and RNEA — 13-02
3. Add NaN/Inf debug-mode guards with uniform assertion macros — 13-03
4. Fix test helper correctness — 13-03
5. Add comprehensive multi-link dynamics tests with external oracles — 13-04
6. CI pipeline (GitHub Actions, 4-matrix, coverage) — 13-05
7. Code quality (namespace, OpenMP, umbrella header, stubs, GTest fallback) — 13-05
8. Formal conventions document — 13-06
9. Gravity invariant/property tests — 13-06
10. Edge case tests with explicit zero-mass contract — 13-07
11. Release-mode stability tests — 13-07
12. README update with gravity API — 13-07

---

## Per-Plan Scores

| Plan | Tasks | Files | Wave | Score | Notes |
|------|-------|-------|------|-------|-------|
| 13-01 | 3 | 8 | 1A | **PASS** | Well-scoped; grep-based call-site validation before/after removal |
| 13-02 | 2 | 4 | 1A | **PASS** | Frame convention documented; backward compatibility verified |
| 13-03 | 2 | 6 | 1B | **PASS** | Creates spatial_assert.h; replaces mixed pattern with uniform macros |
| 13-04 | 3 | 4 | 2 | **FLAG** | WARNING: Template-heavy test code in action; boundary with 13-06 may still blur |
| 13-05 | 3 | 6 | 3 | **FLAG** | WARNING: Still broad despite Group A/B split; compile_smoke_test.cpp not in files_modified |
| 13-06 | 3 | 7 | 2b | **PASS** | Clear boundary with 13-04; convention lock and worked examples |
| 13-07 | 3 | 4 | 3 | **FLAG** | WARNING: Zero-mass ABA test uses std::runtime_error — verify code actually throws this type |

**Overall: FLAG** — Plans are fundamentally sound. 4 warnings, 2 infos, 0 blockers.

---

## Dimension 1: Requirement Coverage

**Status: ✅ PASS**

| Requirement | Plans Covering | Status |
|-------------|---------------|--------|
| VEC-01 (Spatial Vector) | 13-01, 13-03, 13-05, 13-07 | ✅ Covered |
| UTL-03 (Spatial Utils) | 13-01, 13-03, 13-06, 13-07 | ✅ Covered |
| INR-01 (Inertia) | 13-01 | ✅ Covered |
| PLX-04 (Plücker Transform) | 13-01 | ✅ Covered |
| ABA-01 (Forward Dynamics) | 13-02, 13-04, 13-06, 13-07 | ✅ Covered |
| ABA-02 (Forward Dynamics) | 13-04, 13-06, 13-07 | ✅ Covered |
| TST-07 (Integration Tests) | 13-04, 13-06, 13-07 | ✅ Covered |

All requirement IDs from ROADMAP.md that are relevant to Phase 13 appear in at least one plan's `requirements` field. The phase does not introduce new requirement IDs — it extends existing ones.

**No PROJECT.md relevant requirements silently dropped.** Checked against v1.1 milestone docs and all requirements mapped to Phase 13 are accounted for.

---

## Dimension 2: Task Completeness

**Status: ✅ PASS**

All 19 tasks across 7 plans have:
- `<type>` — all `auto` ✅
- `<files>` — all present ✅
- `<read_first>` — all present (Task 1 of 13-01 has valid "discovery-only" explanation) ✅
- `<action>` — all specific with file:line references and code snippets ✅
- `<verify>` — all have `<automated>` sub-element ✅
- `<acceptance_criteria>` — all present with measurable conditions ✅

**No missing verify/automated elements.** Every task has a runnable verification command. Build-and-test commands are used for implementation tasks; grep commands are used for validation tasks.

Some tasks have `<manual-check>` or `<human-check>` as supplementary verify steps (13-03 Task 1, 13-04 Task 2, 13-05 Tasks 2-3, 13-06 Task 1, 13-07 Task 3). These are documented as additional verification on top of automated commands, not replacements — acceptable practice.

---

## Dimension 3: Dependency Correctness

**Status: ✅ PASS**

```
Wave 1A ──┬── 13-01 (no deps)
           └── 13-02 (no deps)
Wave 1B ──┬── 13-03 (depends_on: [13-01])
Wave 2  ──┴── 13-04 (depends_on: [13-01, 13-02, 13-03])
Wave 2b ──┴── 13-06 (depends_on: [13-04])
Wave 3  ──┬── 13-05 (depends_on: [13-01, 13-04])
           └── 13-07 (depends_on: [13-06])
```

- **No cycles detected** ✅
- **All referenced plans exist** ✅
- **Wave ordering respects dependency depth** ✅
- **13-03 (helper fixes) runs before 13-04 (tests that depend on correct helpers)** ✅ — addresses cross-plan HIGH concern
- **13-04 runs before 13-06** ✅ — ensures numerical regression tests exist before invariant tests

**Note:** Wave labels use mixed case (1A, 2b, 2, 3). Wave 2b (lowercase 'b') vs 1A (uppercase 'A'). This is cosmetic — `depends_on` determines execution order, wave labels are hints. If the orchestrator does case-sensitive wave sorting, 2b and 2 may sort differently than intended. Suggest keeping wave labels case-consistent (all uppercase: 1A, 1B, 2, 2A, 3).

---

## Dimension 4: Key Links Planned

**Status: ✅ PASS**

Key links across all plans connect dependent artifacts with explicit patterns:

| Plan | Link | Pattern | Status |
|------|------|---------|--------|
| 13-01 | `SpatialVector.cpp::crossForce` → `SpatialUtils.h::cross()` | `cross(ForceVector(*this), ForceVector(other))` | ✅ Matches action |
| 13-01 | `ForceVector.cpp::crossForce` → `SpatialUtils.h::cross()` | `cross(*this, other)` | ✅ Matches action |
| 13-01 | `ArticulatedBodyInertia.h::operator+` → `skew()` | `skew(other.getCom())` | ✅ Matches action |
| 13-02 | `ForwardDynamics::outwardPass` → `gravity member` | `c = MotionVector(..., -gravity)` | ✅ Matches action |
| 13-02 | `InverseDynamics::outwardPass` → `gravity member` | `- MotionVector(...)` | ✅ Matches action |
| 13-03 | `spatial_assert.h` → All NaN call sites | `#include "spatial_assert.h"` | ✅ |
| 13-05 | `CMakeLists.txt` → `compile_smoke_test.cpp` | `add_executable(SpatialAlgebraCompileSmoke ...)` | ✅ |
| 13-06 | `MATHEMATICAL_CONVENTIONS.md` → Test files | Cross-references to `SingleLinkGravityProportionality` | ✅ |
| 13-07 | EdgeCase tests → `MATHEMATICAL_CONVENTIONS.md` | `degenerate` contract reference | ✅ |

**No missing wiring detected.** Each plan's key links correspond to concrete action statements.

---

## Dimension 5: Scope Sanity

**Status: ⚠️ FLAG (2 warnings)**

### Metrics

| Plan | Tasks | Files | Status |
|------|-------|-------|--------|
| 13-01 | 3 | 8 | ✅ OK |
| 13-02 | 2 | 4 | ✅ OK |
| 13-03 | 2 | 6 | ✅ OK |
| 13-04 | 3 | 4 | ✅ OK (large test code volume — acceptable for test plan) |
| 13-05 | 3 | 6 | ⚠️ **WARNING** |
| 13-06 | 3 | 7 | ✅ OK (document + tests + header cross-links) |
| 13-07 | 3 | 4 | ⚠️ **WARNING** |

**WARNING (13-05):** Plan acknowledges review concern about being "oversized" and splits work into Group A (code hygiene) and Group B (build/CI infrastructure) with explicit acceptance checkpoints. However, Task 2 alone packs 4 distinct CMake changes (GTest fallback, stub removal, coverage option, compile smoke test + .gitignore). If any of these breaks unexpectedly, the executor must backtrack across unrelated concerns in the same task. Consider splitting into two plans: 13-05 (code hygiene: namespace, OpenMP, umbrella header, stubs) and a new 13-05b (build/CI: CMake changes, CI workflow, coverage).

**WARNING (13-07):** Task 3 (README update) is coupled with Task 1-2 (edge case tests). README updates belong in a separate wave since they touch docs not tests. Currently 13-07 is Wave 3, and the README work is in the same plan as test changes. If README review reveals issues, the test changes are held up too. Consider splitting README into a separate plan or making it the 4th task of 13-06.

### Context Budget Assessment
- Total: 19 tasks across 7 plans — reasonable (avg 2.7/plan)
- No plan exceeds 3 tasks ✅
- No plan exceeds 15 files ✅
- Highest file count: 13-01 (8 files) — within bounds ✅
- Estimated total context: ~45% — leaves room for execution ✅

---

## Dimension 6: Verification Derivation (must_haves)

**Status: ✅ PASS**

All 7 plans have `must_haves` with `truths`, `artifacts`, and `key_links`.

**Sample truth quality check:**

| Plan | Truth | User-Observable? | Verifiable? |
|------|-------|-------------------|-------------|
| 13-01 | "Force×force cross product is correct and consistent across all call sites" | ✅ | ✅ (tests + delegation pattern) |
| 13-01 | "PluckerTransform::apply() returns explicit types, not auto-deduced" | ✅ | ✅ (grep) |
| 13-02 | "Gravity is optional: default zero-vector preserves backward compatibility" | ✅ | ✅ (tests pass unchanged) |
| 13-03 | "NaN/Inf detection uses the same assertion macro everywhere" | ⚠️ Implementation-focused | ✅ (grep for uniformity) |
| 13-04 | "Numerical regression tests use external oracles" | ✅ | ✅ (grep for oracle reference) |
| 13-05 | "Library builds without system-installed GTest via FetchContent fallback" | ✅ | ⚠️ Partial (only testable in GTest-free env) |
| 13-06 | "Library has a formal conventions document with a 'convention lock'" | ✅ | ✅ (file exists + grep sections) |
| 13-07 | "Zero-mass/inertia edge case behavior is covered by an explicit contract" | ✅ | ✅ (test throws/nothrows) |

**Minor note:** 13-05 truth about "Library builds without system-installed GTest" is only partially verifiable on the development machine (where GTest is installed). The true test will happen on CI or a GTest-free environment. The plan acknowledges this with the `QUIET` mode gating.

---

## Dimension 7: Context Compliance (CONTEXT.md)

**Status: ✅ PASS — All 28 locked decisions addressed**

### Decision Coverage Map

| Decision | Description | Plan Task | Status |
|----------|-------------|-----------|--------|
| D-01 | Unify force×force cross product | 13-01 Task 2 | ✅ |
| D-02 | Correct formula (τ1×f2 − τ2×f1) | 13-01 Task 2 | ✅ |
| D-03 | Mixed torque+force tests | 13-03 Task 2 (Part C) | ✅ |
| D-04 | Remove MotionVector::crossForce, ForceVector::crossMotion | 13-01 Task 2 | ✅ |
| D-05 | Two-phase ABA inwardPass | Already done (Phase 12) | ✅ |
| D-06 | Known numerical values for validation | 13-04 Tasks 1-3 | ✅ |
| D-07 | Optional gravity parameter | 13-02 Tasks 1-2 | ✅ |
| D-08 | Gravity via base acceleration | 13-02 Tasks 1-2 | ✅ |
| D-09 | GitHub Actions single workflow | 13-05 Task 3 | ✅ |
| D-10 | Matrix build (ubuntu/macos × g++/clang++) | 13-05 Task 3 | ✅ |
| D-11 | Coverage tracking | 13-05 Task 2 | ✅ |
| D-12 | Tests on every push/PR | 13-05 Task 3 | ✅ |
| D-13 | Debug-mode NaN assertions | 13-03 Task 2 | ✅ |
| D-14 | No production overhead | 13-03 Task 1 (release no-op stubs) | ✅ |
| D-15 | Fix test helpers | 13-03 Task 2 (Part B) | ✅ |
| D-16 | RNEA non-zero velocity tests | 13-04 Task 2 | ✅ |
| D-17 | Multi-link numerical validation | 13-04 Task 1 | ✅ |
| D-18 | Remove Eigen version pin | Already done (PATTERNS.md says unpinned) | ✅ |
| D-19 | Fix bugs first, then tests, then CI | Wave structure (1A → 1B → 2 → 3) | ✅ |
| D-20 | ABI constructor arguments | 13-01 Task 3 | ✅ |
| D-21 | Plücker auto return type | 13-01 Task 3 | ✅ |
| D-22 | SpatialOperations unsafe downcasts | No explicit task — PATTERNS.md says already correct | ⚠️ INFO |
| D-23 | Global namespace pollution | 13-05 Task 1 | ✅ |
| D-24 | Include guard consistency | 13-05 Task 1 (PATTERNS.md says already correct) | ✅ |
| D-25 | GTest FetchContent fallback | 13-05 Task 2 | ✅ |
| D-26 | Remove OpenMP from LowerTriangular | 13-05 Task 1 | ✅ |
| D-27 | Remove empty stub .cpp files | 13-05 Task 2 | ✅ |
| D-28 | Umbrella header | 13-05 Task 1 | ✅ |

### Deferred Ideas Check
No tasks implement anything from the Deferred Ideas section (Python bindings, joint limits, branching tree ABA, LowerTriangular threshold, condition number estimation). ✅

### Discretion Areas
Discretion items (test values, cleanup order, namespace host header) are handled appropriately across plans. ✅

---

## Dimension 7b: Scope Reduction Detection

**Status: ✅ PASS — No scope reduction detected**

Scanned all plan action and acceptance criteria for scope reduction language:
- No "v1/v2" versioning patterns
- No "future enhancement" or "will be wired later" deferrals
- No "simplified for now" or "basic version" reductions
- No time-estimate-based scope justification
- No "skip for now" or "stub" patterns

All locked decisions are implemented at full scope as specified in CONTEXT.md. Where a decision references a contract boundary (zero-mass "REJECTED"), it is explicitly documented as such rather than silently reduced.

---

## Dimension 7c: Architectural Tier Compliance

**Status: ✅ PASS**

RESEARCH.md has `## Architectural Responsibility Map` at line 15 with 9 capability-to-tier mappings.

Cross-referencing plan tasks against the responsibility map:

| Capability | Expected Tier | Plan Assigning | Actual Tier | Status |
|------------|---------------|----------------|-------------|--------|
| Cross-product correctness | Library (SpatialUtils.h) | 13-01 | Library (SpatialUtils.h delegation) | ✅ Match |
| ABA inward pass | Library (ForwardDynamics.cpp) | 13-02 | Library (ForwardDynamics.cpp) | ✅ Match |
| Gravity term | Library (Forward/Inverse Dynamics) | 13-02 | Library (Forward/Inverse Dynamics) | ✅ Match |
| NaN/Inf guards | Library (core headers) | 13-03 | Library (spatial_assert.h + headers) | ✅ Match |
| Test helper fixes | Tests | 13-03 | Tests (TestSpatialOperations.cpp) | ✅ Match |
| Multi-link tests | Tests | 13-04, 13-06, 13-07 | Tests (test .cpp files) | ✅ Match |
| CI pipeline | Infrastructure (.github/) | 13-05 | Infrastructure (.github/workflows/ci.yml) | ✅ Match |
| Eigen 5.x compat | Build (CMakeLists.txt) | 13-05 | Build (CMakeLists.txt + CI job) | ✅ Match |
| Code quality fixes | Library headers/source | 13-05 | Library headers/source | ✅ Match |

**No tier mismatches detected.** Security-sensitive capabilities (cross-product correctness, NaN guards) are correctly placed in the Library tier, not in less-trusted tiers.

---

## Dimension 8: Nyquist Compliance

**Status: ⚠️ FLAG (1 warning)**

### Check 8e — VALIDATION.md Existence
VALIDATION.md exists at `13-VALIDATION.md`. ✅

### Check 8a — Automated Verify Presence
- All 19 tasks across 7 plans have `<verify>` with `<automated>` command ✅
- No `<automated>MISSING</automated>` references found ✅

### Check 8b — Feedback Latency Assessment
- Primary verify commands: `cmake --build build && ctest --output-on-failure` — estimated 30-60s ✅
- Filtered test commands: `ctest -R "TestForwardDynamics|TestInverseDynamics"` — faster (~10-15s) ✅
- No watch-mode flags (`--watchAll`) detected ✅
- Grep-based verify commands: near-instantaneous ✅

### Check 8c — Sampling Continuity
| Wave | Plans | Tasks | Automated Verify | Window Check |
|------|-------|-------|-----------------|--------------|
| 1A | 13-01, 13-02 | 5 | 5/5 | ✅ ≥2/3 in any window |
| 1B | 13-03 | 2 | 2/2 | ✅ (only 2 tasks) |
| 2 | 13-04 | 3 | 3/3 | ✅ |
| 2b | 13-06 | 3 | 3/3 | ✅ |
| 3 | 13-05, 13-07 | 6 | 6/6 | ✅ |

No 3 consecutive tasks without automated verify. ✅

### Check 8d — Wave 0 Completeness
No `MISSING` references found in any automated verify command. ✅

### VALIDATION.md Metadata
The VALIDATION.md file has:
- `nyquist_compliant: false` — accurate (stub table not filled)
- `wave_0_complete: false` — claims Wave 0 not complete
- Task-to-validation mapping table uses "TBD" placeholders

**WARNING:** The VALIDATION.md file has not been updated from its initial template state. The task IDs in the validation table are all "TBD" — none map to the actual 19 tasks across 7 plans. While the plans themselves define proper `<automated>` verify commands, the VALIDATION.md document is misleadingly incomplete. This should be updated by the planner to reflect the actual task structure before execution, or updated by the executor as part of execution.

---

## Dimension 9: Cross-Plan Data Contracts

**Status: ✅ PASS**

### Shared File Analysis

| File | Modified By | Waves | Conflict Risk |
|------|-------------|-------|--------------|
| `include/SpatialUtils.h` | 13-01 (Wave 1A), 13-03 (Wave 1B) | Sequential | ✅ No conflict (13-03 after 13-01) |
| `include/ArticulatedBodyInertia.h` | 13-01 (Wave 1A), 13-03 (Wave 1B) | Sequential | ✅ No conflict (fix then add guards) |
| `src/SpatialVector.cpp` | 13-01 (Wave 1A), 13-03 (Wave 1B) | Sequential | ✅ No conflict (delegate then add guards) |
| `include/ForceVector.h` | 13-01 (Wave 1A) | Single | ✅ |
| `tests/TestForwardDynamics.cpp` | 13-04 (Wave 2), 13-06 (Wave 2b), 13-07 (Wave 3) | Sequential | ✅ Append-only additions |
| `tests/TestInverseDynamics.cpp` | 13-04 (Wave 2), 13-06 (Wave 2b), 13-07 (Wave 3) | Sequential | ✅ Append-only additions |
| `tests/TestSpatialOperations.cpp` | 13-03 (Wave 1B), 13-07 (Wave 3) | Sequential | ✅ Append-only additions |
| `tests/TestDynamicsConsistency.cpp` | 13-04 (Wave 2), 13-06 (Wave 2b) | Sequential | ✅ Append-only additions |

**No same-wave file conflicts** — all shared files are modified in different waves with dependency ordering.

**No conflicting data transformations** — all test modifications are append-only (adding new TEST() blocks). No plan strips or transforms data that another plan needs in original form.

---

## Dimension 10: AGENTS.md Compliance

**Status: ✅ PASS**

| AGENTS.md Directive | Plan Compliance | Status |
|--------------------|----------------|--------|
| Build: `cmake -B build && cmake --build build` | All plans use these commands | ✅ |
| Test: `cd build && ctest --output-on-failure` | All test verify commands use this | ✅ |
| Doxygen `@brief`/`@details` on declarations | Plans add/update Doxygen (13-02, 13-06) | ✅ |
| Include guards: `#ifndef`/`#define` | spatial_assert.h: `#ifndef SPATIAL_ASSERT_H`; SpatialAlgebra.h: `#ifndef SPATIAL_ALGEBRA_H` | ✅ |
| Type aliases at file level | Plans follow existing pattern | ✅ |
| Eigen requires `find_package(Eigen3 ...)` | CMakeLists.txt changes preserve this | ✅ |
| Test executables registered in CMakeLists.txt | Plans add new test executables via `add_executable()` + `add_test()` | ✅ |

**No forbidden patterns detected.** Plans introduce no patterns that AGENTS.md forbids. Plans do not skip steps AGENTS.md requires.

---

## Dimension 11: Research Resolution

**Status: ✅ PASS**

RESEARCH.md at `13-RESEARCH.md`:
- `## Open Questions (RESOLVED)` at line 481 ✅
- All 3 questions resolved:
  1. Featherstone reference values → "derive from round-trip consistency" addressed in 13-04 ✅
  2. Vector3d host header → "SpatialVector.h" addressed in 13-05 ✅
  3. CodeCov upload method → "codecov-action@v4" addressed in 13-05 ✅

No unresolved open questions remain. Research confidence is HIGH per metadata.

---

## Dimension 12: Pattern Compliance

**Status: ✅ PASS**

PATTERNS.md at `13-PATTERNS.md` (1020 lines) maps 31 file analogs.

### Analog Reference Check

| Plan | Files Modified | Pattern Section Referenced | Status |
|------|---------------|---------------------------|--------|
| 13-01 | SpatialUtils.h, SpatialVector.cpp, ForceVector.cpp, ForceVector.h, MotionVector.cpp, MotionVector.h, ArticulatedBodyInertia.h, PluckerTransform.h | `src/SpatialVector.cpp`, `src/ForceVector.cpp`, `include/ArticulatedBodyInertia.h`, `Cross-Product Delegation Pattern` | ✅ |
| 13-02 | ForwardDynamics.h/cpp, InverseDynamics.h/cpp | `include/ForwardDynamics.h`, `src/ForwardDynamics.cpp`, `include/InverseDynamics.h`, `src/InverseDynamics.cpp` | ✅ |
| 13-03 | spatial_assert.h, SpatialUtils.h, SpatialVector.cpp, RigidBodyInertia.h, ArticulatedBodyInertia.h, TestSpatialOperations.cpp | `Debug-mode NaN/Inf Assertion Pattern`, `include/SpatialVector.h`, `include/RigidBodyInertia.h` | ✅ |
| 13-05 | SpatialVector.h, LowerTriangular.h, SpatialAlgebra.h, CMakeLists.txt, .github/workflows/ci.yml | `Shared Patterns > Include Guard Convention`, `CMakeLists.txt`, `.github/workflows/ci.yml` | ✅ |
| 13-06 | MATHEMATICAL_CONVENTIONS.md (new), test files | `GTest Test Suite Pattern` | ✅ |
| 13-07 | test files, README.md | `GTest Test Suite Pattern` | ✅ |

### Shared Patterns Coverage
- Include Guard Convention: SpatialAlgebra.h uses `#ifndef SPATIAL_ALGEBRA_H` ✅
- Doxygen Documentation Convention: Plan 13-02 adds Doxygen; Plan 13-06 adds cross-links ✅
- Debug-mode NaN/Inf Assertion Pattern: Created as spatial_assert.h macros ✅
- Exception-based Error Handling: 13-07 test expects `std::runtime_error` ✅
- GTest Test Suite Pattern: All test additions follow the pattern ✅

### No Analog Found Files
- `.github/workflows/ci.yml` — Plan 13-05 uses RESEARCH.md patterns correctly ✅
- `include/SpatialAlgebra.h` — Plan 13-05 uses PATTERNS.md umbrella header pattern ✅

---

## 10 HIGH Review Concerns from 13-REVIEWS.md — Resolution Status

| # | Concern | Plan | Resolution | Status |
|---|---------|------|------------|--------|
| 1 | Removed methods breaking change | 13-01 | Pre-removal grep scan (Task 1) + post-removal validation (Task 3) | ✅ Resolved |
| 2 | Sign convention errors in gravity | 13-02 | World-frame convention explicitly documented; g=0 backward compat verified | ✅ Resolved |
| 3 | Frame ambiguity under-specified | 13-02 | `@note` specifying world/base coordinates in both ForwardDynamics.h and InverseDynamics.h | ✅ Resolved |
| 4 | Test helper fixes prerequisite | 13-03 | Wave 1B (13-03) before Wave 2 (13-04); depends_on captures this | ✅ Resolved |
| 5 | Round-trip not independent oracle | 13-04 | Tests labeled as "self-consistency only"; external oracles (hand-derived, Python) provided | ✅ Resolved |
| 6 | Numerical validation underspecified | 13-04 | Hand-derived single-link oracle (1.0 rad/s²) + Python RNEA cross-validation | ✅ Resolved |
| 7 | Plan 13-05 oversized | 13-05 | Split into Group A/B with acceptance checkpoints; still broad | ⚠️ Addressed but not fully solved |
| 8 | FetchContent conflicts with Homebrew | 13-05 | `find_package(GTest QUIET)` → system install takes priority; FetchContent is fallback | ✅ Resolved |
| 9 | Zero-mass contract ambiguous | 13-07 | Explicit contract: rejected (throws) vs supported (finite output for near-zero) | ✅ Resolved |
| 10 | 13-04/13-06 boundary needed | 13-04, 13-06 | Both plans explicitly document boundary in objectives | ✅ Resolved |

---

## Issues Summary

### Blockers: 0
None — all dimensions pass at the blocking level.

### Warnings: 4

**W-1 [scope_sanity] Plan 13-05 still broad despite group split**
- **File:** `.planning/phases/13-production-readiness/13-05-PLAN.md`
- **Finding:** Plan acknowledges HIGH concern #7 about being "oversized" and splits work into Group A (Task 1) and Group B (Tasks 2-3), but the breadth of changes remains high: namespace cleanup + OpenMP + umbrella header (Task 1) combined with GTest fallback + stub removal + coverage option + compile smoke test + .gitignore (Task 2) + CI workflow (Task 3) — 3 tasks covering 7 distinct engineering concerns. If the CI workflow YAML needs review or the FetchContent approach breaks, it blocks the unrelated namespace cleanup from being verified independently.
- **Fix hint:** Split into two plans: 13-05a (code hygiene: namespace, OpenMP, umbrella header, stubs) and 13-05b (build/CI: CMake changes, CI workflow, coverage). Both at Wave 3 since neither depends on 13-06/13-07.

**W-2 [scope_sanity] Plan 13-07 couples README updates with edge case tests**
- **File:** `.planning/phases/13-production-readiness/13-07-PLAN.md`
- **Finding:** Task 3 (README.md update with gravity example, migration guide, API changes) shares a plan with Tasks 1-2 (edge case tests, release-mode stability). If README review reveals issues, the test work is blocked from verification. README updates belong in a documentation-focused plan or as a separate task in the conventions plan.
- **Fix hint:** Move README update to a separate Wave 4 plan, or add it as Task 4 of 13-06 (conventions + documentation).

**W-3 [nyquist_compliance] VALIDATION.md has TBD placeholders and stub status**
- **File:** `.planning/phases/13-production-readiness/13-VALIDATION.md`
- **Finding:** The validation mapping table still uses "TBD" as task IDs for all entries. `nyquist_compliant: false` and `wave_0_complete: false`. While the plans themselves define proper `<automated>` verify commands, the VALIDATION.md has not been updated to reflect the actual task-to-validation mapping. This document is misleading if read in isolation.
- **Fix hint:** Either update VALIDATION.md with actual task IDs and automated commands (mapping each of the 19 tasks to its verify command), or set `nyquist_compliant: true` after confirming all tasks have automated verify and set `wave_0_complete: true` after confirming no MISSING references.

**W-4 [dependency_correctness] Wave label case inconsistency**
- **File:** Multiple PLAN.md files
- **Finding:** Waves use mixed case: 13-01 `wave: 1A`, 13-02 `wave: 1A`, 13-03 `wave: 1B`, 13-06 `wave: 2b` (lowercase 'b'). If the orchestrator does case-sensitive alphabetical sorting, `2b` may sort differently from `2` and `1B`. The dependency graph (`depends_on`) is the true ordering mechanism, but inconsistent wave labels create confusion and potential execution ordering bugs.
- **Fix hint:** Standardize all wave labels: use `1A`, `1B`, `2`, `2A` (not `2b`), `3`.

### Infos: 2

**I-1 [context_compliance] D-22 not explicitly verified in any plan**
- **File:** `.planning/phases/13-production-readiness/13-01-PLAN.md` through `13-07-PLAN.md`
- **Finding:** D-22 (SpatialOperations unsafe downcasts — change signatures to accept `const MotionVector&` and `const ForceVector&` directly) is a locked decision from CONTEXT.md. PATTERNS.md says the signatures "already accept the correct types" and "no change needed." However, no plan task explicitly verifies this or mentions D-22. While likely already done, a locked decision with zero verification creates risk that it was missed or partially done.
- **Fix hint:** Add a grep verification to any plan: `grep "static_cast.*SpatialVector" include/SpatialOperations.h src/SpatialOperations.cpp` should return zero results. Plan 13-05 (code quality) is a natural home for this.

**I-2 [task_completeness] compile_smoke_test.cpp not listed in 13-05 files_modified**
- **File:** `.planning/phases/13-production-readiness/13-05-PLAN.md`
- **Finding:** Plan 13-05 Task 2 creates `tests/compile_smoke_test.cpp` in its action, but this file is not listed in the plan's `files_modified` frontmatter. This is a metadata gap — the plan will create a file not declared in its manifest. While the acceptance criteria verify its existence, the missing metadata means the orchestrator's file-change tracking won't know about it.
- **Fix hint:** Add `- tests/compile_smoke_test.cpp` to plan 13-05's `files_modified` list.

---

## Recommendations

1. **No blockers — execution can proceed** with the noted warnings and infos.
2. Consider splitting 13-05 into two plans (code hygiene + build/CI) for better isolation (W-1).
3. Consider moving README updates to a separate plan or merging with 13-06 (W-2).
4. Update VALIDATION.md with actual task IDs before marking phase complete (W-3).
5. Standardize wave label casing (W-4).
6. Add D-22 grep verification to 13-05 (I-1).
7. Add `tests/compile_smoke_test.cpp` to 13-05's files_modified (I-2).

---

## Summary

| Dimension | Status |
|-----------|--------|
| 1. Requirement Coverage | ✅ PASS |
| 2. Task Completeness | ✅ PASS |
| 3. Dependency Correctness | ✅ PASS (1 info) |
| 4. Key Links Planned | ✅ PASS |
| 5. Scope Sanity | ⚠️ FLAG (2 warnings) |
| 6. Verification Derivation | ✅ PASS |
| 7. Context Compliance | ✅ PASS (1 info) |
| 7b. Scope Reduction Detection | ✅ PASS |
| 7c. Architectural Tier Compliance | ✅ PASS |
| 8. Nyquist Compliance | ⚠️ FLAG (1 warning) |
| 9. Cross-Plan Data Contracts | ✅ PASS |
| 10. AGENTS.md Compliance | ✅ PASS |
| 11. Research Resolution | ✅ PASS |
| 12. Pattern Compliance | ✅ PASS |

**Plans are fundamentally sound and WILL achieve the phase goal.** 0 blockers, 4 warnings, 2 infos. The 10 HIGH review concerns from 13-REVIEWS.md are addressed across the plans with concrete mitigations. Recommendations are optional improvements — the revision loop COULD iterate on the warnings but is not required given no blocking issues exist.

---

## ISSUES FOUND

**Phase:** 13-production-readiness
**Plans checked:** 7
**Issues:** 0 blocker(s), 4 warning(s), 2 info(s)

### Structured Issues

```yaml
issues:
  - plan: "13-05"
    dimension: scope_sanity
    severity: warning
    description: "Plan still broad despite Group A/B split — 3 tasks cover 7 distinct engineering concerns (namespace, OpenMP, umbrella header, GTest fallback, stub removal, coverage option, compile smoke test, .gitignore, CI workflow)"
    fix_hint: "Split into two plans: 13-05a (code hygiene) and 13-05b (build/CI infrastructure). Both at Wave 3."

  - plan: "13-07"
    dimension: scope_sanity
    severity: warning
    description: "READ.me update (Task 3) couples documentation with edge case tests (Tasks 1-2) — if README review reveals issues, test work is blocked"
    fix_hint: "Move README update to a separate Wave 4 plan or add as Task 4 of 13-06."

  - plan: null
    dimension: nyquist_compliance
    severity: warning
    description: "VALIDATION.md task-to-validation table uses TBD placeholders for all entries; nyquist_compliant: false and wave_0_complete: false despite plans having proper automated verify commands"
    file: "13-VALIDATION.md"
    fix_hint: "Update VALIDATION.md with actual task IDs and verify commands from plans, or mark nyquist_compliant: true after confirming coverage."

  - plan: null
    dimension: dependency_correctness
    severity: warning
    description: "Wave labels use mixed case (1A, 1B, 2b) — case-sensitive sorting may cause execution ordering bugs"
    fix_hint: "Standardize to 1A, 1B, 2, 2A, 3 (uppercase letters)."

  - plan: "13-05"
    dimension: context_compliance
    severity: info
    description: "D-22 (SpatialOperations downcast fix) not explicitly verified in any plan — PATTERNS.md says already correct but no task confirms"
    fix_hint: "Add grep verification to 13-05: check static_cast<SpatialVector> is not used in SpatialOperations."

  - plan: "13-05"
    dimension: task_completeness
    severity: info
    description: "tests/compile_smoke_test.cpp created by Task 2 action but not listed in files_modified frontmatter"
    fix_hint: "Add - tests/compile_smoke_test.cpp to 13-05 files_modified."
```

---

## CHECK COMPLETE

**Overall status: ISSUES FOUND (0 blockers, 4 warnings, 2 infos)**

Plans can execute with current structure. Warnings are actionable but not blocking. Revision loop may iterate at planner's discretion.
