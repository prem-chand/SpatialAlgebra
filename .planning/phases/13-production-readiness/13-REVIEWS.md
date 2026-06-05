---
phase: 13
reviewers: [codex]
reviewed_at: 2026-05-27T13:39:37Z
plans_reviewed: [13-01-PLAN.md, 13-02-PLAN.md, 13-03-PLAN.md, 13-04-PLAN.md, 13-05-PLAN.md, 13-06-PLAN.md, 13-07-PLAN.md]
---

# Cross-AI Plan Review — Phase 13

## Gemini Review

**UNAVAILABLE** — Gemini CLI was detected but not authenticated. Set `GEMINI_API_KEY` or configure `~/.gemini/settings.json` to enable Gemini review.

---

## Codex Review

### Overall Summary

The phase is generally well-shaped: the plans are split along sensible fault lines, the wave ordering mostly matches the true dependency graph, and the scope is aligned with a production-readiness objective rather than feature expansion. The strongest part is that mathematical correctness, testing, CI, and API stability are all represented. The main risks are not missing functionality so much as hidden coupling between plans 13-01/13-02 and the downstream test/documentation/CI work, plus the possibility that "code quality + CI" has become a catch-all bucket with enough breadth to destabilize the phase if not tightly controlled. If executed with strict interfaces and clear acceptance criteria per plan, this phase should achieve its stated goals.

---

### 13-01: Fix 3 Critical Bugs

**Summary:** This is the right place to start. It targets correctness defects at the semantic layer before downstream validation work depends on them. The plan appears appropriately focused on canonicalizing cross-product behavior, correcting a constructor bug, resolving return-type ambiguity, and removing invalid APIs. Those changes directly reduce undefined behavior and conceptual drift in a math library.

**Strengths:**
- Fixes root-cause correctness issues before expanding test coverage.
- Removes physically meaningless APIs instead of preserving misleading surface area.
- Canonical implementation choice for force×force reduces duplication and divergence risk.
- Explicit return types in Plücker code improve readability and compiler portability.

**Concerns:**
- **HIGH:** Removing deprecated methods may still be a breaking change if any internal tests, examples, or downstream code still use them.
- **MEDIUM:** "Single canonical implementation" in `SpatialUtils.h` can create accidental semantic mismatch if old call sites relied on subtle operand ordering or frame conventions.
- **MEDIUM:** ABI constructor arg swap is a classic silent-correctness bug; existing tests may be insufficient to prove no remaining call-site confusion.
- **LOW:** Explicit return type changes may expose previously hidden narrowing or copy behavior differences across compilers.

**Suggestions:**
- Add compile-fail or grep-based validation that removed methods are no longer referenced anywhere in tree.
- Define exact algebraic sign/frame convention in the plan acceptance criteria, not just "unify implementation."
- Add targeted regression tests for the ABI constructor with asymmetric test data so swapped arguments cannot accidentally pass.
- Require before/after API notes in changelog or README if public methods are removed.

**Risk Assessment:** MEDIUM — Correct and necessary, but touches foundational math semantics and public API surface.

---

### 13-02: Add Gravity Support

**Summary:** The gravity plan is well-conceived and mathematically grounded. Using base acceleration for both ABA and RNEA is the right formulation and preserves backward compatibility through default-zero gravity. The main risk is not the concept, but ensuring consistent sign conventions, frame assumptions, and exact parity between inverse and forward dynamics paths.

**Strengths:**
- Backward-compatible API via optional parameter defaulting to zero.
- Uses standard rigid-body dynamics formulation rather than ad hoc force injection.
- Symmetric treatment across ABA and RNEA reduces conceptual fragmentation.
- Gravity support is isolated enough to be testable independently.

**Concerns:**
- **HIGH:** Sign convention errors are easy here and may pass simple single-link tests while failing multi-link chains.
- **HIGH:** Frame ambiguity is under-specified: is gravity expressed in world/base coordinates only, or arbitrary local frame coordinates?
- **MEDIUM:** Default argument changes in headers can create subtle overload/resolution issues if there are existing wrappers or bindings.
- **MEDIUM:** If gravity is applied through base acceleration, plans must verify consistency with existing bias-force terms under non-zero joint velocities.
- **LOW:** Documentation risk if examples omit units and frame assumptions.

**Suggestions:**
- State explicitly in the plan that gravity is expressed in base/world coordinates and document transformation assumptions.
- Add single-link analytical oracle tests in addition to multi-link propagation tests.
- Verify exact equivalence between `g = 0` new path and pre-change dynamics outputs.
- Add tests combining gravity and non-zero velocity/Coriolis terms, not only static gravity cases.

**Risk Assessment:** MEDIUM-HIGH — Mathematically sound approach, but easy to get subtly wrong in ways that only show up in more complex chains.

---

### 13-03: NaN/Inf Guards + Test Helpers

**Summary:** This plan strengthens diagnostics and test trustworthiness. The debug-only guard approach is appropriate for a performance-sensitive numerical library, and fixing broken test helpers is necessary before relying on new tests. The risk is that this plan currently mixes two very different concerns: runtime defensive checks and test infrastructure correctness.

**Strengths:**
- Debug-only assertions avoid release overhead.
- Fixing zero-matrix helper bugs improves the credibility of all subsequent tests.
- Cross-force mixed-input tests directly support the earlier semantic fixes.
- Guards are well-targeted for production-readiness without changing runtime semantics.

**Concerns:**
- **HIGH:** If test helpers are currently wrong, prior test evidence may be unreliable; downstream plans depending on them should not proceed until helper fixes are proven.
- **MEDIUM:** NaN/Inf assertions can be too shallow if only added at API boundaries and not at critical intermediate computations.
- **MEDIUM:** Assertions may change behavior in debug builds enough to mask sequencing issues or make existing tests fail nondeterministically.
- **LOW:** Combining helper fixes and guard insertion in one plan may blur root cause if tests start failing.

**Suggestions:**
- Split acceptance criteria into two checkpoints: helper correctness first, then assertion rollout.
- Add explicit tests that assertions trigger in debug for constructed NaN/Inf inputs if the test framework/build setup allows it.
- Define which "core ops" must be guarded so coverage is not left subjective.
- Re-run all dynamics tests after helper fixes before attributing failures to guard changes.

**Risk Assessment:** MEDIUM — Valuable plan, but it affects the trust model of the entire test suite.

---

### 13-04: Comprehensive Dynamics Tests

**Summary:** This is the highest-value validation plan in the phase. It directly targets the project's remaining correctness risk: whether dynamics algorithms behave properly beyond trivial cases. The test categories are strong, especially multi-link ABA validation and RNEA↔ABA consistency. The main risk is oracle quality: if expected values come from the same flawed assumptions as the implementation, the tests can look comprehensive without being truly independent.

**Strengths:**
- Focuses on multi-link systems where real bugs typically emerge.
- Includes both numerical validation and round-trip consistency.
- Non-zero velocity tests cover Coriolis/centrifugal paths often missed by static tests.
- Gravity propagation tests directly validate the new API behavior.

**Concerns:**
- **HIGH:** Round-trip consistency alone is not an independent oracle; mutually wrong ABA and RNEA implementations can still agree.
- **HIGH:** "Numerical validation" is underspecified; if expected values are generated from the library itself, the tests are weak.
- **MEDIUM:** Multi-link tests for 3+ link chains may still miss branched trees, mixed inertias, or offset COM configurations.
- **MEDIUM:** Floating-point tolerances need explicit policy or CI/compiler variation may create flakiness.
- **LOW:** Large hand-derived expected values can become brittle and hard to maintain without helper structure.

**Suggestions:**
- Require at least one external oracle source per major scenario: hand-derived small systems, symbolic derivation, or cross-check against trusted reference implementation.
- Include asymmetric inertias, non-zero COM offsets, and non-axis-aligned gravity vectors.
- Define tolerance bands per quantity type, not one global epsilon.
- Add separate tests for chain topology growth: 1-link, 2-link, 3-link, and a deeper smoke case.
- Treat known CR-02 cases as explicit regression tests with issue references.

**Risk Assessment:** MEDIUM-HIGH — Essential plan, but only if test oracles are genuinely independent.

---

### 13-05: Code Quality + CI

**Summary:** This plan covers important production infrastructure, but it is the broadest and therefore the most likely to accumulate incidental complexity. CI, build matrix updates, dependency handling, and code cleanup all matter, but mixing them with namespace cleanup, OpenMP removal, umbrella header work, and CMake modernization makes this plan vulnerable to scope creep and hidden regressions.

**Strengths:**
- CI matrix directly supports portability and production-readiness claims.
- Eigen 5.x compatibility and GTest dependency cleanup address real build friction.
- Removing empty test stubs and build dead weight reduces maintenance overhead.
- Coverage integration is useful for preventing future blind spots.

**Concerns:**
- **HIGH:** This plan is doing too many different things: build system evolution, dependency strategy, CI, public header packaging, and code cleanup.
- **HIGH:** `FetchContent` for GTest may conflict with Homebrew/system installs or existing local developer workflows if not carefully gated.
- **MEDIUM:** Umbrella header addition can create include-order problems, symbol collisions, or circular dependencies if not tested widely.
- **MEDIUM:** Coverage flags often differ by compiler/platform; mixing them into a 4-matrix CI can create fragile config logic.
- **MEDIUM:** OpenMP removal is low-risk only if it is truly unused; otherwise performance or build semantics may change unexpectedly.
- **LOW:** Namespace cleanup can create noisy diffs that obscure more meaningful changes during review.

**Suggestions:**
- Split this into at least two logical acceptance buckets internally: build/CI infrastructure and codebase hygiene/public packaging.
- Make `FetchContent` optional or fallback-based rather than mandatory if system GTest is already available.
- Keep coverage to a single dedicated CI job instead of all matrix jobs.
- Add a public-header compile smoke test: include umbrella header alone in a tiny target.
- Gate Eigen 5 compatibility with one explicit CI job and one local documented path, not diffuse support assumptions.

**Risk Assessment:** HIGH — Important work, but it is the plan most likely to overrun or introduce avoidable instability.

---

### 13-06: Formal Conventions + Gravity Invariant Tests

**Summary:** This is a strong supporting plan because it couples documentation with invariant-based validation. In a math-heavy library, formalizing conventions is not optional. The gravity test oracle ideas are especially good because they validate behavior through physics properties rather than implementation mirroring. The main gap is that these invariants are necessary but not sufficient.

**Strengths:**
- Documents mathematical conventions explicitly, reducing future ambiguity.
- Invariant-based tests are more robust than pure golden-value tests for some scenarios.
- Gravity-specific oracles are well chosen for sanity checking physical behavior.
- Independent verifiability improves confidence in the gravity implementation.

**Concerns:**
- **MEDIUM:** Invariants like proportionality and mass-scaling can hold even when sign or frame conventions are wrong in some configurations.
- **MEDIUM:** Documentation can diverge from implementation if not linked to specific tests/examples.
- **LOW:** A standalone conventions doc may be overlooked unless surfaced from README/API docs.
- **LOW:** This plan may duplicate some of 13-04 unless responsibilities are sharply defined.

**Suggestions:**
- Define this plan's test scope as invariant/property tests, while 13-04 owns explicit numerical regression tests.
- Include worked examples in the conventions doc that correspond to test cases.
- Add one "convention lock" section covering frame notation, gravity sign, motion/force ordering, and transform direction.
- Cross-link the document from README and relevant headers.

**Risk Assessment:** LOW-MEDIUM — Good leverage, low implementation risk, but should not be mistaken for full dynamics validation.

---

### 13-07: Edge Cases + Documentation

**Summary:** This is a sensible closing plan for hardening and usability. Zero-mass/inertia cases, release-mode stability, and README updates all fit production readiness. The main question is policy: are zero-mass or singular inertia values supported physical edge cases, invalid inputs, or debug-only traps? Without that decision, tests can encode the wrong contract.

**Strengths:**
- Addresses robustness beyond the "happy path."
- Release-mode stability coverage complements the debug guard work.
- README gravity API updates help prevent misuse after signature changes.
- Good final-pass plan once correctness and CI are in place.

**Concerns:**
- **HIGH:** Zero-mass/inertia behavior is ambiguous unless the library contract explicitly defines whether these inputs are supported, tolerated, or rejected.
- **MEDIUM:** Release-mode stability tests can be weak if they only check "no crash" instead of output finiteness or documented behavior.
- **MEDIUM:** Edge-case numerical behavior may differ across compilers/platforms and become flaky without precise expectations.
- **LOW:** Documentation updates may lag behind actual API names/signatures if done late.

**Suggestions:**
- Add an explicit contract decision: invalid-input debug assert, runtime exception, or mathematically defined behavior.
- For release-mode tests, check finite outputs and invariant preservation, not only process survival.
- Include near-singular cases in addition to exact zero values.
- Update README with one minimal gravity example and one migration note for prior API users.

**Risk Assessment:** MEDIUM — Useful hardening plan, but only if edge-case semantics are defined precisely.

---

### Cross-Plan Concerns

**Strengths:**
- Wave structure is mostly correct: semantic fixes first, validation and infrastructure second.
- The phase covers correctness, diagnostics, tests, docs, and CI together.
- Backward compatibility is considered explicitly in gravity API design.
- Mathematical conventions are treated as first-class, which is appropriate for this domain.

**Concerns:**
- **HIGH:** 13-05 is oversized and risks becoming a "miscellaneous changes" plan that weakens review quality.
- **HIGH:** 13-04 and 13-06 need cleaner boundaries to avoid duplicated or non-independent testing.
- **HIGH:** The dependency graph likely understates that 13-03 helper fixes are a prerequisite for trusting 13-04/13-06 test results.
- **MEDIUM:** Removed methods in 13-01 may require synchronized updates in docs/examples/tests owned by 13-05 or 13-07.
- **MEDIUM:** Edge-case semantics in 13-07 should inform guard/assert behavior in 13-03; there is some reverse coupling there.
- **LOW:** Security concerns are limited for a local C++ math library, but CI supply-chain choices in 13-05 deserve normal dependency pinning and action version pinning.

**Suggestions:**
- Reframe dependencies as: Wave 1A (13-01, 13-02), Wave 1B (13-03 helper fixes), Wave 2 (13-04, 13-06), Wave 3 (13-05, 13-07).
- Narrow 13-05 to build/CI only, and move purely cosmetic/code hygiene items elsewhere if needed.
- Add explicit acceptance criteria for each plan with objective evidence: g=0 parity preserved, multi-link oracle independence, debug assertions covered, umbrella header compile-tested, edge-case contract documented.
- Require one regression ledger mapping each critical bug/decision to a specific test.

**Overall Risk Assessment:** MEDIUM — The phase design is fundamentally sound and should achieve the stated production-readiness goals, but it carries execution risk from one oversized infrastructure plan, some understated test dependencies, and the need for stronger oracle independence in dynamics validation.

---

## Consensus Summary

Only one reviewer (Codex CLI) participated in this review. Gemini CLI was unavailable due to missing authentication configuration.

### Agreed Strengths
- Wave structure (semantic fixes → validation → CI/doc) is fundamentally correct.
- Backward compatibility preserved via default-zero gravity parameter.
- Mathematical conventions treated as first-class concern.
- Proper focus on fixing root causes before expanding test coverage.

### Agreed Concerns
- Plan 13-05 (Code Quality + CI) is oversized — mixes build, CI, namespace, OpenMP, umbrella header, and stub removal into one plan.
- Plans 13-04 and 13-06 have weak boundary separation (both test dynamics, risk of duplication).
- Test oracle independence is not guaranteed — round-trip consistency can mask mutually-wrong implementations.
- Zero-mass edge case semantics are underspecified (supported, tolerated, or rejected).
- Gravity frame conventions need explicit documentation (world/base coordinates).
- Dependency graph understates that 13-03 helper fixes are a prerequisite for trusting downstream test results.

### Divergent Views
- Only one reviewer participated, so no divergent views to report.

---

*Reviewed by Codex CLI. Gemini CLI unavailable (auth not configured).*
