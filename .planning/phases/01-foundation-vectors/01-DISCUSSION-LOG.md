# Phase 1: Foundation Vectors - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions captured in CONTEXT.md — this log preserves the discussion.

**Date:** 2026-05-15
**Phase:** 01-foundation-vectors
**Mode:** discuss
**Areas analyzed:** Test framework approach, Implementation completeness, API design patterns, Verification strategy

## Gray Areas Discussed

### Test Framework Approach
**Question:** Current tests use assert() but Phase 8 is dedicated to test infrastructure with GTest. Should Phase 1 use GTest now or defer?

**Options presented:**
1. Defer to Phase 8 — Keep simple assert() tests for Phase 1, migrate to GTest in Phase 8
2. Use GTest now (Recommended) — Upgrade to GTest now in Phase 1 for better test reporting and assertions
3. Hybrid approach — Use both — assert() for basic checks, add GTest for detailed verification

**Decision:** Use GTest now (Recommended)

**Rationale:** Better test reporting and assertions from the start, even though Phase 8 is dedicated to test infrastructure.

---

### Implementation Completeness
**Question:** MotionVector::crossMotion uses (ω1×ω2, v1×v2) but SpatialVector uses (ω1×ω2, ω1×v2+v1×ω2). Which is correct per Featherstone?

**Options presented:**
1. Fix to match base class (Recommended) — MotionVector should use the same formula as SpatialVector
2. Keep MotionVector version — SpatialVector base class is wrong — MotionVector has the correct specialized formula
3. Document the difference — Both are intentional — different operations for different use cases

**Decision:** Fix to match base class (Recommended)

**Rationale:** The base class SpatialVector has the correct Featherstone formula. MotionVector implementation is buggy.

---

### API Design Patterns
**Question:** Should Phase 1 add convenience methods (Zero(), Random(), factories) or keep the minimal constructor-focused API?

**Options presented:**
1. Keep minimal (Recommended) — Minimal API — only essential constructors and operations (current state)
2. Add convenience factories — Add factory methods like MotionVector::Zero(), MotionVector::Random(), fromScrewAxis()
3. Maximize ergonomics — Add implicit conversions and more operators for ergonomic usage

**Decision:** Keep minimal (Recommended)

**Rationale:** Focus on core functionality first. Convenience methods can be added in later phases if needed.

---

### Verification Strategy
**Question:** For verifying spatial vector operations, what testing approach should we use?

**Options presented:**
1. Textbook examples (Recommended) — Unit tests with known numerical values from Featherstone textbook examples
2. Property-based testing — Property-based tests — verify invariants like R·R^T = I for rotations
3. Combined approach — Both — textbook examples for correctness, properties for robustness

**Decision:** Combined approach

**Rationale:** Mathematical correctness requires both specific numerical verification and general property validation.

---

## Corrections Made

No corrections — all assumptions confirmed through discussion.

## External Research

No external research performed — decisions made from codebase analysis and user discussion.

---

*Discussion completed: 2026-05-15*
