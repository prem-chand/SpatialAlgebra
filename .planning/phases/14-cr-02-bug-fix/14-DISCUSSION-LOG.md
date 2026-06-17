# Phase 14: CR-02 Bug Fix - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-17
**Phase:** 14-cr-02-bug-fix
**Areas discussed:** Test strategy

---

## Test Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Update context and replan | Update CONTEXT.md with the finding, add non-zero COM test requirements, then replan | |
| Run Phase 14 fix as-is | The fix is still valid, just the test baseline was wrong. Execute existing plan, verify with robot examples | |
| Update tests first | TDD: add non-zero COM multi-link tests that fail (red), then fix the code (green), then verify | ✓ |

**User's choice:** Update tests first (TDD approach)
**Notes:** User chose to add failing tests with non-zero COM before implementing the fix. This ensures the fix is validated against tests that actually exercise the bug (existing tests all use zero COM and pass with the buggy code).

---

## Key Finding (pre-discussion)

All 22 existing ABA/FD/consistency tests pass with the current buggy 3-phase code because every multi-link test sets COM to `Vector3d::Zero()`. With zero COM, `skew(com)*mass = 0`, the inertia is block-diagonal, and the Phase 3 correction double-counting bug is not triggered.

The v1.2 robot examples (Phase 18) use non-zero COM (UR5-derived parameters) and correctly demonstrate the bug on 3-link chains.

## the agent's Discretion

- Exact non-zero COM values for new test models
- ThreeLinkNumericalValidation: use round-trip ID→FD check instead of hardcoded expected values
- Doxygen wording specifics

## Deferred Ideas

None.
