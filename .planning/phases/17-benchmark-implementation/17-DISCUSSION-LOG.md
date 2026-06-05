# Phase 17: Benchmark Implementation - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-05
**Phase:** 17-benchmark-implementation
**Areas discussed:** ABA torque input, RNEA acceleration strategy, Microbenchmark scope, State regeneration frequency

---

## ABA Torque Input

| Option | Description | Selected |
|--------|-------------|----------|
| All-ones vector | Eigen::VectorXd::Ones(nDOF) — trivially simple, zero code, fully deterministic | |
| Random via RandomState | Extend RandomState with randomTorques() method — more realistic profile | ✓ |

**User's choice:** Random via RandomState with range [-10, 10] Nm
**Notes:** User chose against the researcher's recommendation (all-ones) for realistic torque profiles. Range [-10, 10] selected to provide wider excitation than default velocity range.

---

## RNEA Acceleration Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Random qddot via RandomState | Extend RandomState with randomAccelerations() — exercises all RNEA code paths | ✓ |
| Zero qddot | Zero vector — simpler but misses S·q̈ term | |

**User's choice:** Random qddot via RandomState with range [-5, 5] rad/s² (matching D-10)
**Notes:** Follows recommendation. [-5, 5] range consistent with velocity range (D-10).

---

## Microbenchmark Scope

| Option | Description | Selected |
|--------|-------------|----------|
| Algorithm-critical subset | transformMotion(mv), inverseTransformForce(fv), cross(mv,mv), cross(mv,fv) | |
| Full API surface | All 8 transform ops + all 4 cross-product variants + inertia transforms | ✓ |

**User's choice:** Full API surface including inertia transforms (20 operations total)
**Notes:** User chose full API for complete performance profiling despite researcher's recommendation for smaller subset.

---

## State Regeneration Frequency

| Option | Description | Selected |
|--------|-------------|----------|
| Per-iteration (pre-allocated) | Fresh random state every iteration using pre-allocated buffers + fill methods | ✓ |
| Once per DOF | Regenerate once per DOF value — lower variance but warm-state bias | |

**User's choice:** Per-iteration with pre-allocated buffers and zero-alloc fill methods
**Notes:** Follows recommendation. Requires adding fillPositions, fillVelocities, fillTorques, fillAccelerations to RandomState.

---

## Claude's Discretion

- Exact benchmark function signatures (within `void BM_*(benchmark::State&, int nDOF)` pattern)
- DOF iteration step size (1 or configurable) within the 1..20 range
- Google Benchmark `MinTime` / iterations configuration
- Specific `DoNotOptimize` / `ClobberMemory` placement in microbenchmark loops

## Deferred Ideas

None.
