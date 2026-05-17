# Phase 13: Production Readiness — Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-05-17
**Phase:** 13-production-readiness
**Areas discussed:** cross-product fix, ABA multi-link, gravity, CI, NaN/Inf guards, test helpers, Eigen 5.x, execution ordering, ABI args, auto return type, downcasts, namespace pollution, include guards, GTest fallback, OpenMP, stub files, umbrella header

---

## Cross-Product Fix Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Single canonical impl | Put correct formula in SpatialUtils.h, delegate from others | ✓ |
| Fix each in place | Fix all 3 independently | |

**User's choice:** Single canonical impl
**Notes:** Mixed torque+force test for verification, remove deprecated overloads.

## ABA Multi-Link Fix

| Option | Description | Selected |
|--------|-------------|----------|
| Two-phase restructure | Phase 1: initialize from RBI. Phase 2: accumulate without re-init | ✓ |
| Surgical fix | Move initialization outside tip→base loop | |

**User's choice:** Two-phase restructure
**Notes:** Verify with Featherstone textbook example numerical values.

## Gravity Term Design

| Option | Description | Selected |
|--------|-------------|----------|
| Parameter to compute methods | Optional gravity vector arg to computeAccelerations/computeTorques | ✓ |
| Set on solver object | add setGravity() method | |

**User's choice:** Parameter to compute methods
**Notes:** Implement via base link acceleration (a₀ = -g). Propagates through outward pass automatically.

## CI Pipeline Choice

| Option | Description | Selected |
|--------|-------------|----------|
| GitHub Actions | Ubuntu + macOS matrix, g++ + clang++ | ✓ |
| Pre-commit hooks only | Git hooks, no PR protection | |

**User's choice:** GitHub Actions
**Notes:** Add coverage tracking with gcov + CodeCov.

## NaN/Inf Guard Pattern

| Option | Description | Selected |
|--------|-------------|----------|
| Debug-mode assertions | assert() on constructors, arithmetic, inertia apply | ✓ |
| Production exceptions | throw std::invalid_argument always | |

**User's choice:** Debug-mode assertions
**Notes:** Covers constructors + arithmetic + inertia apply only, not every method.

## Test Helper Fix + Multi-Link Tests

| Option | Description | Selected |
|--------|-------------|----------|
| Fix in place | Correct createIdentityInertia, add inertia assertions | ✓ |
| Replace with inline constants | Remove helpers, use local constants | |

**User's choice:** Fix in place
**Notes:** RNEA tests with known qdot values and expected tau. Multi-link ABA with Featherstone reference values.

## Eigen 5.x Compatibility

| Option | Description | Selected |
|--------|-------------|----------|
| Remove version pin | find_package(Eigen3 REQUIRED NO_MODULE) | ✓ |
| Add version range | Explicit version range syntax | |

**User's choice:** Remove version pin

## Execution Ordering

| Option | Description | Selected |
|--------|-------------|----------|
| Bugs first, then tests, then CI | Natural dependency order | ✓ |
| CI first | Get test results visible first | |
| Parallel: CI + test helpers | Maximize parallelism | |

**User's choice:** Bugs first, then tests, then CI

## ABI Constructor Args

| Option | Description | Selected |
|--------|-------------|----------|
| Leave to agent discretion | Bug well-documented in CONCERNS.md | ✓ |
| Let me specify | User has specific preferences | |

**User's choice:** Agent discretion

## PluckerTransform auto Return Type

| Option | Description | Selected |
|--------|-------------|----------|
| Fix it | Change auto to PluckerTransform explicitly | ✓ |
| Leave as is | Works currently (only called from apply()) | |

**User's choice:** Fix it

## Cleanup Items (Downcasts, Namespace, Guards)

| Option | Description | Selected |
|--------|-------------|----------|
| Fix all | Signatures, namespace, include guards | ✓ |
| Fix only functional | Downcasts + namespace, skip guards | |
| Defer all | | |

**User's choice:** Fix all

## Batch Cleanup (GTest Fallback, OpenMP, Stubs, Umbrella)

| Option | Description | Selected |
|--------|-------------|----------|
| Fix all | FetchContent, remove OpenMP, remove stubs, add umbrella | ✓ |
| Fix only essential | GTest fallback + stub removal | |
| Defer all | | |

**User's choice:** Fix all

---

## the agent's Discretion

- ABI constructor argument fix implementation details
- Specific test values for multi-link reference tests
- Code cleanup order (auto return, downcasts, namespace, guards, GTest fallback, OpenMP, stubs, umbrella)
- Which header hosts the canonical `using Vector3d` declaration

## Deferred Ideas

- Python bindings (pybind11) — v2 milestone
- Joint limit / singularity handling — future phase
- Branching tree ABA with child index lists — future optimization
- LowerTriangular inverse threshold parameterization — future
- Condition number estimation — not currently needed

---

*Discussion log: 2026-05-17*
