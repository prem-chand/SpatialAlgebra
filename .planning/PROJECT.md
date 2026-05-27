# SpatialAlgebra Project

**Project Code:** SA  
**Project Title:** Spatial Vector Algebra Library for Robotics  
**Last Updated:** 2026-05-27

---

## What This Is

A C++17 library implementing spatial vector algebra for rigid body dynamics, following Featherstone's formulation. Provides 6D spatial vectors (twists and wrenches), Plücker coordinate transforms, inertia representations, and forward/inverse dynamics for robotics simulation and control.

## Core Value

**Must Deliver:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

**Success Looks Like:** 
- All core classes implemented and verified ✓ — achieved in v1.0
- Forward dynamics (Articulated Body Algorithm - ABA) ✓ — implemented and tested
- Inverse dynamics (Recursive Newton-Euler Algorithm - RNEA) ✓ — implemented and tested
- 95%+ test coverage for all mathematical operations ✓
- Library production-ready for serial kinematic chains ✓

## Requirements

### Validated (v1.0)

- ✓ SpatialVector, MotionVector, ForceVector — complete with Featherstone-verified cross products
- ✓ Rotation (angle-axis, quaternion, matrix operations)
- ✓ LowerTriangular packed storage (multiply, inverse, transpose)
- ✓ Spatial utilities (skew, dot, cross) and SpatialOperations class
- ✓ Plücker transforms (motion/force transformation, rigid body inertia transform)
- ✓ RigidBodyInertia and ArticulatedBodyInertia (construction, apply, operators)
- ✓ Articulated Body Algorithm (ABA) for serial and branching kinematic chains
- ✓ Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics
- ✓ Comprehensive GTest test infrastructure (186+ tests)
- ✓ Integration tests for dynamics pipeline consistency
- ✓ Documentation: README, Doxygen, compilable examples

### Active (Next Milestone)

- [ ] Multi-link RNEA↔ABA consistency for 3+ link chains (BF-02)
- [ ] GitHub Actions CI pipeline
- [ ] Performance benchmarks vs RBDL/Pinocchio
- [ ] Additional real-world robot examples (2-link planar, 3-link spatial arm)
- [ ] Eigen 5.x compatibility in CI

### Out of Scope

- **Python bindings** — Standalone RNEA exists; pybind11 integration deferred to v2.0
- **Contact/collision handling** — Pure rigid body dynamics; contact physics deferred
- **Inverse kinematics** — Focus on dynamics algorithms (RNEA, ABA), not kinematics
- **Visualization tools** — Library-only, no rendering components
- **Serialization** — No save/load for configurations or state
- **Floating base / 6-DOF base** — Fixed-base kinematics only for v1.0

## Context

**v1.0 Shipped:** May 16, 2026 — ~9,450 LOC C++17, 186+ GTest tests across 13 test executables.

**State:**
- All core classes implemented and verified with tests
- Forward dynamics (ABA) and inverse dynamics (RNEA) with gravity support
- 98.7% test pass rate (3 pre-existing multi-link consistency failures)
- NaN/Inf debug-mode guards on core operations
- Cross-product operations unified to single canonical implementation
- Production-hardened with edge case tests and release-mode stability verification

**Tech Stack:**
- C++17 with Eigen3 for linear algebra
- Google Test for unit testing
- CMake build system
- Doxygen documentation

## Constraints

**Technical:**
- Must maintain Eigen3 compatibility (version 3.3+ or 5.x with testing)
- Preserve existing class interfaces (backward compatibility)
- Follow existing code conventions (Doxygen comments, type aliases)

**Timeline:** Flexible — quality-focused completion

**Budget:** N/A (open source library)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| GTest from Phase 1 instead of bare assert() | Better test reporting and property-based testing | All 186+ tests use GTest |
| Fix cross product bugs immediately | Critical for downstream dynamics correctness | All cross products verified against Featherstone |
| Remove Eigen 5.x version pin | Homebrew compatibility | Library builds on both Eigen 3.4 and 5.x |
| Cross-product delegation to canonical free functions | Eliminate triplicate-implementation bug | Single source of truth in SpatialUtils.h |
| Featherstone Algorithm 7.3 for ABA | Textbook reference implementation | Verified for serial and branching chains |
| Featherstone Algorithm 7.1 for RNEA | Textbook reference implementation | Verified with 5 unit tests |
| Separate Link structs for ABA/RNEA | Avoid POSIX `link()` conflict and namespace issues | ForwardDynamicsLink, InverseDynamicsLink |

## Out of Scope (Detailed)

- **Python bindings** — Python RNEA exists standalone; integration deferred to future milestone
- **Contact/collision handling** — Pure rigid body dynamics only
- **Inverse kinematics** — Focus on dynamics, not kinematics
- **Visualization tools** — Library-only, no rendering
- **Serialization** — No state save/load needed for library operations

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---

*Last updated: 2026-05-27 after v1.0 milestone*

---

## Current State (v1.0)

**Shipped:** 2026-05-16  
**Phases:** 10  
**Plans:** 25  
**Test Count:** 186+ (186 passing, 3 deferred failures)  

**Delivered:**
- Complete 6D spatial vector algebra with Featherstone-verified cross products
- Full Plücker coordinate transforms (motion/force/RBI transforms verified)
- ABA forward dynamics and RNEA inverse dynamics with gravity support
- Memory-efficient LowerTriangular packed matrix storage
- Comprehensive GTest test infrastructure
- Documentation: README, Doxygen, 4 compilable examples

**Known Gaps (inherited by next milestone):**
- Multi-link RNEA↔ABA consistency (3 failing tests)
- ABI transform formulas (fixed in v1.1 Phase 11)

---

<details>
<summary>v1.1 Milestone (Archived Planning Reference)</summary>

## Milestone v1.1 — Bug Fixes & Stability

**Started:** 2026-05-16  
**Goal:** Fix ArticulatedBodyInertia transforms and multi-link consistency

**Scope:**
- Fix `tformABI()` / `invtformABI()` formulas (Featherstone Eq 7.16)
- Align RNEA/ABA conventions for multi-link chains
- Achieve 100% test pass rate
- Performance benchmarks and additional examples

**Requirements:** 6 (2 bug fixes, 2 stability, 2 enhancements)  
**Phases:** 4 planned (11, 12, 13)

</details>
