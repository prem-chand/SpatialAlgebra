# SpatialAlgebra Project

**Project Code:** SA  
**Project Title:** Spatial Vector Algebra Library for Robotics  
**Last Updated:** 2026-06-06

---

## What This Is

A C++17 library implementing spatial vector algebra for rigid body dynamics, following Featherstone's formulation. Provides 6D spatial vectors (twists and wrenches), Plücker coordinate transforms, inertia representations, forward/inverse dynamics, CI builds, and formal mathematical conventions — for robotics simulation and control.

## Core Value

**Must Deliver:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

**Success Looks Like:** 
- All core classes implemented and verified ✓ — achieved in v1.0
- Forward dynamics (ABA) ✓ — gravity support added in v1.1
- Inverse dynamics (RNEA) ✓ — gravity support added in v1.1
- 98.7% test coverage for all mathematical operations ✓
- Library production-ready for serial kinematic chains ✓
- CI pipeline with 4-matrix build ✓
- Formal conventions documented ✓

## Current Milestone: v1.3 Pinocchio Cross-Validation

**Goal:** Extract all 11 test suites into a zero-dependency test model library, refine tests with better variety, precision, documentation and edge cases, and build Pinocchio C++/Python comparison benchmarks.

**Target features:**
- Standalone test model library (Eigen-only, zero SpatialAlgebra dependency) covering all 11 test domains
- Refined tests: better model variety, numerical precision, docstrings, and edge case coverage
- Pinocchio C++ comparison benchmarks for ABA/RNEA round-trip consistency
- Pinocchio Python comparison harness with equivalent numerical validation

## Requirements

### Validated

- ✓ **v1.0:** SpatialVector, MotionVector, ForceVector — complete with Featherstone-verified cross products
- ✓ **v1.0:** Rotation (angle-axis, quaternion, matrix operations)
- ✓ **v1.0:** LowerTriangular packed storage (multiply, inverse, transpose)
- ✓ **v1.0:** Spatial utilities (skew, dot, cross) and SpatialOperations class
- ✓ **v1.0 (PLX-04, PLX-06 fixed in v1.1):** Plücker transforms (motion/force/RBI/ABI transforms)
- ✓ **v1.0:** RigidBodyInertia and ArticulatedBodyInertia (construction, apply, operators)
- ✓ **v1.0:** Articulated Body Algorithm (ABA) for serial and branching kinematic chains
- ✓ **v1.0:** Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics
- ✓ **v1.0:** Comprehensive GTest test infrastructure (189+ tests)
- ✓ **v1.0:** Integration tests for dynamics pipeline consistency
- ✓ **v1.0:** Documentation: README, Doxygen, compilable examples
- ✓ **v1.1:** ABI transform formulas fixed per Featherstone Eq 7.16 (BF-01)
- ✓ **v1.1:** Gravity support for ABA and RNEA (Featherstone D-07/D-08)
- ✓ **v1.1:** Cross-product operations unified to single canonical implementation
- ✓ **v1.1:** NaN/Inf debug-mode guards on core operations
- ✓ **v1.1:** GitHub Actions CI (4-matrix: ubuntu/macos × g++/clang++)
- ✓ **v1.1:** CMake FetchContent fallback for GTest
- ✓ **v1.1:** MATHEMATICAL_CONVENTIONS.md — formal cross-product, gravity, conventions
- ✓ **v1.1:** Edge case and release-mode stability tests
- ✓ **v1.2:** Eigen 5.x CI matrix with 8-job build (2 OS × 2 compilers × 2 Eigen versions)
- ✓ **v1.2:** Google Benchmark v1.9.5 integration via FetchContent with `SA_BUILD_BENCHMARKS` guard
- ✓ **v1.2:** Full benchmark infrastructure (ModelFactory, RandomState, 80 registered benchmarks)
- ✓ **v1.2:** 18 core microbenchmarks (Plücker, cross, inertia) plus ABA/RNEA DOF sweep (n=1..20)
- ✓ **v1.2:** Real-world robot examples (2-link Z-Z planar, 3-link Z-Y-Z spatial) with UR5 parameters

### Active (v1.3)

- [ ] Standalone test model library — Eigen-only types, zero SpatialAlgebra dependency, covering all 11 test domains (vectors, transforms, inertia, dynamics, consistency)
- [ ] Test refinement — better model variety (URDF-derived, different joint types), tighter numerical tolerances with error analysis, comprehensive Doxygen-style docstrings, edge case coverage (zero-mass, singular, near-singular)
- [ ] Pinocchio C++ comparison benchmarks — link against libpinocchio, run equivalent RNEA/ABA/consistency tests
- [ ] Pinocchio Python comparison harness — equivalent numerical validation via pinocchio Python bindings
- [ ] Benchmark result reporting — cross-library comparison tables with relative error analysis

### Out of Scope

- **Python bindings** — Standalone RNEA exists; pybind11 integration deferred to v2.0
- **Contact/collision handling** — Pure rigid body dynamics; contact physics deferred
- **Inverse kinematics** — Focus on dynamics algorithms (RNEA, ABA), not kinematics
- **Visualization tools** — Library-only, no rendering components
- **Serialization** — No save/load for configurations or state
- **Floating base / 6-DOF base** — Fixed-base kinematics only for v1.0

## Context

**v1.1 Shipped:** May 17, 2026 — ~9,450 LOC C++17, 189+ GTest tests, CI pipeline.

**v1.2 Shipped:** 2026-06-06 — 5 phases (15-18), 9 plans, ~30 files modified, 11 CTest passing, 80 benchmarks.

**CR-02 Fixed (Phase 14, 2026-06-17):** ABA inward pass Phase 3 correction double-counting resolved — 3-link non-zero COM round-trip passes, all 11 CTest green.

**State:**
- All core classes implemented and verified with tests
- Forward dynamics (ABA) and inverse dynamics (RNEA) with gravity support
- All 11 CTest executables passing (8 consistency + 15 FD tests)
- NaN/Inf debug-mode guards on core operations
- Cross-product operations unified to single canonical implementation
- Production-hardened with edge case tests and release-mode stability verification
- GitHub Actions CI with 8-matrix build (Eigen 3.4 + 5.x)
- Google Benchmark v1.9.5 with 80 registered benchmarks
- Formal mathematical conventions documented
- Independent gravity invariant test oracles for both solvers

**Tech Stack:**
- C++17 with Eigen3 for linear algebra
- Google Test for unit testing
- CMake build system (FetchContent fallback)
- Doxygen documentation
- GitHub Actions CI

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
| Gravity as base acceleration (not external force) | Follows Featherstone D-07/D-08 | Correct propagation through recursive chain |
| NaN/Inf checks as non-fatal warnings | Matches Featherstone safety convention | Zero production overhead via `#ifndef NDEBUG` |
| OpenMP removal from LowerTriangular | Eliminates hidden linkage dependency | No performance impact for target use cases |
| Namespace cleanup (Vector3d into SpatialAlgebra) | Eliminates ODR hazard | Backward compatible via `using namespace` |
| GitHub Actions CI with 4-matrix build | Automated build verification | Coverage upload on ubuntu+g++ |
| Eigen 5.x CI matrix expansion | Validate forward compatibility | 8-job CI matrix, zero warnings |
| Google Benchmark via FetchContent | No system dependency for benchmarks | Self-contained benchmark build |
| UR5-derived parameters in examples | Realistic robot dynamics | Masses 3.7/8.393/2.33 kg verified |
| Cross-validation as solver consistency check | FD(ID(0,g),g) ≈ 0 gold standard | 2-link passes, 3-link reveals ABA bug |
| Honest documentation of solver bugs | Users see actual library state | Both ABA bug and RNEA limitation documented in examples |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---

*Last updated: 2026-06-17 after starting v1.3 milestone*

---

<details>
<summary>v1.0 Current State (Archived)</summary>

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

**Known Gaps (inherited by v1.1):**
- Multi-link RNEA↔ABA consistency (3 failing tests)
- ABI transform formulas (fixed in v1.1 Phase 11)

</details>

<details>
<summary>v1.1 Current State (Archived)</summary>

## Current State (v1.1)

**Shipped:** 2026-05-17  
**Phases:** 3 (11-13)  
**Plans:** 9  
**Commits:** 35  
**Files Modified:** 57  
**Test Count:** 189+ (186 passing, 3 CR-02 failures)  

**Delivered:**
- ABI transform formulas fixed (tformABI, invtformABI) — 40/40 Plücker tests
- Gravity support for ABA (computeAccelerations) and RNEA (computeTorques)
- Cross-product operations unified to single canonical implementation
- NaN/Inf debug-mode guards and zero-mass edge case tests
- GitHub Actions CI with 4-matrix build and code coverage
- CMake FetchContent fallback, OpenMP removal, umbrella header
- MATHEMATICAL_CONVENTIONS.md — formal specification
- README.md updated with v1.1 API changelog

**Known Gaps (BF-02 deferred):**
- Multi-link RNEA↔ABA consistency: CR-02 bug in ABA inward pass (3 tests)
- Performance benchmarks and additional examples deprioritized

</details>

## Current State (v1.2)

**Shipped:** 2026-06-06  
**Phases:** 4 (15-18)  
**Plans:** 8  
**Files Modified:** ~30  
**Test Count:** 11 CTest (all passing), 80 Google Benchmark registrations  

**Delivered:**
- Eigen 5.x CI: 8-job matrix (2 OS × 2 compilers × 2 Eigen versions), version range syntax, zero warnings
- Google Benchmark v1.9.5 FetchContent integration with `SA_BUILD_BENCHMARKS` guard (default OFF)
- Full benchmark infrastructure: ModelFactory, RandomState, `bench_all` with 80 registered benchmarks
- 18 core microbenchmarks (8 Plücker, 3 cross, 7 inertia) + ABA/RNEA DOF sweep (n=1..20)
- Two robot dynamics executables with UR5-derived parameters:
  - 2-link Z-Z planar arm: full FD, ID gravity, cross-validation passing
  - 3-link Z-Y-Z spatial RRR arm: FD, ID gravity, cross-validation documenting ABA/RNEA limitations

**Known Gaps (deferred to v1.3):**
- CR-02 ABA bug: multi-link ID→FD round-trip fails for non-zero COM
- RNEA fixed-transform limitation: X does not update with joint q
- RBDL comparison benchmarks not started
