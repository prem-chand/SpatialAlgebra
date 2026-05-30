# SpatialAlgebra Roadmap

**Last Updated:** 2026-05-30  
**Granularity:** Fine

---

## Milestones

- ✅ **v1.0 MVP / Core Library** — Phases 1-10 (shipped 2026-05-16)
- ✅ **v1.1 Bug Fixes & Stability** — Phases 11-13 (shipped 2026-05-17)
- 🚧 **v1.2 Production Quality** — Phases 14-19 (in planning)

## Phases

- [x] **Phase 1: Foundation Vectors** — 6D spatial vector base classes
- [x] **Phase 2: Rotation & Math** — 3D rotation, angle-axis, quaternion
- [x] **Phase 3: Packed Matrix** — Lower-triangular packed storage
- [x] **Phase 4: Spatial Utilities** — skew, dot, cross, SpatialOperations
- [x] **Phase 5: Plücker Transforms** — 6×6 coordinate transforms
- [x] **Phase 6: Inertia Properties** — RigidBodyInertia, ArticulatedBodyInertia
- [x] **Phase 7: Forward Dynamics** — ABA algorithm
- [x] **Phase 8: Test Infrastructure** — GTest, CMake FetchContent
- [x] **Phase 9: Integration Tests** — Dynamics consistency, round-trip
- [x] **Phase 10: Documentation** — README, Doxygen, examples
- [x] **Phase 11: ABI Transform Fixes** — Fixed ABI transform formulas
- [x] **Phase 12: Dynamics Consistency** — Cross-product unification, NaN guards
- [x] **Phase 13: Production Readiness** — CI, CMake FetchContent, conventions doc
- [ ] **Phase 14: CR-02 Bug Fix** — Restructure ABA inward pass per Featherstone Algorithm 7.3
- [ ] **Phase 15: Eigen 5.x CI** — Version range syntax + CI matrix expansion
- [ ] **Phase 16: Benchmark Infrastructure** — Google Benchmarks FetchContent, `benchmarks/` dir
- [ ] **Phase 17: Benchmark Implementation** — ABA/RNEA timing with DOF sweep, microbenchmarks
- [ ] **Phase 18: Robot Examples** — 2-link planar and 3-link spatial arm examples
- [ ] **Phase 19: RBDL Comparison** — Optional RBDL dependency and comparison benchmarks

<details>
<summary>✅ v1.0 MVP / Core Library (Phases 1-10) — SHIPPED 2026-05-16</summary>

- [x] Phase 1: Foundation Vectors (4/4 plans) — completed 2026-05-15
- [x] Phase 2: Rotation & Math (1/1 plan) — completed 2026-05-15
- [x] Phase 3: Packed Matrix (2/2 plans) — completed 2026-05-15
- [x] Phase 4: Spatial Utilities (3/3 plans) — completed 2026-05-16
- [x] Phase 5: Plücker Transforms (3/3 plans) — completed 2026-05-16
- [x] Phase 6: Inertia Properties (3/3 plans) — completed 2026-05-16
- [x] Phase 7: Forward Dynamics (2/2 plans) — completed 2026-05-16
- [x] Phase 8: Test Infrastructure (2/2 plans) — completed 2026-05-16
- [x] Phase 9: Integration Tests (2/2 plans) — completed 2026-05-16
- [x] Phase 10: Documentation (3/3 plans) — completed 2026-05-16

</details>

<details>
<summary>✅ v1.1 Bug Fixes & Stability (Phases 11-13) — SHIPPED 2026-05-17</summary>

- [x] Phase 11: ABI Transform Fixes (1/1 plan) — completed 2026-05-16
- [x] Phase 12: Dynamics Consistency Fixes (1/1 plan) — completed 2026-05-16 (partial)
- [x] Phase 13: Production Readiness (7/7 plans) — completed 2026-05-17

</details>

---

## Phase Details

### 🚧 v1.2 Production Quality (In Planning)

**Milestone Goal:** Close all remaining gaps — fix multi-link dynamics consistency (CR-02), establish performance benchmarks vs RBDL, add Eigen 5.x to CI matrix, and provide real-world robot examples.

---

### Phase 14: CR-02 Bug Fix
**Goal**: Correct forward dynamics for chains with 3+ joints — all 4 multi-link consistency tests pass
**Depends on**: Phase 13 (previous milestone; independent of Phase 15)
**Requirements**: BFIX-01
**Success Criteria** (what must be TRUE):
   1. All 4 multi-link consistency tests pass: ThreeLinkSerialChain, BranchingYConfiguration, TwoLinkRoundTrip, ThreeLinkNumericalValidation
   2. ABA forward dynamics produces bitwise-consistent results with RNEA∘Inverse round-trip for 3+ link chains
   3. All 156 existing tests continue to pass (no regressions)
   4. Bias forces (pa) computed after child inertia accumulation per Featherstone Algorithm 7.3
**Plans**: 1 plan

Plans:
- [ ] 14-01-PLAN.md — Restructure ABA inward pass with condensation, update tests

### Phase 15: Eigen 5.x CI
**Goal**: Library builds and tests pass under Eigen 5.x in CI matrix
**Depends on**: Phase 13 (previous milestone; independent of Phase 14)
**Requirements**: CI-01, CI-02
**Success Criteria** (what must be TRUE):
  1. CMake `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` version range syntax accepted by Eigen 3.4.x and 5.x
  2. Library compiles without warnings under Eigen 5.0.1 on both Ubuntu (g++) and macOS (clang++)
  3. CI matrix includes Eigen 5.0.1 entries alongside existing 3.4.x builds
  4. All tests pass under Eigen 5.x with no failures
**Plans**: TBD

### Phase 16: Benchmark Infrastructure
**Goal**: Build system and shared utilities for performance benchmarks
**Depends on**: Phase 14 (requires correct solver)
**Requirements**: BINF-01, BINF-02, BINF-03
**Success Criteria** (what must be TRUE):
  1. Google Benchmark v1.9.5 integrated via FetchContent in `benchmarks/CMakeLists.txt`
  2. `benchmarks/` directory exists with separate CMakeLists.txt and `SA_BUILD_BENCHMARKS` guard (default OFF)
  3. Shared benchmark utilities available: parameterized model factory for arbitrary n-DOF chains, random joint state generator
  4. Benchmarks compile at `-O3 -DNDEBUG` and link correctly against libSpatialAlgebra.a
**Plans**: TBD

### Phase 17: Benchmark Implementation
**Goal**: Executable microbenchmarks for ABA, RNEA, and core operations
**Depends on**: Phase 16 (requires benchmark infrastructure)
**Requirements**: BENCH-01, BENCH-02, BENCH-03
**Success Criteria** (what must be TRUE):
  1. `bench_aba` runs with parameterized DOF sweep n=1..20 and outputs timing results with statistical rigor
  2. `bench_rnea` runs with parameterized DOF sweep n=1..20 and outputs timing results
  3. `bench_plucker` and `bench_cross_product` microbenchmark executables produce stable per-operation timing
  4. Each iteration generates fresh random joint states to avoid warm-state bias
**Plans**: TBD

### Phase 18: Robot Examples
**Goal**: Real-world robot examples demonstrating correct physics
**Depends on**: Phase 14 (requires correct solver; independent of Phases 16-17)
**Requirements**: EX-01, EX-02
**Success Criteria** (what must be TRUE):
  1. `example_robot_2link` compiles, runs, and prints physically correct forward/inverse dynamics for a Z-Z planar arm
  2. `example_robot_3link` compiles, runs, and prints physically correct dynamics for a Z-Y-Z spatial RRR arm
  3. Both examples demonstrate gravity-compensated torque output that matches expected static equilibrium
**Plans**: TBD

### Phase 19: RBDL Comparison
**Goal**: Optional RBDL-based comparison benchmarks with numerical identity verification
**Depends on**: Phase 17, Phase 18 (requires benchmarks and examples)
**Requirements**: RBDL-01, RBDL-02, RBDL-03
**Success Criteria** (what must be TRUE):
  1. `cmake/FindRBDL.cmake` locates RBDL v3.3.1 installation; comparison benchmarks guarded by `SA_BUILD_COMPARISON_BENCHMARKS` (default OFF)
  2. Numerical identity verified between SpatialAlgebra and RBDL for single-link dynamics to 1e-12
  3. ABA/RNEA timing comparison executable runs for 3+ link serial chains
  4. Frame convention mapping documented (transform directions, joint screw conventions)
**Plans**: TBD

---

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|---------------|--------|-----------|
| 1. Foundation Vectors | v1.0 | 4/4 | Complete | 2026-05-15 |
| 2. Rotation & Math | v1.0 | 1/1 | Complete | 2026-05-15 |
| 3. Packed Matrix | v1.0 | 2/2 | Complete | 2026-05-15 |
| 4. Spatial Utilities | v1.0 | 3/3 | Complete | 2026-05-16 |
| 5. Plücker Transforms | v1.0 | 3/3 | Complete | 2026-05-16 |
| 6. Inertia Properties | v1.0 | 3/3 | Complete | 2026-05-16 |
| 7. Forward Dynamics | v1.0 | 2/2 | Complete | 2026-05-16 |
| 8. Test Infrastructure | v1.0 | 2/2 | Complete | 2026-05-16 |
| 9. Integration Tests | v1.0 | 2/2 | Complete | 2026-05-16 |
| 10. Documentation | v1.0 | 3/3 | Complete | 2026-05-16 |
| 11. ABI Transform Fixes | v1.1 | 1/1 | Complete | 2026-05-16 |
| 12. Dynamics Consistency | v1.1 | 1/1 | Partial | 2026-05-16 |
| 13. Production Readiness | v1.1 | 7/7 | Complete | 2026-05-17 |
| 14. CR-02 Bug Fix | v1.2 | 0/1 | Not started | - |
| 15. Eigen 5.x CI | v1.2 | 0/0 | Not started | - |
| 16. Benchmark Infrastructure | v1.2 | 0/0 | Not started | - |
| 17. Benchmark Implementation | v1.2 | 0/0 | Not started | - |
| 18. Robot Examples | v1.2 | 0/0 | Not started | - |
| 19. RBDL Comparison | v1.2 | 0/0 | Not started | - |

---

## Requirement Coverage

**Total v1 requirements:** 55 (41 v1.0/v1.1 + 14 v1.2)  
**Mapped:** 55/55 ✓

| Requirement | Phase | Status |
|-------------|-------|--------|
| VEC-01 | Phase 1 | Complete |
| VEC-02 | Phase 1 | Complete |
| VEC-03 | Phase 1 | Complete |
| VEC-04 | Phase 1 | Complete |
| ROT-01 | Phase 2 | Complete |
| ROT-02 | Phase 2 | Complete |
| ROT-03 | Phase 2 | Complete |
| ROT-04 | Phase 2 | Complete |
| LTR-01 | Phase 3 | Complete |
| LTR-02 | Phase 3 | Complete |
| LTR-03 | Phase 3 | Complete |
| LTR-04 | Phase 3 | Complete |
| UTL-01 | Phase 4 | Complete |
| UTL-02 | Phase 4 | Complete |
| UTL-03 | Phase 4 | Complete |
| UTL-04 | Phase 4 | Complete |
| PLX-01 | Phase 5 | Complete |
| PLX-02 | Phase 5 | Complete |
| PLX-03 | Phase 5 | Complete |
| PLX-04 | Phase 5 | Complete |
| PLX-05 | Phase 5 | Complete |
| PLX-06 | Phase 5 | Complete |
| INR-01 | Phase 6 | Complete |
| INR-02 | Phase 6 | Complete |
| INR-03 | Phase 6 | Complete |
| INR-04 | Phase 6 | Complete |
| ABA-01 | Phase 7 | Complete |
| ABA-02 | Phase 7 | Complete |
| ABA-03 | Phase 7 | Complete |
| ABA-04 | Phase 7 | Complete |
| TST-01 | Phase 8 | Complete |
| TST-02 | Phase 6 | Complete |
| TST-03 | Phase 6 | Complete |
| TST-04 | Phase 8 | Complete |
| TST-05 | Phase 3 | Complete |
| TST-06 | Phase 8 | Complete |
| TST-07 | Phase 9 | Complete |
| DOC-01 | Phase 10 | Complete |
| DOC-02 | Phase 10 | Complete |
| DOC-03 | Phase 10 | Complete |
| BFIX-01 | Phase 14 | Pending |
| CI-01 | Phase 15 | Pending |
| CI-02 | Phase 15 | Pending |
| BINF-01 | Phase 16 | Pending |
| BINF-02 | Phase 16 | Pending |
| BINF-03 | Phase 16 | Pending |
| BENCH-01 | Phase 17 | Pending |
| BENCH-02 | Phase 17 | Pending |
| BENCH-03 | Phase 17 | Pending |
| EX-01 | Phase 18 | Pending |
| EX-02 | Phase 18 | Pending |
| RBDL-01 | Phase 19 | Pending |
| RBDL-02 | Phase 19 | Pending |
| RBDL-03 | Phase 19 | Pending |

---

*See `.planning/milestones/v1.0-ROADMAP.md` for v1.0 details.*
*See `.planning/milestones/v1.1-ROADMAP.md` for v1.1 details.*
