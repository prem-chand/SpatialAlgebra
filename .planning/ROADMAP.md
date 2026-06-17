# SpatialAlgebra Roadmap

**Last Updated:** 2026-06-06  
**Granularity:** Milestone

---

## Milestones

- ✅ **v1.0 MVP / Core Library** — Phases 1-10 (shipped 2026-05-16)
- ✅ **v1.1 Bug Fixes & Stability** — Phases 11-13 (shipped 2026-05-17)
- ✅ **v1.2 Production Quality** — Phases 15-18 (shipped 2026-06-06); Phase 14 completed 2026-06-17
- 🚧 **v1.3 Pinocchio Cross-Validation** — Phases 20-24 (active)

## Phases

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

<details>
<summary>✅ v1.2 Production Quality (Phases 15-18) — SHIPPED 2026-06-06</summary>

- [x] Phase 15: Eigen 5.x CI (1/1 plan) — completed 2026-06-05
- [x] Phase 16: Benchmark Infrastructure (3/3 plans) — completed 2026-06-05
- [x] Phase 17: Benchmark Implementation (3/3 plans) — completed 2026-06-05
- [x] Phase 18: Robot Examples (1/1 plan) — completed 2026-06-06

</details>

### 🚧 v1.3 TBD (Planned)

- [ ] Phase 14: CR-02 Bug Fix — gap closure: fix ABA multi-link consistency (2 plans)
- [ ] Phase 19: RBDL Comparison — comparison benchmarks vs RBDL (TBD plans)

---

## Phase Details

### Phase 14: CR-02 Bug Fix

**Goal**: Correct forward dynamics for chains with 3+ joints and non-zero COM — TDD with Red→Green→Verify cycle
**Depends on**: Phase 13 (previous milestone; Phase 14-01 cross product fix already completed); independent of Phase 15
**Requirements**: BFIX-01
**Success Criteria** (what must be TRUE):

   1. All multi-link consistency tests with non-zero COM pass: ThreeLinkSerialChainNonZeroCOM, BranchingYNonZeroCOM, TwoLinkGravityNonZeroCOM
   2. ThreeLinkNumericalValidation updated with non-zero COM round-trip check — passes within EPSILON=1e-8
   3. ABA inward pass restructured to single tip-to-base sweep per Featherstone Algorithm 7.3 (no Phase 3 correction)
   4. All 11 test executables pass with zero regressions (existing zero-COM tests + new non-zero COM tests)
   5. Bias forces (pa) computed after child inertia accumulation per Featherstone Algorithm 7.3

**Plans**: 3 plans (1 completed cross product fix + 2 gap closure)

Plans:
**Wave 1**

- [x] 14-01-PLAN.md — Cross product fixes (SpatialUtils.h) — completed 2026-06-04 (partial: tests pass but only with zero COM)
- [ ] 14-02-PLAN.md — RED: Add non-zero COM tests that FAIL (closes SC1, SC2)

**Wave 2** *(blocked on Wave 1 completion)*

- [ ] 14-03-PLAN.md — GREEN+VERIFY: Restructure inwardPass + Doxygen + full regression (closes SC3, SC4, SC5)

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
| 12. Dynamics Consistency | v1.1 | 1/1 | Complete | 2026-05-16 |
| 13. Production Readiness | v1.1 | 7/7 | Complete | 2026-05-17 |
| 14. CR-02 Bug Fix | v1.3 | 0/3 | Gap Closure | - |
| 15. Eigen 5.x CI | v1.2 | 1/1 | Complete | 2026-06-05 |
| 16. Benchmark Infrastructure | v1.2 | 3/3 | Complete | 2026-06-05 |
| 17. Benchmark Implementation | v1.2 | 3/3 | Complete | 2026-06-05 |
| 18. Robot Examples | v1.2 | 1/1 | Complete | 2026-06-06 |
| 19. RBDL Comparison | v1.3 | 0/0 | Deferred | - |

---

## Requirement Coverage

**Total shipped requirements:** 51 (41 v1.0/v1.1 + 10 v1.2)  
**Deferred to v1.3:** 4 (BFIX-01, RBDL-01, RBDL-02, RBDL-03)

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
| BFIX-01 | Phase 14 | Deferred to v1.3 |
| CI-01 | Phase 15 | Complete |
| CI-02 | Phase 15 | Complete |
| BINF-01 | Phase 16 | Complete |
| BINF-02 | Phase 16 | Complete |
| BINF-03 | Phase 16 | Complete |
| BENCH-01 | Phase 17 | Complete |
| BENCH-02 | Phase 17 | Complete |
| BENCH-03 | Phase 17 | Complete |
| EX-01 | Phase 18 | Complete |
| EX-02 | Phase 18 | Complete |
| RBDL-01 | Phase 19 | Deferred to v1.3 |
| RBDL-02 | Phase 19 | Deferred to v1.3 |
| RBDL-03 | Phase 19 | Deferred to v1.3 |

---

*See `.planning/milestones/v1.0-ROADMAP.md` for v1.0 details.*
*See `.planning/milestones/v1.1-ROADMAP.md` for v1.1 details.*
*See `.planning/milestones/v1.2-ROADMAP.md` for v1.2 details.*

---

## 🚧 v1.3 Pinocchio Cross-Validation (Phases 20-24)

**Goal:** Extract all test models into a zero-dependency Eigen-only library with adapter interfaces, refine test coverage and precision, and build Pinocchio C++/Python comparison benchmarks with result reporting.

| # | Phase | Goal | Requirements | Status |
|---|-------|------|--------------|--------|
| 20 | Test Model Library | Extract all 11 test domains into header-only Eigen-only library with adapter interface | TML-01..05 | Planned (4 plans) |
| 21 | Test Refinement | Better variety, precision, documentation, edge cases | TST-01..04 | Pending |
| 22 | Pinocchio C++ Comparison | C++ adapter and comparison benchmarks | PCC-01..04 | Pending |
| 23 | Pinocchio Python Harness | Python comparison script and JSON output | PCP-01..04 | Pending |
| 24 | Result Reporting | Cross-library comparison tables with regression tracking | RPT-01..03 | Pending |

### Phase Details

**Phase 20: Test Model Library**
- Goal: Create standalone header-only library with zero SA dependency, adapter interface, and all kinematic model definitions
- Requirements: TML-01, TML-02, TML-03, TML-04, TML-05
- **Plans:** 4 plans
- Success criteria:
  1. `tests/test-models/` directory with `<test-models/*.h>` headers using only Eigen types
  2. Adapter interface `RobotSolver` with pure virtual `computeTorques()` and `computeAccelerations()` (per D-06)
  3. SpatialAlgebra adapter wrapping existing InverseDynamics/ForwardDynamics classes
  4. CMake INTERFACE library target `test_models` compiles without any SpatialAlgebra headers
  5. All 11 test domains represented (vectors, transforms, rotation, LT, RBI, ABI, utils, ID, FD, consistency, spatial ops)

Plans:
- [ ] 20-01-PLAN.md — Core data types (RobotModel, RobotSolver) + CMake infrastructure (INTERFACE test_models)
- [ ] 20-02-PLAN.md — Model factory functions: 7 simple-domain chains (spatial_vectors through spatial_utils)
- [ ] 20-03-PLAN.md — Model factory functions: 4 dynamics-domain chains (inverse_dynamics, forward_dynamics, consistency, spatial_operations)
- [ ] 20-04-PLAN.md — SpatialAlgebra adapter (PIMPL) + compile_smoke_test zero-dependency verification

**Phase 21: Test Refinement**
- Goal: Improve test model variety, numerical precision analysis, docstrings, and edge case coverage
- Requirements: TST-01, TST-02, TST-03, TST-04
- Success criteria:
  1. 3+ new kinematic configurations (prismatic joints, mixed types, high-DOF chains)
  2. Relative error reporting alongside absolute EXPECT_NEAR assertions
  3. Every test model struct/class has Doxygen docstring with Featherstone reference
  4. Edge case tests pass: near-zero mass, singular configs, n>10 DOF, non-identity rotations

**Phase 22: Pinocchio C++ Comparison**
- Goal: Build Pinocchio adapter and comparison benchmark executable
- Requirements: PCC-01, PCC-02, PCC-03, PCC-04
- Success criteria:
  1. Pinocchio C++ adapter passes all test models through adapter interface
  2. CMake `SA_BUILD_PINOCCHIO_BENCHMARKS` guard works (default OFF)
  3. Comparison executable runs all models through both adapters
  4. Per-joint relative error < 1e-6 for ABA/RNEA equivalence

**Phase 23: Pinocchio Python Harness**
- Goal: Python comparison script using pinocchio Python bindings
- Requirements: PCP-01, PCP-02, PCP-03, PCP-04
- Success criteria:
  1. `pip install pin` works and Python harness runs successfully
  2. Python models mirror all C++ kinematic chain definitions
  3. Round-trip RNEA↔ABA consistency within 1e-6 for all models
  4. JSON output file with comparison results generated

**Phase 24: Result Reporting**
- Goal: Cross-library comparison tables, error analysis, CI integration
- Requirements: RPT-01, RPT-02, RPT-03
- Success criteria:
  1. Cross-library comparison table (SA C++ vs Pinocchio C++ vs Pinocchio Python)
  2. Per-joint relative error analysis with max/mean error per model
  3. CI test passes/fails based on tolerance thresholds

### Requirement Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| TML-01 | Phase 20 | Pending |
| TML-02 | Phase 20 | Pending |
| TML-03 | Phase 20 | Pending |
| TML-04 | Phase 20 | Pending |
| TML-05 | Phase 20 | Pending |
| TST-01 | Phase 21 | Pending |
| TST-02 | Phase 21 | Pending |
| TST-03 | Phase 21 | Pending |
| TST-04 | Phase 21 | Pending |
| PCC-01 | Phase 22 | Pending |
| PCC-02 | Phase 22 | Pending |
| PCC-03 | Phase 22 | Pending |
| PCC-04 | Phase 22 | Pending |
| PCP-01 | Phase 23 | Pending |
| PCP-02 | Phase 23 | Pending |
| PCP-03 | Phase 23 | Pending |
| PCP-04 | Phase 23 | Pending |
| RPT-01 | Phase 24 | Pending |
| RPT-02 | Phase 24 | Pending |
| RPT-03 | Phase 24 | Pending |
