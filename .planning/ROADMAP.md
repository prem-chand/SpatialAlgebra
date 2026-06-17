# SpatialAlgebra Roadmap

**Last Updated:** 2026-06-06  
**Granularity:** Milestone

---

## Milestones

- ✅ **v1.0 MVP / Core Library** — Phases 1-10 (shipped 2026-05-16)
- ✅ **v1.1 Bug Fixes & Stability** — Phases 11-13 (shipped 2026-05-17)
- ✅ **v1.2 Production Quality** — Phases 15-18 (shipped 2026-06-06); Phases 14, 19 deferred
- 🚧 **v1.3 TBD** — Phases 14, 19 (planned)

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
- [x] 14-01-PLAN.md — Cross product fixes (SpatialUtils.h) — completed 2026-06-04 (partial: tests pass but only with zero COM)
- [ ] 14-02-PLAN.md — RED: Add non-zero COM tests that FAIL (closes SC1, SC2)
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
