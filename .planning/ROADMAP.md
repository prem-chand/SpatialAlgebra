# SpatialAlgebra Roadmap

**Last Updated:** 2026-06-18  
**Granularity:** Milestone

---

## Milestones

:- ✅ **v1.0 MVP / Core Library** — Phases 1-10 (shipped 2026-05-16)
:- ✅ **v1.1 Bug Fixes & Stability** — Phases 11-13 (shipped 2026-05-17)
:- ✅ **v1.2 Production Quality** — Phases 15-18 (shipped 2026-06-06)
:- ✅ **v1.3 Pinocchio Cross-Validation** — Phases 14, 20-24 (shipped 2026-06-18)

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

---

<details>
<summary>✅ v1.3 Pinocchio Cross-Validation (Phases 14, 20-24) — SHIPPED 2026-06-18</summary>

- [x] Phase 14: CR-02 Bug Fix (3/3 plans) — completed 2026-06-17
- [x] Phase 20: Test Model Library (4/4 plans) — completed 2026-06-17
- [x] Phase 21: Test Refinement — completed 2026-06-17 (code delivered)
- [x] Phase 22: Pinocchio C++ Comparison — completed 2026-06-17 (code delivered)
- [x] Phase 23: Pinocchio Python Harness — completed 2026-06-17 (code delivered)
- [x] Phase 24: Result Reporting — completed 2026-06-17 (code delivered)

</details>

## Backlog

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
| 14. CR-02 Bug Fix | v1.3 | 3/3 | Complete | 2026-06-17 |
| 15. Eigen 5.x CI | v1.2 | 1/1 | Complete | 2026-06-05 |
| 16. Benchmark Infrastructure | v1.2 | 3/3 | Complete | 2026-06-05 |
| 17. Benchmark Implementation | v1.2 | 3/3 | Complete | 2026-06-05 |
| 18. Robot Examples | v1.2 | 1/1 | Complete | 2026-06-06 |
| 19. RBDL Comparison | v1.3 | 0/0 | Deferred | - |
| 20. Test Model Library | v1.3 | 4/4 | Complete | 2026-06-17 |
| 21. Test Refinement | v1.3 | — | Complete | 2026-06-17 |
| 22. Pinocchio C++ Comparison | v1.3 | — | Complete | 2026-06-17 |
| 23. Pinocchio Python Harness | v1.3 | — | Complete | 2026-06-17 |
| 24. Result Reporting | v1.3 | — | Complete | 2026-06-17 |

---

## Requirement Coverage

**Total shipped requirements:** 67 (41 v1.0/v1.1 + 10 v1.2 + 16 v1.3)  
**Deferred:** 3 (RBDL-01, RBDL-02, RBDL-03)

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
| BFIX-01 | Phase 14 | Complete |
| CI-01 | Phase 15 | Complete |
| CI-02 | Phase 15 | Complete |
| BINF-01 | Phase 16 | Complete |
| BINF-02 | Phase 16 | Complete |
| BINF-03 | Phase 16 | Complete |
| BENCH-01 | Phase 17 | Complete |
| BENCH-02 | Phase 17 | Complete |
| BENCH-03 | Phase 17 | Complete |
| EX-01 | Phase 18 | Complete |
| RBDL-01 | Phase 19 | Deferred |
| RBDL-02 | Phase 19 | Deferred |
| RBDL-03 | Phase 19 | Deferred |
| TML-01 | Phase 20 | Complete |
| TML-02 | Phase 20 | Complete |
| TML-03 | Phase 20 | Complete |
| TML-04 | Phase 20 | Complete |
| TML-05 | Phase 20 | Complete |
| TST-01 | Phase 21 | Complete |
| TST-02 | Phase 21 | Complete |
| TST-03 | Phase 21 | Complete |
| TST-04 | Phase 21 | Complete |
| PCC-01 | Phase 22 | Complete |
| PCC-02 | Phase 22 | Complete |
| PCC-03 | Phase 22 | Complete |
| PCC-04 | Phase 22 | Complete |
| PCP-01 | Phase 23 | Complete |
| PCP-02 | Phase 23 | Complete |
| PCP-03 | Phase 23 | Complete |
| PCP-04 | Phase 23 | Complete |
| RPT-01 | Phase 24 | Complete |
| RPT-02 | Phase 24 | Complete |
| RPT-03 | Phase 24 | Complete |

---

*See `.planning/milestones/v1.0-ROADMAP.md` for v1.0 details.*
*See `.planning/milestones/v1.1-ROADMAP.md` for v1.1 details.*
*See `.planning/milestones/v1.2-ROADMAP.md` for v1.2 details.*
*See `.planning/milestones/v1.3-ROADMAP.md` for v1.3 details.*

---

