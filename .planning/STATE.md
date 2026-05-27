---
gsd_state_version: 1.1
milestone: v1.1
milestone_name: v1.1-bug-fixes-and-stability
status: completed
last_updated: "2026-05-27T14:00:00.000Z"
progress:
  total_phases: 13
  completed_phases: 13
  total_plans: 34
  completed_plans: 34
  percent: 100
---

# SpatialAlgebra State

**Last Updated:** 2026-05-27  
**Current Milestone:** v1.1 ✅ ARCHIVED

---

## v1.1 Summary

**Shipped:** 2026-05-17 | **Archived:** 2026-05-27

**Delivered:**
- ABI transform formulas fixed (tformABI, invtformABI) per Featherstone Eq 7.16
- Gravity support for ForwardDynamics (ABA) and InverseDynamics (RNEA)
- Cross-product operations unified to single canonical implementation
- NaN/Inf debug-mode guards on core operations
- GitHub Actions CI with 4-matrix build and code coverage
- CMake FetchContent fallback, OpenMP removal, umbrella header
- MATHEMATICAL_CONVENTIONS.md formal specification
- Edge case tests (zero-mass, release-mode stability)
- README.md updated with v1.1 API changelog

**Known Gaps (deferred to next milestone):**
- BF-02: Multi-link RNEA↔ABA consistency (3 failing tests, CR-02 bug)
- Performance benchmarks and additional examples deprioritized

---

## Project Reference

**See:** .planning/PROJECT.md (updated 2026-05-27)

**Core value:** Complete, well-tested spatial algebra library — achieved with both v1.0 and v1.1
**Current focus:** Planning next milestone

---

## Archived Artifacts

- `.planning/milestones/v1.0-ROADMAP.md` — v1.0 milestone phase details
- `.planning/milestones/v1.0-REQUIREMENTS.md` — v1.0 requirements with outcomes
- `.planning/milestones/v1.1-ROADMAP.md` — v1.1 milestone phase details
- `.planning/milestones/v1.1-REQUIREMENTS.md` — v1.1 requirements with outcomes
- `.planning/ROADMAP.md` — Updated with milestone grouping (v1.0, v1.1 collapsed)
- `.planning/PROJECT.md` — Current state with Validated requirements
- `.planning/MILESTONES.md` — Milestone summary entries

---

*Milestones v1.0 and v1.1 archived: 2026-05-27*
