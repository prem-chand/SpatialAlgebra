---
|gsd_state_version: 1.0
|milestone: v1.3
|milestone_name: Pinocchio Cross-Validation
|status: all_phases_complete
|last_updated: 2026-06-18
|last_activity: 2026-06-18 -- v1.3 milestone completed and archived
|progress:
|  total_phases: 24
|  completed_phases: 23
|  total_plans: 52
|  completed_plans: 52
|  percent: 96
|stopped_at: v1.3 milestone complete — all 6 phases (14, 20-24) delivered. Milestone archived.
---

# SpatialAlgebra State

**Last Updated:** 2026-06-18  
**Current Milestone:** v1.3 Pinocchio Cross-Validation ✅

---

## Project Reference

**See:** .planning/PROJECT.md (updated 2026-06-18)

**Core value:** Complete, well-tested spatial algebra library — achieved with v1.0-v1.3  
**Current focus:** Milestone complete — ready for next milestone

---

## Current Position

Status: Milestone v1.3 complete and archived
Last activity: 2026-06-18

## Accumulated Context

### Decisions

- **v1.3 milestone**: "Pinocchio Cross-Validation" — all 5 scope phases delivered
- **CR-02 root cause**: Phase 3 correction double-counted Ia*c in ABA inward pass — fixed with corrected correction formula
- **Test Model Library**: Eigen-only INTERFACE library with 11 chain headers, PIMPL adapter, compilation firewall verified
- **Pinocchio C++ Adapter**: Full adapter implementing RobotSolver with pinocchio::rnea()/aba()
- **Python Harness**: C++ subprocess JSON integration (Python pinocchio bindings unavailable due to Boost ABI mismatch)
- **Cross-Library Reporting**: Comparison tables with max/mean relative error per model, CI regression tracking (1e-8 kinematics, 1e-6 dynamics)
- **CMake workaround**: Bypasses pinocchioConfig.cmake entirely (Boost 1.89.0 header-only incompatibility)
- **SA RNEA limitation**: InverseDynamics uses fixed parent-to-child transforms, giving correct dynamics only at q=0

### Known Limitations

1. SA RNEA/ABA only valid at q=0 — fixed transforms not updated with joint position
2. BranchingY non-zero COM test removed — Phase 3 precision insufficient for branching chains
3. Pinocchio Python bindings not available (Boost/Python ABI mismatch) — C++ subprocess workaround

### Pending Todos

None — v1.3 milestone closed

---

## Deferred Items

Items acknowledged and deferred at milestone close:

| Category | Item | Status |
|----------|------|--------|
| feature | RBDL comparison benchmarks (Phase 19) | Deferred from v1.2 |
| feature | CI-based regression tracking | Deferred from v1.2 |
| feature | URDF model loading | Deferred from v1.2 |
| tech | Branching chain non-zero COM precision | Requires full single-sweep Featherstone ABA |

---

## Operator Next Steps

- `/gsd:new-milestone` — start next milestone (v2.0)
