# SpatialAlgebra — Retrospective

## Milestone: v1.2 — Production Quality

**Shipped:** 2026-06-06
**Phases:** 4 (15-18) | **Plans:** 8 | **Tasks:** 13

### What Was Built

- Expanded CI matrix from 4 to 8 jobs with Eigen version range syntax, verified zero-error compilation under Eigen 5.0.1
- Google Benchmark v1.9.5 FetchContent integration, SA_BUILD_BENCHMARKS guard, and full benchmarks/ directory
- Shared benchmark utilities: ModelFactory for n-DOF solver construction, RandomState for deterministic joint state generation
- 18 core microbenchmarks (Plücker, cross, inertia) and ABA/RNEA DOF sweep (n=1..20)
- Two robot dynamics executables with UR5-derived parameters, cross-validation revealing pre-existing solver limitations

### What Worked

- Parallel phase execution (benchmarks and Eigen 5.x) reduced total wall time
- UR5 parameter research enabled realistic physics in examples without real robot data
- Cross-validation pattern caught ABA bug that existing tests missed (degenerate test parameters)
- Plan-checker gate caught missing VALIDATION.md before execution

### What Was Inefficient

- Phase 14 (CR-02 fix) never started — deferred through entire milestone; should have been prioritized first
- Phase 19 (RBDL comparison) also deferred — both incomplete phases should be scoped to v1.3 explicitly
- Phase 14 dependency on benchmark phases was incorrectly stated in ROADMAP.md (benchmarks don't depend on correct multi-link ABA)

### Patterns Established

- Benchmark infrastructure as a separate build target (SA_BUILD_BENCHMARKS guard)
- Cross-validation (FD↔ID round-trip) as gold standard for solver correctness
- Honest bug documentation in example output rather than hiding limitations
- Deterministic RNG with fixed seed for reproducible benchmark inputs

### Key Lessons

- Test degeneracy can hide real bugs: zero-COM tests pass while non-zero-COM chains fail
- Phase dependency ordering must be verified against actual code coupling, not assumed
- Robot examples with realistic parameters are more valuable than synthetic tests for catching bugs

### Cost Observations

- Model mix: ~70% sonnet, ~30% opus (planning/checking)
- Sessions: 5 (Phase 15: 1, Phase 16: 1, Phase 17: 1, Phase 18: 1, Phase 18 execution: 1)
- Notable: Plan-checker caught 1 blocker before execution, preventing invalid state commit

---

## Cross-Milestone Trends

*(No previous milestone retrospective data — this is the first entry.)*
