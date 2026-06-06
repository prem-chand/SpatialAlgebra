# SpatialAlgebra Milestones

## v1.2 Production Quality (Shipped: 2026-06-06)

**Phases completed:** 5 phases, 9 plans, 13 tasks

**Key accomplishments:**

- Plan:
- Expanded CI matrix from 4 to 8 jobs with Eigen version range syntax, verified zero-error compilation and all 11 tests passing under Eigen 5.0.1
- Google Benchmark v1.9.5 FetchContent integration, SA_BUILD_BENCHMARKS guard (default OFF), and full benchmarks/ directory structure with per-domain CMakeLists.txt files
- Shared benchmark utilities: ModelFactory for n-DOF solver construction and RandomState for deterministic random joint state generation
- Programmatic RegisterBenchmark DOF sweep stubs for ABA, RNEA, and core microbenchmarks, producing a linkable bench_all executable with 80 registered benchmarks
- Two robot dynamics executables (2-link Z-Z and 3-link Z-Y-Z arms) with UR5-derived parameters, cross-validation revealing two pre-existing solver limitations: ABA multi-link COM propagation bug and RNEA fixed-transform limitation

---

## v1.1 — Bug Fixes & Stability

**Shipped:** 2026-05-17
**Phases:** 3 (11-13) | **Plans:** 9 | **Tasks:** 18 | **Commits:** 35
**Files modified:** 57 | **Timeline:** 1 day

### Accomplishments

1. Fixed ABI transform formulas and Plücker inverse — 40/40 Plücker tests passing
2. Added gravity support to ForwardDynamics (ABA) and InverseDynamics (RNEA)
3. Production hardening: NaN/Inf guards, zero-mass edge case tests, release-mode stability
4. Build & CI: GitHub Actions workflow (4-matrix), CMake FetchContent, OpenMP removal, umbrella header
5. Mathematical conventions documented in MATHEMATICAL_CONVENTIONS.md
6. Comprehensive dynamics tests: multi-link ABA, gravity invariants, Coriolis, round-trip consistency

### Known Gaps

- BF-02: Multi-link RNEA↔ABA consistency — 2/4 tests failing (CR-02 bug in ABA inward pass)
- EN-01, EN-02: Performance benchmarks and additional examples deprioritized

### Deferred Items

No deferred items at close.

---

## v1.0 — MVP / Core Library

**Shipped:** 2026-05-16
**Phases:** 10 (1-10) | **Plans:** 25 | **Tasks:** 41

*See `.planning/milestones/v1.0-ROADMAP.md` for full details.*
