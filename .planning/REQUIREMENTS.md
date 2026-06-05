# Requirements: SpatialAlgebra

**Defined:** 2026-05-30
**Core Value:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

## v1.2 Requirements

Requirements for v1.2 "Production Quality" milestone. Each maps to roadmap phases.

### Bug Fixes

- [ ] **BFIX-01**: CR-02 ABA inward pass restructured per Featherstone Algorithm 7.3 — all 4 multi-link consistency tests pass

### Build & CI

- [x] **CI-01**: Eigen 5.x compatible via version range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)`
- [x] **CI-02**: Eigen 5.x added to CI build matrix

### Benchmark Infrastructure

- [ ] **BINF-01**: Google Benchmark v1.9.5 integrated via FetchContent
- [ ] **BINF-02**: `benchmarks/` directory with separate CMakeLists.txt guarded by `SA_BUILD_BENCHMARKS` (default OFF)
- [ ] **BINF-03**: Shared benchmark utilities (model factory, random state generator)

### Benchmark Implementation

- [ ] **BENCH-01**: ABA forward dynamics timing with parameterized DOF sweep (n=1..20)
- [ ] **BENCH-02**: RNEA inverse dynamics timing with parameterized DOF sweep (n=1..20)
- [ ] **BENCH-03**: Plücker transform and cross-product microbenchmarks

### Robot Examples

- [ ] **EX-01**: 2-link planar robot dynamics example
- [ ] **EX-02**: 3-link spatial arm dynamics example

### RBDL Comparison

- [ ] **RBDL-01**: RBDL v3.3.1 optional dependency with `FindRBDL.cmake` module
- [ ] **RBDL-02**: Numerical identity verification against RBDL for single-link dynamics
- [ ] **RBDL-03**: ABA/RNEA timing comparison vs RBDL for 3+ link serial chains

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Advanced Features

- **ADV-01**: Performance regression tracking in CI
- **ADV-02**: Memory profiling of dynamics solvers
- **ADV-03**: URDF-based model loading
- **ADV-04**: Python bindings via pybind11
- **ADV-05**: Pinocchio comparison benchmarks

## Out of Scope

| Feature | Reason |
|---------|--------|
| Python benchmark harness | Python is not on the performance-critical path |
| CI benchmark regression tests | Benchmarks are stochastic and machine-dependent |
| URDF parsing | Hand-constructed chains sufficient for v1.2 |
| Floating base dynamics | Fixed-base only (matched to existing library scope) |
| OpenMP threading | Already removed; single-threaded benchmarks only |
| nanobench / custom timers | Google Benchmark is the standard for C++ microbenchmarking |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| BFIX-01 | Phase 14 | Pending |
| CI-01 | Phase 15 | Complete |
| CI-02 | Phase 15 | Complete |
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

**Coverage:**
- v1.2 requirements: 14 total
- Mapped to phases: 14
- Unmapped: 0 ✓

---
*Requirements defined: 2026-05-30*
*Last updated: 2026-05-30 after v1.2 milestone definition*
