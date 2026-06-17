# Requirements: SpatialAlgebra

**Defined:** 2026-06-17
**Core Value:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

## v1.3 Requirements

Requirements for Pinocchio cross-validation milestone. Each maps to roadmap phases.

### TML: Test Model Library

- [ ] **TML-01**: Header-only test model library using only Eigen types (Vector3d, Matrix3d, MatrixXd) — zero SpatialAlgebra dependency
- [ ] **TML-02**: Kinematic chain definitions covering all 11 test domains (spatial vectors, Plücker transforms, rotation, lower triangular, rigid body inertia, articulated body inertia, spatial utils, inverse dynamics, forward dynamics, consistency, spatial operations)
- [ ] **TML-03**: Adapter interface (abstract base class) with `computeTorques(q, qdot, qddot)` and `computeAccelerations(tau)` pure virtual methods
- [ ] **TML-04**: SpatialAlgebra adapter implementing the adapter interface as a thin wrapper around existing solver classes
- [ ] **TML-05**: CMake library target `TestModels` with Eigen3 as sole dependency, installable as a standalone library

### TST: Test Refinement

- [ ] **TST-01**: 3+ new kinematic chain configurations beyond existing (prismatic joints, mixed revolute+prismatic, URDF-inspired multi-DOF chains)
- [ ] **TST-02**: Tighter numerical tolerances with relative error analysis — report relative error |computed - expected| / |expected| alongside absolute error
- [ ] **TST-03**: Comprehensive Doxygen docstrings on every test model class/struct explaining the kinematic layout, expected mathematical behavior, and Featherstone reference
- [ ] **TST-04**: Edge case test coverage (near-zero mass inertia, singular joint configurations, near-singular, high-DOF chains n>10, non-identity rotations)

### PCC: Pinocchio C++ Comparison

- [ ] **PCC-01**: Pinocchio C++ adapter implementing `computeTorques()` via `pinocchio::rnea()` and `computeAccelerations()` via `pinocchio::aba()`
- [ ] **PCC-02**: CMake integration — `find_package(pinocchio REQUIRED)` with optional `SA_BUILD_PINOCCHIO_BENCHMARKS` guard (default OFF)
- [ ] **PCC-03**: Comparison benchmark executable running all test models through both SpatialAlgebra and Pinocchio adapters, producing per-joint numerical comparison
- [ ] **PCC-04**: Model-to-Pinocchio conversion utility translating Eigen-only test model definitions into Pinocchio Model+Data structures

### PCP: Pinocchio Python Harness

- [ ] **PCP-01**: Python comparison script using `pinocchio` Python bindings (`pip install pin`)
- [ ] **PCP-02**: Equivalent kinematic chain models defined in Python mirroring C++ test model library (3-link serial, branching Y, 2-link gravity)
- [ ] **PCP-03**: Round-trip RNEA(rnea(ABA(tau))) ≈ tau consistency check against SpatialAlgebra golden values
- [ ] **PCP-04**: Serialize comparison results to JSON for cross-language aggregation with C++ benchmark results

### RPT: Result Reporting

- [ ] **RPT-01**: Cross-library comparison tables (SpatialAlgebra C++ vs Pinocchio C++ vs Pinocchio Python) with per-model, per-joint error metrics
- [ ] **RPT-02**: Relative error analysis per joint, per solver, per test model — report max relative error, mean relative error, pass/fail per threshold
- [ ] **RPT-03**: CI-integrated regression tracking — comparison results captured as test assertions, pass/fail based on tolerance thresholds (1e-8 for kinematics, 1e-6 for dynamics)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Pinocchio derivative algorithms (ABA derivatives, RNEA derivatives) | Beyond numerical cross-validation scope; deferred to v2.0 |
| RBDL comparison benchmarks | Previously deferred from v1.2; Pinocchio takes priority |
| Real-time comparison profiling (Google Benchmark) | Profiling already exists in benchmarks/; cross-library timing is separate concern |
| URDF model loading in test library | Test models are programmatically constructed; URDF support would add large dependency |
| Floating-base robot models | SpatialAlgebra is fixed-base only; floating-base cross-validation deferred |

## Traceability

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

**Coverage:**
- v1.3 requirements: 20 total
- Mapped to phases: 20
- Unmapped: 0 ✓

---
*Requirements defined: 2026-06-17*
*Last updated: 2026-06-17 after initial definition*
