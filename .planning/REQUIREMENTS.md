# Requirements

**Version:** 1.0  
**Last Updated:** 2026-05-15

---

## v1 Requirements (Active)

### Core Vector Operations

- [x] **VEC-01**: SpatialVector base class fully functional with angular/linear components
- [x] **VEC-02**: MotionVector (twist) operations complete and tested
- [x] **VEC-03**: ForceVector (wrench) operations complete and tested
- [x] **VEC-04**: Vector arithmetic (add, subtract, scale) for all spatial vectors

### Rotation Operations

- [x] **ROT-01**: Rotation matrix operations (multiply, transpose, inverse)
- [x] **ROT-02**: Angle-axis to rotation matrix conversion
- [x] **ROT-03**: Quaternion to rotation matrix conversion
- [x] **ROT-04**: Rotation matrix orthogonality preservation

### Plücker Transform Operations

- [ ] **PLX-01**: Motion vector transformation `transformMotion()` complete
- [ ] **PLX-02**: Force vector transformation `transformForce()` complete
- [ ] **PLX-03**: Rigid body inertia transformation `tformRBI()` complete
- [ ] **PLX-04**: Articulated body inertia transformation `tformABI()` complete
- [ ] **PLX-05**: Inverse transform `inverse()` complete and verified
- [ ] **PLX-06**: Inverse articulated body inertia transform `invtformABI()` implemented

### Inertia Operations

- [ ] **INR-01**: RigidBodyInertia construction with mass, COM, inertia tensor
- [ ] **INR-02**: RigidBodyInertia `apply(MotionVector)` returns correct ForceVector
- [ ] **INR-03**: ArticulatedBodyInertia construction complete
- [ ] **INR-04**: ArticulatedBodyInertia `apply(MotionVector)` returns correct ForceVector

### Lower Triangular Matrix

- [ ] **LTR-01**: Packed storage indexing correct for all operations
- [ ] **LTR-02**: Matrix multiplication with dense matrices
- [ ] **LTR-03**: Matrix multiplication with vectors
- [ ] **LTR-04**: Transpose and inverse operations

### Spatial Utilities

- [ ] **UTL-01**: `skew()` operator for cross product matrix
- [ ] **UTL-02**: `dot()` products for spatial vectors
- [ ] **UTL-03**: `cross()` products for spatial vectors
- [ ] **UTL-04**: SpatialOperations utility class complete

### Forward Dynamics (ABA)

- [ ] **ABA-01**: Articulated Body Algorithm forward dynamics implementation
- [ ] **ABA-02**: ABA handles serial kinematic chains
- [ ] **ABA-03**: ABA handles branching kinematic trees
- [ ] **ABA-04**: ABA integrates with existing PluckerTransform operations

### Test Coverage

- [ ] **TST-01**: All test stub files implemented with GTest
- [ ] **TST-02**: 100% coverage of RigidBodyInertia operations
- [ ] **TST-03**: 100% coverage of ArticulatedBodyInertia operations
- [ ] **TST-04**: 100% coverage of SpatialOperations utilities
- [ ] **TST-05**: 100% coverage of LowerTriangular operations
- [ ] **TST-06**: All PluckerTransform methods tested
- [x] **TST-07**: Integration tests for complete dynamics pipeline

### Documentation

- [ ] **DOC-01**: README.md with build instructions and basic usage
- [ ] **DOC-02**: Doxygen documentation up to date
- [ ] **DOC-03**: Usage examples for core operations

---

## v2 Requirements (Deferred)

### Python Integration

- [ ] **PYB-01**: pybind11 bindings for C++ library
- [ ] **PYB-02**: Python package with pip install support
- [ ] **PYB-03**: Integrate existing RNEA Python implementation with C++ backend

### Advanced Features

- [ ] **JAC-01**: Jacobian computation for kinematic chains
- [ ] **IK-01**: Inverse kinematics solver
- [ ] **CNT-01**: Contact and collision handling

### Performance

- [ ] **PERF-01**: Benchmark suite for core operations
- [ ] **PERF-02**: OpenMP optimization for parallel operations
- [ ] **PERF-03**: Memory layout optimization

---

## Out of Scope

- **Visualization tools** — Library-only, no rendering components
- **Contact/collision handling** — Pure rigid body dynamics; contact physics deferred
- **Inverse kinematics** — Focus on dynamics algorithms (RNEA, ABA), not kinematics
- **Python bindings (v1)** — Standalone Python RNEA exists; C++ integration deferred to v2
- **Serialization** — No save/load for configurations or state

---

## Traceability

*This section is populated by the roadmap. Each phase will map to the requirements it delivers.*

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
| LTR-01 | Phase 3 | Pending |
| LTR-02 | Phase 3 | Pending |
| LTR-03 | Phase 3 | Pending |
| LTR-04 | Phase 3 | Pending |
| UTL-01 | Phase 4 | Pending |
| UTL-02 | Phase 4 | Pending |
| UTL-03 | Phase 4 | Pending |
| UTL-04 | Phase 4 | Pending |
| PLX-01 | Phase 5 | Pending |
| PLX-02 | Phase 5 | Pending |
| PLX-03 | Phase 5 | Pending |
| PLX-04 | Phase 5 | Pending |
| PLX-05 | Phase 5 | Pending |
| PLX-06 | Phase 5 | Pending |
| INR-01 | Phase 6 | Pending |
| INR-02 | Phase 6 | Pending |
| INR-03 | Phase 6 | Pending |
| INR-04 | Phase 6 | Pending |
| ABA-01 | Phase 7 | Pending |
| ABA-02 | Phase 7 | Pending |
| ABA-03 | Phase 7 | Pending |
| ABA-04 | Phase 7 | Pending |
| TST-01 | Phase 8 | Pending |
| TST-02 | Phase 8 | Pending |
| TST-03 | Phase 8 | Pending |
| TST-04 | Phase 8 | Pending |
| TST-05 | Phase 8 | Pending |
| TST-06 | Phase 8 | Pending |
| TST-07 | Phase 9 | Complete |
| DOC-01 | Phase 10 | Pending |
| DOC-02 | Phase 10 | Pending |
| DOC-03 | Phase 10 | Pending |

---

*Requirements defined: 2026-05-15*
