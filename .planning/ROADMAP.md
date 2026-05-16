# SpatialAlgebra Roadmap

**Version:** 1.0  
**Last Updated:** 2026-05-15  
**Granularity:** Fine

---

## Phases

- [ ] **Phase 1: Foundation Vectors** - SpatialVector base class and motion/force vectors
- [x] **Phase 2: Rotation & Math** - Rotation matrix operations and conversions (completed 2026-05-15)
- [x] **Phase 3: Packed Matrix** - LowerTriangular matrix with packed storage (completed 2026-05-15)
- [ ] **Phase 4: Spatial Utilities** - Cross product, dot product, skew operators
- [ ] **Phase 5: Plücker Transforms** - 6D coordinate transformations
- [x] **Phase 6: Inertia Properties** - RigidBodyInertia and ArticulatedBodyInertia (completed 2026-05-16)
- [ ] **Phase 7: Forward Dynamics** - Articulated Body Algorithm implementation
- [ ] **Phase 8: Test Infrastructure** - Complete GTest coverage for all classes
- [ ] **Phase 9: Integration Tests** - End-to-end dynamics pipeline verification
- [ ] **Phase 10: Documentation** - README, examples, Doxygen

---

## Phase Details

### Phase 1: Foundation Vectors
**Goal**: Users can create and manipulate 6D spatial vectors (twists and wrenches)
**Depends on**: Nothing
**Requirements**: VEC-01, VEC-02, VEC-03, VEC-04
**Success Criteria** (what must be TRUE):
  1. User can create SpatialVector with angular and linear components
  2. User can create MotionVector (twist) and access both components
  3. User can create ForceVector (wrench) and access both components
  4. User can add, subtract, and scale spatial vectors
**Plans**: 4 plans

Plans:
- [x] 01-01-PLAN.md — GTest test scaffolding and comprehensive test suite
- [ ] 01-02-PLAN.md — Fix MotionVector::crossMotion bug per Featherstone
- [x] 01-03-PLAN.md — Complete ForceVector and SpatialVector::crossForce
- [ ] 01-04-PLAN.md — Property-based tests and textbook examples

### Phase 2: Rotation & Math
**Goal**: Users can perform rotation operations and convert between representations
**Depends on**: Phase 1
**Requirements**: ROT-01, ROT-02, ROT-03, ROT-04
**Success Criteria** (what must be TRUE):
  1. User can multiply, transpose, and invert rotation matrices
  2. User can convert angle-axis representation to rotation matrix
  3. User can convert quaternion to rotation matrix
  4. Rotation matrices maintain orthogonality through operations
**Plans**: 1 plan

Plans:
- [x] 02-01-PLAN.md — Comprehensive GTest test suite for Rotation class

### Phase 3: Packed Matrix
**Goal**: Users can use memory-efficient lower triangular matrices
**Depends on**: Phase 1
**Requirements**: LTR-01, LTR-02, LTR-03, LTR-04
**Success Criteria** (what must be TRUE):
  1. User can create LowerTriangular matrix with correct packed storage indexing
  2. User can multiply LowerTriangular with dense matrices
  3. User can multiply LowerTriangular with vectors
  4. User can compute transpose and inverse of LowerTriangular matrix
**Plans**: 2 plans

Plans:
- [x] 03-01-PLAN.md — Implement missing LowerTriangular methods (operator*, inverse)
- [x] 03-02-PLAN.md — Create GTest test suite for LowerTriangular class

### Phase 4: Spatial Utilities
**Goal**: Users can compute cross products and dot products for spatial vectors
**Depends on**: Phase 1
**Requirements**: UTL-01, UTL-02, UTL-03, UTL-04
**Success Criteria** (what must be TRUE):
  1. User can create 3x3 skew-symmetric matrix from Vector3d
  2. User can compute dot products for spatial vectors
  3. User can compute cross products for spatial vectors
  4. User can use SpatialOperations utility class for common operations
**Plans**: 3 plans

Plans:
- [ ] 04-01-PLAN.md — Complete SpatialUtils.h free functions (skew, dot, cross)
- [ ] 04-02-PLAN.md — Implement SpatialOperations static class methods
- [ ] 04-03-PLAN.md — Create comprehensive GTest test suite

### Phase 5: Plücker Transforms
**Goal**: Users can transform spatial vectors and inertias between coordinate frames
**Depends on**: Phase 1, Phase 2, Phase 3
**Requirements**: PLX-01, PLX-02, PLX-03, PLX-04, PLX-05, PLX-06
**Success Criteria** (what must be TRUE):
  1. User can transform motion vectors between frames using `transformMotion()`
  2. User can transform force vectors between frames using `transformForce()`
  3. User can transform rigid body inertia using `tformRBI()`
  4. User can transform articulated body inertia using `tformABI()`
  5. User can compute inverse Plücker transform
  6. User can compute inverse articulated body inertia transform using `invtformABI()`
**Plans**: 3 plans

Plans:
- [ ] 05-01-PLAN.md — Verify and test transformMotion/transformForce (PLX-01, PLX-02)
- [ ] 05-02-PLAN.md — Verify and test tformRBI/invtformRBI/inverse (PLX-03, PLX-05)
- [ ] 05-03-PLAN.md — Implement invtformABI and test tformABI (PLX-04, PLX-06)

### Phase 6: Inertia Properties
**Goal**: Users can define and manipulate rigid body and articulated body inertias
**Depends on**: Phase 1, Phase 3
**Requirements**: INR-01, INR-02, INR-03, INR-04
**Success Criteria** (what must be TRUE):
  1. User can create RigidBodyInertia with mass, COM, and inertia tensor
  2. User can apply MotionVector to RigidBodyInertia to get ForceVector
  3. User can create ArticulatedBodyInertia with full parameterization
  4. User can apply MotionVector to ArticulatedBodyInertia to get ForceVector
**Plans**: 3 plans

Plans:
- [x] 06-01-PLAN.md — Implement ArticulatedBodyInertia::apply() and print() methods
- [x] 06-02-PLAN.md — Create GTest test suite for RigidBodyInertia
- [x] 06-03-PLAN.md — Create GTest test suite for ArticulatedBodyInertia

### Phase 7: Forward Dynamics
**Goal**: Users can compute forward dynamics for articulated rigid body systems
**Depends on**: Phase 5, Phase 6, Phase 4
**Requirements**: ABA-01, ABA-02, ABA-03, ABA-04
**Success Criteria** (what must be TRUE):
  1. User can run Articulated Body Algorithm for forward dynamics
  2. User can compute accelerations for serial kinematic chains
  3. User can compute accelerations for branching kinematic trees
  4. ABA correctly uses PluckerTransform operations for coordinate transformations
**Plans**: 2 plans

Plans:
- [ ] 07-01-PLAN.md — Implement ABA algorithm (outward/inward pass, computeAccelerations)
- [ ] 07-02-PLAN.md — Create comprehensive GTest test suite for ABA

### Phase 8: Test Infrastructure
**Goal**: All core classes have comprehensive GTest test suites
**Depends on**: Phase 1, Phase 2, Phase 3, Phase 4, Phase 5, Phase 6
**Requirements**: TST-01, TST-02, TST-03, TST-04, TST-05, TST-06
**Success Criteria** (what must be TRUE):
  1. All empty test stub files are implemented with GTest
  2. RigidBodyInertia has 100% test coverage for all operations
  3. ArticulatedBodyInertia has 100% test coverage for all operations
  4. SpatialOperations utilities have 100% test coverage
  5. LowerTriangular has 100% test coverage for all operations
  6. All PluckerTransform methods have tests
**Plans**: TBD

### Phase 9: Integration Tests
**Goal**: Complete dynamics pipeline works end-to-end
**Depends on**: Phase 7, Phase 8
**Requirements**: TST-07
**Success Criteria** (what must be TRUE):
  1. Integration tests verify complete dynamics pipeline (RNEA + ABA)
  2. Tests verify consistency between inverse and forward dynamics
  3. Tests pass for multi-body systems with multiple links
**Plans**: TBD

### Phase 10: Documentation
**Goal**: Users can learn and use the library from documentation
**Depends on**: Phase 1, Phase 2, Phase 3, Phase 4, Phase 5, Phase 6, Phase 7
**Requirements**: DOC-01, DOC-02, DOC-03
**Success Criteria** (what must be TRUE):
  1. README.md contains build instructions and basic usage examples
  2. Doxygen documentation is generated and up to date
  3. Usage examples demonstrate core operations (vectors, transforms, inertia, ABA)
**Plans**: TBD

---

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation Vectors | 1/4 | In Progress|  |
| 2. Rotation & Math | 1/1 | Complete   | 2026-05-15 |
| 3. Packed Matrix | 0/4 | Not started | - |
| 4. Spatial Utilities | 0/4 | Not started | - |
| 5. Plücker Transforms | 0/6 | Not started | - |
| 6. Inertia Properties | 0/4 | Not started | - |
| 7. Forward Dynamics | 0/4 | Not started | - |
| 8. Test Infrastructure | 0/6 | Not started | - |
| 9. Integration Tests | 0/3 | Not started | - |
| 10. Documentation | 0/3 | Not started | - |

---

## Requirement Coverage

**Total v1 requirements:** 37  
**Mapped:** 37/37 ✓

| Requirement | Phase | Status |
|-------------|-------|--------|
| VEC-01 | Phase 1 | Pending |
| VEC-02 | Phase 1 | Pending |
| VEC-03 | Phase 1 | Pending |
| VEC-04 | Phase 1 | Pending |
| ROT-01 | Phase 2 | Pending |
| ROT-02 | Phase 2 | Pending |
| ROT-03 | Phase 2 | Pending |
| ROT-04 | Phase 2 | Pending |
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
| TST-07 | Phase 9 | Pending |
| DOC-01 | Phase 10 | Pending |
| DOC-02 | Phase 10 | Pending |
| DOC-03 | Phase 10 | Pending |

---

*Roadmap created: 2026-05-15*
