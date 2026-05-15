# SpatialAlgebra Project

**Project Code:** SA  
**Project Title:** Spatial Vector Algebra Library for Robotics  
**Last Updated:** 2026-05-15

---

## What This Is

A C++17 library implementing spatial vector algebra for rigid body dynamics, following Featherstone's formulation. Provides 6D spatial vectors (twists and wrenches), Plücker coordinate transforms, and inertia representations for robotics simulation and control.

## Core Value

**Must Deliver:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

**Success Looks Like:** 
- All incomplete implementations finished (ArticulatedBodyInertia, RigidBodyInertia, SpatialOperations)
- Forward dynamics (Articulated Body Algorithm - ABA) implemented
- 100% test coverage for all mathematical operations
- Library ready for integration into robotics projects

## Context

**Brownfield Project:** Existing codebase with partial implementations.

**Current State:**
- Core classes implemented: SpatialVector, MotionVector, ForceVector, Rotation, PluckerTransform
- Incomplete: RigidBodyInertia, ArticulatedBodyInertia (stub tests, incomplete inverse transform)
- Missing: Forward dynamics (ABA), no Python integration
- Test coverage: 2/5 test files implemented, 3 empty stubs
- Standalone Python RNEA exists but not integrated with C++ library

**Tech Stack:**
- C++17 with Eigen3 for linear algebra
- Google Test for unit testing
- CMake build system
- Doxygen documentation

## Constraints

**Technical:**
- Must maintain Eigen3 compatibility (version 3.3+ or 5.x with testing)
- Preserve existing class interfaces (backward compatibility)
- Follow existing code conventions (Doxygen comments, type aliases)

**Timeline:** Flexible — quality-focused completion

**Budget:** N/A (open source library)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Complete all incomplete code before adding features | Foundation must be solid before extending | All stub files will be implemented |
| Add forward dynamics (ABA) as next feature | Natural complement to existing RNEA (inverse dynamics) | ABA will be implemented after completing existing code |
| Use fine-grained phases | Allows focused verification of each component | 8-12 phases with 5-10 plans each |
| Parallel execution | Independent components can be developed simultaneously | Faster iteration |
| Full verification workflow | Mathematical correctness is critical | Research, plan check, and verifier enabled |

## Out of Scope

- **Python bindings** — Python RNEA exists standalone; integration deferred to future milestone
- **Contact/collision handling** — Pure rigid body dynamics only
- **Inverse kinematics** — Focus on dynamics, not kinematics
- **Visualization tools** — Library-only, no rendering

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---

*Last updated: 2026-05-15 after initialization*
