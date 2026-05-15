# Phase 1: Foundation Vectors - Context

**Gathered:** 2026-05-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Phase 1 delivers complete and verified SpatialVector, MotionVector, and ForceVector implementations. This includes fixing the crossMotion bug in MotionVector, upgrading tests to GTest, and implementing comprehensive verification using both textbook examples and property-based tests.

</domain>

<decisions>
## Implementation Decisions

### Test Framework
- **D-01:** Use GTest for Phase 1 tests instead of deferring to Phase 8 — provides better test reporting and assertions from the start

### Implementation Completeness
- **D-02:** Fix MotionVector::crossMotion to match SpatialVector base class formula: `[ω1×ω2; ω1×v2 + v1×ω2]` — current implementation `[ω1×ω2; v1×v2]` is incorrect per Featherstone

### API Design Patterns
- **D-03:** Keep minimal API — only essential constructors and operations (current state). No convenience factories (Zero(), Random()) or additional ergonomic operators in Phase 1

### Verification Strategy
- **D-04:** Use combined verification approach:
  - Unit tests with known numerical values from Featherstone textbook examples
  - Property-based tests verifying invariants (e.g., cross product anti-commutativity, distributivity)

### the agent's Discretion
- Test file structure and organization
- Specific test case names and grouping
- Order of implementation (fix bug first vs write tests first)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Spatial Algebra Formulation
- `include/SpatialVector.h` — Base class definition and cross product formulas
- `include/MotionVector.h` — Motion vector specialization (has bug in crossMotion)
- `include/ForceVector.h` — Force vector specialization
- `src/SpatialVector.cpp:42-46` — Correct crossMotion implementation: `[angular.cross(other.angular), linear.cross(other.angular) + angular.cross(other.linear)]`
- `src/MotionVector.cpp:42-46` — Bug location: uses `linear.cross(other.linear)` instead of correct formula

### Testing
- `tests/TestSpatialVector.cpp` — Current minimal test file (uses assert())
- `tests/TestPluckerTransform.cpp` — Example of GTest usage in project

### Mathematical Reference
- Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2 — Spatial vector algebra formulation

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **SpatialVector base class** — Fully implemented with correct crossMotion and crossForce formulas
- **ForceVector** — Implementation appears complete, no bugs detected
- **Eigen3 integration** — Vector3d type alias established for 3D vector operations

### Established Patterns
- **Inheritance pattern** — MotionVector and ForceVector extend SpatialVector
- **Type aliases** — `mv` for MotionVector, `fv` for ForceVector used throughout codebase
- **Doxygen documentation** — Every class and method has detailed @brief and @details
- **Include guards** — `#ifndef`/`#define`/`#endif` pattern (except LowerTriangular.h uses #pragma once)

### Integration Points
- **MotionVector::crossForce** — Used by PluckerTransform for coordinate transformations
- **ForceVector operations** — Used by RigidBodyInertia::apply() to compute wrenches from twists
- **Test executables** — CMakeLists.txt registers TestSpatialVector as separate test executable

</code_context>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches for test organization and implementation order.

</specifics>

<deferred>
## Deferred Ideas

### Reviewed Todos (not folded)
None — discussion stayed within phase scope

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-foundation-vectors*
*Context gathered: 2026-05-15*
