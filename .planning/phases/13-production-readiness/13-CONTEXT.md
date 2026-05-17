# Phase 13: Production Readiness — Context

**Gathered:** 2026-05-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Phase 13 fixes all P0/P1 critical bugs, fills testing gaps, adds gravity support, adds CI pipeline, and fixes code quality issues — making the library safe for real robotics use. This is an implementation-quality phase, not a new-feature phase.

</domain>

<decisions>
## Implementation Decisions

### Cross-Product Fix Strategy
- **D-01:** Unify force×force cross product to a single canonical implementation in `include/SpatialUtils.h` as a free function `cross(ForceVector, ForceVector)`. `SpatialVector::crossForce`, `ForceVector::crossForce`, and the existing free function must all delegate to this single implementation.
- **D-02:** Correct formula: `angular = τ1×τ2 + f1×f2`, `linear = τ1×f2 - τ2×f1` (matches Featherstone, anti-commutative).
- **D-03:** Test with mixed torque+force inputs (e.g., τ=(1,0,0), f=(0,1,0)) — triggers all 3 terms. Verify anti-commutativity property.
- **D-04:** Remove deprecated/useless overloads: `MotionVector::crossForce` (motion×motion has no physical meaning) and `ForceVector::crossMotion` (force×force with motion-cross formula is a category error). Do not mark deprecated — remove.

### ABA Multi-Link Fix
- **D-05:** Restructure `ForwardDynamics::inwardPass` into a clean two-phase approach:
  - Phase 1: Initialize all `Ia`/`pa` from rigid body inertia only (base→tip or separate init).
  - Phase 2: Accumulate child contributions into parent (tip→base) without re-initialization.
- **D-06:** Verify with known numerical acceleration values from Featherstone textbook examples (Chapter 7), not just "positive and finite" checks.

### Gravity Term Design
- **D-07:** Add gravity vector as optional parameter to `computeAccelerations()` and `computeTorques()` (defaults to zero vector for backward compatibility).
- **D-08:** Implement via base link acceleration: set `a₀ = -g` (gravity as spatial acceleration of base). Propagates through the outward pass automatically per Featherstone formulation.

### CI Pipeline
- **D-09:** Use GitHub Actions with a single workflow file.
- **D-10:** Matrix build: `ubuntu-latest` and `macos-latest`, each with `g++` and `clang++`.
- **D-11:** Add coverage tracking: compile with `--coverage`, generate gcov/lcov reports, upload to CodeCov.
- **D-12:** Tests run on every push and PR.

### NaN/Inf Validation
- **D-13:** Add debug-mode assertions (`#ifndef NDEBUG`) for NaN/Inf in: SpatialVector constructors, operator+, operator-, operator*, cross() methods, RigidBodyInertia::apply(), ArticulatedBodyInertia::apply().
- **D-14:** No production runtime overhead. Existing NaN/Inf checks at dynamics solver entry points remain as-is.

### Test Helper and Multi-Link Tests
- **D-15:** Fix `createIdentityInertia()` and `createDiagonalInertia()` in `tests/TestSpatialOperations.cpp` to return actual identity/diagonal matrices. Add assertions in existing tests to verify inertia matrix values (not just getMass/getCom).
- **D-16:** Add RNEA tests with known non-zero `qdot` values and expected `tau` values — exercises Coriolis/centrifugal computation.
- **D-17:** Multi-link ABA tests validate numerical acceleration values from Featherstone textbook examples.

### Eigen 5.x Compatibility
- **D-18:** Remove the version pin from `CMakeLists.txt` — change `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` to `find_package(Eigen3 REQUIRED NO_MODULE)`. Works with both Eigen 3.x and 5.x.

### Execution Ordering
- **D-19:** Fix bugs first (cross-product, ABI args, ABA, gravity), then add rigorous tests to verify fixes, then set up CI.

### ABI Constructor Arguments
- **D-20:** Leave to agent discretion — bug is well-documented in CONCERNS.md (CR-01). Fix is clear: swap first/third arguments in `operator+(RigidBodyInertia)`, add mass multiplier to coupling term.

### PluckerTransform auto Return Type
- **D-21:** Change `multiply()` return type from `auto` to `PluckerTransform` explicitly in the header. Ensure the definition is visible at call sites or move it to the header.

### SpatialOperations Unsafe Downcasts
- **D-22:** Change `crossProductMotion` and `crossProductForce` signatures to accept `const MotionVector&` and `const ForceVector&` directly, eliminating the `static_cast`.

### Global Namespace Pollution
- **D-23:** Move `using Vector3d` from global scope into `namespace SpatialAlgebra` in a single header. Remove the duplicate declarations from the other 5 headers.

### Include Guard Consistency
- **D-24:** Change `LowerTriangular.h` from `#pragma once` to `#ifndef LOWER_TRIANGULAR_H` / `#define LOWER_TRIANGULAR_H` / `#endif` to match the rest of the codebase.

### GTest Dependency Fallback
- **D-25:** Add `FetchContent` fallback for Google Test in `CMakeLists.txt` so the library builds even when GTest is not system-installed.

### LowerTriangular OpenMP
- **D-26:** Remove `#pragma omp parallel for collapse(2)` from `LowerTriangular` — it is silently ignored on macOS and provides no benefit for 3×3 matrices where thread spawning overhead exceeds the parallel gain.

### Empty Stub Source Files
- **D-27:** Remove `src/RigidBodyInertia.cpp` and `src/ArticulatedBodyInertia.cpp` (empty stubs). Update CMake to not expect them.

### Umbrella Header
- **D-28:** Create `include/SpatialAlgebra.h` that includes all public headers in documented dependency order.

### the agent's Discretion
- Specific test case values for multi-link reference tests (within Featherstone examples)
- Code cleanup order among the smaller items (PluckerTransform auto, downcasts, namespace, guards, GTest fallback, OpenMP, stubs, umbrella header)
- Namespace cleanup — which header hosts the canonical `using Vector3d` declaration

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Production Readiness Analysis
- `.planning/codebase/CONCERNS.md` — Full catalog of known bugs, tech debt, and concerns (CR-01 through CR-03, ABI transform errors, multi-link dynamics gaps)
- `.planning/codebase/TESTING.md` — Current testing patterns, framework usage, known gaps

### Mathematical Reference
- Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2 (spatial vector algebra), Chapter 7 (forward/inverse dynamics, gravity formulation).
- `include/SpatialUtils.h` — Location for canonical cross product implementation (D-01)
- `src/SpatialVector.cpp:48-55` — Current (buggy) crossForce implementation
- `src/ForceVector.cpp:40-46` — Current (buggy) crossForce override
- `include/SpatialUtils.h:115-125` — Current (buggy) free function cross

### ABA / Forward Dynamics
- `src/ForwardDynamics.cpp:54-107` — Inward pass with the child-inertia overwrite bug (CR-02)
- `include/ForwardDynamics.h` — Link struct and ForwardDynamics class interface
- `tests/TestForwardDynamics.cpp` — Current tests (single-link only, "positive and finite")

### RNEA / Inverse Dynamics
- `src/InverseDynamics.cpp` — Current implementation (no gravity support)
- `include/InverseDynamics.h` — InverseDynamicsLink and class interface
- `tests/TestInverseDynamics.cpp` — Current tests (all zero-velocity)

### ABI / Inertia
- `include/ArticulatedBodyInertia.h:150-155` — Swapped constructor arguments (CR-01)
- `include/ArticulatedBodyInertia.h:95-96` — Same bug, identified but unfixed

### Build System
- `CMakeLists.txt` — Eigen version pin (D-18), GTest dependency (D-25), OpenMP (D-26)

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **SpatialUtils.h free functions** — Existing `cross()` overload pattern to extend with correct force×force implementation
- **ForwardDynamics Link struct** — Contains all fields needed for two-phase restructure (Ia, pa, I, etc.)
- **Test patterns in TestPluckerTransform.cpp** — Well-structured GTest with property-based tests, good model for new tests

### Established Patterns
- **Exception-based error handling** — `std::invalid_argument`, `std::runtime_error` used across library
- **Debug-mode assertions** — `#ifndef NDEBUG` pattern already used in LowerTriangular bounds checking
- **Inline-heavy design** — Inertia classes are all-inline in headers; new gravity parameter follows that pattern
- **Algorithim three-phase pattern** — outward→inward→solve in both ABA and RNEA

### Integration Points
- **`computeAccelerations()` signature change** — Adding gravity parameter requires updating all call sites
- **`computeTorques()` signature change** — Same for inverse dynamics
- **CMakeLists.txt** — Multiple changes: version pin, GTest fallback, coverage flags, remove stub .cpp from glob or list

### Creative Options
- Gravity as base acceleration is minimal-change approach but RNEA needs base acceleration in outward pass too — same pattern works for both
- Two-phase ABA restructure can reuse the existing loop structure with an added initialization pass

</code_context>

<specifics>
## Specific Ideas

No specific external examples cited — open to standard Featherstone-based approaches for multi-link reference values.

</specifics>

<deferred>
## Deferred Ideas

### Future Concerns (not in Phase 13 scope)
- **Python bindings** — Pybind11 integration deferred to v2 milestone (per PROJECT.md)
- **Joint limit / singularity handling** — Damped least-squares regularization for near-singular ABA configurations
- **Branching tree ABA scaling** — Explicit child-link index lists for O(n) branching tree traversal
- **LowerTriangular inverse threshold** — Hardcoded 1e-15 threshold could be parameterized in a future phase
- **Condition number estimation** — Not needed for current use cases, deferred

None — discussion stayed within phase scope

</deferred>

---

*Phase: 13-production-readiness*
*Context gathered: 2026-05-17*
