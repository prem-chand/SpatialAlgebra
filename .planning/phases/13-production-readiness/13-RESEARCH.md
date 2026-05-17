# Phase 13: Production Readiness — Research

**Researched:** 2026-05-17
**Domain:** C++17 spatial algebra — bug fixing, testing, CI, code quality
**Confidence:** HIGH

## Summary

Phase 13 transforms the SpatialAlgebra library from a research prototype with known correctness issues into a production-ready library safe for real robotics use. The phase addresses 3 critical bugs (cross-product inconsistency, ABI constructor argument swap, ABA child-inertia overwrite), adds gravity support, sets up CI/CD, and repairs multiple code-quality issues.

**Current state:** 151/158 tests pass (95.6%), 7 known failures from v1.0 verification plus additional concerns documented in CONCERNS.md. The 7 failures break down as: 5 ABI transform tests (now fixed in a prior phase, currently 0 failures from PluckerTransform), and 2 dynamics consistency tests (still failing). Additional bugs exist that are masked by limited test coverage (cross-product with mixed inputs, ABA multi-link, test helper returning zero matrices).

**Primary recommendation:** Fix bugs first (cross-product → ABI args → ABA → gravity), then add rigorous tests to verify fixes, then set up CI. This ordering ensures that new tests validate correct behavior and CI catches regressions.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Cross-product correctness | Library (SpatialUtils.h) | — | Single canonical implementation in free functions, other callers delegate |
| ABA inward pass | Library (ForwardDynamics.cpp) | — | Algorithm implementation; no UI or external service involvement |
| Gravity term | Library (ForwardDynamics, InverseDynamics) | — | Added as parameter to existing API; propagates through algorithm passes |
| NaN/Inf guards | Library (core vector/inertia headers/cpp) | — | Debug-mode assertions in math operations; no runtime config |
| Test helper fixes | Tests | — | Test-local helper functions return correct values |
| Multi-link tests | Tests | — | GTest-based verification against known Featherstone values |
| CI pipeline | Infrastructure (GitHub Actions) | — | Build + test + coverage on push/PR |
| Eigen 5.x compat | Build (CMakeLists.txt) | — | Remove version pin |
| Code quality fixes | Library headers/source | — | Individual fixes per D-20 through D-28 |

## User Constraints (from CONTEXT.md)

### Locked Decisions

**D-01:** Unify force×force cross product to a single canonical implementation in `include/SpatialUtils.h` as a free function `cross(ForceVector, ForceVector)`. `SpatialVector::crossForce`, `ForceVector::crossForce`, and the existing free function must all delegate to this single implementation.

**D-02:** Correct formula: `angular = τ1×τ2 + f1×f2`, `linear = τ1×f2 - τ2×f1` (matches Featherstone, anti-commutative).

**D-03:** Test with mixed torque+force inputs (e.g., τ=(1,0,0), f=(0,1,0)) — triggers all 3 terms. Verify anti-commutativity property.

**D-04:** Remove deprecated/useless overloads: `MotionVector::crossForce` (motion×motion has no physical meaning) and `ForceVector::crossMotion` (force×force with motion-cross formula is a category error). Do not mark deprecated — remove.

**D-05:** Restructure `ForwardDynamics::inwardPass` into a clean two-phase approach: Phase 1 initializes Ia/pa from rigid body inertia only, Phase 2 accumulates child contributions into parent without re-initialization.

**D-06:** Verify with known numerical acceleration values from Featherstone textbook examples (Chapter 7), not just "positive and finite" checks.

**D-07:** Add gravity vector as optional parameter to `computeAccelerations()` and `computeTorques()` (defaults to zero vector for backward compatibility).

**D-08:** Implement via base link acceleration: set `a₀ = -g` (gravity as spatial acceleration of base). Propagates through the outward pass automatically per Featherstone formulation.

**D-09:** Use GitHub Actions with a single workflow file.

**D-10:** Matrix build: `ubuntu-latest` and `macos-latest`, each with `g++` and `clang++`.

**D-11:** Add coverage tracking: compile with `--coverage`, generate gcov/lcov reports, upload to CodeCov.

**D-12:** Tests run on every push and PR.

**D-13:** Add debug-mode assertions (`#ifndef NDEBUG`) for NaN/Inf in: SpatialVector constructors, operator+, operator-, operator*, cross() methods, RigidBodyInertia::apply(), ArticulatedBodyInertia::apply().

**D-14:** No production runtime overhead. Existing NaN/Inf checks at dynamics solver entry points remain as-is.

**D-15:** Fix `createIdentityInertia()` and `createDiagonalInertia()` in `tests/TestSpatialOperations.cpp` to return actual identity/diagonal matrices. Add assertions in existing tests to verify inertia matrix values (not just getMass/getCom).

**D-16:** Add RNEA tests with known non-zero `qdot` values and expected `tau` values — exercises Coriolis/centrifugal computation.

**D-17:** Multi-link ABA tests validate numerical acceleration values from Featherstone textbook examples.

**D-18:** Remove the version pin from `CMakeLists.txt` — change `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` to `find_package(Eigen3 REQUIRED NO_MODULE)`. Works with both Eigen 3.x and 5.x.

**D-19:** Fix bugs first (cross-product, ABI args, ABA, gravity), then add rigorous tests to verify fixes, then set up CI.

**D-20:** Fix ABI constructor argument swap — swap first/third arguments in `operator+(RigidBodyInertia)`, add mass multiplier to coupling term.

**D-21:** Change `multiply()` return type from `auto` to `PluckerTransform` explicitly in the header.

**D-22:** Change `crossProductMotion` and `crossProductForce` signatures to accept `const MotionVector&` and `const ForceVector&` directly.

**D-23:** Move `using Vector3d` from global scope into `namespace SpatialAlgebra` in a single header. Remove the duplicate from `LowerTriangular.h`.

**D-24:** Change `LowerTriangular.h` from `#pragma once` to `#ifndef LOWER_TRIANGULAR_H` / `#define LOWER_TRIANGULAR_H` / `#endif`.

**D-25:** Add `FetchContent` fallback for Google Test in `CMakeLists.txt`.

**D-26:** Remove `#pragma omp parallel for collapse(2)` from `LowerTriangular`.

**D-27:** Remove `src/RigidBodyInertia.cpp` and `src/ArticulatedBodyInertia.cpp` (empty stubs). Update CMake.

**D-28:** Create `include/SpatialAlgebra.h` umbrella header.

### the agent's Discretion
- Specific test case values for multi-link reference tests (within Featherstone examples)
- Code cleanup order among smaller items (D-21 through D-28)
- Namespace cleanup — which header hosts the canonical `using Vector3d` declaration

### Deferred Ideas (OUT OF SCOPE)
- Python bindings (pybind11)
- Joint limit / singularity handling
- Branching tree ABA scaling
- LowerTriangular inverse threshold parameterization
- Condition number estimation

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| VEC-01 | SpatialVector base class fully functional | Currently 26/26 tests pass. Cross-product bug (CR-03) affects correctness of force×force but doesn't cause test failures due to weak test coverage. Fix requires unifying 3 buggy implementations into one canonical free function in SpatialUtils.h. |
| UTL-03 | cross() products for spatial vectors | Buggy free function `cross(ForceVector,ForceVector)` in SpatialUtils.h:115-125 has correct angular term (`τ1×τ2 + f1×f2`) but missing the `-τ2×f1` term in the linear component. D-01 through D-04 specify the fix. |
| PLX-04 | Articulated body inertia transformation | 5 ABI tests were failing but have been fixed in a prior phase (PluckerTransform now passes all 27+9=36 tests). No additional work needed unless the fix regresses. |
| INR-01 | RigidBodyInertia construction | Currently passes. The ABI constructor argument swap (CR-01) affects `operator+(RigidBodyInertia)` — the inertia tensor and mass matrix arguments are swapped, and the coupling term lacks a mass multiplier. D-20 addresses this. |
| ABA-01 | Articulated Body Algorithm forward dynamics | CR-02 bug: inward pass re-initializes Ia/pa for each link, overwriting child contributions. For a 3-link chain (0→1→2): link 2's contribution is added to link 1, then link 1 re-initializes itself (destroying link 2's contribution). D-05/D-06 restructure to two-phase approach. |
| ABA-02 | ABA handles serial kinematic chains | Two-phase restructure (D-05) is required before multi-link works. Current tests only pass for single-link (no children to accumulate). |
| TST-07 | Integration tests for complete dynamics pipeline | 2 consistency tests fail (ThreeLinkSerialChain, BranchingYConfiguration). After fixing CR-02 and CR-01, these should pass. D-16/D-17 add more rigorous tests. |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Eigen3 | 3.3+ / 5.x | Linear algebra backend | Existing dependency, universally used in robotics C++. Version pin will be removed (D-18). |
| Google Test | latest (brew) | Unit test framework | Existing, used in all test files. D-25 adds FetchContent fallback for CI. |
| gcov/lcov | system | Code coverage | Standard GCC/Clang toolchain. D-11 specifies --coverage flags + CodeCov upload. |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| GitHub Actions | N/A | CI pipeline | CI config only (D-09 through D-12). No library dependency. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| gcov/lcov | CodeCov's bash uploader + kcov | gcov is free, works with both GCC and Clang. kcov requires Ruby. |
| FetchContent GTest | vcpkg / Conan | FetchContent is CMake-native, no extra package manager needed for CI. |

**Installation:**
```bash
# No new dependencies. Existing deps:
brew install eigen
brew install googletest
```

**Version verification:** Eigen 3.3+ and 5.x are verified to work. The `using Vector3d` at global scope is the only potential compatibility issue with newer Eigen versions — fix per D-23.

## Package Legitimacy Audit

> No external packages are being added in this phase. All changes are to existing code, tests, CMakeLists.txt, and new CI workflow YAML. The only new file is `include/SpatialAlgebra.h` (umbrella header) and `.github/workflows/ci.yml`. No `npm`, `pip`, or `cargo` package installation is required.

## Architecture Patterns

### System Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                        Phase 13 Changes                             │
│                                                                     │
│  ┌──────────────┐   ┌───────────────┐   ┌────────────────────────┐ │
│  │  Bug Fixes    │   │  New Features │   │  Code Quality / CI     │ │
│  │              │   │               │   │                        │ │
│  │• Cross-force │   │• Gravity      │   │• NaN/Inf assertions    │ │
│  │  unification │   │  parameter    │   │• Umbrella header       │ │
│  │• ABI args    │   │  (a₀ = -g)    │   │• Namespace cleanup     │ │
│  │• ABA inward  │   │               │   │• GTest FetchContent    │ │
│  │  pass split  │   │               │   │• OpenMP removal        │ │
│  └──────┬───────┘   └───────┬───────┘   • auto→explicit types    │ │
│         │                   │             • Remove empty stubs   │ │
│         │                   │             • Include guard fix    │ │
│         │                   │             • CI pipeline          │ │
│         │                   │             └─────────┬────────────┘ │
│         │                   │                       │              │
└─────────┼───────────────────┼───────────────────────┼──────────────┘
          │                   │                       │
          ▼                   ▼                       ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    Existing Library Surface                          │
│                                                                     │
│  SpatialUtils.h    ForwardDynamics.cpp   InverseDynamics.cpp        │
│  SpatialVector.cpp ForceVector.cpp       ArticulatedBodyInertia.h   │
│  LowerTriangular.h PluckerTransform.h/.cpp  SpatialOperations.cpp   │
│                                                                     │
│  Tests: TestSpatialOperations.cpp  TestForwardDynamics.cpp          │
│         TestInverseDynamics.cpp    TestDynamicsConsistency.cpp      │
│         TestPluckerTransform.cpp   TestSpatialVector.cpp            │
│         TestRigidBodyInertia.cpp   TestArticulatedBodyInertia.cpp   │
└──────────┬──────────────────────────────────────────────────────────┘
           │
           ▼
┌──────────────────────┐
│  CI Pipeline (new)   │
│                      │
│  GitHub Actions      │
│  ┌────────────────┐  │
│  │ ubuntu g++     │  │
│  │ ubuntu clang++ │  │
│  │ macos g++      │  │
│  │ macos clang++  │  │
│  │                │  │
│  │ gcov → CodeCov │  │
│  └────────────────┘  │
└──────────────────────┘
```

### Recommended Project Structure
```
SpatialAlgebra/
├── include/
│   ├── SpatialAlgebra.h          ★ NEW: Umbrella header (D-28)
│   ├── SpatialVector.h           ※ Fix: move using Vector3d into namespace (D-23)
│   ├── LowerTriangular.h         ※ Fix: include guard, remove OpenMP (D-24, D-26)
│   │                                Fix: add NaN/Inf assertions (D-13)
│   ├── SpatialUtils.h            ★ Fix: canonical cross() (D-01), add NaN assertions
│   ├── ForceVector.h             ※ Remove crossForce (partial, D-04)
│   ├── MotionVector.h            ※ Remove crossForce (D-04)
│   ├── PluckerTransform.h        ※ Fix: auto→explicit multiply() (D-21)
│   ├── ArticulatedBodyInertia.h  ★ Fix: operator+ arg order + mass multiplier (CR-01)
│   ├── RigidBodyInertia.h        ※ Add NaN assertions
│   ├── ForwardDynamics.h         ★ Add gravity param (D-07)
│   ├── InverseDynamics.h         ★ Add gravity param (D-07)
│   └── ... (Rotation.h, SpatialOperations.h, etc.)
├── src/
│   ├── ForwardDynamics.cpp       ★ Fix: two-phase inward pass (CR-02), gravity (D-08)
│   ├── InverseDynamics.cpp       ★ Add gravity (D-08)
│   ├── SpatialVector.cpp         ※ Fix: delegate crossForce to canonical
│   ├── ForceVector.cpp           ※ Fix: delegate crossForce, remove crossMotion
│   ├── MotionVector.cpp          ※ Remove crossForce
│   ├── PluckerTransform.cpp      ※ Fix: auto→explicit multiply()
│   ├── SpatialOperations.cpp     ※ Fix: unsafe downcasts (D-22)
│   ├── LowerTriangular.cpp       ※ Fix: remove OMP pragma
│   ├── RigidBodyInertia.cpp      ✗ REMOVE empty stub (D-27)
│   └── ArticulatedBodyInertia.cpp ✗ REMOVE empty stub (D-27)
├── tests/
│   ├── TestSpatialOperations.cpp ★ Fix: helper functions, add cross-force mixed tests
│   ├── TestForwardDynamics.cpp   ★ Add multi-link numerical validation (D-17)
│   ├── TestInverseDynamics.cpp   ★ Add non-zero velocity tests (D-16)
│   ├── TestDynamicsConsistency.cpp (modified to test gravity)
│   └── ... (add NaN propagation tests for SpatialVector, RigidBodyInertia, etc.)
├── CMakeLists.txt                ※ Fix: Eigen version pin, GTest FetchContent, coverage flags, remove stubs
├── .github/workflows/
│   └── ci.yml                    ★ NEW: GitHub Actions CI (D-09 through D-12)
└── ... (existing docs, examples, etc.)
```

### Pattern 1: Cross-Product Unification via Delegation
**What:** All 3 buggy force×force cross product implementations are replaced with delegation to a single canonical free function in SpatialUtils.h. This eliminates the inconsistency problem at its root.
**When to use:** Any time the same mathematical operation has multiple implementations that must produce identical results.
**Details:**
- Canonical implementation in `include/SpatialUtils.h` as free function `cross(ForceVector, ForceVector)`
- `SpatialVector::crossForce()` delegates: `return ForceVector(cross(ForceVector(*this), other));`
- `ForceVector::crossForce()` delegates: `return ForceVector(cross(*this, other));`
- The existing free function `cross(ForceVector, ForceVector)` in SpatialUtils.h is updated in-place to the correct formula

**Canonical formula (D-02):**
```cpp
// angular = τ1×τ2 + f1×f2
// linear  = τ1×f2 - τ2×f1
inline ForceVector cross(const ForceVector& v1, const ForceVector& v2) noexcept {
    const Vector3d& t1 = v1.getAngular();
    const Vector3d& f1 = v1.getLinear();
    const Vector3d& t2 = v2.getAngular();
    const Vector3d& f2 = v2.getLinear();
    return ForceVector(
        t1.cross(t2) + f1.cross(f2),      // angular: τ1×τ2 + f1×f2
        t1.cross(f2) - t2.cross(f1)        // linear:  τ1×f2 - τ2×f1
    );
}
```

### Pattern 2: Two-Phase ABA Inward Pass
**What:** The inward pass is split into a clear initialization phase (base→tip or separate init) and an accumulation/solve phase (tip→base).
**When to use:** Multi-link ABA where child inertias must be accumulated into parents without being overwritten.
**Details (D-05):**
```
Phase 1 (init, any order): For each link i:
    links[i].Ia = ArticulatedBodyInertia(I_i, skew(c_i)*m_i, m_i*I)
    links[i].pa = Ia·c_i + v_i × Ia·v_i + f_i

Phase 2 (accumulate, tip→base): For each link i from n-1 down to 0:
    If i has a parent p:
        links[p].Ia += X_i^{-1} · Ia_i · X_i^{-T}
        links[p].pa += X_i^{-T} · pa_i

Phase 3 (solve, tip→base): For each link i from n-1 down to 0:
    qddot_i = (tau_i - S_i^T · pa_i) / (S_i^T · Ia_i · S_i)
```

### Pattern 3: Gravity via Base Acceleration
**What:** Gravity is implemented by setting the base link's spatial acceleration to `a₀ = -g` (gravity vector as a spatial acceleration). This propagates through the outward pass automatically.
**When to use:** RNEA and ABA both need gravity. Same pattern works for both.
**Details (D-08):**
```
For RNEA outward pass:
    if (parent == -1):
        a₀ = S₀·q̈₀ - g  [gravity as base acceleration]
    else:
        aᵢ = Xᵢ·a_parent + Sᵢ·q̈ᵢ + vᵢ × Sᵢ·q̇ᵢ

For ABA outward pass:
    if (parent == -1):
        v₀ = S₀·q̇₀
        c₀ = -g  [gravity sets bias acceleration for base]
    else:
        vᵢ = Xᵢ·v_parent + Sᵢ·q̇ᵢ
        cᵢ = Xᵢ·c_parent + vᵢ × Sᵢ·q̇ᵢ
```

### Anti-Patterns to Avoid
- **Three mutally-inconsistent implementations of the same formula:** Never maintain separate implementations of the same mathematical operation. Always delegate to one canonical implementation.
- **Overwriting accumulated state during iteration:** In the ABA inward pass, initializing Ia/pa within the same loop that accumulates child contributions destroys previously accumulated values. Separate initialization from accumulation.
- **auto return type in header with definition in .cpp:** C++17 requires the definition to be visible for auto return type deduction. Headers with auto-declared functions whose definitions are in .cpp files cause link errors for external callers.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Test framework | Custom test runner | Google Test (GTest) | Already established. FetchContent fallback ensures CI can build tests without system-installed GTest. |
| Coverage tooling | Custom coverage | gcov + lcov + CodeCov | GCC/Clang native coverage generation. Upload to CodeCov for PR integration. |
| CI pipeline | Self-hosted CI | GitHub Actions | Free for public repos, matrix build support, CodeCov integration built-in, no infrastructure to maintain. |

**Key insight:** This phase is about fixing existing code and establishing quality infrastructure, not about adding new dependencies. The only new infrastructure is GitHub Actions (YAML configuration) and CodeCov (SaaS, token-based upload).

## Common Pitfalls

### Pitfall 1: Cross-Product Tests Pass Despite Bugs
**What goes wrong:** The 3 buggy force×force cross product implementations happen to produce correct results for the test cases used (pure torque or pure force inputs). Mixed torque+force inputs are required to trigger the inconsistencies.
**Why it happens:** With pure torque inputs (f=0), the `f1×f2` term in the angular component is zero, so "missing it" has no effect. With pure force inputs (τ=0), the `τ1×f2` and `τ2×f1` terms are zero, masking the missing term.
**How to avoid:** Always test with mixed torque+force inputs (e.g., τ=(1,0,0), f=(0,1,0)) to exercise all 3 terms in the formula.
**Warning signs:** Cross-product property tests (anti-commutativity, scalar multiplication) passing but a user reports wrong results for general wrenches.

### Pitfall 2: ABA Inward Pass Initialization Overwrites Accumulation
**What goes wrong:** When Ia and pa are initialized inside the tip→base loop, child contributions added in a previous iteration are overwritten by the next link's self-initialization.
**Why it happens:** The inward pass must first initialize all Ia/pa from the link's own rigid body inertia, then accumulate child contributions. Doing both in one loop or re-initializing after accumulation destroys the accumulated values.
**How to avoid:** Use explicit two-phase approach (D-05). Phase 1 initializes all links. Phase 2 accumulates without re-initializing.
**Warning signs:** Single-link tests pass but multi-link tests fail or show inconsistent results.

### Pitfall 3: Gravity Implementation Breaking Backward Compatibility
**What goes wrong:** Adding a gravity parameter to `computeAccelerations()` and `computeTorques()` changes the public API. Existing code that calls these without the parameter should continue to work (gravity defaults to zero).
**Why it happens:** A required (non-defaulted) parameter breaks all existing call sites.
**How to avoid:** Use an optional parameter with a default value of `Vector3d::Zero()` (D-07). This maintains backward compatibility.
**Warning signs:** Build breaks in existing examples or user code that calls these methods.

### Pitfall 4: Eigen 5.x CMake Compatibility
**What goes wrong:** `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` fails on systems with Eigen 5.x because Eigen 5 changed its CMake version compatibility range.
**Why it happens:** The version pin `3.3` restricts the acceptable range. Eigen 5.x reports a version outside this range.
**How to avoid:** Remove the version pin entirely (`find_package(Eigen3 REQUIRED NO_MODULE)`). The library uses only stable Eigen API (Matrix3d, cross products, etc.) that hasn't changed between 3.x and 5.x.
**Warning signs:** CMake configure fails with "Could not find a configuration file for package Eigen3" on Homebrew systems.

## Code Examples

### Canonical Cross-Product Implementation (SpatialUtils.h)
```cpp
// Source: D-02 derived from Featherstone (2008) Chapter 2
inline ForceVector cross(const ForceVector& v1, const ForceVector& v2) noexcept {
    const Vector3d& t1 = v1.getAngular();
    const Vector3d& f1 = v1.getLinear();
    const Vector3d& t2 = v2.getAngular();
    const Vector3d& f2 = v2.getLinear();
    return ForceVector(
        t1.cross(t2) + f1.cross(f2),   // angular = τ1×τ2 + f1×f2
        t1.cross(f2) - t2.cross(f1)    // linear  = τ1×f2 - τ2×f1
    );
}
```

### Delegation Pattern (SpatialVector.cpp)
```cpp
SpatialVector SpatialVector::crossForce(const SpatialVector &other) const {
    // Delegate to canonical free function in SpatialUtils.h
    return ForceVector(cross(ForceVector(*this), ForceVector(other)));
}
```

### Two-Phase ABA Inward Pass (ForwardDynamics.cpp)
```cpp
void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau) {
    // Phase 1: Initialize Ia and pa from rigid body inertia (all links)
    for (int i = 0; i < static_cast<int>(links.size()); i++) {
        double mass = links[i].I.getMass();
        Vector3d com = links[i].I.getCom();
        links[i].Ia = ArticulatedBodyInertia(
            links[i].I.getInertiaMatrixLT(),
            skew(com) * mass,
            lt::Identity(3) * mass
        );
        MotionVector IaV = links[i].Ia.apply(links[i].v);
        links[i].pa = links[i].Ia.apply(links[i].c)
                    + cross(links[i].v, IaV)
                    + links[i].f;
    }

    // Phase 2: Accumulate child contributions (tip→base, no re-init)
    for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--) {
        int parent = links[i].parent;
        if (parent != -1) {
            ArticulatedBodyInertia IaTransformed =
                links[i].X.invtformABI(links[i].Ia);
            links[parent].Ia = links[parent].Ia + IaTransformed;

            ForceVector paTransformed =
                links[i].X.inverseTransformForce(links[i].pa);
            links[parent].pa = ForceVector(
                links[parent].pa.getAngular() + paTransformed.getAngular(),
                links[parent].pa.getLinear() + paTransformed.getLinear()
            );
        }
    }

    // Phase 3: Solve for joint accelerations (tip→base, existing code)
    // ... unchanged from current implementation
}
```

### Gravity as Base Acceleration in RNEA (InverseDynamics.cpp)
```cpp
void InverseDynamics::outwardPass(const Vector3d& gravity) {
    for (int i = 0; i < static_cast<int>(links.size()); i++) {
        int parent = links[i].parent;
        if (parent == -1) {
            links[i].v = links[i].S * links[i].qdot;
            // Gravity: a₀ = S₀·q̈₀ - g (spatial acceleration of base)
            links[i].a = links[i].S * links[i].qddot
                       - MotionVector(Vector3d::Zero(), gravity);
        } else {
            MotionVector vParent = links[parent].v;
            links[i].v = links[i].X.transformMotion(vParent)
                       + links[i].S * links[i].qdot;
            MotionVector aParent = links[parent].a;
            links[i].a = links[i].X.transformMotion(aParent)
                       + links[i].S * links[i].qddot
                       + cross(links[i].v, links[i].S) * links[i].qdot;
        }
    }
}
```

### NaN/Inf Debug Assertion Pattern
```cpp
// Source: D-13, D-14
SpatialVector::SpatialVector(const Vector3d &a, const Vector3d &l)
    : angular(a), linear(l)
{
#ifndef NDEBUG
    if (a.hasNaN() || l.hasNaN() || 
        a.array().isInf().any() || l.array().isInf().any()) {
        std::cerr << "WARNING: NaN or Inf detected in SpatialVector constructor\n";
    }
#endif
}
```

### Featherstone Multi-Link Reference Values (for D-06, D-17)
The Featherstone textbook (Chapter 7) provides worked examples for 2-link and 3-link serial chains. Specific numerical values allow verifying that the ABA produces correct accelerations. The test should use the simple chain:
- 3 links with identity masses, transforms along X, Z-axis revolute joints
- Known torque inputs → expected acceleration outputs
- Reference: Featherstone (2008) Section 7.2, Algorithm 7.3

A simple verification approach: `ABA(RNEA(qddot)) ≈ qddot` for multi-link chains. After fixing CR-02 (ABA inward pass) and CR-01 (ABI args), the round-trip should be exact.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| CMake `find_package` with Eigen version pin | Version-pin removed | This phase (D-18) | Builds with Eigen 5.x without workarounds |
| No CI pipeline | GitHub Actions with 4× matrix | This phase (D-09) | Every push builds and tests on ubuntu/macos × g++/clang++ |
| No coverage tracking | gcov + CodeCov upload | This phase (D-11) | Coverage visible in PRs, unknown areas exposed |
| No NaN/Inf guards in core ops | Debug-mode assertions | This phase (D-13) | Catches NaN propagation early in development |
| No umbrella header | `include/SpatialAlgebra.h` | This phase (D-28) | Users include one header instead of 12 |
| System-installed GTest required | FetchContent fallback | This phase (D-25) | Builds on bare CI runners without brew |

**Deprecated/outdated:**
- `MotionVector::crossForce()` — physically meaningless, being removed (D-04)
- `ForceVector::crossMotion()` — category error, being removed (D-04)
- `#pragma omp parallel for` in LowerTriangular — silently ignored on macOS, overhead > benefit for 3×3 matrices (D-26)
- Empty stub `.cpp` files — waste compilation time (D-27)

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | ABI transform (tformABI/invtformABI) 5 failing tests have already been fixed in a prior phase | Phase Requirements | LOW — current TestPluckerTransform passes. Verify during execution. |
| A2 | Eigen API used is compatible between 3.4.x and 5.x | Standard Stack | LOW — only basic MatrixXd/Vector3d/cross operations used. Verify during CI testing on both. |
| A3 | Featherstone Chapter 7 provides reference numerical values for multi-link verification | Code Examples | MEDIUM — if not, we may need to derive them manually or use cross-validation (ABA↔RNEA round-trip) as the test oracle. |
| A4 | CodeCov service is available and the token can be obtained for this repo | Architecture Patterns | LOW — coverage upload can be skipped if CodeCov not available; CI still runs tests. |

## Open Questions (RESOLVED)

1. **What are the exact Featherstone textbook example values for multi-link verification?**
   - What we know: Chapter 7 of Featherstone (2008) provides Algorithm 7.3 with worked examples
   - What's unclear: The specific numerical values for the 3-link chain with unit masses
   - Recommendation: For "validation against known values", derive from round-trip consistency (ABA(RNEA(qddot)) ≈ qddot) which is the mathematically correct relationship. If exact Featherstone values are available, use those; otherwise, round-trip serves as an equally valid test oracle.

2. **Which header should host the canonical `using Vector3d` declaration?**
   - What we know: Must be in `namespace SpatialAlgebra`, not global scope. Currently in SpatialVector.h:44 and LowerTriangular.h:55 (both global).
   - Recommendation: SpatialVector.h (it's the base class, included by everything). Move the declaration into the namespace and remove the duplicate from LowerTriangular.h.

3. **Should coverage upload use CodeCov's bash uploader or GitHub Action?**
   - What we know: CodeCov provides both `codecov/codecov-action@v4` and a bash uploader
   - Recommendation: Use `codecov/codecov-action@v4` — tracks build, fail-on-error, integrates with PR checks. The bash uploader is simpler but offers less integration.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| CMake | Build system | ✓ | 4.2.0 | — |
| Apple Clang | C++ compilation | ✓ | 17.0.0 | — |
| GNU G++ (via Homebrew) | CI matrix (linux) | ✗ (macOS) | — | Use Apple Clang on macOS; CI will provide g++ on ubuntu |
| Eigen3 | Library build | ✓ | 3.4.0_1 (brew) | — |
| Google Test | Test build | ✓ | brew-installed | FetchContent fallback (D-25) |
| gcov/lcov | Coverage | ✓ | system | — |
| Docker | CI-local testing | ✓ | 24.0.7 | Not required — CI runs natively |

**Missing dependencies with no fallback:** None

**Missing dependencies with fallback:** Google Test — FetchContent fallback will be added (D-25)

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Google Test (GTest) |
| Config file | none — GTest test executables defined in CMakeLists.txt |
| Quick run command | `cmake --build build && cd build && ctest --output-on-failure` |
| Full suite command | `cmake --build build && cd build && ctest --output-on-failure` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| VEC-01 | Cross-force correctness with mixed inputs | unit | `build/TestSpatialVector` | ✅ existing, update |
| VEC-01 | NaN/Inf guard in SpatialVector | unit | `build/TestSpatialVector` | ➕ new tests |
| UTL-03 | Canonical cross() anti-commutative | unit | `build/TestSpatialUtils` | ✅ existing, update |
| PLX-04 | ABI transform tests pass | unit | `build/TestPluckerTransform` | ✅ already passes |
| INR-01 | ABI+RBI operator+ correct args | unit | `build/TestArticulatedBodyInertia` | ✅ existing |
| ABA-01 | Multi-link ABA numerical check | integration | `build/TestForwardDynamics` | ➕ new tests |
| ABA-02 | 3-link serial chain dynamics | integration | `build/TestForwardDynamics` | ➕ new tests |
| ABA-02 | ABA with gravity | integration | `build/TestForwardDynamics` | ➕ new tests |
| TST-07 | RNEA↔ABA round-trip multi-link | integration | `build/TestDynamicsConsistency` | ✅ existing, fix |
| TST-07 | RNEA with non-zero velocity | integration | `build/TestInverseDynamics` | ➕ new tests |
| TST-07 | RNEA with gravity | integration | `build/TestInverseDynamics` | ➕ new tests |
| — | Test helper returns correct inertia | unit | `build/TestSpatialOperations` | ✅ existing, fix |

### Sampling Rate
- **Per task commit:** `cmake --build build && cd build && ctest --output-on-failure`
- **Per wave merge:** Full suite green
- **Phase gate:** Full suite green, all 11 test executables passing, CI workflow verified

### Wave 0 Gaps
- [ ] Enhanced NaN test cases in TestSpatialVector.cpp / TestRigidBodyInertia.cpp
- [ ] Cross-force mixed-input tests in TestSpatialUtils.cpp
- [ ] Multi-link ABA reference-value tests in TestForwardDynamics.cpp
- [ ] Non-zero-velocity RNEA tests in TestInverseDynamics.cpp
- [ ] Gravity propagation tests in TestForwardDynamics.cpp / TestInverseDynamics.cpp

## Security Domain

> **security_enforcement:** false — No network-facing, auth, or data-access concerns in this C++ spatial algebra library.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V5 Input Validation | partial | NaN/Inf checks in debug mode (D-13, D-14) |
| V6 Cryptography | no | Not applicable |

### Known Threat Patterns for {stack}
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| NaN/Inf propagation → unsafe robot motion | Denial of Service | Debug-mode assertions on core ops (D-13) + existing entry-point checks (D-14) |

The primary security concern is numerical safety — NaN/Inf propagation causing silent failures in downstream applications (e.g., robot controllers making unsafe movements). The fix (D-13/D-14) is debug-mode assertions that catch propagation early without production overhead.

## Sources

### Primary (HIGH confidence)
- Codebase examination of all 11 `include/*.h` and 12 `src/*.cpp` files — verified current buggy implementations
- Current test run output (156/158 passing) — confirmed 2 failing dynamics consistency tests
- v1.0-VERIFICATION.md — documented 7 known failures (5 ABI + 2 consistency)
- CONTEXT.md — 28 locked decisions from discuss-phase
- CONCERNS.md — CR-01 through CR-03 detailed analysis with file:line references
- TESTING.md — GTest patterns, test structure documentation

### Secondary (MEDIUM confidence)
- [CITED: Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2, 7] — Cross-product formulas, ABA algorithm, gravity via base acceleration
- [CITED: CodeCov GitHub Action docs (https://github.com/codecov/codecov-action)] — Coverage upload pattern
- [CITED: GitHub Actions matrix build docs (https://docs.github.com/en/actions/using-jobs/using-a-matrix-for-your-jobs)] — CI matrix strategy

### Tertiary (LOW confidence)
- No tertiary sources used. All findings verified against existing codebase or documented decisions.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — No new dependencies added
- Architecture: HIGH — All 28 decisions locked in CONTEXT.md
- Pitfalls: HIGH — Verified against actual codebase state and known failure patterns
- CI patterns: MEDIUM — Standard GitHub Actions but not tested on this repo yet

**Research date:** 2026-05-17
**Valid until:** 2026-06-17 (30 days — stable environment)
