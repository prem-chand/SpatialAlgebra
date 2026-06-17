# Phase 20: Test Model Library - Research

**Researched:** 2026-06-17
**Domain:** Header-only C++ library architecture, CMake INTERFACE targets, robotics solver abstraction, Eigen3 workflows
**Confidence:** HIGH

## Summary

Phase 20 creates a standalone, header-only test model library (`tests/test-models/`) with zero SpatialAlgebra dependency. The library defines pure-data kinematic chain descriptions (`RobotModel` + `JointSpec`) using only Eigen3 types, an abstract solver interface (`RobotSolver`), and a SpatialAlgebra adapter that bridges the two.

The key architectural challenge is enforcing the zero-dependency guarantee: `robot_model.h` and `robot_solver.h` must never transitively include any SpatialAlgebra header. This is achieved through an INTERFACE CMake target (`test_models`) that links only `Eigen3::Eigen`, with a separate compiled library (`sa_test_adapter`) that links `SpatialAlgebra` and implements the adapter.

The research confirms all 11 test domains can be represented, the JointSpec design is compatible with existing SA types, and the CMake strategy creates a clean compilation firewall. The library is designed to be forward-compatible with Phase 22 (Pinocchio adapter) and Phase 23 (Python model mirrors).

**Primary recommendation:** Implement as an INTERFACE CMake library with 11 `chains/` headers (one per test domain), with `sa_adapter.cpp` as the sole file that includes SpatialAlgebra headers. The existing test files are NOT modified in Phase 20 (deferred to Phase 21).

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Kinematic chain definition (RobotModel) | Library (header-only) | — | Pure data POD structs, zero computation |
| Solver abstraction (RobotSolver) | Library (header-only) | — | Abstract interface, no implementation binding |
| SpatialAlgebra solver binding | Adapter (compiled) | — | sa_adapter.cpp is the sole bridge between test-models and libSpatialAlgebra |
| Model factory functions (chains/) | Library (header-only) | — | Pure Eigen arithmetic to produce RobotModel instances |
| Compilation firewall enforcement | Build System (CMake) | — | INTERFACE target test_models links only Eigen3; sa_test_adapter links SpatialAlgebra |
| Smoke test verification | Test executable | — | compile_smoke_test verifies zero-dependency guarantee compiles |

## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Abstract base class `RobotSolver` with pure virtual methods. No CRTP, no C++20 concepts. Standard polymorphism for GTest-friendly mocking and solver swapping.
- **D-02:** Full robotics solver API: `setState(q, qdot)`, `computeTorques(qddot, gravity)` → tau, `computeAccelerations(tau, gravity)` → qddot, `forwardKinematics()`, `getJointTransform(idx)`, `computeMassMatrix()` → MatrixXd, `computeGravityTorques(gravity)`, `computeJointSpaceJacobian(idx)`, `getLinkCOM(idx)`, `getDOF()`.
- **D-03:** All public methods use dynamic Eigen types (`Eigen::VectorXd`, `Eigen::MatrixXd`). NO templates on scalar or DOF count.
- **D-04:** Hybrid state model: `RobotModel` holds reference joint config; adapter caches internal state. `setState()` copies values and rebuilds internal cache.
- **D-05:** Builder pattern: `RobotModel` is pure data; adapters provide `build(const RobotModel&)` → ready-to-use solver.
- **D-06:** Naming: `RobotSolver` (abstract), `RobotModel` (chain description POD), `JointSpec` (per-link params).
- **D-07:** Pure data POD structs with only Eigen3 types (`Vector3d`, `Matrix3d`, `Matrix4d`). Zero methods beyond accessors, no SpatialAlgebra headers.
- **D-08:** `RobotModel` = array of `JointSpec`. Each `JointSpec`: `parent`, `parentToJoint` (Matrix4d), `jointAxis` (Vector3d), `type` (JointType), `mass`, `com` (Vector3d), `inertia` (Matrix3d dense), `name` (string).
- **D-09:** Joint types: `enum class JointType { REVOLUTE, PRISMATIC, FIXED }`. Revolute: `jointAxis` = rotation axis (ω). Prismatic: `jointAxis` = translation direction (v).
- **D-10:** Library at `tests/test-models/`. Layout: `CMakeLists.txt`, `robot_model.h`, `robot_solver.h`, `sa_adapter.h`, `sa_adapter.cpp`, `chains/` with 11 domain headers.
- **D-11:** One header per test domain in `chains/`. Model-building functions return `RobotModel`. Naming: `makeThreeLinkSerialChain()`, `makeBranchingYConfiguration()`, etc.
- **D-12:** `sa_adapter.h/.cpp` is the ONLY file that includes SpatialAlgebra headers. INTERFACE target does NOT link SpatialAlgebra.
- **D-13:** INTERFACE library target named `test_models`. Links only `Eigen3::Eigen`. Header-only — no compiled source.
- **D-14:** SpatialAlgebra adapter compiled separately as `sa_test_adapter`, linking `SpatialAlgebra`. Tests link `test_models` + optionally `sa_test_adapter`.
- **D-15:** Existing tests are NOT modified in Phase 20. Phase 21 converts existing tests. Phase 20 produces the library + smoke test.

### the agent's Discretion

- Exact `RobotSolver` method signatures (const-correctness, return-by-value vs return-by-reference for large matrices)
- Internal caching strategy for the SA adapter (what gets cached, when to invalidate)
- Whether `forwardKinematics()` is called automatically before `getJointTransform()` or must be explicit
- Error handling convention (exceptions vs error codes for invalid state)
- Whether `RobotModel` stores DOF count explicitly or derives it from `JointSpec` array size

### Deferred Ideas (OUT OF SCOPE)

- Converting existing tests to use the new library — Phase 21
- Pinocchio and RBDL adapters — Phase 22+
- Python bindings for model definitions — out of scope
- URDF/SDF/MJCF parser integration — explicitly excluded

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| TML-01 | Header-only test model library using only Eigen types — zero SpatialAlgebra dependency | INTERFACE CMake target + Eigen-only includes. Verified: INTERFACE targets prevent transitive link pollution. |
| TML-02 | Kinematic chain definitions covering all 11 test domains | Each domain mapped to chains/ header; 2-4 model functions per domain covering the test scenarios found in existing test files. |
| TML-03 | Adapter interface with computeTorques() and computeAccelerations() pure virtual methods | RobotSolver abstract base class mapping to ID+FD solver API. Verified: matches existing ForwardDynamics/InverseDynamics signatures. |
| TML-04 | SpatialAlgebra adapter as thin wrapper around existing solver classes | sa_adapter.cpp converts RobotModel → ForwardDynamics::Link/InverseDynamicsLink arrays. Verified: conversion pattern matches model_factory.cpp. |
| TML-05 | CMake library target TestModels with Eigen3 as sole dependency | Verified: INTERFACE target with `target_link_libraries(test_models INTERFACE Eigen3::Eigen)`. No compiled objects = zero SA dependency. |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Eigen3 | 3.4+ (project requirement) | All matrix/vector types (Vector3d, Matrix3d, Matrix4d, MatrixXd, VectorXd) | Already project dependency; single source of truth for linear algebra [VERIFIED: project CMakeLists.txt] |
| Google Test | 1.12.1 (via FetchContent) | Smoke test verification | Already project dependency; smoke test follows existing test patterns [VERIFIED: project CMakeLists.txt] |
| C++17 STL | — | `std::vector`, `std::string`, `std::unique_ptr`, `std::runtime_error` | Standard library; no additional dependencies needed |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| SpatialAlgebra | (project-local) | SA adapter implementation | ONLY in `sa_adapter.cpp` and tests that link `sa_test_adapter` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Eigen::Matrix4d for transforms | Custom Transform struct | Matrix4d is universally portable to Pinocchio, RBDL, URDF; custom struct would need per-solver conversion |
| Eigen::Matrix3d for dense inertia | LowerTriangular packed storage | Dense 3×3 is portable to all solvers; SA internally uses LowerTriangular but conversion is cheap |
| JointType enum | String-based joint types | Enum is type-safe, compiler-checked; strings are fragile and error-prone |

**Installation:**
```cmake
# In tests/test-models/CMakeLists.txt:
add_library(test_models INTERFACE)
target_include_directories(test_models INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_models INTERFACE Eigen3::Eigen)

add_library(sa_test_adapter sa_adapter.cpp)
target_link_libraries(sa_test_adapter PUBLIC test_models SpatialAlgebra)
```

**Version verification:** Eigen3 version is enforced by root CMakeLists.txt: `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` [VERIFIED: project CMakeLists.txt line 12]. No other external packages are needed for this phase.

## Package Legitimacy Audit

> This phase introduces zero new external packages. The only dependency is Eigen3, which is already a verified project dependency. No `npm install`, `pip install`, or `cargo` operations are required.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| Eigen3 | Homebrew (system) | 10+ yrs | N/A (C++ header lib) | gitlab.com/libeigen/eigen | N/A (not applicable) | Pre-verified project dependency |

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

*No external packages are added in this phase. The library uses only project-internal headers and the pre-existing Eigen3 dependency.*

## Architecture Patterns

### System Architecture Diagram

```
┌──────────────────────────────────┐
│   tests/test-models/             │
│                                  │
│  ┌──────────────────────────┐    │
│  │ robot_model.h            │    │
│  │ JointSpec, RobotModel    │    │
│  │ JointType enum           │    │
│  │ (Eigen3 types ONLY)      │    │
│  └─────────────┬────────────┘    │
│                │                 │
│  ┌─────────────▼────────────┐    │
│  │ robot_solver.h           │    │
│  │ RobotSolver (abstract)   │    │
│  │ (Eigen3 types ONLY)      │    │
│  └─────────────┬────────────┘    │
│                │                 │
│  ┌─────────────▼────────────┐    │
│  │ chains/ (11 headers)     │    │
│  │ make…() → RobotModel     │    │
│  │ (Eigen3 types ONLY)      │    │
│  └──────────────────────────┘    │
│                                  │
│  ┌──────────────────────────┐    │
│  │ sa_adapter.h + .cpp      │◄───┼── SpatialAlgebra headers
│  │ (ONLY file with SA deps) │    │
│  └──────────────────────────┘    │
│                                  │
└──────────────────────────────────┘
         │                    │
         │ test_models        │ sa_test_adapter
         │ (INTERFACE)        │ (STATIC library)
         │ links: Eigen3 only │ links: SpatialAlgebra
         │                    │
         ▼                    ▼
    ┌──────────────────────────────┐
    │  CMake build targets          │
    │                              │
    │  compile_smoke_test          │
    │    links: test_models        │
    │    (verifies zero SA dep)    │
    │                              │
    │  Future test executables     │
    │    links: test_models +      │
    │           sa_test_adapter    │
    └──────────────────────────────┘
```

**Data flow for test model usage (Phase 21+):**

```
chains/makeThreeLinkSerial()
  │
  ▼
RobotModel { joints: [JointSpec×3] }
  │
  ├──► SAAdapter::build(model) ──► SAAdapter (cached links)
  │      │
  │      ├── setState(q, qdot)
  │      ├── computeTorques(qddot, gravity) ──► InverseDynamics::computeTorques()
  │      ├── computeAccelerations(tau, gravity) ──► ForwardDynamics::computeAccelerations()
  │      └── computeMassMatrix() ──► internal assembly from links
  │
  └──► (Future) PinocchioAdapter::build(model) ──► pinocchio::Model + pinocchio::Data
```

### Recommended Project Structure

```
tests/test-models/
├── CMakeLists.txt                # INTERFACE target test_models + sa_test_adapter
├── robot_model.h                 # JointSpec, RobotModel, JointType enum
├── robot_solver.h                # RobotSolver abstract base class
├── sa_adapter.h                  # SpatialAlgebraAdapter declaration
├── sa_adapter.cpp                # SpatialAlgebraAdapter implementation
└── chains/                       # Model definitions, one header per test domain
    ├── spatial_vectors.h         # Single-link chain, motion vector fixtures
    ├── plucker_transforms.h      # Transform fixtures (identity, 90° Z, translation)
    ├── rotation.h                # Rotation angle-axis fixtures
    ├── lower_triangular.h        # LT matrix fixtures (identity, diagonal)
    ├── rigid_body_inertia.h      # RBI configurations (mass, COM, inertia variants)
    ├── articulated_body.h        # ABI configurations (I, H, M variants)
    ├── spatial_utils.h           # Vector pairs for cross/dot product tests
    ├── inverse_dynamics.h        # Serial chains (1, 2, 3 link), branching Y
    ├── forward_dynamics.h        # Same chains + rotated, gravity variants
    ├── consistency.h             # Round-trip chain pairs (ID+FD equivalent)
    └── spatial_operations.h      # Transform+RBI pairs for spatial operations
```

### Pattern 1: INTERFACE CMake Target for Header-Only Library

**What:** An INTERFACE library target does not compile any source files. It only specifies include directories and linked dependencies. Including its headers cannot transitively pull in anything not in its `INTERFACE` link libraries.

**When to use:** When the zero-dependency guarantee must be enforced by the build system, not just by convention.

**Example:**
```cmake
# tests/test-models/CMakeLists.txt
add_library(test_models INTERFACE)
target_include_directories(test_models INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_models INTERFACE Eigen3::Eigen)
# NOTE: SpatialAlgebra is NOT linked here — enforcing zero-dependency guarantee.

# Compiled adapter (links SA)
add_library(sa_test_adapter sa_adapter.cpp)
target_link_libraries(sa_test_adapter PUBLIC test_models SpatialAlgebra)
```
[CITED: tests/test-models/ directory structure per D-10; CMake INTERFACE target pattern per D-13]

### Pattern 2: Model-to-Solver Conversion (Adapter)

**What:** The SA adapter converts pure-Eigen `RobotModel` → SA-specific link arrays. This follows the exact pattern from `benchmarks/common/model_factory.cpp`.

**Example:**
```cpp
// Source: benchmarks/common/model_factory.cpp lines 19-50 (createFD)
// The SA adapter follows the same conversion pattern:

ForwardDynamics SAAdapter::buildFD(const RobotModel& model) {
    ForwardDynamics fd;
    fd.links.reserve(model.joints.size());
    for (const auto& js : model.joints) {
        Link link;
        link.parent = js.parent;
        
        // Convert Matrix4d → PluckerTransform
        link.X = PluckerTransform(
            Rotation(js.parentToJoint.topLeftCorner<3,3>()),
            js.parentToJoint.topRightCorner<3,1>()
        );
        
        // Convert mass/com/dense-inertia → RigidBodyInertia
        link.I = RigidBodyInertia(
            js.mass, js.com, 
            LowerTriangular::fromFullMatrix(js.inertia)
        );
        
        // Convert Vector3d + JointType → MotionVector (screw axis)
        switch (js.type) {
            case JointType::REVOLUTE:
                link.S = MotionVector(js.jointAxis, Vector3d::Zero());
                break;
            case JointType::PRISMATIC:
                link.S = MotionVector(Vector3d::Zero(), js.jointAxis);
                break;
            case JointType::FIXED:
                link.S = MotionVector(Vector3d::Zero(), Vector3d::Zero());
                break;
        }
        fd.links.push_back(std::move(link));
    }
    return fd;
}
```
[VERIFIED: model_factory.cpp conversion pattern; CONTEXT.md D-08 JointSpec design]

### Pattern 3: Model Factory Functions (chains/)

**What:** Each `chains/<domain>.h` header declares 2-4 functions that return `RobotModel` instances with specific kinematic configurations. These are pure data factories — no solver logic, no SA types.

**Example (inverse_dynamics.h):**
```cpp
#pragma once
#include "robot_model.h"

namespace test_models {

/**
 * @brief Single-link Z-revolute chain for basic ID tests
 * @details 1-DOF chain: mass=1, COM=origin, identity inertia,
 *          Z-axis revolute joint.
 * @return RobotModel with one JointSpec (parent=-1, parentToJoint=I)
 */
inline RobotModel makeSingleLinkChain() {
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "base";
    model.joints.push_back(link);
    return model;
}

// ... more functions

} // namespace test_models
```
[CITED: CONTEXT.md D-11, "one header per test domain in chains/"]

### Anti-Patterns to Avoid

- **Accidental SA inclusion:** Including any header from `include/` in `robot_model.h` or `robot_solver.h`. The compilation firewall depends on these headers never transitively pulling in SA types. Solution: use forward declarations or Eigen-only types in the interface.
- **Template-heavy adapter:** Templating the SA adapter on DOF count or scalar type. D-03 explicitly forbids this — use dynamic Eigen types (`VectorXd`, `MatrixXd`) everywhere.
- **Computation in RobotModel:** Adding methods like `computeMass()` or `getTotalInertia()` to `RobotModel`. D-07 defines it as pure data POD. All computation belongs in adapters.
- **Mixing adapter concerns:** Putting Pinocchio-specific or RBDL-specific logic in `sa_adapter.cpp`. Each adapter is a separate compilation unit — SA adapter handles only SA types.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Matrix/vector types | Custom linear algebra structs | `Eigen::Vector3d`, `Eigen::Matrix3d`, `Eigen::Matrix4d`, `Eigen::VectorXd`, `Eigen::MatrixXd` | Eigen is the project standard; portable to all robotics libraries |
| Kinematic transform math | Custom 4×4 transform | `Eigen::Matrix4d` (homogeneous transform) | Standard robotics convention; trivial conversion to PluckerTransform, SE3, URDF transforms |
| Rotation representation | Custom rotation class | `Eigen::Matrix4d::topLeftCorner<3,3>()` for rotation part | Extracted from homogeneous matrix on demand; avoids duplicating rotation storage |
| Build-system dependency firewall | Manual header discipline | CMake INTERFACE target | Build system enforces the boundary; accidental includes become compile errors, not runtime bugs |
| Virtual dispatch | Custom dispatch/vtable | C++ virtual methods (abstract base class) | D-01 explicitly requires standard polymorphism; GTest mock-compatible |

**Key insight:** The zero-dependency guarantee is enforced by the build system, not by convention. If anyone accidentally adds `#include "PluckerTransform.h"` to `robot_model.h`, the `test_models` INTERFACE target would fail to compile because `PluckerTransform.h` depends on `SpatialVector.h` which is not in the include path chain. This makes the boundary self-policing.

## Runtime State Inventory

> **SKIPPED:** Phase 20 is a greenfield library creation phase. It introduces new files into `tests/test-models/` but does not rename, refactor, or migrate any existing code, data, or state. No runtime state is affected.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — greenfield library, no data migration | None |
| Live service config | None — no services affected | None |
| OS-registered state | None — no OS registrations | None |
| Secrets/env vars | None — no secrets or env vars modified | None |
| Build artifacts | None — new library target, no existing artifacts to migrate | None |

## Common Pitfalls

### Pitfall 1: Accidental Transitive Include Chain

**What goes wrong:** A `chains/` header includes `robot_model.h` (ok), but `robot_model.h` includes `<Eigen/Dense>` which pulls in Eigen's internal headers. Someone later adds `#include "SpatialVector.h"` to `robot_model.h` thinking Eigen is already available, breaking the zero-dependency guarantee.

**Why it happens:** Developers assume "Eigen is already included" means "SA headers are safe to add." The include path isn't the issue — it's the semantic dependency.

**How to avoid:** The compile smoke test (`compile_smoke_test.cpp`) specifically tests that `robot_model.h` and `robot_solver.h` compile without any SA headers. If any SA type appears in these headers, the smoke test fails. Check this in CI.

**Warning signs:** Smoke test fails with "unknown type 'PluckerTransform'" or similar. Indicates SA type leaked into test_models headers.

### Pitfall 2: LowerTriangular Conversion Loss

**What goes wrong:** `JointSpec::inertia` stores a dense 3×3 `Eigen::Matrix3d`. Converting to `LowerTriangular` via `fromFullMatrix()` extracts only the lower triangle. If the input matrix is not symmetric (e.g., due to numerical noise), information is silently lost.

**Why it happens:** Dense matrix → packed storage conversion discards the upper triangle. For a truly symmetric inertia tensor this is correct, but if the matrix has accumulated floating-point asymmetry from upstream conversions, the result may be subtly wrong.

**How to avoid:** In the SA adapter, symmetrize the dense inertia before conversion: `0.5 * (M + M.transpose())`. This is standard practice and cheap (3×3 matrix). Document the symmetrization in a comment.

**Warning signs:** Non-zero `(inertia - inertia.transpose()).norm()` in debug output. Asymmetric inertia tensors produce physically incorrect dynamics.

### Pitfall 3: PluckerTransform vs Homogeneous Matrix Convention Mismatch

**What goes wrong:** Featherstone's PluckerTransform uses `[R, 0; -R*skew(r), R]` for motion transforms, where `r` is the translation FROM parent TO child expressed IN the parent frame. The `JointSpec::parentToJoint` is a standard homogeneous matrix `[R | t; 0 | 1]` where `t` is the translation in the parent frame. These match exactly for the rotation part, but the sign convention for translation must be verified.

**Why it happens:** Different robotics libraries use different conventions for the translation direction in spatial transforms.

**How to avoid:** Follow the exact conversion pattern used in `model_factory.cpp`: `PluckerTransform(Rotation(matrix.topLeftCorner<3,3>()), matrix.topRightCorner<3,1>())`. The `model_factory.cpp` at lines 28-31 shows the canonical conversion pattern.

**Warning signs:** Sign errors in dynamics results when comparing SA adapter output with direct SA usage. If `computeTorques()` via adapter produces different signs from direct `InverseDynamics` calls with equivalent parameters, the transform conversion is suspect.

### Pitfall 4: INTERFACE Target and add_subdirectory Order

**What goes wrong:** `test-models/CMakeLists.txt` is added via `add_subdirectory()` from the root `CMakeLists.txt`. If `add_subdirectory(tests/test-models)` is placed before `find_package(Eigen3)`, the INTERFACE target cannot link `Eigen3::Eigen`.

**Why it happens:** CMake processes subdirectories in order. The INTERFACE target needs `Eigen3::Eigen` to be available as an imported target.

**How to avoid:** Place `add_subdirectory(tests/test-models)` AFTER the existing `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` call at line 12 of the root CMakeLists.txt. Current file order is already correct — add the subdirectory call near the existing test registration block (around line 58, after `enable_testing()`).

**Warning signs:** CMake error: "Target 'test_models' links to 'Eigen3::Eigen' but the target was not found."

## Test Domain Extraction: Model Inventory

Research into the 11 existing test files reveals the following model requirements per domain:

### Domain 1: Spatial Vectors (spatial_vectors.h)
**Source:** `tests/TestSpatialVector.cpp` (609 lines, GTest)  
**Model patterns:** Tests use `MotionVector`/`ForceVector` constructors directly with `Vector3d` components. No kinematic chains. Tests verify arithmetic (add, subtract, scale), dot products, cross products, and property invariants (anti-commutativity, distributivity).  
**Chain models needed:**
- Single-link Z-revolute chain for transform tests (e.g., Featherstone Example 2.1 uses spatial vectors with transforms)
- `makeSingleLink()` — 1-DOF, mass=1, Z-axis revolute, COM=origin, identity inertia

### Domain 2: Plücker Transforms (plucker_transforms.h)
**Source:** `tests/TestPluckerTransform.cpp` (940 lines, GTest)  
**Model patterns:** Tests construct `PluckerTransform` directly from `Rotation` + `Vector3d`. Test fixtures: identity, 90° Z rotation, pure translation [1,0,0], combined rotation+translation, non-orthogonal axes.  
**Chain models needed:**
- `makeIdentityTransform()` — identity rotation, zero translation
- `make90DegreeZRotation()` — 90° Z, zero translation
- `makePureTranslation()` — identity rotation, [1,0,0] translation
- `makeCombinedTransform()` — 90° Z + [1,0,0] translation

### Domain 3: Rotation (rotation.h)
**Source:** `tests/TestRotation.cpp` (414 lines, GTest)  
**Model patterns:** Tests use `Rotation` directly from `AngleAxisd`, `Quaterniond`, `Matrix3d`.  
**Chain models needed:**
- `makeIdentityRotation()` — zero rotation
- `make45DegreeXRotation()` — 45° around X
- `make90DegreeZRotationFixture()` — 90° around Z

### Domain 4: Lower Triangular (lower_triangular.h)
**Source:** `tests/TestLowerTriangular.cpp` (451 lines, GTest)  
**Model patterns:** Tests use `LowerTriangular` directly with packed data or `fromFullMatrix`.  
**Chain models needed:**
- `makeIdentityLT(size=3)` — identity matrix
- `makeDiagonalLT(size=3, value=2.0)` — scaled identity 
- `makeArbitraryLT()` — non-trivial lower triangular values

### Domain 5: Rigid Body Inertia (rigid_body_inertia.h)
**Source:** `tests/TestRigidBodyInertia.cpp` (376 lines, GTest)  
**Model patterns:** Tests construct `RigidBodyInertia(mass, com, lt)` directly. Tested: default zero, parameterized (mass=5, COM=[1,2,3]), apply() with pure rotation/translation/combined.  
**Chain models needed:**
- `makeRBIDefault()` — mass=1, COM=origin, identity LT
- `makeRBIOffsetCOM(mass, comX, comY, comZ)` — configurable COM
- `makeRBIDiagonalInertia(mass, com, diag)` — configurable inertia diagonal

### Domain 6: Articulated Body Inertia (articulated_body.h)
**Source:** `tests/TestArticulatedBodyInertia.cpp` (531 lines, GTest)  
**Model patterns:** Tests construct `ArticulatedBodyInertia(I, H, M)` with 3×3 matrices and `LowerTriangular`. Tested: default, parameterized, addition, scaling, apply(), reduced case matching RBI.  
**Chain models needed:**
- `makeABIIdentity()` — I=M=identity, H=identity
- `makeABIReduced()` — H=zero, M=identity (matches RBI behavior)

### Domain 7: Spatial Utils (spatial_utils.h)
**Source:** `tests/TestSpatialUtils.cpp` (335 lines, GTest)  
**Model patterns:** Tests use `skew()`, `dot()`, `cross()` free functions. Test fixtures: pairs of `MotionVector` and `ForceVector` for cross/dot tests.  
**Chain models needed:**
- `makeSkewFixture(v)` — converts Vector3d to skew-symmetric matrix
- `makeMotionVectorPair()` — two motion vectors for cross/dot tests
- `makeForceVectorPair()` — two force vectors for cross/dot tests

### Domain 8: Inverse Dynamics (inverse_dynamics.h)
**Source:** `tests/TestInverseDynamics.cpp` (550 lines, GTest)  
**Model patterns:** Builds `InverseDynamicsLink` arrays manually. Tested chains:
- Single-link Z-revolute (mass=1, COM=origin, identity inertia)
- Two-link serial Z-revolute (1m X spacing)
- Branching Y (base + two children at ±1m X)
- Non-zero velocity (qdot=2.0, 1.0) with COM offset [0,0.1,0]
- Gravity variants (single-link with -Z gravity, 2-link combined)
- Static gravity proportionality (X-axis joint, COM=[0,0.5,0], g={0,5,10})
- Zero-torque-at-vertical (collinear COM and gravity)
- Release-mode stability  
**Chain models needed:**
- `makeSingleLinkChain()` — 1-DOF, Z-revolute, mass=1, COM=zero
- `makeTwoLinkSerialChain()` — 2-DOF, Z-revolute, 1m spacing, mass=1
- `makeBranchingYConfiguration()` — 3-DOF, base + two symmetric children
- `makeTwoLinkWithCOMOffset()` — 2-DOF with COM=[0,0.1,0] on each link
- `makeSingleLinkXAxisWithCOMOffset()` — X-axis revolute, COM=[0,0.5,0]

### Domain 9: Forward Dynamics (forward_dynamics.h)
**Source:** `tests/TestForwardDynamics.cpp` (738 lines, GTest)  
**Model patterns:** Same chains as ID but with `ForwardDynamics::Link`. Additional tests:
- PluckerTransformUsage (90° Z rotation on link 1)
- ThreeLinkNumericalValidation (3-link with COM=[0.1,0,0])
- CondensationReducesInertiaNorm (property test on Ia)
- ThreeLinkSingleTorque (tau=[1,0,0])
- SingleLinkWithGravity, TwoLinkWithGravity
- GravityEffectScalesWithMass (m=1 → m=2, ratio=2.0)
- GravityProportionalityInvariant (diff_g1/diff_g2 = g1/g2)
- ZeroMassEdgeCase (degenerate inertia)
- ReleaseModeStability  
**Chain models needed:** All ID chains plus:
- `makeThreeLinkSerialChain()` — 3-DOF Z-revolute, 1m spacing, mass=1
- `makeThreeLinkSerialChainNonZeroCOM()` — 3-DOF with COM=[0.1,0,0]
- `makeTwoLinkWith90DegreeRotation()` — link 1 has 90° Z rotation

### Domain 10: Consistency (consistency.h)
**Source:** `tests/TestDynamicsConsistency.cpp` (526 lines, GTest)  
**Model patterns:** Duplicate chains for ID and FD — both constructed from equivalent parameters. Tests:
- RoundTripABARNEA (1-link)
- RoundTripRNEAABA (1-link)
- ThreeLinkSerialChain (3-link consistency)
- ThreeLinkSerialChainNonZeroCOM (3-link with COM offset, CR-02 test)
- BranchingYConfiguration (3-link Y tree)
- TwoLinkRoundTripWithGravity
- TwoLinkGravityNonZeroCOM (CR-02 test)
- RoundTripABA_RNEA_DirectComparison (X-axis joint, single-link)  
**Chain models needed:** Same chains as ID+FD — the consistency models are parameter-equivalent copies of the ID and FD chains. The `RobotModel` is the single source of truth; both adapters construct from the same model.

### Domain 11: Spatial Operations (spatial_operations.h)
**Source:** `tests/TestSpatialOperations.cpp` (424 lines, GTest)  
**Model patterns:** Tests `SpatialOperations::crossProductMotion`, `crossProductForce`, `transformInertia`. Uses `MotionVector`, `ForceVector`, `RigidBodyInertia`, `PluckerTransform` directly.  
**Chain models needed:**
- `makeTransformRBIPair()` — PluckerTransform + RigidBodyInertia for transformInertia tests
- `makeCrossProductMotionFixtures()` — MotionVector pairs for crossMotion tests
- `makeCrossProductForceFixtures()` — ForceVector/MotionVector pairs for crossForce tests

## RobotModel Struct Design Validation

### Matrix4d for parentToJoint — CONFIRMED

`Eigen::Matrix4d` is the correct type for homogeneous transforms. This is the universal robotics convention and is directly convertible to:
- PluckerTransform: `Rotation(M.topLeftCorner<3,3>())` + `M.topRightCorner<3,1>()` [VERIFIED: model_factory.cpp lines 28-31]
- Pinocchio SE3: `pinocchio::SE3(M.topLeftCorner<3,3>(), M.topRightCorner<3,1>())` [ASSUMED standard Pinocchio API]
- RBDL Math::SpatialTransform: similar pattern [ASSUMED]
- URDF origin: `origin.xyz = t`, `origin.rpy = euler angles from R` [ASSUMED]

The bottom row `[0 0 0 1]` convention is assumed (verified by construction in factory functions).

### Matrix3d Dense Inertia — CONFIRMED

The 3×3 dense `Eigen::Matrix3d` for rotational inertia at COM is adequate and portable:
- SpatialAlgebra: `LowerTriangular::fromFullMatrix(inertia)` trivial conversion [VERIFIED: LowerTriangular API]
- Pinocchio: `pinocchio::Inertia(mass, com, inertia)` takes dense 3×3 [ASSUMED]
- The 6 unique elements of a symmetric 3×3 are fully captured in a 3×3 dense matrix.

### JointType Enum — CONFIRMED

Three types (`REVOLUTE`, `PRISMATIC`, `FIXED`) match the existing joint configurations in the test suite:
- All 11 test domains use REVOLUTE joints (Z-axis or X-axis)
- PRISMATIC is needed for Phase 21 (TST-01: "prismatic joints, mixed revolute+prismatic")
- FIXED is used for static transforms (non-actuated links)

Conversion to SpatialAlgebra `MotionVector`:
- REVOLUTE: `MotionVector(axis, Vector3d::Zero())` — screw axis with pure angular component
- PRISMATIC: `MotionVector(Vector3d::Zero(), axis)` — screw axis with pure linear component
- FIXED: `MotionVector(Vector3d::Zero(), Vector3d::Zero())` — zero motion subspace
[VERIFIED: model_factory.cpp uses MotionVector for all joint axes; Featherstone convention]

### RobotModel DOF Count — DISCRETION

Whether to store `dof` explicitly or derive from `joints.size()`:
- **Option A:** `int dof` field (explicit). Pro: self-documenting, allows non-1:1 joint-to-DOF mapping. Con: risk of inconsistency with `joints.size()`.
- **Option B:** `int getDOF() const { return static_cast<int>(joints.size()); }` (derived). Pro: cannot be inconsistent. Con: assumes exactly 1 DOF per joint.
- **Recommendation:** Use derived `getDOF()`. All existing test models use 1-DOF-per-joint exclusively, and Phase 21 multi-DOF joints can extend the model if needed. The planner should use Option B for simplicity, adding an explicit `dof` field only if Phase 21 requirements demand it. [ASSUMED based on test file analysis; all 11 domains use 1-DOF-per-joint]

## CMake Integration Strategy

### Placement of add_subdirectory

The `add_subdirectory(tests/test-models)` should be inserted in the root `CMakeLists.txt`:

```cmake
# Root CMakeLists.txt (after line 55: enable_testing())
enable_testing()

# NEW: Test model library (zero SA dependency)
add_subdirectory(tests/test-models)

# Existing test executables (unchanged for Phase 20)
add_executable(TestSpatialVector tests/TestSpatialVector.cpp)
# ...
```

**Why here:** `find_package(Eigen3)` is already resolved (line 12). `enable_testing()` is already called. The `tests/test-models/CMakeLists.txt` can add smoke tests via `add_test()`. Existing test targets are not modified — they don't link `test_models` or `sa_test_adapter` in Phase 20.

### INTERFACE Target Mechanics

```cmake
# tests/test-models/CMakeLists.txt
add_library(test_models INTERFACE)
target_include_directories(test_models INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_models INTERFACE Eigen3::Eigen)
```

The INTERFACE target ensures:
1. No compiled objects — headers only, can't accidentally link SA `.o` files
2. `target_include_directories` adds `tests/test-models/` to the include path, so `#include "robot_model.h"` resolves
3. `target_link_libraries(... Eigen3::Eigen)` propagates Eigen include dirs transitively

**Critical:** `test_models` does NOT link `SpatialAlgebra`. If any test-models header transitively includes an SA header, the compiler will fail because SA headers are not on the include path for this target. This is a build-system-enforced guarantee.

### sa_test_adapter Compilation

```cmake
add_library(sa_test_adapter sa_adapter.cpp)
target_link_libraries(sa_test_adapter PUBLIC test_models SpatialAlgebra)
target_include_directories(sa_test_adapter PRIVATE ${CMAKE_SOURCE_DIR}/include)
```

The adapter links BOTH `test_models` (for `RobotModel`, `JointSpec`, `RobotSolver` types) and `SpatialAlgebra` (for `ForwardDynamics`, `InverseDynamics`, `PluckerTransform`, etc.). This is the ONLY compilation unit that touches both worlds.

**Installation (TML-05):** The `test_models` INTERFACE target can be made installable via:
```cmake
install(TARGETS test_models EXPORT SpatialAlgebraTargets)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ 
    DESTINATION include/test-models
    FILES_MATCHING PATTERN "*.h"
)
```
This satisfies TML-05's requirement that `TestModels` be "installable as a standalone library."

### Smoke Test Integration

```cmake
# Add to tests/test-models/CMakeLists.txt
add_executable(compile_smoke_test 
    ${CMAKE_SOURCE_DIR}/tests/compile_smoke_test.cpp
)
target_link_libraries(compile_smoke_test PRIVATE test_models)
# NOTE: Does NOT link SpatialAlgebra — this is the zero-dependency verification
add_test(NAME compile_smoke_test COMMAND compile_smoke_test)
```

The smoke test should include `robot_model.h` and `robot_solver.h` and instantiate types to verify compilation. It does NOT include any SA headers. [CITED: CONTEXT.md line 147, "compile smoke test that #include 'robot_model.h' and #include 'robot_solver.h' compile WITHOUT including any SpatialAlgebra headers"]

## Code Examples

### RobotModel and JointSpec Definition

```cpp
// Source: CONTEXT.md D-07, D-08, D-09
// File: tests/test-models/robot_model.h

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace test_models {

/**
 * @brief Joint type enumeration
 * @details REVOLUTE: rotational joint about jointAxis
 *          PRISMATIC: translational joint along jointAxis
 *          FIXED: rigid connection (0-DOF)
 */
enum class JointType { REVOLUTE, PRISMATIC, FIXED };

/**
 * @brief Per-link specification for a kinematic chain
 * @details Pure data POD struct using only Eigen3 types.
 *          No SpatialAlgebra dependency. Designed to be
 *          compatible with Featherstone's spatial vector
 *          algebra (Featherstone 2008).
 * 
 *          The parent-to-joint transform uses a 4×4 homogeneous
 *          matrix [R | t; 0 | 1] in the parent frame.
 *          The inertia matrix is a 3×3 dense rotational inertia
 *          at the center of mass.
 */
struct JointSpec {
    int parent = -1;                     ///< Parent link index (-1 for base)
    Eigen::Matrix4d parentToJoint =      ///< Homogeneous transform parent→link
        Eigen::Matrix4d::Identity();
    Eigen::Vector3d jointAxis =          ///< Joint axis direction
        Eigen::Vector3d::UnitZ();
    JointType type = JointType::REVOLUTE;///< Joint type
    double mass = 1.0;                   ///< Link mass (kg)
    Eigen::Vector3d com =                ///< Center of mass in link frame
        Eigen::Vector3d::Zero();
    Eigen::Matrix3d inertia =            ///< 3×3 rotational inertia at COM
        Eigen::Matrix3d::Identity();
    std::string name;                    ///< Optional link name
};

/**
 * @brief Complete kinematic chain description
 * @details An ordered collection of JointSpec entries.
 *          Links must be in topological order (parents before children).
 *          For branching trees, children of the same parent must have
 *          consecutive indices.
 * 
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */
struct RobotModel {
    std::vector<JointSpec> joints;       ///< Link specifications
    
    /**
     * @brief Number of degrees of freedom
     * @return Number of joints (1 DOF per joint)
     */
    int getDOF() const {
        return static_cast<int>(joints.size());
    }
};

} // namespace test_models
```

### RobotSolver Abstract Interface

```cpp
// Source: CONTEXT.md D-01, D-02, D-03
// File: tests/test-models/robot_solver.h

#pragma once
#include <Eigen/Dense>
#include <memory>
#include "robot_model.h"

namespace test_models {

/**
 * @brief Abstract solver interface for robotics dynamics
 * @details Pure virtual interface enabling solver-agnostic
 *          test models. Implementations wrap specific solvers
 *          (SpatialAlgebra, Pinocchio, RBDL) behind a uniform API.
 * 
 *          State management: setState() loads joint configuration;
 *          computation methods read from cached internal state.
 * 
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */
class RobotSolver {
public:
    virtual ~RobotSolver() = default;

    /**
     * @brief Set joint configuration and velocity
     * @param q Joint positions (size = DOF)
     * @param qdot Joint velocities (size = DOF)
     * @details Copies state and rebuilds internal cache as needed.
     */
    virtual void setState(const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& qdot) = 0;

    /**
     * @brief Compute joint torques via inverse dynamics (RNEA)
     * @param qddot Joint accelerations (size = DOF)
     * @param gravity Gravity vector in world frame (default zero)
     * @return Joint torques (size = DOF)
     */
    virtual Eigen::VectorXd computeTorques(
        const Eigen::VectorXd& qddot,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) = 0;

    /**
     * @brief Compute joint accelerations via forward dynamics (ABA)
     * @param tau Joint torques (size = DOF)
     * @param gravity Gravity vector in world frame (default zero)
     * @return Joint accelerations (size = DOF)
     */
    virtual Eigen::VectorXd computeAccelerations(
        const Eigen::VectorXd& tau,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) = 0;

    /**
     * @brief Update internal joint transforms
     * @details Recomputes link-to-world transforms from current
     *          joint positions. Must be called before getJointTransform().
     */
    virtual void forwardKinematics() = 0;

    /**
     * @brief Get world-frame transform for a link
     * @param idx Link index (0 to DOF-1)
     * @return 4×4 homogeneous transform
     */
    virtual Eigen::Matrix4d getJointTransform(int idx) const = 0;

    /**
     * @brief Compute joint-space inertia matrix H(q)
     * @return DOF×DOF symmetric matrix
     */
    virtual Eigen::MatrixXd computeMassMatrix() = 0;

    /**
     * @brief Compute gravity compensation torques
     * @param gravity Gravity vector in world frame
     * @return Static joint torques (size = DOF)
     */
    virtual Eigen::VectorXd computeGravityTorques(
        const Eigen::Vector3d& gravity) = 0;

    /**
     * @brief Compute Jacobian for a specific link
     * @param idx Link index
     * @return 6×DOF Jacobian matrix (spatial velocity)
     */
    virtual Eigen::MatrixXd computeJointSpaceJacobian(int idx) = 0;

    /**
     * @brief Get link COM in world frame
     * @param idx Link index
     * @return World-frame COM position
     */
    virtual Eigen::Vector3d getLinkCOM(int idx) const = 0;

    /**
     * @brief Number of degrees of freedom
     * @return Total joint DOF
     */
    virtual int getDOF() const = 0;
};

} // namespace test_models
```

### SpatialAlgebra Adapter Declaration

```cpp
// Source: CONTEXT.md D-04, D-05, D-12
// File: tests/test-models/sa_adapter.h

#pragma once
#include "robot_solver.h"
#include "robot_model.h"

namespace test_models {

/**
 * @brief SpatialAlgebra adapter implementing RobotSolver
 * @details Wraps ForwardDynamics and InverseDynamics classes
 *          behind the solver-agnostic RobotSolver interface.
 *          Constructed via build() factory method from a
 *          RobotModel description.
 * 
 *          Internal state: maintains ForwardDynamics and
 *          InverseDynamics solver objects as cached internal
 *          representations. setState() updates joint positions
 *          and velocities on both. Computation methods delegate
 *          to the appropriate solver.
 * 
 *          This is the ONLY file in tests/test-models/ that
 *          includes SpatialAlgebra headers (per D-12).
 */
class SpatialAlgebraAdapter : public RobotSolver {
public:
    /**
     * @brief Build an SA adapter from a RobotModel
     * @param model Kinematic chain description
     * @return Fully initialized adapter
     */
    static std::unique_ptr<SpatialAlgebraAdapter> build(
        const RobotModel& model);

    // RobotSolver interface implementation
    void setState(const Eigen::VectorXd& q,
                  const Eigen::VectorXd& qdot) override;
    Eigen::VectorXd computeTorques(
        const Eigen::VectorXd& qddot,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;
    Eigen::VectorXd computeAccelerations(
        const Eigen::VectorXd& tau,
        const Eigen::Vector3d& gravity = Eigen::Vector3d::Zero()) override;
    void forwardKinematics() override;
    Eigen::Matrix4d getJointTransform(int idx) const override;
    Eigen::MatrixXd computeMassMatrix() override;
    Eigen::VectorXd computeGravityTorques(
        const Eigen::Vector3d& gravity) override;
    Eigen::MatrixXd computeJointSpaceJacobian(int idx) override;
    Eigen::Vector3d getLinkCOM(int idx) const override;
    int getDOF() const override;

private:
    SpatialAlgebraAdapter() = default;  // Private; use build()
    
    // Internal SA solver objects (implementation detail)
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace test_models
```

**Note on Impl pattern:** Uses PIMPL (`struct Impl`) to hide SpatialAlgebra types (`ForwardDynamics`, `InverseDynamics`) from the header. This is critical: `sa_adapter.h` must NOT expose SA types in its public interface, or it would defeat the zero-dependency guarantee for consumers that include it. The Impl struct is defined in `sa_adapter.cpp` where SA headers are included. [CITED: CONTEXT.md D-12, "ONLY file that includes SpatialAlgebra headers"]

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Hardcoded model construction in each test | Declarative RobotModel + factory functions | Phase 20 | Single source of truth; enables cross-solver comparison |
| Direct SA dependency in test definitions | Zero-dependency Eigen-only models | Phase 20 | Test models reusable across SA, Pinocchio, RBDL |
| Manual ID+FD duplication for consistency tests | Single RobotModel → two adapters | Phase 20 | Eliminates parameter drift between ID and FD chains |
| JointConfig (benchmarks) with SA MotionVector | JointSpec with Eigen::Vector3d + JointType | Phase 20 | Solver-agnostic, portable across libraries |

**Deprecated/outdated:**
- `benchmarks/common/model_factory.h:JointConfig` — replaced by `test_models::JointSpec`. The `JointConfig` uses SA types (`MotionVector`) and is solver-specific. The new `JointSpec` is solver-agnostic. (ModelFactory itself remains useful for benchmarks but is not the canonical model source.)

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Pinocchio SE3 constructor accepts `SE3(R.topLeftCorner<3,3>(), R.topRightCorner<3,1>())` | RobotModel Struct Design | Phase 22 Pinocchio adapter conversion formula changes; minor rework |
| A2 | All 11 test domains' model requirements are adequately captured by 2-4 factory functions each | Test Domain Extraction | Some edge cases in existing tests won't be representable; requires adding factory functions in Phase 21 |
| A3 | Derived `getDOF()` from `joints.size()` is sufficient (1 DOF per joint) | RobotModel Struct Design | Multi-DOF joints in Phase 21 would require model extension; low effort to add explicit `dof` field |
| A4 | `Eigen::Matrix4d` bottom row `[0 0 0 1]` convention is standard and doesn't need explicit verification | RobotModel Struct Design | Non-homogeneous transforms would produce incorrect Plücker conversions; caught by smoke test |
| A5 | PIMPL pattern in `sa_adapter.h` successfully hides all SA types from the public header | SA Adapter Design | SA types would leak into `test_models` consumers, violating zero-dependency guarantee; caught by compile smoke test |

**If this table is empty:** All claims in this research were verified or cited — no user confirmation needed.

## Open Questions

1. **`forwardKinematics()` auto-call vs explicit call**
   - What we know: D-02 specifies `forwardKinematics()` as part of the API. Pinocchio recomputes transforms on every `computeTorques`/`computeAccelerations` call. SA currently does NOT recompute transforms (known limitation).
   - What's unclear: Should `getJointTransform()`, `computeMassMatrix()`, `computeGravityTorques()` automatically call `forwardKinematics()` if state is dirty? Or should the user call `forwardKinematics()` explicitly before these queries?
   - Recommendation: Auto-call `forwardKinematics()` internally with a dirty flag. `setState()` marks state as dirty; first call to any query method triggers `forwardKinematics()`. This matches Pinocchio's behavior and prevents stale transform bugs. The planner should implement a `bool _fkDirty` flag pattern.

2. **`computeMassMatrix()` implementation approach**
   - What we know: SA has no built-in mass matrix computation. The composite rigid body algorithm (CRBA) or unit-vector method could compute H(q). RNEA-based column extraction uses N+1 calls to `computeTorques()`.
   - What's unclear: Whether to implement CRBA in the adapter (complex, but O(n) vs O(n²)) or use the simpler RNEA-column method.
   - Recommendation: Use the RNEA-column method for Phase 20 (simple, correct, easy to verify). Call `computeTorques(e_i, zero_gravity)` for each unit vector e_i to extract column i of H(q). This is O(n²) but for n ≤ 10 (test model purpose), it's negligible. The planner can optimize later if needed.

3. **`computeJointSpaceJacobian()` implementation for SA**
   - What we know: SA has no built-in Jacobian computation. The Jacobian J_i(q) maps joint velocities to end-effector spatial velocity: v_i = J_i * qdot.
   - What's unclear: Whether to implement geometric Jacobian (column-by-column using Plücker transforms) or analytical Jacobian. SA lacks automatic differentiation.
   - Recommendation: Implement geometric Jacobian column-by-column. Column j is the world-frame joint screw axis for link idx: J_col(j) = X_world * S_j if j affects link i, else zero. This is standard Featherstone (2008, §3.6). The planner should implement this directly in the adapter using the cached link transforms.

4. **Error handling for invalid state**
   - What we know: SA classes throw `std::invalid_argument` for size mismatches and `std::runtime_error` for singular configurations.
   - What's unclear: Whether the adapter should propagate SA exceptions, wrap them, or use error codes.
   - Recommendation: Propagate SA exceptions directly. The adapter is a thin wrapper — it doesn't add new error semantics. Document that `computeTorques()` and `computeAccelerations()` may throw `std::invalid_argument` (size mismatch) or `std::runtime_error` (singular config). This is consistent with the existing codebase conventions [VERIFIED: ForwardDynamics.h throws std::invalid_argument, InverseDynamics.h throws std::invalid_argument].

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| CMake | Build system | ✓ | 3.19+ (project requirement) | — |
| C++17 compiler (g++/clang++) | Compilation | ✓ | System default | — |
| Eigen3 | test_models (INTERFACE) | ✓ | 3.4+ via Homebrew | — |
| Google Test | Smoke test compilation | ✓ | 1.12.1 via FetchContent | — |
| SpatialAlgebra (lib) | sa_test_adapter | ✓ | Project-local build | — |

**Missing dependencies with no fallback:** none
**Missing dependencies with fallback:** none

*All required dependencies are available on the development machine. The phase is self-contained — no new tools or libraries are needed.*

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Google Test 1.12.1 (GTest) |
| Config file | none — GTest is auto-detected via CMake `find_package(GTest)` |
| Quick run command | `cmake --build build && cd build && ./compile_smoke_test` |
| Full suite command | `cmake --build build && cd build && ctest --output-on-failure` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| TML-01 | `robot_model.h` + `robot_solver.h` compile without SA headers | compile | `cmake --build build --target compile_smoke_test` | ❌ Wave 0 |
| TML-02 | All 11 chains/ headers produce valid RobotModel instances | unit | `cd build && ./compile_smoke_test` (verifies includes, basic instantiation) | ❌ Wave 0 |
| TML-03 | RobotSolver pure virtual class compiles and is mockable | compile | `cmake --build build --target compile_smoke_test` | ❌ Wave 0 |
| TML-04 | SA adapter wraps ForwardDynamics/InverseDynamics correctly | integration | `cd build && ctest -R smoke` (Phase 20); full round-trip in Phase 21 | ❌ Wave 0 |
| TML-05 | test_models INTERFACE target links only Eigen3 | build | `cmake --build build --target test_models` (no-op for INTERFACE) | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `cmake --build build && cd build && ctest -R smoke`
- **Per wave merge:** `cmake --build build && cd build && ctest --output-on-failure`
- **Phase gate:** All 11 existing test executables pass (no regression) + smoke test passes

### Wave 0 Gaps
- [ ] `tests/compile_smoke_test.cpp` — update to include `robot_model.h` and `robot_solver.h`, instantiate `RobotModel` and `JointSpec` types (currently tests SA umbrella header)
- [ ] `tests/test-models/CMakeLists.txt` — new file defining `test_models` INTERFACE target, `sa_test_adapter` library, and `compile_smoke_test` executable
- [ ] `tests/test-models/robot_model.h` — new file: JointSpec, RobotModel, JointType enum
- [ ] `tests/test-models/robot_solver.h` — new file: RobotSolver abstract base class
- [ ] `tests/test-models/sa_adapter.h` — new file: SpatialAlgebraAdapter declaration
- [ ] `tests/test-models/sa_adapter.cpp` — new file: SA adapter implementation
- [ ] `tests/test-models/chains/*.h` — 11 new files: model factory functions per test domain
- [ ] Root `CMakeLists.txt` — add `add_subdirectory(tests/test-models)` after `enable_testing()`

*(Major gaps: the entire `tests/test-models/` directory is new. All files are Wave 0.)*

## Security Domain

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|------------------|
| V2 Authentication | No | N/A — library has no authentication surface |
| V3 Session Management | No | N/A — no sessions |
| V4 Access Control | No | N/A — no access controls |
| V5 Input Validation | Yes | Parameter validation: size checks on `q`, `qdot`, `qddot`, `tau` vectors vs `getDOF()`; bounds checks on `idx` parameters; NaN/Inf detection on inputs |
| V6 Cryptography | No | N/A — no cryptographic operations |

### Known Threat Patterns for C++ Header-Only Library

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Size mismatch: `tau.size() != getDOF()` | Denial of Service | `std::invalid_argument` exception (propagated from underlying SA solver) |
| Index out-of-bounds: `getJointTransform(n)` where `n >= getDOF()` | Tampering | `std::out_of_range` exception (matching SA conventions) |
| NaN/Inf propagation through Eigen operations | Denial of Service | Default Eigen behavior propagates NaN/Inf to output; test assertions catch these at verification boundaries |
| Memory exhaustion from large RobotModel | Denial of Service | `std::vector` raises `std::bad_alloc` on OOM; caller controls model size |
| Use-after-move of adapter after `build()` | Tampering | `std::unique_ptr` ownership prevents accidental double-use; use-after-move is undefined behavior but rare in test usage |

**Mitigation note:** The SA adapter delegates all computation to the existing, tested `ForwardDynamics` and `InverseDynamics` classes. Security concerns are primarily input validation (ensuring size consistency) and numerical stability (ensuring finite inputs produce finite outputs). The SA classes already handle singular configuration detection (throws `std::runtime_error`). The adapter adds only size-checking validation.

## Sources

### Primary (HIGH confidence)
- `benchmarks/common/model_factory.h` — JointConfig struct design, ModelFactory API (conversion pattern reference) [CITED]
- `benchmarks/common/model_factory.cpp` — SA link construction pattern (PluckerTransform from rotation+translation, RigidBodyInertia from mass+com+lt) [VERIFIED]
- `include/ForwardDynamics.h` — ForwardDynamicsLink struct (parent, X, I, S, q, qdot, qddot, v, c, f, Ia, pa fields) [VERIFIED]
- `include/InverseDynamics.h` — InverseDynamicsLink struct (parent, X, I, S, q, qdot, qddot, v, a fields) [VERIFIED]
- `include/RigidBodyInertia.h` — Constructor: `RigidBodyInertia(mass, com, lt)` [VERIFIED]
- `CMakeLists.txt` (root) — Eigen3 find_package, GTest FetchContent, test executable registration [VERIFIED]
- `.planning/phases/20-test-model-library/20-CONTEXT.md` — Locked decisions D-01 through D-15 [CITED]
- `.planning/PROJECT.md` — Eigen3 compatibility requirement, C++17 standard [VERIFIED]

### Secondary (MEDIUM confidence)
- `tests/TestForwardDynamics.cpp` — Model construction patterns (single-link, 2-link, 3-link, branching, PluckerTransform, gravity, edge cases) [VERIFIED]
- `tests/TestInverseDynamics.cpp` — RNEA model construction (singles, serial, branching, velocity, gravity, edge cases) [VERIFIED]
- `tests/TestDynamicsConsistency.cpp` — Round-trip model patterns (ID+FD duplication, non-zero COM, gravity) [VERIFIED]
- `tests/TestPluckerTransform.cpp` — Transform fixture patterns (identity, 90° Z, translation, combined) [VERIFIED]
- `tests/TestSpatialOperations.cpp` — crossProductMotion, crossProductForce, transformInertia fixture patterns [VERIFIED]
- `tests/TestSpatialVector.cpp` — MotionVector/ForceVector test patterns [VERIFIED]
- `tests/TestRigidBodyInertia.cpp` — RBI construction and apply() patterns [VERIFIED]
- `tests/TestArticulatedBodyInertia.cpp` — ABI construction and apply() patterns [VERIFIED]
- `tests/TestLowerTriangular.cpp` — LT packed storage and matrix operations [VERIFIED]
- `tests/TestRotation.cpp` — Rotation constructors and operations [VERIFIED]
- `tests/TestSpatialUtils.cpp` — skew(), dot(), cross() usage patterns [VERIFIED]
- `tests/compile_smoke_test.cpp` — Existing smoke test pattern (self-contained include test) [VERIFIED]
- `examples/dynamics.cpp` — Manual link construction example [VERIFIED]

### Tertiary (LOW confidence)
- Pinocchio API assumptions (SE3 constructor, Model+Data structures) — marked [ASSUMED], verification needed in Phase 22
- RBDL API assumptions — not used in this phase, out of scope

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — Eigen3 is the only external dependency, already verified in the project build
- Architecture: HIGH — CMake INTERFACE target pattern is well-understood; conversion pattern is verified from model_factory.cpp
- Pitfalls: MEDIUM — Pitfalls identified from code review of existing conversion patterns; runtime edge cases (branching tree handling, NaN propagation) are partially validated by existing tests
- Test domain mapping: HIGH — All 11 test files were read and analyzed; model requirements directly extracted from test code
- SA adapter design: HIGH — Conversion formulas verified against model_factory.cpp patterns; PIMPL boundary verified as standard C++ pattern

**Research date:** 2026-06-17
**Valid until:** 2026-07-17 (30 days — stable domain, no fast-moving ecosystem dependencies)
