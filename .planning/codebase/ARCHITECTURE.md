<!-- refreshed: 2026-06-05 -->
# Architecture

**Analysis Date:** 2026-06-05

## System Overview

```text
┌─────────────────────────────────────────────────────────────────────┐
│                        ALGORITHM LAYER                              │
│  ┌──────────────────────┐  ┌──────────────────────────────────────┐ │
│  │   ForwardDynamics    │  │        InverseDynamics                │ │
│  │   (Articulated Body  │  │   (Recursive Newton-Euler)           │ │
│  │    Algorithm)        │  │                                      │ │
│  │  `src/ForwardDynamics` │  │  `src/InverseDynamics.cpp`          │ │
│  └──────────┬───────────┘  └──────────────┬───────────────────────┘ │
│             │                              │                         │
│  ┌──────────▼──────────────────────────────▼───────────────────────┐ │
│  │                  SpatialOperations (static utility)             │ │
│  │                  `include/SpatialOperations.h`                  │ │
│  └─────────────────────────────────────────────────────────────────┘ │
├─────────────────────────────────────────────────────────────────────┤
│                       GEOMETRY / INERTIA LAYER                       │
│  ┌────────────────┐  ┌──────────────────┐  ┌──────────────────────┐ │
│  │  PluckerTransform│  │  RigidBodyInertia│  │ArticulatedBodyInertia│ │
│  │  `include/...`  │  │  `include/...`   │  │  `include/...`      │ │
│  └───────┬────────┘  └────────┬─────────┘  └──────────┬───────────┘ │
│          │                    │                         │             │
│  ┌───────▼────────────────────▼─────────────────────────▼───────────┐ │
│  │               LowerTriangular (packed storage)                    │ │
│  │               `include/LowerTriangular.h`                         │ │
│  └───────────────────────────────────────────────────────────────────┘ │
├─────────────────────────────────────────────────────────────────────┤
│                       SPATIAL VECTOR LAYER                           │
│  ┌────────────────┐  ┌──────────────────────┐  ┌──────────────────┐ │
│  │   MotionVector │  │     SpatialUtils     │  │    ForceVector   │ │
│  │    (twist)     │  │   cross(), dot(),     │  │    (wrench)      │ │
│  │                │  │   skew()             │  │                   │ │
│  └───────┬────────┘  └──────────────────────┘  └─────────┬─────────┘ │
│          │                                                │           │
│  ┌───────▼────────────────────────────────────────────────▼─────────┐ │
│  │                    SpatialVector (6D base class)                   │ │
│  │                    `include/SpatialVector.h`                       │ │
│  └───────────────────────────────────────────────────────────────────┘ │
├─────────────────────────────────────────────────────────────────────┤
│                        FOUNDATION LAYER                              │
│  ┌────────────────────┐  ┌─────────────────────────────────────────┐ │
│  │  Rotation          │  │  Eigen3 (all linear algebra)            │ │
│  │  extends Matrix3d  │  │  Used by every class above              │ │
│  │  `include/Rotation.h│  │  External dependency                   │ │
│  └────────────────────┘  └─────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| SpatialVector | 6D base: angular + linear components, arithmetic, cross/dot products | `include/SpatialVector.h`, `src/SpatialVector.cpp` |
| MotionVector | Type-safe twist (motion) representation | `include/MotionVector.h`, `src/MotionVector.cpp` |
| ForceVector | Type-safe wrench (force) representation | `include/ForceVector.h`, `src/ForceVector.cpp` |
| Rotation | 3D rotation matrix extending Eigen::Matrix3d | `include/Rotation.h`, `src/Rotation.cpp` |
| PluckerTransform | 6x6 Plücker coordinate transform (rotation+translation) | `include/PluckerTransform.h`, `src/PluckerTransform.cpp` |
| RigidBodyInertia | Body mass, COM, inertia tensor (fully inline) | `include/RigidBodyInertia.h` |
| ArticulatedBodyInertia | Articulated body [I, H; H^T, M] block inertia | `include/ArticulatedBodyInertia.h` |
| LowerTriangular | Packed-storage lower triangular matrix (custom, not Eigen) | `include/LowerTriangular.h`, `src/LowerTriangular.cpp` |
| SpatialUtils | Free functions: skew(), dot(), cross() (forced `noexcept`) | `include/SpatialUtils.h` |
| SpatialOperations | Static utility: crossProductMotion, crossProductForce, transformInertia | `include/SpatialOperations.h`, `src/SpatialOperations.cpp` |
| InverseDynamics | Recursive Newton-Euler Algorithm solver | `include/InverseDynamics.h`, `src/InverseDynamics.cpp` |
| ForwardDynamics | Articulated Body Algorithm solver | `include/ForwardDynamics.h`, `src/ForwardDynamics.cpp` |
| SpatialAlgebra.h | Umbrella header (includes all public headers in dependency order) | `include/SpatialAlgebra.h` |

## Pattern Overview

**Overall:** Layered architecture with inheritance for spatial vector types, composition for transforms, and protocol-based link structures for dynamics algorithms.

**Key Characteristics:**
- `SpatialVector` → `MotionVector` / `ForceVector`: inheritance used for type safety (not polymorphism; no virtual methods except destructor)
- `Rotation` extends `Eigen::Matrix3d`: inheritance used for seamless Eigen integration
- `PluckerTransform` uses composition: stores `Rotation` + `Vector3d` translation separately (not a 6x6 matrix)
- `RigidBodyInertia` and `ArticulatedBodyInertia` use `LowerTriangular` for packed symmetric storage (3x3 inertia tensor)
- `LowerTriangular` is a custom packed-storage class (not Eigen-based internally, but uses Eigen vectors for data)
- Dynamics algorithms (`ForwardDynamics`, `InverseDynamics`) are struct-of-arrays-like, storing vectors of `Link`/`InverseDynamicsLink`
- All types use value semantics (pass by value/const reference, no heap allocation)
- All code lives in `namespace SpatialAlgebra`
- Mathematical notation follows Featherstone (2008) throughout

## Layers

**Foundation Layer:**
- Purpose: Fundamental 3D math and linear algebra primitives
- Location: `include/Rotation.h`, `include/SpatialVector.h`, `include/SpatialUtils.h`, `include/LowerTriangular.h`
- Contains: `Rotation` (3x3 matrix extending Eigen), `SpatialVector` (6D base), `SpatialUtils` (free function cross/dot/skew), `LowerTriangular` (packed matrix)
- Depends on: Eigen3 (`Eigen/Dense`, `Eigen/Geometry`)
- Used by: All higher layers

**Spatial Vector Layer:**
- Purpose: Type-safe motion (twist) and force (wrench) 6D vectors
- Location: `include/MotionVector.h`, `include/ForceVector.h`, `src/MotionVector.cpp`, `src/ForceVector.cpp`
- Contains: `MotionVector` (angular→ω, linear→v), `ForceVector` (angular→τ, linear→f), each with typed arithmetic and cross/dot operations
- Depends on: `SpatialVector`
- Used by: `PluckerTransform`, inertia classes, dynamics algorithms

**Geometry/Inertia Layer:**
- Purpose: Coordinate frame transforms, mass property representations
- Location: `include/PluckerTransform.h`, `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`, `src/PluckerTransform.cpp`
- Contains: `PluckerTransform` (motion/force transforms + inertia transformation), `RigidBodyInertia` (mass, COM, inertia), `ArticulatedBodyInertia` (block inertia for ABA)
- Depends on: `SpatialVector`, `MotionVector`, `ForceVector`, `Rotation`, `LowerTriangular`
- Used by: Dynamics algorithms

**Algorithm Layer:**
- Purpose: Rigid body dynamics solvers
- Location: `include/ForwardDynamics.h`, `include/InverseDynamics.h`, `src/ForwardDynamics.cpp`, `src/InverseDynamics.cpp`, `include/SpatialOperations.h`, `src/SpatialOperations.cpp`
- Contains: `ForwardDynamics` (ABA for joint accelerations), `InverseDynamics` (RNEA for joint torques), `SpatialOperations` (static wrappers)
- Depends on: All lower layers
- Used by: `src/main.cpp` (demo), `examples/` (usage demonstrations)

## Data Flow

### Primary Request Path — Inverse Dynamics (RNEA)

1. User populates `InverseDynamics::links` with `InverseDynamicsLink` structs (parent, X, I, S, q, qdot, qddot) (`include/InverseDynamics.h:75-100`)
2. User calls `computeTorques(qddot, gravity)` (`src/InverseDynamics.cpp:96-136`)
3. **Outward pass** (base→tip): Propagates velocities and accelerations through kinematic chain (`src/InverseDynamics.cpp:19-52`)
   - Base: `v = S·q̇`, `a = S·q̈ - g`
   - Children: `v = X·v_parent + S·q̇`, `a = X·a_parent + S·q̈ + v×S·q̇`
4. **Inward pass** (tip→base): Propagates forces, computes joint torques (`src/InverseDynamics.cpp:55-94`)
   - `f = I·a + v×I·v`
   - `τ = S·f`
   - Propagate to parent: `f_parent += X⁻ᵀ·f`
5. Returns `Eigen::VectorXd tau` of joint torques

### Primary Request Path — Forward Dynamics (ABA)

1. User populates `ForwardDynamics::links` with `Link` structs (parent, X, I, S, q, qdot) (`include/ForwardDynamics.h:79-111`)
2. User calls `computeAccelerations(tau, gravity)` (`src/ForwardDynamics.cpp:182-209`)
3. **Outward pass** (base→tip): Propagates velocities, computes bias accelerations (`src/ForwardDynamics.cpp:57-80`)
   - Base: `v = S·q̇`, `c = -g`
   - Children: `v = X·v_parent + S·q̇`, `c = X·c_parent + v×S·q̇`
4. **Inward pass** (tip→base): Accumulates articulated inertias, computes partial accelerations (`src/ForwardDynamics.cpp:82-180`, 3 phases)
   - Phase 1: Initialize `Ia`, `pa` from rigid body inertia + bias forces
   - Phase 2: Condense `Ia`/`pa` along joint axis, propagate to parent
   - Phase 3: Correct partial `qddot` using parent acceleration
5. Results stored in `links[i].qddot`

### Transform Pipeline (PluckerTransform)

1. Motion vector transform: `v' = X·v = [R·ω; R·(v - r×ω)]` (`src/PluckerTransform.cpp:15-26`)
2. Force vector transform: `f' = X⁻ᵀ·f = [R·(τ - r×f); R·f]` (`src/PluckerTransform.cpp:28-40`)
3. Inverse motion: `v = X⁻¹·v' = [Rᵀ·ω'; Rᵀ·v' + r×(Rᵀ·ω')]` (`src/PluckerTransform.cpp:42-51`)
4. Inverse force: `f = Xᵀ·f' = [Rᵀ·τ' + r×(Rᵀ·f'); Rᵀ·f']` (`src/PluckerTransform.cpp:53-64`)
5. Composition: `X₁·X₂` uses `R = R₁·R₂`, `t = t₂ + R₂ᵀ·t₁` (`src/PluckerTransform.cpp:229-240`)
6. Inverse: `X⁻¹` uses `R_inv = Rᵀ`, `t_inv = -R·t` (`src/PluckerTransform.cpp:217-227`)
7. Inertia transform (`tformRBI`): `I' = X·I·Xᵀ` (`src/PluckerTransform.cpp:66-83`)
8. ABI transform (`tformABI`): uses 6x6 block matrix construction (`src/PluckerTransform.cpp:106-163`)

**State Management:**
- Pure value semantics: no shared pointers, no global state
- Dynamics solvers store link vectors as mutable state; user configures before each call
- No persistent state between dynamics solves (caller must re-populate or mutate links)
- No lazy evaluation, no caching

## Key Abstractions

**SpatialVector (6D base):**
- Purpose: 6D vector combining angular + linear 3D components
- Examples: `include/SpatialVector.h:67-174`
- Pattern: Base class with protected `Vector3d angular, linear` members; virtual destructor for type safety
- Provides: arithmetic (+/-/*), `crossMotion`, `crossForce`, `dot`, `print`

**MotionVector / ForceVector (typed spatial vectors):**
- Purpose: Type-safe wrappers with correct physical interpretation
- Examples: `include/MotionVector.h:69-146`, `include/ForceVector.h:71-148`
- Pattern: Inherit `SpatialVector`, add typed operators returning derived type
- Type aliases: `using mv = MotionVector`, `using fv = ForceVector`

**Rotation (3D rotation):**
- Purpose: 3x3 rotation matrix extending Eigen
- Examples: `include/Rotation.h:57-171`
- Pattern: Inherits `Eigen::Matrix3d`, adds constructors from AngleAxis/Quaternion, `inverse()`, `transpose()`, `toAngleAxis()`, `toQuaternion()`

**PluckerTransform (6x6 spatial transform):**
- Purpose: Rigid body coordinate transformation in Plücker coordinates
- Examples: `include/PluckerTransform.h:77-216`
- Pattern: Stores `Rotation` and `Vector3d translation` separately; provides `transformMotion`, `transformForce`, inverse variants, and inertia transformations (`tformRBI`, `tformABI`)
- Type alias: `using plux = PluckerTransform`

**RigidBodyInertia:**
- Purpose: Body mass properties (mass, COM, inertia tensor)
- Examples: `include/RigidBodyInertia.h:27-116`
- Pattern: Stores `double mass`, `Vector3d com`, `LowerTriangular inertiaMatrixLT`; fully inline implementation
- Key method: `apply(mv)` → returns `fv = [Iω + com×v; m·v - com×ω]`
- Type alias: `using rbi = RigidBodyInertia`

**ArticulatedBodyInertia:**
- Purpose: Block inertia for articulated bodies: [I, H; Hᵀ, M]
- Examples: `include/ArticulatedBodyInertia.h:77-213`
- Pattern: Stores `lt Inertia` (rotational), `Matrix3d H` (coupling), `lt M` (mass); inline implementation
- Key method: `apply(mv)` → returns `fv = [Iω + H·v; Hᵀω + M·v]`
- Type alias: `using abi = ArticulatedBodyInertia`

**LowerTriangular (packed storage):**
- Purpose: Memory-efficient lower triangular matrix for symmetric storage
- Examples: `include/LowerTriangular.h:75-564`
- Pattern: 1D `Eigen::VectorXd` array with index mapping `idx = i*(i+1)/2 + j`; provides multiply, inverse, transpose, `getSymmetricMatrix()`, `fromFullMatrix()`, `multiplySymmetric()`
- OpenMP parallelization commented in docs, not observed in implementation
- Type alias: `using lt = LowerTriangular`

**Link / InverseDynamicsLink:**
- Purpose: Protocol struct for kinematic chain configuration
- Examples: `include/ForwardDynamics.h:79-111`, `include/InverseDynamics.h:75-100`
- Pattern: Struct with parent index, X (PluckerTransform), I (RigidBodyInertia), S (joint axis MotionVector), q/qdot/qddot, and intermediate computed quantities (v, c, a, f, Ia, pa)

## Entry Points

**Public API Headers:**
- Location: `include/*.h` (13 headers)
- Triggers: User includes one or more headers and instantiates classes
- Responsibilities: Provide full spatial algebra API
- Umbrella header: `include/SpatialAlgebra.h` includes all public headers in order

**Dynamics Solvers:**
- `InverseDynamics::computeTorques(qddot, gravity)` — entry for inverse dynamics (`src/InverseDynamics.cpp:96`)
- `ForwardDynamics::computeAccelerations(tau, gravity)` — entry for forward dynamics (`src/ForwardDynamics.cpp:182`)

**Demo Executable:**
- `src/main.cpp`: Quick demonstration of vectors, transforms, inertia
- `examples/basic_vectors.cpp`, `examples/transforms.cpp`, `examples/inertia.cpp`, `examples/dynamics.cpp`: Usage examples built as separate executables

**Test Executables:**
- 10 test executables registered in `CMakeLists.txt` (lines 52-154): `TestSpatialVector`, `TestPluckerTransform`, `TestRotation`, `TestLowerTriangular`, `TestSpatialUtils`, `TestRigidBodyInertia`, `TestArticulatedBodyInertia`, `TestForwardDynamics`, `TestSpatialOperations`, `TestInverseDynamics`
- Plus: `TestDynamicsConsistency` (cross-check between inverse and forward dynamics)

**Python Side:**
- `robot_dynamics/rnea.py`: Standalone RNEA implementation using NumPy (not integrated with C++ library)

## Architectural Constraints

- **Threading:** Single-threaded. OpenMP is mentioned in `LowerTriangular.h` documentation but not used in the actual implementation. No thread safety considerations.
- **Global state:** None. No module-level singletons or shared mutable state. Each dynamics solver has its own `links` vector.
- **Circular imports:** None detected. Header dependency graph is a DAG with `SpatialVector.h` at root. `SpatialAlgebra.h` includes in strict dependency order.
- **Memory model:** All objects are value types. `LowerTriangular` owns a heap-allocated `Eigen::VectorXd` internally (RAII). Dynamics solvers own `std::vector<Link>` by value.
- **No runtime polymorphism:** The virtual destructor on `SpatialVector` is the only virtual method. All operations use static dispatch through typed derived classes.
- **`noexcept` correctness:** Only `SpatialUtils.h` free functions (`skew`, `dot`, `cross` variants) are marked `noexcept`. Class methods are not.

## Anti-Patterns

### Inconsistency: `noexcept` on free functions but not on class methods

**What happens:** `SpatialUtils.h` marks all free functions `noexcept`, but identical logic in `SpatialVector` (e.g., `dot()`, `crossMotion()`) and derived classes is not marked `noexcept`.
**Why it's wrong:** Inconsistent exception specification makes it unclear which operations can throw. Most vector operations are pure math and cannot throw.
**Do this instead:** Mirror the `noexcept` annotation from `SpatialUtils.h` onto equivalent `SpatialVector`/`MotionVector`/`ForceVector` methods, as seen in `include/SpatialUtils.h:24`.

### SpatialOperations is a thin wrapper

**What happens:** `SpatialOperations` static methods (`crossProductMotion`, `crossProductForce`, `transformInertia`) simply delegate to free functions or `PluckerTransform` methods (`src/SpatialOperations.cpp:11-24`).
**Why it's wrong:** The class adds no value — it duplicates existing API surface without additional logic or abstraction.
**Do this instead:** Eliminate `SpatialOperations` and call `cross()` / `tformRBI()` directly.

### Mixed test style

**What happens:** `TestSpatialVector.cpp` originally used bare `assert()`, while `TestPluckerTransform.cpp` uses GTest (`TEST()` / `EXPECT_DOUBLE_EQ`). The test file now has been partially converted but bare asserts may remain in older revisions.
**Why it's wrong:** Inconsistent test patterns make it harder to run all tests uniformly (bare `assert()` aborts on failure instead of reporting).
**Do this instead:** Use GTest exclusively across all test files.

## Error Handling

**Strategy:** Exceptions for runtime errors, debug-only assertions for preconditions.

**Patterns:**
- `std::invalid_argument` for dimension mismatches (`LowerTriangular` operations, `computeAccelerations`, `computeTorques`)
- `std::out_of_range` for index bounds (debug mode only in `LowerTriangular`)
- `std::runtime_error` for singular matrix / near-zero inertia in ABA (`ForwardDynamics::inwardPass`)
- No error codes, no `std::optional`, no custom exception types
- NaN/Inf detection in debug mode with `std::cerr` warnings (`RigidBodyInertia::apply`, `ArticulatedBodyInertia::apply`, `SpatialVector` constructor)
- Input validation: dynamics solvers check for NaN/Inf in input vectors before computation

## Cross-Cutting Concerns

**Logging:** All classes have `print()` method outputting to `std::cout`. No structured logging, no log levels.

**Documentation:** Every header uses Doxygen `@brief`/`@details`/`@param`/`@return`/`@note`/`@warning`/`@see` on every declaration. Config: `Doxyfile` (Doxygen 1.12.0).

**Validation:** Debug-mode NaN/Inf checks in constructors and critical methods, guarded by `#ifndef NDEBUG`.

**Testing:** 10 GTest-based test executables, 1 `TestDynamicsConsistency` cross-check, 1 `compile_smoke_test.cpp` for compilation verification. Results registered with CTest.

---

*Architecture analysis: 2026-06-05*
