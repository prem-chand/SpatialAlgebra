<!-- refreshed: 2026-05-17 -->
# Architecture

**Analysis Date:** 2026-05-17

## System Overview

```text
┌─────────────────────────────────────────────────────────────────┐
│                     APPLICATION LAYER                           │
│  `examples/`  `src/main.cpp`  `robot_dynamics/rnea.py`          │
├─────────────────────────────────────────────────────────────────┤
│                     DYNAMICS ALGORITHMS                         │
│  ┌──────────────────────┐  ┌──────────────────────────────────┐ │
│  │  ForwardDynamics     │  │  InverseDynamics                │ │
│  │  (ABA, Algorithm7.3) │  │  (RNEA, Algorithm 7.1)          │ │
│  │  `src/FD.cpp`         │  │  `src/ID.cpp`                  │ │
│  │  `include/FD.h`      │  │  `include/ID.h`                │ │
│  ├──────────────────────┤  ├──────────────────────────────────┤ │
│  │ struct Link (ABA)    │  │ struct InverseDynamicsLink (RNEA)│ │
│  │  parent, X, I, S,    │  │  parent, X, I, S,               │ │
│  │  v, c, f, Ia, pa     │  │  v, a                           │ │
│  └──────────────────────┘  └──────────────────────────────────┘ │
├─────────────────────────────────────────────────────────────────┤
│                     TRANSFORMS & INERTIA                         │
│  ┌──────────────────┐  ┌──────────────────┐  ┌───────────────┐  │
│  │  PluckerTransform│  │  RigidBodyInertia│  │ ArticBodyIner │  │
│  │  6×6 motion/force│  │  m, com, I_lt    │  │ I, H, M (ABI) │  │
│  │  transform       │  │  apply(mv)->fv   │  │ apply(mv)->fv │  │
│  └────────┬─────────┘  └──────────────────┘  └───────┬───────┘  │
│           │                                          │          │
├───────────┼──────────────────────────────────────────┼──────────┤
│           ▼                                          ▼          │
│                    CORE SPATIAL VECTORS                          │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │  Rotation ──┬── extends Eigen::Matrix3d                 │   │
│  │  SpatialVector (base)         ┌──────────────────┐      │   │
│  │   ├── motion/force duality    │ LowerTriangular  │      │   │
│  │   ├── 6D [angular; linear]    │ packed storage   │      │   │
│  │   ├── crossMotion / crossForce│ O(n²/2) memory   │      │   │
│  │   └── dot()                   └──────────────────┘      │   │
│  │    ↳ MotionVector (twist) ──── cross/ +/dot             │   │
│  │    ↳ ForceVector (wrench) ──── cross/ +/dot             │   │
│  └─────────────────────────────────────────────────────────┘   │
├─────────────────────────────────────────────────────────────────┤
│                       LINEAR ALGEBRA BACKEND                    │
│                         Eigen3 3.3+                             │
│                    `#include <Eigen/Dense>`                     │
└─────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| SpatialVector | 6D vector base: angular+linear components, cross/dot ops | `include/SpatialVector.h`, `src/SpatialVector.cpp` |
| MotionVector | Twist type: type-safe motion vectors, covariant transforms | `include/MotionVector.h`, `src/MotionVector.cpp` |
| ForceVector | Wrench type: type-safe force vectors, contravariant transforms | `include/ForceVector.h`, `src/ForceVector.cpp` |
| Rotation | 3D rotation: angle-axis/quaternion/3×3 matrix, extends Eigen::Matrix3d | `include/Rotation.h`, `src/Rotation.cpp` |
| PluckerTransform | 6×6 spatial transform: motion/force/inertia across coordinate frames | `include/PluckerTransform.h`, `src/PluckerTransform.cpp` |
| RigidBodyInertia | Body mass properties: mass, COM, inertia tensor (LT storage) | `include/RigidBodyInertia.h` (all inline) |
| ArticulatedBodyInertia | Composite inertia for ABA: I, H, M block structure (LT storage) | `include/ArticulatedBodyInertia.h` (all inline) |
| LowerTriangular | Packed lower-triangular matrix: O(n²/2) storage, symmetric mult | `include/LowerTriangular.h`, `src/LowerTriangular.cpp` |
| SpatialUtils | Free functions: skew(), dot(), cross() overloads | `include/SpatialUtils.h` |
| SpatialOperations | Static utility class: wraps cross ops, transformInertia | `include/SpatialOperations.h`, `src/SpatialOperations.cpp` |
| ForwardDynamics | ABA solver: outward→inward pass, resolves q̈ from τ | `include/ForwardDynamics.h`, `src/ForwardDynamics.cpp` |
| InverseDynamics | RNEA solver: outward→inward pass, resolves τ from q̈ | `include/InverseDynamics.h`, `src/InverseDynamics.cpp` |
| `struct Link` | Forward dynamics link data: parent, X, I, S, v, c, f, Ia, pa | `include/ForwardDynamics.h:79-111` |
| `struct InverseDynamicsLink` | Inverse dynamics link data: parent, X, I, S, v, a | `include/InverseDynamics.h:75-100` |

## Pattern Overview

**Overall:** Object-oriented spatial algebra library following Featherstone's formulation. Core pattern is **base-class specialization** for motion/force duality with **composition** for transforms (PluckerTransform contains Rotation + translation bound together in a 6×6 transformation matrix).

**Key Characteristics:**
- **Class inheritance vs composition hybrid**: `Rotation` extends `Eigen::Matrix3d` (inheritance). `SpatialVector` is the base class, `MotionVector` and `ForceVector` inherit from it. `PluckerTransform` composes `Rotation` + `Vector3d translation`.
- **Value semantics**: All objects passed by value or const reference. No pointer-based ownership. Small objects returned by value.
- **Packed storage specialization**: `LowerTriangular` uses a custom 1D packed array instead of Eigen's dense storage for inertia tensors, with O(n²/2) memory footprint.
- **Three-phase algorithm pattern**: Both dynamics solvers (ABA, RNEA) follow outward-pass → inward-pass → solve loop structure.
- **Inline-heavy design**: `RigidBodyInertia` and `ArticulatedBodyInertia` are entirely inline in their headers. Other classes split declarations in headers and implementations in `.cpp` files.

## Layers

**Core Layer (Spatial Vectors & Utils):**
- Purpose: Fundamental 6D vector types and linear algebra utilities
- Location: `include/SpatialVector.h`, `include/MotionVector.h`, `include/ForceVector.h`, `include/SpatialUtils.h`
- Contains: SpatialVector (base), MotionVector (twist), ForceVector (wrench), skew/dot/cross free functions, LowerTriangular (packed storage)
- Depends on: Eigen3 (`<Eigen/Dense>`)
- Used by: All higher-level layers

**Transforms & Inertia Layer:**
- Purpose: Coordinate transformations and mass property representations
- Location: `include/Rotation.h`, `include/PluckerTransform.h`, `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`, `include/LowerTriangular.h`
- Contains: Rotation (3D rotations extending Eigen::Matrix3d), PluckerTransform (6×6 Plücker transforms), RigidBodyInertia (mass+COM+inertia_tensor), ArticulatedBodyInertia (I/H/M block-matrix form)
- Depends on: Core Layer (SpatialVector, MotionVector, ForceVector, SpatialUtils)
- Used by: Dynamics Algorithms Layer, Application Layer

**Dynamics Algorithms Layer:**
- Purpose: O(n) recursive dynamics solvers for kinematic chains
- Location: `include/ForwardDynamics.h`, `include/InverseDynamics.h`, `src/ForwardDynamics.cpp`, `src/InverseDynamics.cpp`
- Contains: ForwardDynamics class (ABA algorithm), InverseDynamics class (RNEA algorithm), Link and InverseDynamicsLink structs
- Depends on: Transforms & Inertia Layer (PluckerTransform, RigidBodyInertia, ArticulatedBodyInertia, SpatialVector family)
- Used by: Application Layer

**Application Layer:**
- Purpose: Library demos, testing, and standalone Python implementation
- Location: `src/main.cpp`, `examples/`, `tests/`, `robot_dynamics/rnea.py`
- Contains: Demonstration executables (basic_vectors, transforms, inertia, dynamics), all test suites, standalone Python RNEA
- Depends on: All library layers

## Data Flow

### Primary Request Path — Forward Dynamics (ABA)

1. **Setup**: User populates `ForwardDynamics::links` vector with `Link` structs (parent, X, I, S, q, qdot) (`include/ForwardDynamics.h:137`)
2. **Entry**: `fd.computeAccelerations(tau)` called with joint torque vector (`src/ForwardDynamics.cpp:128`)
3. **Outward Pass (base→tip)**: Propagate spatial velocities from parent to child links, compute bias accelerations (`src/ForwardDynamics.cpp:19-51`):
   - Base: `v₀ = S₀·q̇₀`, `c₀ = 0`
   - Child: `vᵢ = Xᵢ·v_parent + Sᵢ·q̇ᵢ`, `cᵢ = Xᵢ·c_parent + vᵢ × Sᵢ·q̇ᵢ`
4. **Inward Pass (tip→base)**: Initialize articulated inertias Iₐ from rigid body inertias, accumulate child contributions (`src/ForwardDynamics.cpp:53-126`):
   - Initialize: `Iₐᵢ = rbi_to_abi(Iᵢ)`, `pₐᵢ = Iₐᵢ·cᵢ + vᵢ × Iₐᵢ·vᵢ + fᵢ`
   - Accumulate: `Iₐ_parent += X⁻¹·Iₐ_child·X⁻ᵀ`, `pₐ_parent += X⁻ᵀ·pₐ_child`
5. **Solve**: Compute joint accelerations: `q̈ᵢ = (τᵢ - Sᵢᵀ·pₐᵢ) / (Sᵢᵀ·Iₐᵢ·Sᵢ)` (`src/ForwardDynamics.cpp:96-125`)
6. **Output**: `Link::qddot` populated for each link

### Secondary Flow — Inverse Dynamics (RNEA)

1. **Setup**: User populates `InverseDynamics::links` vector with `InverseDynamicsLink` structs (`include/InverseDynamics.h:123`)
2. **Entry**: `id.computeTorques(qddot)` called with joint acceleration vector (`src/InverseDynamics.cpp:103`)
3. **Outward Pass (base→tip)**: Propagate velocities and accelerations (`src/InverseDynamics.cpp:19-52`):
   - Base: `v₀ = S₀·q̇₀`, `a₀ = S₀·q̈₀`
   - Child: `vᵢ = Xᵢ·v_parent + Sᵢ·q̇ᵢ`, `aᵢ = Xᵢ·a_parent + Sᵢ·q̈ᵢ + vᵢ × Sᵢ·q̇ᵢ`
4. **Inward Pass (tip→base)**: Compute spatial forces, accumulate child forces, project onto joint axes (`src/InverseDynamics.cpp:54-101`):
   - `fᵢ = Iᵢ·aᵢ + vᵢ × Iᵢ·vᵢ + Σ Xⱼ⁻ᵀ·fⱼ (children)`
   - `τᵢ = fᵢ·Sᵢ`
5. **Output**: Returns `Eigen::VectorXd tau` — joint torques

### Transform Flow — Plücker Coordinates

1. Create `PluckerTransform(rotation, translation)` — stores rotation and translation separately (`src/PluckerTransform.cpp:12-13`)
2. `transformMotion(mv)`: applies `X_m = [R, 0; -R·[t]×, R]` to a motion vector (`src/PluckerTransform.cpp:15-26`)
3. `transformForce(fv)`: applies `X_f = X_m⁻ᵀ` to a force vector (`src/PluckerTransform.cpp:28-38`)
4. `tformRBI(rbi)`: applies `I' = X·I·Xᵀ` for rigid body inertias (`src/PluckerTransform.cpp:64-81`)
5. `tformABI(abi)`: applies 6×6 block transform `Iₐ' = X·Iₐ·Xᵀ` for articulated body inertias (`src/PluckerTransform.cpp:104-161`)

**State Management:**
- No global state. All computation is local to the class instance.
- `ForwardDynamics` and `InverseDynamics` store link state as `std::vector<Link>` public members, mutated in-place during passes.
- Dynamics classes are reusable — callers update link state and re-invoke.

## Key Abstractions

**SpatialVector (6D vector):**
- Purpose: Base class representing a 6D spatial vector with angular and linear components
- Examples: `include/SpatialVector.h:68`
- Pattern: Concrete base class with protected `Vector3d angular` and `Vector3d linear` members. Provides `crossMotion()`, `crossForce()`, and `dot()` operations. Not intended for direct use — use `MotionVector` or `ForceVector` for type safety.

**MotionVector / ForceVector (physical specialization):**
- Purpose: Type-safe specialization for motion (twist) and force (wrench) duality
- Examples: `include/MotionVector.h:69`, `include/ForceVector.h:71`
- Pattern: Inheritance from `SpatialVector`. Each overrides operators to return their own type. `MotionVector` transforms covariantly under Plücker transforms; `ForceVector` transforms contravariantly.

**PluckerTransform (6×6 spatial transform):**
- Purpose: Rigid body coordinate transformation in Plücker coordinates
- Examples: `include/PluckerTransform.h:77`
- Pattern: Composes `Rotation rotation` and `Vector3d translation`. Provides `transformMotion()`, `transformForce()`, `inverse()`, `multiply()`, `tformRBI()`, `tformABI()`. All transforms use the 6×6 block matrix formulation from Featherstone.

**LowerTriangular (packed matrix):**
- Purpose: Memory-efficient storage for symmetric matrices (inertia tensors)
- Examples: `include/LowerTriangular.h:75`
- Pattern: 1D `Eigen::VectorXd` array with index mapping `idx = i*(i+1)/2 + j`. Provides `multiplySymmetric()` for symmetric matrix-vector multiplication, `getSymmetricMatrix()` for full reconstruction, OpenMP parallelization for matrix-matrix multiply.

**Link (kinematic chain node):**
- Purpose: Represents a single link in a kinematic tree for dynamics algorithms
- Examples: `include/ForwardDynamics.h:79` (ABA), `include/InverseDynamics.h:75` (RNEA)
- Pattern: Public struct with parent index, Plücker transform, rigid body inertia, joint motion axis, state (q, qdot, qddot), and intermediate quantities (v, a/c, f, Ia, pa). Parent index -1 indicates base link.

**Rotation (3D rotation):**
- Purpose: 3D rotation matrix with multiple representation support
- Examples: `include/Rotation.h:57`
- Pattern: Extends `Eigen::Matrix3d` via inheritance. Constructs from `AngleAxisd`, `Quaterniond`, or `Matrix3d`. Provides `inverse()`, `transpose()`, `toAngleAxis()`, `toQuaternion()`.

## Entry Points

**Library Headers (Primary):**
- Location: `include/*.h` (12 header files)
- Triggers: User `#include`s the desired header(s)
- Responsibilities: Declare all class interfaces. The library has **no umbrella header** — users include what they need individually.

**Library Source (Compilation):**
- Location: `src/*.cpp` (12 source files, plus `main.cpp`)
- Triggers: CMake compiles all `src/*.cpp` into `libSpatialAlgebra.a`
- Responsibilities: Define non-inline methods. Note: `RigidBodyInertia.cpp` and `ArticulatedBodyInertia.cpp` are empty stubs (all inline in headers).

**Test Executables:**
- Location: `tests/Test*.cpp` (9 test files, all with `main()` calling `RUN_ALL_TESTS()`)
- Triggers: Built and registered via CTest when `cmake --build build && ctest` is run
- Responsibilities: Verify correctness of all library classes

**Example Executables:**
- Location: `examples/basic_vectors.cpp`, `transforms.cpp`, `inertia.cpp`, `dynamics.cpp`
- Triggers: Built when cmake processes `examples/CMakeLists.txt`
- Responsibilities: Demonstrate library usage patterns

**Python RNEA:**
- Location: `robot_dynamics/rnea.py`
- Triggers: Python import or direct execution
- Responsibilities: Standalone inverse dynamics (not integrated with C++ library)

## Architectural Constraints

- **Threading:** Single-threaded by default. `LowerTriangular::operator*` uses `#pragma omp parallel for` for matrix-matrix multiplication if OpenMP is available. No other parallel constructs.
- **Global state:** None. All state is instance-local. The `using` type aliases (`mv`, `fv`, `plux`, `rbi`, `abi`, `lt`) at namespace scope in `SpatialAlgebra` are the only file-scope declarations.
- **Circular imports:** None detected. The dependency graph is a DAG: SpatialVector → MotionVector/ForceVector → PluckerTransform → RigidBodyInertia/ArticulatedBodyInertia → ForwardDynamics/InverseDynamics.
- **Memory model:** Value semantics throughout. No raw `new`/`delete`. No `std::shared_ptr` or `std::unique_ptr`. All vectors use `std::vector<T>` with value types.
- **No virtual methods:** The class hierarchy (`SpatialVector → MotionVector/ForceVector`) uses no virtual functions. Specialization is achieved via method overriding (non-virtual) and type-specific return types. This is a deliberate design choice for performance (no vtable overhead).

## Anti-Patterns

### Stub .cpp files
**What happens:** `src/RigidBodyInertia.cpp` and `src/ArticulatedBodyInertia.cpp` exist but contain only comments saying "No implementation needed — all methods are inline in header."
**Why it's wrong:** These files are compiled as part of `file(GLOB SOURCES "src/*.cpp")` in CMake, so they produce empty translation units. While harmless, they're misleading and unnecessary.
**Do this instead:** Remove these stub `.cpp` files and exclude them from the glob, or remove the glob pattern and list sources explicitly.

### Inconsistent include guard style
**What happens:** `LowerTriangular.h` uses `#pragma once` while all other headers use traditional `#ifndef`/`#define`/`#endif` guards.
**Why it's wrong:** Inconsistent guard style across the codebase. `#pragma once` is non-standard (though widely supported).
**Do this instead:** Standardize on `#ifndef`/`#define`/`#endif` guards matching the filename pattern (as done in all other headers).

### Non-virtual inheritance in class hierarchy
**What happens:** `MotionVector` and `ForceVector` inherit from `SpatialVector`, but no methods are virtual. Passing by base pointer/reference leads to static dispatch, not dynamic.
**Why it's wrong:** This breaks polymorphism if someone tries to use `SpatialVector&` to hold either type, but it's a conscious performance trade-off.
**Do this instead:** This is intentional — the hierarchy is for type safety and code reuse, not runtime polymorphism. Methods are overridden to return concrete types (`MotionVector`, `ForceVector`). Keep as-is but document the intent.

## Error Handling

**Strategy:** Exception-based error handling with input validation at public method boundaries.

**Patterns:**
- `std::invalid_argument` for dimension mismatches (`LowerTriangular` operations, `ForwardDynamics::computeAccelerations`, `InverseDynamics::computeTorques`)
- `std::out_of_range` for index bounds (debug mode only in `LowerTriangular`)
- `std::runtime_error` for singular matrix conditions (`LowerTriangular::inverse()`, `ForwardDynamics::inwardPass` near-zero inertia)
- `std::isnan`/`std::isinf` validation in dynamics solvers (`src/ForwardDynamics.cpp:143`, `src/InverseDynamics.cpp:117`)
- `assert()` used in `TestSpatialVector.cpp` (not in library code)

## Cross-Cutting Concerns

**Logging:** Every class has a `print()` method that writes to `std::cout`. No structured logging, no log levels, no log file output.

**Validation:** Minimal. Relies on caller correctness for most operations. Bounds checking only in debug mode via `#ifndef NDEBUG`. Input validation (NaN/Inf, dimension mismatch) exists in dynamics solver public entry points.

**Authentication:** Not applicable — this is a native compiled library with no network communication.

**Documentation:** All classes and methods have Doxygen `@brief`/`@details` comments. Generated via `doxygen Doxyfile` to `docs/html/` and `docs/latex/`.

---

*Architecture analysis: 2026-05-17*
