# Architecture

**Analysis Date:** 2026-05-15

## Pattern Overview

**Overall:** Object-oriented spatial algebra library with inheritance-based type hierarchy

**Key Characteristics:**
- CRTP-like inheritance for spatial vectors (SpatialVector → MotionVector/ForceVector)
- Composition over inheritance for transforms (PluckerTransform contains Rotation + translation)
- Value semantics (objects passed by value/const reference)
- Eigen3 integration via inheritance (Rotation extends Eigen::Matrix3d)

## Class Hierarchy

**SpatialVector Base:**
```
SpatialVector (include/SpatialVector.h)
├── MotionVector (include/MotionVector.h) - alias: mv
└── ForceVector (include/ForceVector.h) - alias: fv
```

**Transform & Inertia:**
```
PluckerTransform (include/PluckerTransform.h) - alias: plux
├── Contains: Rotation, Vector3d (translation)
└── Transforms: MotionVector, ForceVector, RigidBodyInertia, ArticulatedBodyInertia

RigidBodyInertia (include/RigidBodyInertia.h) - alias: rbi
├── Contains: mass (double), com (Vector3d), inertiaMatrixLT (LowerTriangular)
└── apply(MotionVector) → ForceVector

ArticulatedBodyInertia (include/ArticulatedBodyInertia.h) - alias: abi
├── Contains: Inertia (LowerTriangular), H (Matrix3d), M (LowerTriangular)
└── apply(MotionVector) → ForceVector
```

**Utilities:**
```
Rotation (include/Rotation.h) - extends Eigen::Matrix3d
LowerTriangular (include/LowerTriangular.h) - packed storage matrix
SpatialOperations (include/SpatialOperations.h) - static utility class
SpatialUtils (include/SpatialUtils.h) - free functions (skew, dot, cross)
```

## Layers

**Core Algebra Layer:**
- Purpose: Fundamental spatial vector operations
- Location: `include/SpatialVector.h`, `src/SpatialVector.cpp`
- Contains: SpatialVector base class with angular/linear components
- Depends on: Eigen3
- Used by: All higher-level classes

**Motion/Force Layer:**
- Purpose: Type-safe motion (twist) and force (wrench) representations
- Location: `include/MotionVector.h`, `include/ForceVector.h`
- Contains: Specialized spatial vectors with physical interpretation
- Depends on: SpatialVector
- Used by: PluckerTransform, inertia classes

**Transform Layer:**
- Purpose: Coordinate frame transformations in Plücker coordinates
- Location: `include/PluckerTransform.h`, `src/PluckerTransform.cpp`
- Contains: 6x6 spatial transforms (rotation + translation)
- Depends on: Rotation, SpatialVector, RigidBodyInertia, ArticulatedBodyInertia
- Used by: Dynamics algorithms

**Inertia Layer:**
- Purpose: Mass property representations
- Location: `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`
- Contains: Mass, COM, inertia tensor (packed storage)
- Depends on: LowerTriangular, SpatialVector
- Used by: Dynamics algorithms

**Utility Layer:**
- Purpose: Helper functions and specialized data structures
- Location: `include/SpatialUtils.h`, `include/LowerTriangular.h`, `include/Rotation.h`
- Contains: skew(), dot(), cross(), packed matrix storage
- Depends on: Eigen3
- Used by: All layers

## Data Flow

**Motion Transformation:**
1. Create MotionVector with angular/linear velocity
2. Create PluckerTransform with Rotation + translation
3. Call `transformMotion()` → returns transformed MotionVector
4. Formula: `v' = [E    0; -E*r̂  E] * [ω; v]`

**Force Transformation:**
1. Create ForceVector with torque/force
2. Create PluckerTransform
3. Call `transformForce()` → returns transformed ForceVector
4. Formula: `f' = [E    E*r̂; 0    E] * [τ; f]`

**Inertia Transformation:**
1. Create RigidBodyInertia with mass, COM, inertia tensor
2. Create PluckerTransform
3. Call `tformRBI()` → returns transformed inertia
4. Formula: `I' = X * I * X^T`

**RNEA (Python):**
1. Forward recursion: propagate velocities/accelerations
2. Backward recursion: propagate forces/moments
3. Project onto joint axes → joint torques

## Key Abstractions

**SpatialVector:**
- Purpose: 6D vector combining angular + linear components
- Examples: `include/SpatialVector.h:68`
- Pattern: Base class with protected angular/linear Vector3d members

**PluckerTransform:**
- Purpose: Rigid body coordinate transformation
- Examples: `include/PluckerTransform.h:79`
- Pattern: Stores rotation and translation separately, provides motion/force transform methods

**LowerTriangular:**
- Purpose: Memory-efficient packed storage for symmetric matrices
- Examples: `include/LowerTriangular.h:70`
- Pattern: 1D array storage with index mapping: `idx = i*(i+1)/2 + j`

## Entry Points

**Library Usage:**
- Location: `include/*.h` (header files)
- Triggers: User includes headers and instantiates classes
- Responsibilities: Provide spatial algebra operations

**Test Executables:**
- Location: `tests/TestSpatialVector.cpp`, `tests/TestPluckerTransform.cpp`
- Triggers: Manual execution via `build/TestSpatialVector`
- Responsibilities: Verify core functionality

**Main Demo:**
- Location: `src/main.cpp`
- Triggers: `cmake --build build` produces executable
- Responsibilities: Demonstrate library usage

## Error Handling

**Strategy:** Exception-based with debug-mode assertions

**Patterns:**
- `std::invalid_argument` for dimension mismatches (LowerTriangular operations)
- `std::out_of_range` for index bounds (debug mode only)
- `assert()` for basic validation (TestSpatialVector.cpp)
- No error codes or result types

## Cross-Cutting Concerns

**Logging:** `print()` methods on all classes output to `std::cout`

**Validation:** Minimal runtime validation; relies on caller correctness

**Documentation:** Doxygen comments on all declarations (`@brief`, `@details`, `@param`, `@return`)

**Type Safety:** Separate MotionVector/ForceVector classes prevent mixing motion and force operations

---

*Architecture analysis: 2026-05-15*
