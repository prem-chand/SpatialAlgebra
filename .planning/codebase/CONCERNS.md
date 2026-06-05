# Codebase Concerns

**Analysis Date:** 2026-06-05

## Tech Debt

### Empty stub `.cpp` files excluded from build

**Issue:** `src/RigidBodyInertia.cpp` and `src/ArticulatedBodyInertia.cpp` are empty stubs explicitly excluded from compilation in `CMakeLists.txt:30-33`. The entire `RigidBodyInertia` and `ArticulatedBodyInertia` classes are defined inline in their headers, making the `.cpp` files dead code that invites confusion.

**Files:** `src/RigidBodyInertia.cpp`, `src/ArticulatedBodyInertia.cpp`, `CMakeLists.txt:30-33`

**Impact:** Developers may mistakenly add method definitions to these files without updating CMakeLists.txt, resulting in linker errors. The `list(REMOVE_ITEM ...)` hack in CMakeLists.txt is a maintenance trap.

**Fix approach:** Either delete the stub files and remove the exclusion from CMakeLists.txt, or move the class implementations out of the headers into these `.cpp` files.

---

### SpatialOperations is a near-empty delegation layer

**Issue:** `SpatialOperations` (declared in `include/SpatialOperations.h`, 35 lines) provides only 3 static methods that are thin wrappers:
- `crossProductMotion()` delegates to free function `cross()` in `SpatialUtils.h`
- `crossProductForce()` delegates to free function `cross()` in `SpatialUtils.h`
- `transformInertia()` delegates to `PluckerTransform::tformRBI()`

The class serves no purpose — it's an unnecessary indirection between callers and the actual implementations.

**Files:** `include/SpatialOperations.h:19-31`, `src/SpatialOperations.cpp:11-25`

**Impact:** Callers have two equivalent ways to call cross product operations (`SpatialOperations::crossProductMotion(a,b)` vs `cross(a,b)` from `SpatialUtils.h`), creating API ambiguity. The class adds conceptual weight with no value.

**Fix approach:** Deprecate or remove `SpatialOperations` and let callers use `SpatialUtils` free functions directly.

---

### Include guard inconsistency

**Issue:** Most headers use `#ifndef`/`#define`/`#endif` include guards (e.g., `SPATIAL_VECTOR_H`, `PLUCKER_TRANSFORM_H`). However, `LowerTriangular.h` uses a different style (`LOWER_TRIANGULAR_H` — no `#pragma once` actually, but the guard naming convention is inconsistent with the file naming pattern). The `AGENTS.md` mentions `LowerTriangular.h` previously used `#pragma once`, which may have been changed.

**Files:** `include/*.h`

**Impact:** Cosmetic inconsistency. No functional impact but suggests the codebase has been edited by multiple contributors with different conventions.

---

### ForwardDynamics uses raw 6x6 matrix conversion instead of API methods

**Issue:** The anonymous namespace helper `transformInertiaToParent()` in `src/ForwardDynamics.cpp:22-53` manually converts `ArticulatedBodyInertia` to a dense 6x6 matrix, constructs the Plücker transform as a 6x6 matrix, and performs the triple product `X^T * Ia * X`. It then extracts the blocks back into `ArticulatedBodyInertia` components. This duplicates the logic already present in `PluckerTransform::tformABI()` and `PluckerTransform::invtformABI()`.

**Files:** `src/ForwardDynamics.cpp:22-53`, `include/PluckerTransform.h:106-163`

**Impact:** Maintains two implementations of the same ABI transform logic. Any fix to one must be manually mirrored in the other. The local helper uses `X^T * Ia * X` while `tformABI()` uses `X * Ia * X^T` — the direction difference suggests one may be incorrect.

**Fix approach:** Replace the local helper with a call to the appropriate `PluckerTransform` method (likely `invtformABI()`), and delete the duplicated matrix construction code.

---

### Python RNEA not integrated with C++ library

**Issue:** `robot_dynamics/rnea.py` contains a standalone NumPy-based RNEA implementation that is completely disconnected from the C++ library. It defines its own `RigidBodyParams` and `MultiBodySystem` classes instead of using the C++ `RigidBodyInertia` and `Link` structures. The example usage at the bottom of the file is broken (`[link1, ...]` is not valid Python).

**Files:** `robot_dynamics/rnea.py:1-120`

**Impact:** Two separate RNEA implementations to maintain. The Python version won't benefit from C++ bug fixes and vice versa. The broken example code erodes confidence.

**Fix approach:** Either add Python bindings (pybind11) to the C++ library, or remove the orphaned Python file.

---

### Redundant definitions of the same arithmetic operators

**Issue:** `SpatialVector` defines `operator+`, `operator-`, `operator*`, `crossMotion`, `crossForce`, and `dot`. `MotionVector` and `ForceVector` each redundantly re-define the same arithmetic operators instead of using the base class. Many of these just call similarly-named methods on the base or are duplicates of the free functions in `SpatialUtils.h`.

**Files:** `src/MotionVector.cpp:16-29`, `src/ForceVector.cpp:16-29`, `src/SpatialVector.cpp:36-69`, `include/SpatialUtils.h:38-127`

**Impact:** The codebase has up to 3 different ways to compute the same spatial cross product. This violates DRY and increases the surface area for bugs.

---

## Known Bugs

### ForwardDynamicsLink default Ia is not zero

**Issue:** The `Link` struct default constructor (`include/ForwardDynamics.h:101-110`) initializes `Ia` (articulated body inertia) as:
```cpp
Ia(lt::Identity(3), Eigen::Matrix3d::Zero(), lt::Identity(3))
```
This is a non-zero articulated inertia (identity blocks on the rotation and mass sub-blocks), which is not physically correct for an uninitialized link. The `ArticulatedBodyInertia` default constructor (`include/ArticulatedBodyInertia.h:104-107`) correctly creates a zero inertia. The `Link` default should use the zero ABI constructor.

**Files:** `include/ForwardDynamics.h:109`, `include/ArticulatedBodyInertia.h:104-107`

**Trigger:** Any use of a `Link` without explicitly setting `Ia` will have incorrect initial articulated body inertia. The `inwardPass()` then overwrites `Ia` before use (line 89-93 of `ForwardDynamics.cpp`), which mitigates the issue but makes the default misleading.

**Impact:** Low (Ia is always overwritten in inwardPass before use). However, the incorrect default could mask bugs if code is refactored.

**Workaround:** `inwardPass()` explicitly initializes `Ia` at `ForwardDynamics.cpp:89-93`.

---

### PluckerTransform::multiply() rotation composition without re-orthogonalization

**Issue:** `PluckerTransform::multiply()` has a TODO comment at `src/PluckerTransform.cpp:236`:
```
// TODO: how to ensure product to 2 rotation matrices is still a rotation matrix upto finite precision?
```
Repeated multiplication of rotation matrices can cause the result to drift from SO(3) due to floating-point error accumulation.

**Files:** `src/PluckerTransform.cpp:236-237`

**Trigger:** Extended kinematic chains with many transform compositions (e.g., 100+ links). Each composition introduces O(ε) drift in orthogonality.

**Impact:** Low for typical robot models (6-50 links). Could cause instability in very long chains or iterative algorithms that repeatedly multiply transforms.

**Workaround:** None implemented. A re-orthogonalization step (e.g., SVD projection onto SO(3) or Gram-Schmidt on the rotation matrix columns) would fix this.

---

### InverseDynamics does not use joint position `q`

**Issue:** The `InverseDynamicsLink` struct stores `q` (joint position), but `InverseDynamics::outwardPass()` never reads it. The joint position is only needed in the Plücker transform `X`, which is already stored separately and does not automatically update when `q` changes. The `InverseDynamics::computeTorques()` documentation lists `q` as a state variable, but the algorithm ignores it for kinematics.

**Files:** `include/InverseDynamics.h:86`, `src/InverseDynamics.cpp:19-53`

**Impact:** If the user updates `q` but not the corresponding `X` (Plücker transform), the kinematics are wrong. This is a design flaw — the transform should be a function of `q`, but the library doesn't enforce this relationship.

---

## Security Considerations

### No input validation in release builds

**Risk:** NaN/Inf checking is only active under `#ifndef NDEBUG` (debug builds). Release builds have zero input validation. A NaN or Inf value propagating through spatial algebra computations can cause silent numerical corruption (NaN poisoning) without any warning or error.

**Files:**
- `include/SpatialVector.h` — no validation at all (base class)
- `src/SpatialVector.cpp:15-20` — `#ifndef NDEBUG` guard
- `include/RigidBodyInertia.h:90-94` — `#ifndef NDEBUG` guard
- `include/ArticulatedBodyInertia.h:183-187` — `#ifndef NDEBUG` guard

**Current mitigation:** NaN/Inf warnings printed to stderr in debug builds only.

**Recommendations:** Promote NaN/Inf checks to release builds by either removing the `#ifndef NDEBUG` guard or adding a library-level configuration macro (e.g., `SPATIAL_ALGEBRA_RELEASE_VALIDATION`). The `computeTorques()` and `computeAccelerations()` entry points already validate inputs unconditionally (`src/InverseDynamics.cpp:109-118`, `src/ForwardDynamics.cpp:195-201`), but lower-level class constructors and methods do not.

---

### Zero-mass edge case can throw or silently produce NaNs

**Risk:** `ForwardDynamics::inwardPass()` throws a `std::runtime_error` when the scalar articulated inertia `D = S^T·Ia·S` is near zero (`src/ForwardDynamics.cpp:118-124`). However, the `ZeroMassEdgeCase` test (`tests/TestForwardDynamics.cpp:676-704`) explicitly allows both exceptions and silent finite results — meaning the behavior for massless links is implementation-defined and could change silently.

**Files:** `src/ForwardDynamics.cpp:118-124`, `tests/TestForwardDynamics.cpp:676-704`

**Recommendations:** Define explicit behavior for degenerate inertias (mass = 0). Either always throw, or always handle gracefully with documented output semantics.

---

## Performance Bottlenecks

### Dense 6x6 matrix construction in ABI transforms

**Problem:** `PluckerTransform::tformABI()` (`src/PluckerTransform.cpp:106-163`) and `invtformABI()` (`src/PluckerTransform.cpp:165-215`) construct full 6x6 dense matrices (Eigen::MatrixXd) for both the ABI and the transform, perform a full 6x6 matrix triple product, then extract 3x3 blocks. The same pattern is duplicated in `ForwardDynamics.cpp:22-53`. This is O(6³)=216 multiply-adds when analytic formulas could do the same in ~O(3³) per block.

**Files:** `src/PluckerTransform.cpp:134-149`, `src/PluckerTransform.cpp:186-202`, `src/ForwardDynamics.cpp:31-52`

**Cause:** Laziness — the 6x6 matrix path is simpler to implement but computationally wasteful.

**Improvement path:** Derive and implement block-form formulas for `X * Ia * X^T` and `X^{-1} * Ia * X^{-T}` operating directly on the 3x3 blocks, avoiding the 6x6 construction entirely. This would give approximately a 2-3x speedup for inertia transforms.

---

### LowerTriangular::multiplySymmetric uses conditionals in inner loop

**Problem:** The generic `multiplySymmetric()` (`include/LowerTriangular.h:327-349`) has an `if (i >= j)` branch inside a double loop over `n x n`. For each element, it decides whether to read from the lower triangle or mirror from the upper. This defeats auto-vectorization and is ~2x slower than an explicit unrolled loop for the 3x3 case (which already exists separately at line 359-369).

**Files:** `include/LowerTriangular.h:327-349`

**Cause:** The generic implementation was written for readability, not performance.

**Improvement path:** The 3x3 specialization (lines 359-369) is already efficient — it should be the primary path since all current uses are 3x3. The generic version could be replaced with two explicit loops without branches.

---

### LowerTriangular::operator<< accesses elements through virtual dispatch overhead

**Problem:** The `operator<<` for LowerTriangular (line 447-458) calls `operator()(i,j)` for every element, which has an implicit branch for upper-triangular elements and bound-checking in debug mode. For an `n x n` dense output, this is `n²` function calls.

**Files:** `include/LowerTriangular.h:447-458`

**Improvement path:** Directly index into `data` using the `getIndex()` formula for output, or iterate over the packed array directly with formatted column breaks.

---

## Fragile Areas

### ForwardDynamics::inwardPass() — high complexity

**Files:** `src/ForwardDynamics.cpp:82-180`

**Why fragile:** This single method (~100 lines) handles three distinct phases:
1. Phase 1: Initialize `Ia` and `pa` from rigid body inertia (lines 84-102)
2. Phase 2: Backward pass — condense articulated inertias and propagate to parent (lines 104-158)
3. Phase 3: Forward pass — correct `qddot` and compute spatial accelerations (lines 160-179)

Phase 2 itself is the most complex part — it computes `IaS`, `D`, `u`, partial `qddot`, the condensation corrections (`inertiaCorr`, `HCorr`, `massCorr`), then transforms and propagates to the parent. There are 3 distinct index vectors (`links`, `D_store`, `Ia_unc`) plus intermediate states that must be stored for Phase 3.

**Safe modification:** Any change to the condensation formula or propagation logic must maintain consistency between `Ia_unc[i]` (stored before condensation) and the post-condensation propagated `Ia`. The correction step at line 174 that uses `Ia_unc[i].apply(a_prime)` is particularly sensitive.

**Test coverage:** The test suite covers single-link, two-link, branching, and gravity cases, but does not test the condensation logic with more than 3 links where the difference between pre-condensation `Ia_unc` and post-condensation `Ia` becomes non-trivial.

---

### PluckerTransform inherits includes from multiple layers

**Files:** `include/PluckerTransform.h:46-52`

**Why fragile:** `PluckerTransform.h` includes `RigidBodyInertia.h` and `ArticulatedBodyInertia.h`, but Plucker transforms are conceptually a lower-level concept than inertia. This creates a circular dependency risk (though currently avoided) and means any change to `RigidBodyInertia.h` triggers a recompile of everything that includes `PluckerTransform.h`.

**Details:** The `tformRBI()`, `invtformRBI()`, `tformABI()`, and `invtformABI()` methods are defined in `PluckerTransform.cpp`, so technically they only need the inertia class declarations (not definitions). The includes could be replaced with forward declarations.

---

### MotionVector/ForceVector access protected members of SpatialVector

**Files:** `src/MotionVector.cpp:18,28,35,41`, `src/ForceVector.cpp:18,28,35,41`

**Why fragile:** Both `MotionVector` and `ForceVector` derive from `SpatialVector` and directly access the protected members `angular` and `linear`. If the base class storage format changes (e.g., to a single `Eigen::Matrix<double, 6, 1>` instead of two `Vector3d`), all derived class implementations must be updated.

---

### No unique_ptr or ownership semantics for large objects

**Files:** `include/ForwardDynamics.h:137`, `include/InverseDynamics.h:123`

**Why fragile:** `ForwardDynamics` and `InverseDynamics` store `std::vector<Link>` and `std::vector<InverseDynamicsLink>` as public members. These `Link` structs contain `PluckerTransform`, `RigidBodyInertia`, `ArticulatedBodyInertia`, etc. by value. Copying a `ForwardDynamics` or `InverseDynamics` object performs a deep copy of every link, including all inertia and transform data. While this is fine for small chains, it's a performance trap for larger models and there's no move-semantics optimization anywhere.

---

## Scaling Limits

### Robot model complexity

**Current capacity:** Tested with serial chains up to 3 links and branching trees with 3 links.

**Limit:** The `ForwardDynamics::inwardPass()` uses `std::vector<double>` and `std::vector<ArticulatedBodyInertia>` for temporary storage (lines 106-107), sized to `links.size()`. With the naive O(n) ABA algorithm, complexity scales linearly, but the current implementation uses dense 6x6 matrices for each transform, not block-form optimizations.

**Scaling path:** Implement block-form ABI transforms (see Performance section). Add benchmark tests for 10, 50, 100+ link chains. Consider removing the pre-condensation state copy (`Ia_unc`) if the correction formula can be reformulated.

---

### LowerTriangular matrix dimension

**Current capacity:** Only used for 3x3 matrices throughout the codebase. The generic implementation supports arbitrary `n`.

**Limit:** The `inverse()` method has O(n³) complexity for an n×n matrix, with an O(n²) memory allocation on every call (`LowerTriangular result(n)` constructs a new `Eigen::VectorXd` of size `n(n+1)/2`). For n=3 this is negligible.

---

## Dependencies at Risk

### Eigen 5.0.1 vendored in repository (116 MB)

**Risk:** The entire Eigen 5.0.1 source tree is vendored at `eigen-5.0.1/` (116 MB on disk). This bloats the repository size significantly. The CMakeLists.txt uses `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` which can find either the system-installed Eigen or the vendored copy, but there's no explicit path setup to prefer the vendored copy.

**Files:** `eigen-5.0.1/`, `CMakeLists.txt:12`

**Impact:** Large clone sizes, longer `git clone` times, and potential confusion about which Eigen version is actually used. The `build-eigen5/` directory (present in directory listing) suggests a separate build was made for testing Eigen 5 compatibility, but this build configuration is not documented.

**Migration plan:** Remove the vendored `eigen-5.0.1/` directory, document the required system dependency (`brew install eigen`), or use CMake's `FetchContent` to download Eigen at build time if not found.

---

### No CI pipeline

**Risk:** There are no GitHub Actions workflows for CI (`ci.yml` exists but may not be active or comprehensive). The codebase has no automated build verification, test runner, or static analysis.

**Files:** `.github/workflows/ci.yml`

**Impact:** Bugs can be introduced and only caught on the developer's machine. No regression detection. No enforced coding standards.

---

## Missing Critical Features

### No Joint abstraction

**Problem:** Joint types are represented by a `MotionVector S` (screw axis) and a `double q`/`qdot`/`qddot`. There is no `Joint` class that encapsulates joint type (revolute, prismatic, spherical), limits, friction, or actuation mode. The Plücker transform `X` must be manually set by the user for each joint position — it is not computed automatically from `q`.

**Files:** `include/ForwardDynamics.h:79-111`, `include/InverseDynamics.h:75-100`

**Blocks:** Realistic robot simulation where joint types beyond simple 1-DOF revolute/prismatic are needed.

---

### No robot model loader

**Problem:** No URDF, SDF, or any standard robot model format importer. Users must manually create `Link` structs and set transform/inertia/joint data by hand.

**Blocks:** Using the library with real robot models from standard sources.

---

### No benchmark/performance tests

**Problem:** The test suite has 13 test executables with functional tests, but there are no benchmarks measuring execution time for any operation. This makes it impossible to detect performance regressions or evaluate optimization improvements.

**Files:** `tests/*.cpp`

---

## Test Coverage Gaps

### Forward/Inverse dynamics validated only for simple chains

**What's not tested:** The round-trip consistency tests (`TestDynamicsConsistency.cpp`) only use single-link chains. The branching tree tests only verify positive/finite output, not numerical correctness. Condensation behavior in the ABA inward pass is not explicitly validated.

**Files:** `tests/TestDynamicsConsistency.cpp:16-93`, `tests/TestForwardDynamics.cpp:51-131`

**Risk:** The core ABA algorithm (`inwardPass`) has the highest complexity in the codebase but the weakest validation. Bugs in inertia condensation or parent propagation could go undetected.

**Priority:** High

---

### LowerTriangular::inverse() limited coverage

**What's not tested:** The `inverse()` method in `LowerTriangular.cpp:25-53` has no dedicated test coverage that verifies `L * L^{-1} = I`. The test file `TestLowerTriangular.cpp` focuses on storage, indexing, and matrix-vector multiplication.

**Files:** `src/LowerTriangular.cpp:25-53`, `tests/TestLowerTriangular.cpp`

**Risk:** The inverse computation uses forward substitution with three nested loops. Any off-by-one error in the index arithmetic would produce silently wrong results.

**Priority:** Medium

---

### No tests for `PluckerTransform::tformABI()` and `invtformABI()`

**What's not tested:** `TestPluckerTransform.cpp` does not test the `tformABI()` and `invtformABI()` methods. These are the newest and most complex transforms in the file (`src/PluckerTransform.cpp:106-215`).

**Files:** `src/PluckerTransform.cpp:106-215`, `tests/TestPluckerTransform.cpp`

**Risk:** These methods build full 6x6 matrices and perform triple products. Without dedicated tests, any bug here would only be caught indirectly through `ForwardDynamics` tests.

**Priority:** Medium

---

### No NaN/Inf propagation tests for release builds

**What's not tested:** All NaN/Inf guards are behind `#ifndef NDEBUG`. There are no tests that verify behavior when NaN/Inf inputs enter the library in release mode.

**Files:** `src/SpatialVector.cpp:15-20`, `include/RigidBodyInertia.h:90-94`, `include/ArticulatedBodyInertia.h:183-187`

**Risk:** NaN values can silently propagate through matrix operations in release builds, corrupting downstream results.

**Priority:** Medium

---

### ZeroMassEdgeCase test non-deterministic

**What's not tested:** The test `ForwardDynamicsTest.ZeroMassEdgeCase` (`TestForwardDynamics.cpp:676-704`) accepts two possible outcomes (exception OR finite result), making it a weak test. The specific behavior for degenerate inertias is not defined.

**Files:** `tests/TestForwardDynamics.cpp:676-704`

**Priority:** Low

---

*Concerns audit: 2026-06-05*
