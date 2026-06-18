---
phase: 14-cr-02-bug-fix
reviewed: 2026-06-17T00:00:00Z
depth: standard
files_reviewed: 9
files_reviewed_list:
  - include/SpatialUtils.h
  - include/SpatialVector.h
  - src/SpatialVector.cpp
  - src/InverseDynamics.cpp
  - include/LowerTriangular.h
  - include/Rotation.h
  - tests/TestLowerTriangular.cpp
  - tests/TestRotation.cpp
  - tests/TestSpatialOperations.cpp
findings:
  critical: 1
  warning: 4
  info: 5
  total: 10
status: issues_found
---

# Phase 14: Code Review Report

**Reviewed:** 2026-06-17T00:00:00Z
**Depth:** standard
**Files Reviewed:** 9
**Status:** issues_found

## Summary

Reviewed 6 source files and 3 test files from the SpatialAlgebra C++17 library. The core spatial vector algebra (`SpatialVector`, `SpatialUtils`, `InverseDynamics`) is mathematically correct and well-implemented. However, a significant safety issue was found in the `LowerTriangular` class: release-mode builds skip bounds checking, enabling out-of-bounds memory access with negative or oversized indices. The `LowerTriangular::transpose()` method has a misleading interface that silently discards off-diagonal data. Several quality issues include missing validation, undocumented code, and dead test files that are not compiled by CMakeLists.

## Critical Issues

### CR-01: Debug-only bounds checking allows out-of-bounds memory access in release builds

**File:** `include/LowerTriangular.h:146–153` and `164–173`
**Issue:** The `operator()(int i, int j)` accessors guard index bounds (`i >= n || j >= n || i < 0 || j < 0`) with `#ifndef NDEBUG`. In release builds (`NDEBUG` defined), these checks are removed entirely. Calling `operator()` with negative or out-of-range indices on a release build invokes `getIndex(i, j)` which computes an invalid packed-data offset, leading to out-of-bounds read/write on the internal `Eigen::VectorXd data`. Concrete examples:
- `operator()(0, -1)` → `getIndex(0, -1)` = `(0*1)/2 + (-1)` = **−1** → `data[-1]` (crash/UB)
- `operator()(n+1, 0)` on a 3×3 matrix → `getIndex(4, 0)` = `(4*5)/2 + 0` = 10, but `data.size()` = 6 → `data[10]` (crash/UB)

**Fix:**
```cpp
// Always validate indices — remove the #ifndef NDEBUG guard:
inline double operator()(int i, int j) const
{
    if (i >= n || j >= n || i < 0 || j < 0)
        throw std::out_of_range(
            "LowerTriangular::operator() index out of bounds: (" +
            std::to_string(i) + "," + std::to_string(j) + ") for " +
            std::to_string(n) + "x" + std::to_string(n) + " matrix");
    return (i < j) ? 0.0 : data[getIndex(i, j)];
}

inline double &operator()(int i, int j)
{
    if (i >= n || j >= n || i < 0 || j < 0)
        throw std::out_of_range(
            "LowerTriangular::operator() index out of bounds: (" +
            std::to_string(i) + "," + std::to_string(j) + ") for " +
            std::to_string(n) + "x" + std::to_string(n) + " matrix");
    if (i < j)
        throw std::invalid_argument(
            "Cannot modify upper triangular element (" +
            std::to_string(i) + "," + std::to_string(j) + ")");
    return data[getIndex(i, j)];
}
```
If runtime performance in release builds is an absolute requirement and callers are trusted, consider adding `assert()` macros from `<cassert>` (which behave differently from `#ifndef NDEBUG` guards) or offering a separate unchecked `operator()_unsafe` accessor.

---

## Warnings

### WR-01: Constructor accepts negative size without validation

**File:** `include/LowerTriangular.h:115–119`
**Issue:** The `explicit LowerTriangular(int size)` constructor sets `n = size` and resizes `data` using `size * (size + 1) / 2` without checking if `size` is negative. A negative size produces a negative `n` but a non-negative (or zero) `data` size. For example, `size = -5` → `n = -5`, `data.resize(10)` → `getSize()` returns −5, but `data` has 10 elements. Subsequent `operator()` calls with valid 0-based indices will incorrectly pass the debug bounds check (`0 >= -5` = true) or access `data` at those indices in release builds, producing silent data corruption.

**Fix:**
```cpp
explicit LowerTriangular(int size) : n(size)
{
    if (size < 0)
        throw std::invalid_argument(
            "LowerTriangular size must be non-negative, got " +
            std::to_string(size));
    data.resize(size * (size + 1) / 2);
    data.setZero();
}
```

### WR-02: `transpose()` silently discards off-diagonal data; misleading method name

**File:** `include/LowerTriangular.h:393–404`
**Issue:** The `transpose()` method returns a `LowerTriangular` where off-diagonal elements are all zero. For a lower-triangular matrix L, the mathematical transpose Lᵀ is upper-triangular (all non-trivial elements are above the diagonal). Storing Lᵀ in a `LowerTriangular` representation means only the diagonal elements survive:
```
L = [1 0 0]      Lᵀ = [1 2 4]      lower(Lᵀ) = [1 0 0]
    [2 3 0]           [0 3 5]                  [0 3 0]
    [4 5 6]           [0 0 6]                  [0 0 6]
```
The implementation produces the right-hand result — the mathematically correct lower triangle of Lᵀ — but the method name `transpose()` is misleading because users expect the full transpose. The test (`TestLowerTriangular.cpp:254`) acknowledges this with the comment "The actual transpose would be upper triangular, but we store as lower."

**Fix:** Either:
- **(A)** Return `Eigen::MatrixXd` from `transpose()` so no information is lost:
  ```cpp
  Eigen::MatrixXd transpose() const
  {
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
      for (int i = 0; i < n; ++i)
          for (int j = 0; j <= i; ++j)
              result(j, i) = (*this)(i, j);  // write to upper triangle
      return result;
  }
  ```
- **(B)** Rename the method to `lowerOfTranspose()` or `transposeAsLower()` to accurately reflect semantics.
- **(C)** Add a `@warning` in the Doxygen comment explaining that the returned `LowerTriangular` contains only the lower triangle of Lᵀ (i.e., only diagonal elements from the original matrix survive).

### WR-03: `InverseDynamics::outwardPass()` uses naive ordering — fragile for non-topological link arrays

**File:** `src/InverseDynamics.cpp:22–52`
**Issue:** The outward pass iterates `for (int i = 0; i < ...; i++)` and accesses `links[parent].v` / `links[parent].a` assuming that the parent has already been processed. `InverseDynamics.h:117` documents: *"Links must be ordered such that parents appear before children."* There is no runtime validation of this constraint. If a user violates the ordering, the algorithm reads uninitialized `v`/`a` from the parent link, producing silent garbage results.

**Fix:** Add a pre-pass validation in `outwardPass()` that verifies `parent < i` for every link:
```cpp
void InverseDynamics::outwardPass()
{
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        int parent = links[i].parent;
        if (parent >= i && parent != -1)
            throw std::invalid_argument(
                "InverseDynamics::outwardPass: link " + std::to_string(i) +
                " has parent " + std::to_string(parent) +
                " — links must be topologically ordered (parent index < child index)");
        // ... rest of loop ...
    }
}
```

### WR-04: `InverseDynamics::computeTorques()` iteration variable signed/unsigned mismatch

**File:** `src/InverseDynamics.cpp:109`
**Issue:** The loop `for (int i = 0; i < qddot.size(); i++)` compares signed `int i` against `Eigen::Index` (typically `long`) from `qddot.size()`. This generates a compiler warning (`-Wsign-compare`) and can produce incorrect behavior if `qddot.size()` exceeds `INT_MAX` (theoretical for very large kinematic trees, though unlikely in practice).

**Fix:**
```cpp
for (Eigen::Index i = 0; i < qddot.size(); ++i)
```

---

## Info

### IN-01: Missing Doxygen comment on `dot(MotionVector, ForceVector)` overload

**File:** `include/SpatialUtils.h:65–68`
**Issue:** The `dot(const MotionVector&, const ForceVector&)` overload implements the physically meaningful power computation but lacks any `@brief`/`@details` Doxygen documentation, unlike all other `dot()` and `cross()` overloads in the same file.

**Fix:** Add a Doxygen block:
```cpp
/**
 * @brief Computes the spatial dot product (power) between a motion and force vector
 * @param v1 Motion vector (twist)
 * @param v2 Force vector (wrench)
 * @return Scalar power = ω·τ + v·f
 * @details This is the physically meaningful dot product: motion · force = power.
 */
inline double dot(const MotionVector& v1, const ForceVector& v2) noexcept { ... }
```

### IN-02: Misleading variable naming in `cross(MotionVector, ForceVector)`

**File:** `include/SpatialUtils.h:79–82`
**Issue:** The variable `f1` at line 81 is assigned `v2.getAngular()`. Since `v2` is a `ForceVector`, `getAngular()` returns the torque component τ₂, not a force component. The name `f1` (which denotes the linear force in spatial algebra notation) is misleading — the Doxygen formula uses τ₂ for this value. The math is correct; only the naming is confusing and risks future maintenance errors.

**Fix:** Rename variables to match the Doxygen formula:
```cpp
const Vector3d& w1 = v1.getAngular();
const Vector3d& v1_lin = v1.getLinear();
const Vector3d& tau2 = v2.getAngular();  // was f1
const Vector3d& f2 = v2.getLinear();
return ForceVector(
    w1.cross(tau2) + v1_lin.cross(f2),
    w1.cross(f2)
);
```

### IN-03: Three test files are dead code — not compiled or linked

**Files:** `tests/TestLowerTriangular.cpp`, `tests/TestRotation.cpp`, `tests/TestSpatialOperations.cpp`
**Issue:** `CMakeLists.txt` only registers two test executables: `TestSpatialVector` and `TestPluckerTransform`. These three test files contain comprehensive GTest suites but are never compiled, linked, or executed. Any bugs in `LowerTriangular`, `Rotation`, or `SpatialOperations` will not be caught by CI or manual testing since these tests never run.

**Fix:** Add the test executables to `CMakeLists.txt`:
```cmake
add_executable(TestLowerTriangular tests/TestLowerTriangular.cpp src/LowerTriangular.cpp)
target_link_libraries(TestLowerTriangular GTest::GTest GTest::Main)

add_executable(TestRotation tests/TestRotation.cpp src/Rotation.cpp)
target_link_libraries(TestRotation GTest::GTest GTest::Main)

add_executable(TestSpatialOperations tests/TestSpatialOperations.cpp ...)
target_link_libraries(TestSpatialOperations GTest::GTest GTest::Main)

add_test(NAME TestLowerTriangular COMMAND TestLowerTriangular)
add_test(NAME TestRotation        COMMAND TestRotation)
add_test(NAME TestSpatialOperations COMMAND TestSpatialOperations)
```

### IN-04: Typo in test tolerance constant name

**File:** `tests/TestSpatialOperations.cpp:24`
**Issue:** The tolerance constant is misspelled as `EPSILON` instead of `EPSILON`.

**Fix:** Rename to `EPSILON`.

### IN-05: `crossProductMotion` test loses type information via base class slicing

**File:** `tests/TestSpatialOperations.cpp:62`
**Issue:** `SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);` — the return value is stored as `SpatialVector` rather than `MotionVector`. Since `MotionVector` and `ForceVector` have distinct semantics (they transform differently under coordinate changes and their dot product has physical meaning), losing type information at the call site weakens type safety.

**Fix:** Store result as `MotionVector`:
```cpp
MotionVector result = SpatialOperations::crossProductMotion(v1, v2);
```
(The compiler will use the implicit `MotionVector(const SpatialVector&)` conversion — but better yet, `crossProductMotion` should be declared to return `MotionVector`.)

---

_Reviewed: 2026-06-17T00:00:00Z_
_Reviewer: the agent (gsd-code-reviewer)_
_Depth: standard_
