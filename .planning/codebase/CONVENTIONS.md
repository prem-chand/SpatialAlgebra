# Coding Conventions

**Analysis Date:** 2026-06-05

## Naming Patterns

**Files:**
- PascalCase with descriptive names: `SpatialVector.h`, `PluckerTransform.cpp`, `RigidBodyInertia.h`
- Test files uniformly prefixed with `Test`: `TestSpatialVector.cpp`, `TestRotation.cpp`
- Header guards match filename: `SPATIAL_VECTOR_H` (in `SpatialVector.h`), `PLUCKER_TRANSFORM_H`, `ROTATION_H`
- Exception: `LowerTriangular.h` uses `#pragma once` instead of include guards

**Classes:**
- PascalCase reflecting the mathematical concept: `SpatialVector`, `MotionVector`, `ForceVector`, `PluckerTransform`, `Rotation`, `RigidBodyInertia`, `ArticulatedBodyInertia`, `LowerTriangular`
- `SpatialOperations` is a static utility class (no instances, all static methods)

**Functions:**
- camelCase for methods: `transformMotion`, `crossProductMotion`, `getAngular`, `setFromAngleAxis`, `computeTorques`, `computeAccelerations`
- Verbs for operations: `transform`, `apply`, `inverse`, `transpose`, `dot`, `cross`, `print`
- Getters prefixed with `get`: `getMass`, `getCom`, `getInertia`, `getAngular`, `getLinear`, `getData`, `getH`, `getM`
- `SpatialUtils.h` free functions are lowercase: `skew()`, `dot()`, `cross()`, `crossSpatial()`
- Free function dot product is overloaded: `dot(const MotionVector&, const MotionVector&)`, `dot(const ForceVector&, const ForceVector&)` (see `SpatialUtils.h:64-82`)

**Variables:**
- Private members: lowercase (e.g., `angular`, `linear`, `rotation`, `translation`, `mass`, `com`, `data`, `size`)
- Underscore suffix for some private members: `data_`, `size_`, `rotation_`, `translation_` (inconsistent — `SpatialVector.h:77` uses `angular`, `linear` while `LowerTriangular.h:87` uses `data_`, `size_`)
- Local variables: camelCase (e.g., `transformedAngular`, `newRotation`, `invRotation`, `twist1`, `dotProduct`)
- Test variables: descriptive of the scenario (e.g., `mv1`, `fv1`, `wrench1`, `qddot_input`, `tau_output`)

**Type Aliases:**
- `using Vector3d = Eigen::Matrix<double, 3, 1>` at file scope in headers
- `using Vector6d = Eigen::Matrix<double, 6, 1>` in `SpatialVector.h:13`
- Namespace-level shorthand aliases in `SpatialVector.h:18-21`:
  ```cpp
  using mv = MotionVector;
  using fv = ForceVector;
  using plux = PluckerTransform;
  using rbi = RigidBodyInertia;
  using abi = ArticulatedBodyInertia;
  using lt = LowerTriangular;
  ```
- `lt` is both a type alias and used as a variable name for LowerTriangular instances in test code
- `InverseDynamicsLink` and `ForwardDynamicsLink` are struct types (not aliases), defined in respective headers

## Code Style

**Formatting:**
- Indentation: 4 spaces (no tabs)
- Line length: ~100-120 characters typical
- Braces: K&R style for functions (`namespace {` on same line), Allman style for classes (opening brace on next line with colon)
- No formatter configured — no `.clang-format` in project root (only `eigen-5.0.1/.clang-format` exists for the vendored Eigen)

**Linting:**
- No linter configured
- No `.clang-tidy` file in project root
- No CI workflow that enforces style checks

**Const Correctness:**
- Method parameters passed by `const&` for objects: `const SpatialVector &other`, `const Vector3d &angular`
- Methods marked `const` where appropriate: all getters, `dot()`, `print()` (though some older methods may miss this)
- Return by const reference for internal data: `const Vector3d &getCom() const` (`RigidBodyInertia.h:57`), `const Eigen::VectorXd& getData() const` (`LowerTriangular.h:107`)

**noexcept:**
- Used sparingly. `Rotation.h` methods like `transpose()`, `inverse()` are `noexcept`. `SpatialVector.h` operators (`+`, `-`, `*`, etc.) are `noexcept`.
- Not consistently applied across all classes — `PluckerTransform` constructors and methods omit `noexcept`.

**`explicit`:**
- Used on single-argument constructors: `explicit SpatialVector(const Vector3d &angular)` (`SpatialVector.h:82`), `explicit Rotation(const Eigen::Matrix3d &m)` (`Rotation.h:43`)
- Not used on multi-argument constructors or default constructors

## Import Organization

**Include Guard Pattern:**
```cpp
#ifndef SPATIAL_VECTOR_H
#define SPATIAL_VECTOR_H
// ... header content ...
#endif // SPATIAL_VECTOR_H
```
(Exception: `LowerTriangular.h` uses `#pragma once`.)

**Include Order in Headers:**
1. Standard library / Eigen headers (alphabetically)
2. Project headers (alphabetically)
```cpp
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "ClassName.h"
#include "OtherDependency.h"
```

**Include Order in Source Files:**
1. Corresponding header first (for self-containment check)
2. Standard library / Eigen headers
3. Other project headers
4. GTest headers (in test files)
```cpp
#include "ClassName.h"

#include <Eigen/Dense>
#include <iostream>

#include "Dependency.h"
#include "OtherDependency.h"
```

**Test file includes:**
```cpp
#include "ClassName.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
```

**Path Aliases:**
- None configured. All includes use bare filenames resolved by `include_directories(include)` in CMake

## Error Handling

**Strategy:** C++ exceptions for runtime errors, assertions for debug-mode validation.

**Patterns:**
- `std::invalid_argument` for dimension mismatches: `LowerTriangular` operations throw when matrix dimensions don't match (`LowerTriangular.cpp:45-48`)
- `std::runtime_error` for singular matrices in `Rotation::inverse()` (`Rotation.cpp:37-42`)
- Debug-mode NaN/Inf checking: constructors like `SpatialVector(const Vector3d &angular)` check `hasNaN()` and `hasInf()` via `#ifndef NDEBUG` guards with `std::cerr` warnings (`SpatialVector.cpp:7-19`)
- Bounds checking only in debug mode: `LowerTriangular::operator()` checks index bounds via `assert()` inside `#ifndef NDEBUG` guards
- No custom exception types defined
- No error codes or `std::optional` / `std::expected` return types
- No input validation at public API boundaries (caller is expected to pass valid data)

## Logging

**Approach:** `print()` method on every class for human-readable output.

**Patterns:**
```cpp
void SpatialVector::print() const {
    std::cout << "Angular: [" << angular[0] << ", " << angular[1] << ", " << angular[2] << "]\n";
    std::cout << "Linear:  [" << linear[0] << ", " << linear[1] << ", " << linear[2] << "]\n";
}
```
- Present on `SpatialVector` (`SpatialVector.h:145`), `MotionVector`, `ForceVector`, `PluckerTransform` (`PluckerTransform.h:131`), `Rotation` (`Rotation.h:85`), `LowerTriangular` (`LowerTriangular.h:107`)
- Uses `std::cout` directly — no structured logging, no log levels, no log file output
- `print()` is not `const` in all cases

## Comments

**Doxygen Style:**
- Block comments (`/** ... */`) on every class and method declaration
- Tags used: `@brief`, `@details`, `@param`, `@return`, `@note`, `@warning`, `@see`, `@code`
- Example:
  ```cpp
  /**
   * @brief Transform a motion vector (twist) to a new coordinate frame
   * @details Applies the Plücker transform X = [R, 0; -R[r]x, R] to the
   *          input motion vector. The transformation formula is:
   *          ω' = R·ω
   *          v' = R·(v - r×ω)
   * @param input The motion vector to transform
   * @return Transformed motion vector in the new coordinate frame
   */
  MotionVector transformMotion(const MotionVector &input) const;
  ```
- Test files use Doxygen `@brief` and `@details` blocks above each `TEST()` or `TEST_F()`
- Inline comments (``) explain mathematical formulas and intermediate results
- Mathematical notation preserved in comments: `[ω; v]`, `[τ; f]`, `[R, 0; -R[t]x, R]`
- References to Featherstone textbook: `@see Featherstone 2008, Chapter 2`

**When to Comment:**
- Every class and public method gets a Doxygen block
- Complex formulas get inline explanation with intermediate value calculations
- Test cases explain the scenario and expected result

## Function Design

**Size:** Small to medium functions (10-50 lines typical). Operators are 3-10 lines. Complex transforms are 15-30 lines.

**Parameters/Return:**
- Pass objects by `const&`: `const SpatialVector &other`
- Pass primitives by value: `double scalar`, `int`
- Return by value for small objects (SpatialVector, Rotation, PluckerTransform) — copy elision / RVO expected
- Return by `const&` for internal data members: `const Vector3d &getCom() const`
- No `std::optional` or `std::variant` usage anywhere
- Named parameters not supported (C++ limitation)

**Inline methods:**
- Trivial getters and operators inlined in headers: `getAngular()`, `getLinear()`, `operator+`, `operator-`, `operator*`
- Complex logic in `.cpp` files: `PluckerTransform::transformMotion`, `Rotation::inverse()`
- `RigidBodyInertia` is entirely inline in header (`RigidBodyInertia.h`)
- Free functions in `SpatialUtils.h` are all inline in the header

**Static free functions:**
- `SpatialOperations` is a class with all static methods: `SpatialOperations::crossProductMotion()`, `SpatialOperations::crossProductForce()`, `SpatialOperations::transformInertia()`
- `SpatialUtils.h` provides free functions in `namespace SpatialAlgebra`: `skew()`, `dot()`, `cross()`, `crossSpatial()`

## Module Design

**Namespace:** All classes in `namespace SpatialAlgebra`. Free functions also in `namespace SpatialAlgebra`.

**Exports:**
- No umbrella include pattern (though `SpatialAlgebra.h` exists that includes all public headers in dependency order)
- Users include specific headers or `SpatialAlgebra.h`

**Source/Header Split:**
- Declarations in `.h` files in `include/`
- Definitions in `.cpp` files in `src/`
- Inline methods in headers for trivial operations (getters, operators)
- `RigidBodyInertia` is entirely inline in header

**Type Alias Pattern:**
- Aliases at file scope in `SpatialVector.h:18-21`
- Used consistently throughout the codebase: `mv`, `fv`, `plux`, `rbi`, `abi`, `lt`

## Documentation

**Config:** `Doxyfile` in project root

**Build:** `doxygen Doxyfile` generates HTML docs into `docs/html/` and LaTeX into `docs/latex/`

**Output:** `docs/html/`, `docs/latex/` (both gitignored)

**Tags used:** `@brief`, `@details`, `@param`, `@return`, `@note`, `@warning`, `@see`, `@code`

## Mathematical Notation

- Greek letters transliterated: `omega`, `tau` (not ω, τ) in variable names
- Formulas in code comments with Greek notation: `[ω; v]`, `[τ; f]`
- Code uses `angular` / `linear` for components of spatial vectors
- `crossMotion` formula: `[ω1×ω2; ω1×v2 + v1×ω2]`
- `crossForce` formula: `[τ1×τ2 + f1×f2; τ1×f2 - τ2×f1]`
- Plücker transform: `X = [R, 0; -R[r]x, R]`
- Test files annotate intermediate hand calculations to verify results

---

*Convention analysis: 2026-06-05*
