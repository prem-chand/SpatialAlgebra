# Coding Conventions

**Analysis Date:** 2026-05-15

## Naming Patterns

**Files:**
- PascalCase with descriptive names: `SpatialVector.h`, `PluckerTransform.cpp`
- Test files prefixed with `Test`: `TestSpatialVector.cpp`
- Header guards match filename: `SPATIAL_VECTOR_H`, `PLUCKER_TRANSFORM_H`

**Classes:**
- PascalCase: `SpatialVector`, `MotionVector`, `ForceVector`, `PluckerTransform`
- Descriptive of mathematical concept being modeled

**Functions:**
- camelCase: `transformMotion`, `crossProductMotion`, `getAngular`, `setFromAngleAxis`
- Verbs for operations: `transform`, `apply`, `inverse`, `transpose`
- Getters prefixed with `get`: `getMass`, `getCom`, `getInertia`

**Variables:**
- Private members: lowercase (e.g., `angular`, `linear`, `rotation`, `translation`, `mass`, `com`)
- Local variables: camelCase (e.g., `transformedAngular`, `newRotation`, `invRotation`)
- Eigen types: `Vector3d`, `Matrix3d`, `VectorXd` (using Eigen's typedefs)

**Types:**
- Type aliases at file scope: `using Vector3d = Eigen::Matrix<double, 3, 1>`
- Namespace aliases: `using mv = MotionVector`, `using fv = ForceVector`
- Lowercase abbreviations for concise notation: `mv`, `fv`, `plux`, `rbi`, `abi`, `lt`

## Code Style

**Formatting:**
- Indentation: 4 spaces (no tabs observed)
- Line length: ~100-120 characters typical
- Braces: K&R style for functions, Allman style for classes

**Linting:**
- No linter configured
- No formatter configured (no .clang-format, .prettierrc, etc.)

**Include Order:**
1. Corresponding header first (e.g., `#include "SpatialVector.h"` in `SpatialVector.cpp`)
2. Standard library headers
3. Eigen headers
4. Project headers

## Import Organization

**Header Guards:**
- Traditional `#ifndef`/`#define`/`#endif` pattern
- Example: `#ifndef SPATIAL_VECTOR_H` → `#define SPATIAL_VECTOR_H`
- Exception: `LowerTriangular.h` uses `#pragma once`

**Include Patterns:**
```cpp
// In headers
#include <Eigen/Dense>
#include "Dependency.h"

// In sources
#include "ClassName.h"      // Corresponding header first
#include "OtherDependency.h"
#include <Eigen/Dense>
#include <iostream>
```

**Path Aliases:**
- None configured in CMake
- VSCode config specifies: `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`

## Error Handling

**Patterns:**
- Exceptions for runtime errors:
  ```cpp
  throw std::invalid_argument("Matrix size mismatch");
  throw std::out_of_range("Index out of bounds");
  ```
- Debug-mode assertions via `#ifndef NDEBUG` guards
- No error codes or `std::optional` return types
- No custom exception types

**Validation:**
- Minimal input validation
- Relies on caller correctness for most operations
- Bounds checking only in debug mode

## Logging

**Framework:** None (direct `std::cout` usage)

**Patterns:**
- `print()` method on every class:
  ```cpp
  void SpatialVector::print() const {
      std::cout << "Angular: " << angular.transpose() << std::endl;
      std::cout << "Linear: " << linear.transpose() << std::endl;
  }
  ```
- No structured logging
- No log levels
- No log file output

## Comments

**When to Comment:**
- Every class and method has Doxygen documentation
- Inline comments explain mathematical formulas
- Complex operations have detailed `@details` sections

**JSDoc/TSDoc (Doxygen):**
- Block comments on all declarations:
  ```cpp
  /**
   * @brief Brief description
   * @param name Parameter description
   * @return Return value description
   * @details Extended explanation with formulas
   * @note Important notes
   * @warning Warnings about usage
   */
  ```
- Mathematical notation in comments (e.g., `[ω; v]`, `[τ; f]`)
- References to literature (Featherstone textbook)

## Function Design

**Size:**
- Small to medium functions (10-50 lines typical)
- Single responsibility per function
- Inline methods for trivial getters in headers

**Parameters:**
- Pass by const reference for objects: `const SpatialVector &other`
- Pass by value for primitives: `double scalar`, `int size`
- Named parameters not used (C++ limitation)

**Return Values:**
- Return by value for small objects (SpatialVector, Rotation)
- Return by const reference for internal data: `const Vector3d &getCom() const`
- No `std::optional` or `std::variant` usage

## Module Design

**Exports:**
- All classes in `namespace SpatialAlgebra`
- Type aliases at namespace scope: `using mv = MotionVector`
- Free functions in `namespace SpatialAlgebra` (SpatialUtils.h)

**Barrel Files:**
- None (users include specific headers)
- No `include/SpatialAlgebra.h` umbrella header

**Header-Only vs Separation:**
- Most methods declared in headers, defined in `.cpp` files
- Inline methods in headers for trivial operations (getters, operators)
- `RigidBodyInertia` entirely inline in header

## Documentation

**Doxygen:**
- Config: `Doxyfile`
- Command: `doxygen Doxyfile`
- Output: `docs/html/`, `docs/latex/`
- Tags used: `@brief`, `@details`, `@param`, `@return`, `@note`, `@warning`, `@see`, `@code`

**Example Documentation:**
```cpp
/**
 * @brief Compute the motion cross product
 * @param other The spatial vector to cross with
 * @return The resulting spatial vector
 * @details Implements the motion cross product operation:
 *          [ω1×ω2; ω1×v2 + v1×ω2]
 * @note This operation is specific to motion vectors (twists)
 */
SpatialVector crossMotion(const SpatialVector &other) const;
```

## Mathematical Notation

**In Code:**
- Greek letters transliterated: `omega`, `tau` (not ω, τ)
- Vector notation: `Vector3d`, `Vector6d`
- Matrix notation: `Matrix3d`, `MatrixXd`

**In Comments:**
- Mathematical notation preserved: `[ω; v]`, `[τ; f]`
- Formulas in code comments: `X = [R, 0; -R[t]x, R]`
- References to Featherstone textbook

---

*Convention analysis: 2026-05-15*
