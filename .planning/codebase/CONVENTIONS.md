# Coding Conventions

**Analysis Date:** 2026-05-17

## Naming Patterns

**Files:**
- PascalCase for class header/source files: `SpatialVector.h`, `PluckerTransform.cpp`, `Rotation.h`
- Test files prefixed with `Test`: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`
- Utility headers named descriptively: `SpatialUtils.h`, `SpatialOperations.h`
- Python script: `rnea.py` (snake_case)

**Classes:**
- PascalCase, descriptive of the mathematical concept: `SpatialVector`, `MotionVector`, `ForceVector`, `PluckerTransform`, `RigidBodyInertia`, `ArticulatedBodyInertia`, `LowerTriangular`, `Rotation`

**Functions:**
- camelCase for methods: `transformMotion`, `crossProductMotion`, `getAngular`, `setFromAngleAxis`
- Verbs for operations: `transform`, `apply`, `inverse`, `transpose`, `multiplySymmetric`
- Getters prefixed with `get`: `getMass`, `getCom`, `getInertia`, `getAngular`, `getLinear`
- Static methods also camelCase: `SpatialOperations::crossProductMotion`, `LowerTriangular::fromFullMatrix`

**Variables:**
- Private members: lowercase (e.g., `angular`, `linear`, `rotation`, `translation`, `mass`, `com`, `Inertia`, `H`, `M`, `data`, `n`)
- Local variables: camelCase (e.g., `transformedAngular`, `newRotation`, `invRotation`, `combinedInertia`, `a_cross_b`)
- Parameters: camelCase (e.g., `const SpatialVector &other`, `const Rotation &rotation`)
- Loop counters: `i`, `j`, `k`

**Types:**
- Eigen types use Eigen's typedefs: `Vector3d`, `Matrix3d`, `VectorXd`, `Vector6d`
- File-scope `using` declarations for concise notation:
```cpp
using Vector3d = Eigen::Matrix<double, 3, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
```
- Namespace-level type aliases for shorthand:
```cpp
using mv = MotionVector;
using fv = ForceVector;
using plux = PluckerTransform;
using rbi = RigidBodyInertia;
using abi = ArticulatedBodyInertia;
using lt = LowerTriangular;
```
- Defined at namespace scope in their respective header files (e.g., `using mv = MotionVector` in `include/MotionVector.h:157`)
- `lt` alias defined in `include/RigidBodyInertia.h:19`

**Header Guards:**
- `#ifndef`/`#define`/`#endif` pattern (traditional):
```cpp
#ifndef SPATIAL_VECTOR_H
#define SPATIAL_VECTOR_H
// ...
#endif // SPATIAL_VECTOR_H
```
- Guard names match filename, uppercase with underscores: `SPATIAL_VECTOR_H`, `PLUCKER_TRANSFORM_H`, `MOTION_VECTOR_H`, `RIGID_BODY_INERTIA_H`, `ARTIC_BODY_INERTIA_H`
- Exception: `include/LowerTriangular.h` uses `#pragma once` (line 1)

## Code Style

**Indentation:** 4 spaces (no tabs observed)

**Braces:**
- K&R style for functions (opening brace on same line):
```cpp
void SpatialVector::print() const
{
    std::cout << "Angular: " << angular.transpose() << std::endl;
    std::cout << "Linear: " << linear.transpose() << std::endl;
}
```
- Allman style for class/struct/namespace declarations (opening brace on next line):
```cpp
namespace SpatialAlgebra
{
    class SpatialVector
    {
    protected:
        Vector3d angular;
    public:
        SpatialVector();
    };
}
```
- Space before opening paren in control flow statements:
```cpp
if (n != other.n)
    throw std::invalid_argument("Matrix size mismatch");
```
- No spaces inside parentheses for function calls: `sum.getAngular()[0]`

**Line length:** ~80-120 characters typical, with some longer lines for complex expressions

**Horizontal spacing:**
- Spaces around binary operators: `mass * scalar`, `com.cross(v)`, `i * (i + 1) / 2 + j`
- No space after unary operators: `*this`, `&other`
- Space after comma in parameter lists: `const Vector3d &angular, const Vector3d &linear`

**Linting:**
- No linter configured (no `.clang-tidy`, no `.eslintrc`, no `biome.json`)
- No CI pipeline for lint enforcement

**Formatting:**
- No formatter configured (no `.clang-format`)

## Import Organization

**Include order within `.cpp` files:**
1. Corresponding header first: `#include "SpatialVector.h"`
2. Other library headers: `#include <Eigen/Dense>`
3. Standard library headers: `#include <iostream>`, `#include <cmath>`, `#include <stdexcept>`
4. GTest headers (in test files): `#include <gtest/gtest.h>`

**Include order within `.h` files:**
1. Standard library: `#include <iostream>`, `#include <vector>`, `#include <iomanip>`
2. Eigen headers: `#include <Eigen/Dense>`, `#include <Eigen/Geometry>`
3. Project headers: `#include "SpatialVector.h"`, `#include "Rotation.h"`

**Test file includes follow this order:**
```cpp
#include "ClassUnderTest.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
// Optional: <Eigen/Geometry>, <cmath>, <stdexcept>
```

**Path Aliases:**
- None configured in CMake (no `-I` include path overrides beyond `include/` directory)
- VSCode config in `.vscode/c_cpp_properties.json` specifies: `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`

## Namespace Usage

- All library code in `namespace SpatialAlgebra`:
```cpp
namespace SpatialAlgebra
{
    // class definitions
}
```
- Test files: `using namespace SpatialAlgebra;` at file scope
- `LowerTriangular.h` has a broader export: `using SpatialAlgebra::LowerTriangular;` at file scope (line 569) and `using SpatialAlgebra::Rotation;` in `Rotation.h:175`
- Test files may also use `using namespace Eigen;`

## Comments and Documentation

**Doxygen Tags Used:**
- `@brief` — Required on every class, method, and member declaration
- `@details` — Detailed explanation of mathematical meaning (almost every declaration)
- `@param` — Parameter descriptions with units where applicable (e.g., `(rad/s)`, `(kg⋅m²)`)
- `@return` — Return value description
- `@note` — Implementation notes, warnings about physics interpretation
- `@warning` — Usage warnings (e.g., type-safety, physical validity)
- `@see` — Cross-references to Featherstone textbook chapters, other classes
- `@code{.cpp}` / `@endcode` — Example usage blocks
- `@file` — On header files
- `@throws` — Exception documentation (in headers)

**Comment Style:**
- Block Doxygen comments (`/** ... */`) on all declarations
- Inline `///<` for member variable documentation:
```cpp
Vector3d angular; ///< Angular component (ω for motion vectors, τ for force vectors)
```
- Section separators in test files:
```cpp
// ============================================================================
// TransformRBI Tests
// ============================================================================
```
- Mathematical notation in comments (Greek letters transliterated):
```cpp
// Result: [(0,0,1); (0,-1,0)]
// Linear = (0,-1,0) + (0,0,0) = (0,-1,0)
```
- References to Featherstone textbook with chapter notation:
```cpp
/**
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 2.
 */
```

## Error Handling

**Exception Types Used:**
- `std::invalid_argument` — Dimension mismatches for matrix operations (`include/LowerTriangular.h:134,199,229,251,311,331,381`)
- `std::out_of_range` — Index bounds checking (debug mode only via `#ifndef NDEBUG`):
```cpp
#ifndef NDEBUG
    if (i >= n || j >= n || i < 0 || j < 0)
        throw std::out_of_range("Index out of bounds");
#endif
```
- `std::runtime_error` — Singular matrix in inversion (`include/LowerTriangular.h:410` in the `inverse()` method documentation)

**Pattern:**
- Exceptions used for programmer errors (invalid arguments, out of bounds)
- No error codes or `std::optional` return types
- No custom exception types
- Minimal input validation — relies on caller correctness for most operations
- `InverseDynamics::computeTorques` documents `@throws std::invalid_argument` for NaN/Inf and size mismatch

**Documentation Philosophy:**
- Exceptions are documented with `@throws` Doxygen tag in header files
- `noexcept` specifier used on simple accessors and utility functions:
```cpp
inline const Vector3d &getData() const noexcept { return data; }
inline int getSize() const noexcept { return n; }
```

## Logging

**Framework:** None — uses raw `std::cout`

**Pattern:**
- `print()` method on every class:
```cpp
void print() const;
```
- Implementations write to `std::cout`:
```cpp
void SpatialVector::print() const
{
    std::cout << "Angular: " << angular.transpose() << std::endl;
    std::cout << "Linear: " << linear.transpose() << std::endl;
}
```
- No structured logging
- No log levels
- No log file output

## Function Design

**Size:**
- Small to medium functions (10-50 lines typical)
- Operators implemented inline for small classes (e.g., `RigidBodyInertia.h` entirely inline)

**Parameters:**
- Pass by const reference for objects: `const SpatialVector &other`, `const Eigen::Vector3d &v`
- Pass by value for primitives: `double scalar`, `int size`
- Named parameters not used (C++ limitation)

**Return Values:**
- Return by value for small objects: `SpatialVector`, `Rotation`, `ForceVector`, `MotionVector`
- Return by const reference for internal data: `const Vector3d &getCom() const`, `const Eigen::VectorXd &getData() const`
- No `std::optional` or `std::variant` usage

**Inline Methods:**
- Trivial getters in headers: `inline double getMass() const { return mass; }`
- Simple operators in headers: `RigidBodyInertia` operators are all inline
- Most methods declared in headers, defined in `.cpp` files for non-trivial classes
- `RigidBodyInertia` and `ArticulatedBodyInertia` are entirely inline in their headers

**Method Chaining:**
- Not used — functions return results, not `*this`

## Module Design

**Exports:**
- All classes in `namespace SpatialAlgebra`
- Type aliases at namespace scope in respective headers
- Free functions in `namespace SpatialAlgebra` (`SpatialUtils.h`, `SpatialOperations.h`)
- No umbrella header (`include/SpatialAlgebra.h` not present)

**Barrel Files:**
- None — users include specific headers by class name

**Header Dependencies:**
- `SpatialVector.h` is the foundational header (included by all others)
- `PluckerTransform.h` has the most includes, including forward declarations of `RigidBodyInertia` and `ArticulatedBodyInertia`
- `RigidBodyInertia.h` and `ArticulatedBodyInertia.h` have circular dependency via `PluckerTransform.h` in the transform methods

## Mathematical Notation Conventions

**Variable Names:**
- Greek letters transliterated: `omega`, `tau` (not `ω`, `τ`)
- Vector components as `Vector3d`, `Vector6d`
- Matrix operations as `Matrix3d`, `MatrixXd`
- Cross product: `angular.cross(other.angular)`
- Dot product: `angular.dot(other.angular) + linear.dot(other.linear)`

**Comments:**
- Formulas preserved in mathematical notation with transliteration:
```
// V1 x V2 = [w1 x w2, w1 x v2 + v1 x w2]
```
- Source code comments include manual calculation verification:
```
// ω1·ω2 = 1*2 + 2*3 + 3*4 = 2 + 6 + 12 = 20
// v1·v2 = 4*5 + 5*6 + 6*7 = 20 + 30 + 42 = 92
// Total = 20 + 92 = 112
```

**Units:**
- Documented in `@param` Doxygen tags:
- Motion vectors: `rad/s`, `m/s`
- Force vectors: `N⋅m`, `N`
- Inertia: `kg⋅m²`, `kg⋅m`, `kg`

---

*Convention analysis: 2026-05-17*
