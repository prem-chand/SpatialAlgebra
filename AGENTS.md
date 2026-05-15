# SpatialAlgebra — Agent Instructions

## Build

```sh
cmake -B build
cmake --build build
```

C++17, requires **Eigen3** (`brew install eigen`) and **Google Test** (`brew install googletest`). The library builds to `build/libSpatialAlgebra.a`.

### Eigen version mismatch

Eigen 5.0.1 installed via Homebrew may fail CMake's `find_package(Eigen3 3.3 REQUIRED)` because Eigen 5.x changed its CMake version-compatibility range. Either:

- Remove the version pin: change `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` → `find_package(Eigen3 REQUIRED NO_MODULE)` in `CMakeLists.txt:13`
- Or set `-DEigen3_DIR=$(brew --prefix eigen)/share/eigen3/cmake`

## Run tests

```sh
cmake --build build && cd build && ctest --output-on-failure
# or run individual test executables directly:
build/TestSpatialVector
build/TestPluckerTransform
```

Only 2 test executables are registered in CMakeLists.txt (`TestSpatialVector`, `TestPluckerTransform`). The other 4 test files in `tests/` are empty stubs and are **not** compiled or linked.

Test style is mixed: `TestSpatialVector.cpp` uses bare `assert()`, `TestPluckerTransform.cpp` uses GTest (`TEST()` / `EXPECT_DOUBLE_EQ`).

## Architecture

Everything lives in `namespace SpatialAlgebra`. Core class hierarchy:

- `SpatialVector` (6D base: angular + linear components) → `MotionVector` (twist) / `ForceVector` (wrench)
- `Rotation` — 3×3 matrix extending `Eigen::Matrix3d`
- `PluckerTransform` — 6×6 Plücker coordinate transform (rotation + translation)
- `RigidBodyInertia` / `ArticulatedBodyInertia` — inertia representations
- `LowerTriangular` — packed-storage lower-triangular matrix (custom, not Eigen)
- `SpatialOperations` — static utility class

**Type aliases** (used throughout code): `mv` = MotionVector, `fv` = ForceVector, `plux` = PluckerTransform, `rbi` = RigidBodyInertia, `abi` = ArticulatedBodyInertia, `lt` = LowerTriangular.

**Utility functions** in `SpatialUtils.h`: `skew()`, `dot()`, `cross()` (free functions, not methods).

**Python side**: `robot_dynamics/rnea.py` has a standalone RNEA implementation using NumPy. It is *not* integrated with the C++ library.

## Code conventions

- Doxygen `@brief` / `@details` on every declaration (block comments, not inline)
- Include guards are `#ifndef`/`#define` (except `LowerTriangular.h` uses `#pragma once`)
- Typedefs at file level: `using Vector3d = Eigen::Matrix<double, 3, 1>`
- Codegen: `doxygen Doxyfile` generates HTML docs into `docs/html/`
- No CI workflows, no linter, no formatter, no type checker

<!-- GSD:project-start source:PROJECT.md -->
## Project

**SpatialAlgebra Project**

A C++17 library implementing spatial vector algebra for rigid body dynamics, following Featherstone's formulation. Provides 6D spatial vectors (twists and wrenches), Plücker coordinate transforms, and inertia representations for robotics simulation and control.

**Core Value:** **Must Deliver:** Complete, well-tested spatial algebra library where all core classes are fully implemented and verified with comprehensive tests.

**Success Looks Like:** 
- All incomplete implementations finished (ArticulatedBodyInertia, RigidBodyInertia, SpatialOperations)
- Forward dynamics (Articulated Body Algorithm - ABA) implemented
- 100% test coverage for all mathematical operations
- Library ready for integration into robotics projects

### Constraints

**Technical:**
- Must maintain Eigen3 compatibility (version 3.3+ or 5.x with testing)
- Preserve existing class interfaces (backward compatibility)
- Follow existing code conventions (Doxygen comments, type aliases)

**Timeline:** Flexible — quality-focused completion

**Budget:** N/A (open source library)
<!-- GSD:project-end -->

<!-- GSD:stack-start source:codebase/STACK.md -->
## Technology Stack

## Languages
- C++17 - Core library implementation in `include/` and `src/`
- Python 3 - Standalone RNEA implementation in `robot_dynamics/rnea.py`
## Runtime
- Native compiled C++ library (no runtime)
- macOS (Darwin) development environment
- None (CMake-based build)
- Homebrew for system dependencies (Eigen3, Google Test)
## Frameworks
- Eigen3 3.3+ - Linear algebra backend for all matrix/vector operations
- Google Test (GTest) - Test framework for `TestPluckerTransform.cpp`
- Basic `assert()` - Minimal testing in `TestSpatialVector.cpp`
- CMake 3.10+ - Build system configuration
- g++ (GCC) - C++ compiler (configured in `CMakeLists.txt:7`)
## Key Dependencies
- Eigen3 3.3+ - Required for all linear algebra operations
- Google Test - Unit testing framework
## Configuration
- No runtime environment variables required
- Pure library with compile-time configuration
- `CMakeLists.txt` - Main build configuration
- `.vscode/c_cpp_properties.json` - VSCode IntelliSense configuration
- `Doxyfile` - Doxygen documentation generation
## Platform Requirements
- CMake 3.10+
- C++17 compatible compiler (g++, clang++)
- Eigen3 3.3+ (`brew install eigen`)
- Google Test (`brew install googletest`)
- Compiled static library: `build/libSpatialAlgebra.a`
- Header-only style with separate implementation files
- No runtime dependencies beyond standard C++ library
## Build Artifacts
- `build/libSpatialAlgebra.a` - Static library
- `build/TestSpatialVector` - Spatial vector tests
- `build/TestPluckerTransform` - Plücker transform tests
- `docs/html/` - Generated Doxygen HTML documentation
- `docs/latex/` - Generated Doxygen LaTeX documentation
<!-- GSD:stack-end -->

<!-- GSD:conventions-start source:CONVENTIONS.md -->
## Conventions

## Naming Patterns
- PascalCase with descriptive names: `SpatialVector.h`, `PluckerTransform.cpp`
- Test files prefixed with `Test`: `TestSpatialVector.cpp`
- Header guards match filename: `SPATIAL_VECTOR_H`, `PLUCKER_TRANSFORM_H`
- PascalCase: `SpatialVector`, `MotionVector`, `ForceVector`, `PluckerTransform`
- Descriptive of mathematical concept being modeled
- camelCase: `transformMotion`, `crossProductMotion`, `getAngular`, `setFromAngleAxis`
- Verbs for operations: `transform`, `apply`, `inverse`, `transpose`
- Getters prefixed with `get`: `getMass`, `getCom`, `getInertia`
- Private members: lowercase (e.g., `angular`, `linear`, `rotation`, `translation`, `mass`, `com`)
- Local variables: camelCase (e.g., `transformedAngular`, `newRotation`, `invRotation`)
- Eigen types: `Vector3d`, `Matrix3d`, `VectorXd` (using Eigen's typedefs)
- Type aliases at file scope: `using Vector3d = Eigen::Matrix<double, 3, 1>`
- Namespace aliases: `using mv = MotionVector`, `using fv = ForceVector`
- Lowercase abbreviations for concise notation: `mv`, `fv`, `plux`, `rbi`, `abi`, `lt`
## Code Style
- Indentation: 4 spaces (no tabs observed)
- Line length: ~100-120 characters typical
- Braces: K&R style for functions, Allman style for classes
- No linter configured
- No formatter configured (no .clang-format, .prettierrc, etc.)
## Import Organization
- Traditional `#ifndef`/`#define`/`#endif` pattern
- Example: `#ifndef SPATIAL_VECTOR_H` → `#define SPATIAL_VECTOR_H`
- Exception: `LowerTriangular.h` uses `#pragma once`
#include <Eigen/Dense>
#include "Dependency.h"
#include "ClassName.h"      // Corresponding header first
#include "OtherDependency.h"
#include <Eigen/Dense>
#include <iostream>
- None configured in CMake
- VSCode config specifies: `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`
## Error Handling
- Exceptions for runtime errors:
- Debug-mode assertions via `#ifndef NDEBUG` guards
- No error codes or `std::optional` return types
- No custom exception types
- Minimal input validation
- Relies on caller correctness for most operations
- Bounds checking only in debug mode
## Logging
- `print()` method on every class:
- No structured logging
- No log levels
- No log file output
## Comments
- Every class and method has Doxygen documentation
- Inline comments explain mathematical formulas
- Complex operations have detailed `@details` sections
- Block comments on all declarations:
- Mathematical notation in comments (e.g., `[ω; v]`, `[τ; f]`)
- References to literature (Featherstone textbook)
## Function Design
- Small to medium functions (10-50 lines typical)
- Single responsibility per function
- Inline methods for trivial getters in headers
- Pass by const reference for objects: `const SpatialVector &other`
- Pass by value for primitives: `double scalar`, `int size`
- Named parameters not used (C++ limitation)
- Return by value for small objects (SpatialVector, Rotation)
- Return by const reference for internal data: `const Vector3d &getCom() const`
- No `std::optional` or `std::variant` usage
## Module Design
- All classes in `namespace SpatialAlgebra`
- Type aliases at namespace scope: `using mv = MotionVector`
- Free functions in `namespace SpatialAlgebra` (SpatialUtils.h)
- None (users include specific headers)
- No `include/SpatialAlgebra.h` umbrella header
- Most methods declared in headers, defined in `.cpp` files
- Inline methods in headers for trivial operations (getters, operators)
- `RigidBodyInertia` entirely inline in header
## Documentation
- Config: `Doxyfile`
- Command: `doxygen Doxyfile`
- Output: `docs/html/`, `docs/latex/`
- Tags used: `@brief`, `@details`, `@param`, `@return`, `@note`, `@warning`, `@see`, `@code`
## Mathematical Notation
- Greek letters transliterated: `omega`, `tau` (not ω, τ)
- Vector notation: `Vector3d`, `Vector6d`
- Matrix notation: `Matrix3d`, `MatrixXd`
- Mathematical notation preserved: `[ω; v]`, `[τ; f]`
- Formulas in code comments: `X = [R, 0; -R[t]x, R]`
- References to Featherstone textbook
<!-- GSD:conventions-end -->

<!-- GSD:architecture-start source:ARCHITECTURE.md -->
## Architecture

## Pattern Overview
- CRTP-like inheritance for spatial vectors (SpatialVector → MotionVector/ForceVector)
- Composition over inheritance for transforms (PluckerTransform contains Rotation + translation)
- Value semantics (objects passed by value/const reference)
- Eigen3 integration via inheritance (Rotation extends Eigen::Matrix3d)
## Class Hierarchy
```
```
```
```
```
```
## Layers
- Purpose: Fundamental spatial vector operations
- Location: `include/SpatialVector.h`, `src/SpatialVector.cpp`
- Contains: SpatialVector base class with angular/linear components
- Depends on: Eigen3
- Used by: All higher-level classes
- Purpose: Type-safe motion (twist) and force (wrench) representations
- Location: `include/MotionVector.h`, `include/ForceVector.h`
- Contains: Specialized spatial vectors with physical interpretation
- Depends on: SpatialVector
- Used by: PluckerTransform, inertia classes
- Purpose: Coordinate frame transformations in Plücker coordinates
- Location: `include/PluckerTransform.h`, `src/PluckerTransform.cpp`
- Contains: 6x6 spatial transforms (rotation + translation)
- Depends on: Rotation, SpatialVector, RigidBodyInertia, ArticulatedBodyInertia
- Used by: Dynamics algorithms
- Purpose: Mass property representations
- Location: `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`
- Contains: Mass, COM, inertia tensor (packed storage)
- Depends on: LowerTriangular, SpatialVector
- Used by: Dynamics algorithms
- Purpose: Helper functions and specialized data structures
- Location: `include/SpatialUtils.h`, `include/LowerTriangular.h`, `include/Rotation.h`
- Contains: skew(), dot(), cross(), packed matrix storage
- Depends on: Eigen3
- Used by: All layers
## Data Flow
## Key Abstractions
- Purpose: 6D vector combining angular + linear components
- Examples: `include/SpatialVector.h:68`
- Pattern: Base class with protected angular/linear Vector3d members
- Purpose: Rigid body coordinate transformation
- Examples: `include/PluckerTransform.h:79`
- Pattern: Stores rotation and translation separately, provides motion/force transform methods
- Purpose: Memory-efficient packed storage for symmetric matrices
- Examples: `include/LowerTriangular.h:70`
- Pattern: 1D array storage with index mapping: `idx = i*(i+1)/2 + j`
## Entry Points
- Location: `include/*.h` (header files)
- Triggers: User includes headers and instantiates classes
- Responsibilities: Provide spatial algebra operations
- Location: `tests/TestSpatialVector.cpp`, `tests/TestPluckerTransform.cpp`
- Triggers: Manual execution via `build/TestSpatialVector`
- Responsibilities: Verify core functionality
- Location: `src/main.cpp`
- Triggers: `cmake --build build` produces executable
- Responsibilities: Demonstrate library usage
## Error Handling
- `std::invalid_argument` for dimension mismatches (LowerTriangular operations)
- `std::out_of_range` for index bounds (debug mode only)
- `assert()` for basic validation (TestSpatialVector.cpp)
- No error codes or result types
## Cross-Cutting Concerns
<!-- GSD:architecture-end -->

<!-- GSD:skills-start source:skills/ -->
## Project Skills

No project skills found. Add skills to any of: `.claude/skills/`, `.agents/skills/`, `.cursor/skills/`, `.github/skills/`, or `.codex/skills/` with a `SKILL.md` index file.
<!-- GSD:skills-end -->

<!-- GSD:workflow-start source:GSD defaults -->
## GSD Workflow Enforcement

Before using Edit, Write, or other file-changing tools, start work through a GSD command so planning artifacts and execution context stay in sync.

Use these entry points:
- `/gsd-quick` for small fixes, doc updates, and ad-hoc tasks
- `/gsd-debug` for investigation and bug fixing
- `/gsd-execute-phase` for planned phase work

Do not make direct repo edits outside a GSD workflow unless the user explicitly asks to bypass it.
<!-- GSD:workflow-end -->

<!-- GSD:profile-start -->
## Developer Profile

> Profile not yet configured. Run `/gsd-profile-user` to generate your developer profile.
> This section is managed by `generate-claude-profile` -- do not edit manually.
<!-- GSD:profile-end -->
