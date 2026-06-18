# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```sh
# Configure and build
cmake -B build
cmake --build build

# Run all tests
cd build && ctest --output-on-failure

# Run a single test executable
./build/TestSpatialVector
./build/TestPluckerTransform
./build/TestRotation
./build/TestLowerTriangular
./build/TestRigidBodyInertia
./build/TestArticulatedBodyInertia
./build/TestForwardDynamics
./build/TestInverseDynamics
./build/TestDynamicsConsistency
./build/TestSpatialOperations
./build/TestSpatialUtils

# Build and run examples
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics
./build/examples/example_dynamics

# Generate documentation
doxygen Doxyfile   # outputs to docs/html/ and docs/latex/
```

### Eigen version mismatch

Eigen 5.x (Homebrew) changes CMake version compatibility. If `find_package(Eigen3 REQUIRED)` fails, either remove the version pin in `CMakeLists.txt:13` or pass `-DEigen3_DIR=$(brew --prefix eigen)/share/eigen3/cmake` at configure time.

### Boost / Pinocchio configure failure

Homebrew Boost 1.89.0 makes `boost_system` header-only (no `libboost_system.dylib`). If you enable `-DSA_BUILD_PINOCCHIO_BENCHMARKS=ON` and configure fails with "Could NOT find Boost (missing: system)", the fix is already in `CMakeLists.txt`: pinocchio is located via `find_path`/`find_library` instead of `find_package(pinocchio REQUIRED)`, which bypasses `pinocchioConfig.cmake`'s `find_package(Boost REQUIRED COMPONENTS system)` call. Do not revert to `find_package(pinocchio REQUIRED)`.

### Coverage build

```sh
cmake -B build -DENABLE_COVERAGE=ON
cmake --build build
```

## GSD Workflow

Before making file edits, use a GSD entry point:
- `/gsd-quick` — small fixes, doc updates, ad-hoc tasks
- `/gsd-debug` — investigation and bug fixing
- `/gsd-execute-phase` — planned phase work

Do not make direct repo edits outside a GSD workflow unless explicitly asked to bypass it.

## Architecture

All code lives in `namespace SpatialAlgebra`. The library is a C++17 static library (`build/libSpatialAlgebra.a`) with Eigen3 as its only runtime dependency.

**Class hierarchy and dependency order** (matches `include/SpatialAlgebra.h` include order):

```
Eigen::Matrix3d
  └── Rotation                  (3×3 rotation matrix)

SpatialVector                   (6D base: angular + linear Vector3d)
  ├── MotionVector  (twist)     [ω; v]
  └── ForceVector   (wrench)    [τ; f]

LowerTriangular                 (packed lower-triangular, custom storage: idx = i*(i+1)/2 + j)

RigidBodyInertia                (mass, COM, inertia tensor via LowerTriangular)
ArticulatedBodyInertia          (ABA inertia — three LowerTriangular blocks)

PluckerTransform                (Rotation + translation; transforms motions, forces, inertias)

ForwardDynamics / Link          (ABA: outward pass → inward pass → solve qddot)
InverseDynamics                 (RNEA: outward pass → inward pass → compute tau)
```

**Free functions** in `SpatialUtils.h`: `skew()`, `dot()`, `cross()` (four overloads: mv×mv, mv×fv, fv×mv, fv×fv).

**Type aliases** used throughout the codebase:

```cpp
using mv   = MotionVector;
using fv   = ForceVector;
using plux = PluckerTransform;
using rbi  = RigidBodyInertia;
using abi  = ArticulatedBodyInertia;
using lt   = LowerTriangular;
using Vector3d = Eigen::Matrix<double, 3, 1>;  // inside namespace SpatialAlgebra
```

**Python**: `robot_dynamics/rnea.py` is a standalone NumPy RNEA implementation, not integrated with the C++ library.

## Mathematical Conventions

Follows Featherstone (2008) throughout. Key conventions:

- Spatial vectors: `[angular; linear]` — angular component first
- Transform direction: `X` stored per-link transforms **parent → child**; `v_child = X · v_parent`
- Force transforms contravariantly: `f_child = X.transformForce(f)` = X^{-T} · f
- Inward pass uses `X.inverseTransformForce(f_child)` to propagate forces to parent
- ABI transform: `X.invtformABI(Ia)` = X^{-1} · Ia · X^{-T}
- Gravity via base acceleration: ABA sets `c₀ = -g`; RNEA sets `a₀ = -[0; g]`

**Numerical tolerances** (do not tighten without justification):
- `1e-10` — single-step unit tests (direct algebraic ops)
- `1e-8` — round-trip dynamics tests (accumulated FP error across RNEA + ABA)

## Code Conventions

- Header guards: `#ifndef`/`#define` (except `LowerTriangular.h` uses `#pragma once`)
- Every class and method has Doxygen `@brief` / `@details`
- Inline comments use mathematical notation: `[ω; v]`, `X = [R, 0; -R[t]×, R]`
- Naming: PascalCase classes, camelCase methods, `get`-prefixed accessors, 4-space indentation
- Pass objects by `const&`; return by value for small spatial objects
- No linter, no formatter, no CI
