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
./build/TestCompositeRigidBodyInertia

# Build and run examples
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics example_robot_2link example_robot_3link
./build/examples/example_dynamics

# Generate documentation
doxygen Doxyfile   # outputs to docs/html/ and docs/latex/
```

### Coverage build

```sh
cmake -B build -DENABLE_COVERAGE=ON
cmake --build build
```

## Performance benchmarks (optional)

Google Benchmark is fetched via `FetchContent`. Opt in with `-DSA_BUILD_BENCHMARKS=ON`:

```sh
cmake -B build -DSA_BUILD_BENCHMARKS=ON
cmake --build build --target bench_all
./build/benchmarks/bench_all
```

Benchmark sources live under `benchmarks/{aba,rnea,core}/` with shared model-generation and randomized-state helpers in `benchmarks/common/`. Compiled with `-O3 -DNDEBUG -flto` regardless of the top-level build type.

## MuJoCo cross-validation (optional)

Requires the MuJoCo C library (3.4+): `brew install mujoco`. Build with:

```sh
cmake -B build
cmake --build build
```

If CMake finds MuJoCo, it builds `test_mujoco_adapter` and adds it to CTest. This test compares SpatialAlgebra with MuJoCo on a rotated chain and on Go1. It also compares each parsed Menagerie model with MuJoCo's native model of the same MJCF. The `bench_mujoco_comparison` accuracy and timing benchmark is opt-in with `-DSA_BUILD_MUJOCO_BENCHMARKS=ON`. It is not in CTest.

### Supported robot models

The benchmark suite uses MJCF models from the [MuJoCo Menagerie](https://github.com/google-deepmind/mujoco_menagerie):
- **Unitree Go1** — 12 DOF quadruped
- **Unitree Go2** — 12 DOF quadruped  
- **Unitree H1** — 19 DOF humanoid
- **Unitree G1** — 29 DOF humanoid

Run `scripts/fetch_models.sh` to download them before running benchmarks.

### Floating-base handling

Floating-base robots (all models above) have a 6-DOF free joint between the world and the base body. The parser converts this into 6 sequential 1-DOF joints (prismatic X/Y/Z, then revolute X/Y/Z) that correctly model the floating base dynamics. The base body's mass, COM, and inertia are assigned to the last virtual joint (revolute Z).

The SA library's RNEA and ABA algorithms handle these virtual joints natively — no special floating-base code path is needed. The dynamics are mathematically correct: gravity propagates from the first link, Coriolis terms include the rotational-translational coupling, and the mass matrix includes the full 6×6 floating base inertia block.

Both adapters in the comparison use this same serial representation. So SpatialAlgebra and MuJoCo agree to about 1e-9 on torques, accelerations, gravity torques, and the mass matrix. The parser drops armature, damping, friction loss, and joint limits.

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
CompositeRigidBodyInertia       (CRBA: O(n²) mass matrix H(q); operates on ForwardDynamics::Link)
```

`CompositeRigidBodyInertia.h` is **not** part of the `SpatialAlgebra.h` umbrella include — include it separately.

**Joint kinematics**: `Link` carries `q`/`qdot`/`qddot`, a `JointType` (`FIXED`/`REVOLUTE`/`PRISMATIC`), `jointAxis`, and a reference transform `X0` (pose at `q=0`). `updateKinematics()` on `ForwardDynamics`/`InverseDynamics` recomputes each link's `X` from `X0` and the current `q` (`X(q) = X0 * X_joint(q)`); it runs automatically at the start of `computeAccelerations()`/`computeTorques()`, and must be called manually before `computeCompositeRigidBodyInertia()`.

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

**`tests/test-models/`** is a separate, zero-SpatialAlgebra-dependency layer used to validate the library against ground truth without coupling test model definitions to `SpatialAlgebra` types:
- `test_models` — INTERFACE target, header-only, links only `Eigen3::Eigen`. Model/robot definitions (`robot_model.h`, `robot_solver.h`) live here.
- `sa_test_adapter` (`sa_adapter.cpp`) — the sole bridge from `test_models` to `SpatialAlgebra`.
- `mujoco_test_adapter` / `mjcf_parser` — optional, built only when the MuJoCo C library is found; parse MJCF into `test_models` types for cross-validation.
- `compile_smoke_test` links only `test_models` and fails to build if any `SpatialAlgebra` type leaks into the model headers — this is the zero-dependency guarantee check.

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
