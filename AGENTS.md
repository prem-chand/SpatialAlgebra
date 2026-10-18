# SpatialAlgebra agent instructions

## Build

```sh
cmake -S . -B build
cmake --build build
```

The project requires CMake 3.19+, C++17, Eigen 3.4 through 5, and Google Test. CMake fetches Google Test when it is not installed. The static library is `build/libSpatialAlgebra.a`.

## Test

```sh
ctest --test-dir build --output-on-failure
```

CTest covers the spatial vector, transform, rotation, lower-triangular, utility, inertia, dynamics, CRBA, and test-model smoke-test targets. The MuJoCo adapter smoke test is added when CMake finds the MuJoCo library and headers. The longer comparison benchmark is opt-in with `-DSA_BUILD_MUJOCO_BENCHMARKS=ON`.

## Examples and optional benchmarks

Build examples with `cmake --build build --target example_vectors example_transforms example_inertia example_dynamics example_robot_2link example_robot_3link`.

Google Benchmark targets are disabled by default. Enable them with `-DSA_BUILD_BENCHMARKS=ON`. MuJoCo adapter targets are enabled automatically when MuJoCo is installed. Enable the comparison benchmark with `-DSA_BUILD_MUJOCO_BENCHMARKS=ON`. The model fetch script is `scripts/fetch_models.sh`.

## Architecture

All C++ classes and functions live in `namespace SpatialAlgebra`.

- `SpatialVector` is the six-dimensional base for `MotionVector` and `ForceVector`.
- `Rotation` and `PluckerTransform` represent orientation and spatial frame transforms.
- `RigidBodyInertia` and `ArticulatedBodyInertia` represent body and articulated inertia.
- `ForwardDynamics` implements ABA. `InverseDynamics` implements RNEA.
- `CompositeRigidBodyInertia.h` declares CRBA for computing the joint-space mass matrix. Include it separately from `SpatialAlgebra.h`.
- `LowerTriangular` stores a symmetric lower-triangular matrix in packed form.
- `SpatialUtils.h` provides free functions including `skew()`, `dot()`, and `cross()`.

The standalone Python implementation in `robot_dynamics/rnea.py` is not connected to the C++ library.

## Project layout

- `include/` contains public C++ headers.
- `src/` contains C++ implementations.
- `examples/` contains runnable programs.
- `tests/` contains unit and integration tests.
- `tests/test-models/` contains reusable test models and adapters.
- `benchmarks/` contains optional Google Benchmark programs.
- `models/menagerie/` contains optional MuJoCo model assets.
- `docs/` contains generated Doxygen output.

## Conventions

- Use C++17 and preserve existing public interfaces unless the task requests an API change.
- Document public declarations with Doxygen `@brief` and `@details` block comments.
- Use `#ifndef`/`#define` include guards. `LowerTriangular.h` retains its existing `#pragma once` guard.
- Follow the existing four-space indentation and naming conventions.
- Do not edit generated build output under `build/` or generated Doxygen output under `docs/` by hand.

## Project skill

Use [`.claude/skills/verify/SKILL.md`](.claude/skills/verify/SKILL.md) to build and drive the examples while capturing verification evidence.
