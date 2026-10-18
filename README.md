# SpatialAlgebra

SpatialAlgebra is a C++17 library for spatial vector algebra and rigid-body dynamics. It follows Featherstone's formulation and uses Eigen for its vector and matrix types.

## Features

- Six-dimensional motion vectors and force vectors, with spatial cross and dot products.
- Rotations and Plucker transforms for changing coordinate frames.
- Rigid-body and articulated-body inertia representations.
- Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics.
- Articulated Body Algorithm (ABA) for forward dynamics.
- Joint-position-dependent kinematics for revolute and prismatic joints.
- Composite Rigid Body Algorithm (CRBA) for the joint-space mass matrix. CRBA is declared in `CompositeRigidBodyInertia.h` and is not included by the `SpatialAlgebra.h` umbrella header.
- Optional Google Benchmark and MuJoCo comparison targets.

## Requirements

- CMake 3.19 or newer.
- A C++17 compiler.
- Eigen 3.4 through 5.
- Google Test. CMake uses FetchContent when it is not installed.

MuJoCo is optional. CMake enables the MuJoCo adapter and comparison programs when it finds the MuJoCo C library and headers.

## Build and test

Configure and build the library, examples, and tests:

```sh
cmake -S . -B build
cmake --build build
```

Run the registered CTest tests:

```sh
ctest --test-dir build --output-on-failure
```

When MuJoCo is available, CMake adds the adapter smoke test. The longer comparison benchmark is excluded from the default build and CTest run.

The static library is written to `build/libSpatialAlgebra.a`.

## Examples

Build all example programs:

```sh
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics example_robot_2link example_robot_3link
```

Run a two-link forward-dynamics example:

```sh
./build/examples/example_dynamics
```

The `examples/` directory also contains vector, transform, inertia, and two-link robot examples. Each program uses fixed sample values and prints its results to standard output.

## Basic use

Include the umbrella header for the core algebra types:

```cpp
#include "SpatialAlgebra.h"

using namespace SpatialAlgebra;

MotionVector twist(Vector3d(0, 0, 1), Vector3d::Zero());
ForceVector wrench(Vector3d(0, 0, 2), Vector3d(3, 0, 0));
double power = dot(twist, wrench);
```

Build a solver by adding `Link` values to `ForwardDynamics::links` or `InverseDynamics::links`. Set each link's reference transform `X0`, joint type, axis, and position. `updateKinematics()` derives the current transform from the joint position. The forward and inverse solvers update kinematics when computing accelerations or torques.

Call `computeCompositeRigidBodyInertia()` from `CompositeRigidBodyInertia.h` after updating the link kinematics to calculate the joint-space mass matrix. The current CRBA implementation documents serial chains as its supported and tested topology.

## Benchmarks and MuJoCo comparison

Google Benchmark programs are optional and disabled by default. Enable them with:

```sh
cmake -S . -B build -DSA_BUILD_BENCHMARKS=ON
cmake --build build --target bench_all
./build/benchmarks/bench_all
```

The optional MuJoCo comparison uses MJCF robot models. Fetch the supported Menagerie models before running the comparison:

```sh
./scripts/fetch_models.sh
cmake -S . -B build -DSA_BUILD_MUJOCO_BENCHMARKS=ON
cmake --build build --target test_mujoco_adapter bench_mujoco_comparison
./build/tests/test-models/bench_mujoco_comparison
```

The fetch script downloads model assets into `models/menagerie/`. Keep those third-party model files with their upstream license notices.

## Repository layout

```text
include/                 Public C++ headers
src/                     Library implementation
examples/                Runnable usage examples
tests/                   Unit, dynamics, and integration tests
tests/test-models/       Test model types and solver adapters
benchmarks/              Optional Google Benchmark programs
models/menagerie/        Optional MuJoCo model assets
robot_dynamics/          Standalone Python RNEA implementation
scripts/                 Development and model-fetch scripts
docs/                    Generated Doxygen output
```

## Documentation

Generate the Doxygen API reference with:

```sh
doxygen Doxyfile
```

The HTML output is written to `docs/html/`.

## License

SpatialAlgebra is distributed under the GNU General Public License version 3. See [LICENSE](LICENSE).

## References

- Featherstone, R. (2008). *Rigid Body Dynamics Algorithms*. Springer.
- Featherstone, R. (1987). *Robot Dynamics Algorithms*. Kluwer Academic Publishers.
