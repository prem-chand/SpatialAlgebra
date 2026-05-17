# SpatialAlgebra

A C++17 library implementing spatial vector algebra for rigid body dynamics, following Featherstone's formulation. Provides 6D spatial vectors (twists and wrenches), Plücker coordinate transforms, and inertia representations for robotics simulation and control.

## Features

- **Spatial Vectors**: 6D vectors combining angular and linear components
  - `SpatialVector` — Base class for 6D spatial vectors
  - `MotionVector` (twist) — Spatial velocity [ω; v]
  - `ForceVector` (wrench) — Spatial force [τ; f]

- **Transforms**: Coordinate frame transformations in Plücker coordinates
  - `Rotation` — 3×3 rotation matrices (extends Eigen::Matrix3d)
  - `PluckerTransform` — 6×6 spatial transforms (rotation + translation)

- **Inertia Representations**: Mass property representations
  - `RigidBodyInertia` — Mass, center of mass, inertia tensor
  - `ArticulatedBodyInertia` — Articulated body inertia for forward dynamics

- **Dynamics Algorithms**:
  - `ForwardDynamics` — Articulated Body Algorithm (ABA) for computing accelerations
  - `InverseDynamics` — Recursive Newton-Euler Algorithm (RNEA) for computing torques
  - **Gravity Support** — Optional gravity parameter on `computeAccelerations()` and `computeTorques()` using Featherstone base acceleration formulation (`a₀ = -g`). Defaults to zero for backward compatibility.

- **Utilities**:
  - `LowerTriangular` — Packed-storage lower-triangular matrix
  - `SpatialOperations` — Static utility functions
  - `SpatialUtils` — Free functions: `skew()`, `dot()`, `cross()`

- **Single Header**: Include `SpatialAlgebra.h` to bring in all library headers at once.

## Requirements

- **CMake** 3.10+
- **C++17** compatible compiler (g++, clang++)
- **Eigen3** 3.3+ (`brew install eigen`)
- **Google Test** (`brew install googletest`)

## Build Instructions

```sh
# Configure
cmake -B build

# Build library and tests
cmake --build build

# Run tests
cd build && ctest --output-on-failure
```

The library builds to `build/libSpatialAlgebra.a`.

## Usage Examples

### Creating Motion and Force Vectors

```cpp
#include "MotionVector.h"
#include "ForceVector.h"
#include "SpatialUtils.h"
#include <iostream>

using namespace SpatialAlgebra;

int main() {
    // Create a motion vector (twist) with angular velocity [1, 0, 0] and linear velocity [0, 1, 0]
    MotionVector twist(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    
    // Create a force vector (wrench) with torque [0, 0, 5] and force [10, 0, 0]
    ForceVector wrench(Vector3d(0, 0, 5), Vector3d(10, 0, 0));
    
    // Vector operations
    MotionVector twist2 = twist * 2.0;
    ForceVector total = wrench + wrench;
    
    // Cross products (free functions from SpatialUtils.h)
    ForceVector result = cross(twist, wrench);   // Motion × Force → Force
    MotionVector motion_cross = cross(twist, twist);  // Motion × Motion → Motion
    ForceVector force_cross = cross(wrench, wrench);  // Force × Force → Force
    
    twist.print();
    wrench.print();
    
    return 0;
}
```

### Creating and Using Plücker Transforms

```cpp
#include "PluckerTransform.h"
#include "Rotation.h"
#include <iostream>

using namespace SpatialAlgebra;

int main() {
    // Create a rotation of 90 degrees around Z axis
    Rotation R;
    R.setFromAngleAxis(M_PI_2, Vector3d(0, 0, 1));
    
    // Create a translation of [1, 0, 0]
    Vector3d translation(1, 0, 0);
    
    // Create Plücker transform
    PluckerTransform X(R, translation);
    
    // Create a motion vector in the source frame
    MotionVector v_source(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    
    // Transform motion to target frame
    MotionVector v_target = X.transformMotion(v_source);
    
    // Transform force to target frame
    ForceVector f_source(Vector3d(0, 0, 0), Vector3d(1, 0, 0));
    ForceVector f_target = X.transformForce(f_source);
    
    // Compute inverse transform
    PluckerTransform X_inv = X.inverse();
    
    v_target.print();
    f_target.print();
    
    return 0;
}
```

### Working with Rigid Body Inertia

```cpp
#include "RigidBodyInertia.h"
#include "ArticulatedBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include <iostream>

using namespace SpatialAlgebra;

int main() {
    // Create rigid body inertia: mass=2kg, COM at [0, 0, 0.5], inertia tensor
    Matrix3d I_body;
    I_body << 0.1, 0, 0,
              0, 0.1, 0,
              0, 0, 0.2;
    
    RigidBodyInertia I(2.0, Vector3d(0, 0, 0.5), I_body);
    
    // Apply inertia to motion to get force
    MotionVector v(Vector3d(0, 0, 1), Vector3d(0, 0, 0));
    ForceVector f = I.apply(v);
    
    // Create articulated body inertia from rigid body inertia
    ArticulatedBodyInertia Ia(I);
    
    // Transform inertia by Plücker transform
    Rotation R;
    R.setIdentity();
    PluckerTransform X(R, Vector3d(0, 0, 0.1));
    RigidBodyInertia I_transformed = X.tformRBI(I);
    
    I.print();
    f.print();
    
    return 0;
}
```

### Forward Dynamics with Articulated Body Algorithm

```cpp
#include "ForwardDynamics.h"
#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include <iostream>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
using lt = LowerTriangular;

int main() {
    // Create a 2-link chain using Link structures
    ForwardDynamics fd;
    
    // Link 0 (base): identity transform, Z-axis revolute joint
    Link link0;
    link0.parent = -1;
    link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link0.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link0.q = 0.0;
    link0.qdot = 0.0;
    fd.links.push_back(link0);
    
    // Link 1 (child of link 0, offset along X)
    Link link1;
    link1.parent = 0;
    link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link1.q = 0.0;
    link1.qdot = 0.0;
    fd.links.push_back(link1);
    
    // Joint torques
    Eigen::VectorXd tau(2);
    tau[0] = 1.0;
    tau[1] = 0.5;
    
    // Compute joint accelerations using ABA (with optional gravity)
    fd.computeAccelerations(tau);
    // With gravity: fd.computeAccelerations(tau, Vector3d(0, 0, -9.81));
    
    // Results stored in link.qddot
    std::cout << "Joint accelerations:" << std::endl;
    std::cout << "  Joint 1: " << fd.links[0].qddot << " rad/s²" << std::endl;
    std::cout << "  Joint 2: " << fd.links[1].qddot << " rad/s²" << std::endl;
    
    return 0;
}
```

## Running Examples

The `examples/` directory contains compilable examples demonstrating core operations:

```sh
# Build all examples
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics

# Run individual examples
./build/examples/example_vectors      # Basic vector operations
./build/examples/example_transforms   # Plücker coordinate transforms
./build/examples/example_inertia      # Rigid body inertia
./build/examples/example_dynamics     # Forward dynamics (ABA)
```

## Testing

```sh
# Build and run all tests
cmake --build build
cd build && ctest --output-on-failure

# Run individual test executables
./build/TestSpatialVector
./build/TestPluckerTransform
./build/TestRotation
./build/TestLowerTriangular
```

Test coverage:
- `TestSpatialVector` — Spatial vector operations (MotionVector, ForceVector)
- `TestPluckerTransform` — Plücker transform operations
- `TestRotation` — Rotation matrix operations
- `TestLowerTriangular` — Packed lower-triangular matrix operations

## Documentation

Generate browsable API documentation with Doxygen:

```sh
doxygen Doxyfile
```

Documentation is generated in:
- `docs/html/` — HTML documentation (open `docs/html/index.html` in a browser)
- `docs/latex/` — LaTeX documentation

## Project Structure

```
SpatialAlgebra/
├── include/              # Header files
│   ├── SpatialAlgebra.h  # Umbrella header (includes all public headers)
│   ├── SpatialVector.h
│   ├── MotionVector.h
│   ├── ForceVector.h
│   ├── Rotation.h
│   ├── PluckerTransform.h
│   ├── RigidBodyInertia.h
│   ├── ArticulatedBodyInertia.h
│   ├── ForwardDynamics.h
│   ├── InverseDynamics.h
│   ├── LowerTriangular.h
│   ├── SpatialOperations.h
│   └── SpatialUtils.h
├── src/                  # Implementation files
├── tests/                # Test files
├── examples/             # Usage examples
├── CMakeLists.txt        # Build configuration
├── Doxyfile              # Doxygen configuration
└── README.md             # This file
```

## Type Aliases

For concise notation, the library uses these type aliases:

```cpp
using mv = MotionVector;           // Motion vector (twist)
using fv = ForceVector;            // Force vector (wrench)
using plux = PluckerTransform;     // Plücker transform
using rbi = RigidBodyInertia;      // Rigid body inertia
using abi = ArticulatedBodyInertia; // Articulated body inertia
using lt = LowerTriangular;        // Lower triangular matrix
```

## Mathematical Notation

The library follows Featherstone's spatial vector algebra notation:

- Motion vectors (twists): `[ω; v]` where ω is angular velocity, v is linear velocity
- Force vectors (wrenches): `[τ; f]` where τ is torque, f is force
- Plücker transform: `X = [R, 0; -R[t]×, R]` where R is rotation, t is translation
- Cross products: `×` denotes spatial cross product (different from 3D cross product)

## API Changes in v1.1

The following changes were introduced in version 1.1:

- **Gravity parameter**: `ForwardDynamics::computeAccelerations()` and `InverseDynamics::computeTorques()` now accept an optional `Vector3d gravity` parameter (default `Vector3d::Zero()`) implementing Featherstone's base acceleration formulation (`a₀ = -g`). This enables gravity-aware dynamics without breaking existing code that omits the parameter.

- **Removed cross product overloads**: `MotionVector::crossForce()` and `ForceVector::crossMotion()` member functions have been removed. Use the free functions from `SpatialUtils.h` instead:
  - `cross(MotionVector, ForceVector)` → `ForceVector`
  - `cross(MotionVector, MotionVector)` → `MotionVector`
  - `cross(ForceVector, ForceVector)` → `ForceVector`

- **Namespace change**: `Vector3d` moved from global scope into the `SpatialAlgebra` namespace. Existing code using `Vector3d` without namespace qualification may need `using SpatialAlgebra::Vector3d;` or `using namespace SpatialAlgebra;`.

- **Umbrella header**: New `#include "SpatialAlgebra.h"` includes all public library headers in dependency order, providing a single include point for convenience.

- **Empty stub files removed**: Source stubs (`src/RigidBodyInertia.cpp`, `src/ArticulatedBodyInertia.cpp`, `src/SpatialOperations.cpp`) have been removed. These classes are now fully header-inline or implemented in existing source files.

- **OpenMP dependency removed**: `LowerTriangular` operations no longer require OpenMP. The library builds without special compiler flags for parallel execution.

- **GTest FetchContent**: CMakeLists.txt includes a FetchContent fallback for Google Test when it is not installed system-wide, improving CI and cross-platform build compatibility.

## License

[License information to be added]

## References

- Featherstone, R. (2008). *Rigid Body Dynamics Algorithms*. Springer.
- Featherstone, R. (1987). *Robot Dynamics Algorithms*. Kluwer Academic Publishers.
