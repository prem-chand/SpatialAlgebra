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

- **Utilities**:
  - `LowerTriangular` — Packed-storage lower-triangular matrix
  - `SpatialOperations` — Static utility functions
  - `SpatialUtils` — Free functions: `skew()`, `dot()`, `cross()`

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
#include "SpatialAlgebra/MotionVector.h"
#include "SpatialAlgebra/ForceVector.h"
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
    
    // Cross products
    ForceVector result = twist.crossForce(wrench);
    
    twist.print();
    wrench.print();
    
    return 0;
}
```

### Creating and Using Plücker Transforms

```cpp
#include "SpatialAlgebra/PluckerTransform.h"
#include "SpatialAlgebra/Rotation.h"
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
#include "SpatialAlgebra/RigidBodyInertia.h"
#include "SpatialAlgebra/ArticulatedBodyInertia.h"
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
#include "SpatialAlgebra/ForwardDynamics.h"
#include "SpatialAlgebra/RigidBodyInertia.h"
#include "SpatialAlgebra/PluckerTransform.h"
#include <iostream>

using namespace SpatialAlgebra;

int main() {
    // Create a 2-link chain
    ForwardDynamics fd(2);
    
    // Link 1: inertia, transform to parent, joint axis
    RigidBodyInertia I1(1.0, Vector3d(0, 0, 0.5), Matrix3d::Identity() * 0.1);
    Rotation R1;
    R1.setIdentity();
    PluckerTransform X1(R1, Vector3d(0, 0, 0.5));
    Vector3d S1(0, 0, 1); // Rotation about Z axis
    
    fd.setLink(0, I1, X1, S1);
    
    // Link 2: same properties
    RigidBodyInertia I2(1.0, Vector3d(0, 0, 0.5), Matrix3d::Identity() * 0.1);
    Rotation R2;
    R2.setIdentity();
    PluckerTransform X2(R2, Vector3d(0, 0, 0.5));
    Vector3d S2(0, 0, 1);
    
    fd.setLink(1, I2, X2, S2);
    
    // Joint torques
    VectorXd tau(2);
    tau << 1.0, 0.5;
    
    // Compute joint accelerations using ABA
    VectorXd qdd = fd.computeAccelerations(tau);
    
    std::cout << "Joint accelerations:" << std::endl;
    std::cout << qdd.transpose() << std::endl;
    
    return 0;
}
```

## Running Examples

The `examples/` directory contains compilable examples demonstrating core operations:

```sh
# Build all examples
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics

# Run individual examples
./build/example_vectors      # Basic vector operations
./build/example_transforms   # Plücker coordinate transforms
./build/example_inertia      # Rigid body inertia
./build/example_dynamics     # Forward dynamics (ABA)
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

## License

[License information to be added]

## References

- Featherstone, R. (2008). *Rigid Body Dynamics Algorithms*. Springer.
- Featherstone, R. (1987). *Robot Dynamics Algorithms*. Kluwer Academic Publishers.
