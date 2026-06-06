# Phase 18: Robot Examples - Research

**Researched:** 2026-06-06
**Domain:** Robot dynamics examples — Forward Dynamics (ABA) + Inverse Dynamics (RNEA) with gravity cross-validation
**Confidence:** HIGH

## Summary

Phase 18 creates two self-contained example executables demonstrating real-world robot physics using the SpatialAlgebra library. Both examples follow the existing verbose output pattern (`examples/*.cpp`) and are standalone `.cpp` files registered in `examples/CMakeLists.txt`.

**Primary recommendation:** Use UR5-derived link parameters (published masses, COMs, approximate inertias from the ROS industrial URDF and Cambridge 2024 paper) to build two kinematic chains — a 2-link Z-Z planar arm and a 3-link Z-Y-Z spatial RRR arm. Both demonstrate FD (ABA via `ForwardDynamics::computeAccelerations`) and ID (RNEA via `InverseDynamics::computeTorques`) with gravity cross-validation: `ID(q, 0, 0, g) = tau_g → FD(tau_g, g) = qddot ≈ 0`.

**Key physics insight:** The 2-link Z-Z planar arm with gravity [0,0,-9.81] produces ZERO gravity torques (gravity acts parallel to Z joint axes). The cross-validation still validates solver consistency. The 3-link Z-Y-Z arm produces ~30 N·m gravity torque on joint 2 (shoulder pitch), providing the non-trivial demonstration.

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Forward dynamics (ABA) | API / Backend | — | `ForwardDynamics` solver computes qddot from tau and gravity. Examples are consumers, not implementors. |
| Inverse dynamics (RNEA) | API / Backend | — | `InverseDynamics` solver computes tau from qddot and gravity. Examples are consumers. |
| Link model setup | API / Backend | — | Defining Link structs with transforms, inertia, joint axes is the solver setup pattern. |
| Example output / demonstration | Client | — | Each example is a standalone `main()` printing results to stdout. |
| Build registration | Build System | — | `examples/CMakeLists.txt` adds `add_executable` targets. |

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Use UR5-like parameters (mass, inertia, link length) adapted for Z-Z planar and Z-Y-Z spatial configurations.
- **D-02:** Gravity vector = `Vector3d(0, 0, -9.81)` — standard Z-down robotics convention.
- **D-03:** FD+ID cross-validation pattern in both examples.
- **D-04:** Demonstrate static equilibrium: gravity torques exactly cancel gravitational forces, so no acceleration occurs.
- **D-05:** Fully self-contained .cpp files — each example is standalone.
- **D-06:** File names: `examples/example_robot_2link.cpp` and `examples/example_robot_3link.cpp`.
- **D-07:** Executable names: `example_robot_2link` and `example_robot_3link`, linking `SpatialAlgebra Eigen3::Eigen`.
- **D-08:** Verbose output matching existing examples.
- **D-09:** Each example demonstrates: model setup, FD, ID, cross-validation, physical interpretation.
- **D-10:** 2-link Z-Z planar arm: both joints rotate about local Z axis.
- **D-11:** 3-link Z-Y-Z spatial RRR arm: joint 1 about Z (waist yaw), joint 2 about Y (shoulder pitch), joint 3 about Z (elbow roll).
- **D-12:** Cross-validation: `ID(q, 0, 0, gravity) = tau_gravity` → `FD(tau_gravity, gravity) → qddot ≈ 0` (within numerical precision ~1e-12).

### Agent's Discretion
- Exact UR5 parameter values (masses, inertias, link lengths) to use
- Specific static pose(s) used for gravity demonstration
- Number of test cases per example (1-3 different poses/torques)
- Exact cout output format within the verbose style
- Include guards and Doxygen for example files (following existing pattern)
- CMake `add_executable` details in examples/CMakeLists.txt

### Deferred Ideas (OUT OF SCOPE)
None.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| EX-01 | 2-link planar robot dynamics example | UR5 parameters mapped to 2-link Z-Z configuration (see §Standard Stack - Parameter Tables). Uses `ForwardDynamics` / `InverseDynamics` solvers with gravity cross-validation. |
| EX-02 | 3-link spatial arm dynamics example | UR5 parameters mapped to 3-link Z-Y-Z configuration (see §Standard Stack - Parameter Tables). Joint 2 (Y-axis) bears gravity load for non-trivial cross-validation. |
</phase_requirements>

---

## Standard Stack

### Core API
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `ForwardDynamics` | current (Phase 14 fixed) | ABA solver for joint accelerations | Only FD solver in the library |
| `InverseDynamics` | current | RNEA solver for joint torques | Only ID solver in the library |
| `RigidBodyInertia` | current | Rigid body inertia (mass, COM, inertia tensor) | Only inertia API in the library |
| `PluckerTransform` | current | 6×6 spatial coordinate transform | Only spatial transform in the library |
| `MotionVector` | current | Spatial motion vector (twist / screw axis) | Only motion vector API in the library |

### Supporting
| Library | Purpose | When to Use |
|---------|---------|-------------|
| `Rotation` | 3×3 rotation matrix construction | Create identity or rotated link frames |
| `LowerTriangular` | Packed-storage inertia tensor | Construct link rotational inertias |
| `<Eigen/Dense>` | Eigen vector types (`Vector3d`, `VectorXd`) | Torque vectors, state vectors |

### Installation
No new packages — existing library only. Examples link `SpatialAlgebra Eigen3::Eigen` in `examples/CMakeLists.txt`.

### Verified API Signatures

```
[VERIFIED: source code at include/ForwardDynamics.h:156, src/ForwardDynamics.cpp:182]
void ForwardDynamics::computeAccelerations(
    const Eigen::VectorXd& tau,
    const Vector3d& gravity = Vector3d::Zero()
);
// Results stored in fd.links[i].qddot

[VERIFIED: source code at include/InverseDynamics.h:144, src/InverseDynamics.cpp:96]
Eigen::VectorXd InverseDynamics::computeTorques(
    const Eigen::VectorXd& qddot,
    const Vector3d& gravity = Vector3d::Zero()
);
// Returns tau vector

[VERIFIED: source code at include/RigidBodyInertia.h:41]
RigidBodyInertia::RigidBodyInertia(double mass, const Vector3d &com, const lt &inertiaMatrixLT);

[VERIFIED: source code at include/PluckerTransform.h:91]
PluckerTransform::PluckerTransform(const Rotation &rotation, const Vector3d &translation);

[VERIFIED: source code at include/MotionVector.h:86]
MotionVector::MotionVector(const Vector3d &angular, const Vector3d &linear);
// For joint axes: S = MotionVector(axis_direction, Vector3d::Zero()) for revolute joints

[VERIFIED: source code at include/LowerTriangular.h:431]
static LowerTriangular LowerTriangular::Identity(int size);
// Creates identity matrix: lt::Identity(3)

[VERIFIED: source code at include/Rotation.h:64-69]
Rotation();  // Default constructor = identity
Rotation(const Eigen::Matrix3d &matrix);  // From full matrix
// Also: R.setIdentity() via Eigen::Matrix3d inheritance
```

### UR5 Parameter Values

**Source:** [CITED: https://github.com/ros-industrial/universal_robot — UR5 URDF xacro with masses/COMs] and [CITED: Cambridge University Press 2024 — "Electromechanical Modeling and Identification of the UR5e" (Table II, Table III)]

| Link | Mass (kg) | COM (m) | Rotational Inertia (kg·m²) [Ixx, Iyy, Izz, Ixy, Ixz, Iyz] | DH length (m) |
|------|-----------|---------|------------------------------------------------------------|----------------|
| 1 (shoulder) | 3.7 | [0, -0.02561, 0.00193] | [0.0067, 0.0064, 0.0067, 0, 0, 0] | d1 = 0.089 |
| 2 (upper arm) | 8.393 | [0.2125, 0, 0.11336] | [0.0149, 0.3564, 0.3553, 0, 0, 0] | a2 = 0.425 |
| 3 (forearm) | 2.33 | [0.15, 0, 0.0265] | [0.0025, 0.0551, 0.0546, 0, 0, 0.0034] | a3 = 0.392 |

---

## Package Legitimacy Audit

No external packages are installed for this phase. The phase creates example `.cpp` files using only the existing library API (already built and verified in prior phases). Both `ForwardDynamics` and `InverseDynamics` are compiled from `src/` and linked via `SpatialAlgebra` static library.

**Packages removed due to slopcheck [SLOP] verdict:** none
**Packages flagged as suspicious [SUS]:** none

---

## Architecture Patterns

### System Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                      Example Executables                        │
│                                                                 │
│  ┌─────────────────────┐     ┌─────────────────────────────┐    │
│  │example_robot_2link  │     │  example_robot_3link        │    │
│  │                     │     │                             │    │
│  │ 1. Build Z-Z arm    │     │  1. Build Z-Y-Z arm        │    │
│  │ 2. FD(τ) → qddot    │     │  2. FD(τ) → qddot          │    │
│  │ 3. ID(0, g) → τ_g   │     │  3. ID(0, g) → τ_g         │    │
│  │ 4. FD(τ_g, g)→0     │     │  4. FD(τ_g, g)→0           │    │
│  └────────┬────────────┘     └───────────┬─────────────────┘    │
│           │                              │                      │
└───────────┼──────────────────────────────┼──────────────────────┘
            │                              │
            ▼                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   SpatialAlgebra Library                         │
│                                                                 │
│  ┌─────────────────┐  ┌─────────────────┐                       │
│  │ ForwardDynamics  │  │ InverseDynamics │                       │
│  │ (ABA, O(n))      │  │ (RNEA, O(n))    │                       │
│  │ computeAccels()  │  │ computeTorques() │                      │
│  └────────┬─────────┘  └────────┬────────┘                      │
│           │                     │                                │
│           ▼                     ▼                                │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │  Core Types: Link / InverseDynamicsLink, PluckerTransform │   │
│  │  RigidBodyInertia, MotionVector (joint axes), Rotation    │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
            │                              │
            ▼                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Eigen3 (linear algebra backend)                                │
│  - Vector3d, VectorXd, Matrix3d, MatrixXd                       │
│  - Rotation via Eigen::Matrix3d inheritance                     │
└─────────────────────────────────────────────────────────────────┘
```

Data flow for cross-validation:
```
tau_g = ID(q, qdot=0, qddot=0, gravity)
  ──→ qddot_result = FD(tau_g, gravity)
  ──→ verify |qddot_result| < 1e-10
```

### Recommended Project Structure

No new files beyond the two examples and CMakeLists.txt update:

```
examples/
├── CMakeLists.txt           # Add 2 new add_executable entries (existing file, edit)
├── example_robot_2link.cpp  # NEW: Z-Z planar arm example
├── example_robot_3link.cpp  # NEW: Z-Y-Z spatial arm example
├── dynamics.cpp             # Existing (unchanged)
├── basic_vectors.cpp        # Existing (unchanged)
├── transforms.cpp           # Existing (unchanged)
└── inertia.cpp              # Existing (unchanged)
```

### Pattern 1: Example File Structure

**What:** Self-contained `main()` with Doxygen file comment, includes, type aliases, and verbose output sections. Matches the existing pattern from `examples/dynamics.cpp`.

**When to use:** All new example files.

**Example structure:**

```cpp
/**
 * @file example_robot_2link.cpp
 * @brief Demonstrates forward and inverse dynamics for a 2-link Z-Z planar arm.
 * 
 * This example shows:
 * - Setting up a 2-link kinematic chain with UR5-like parameters
 * - Forward dynamics (ABA) with applied joint torques
 * - Inverse dynamics (RNEA) for gravity compensation
 * - FD+ID cross-validation verifying static equilibrium
 * 
 * Robot Configuration:
 * - Joint 1: Z-axis revolute (waist)
 * - Joint 2: Z-axis revolute (elbow)
 * - Arm operates in XY plane (horizontal)
 * - Gravity: [0, 0, -9.81] m/s²
 */

#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include "LowerTriangular.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
using lt = LowerTriangular;
using mv = MotionVector;
using fv = ForceVector;
using plux = PluckerTransform;
using rbi = RigidBodyInertia;

int main() {
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "=== SpatialAlgebra 2-Link Z-Z Planar Arm Example ===" << std::endl;
    // ...
    return 0;
}
```

### Concrete Parameter Tables

#### 2-Link Z-Z Planar Arm Parameters

| Property | Link 0 (base — parent=-1) | Link 1 (tip — parent=0) |
|----------|--------------------------|-------------------------|
| Mass (kg) | 8.393 | 2.33 |
| COM (m) | [0.2125, 0, 0.02] | [0.15, 0, 0.02] |
| Inertia Ixx | 0.0149 | 0.0025 |
| Inertia Iyy | 0.3564 | 0.0551 |
| Inertia Izz | 0.3553 | 0.0546 |
| Parent | -1 | 0 |
| X from parent | `PluckerTransform(Identity, [0,0,0])` | `PluckerTransform(Identity, [0.425, 0, 0])` |
| Joint axis S | `MotionVector([0,0,1], [0,0,0])` (Z) | `MotionVector([0,0,1], [0,0,0])` (Z) |
| Source UR5 link | Link 2 (upper arm) | Link 3 (forearm) |

#### 3-Link Z-Y-Z Spatial Arm Parameters

| Property | Link 0 (waist) | Link 1 (shoulder) | Link 2 (elbow) |
|----------|----------------|-------------------|----------------|
| Mass (kg) | 3.7 | 8.393 | 2.33 |
| COM (m) | [0, -0.02561, 0.00193] | [0.2125, 0, 0.11336] | [0.15, 0, 0.0265] |
| Inertia Ixx | 0.0067 | 0.0149 | 0.0025 |
| Inertia Iyy | 0.0064 | 0.3564 | 0.0551 |
| Inertia Izz | 0.0067 | 0.3553 | 0.0546 |
| Inertia Ixz | 0 | 0 | 0.0034 |
| Parent | -1 | 0 | 1 |
| X from parent | `PluckerTransform(Identity, [0,0,0])` | `PluckerTransform(Identity, [0,0,0.089])` | `PluckerTransform(Identity, [0.425, 0, 0])` |
| Joint axis S | `[0,0,1]` (Z, waist yaw) | `[0,1,0]` (Y, shoulder pitch) | `[0,0,1]` (Z, elbow roll) |
| Source UR5 link | Link 1 (shoulder) | Link 2 (upper arm) | Link 3 (forearm) |

### Inertia Tensor Construction Pattern

```cpp
// Example: link 0 inertia (from UR5 link 2: Ixx=0.0149, Iyy=0.3564, Izz=0.3553)
lt I(3);
I(0, 0) = 0.0149;  // Ixx
I(1, 0) = 0.0;     // Ixy
I(1, 1) = 0.3564;  // Iyy
I(2, 0) = 0.0;     // Ixz
I(2, 1) = 0.0;     // Iyz
I(2, 2) = 0.3553;  // Izz
link.I = RigidBodyInertia(mass, com, I);
```

For links with non-zero Ixz (link 3: Ixz=0.0034):
```cpp
I(2, 0) = 0.0034;  // Ixz
```

### Cross-Validation Pattern

```cpp
// Step 1: Setup both solvers with identical chain
ForwardDynamics fd;
InverseDynamics id;
// ... populate fd.links and id.links identically ...

// Step 2: ID at static pose — compute gravity compensation torques
Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(nDof);
Vector3d gravity(0, 0, -9.81);
Eigen::VectorXd tau_gravity = id.computeTorques(qddot_zero, gravity);

// Step 3: FD with gravity torques — verify zero acceleration
fd.computeAccelerations(tau_gravity, gravity);

// Step 4: Check
for (int i = 0; i < nDof; i++) {
    std::cout << "  Joint " << i << " qddot: " << fd.links[i].qddot;
    std::cout << "  (expected ~0, tolerance 1e-10)" << std::endl;
}
```

### Anti-Patterns to Avoid
- **Sharing setup code between examples:** D-05 requires fully self-contained `.cpp` files. Do not create a shared `robot_utils.h`.
- **Modifying library code:** All work is in `examples/` — no changes to `include/`, `src/`, or `tests/`.
- **Using ForwardDynamicsLink and InverseDynamicsLink inconsistently:** `ForwardDynamics::links` stores `Link` (which is `ForwardDynamicsLink`), `InverseDynamics::links` stores `InverseDynamicsLink`. These are separate struct types even though they share fields. The solver models must be built independently.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Forward dynamics | Custom O(n²) solver | `ForwardDynamics::computeAccelerations()` | ABA is O(n), tested, gravity-aware |
| Inverse dynamics | Custom RNEA | `InverseDynamics::computeTorques()` | RNEA is O(n), tested, gravity-aware |
| Inertia construction | Manual 6×6 matrix assembly | `RigidBodyInertia(mass, com, lt)` | Handles COM coupling correctly |
| Spatial transforms | Manual 6×6 matrix multiply | `PluckerTransform::transformMotion()` | Verified Plücker convention |

**Key insight:** All the heavy dynamics machinery already exists in the library. The examples are consumers only.

---

## Gravity Cross-Validation Math

### How Gravity is Handled Internally

Both solvers implement the Featherstone gravity convention:

**InverseDynamics** (see `src/InverseDynamics.cpp:33-34`):
```cpp
// Base link: a₀ = S₀·q̈₀ − [0; g]  (gravity subtracts from linear base acceleration)
links[i].a = links[i].S * links[i].qddot
           - MotionVector(Vector3d::Zero(), gravity);
```

**ForwardDynamics** (see `src/ForwardDynamics.cpp:67`):
```cpp
// Base link: c₀ = −[0; g]  (bias acceleration includes equivalent upward acceleration)
links[i].c = MotionVector(Vector3d::Zero(), -this->gravity);
```

Both treat gravity as an upward equivalent acceleration applied to the base link's linear component. This is the standard Featherstone convention (D-08).

### Cross-Validation Identity

For any static configuration (qdot=0, qddot=0):

1. **ID computes gravity torque:** `tau_g = ID(q, 0, 0, g)` = torques needed to hold the arm in static equilibrium against gravity
2. **FD with gravity torque:** `FD(tau_g, g)` should produce qddot ≈ 0

This identity holds because both solvers use the same gravity convention. If they disagree, there is a bug.

### Expected Precision

From existing test infrastructure (`tests/TestDynamicsConsistency.cpp`):
- Single-link round-trip with gravity: `EXPECT_NEAR(..., 1e-8)` [VERIFIED: test code line 303]
- Single-link direct comparison: `EXPECT_NEAR(..., 1e-10)` [VERIFIED: test code line 422]
- Multi-link round-trip: `EXPECT_NEAR(..., 1e-8)` [VERIFIED: test code line 175]

For examples output, typical qddot values after cross-validation should be on the order of **1e-12 to 1e-10**. The planner should check for `abs(qddot) < 1e-8` as a pass criterion.

### 2-Link Z-Z Arm: Zero Gravity Torque

**Key physics finding:** For a Z-Z planar arm operating in the XY plane with gravity [0, 0, -9.81], gravity acts PARALLEL to both joint axes. The cross product r × (m·g) has zero component along Z, so gravity produces NO torque about Z-axis joints. The gravity compensation torques are exactly [0, 0] for any pose.

This is physically correct behavior — the arm is a horizontal SCARA-like configuration. The cross-validation still works:
- `tau_g = ID(q, 0, 0, g)` = [0, 0] (correct: no gravity torque needed)
- `FD([0, 0], g)` → qddot ≈ [0, 0] (correct: arm stays at rest)

The example should explain this physics and note that the interesting gravity case is in the 3-link Z-Y-Z arm.

### 3-Link Z-Y-Z Arm: Approximate Gravity Torque

For the 3-link arm at q=[0, 0, 0] (waist at 0, arm extended horizontally, elbow at 0):

- Joint 1 (Z, waist yaw): τ₁ ≈ 0 N·m (gravity along Z, joint axis Z → zero torque)
- Joint 2 (Y, shoulder pitch): τ₂ ≈ m₂·g·x₂ + m₃·g·(L₂ + x₃)
  - = 8.393 × 9.81 × 0.2125 + 2.33 × 9.81 × (0.425 + 0.15)
  - = 17.5 + 13.2 = **~30.7 N·m**
- Joint 3 (Z, elbow roll): τ₃ ≈ 0 N·m (gravity along Z, joint axis Z → zero torque)

The gravity torque on joint 2 is substantial (~30 N·m for arm extended horizontally, or ~15 N·m for arm at 45°). This demonstrates physically meaningful gravity compensation.

---

## Common Pitfalls

### Pitfall 1: Incorrect Parent-Child Transform Direction
**What goes wrong:** The Plücker transform `X` from parent to child describes how to transform a motion vector from the parent frame to the child frame. Confusing this direction produces wrong dynamics.
**How to avoid:** The transform `X` is the transform FROM parent TO child — a vector expressed in the parent frame is transformed to the child frame via `X.transformMotion(v_parent)`. This is already established in existing examples (`dynamics.cpp:72`). [VERIFIED: src/ForwardDynamics.cpp:72]
**Warning signs:** Joint accelerations are wildly wrong (orders of magnitude off).

### Pitfall 2: Using `Link` vs `InverseDynamicsLink`
**What goes wrong:** `ForwardDynamics::links` holds `Link` (aliased as `ForwardDynamicsLink`), while `InverseDynamics::links` holds `InverseDynamicsLink`. These are separate struct types. Building an FD model and passing it to ID will not compile.
**How to avoid:** Build separate FD and ID link vectors, even with identical parameters. The existing test code (`TestDynamicsConsistency.cpp`) demonstrates this pattern — note how `id_link0` and `fd_link0` are separate variables. [VERIFIED: test code lines 105-111 and 144-150]
**Warning signs:** Compiler error "cannot convert ForwardDynamicsLink to InverseDynamicsLink" or vice versa.

### Pitfall 3: Forgetting the Z-Z Arm Has Zero Gravity Torque
**What goes wrong:** The planner/implementor expects non-zero gravity compensation and writes incorrect physics commentary.
**How to avoid:** The 2-link Z-Z arm with gravity [0,0,-9.81] has gravity parallel to joint axes → zero gravity torque. This is correct. Document this clearly in the output.
**Warning signs:** Comments claiming "gravity torque on Z-axis joint" for a Z-Z arm with Z-down gravity.

### Pitfall 4: Mismatched Link Counts Between FD and ID
**What goes wrong:** FD and ID models must have the **same number of links** with the **same topology** (same parent indices, same transforms). Any mismatch breaks the cross-validation.
**How to avoid:** Build both models in parallel setup sections. The `computeTorques(qddot)` and `computeAccelerations(tau)` methods validate that input vector size matches `links.size()`.
**Warning signs:** `std::invalid_argument` exception at runtime.

---

## Code Examples

### Verified Pattern: Identity Rotation Construction

Two patterns used interchangeably in existing code:

```cpp
// Pattern A (used in examples/dynamics.cpp:53-54):
Rotation R;
R.setIdentity();

// Pattern B (used in tests/TestDynamicsConsistency.cpp:25):
Rotation R(Eigen::Matrix3d::Identity());
```

Both produce the same result. Pattern B is preferred for one-liners.

### Verified Pattern: LowerTriangular Inertia Construction

```cpp
// Source: examples/dynamics.cpp:47-49
lt I_tensor(3);
I_tensor(0, 0) = 0.1; I_tensor(1, 0) = 0.0; I_tensor(1, 1) = 0.1;
I_tensor(2, 0) = 0.0; I_tensor(2, 1) = 0.0; I_tensor(2, 2) = 0.1;
link.I = RigidBodyInertia(mass, Vector3d::Zero(), I_tensor);
```

### Verified Pattern: Joint Axis Construction

```cpp
// Revolute joint about Z axis (used in dynamics.cpp:58):
link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
// Revolute joint about Y axis:
link.S = MotionVector(Vector3d(0, 1, 0), Vector3d::Zero());
```

### Verified Pattern: CMakeLists.txt Registration

```cmake
# Source: examples/CMakeLists.txt:21-23 (existing pattern)
add_executable(example_dynamics dynamics.cpp)
target_link_libraries(example_dynamics SpatialAlgebra Eigen3::Eigen)
```

New entries follow the same pattern:
```cmake
add_executable(example_robot_2link example_robot_2link.cpp)
target_link_libraries(example_robot_2link SpatialAlgebra Eigen3::Eigen)

add_executable(example_robot_3link example_robot_3link.cpp)
target_link_libraries(example_robot_3link SpatialAlgebra Eigen3::Eigen)
```

### Verified Pattern: Gravity Round-Trip (from TestDynamicsConsistency.cpp)

```cpp
// Source: tests/TestDynamicsConsistency.cpp:269-303
Vector3d gravity(0, 0, -9.81);

// ID
InverseDynamics id;
// ... setup links identically for ID and FD ...
Eigen::VectorXd tau = id.computeTorques(qddot_input, gravity);

// FD
ForwardDynamics fd;
// ... setup links identically for ID and FD ...
fd.computeAccelerations(tau, gravity);

EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], 1e-8);
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| No examples with real robot parameters | UR5-based examples (this phase) | Phase 18 | Demonstrates library usefulness |
| Zero-gravity FD only | Gravity-aware FD+ID cross-validation | Phase 14 | Enables this phase |
| Separate `examples/CMakeLists.txt` per example | Batch examples in single CMakeLists.txt | v1.0 | Pattern established |

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | UR5 link masses and COMs from ROS URDF are close approximations for real UR5 | Standard Stack - Parameter Tables | Low — D-01 says "UR5-LIKE", exact precision isn't critical for demo |
| A2 | Inertia values from Cambridge 2024 paper approximated for UR5e are suitable for UR5 examples | Standard Stack - Parameter Tables | Low — inertias affect dynamics quantitatively, but cross-validation tests solver consistency, not absolute accuracy |
| A3 | The 2-link Z-Z arm with gravity [0,0,-9.81] produces zero gravity torque | Gravity Cross-Validation Math | MEDIUM — if the library convention differs from my physics analysis, tau_g may be non-zero. This should be verified by running the compiled example. |

---

## Open Questions (RESOLVED)

1. **What poses to use for each example?** (RESOLVED by 18-01-PLAN.md)
   - What we know: 2-link arm zero position (both q=0) and one rotated pose; 3-link arm with horizontal extension
   - Resolution: 2-link uses q=[0,0] and q=[π/4, π/3]; 3-link uses q=[0,0,0] (horizontal) and q=[0,π/4,0] (45°)

2. **How many torque configurations per example?** (RESOLVED by 18-01-PLAN.md)
   - What we know: existing dynamics.cpp shows 3 torque configurations
   - Resolution: 2 FD tests + 1 ID gravity test + 1 cross-validation test per example, matching the 3 torque configurations from dynamics.cpp plus the dedicated gravity section

---

## Environment Availability

> Skipped — this phase creates `.cpp` source files and updates CMakeLists.txt. No external dependencies beyond the existing build toolchain (cmake, compiler, Eigen3), which are already verified working from prior phases.

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | CMake build + runtime execution |
| Config file | `examples/CMakeLists.txt` |
| Quick run command | `cmake --build build && build/example_robot_2link && build/example_robot_3link` |
| Full suite command | `cmake --build build && build/example_robot_2link > /dev/null && build/example_robot_3link > /dev/null` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| EX-01 | 2-link example compiles and runs without crash | smoke | `cmake --build build && build/example_robot_2link` | ❌ Wave 0 |
| EX-02 | 3-link example compiles and runs without crash | smoke | `cmake --build build && build/example_robot_3link` | ❌ Wave 0 |

### Verification Protocol (not automated tests — examples are manually inspected)

The examples output to stdout with no assertions. Verification is by inspection:
1. **Compilation:** `cmake --build build` succeeds (no link errors)
2. **Execution:** Both executables run to completion (return 0)
3. **Numerical cross-validation:** FD+ID round-trip produces qddot values < 1e-8 in magnitude
4. **Physical interpretation:** Output includes explanatory text matching expected physics
5. **Exit codes:** Both executables `return 0` on success

### Sampling Rate
- **Per task commit:** `cmake --build build` compilation check
- **Per wave merge:** Full runtime test of both executables
- **Phase gate:** Both examples compile, run, and produce physically correct output before `/gsd-verify-work`

### Wave 0 Gaps
- [ ] `examples/example_robot_2link.cpp` — does not exist yet (this phase creates it)
- [ ] `examples/example_robot_3link.cpp` — does not exist yet (this phase creates it)
- [ ] `examples/CMakeLists.txt` — edit to add two new targets

---

## Security Domain

> Not applicable — this phase creates example executable code that runs locally with no network access, no user input, and no persistent state. No authentication, authorization, input validation, or cryptography concerns.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | — |
| V3 Session Management | no | — |
| V4 Access Control | no | — |
| V5 Input Validation | no | — |
| V6 Cryptography | no | — |

---

## Sources

### Primary (HIGH confidence)
- `include/ForwardDynamics.h` + `src/ForwardDynamics.cpp` — API signatures and gravity implementation
- `include/InverseDynamics.h` + `src/InverseDynamics.cpp` — API signatures and gravity implementation
- `include/RigidBodyInertia.h` — Constructor, `apply()`, member access
- `include/PluckerTransform.h` — Constructor, transform methods
- `include/MotionVector.h` — Constructor for joint screw axes
- `include/LowerTriangular.h` — `Identity()`, element access, packed storage
- `include/Rotation.h` — Default constructor, `setIdentity()`, matrix construction
- `examples/dynamics.cpp` — Example file structure, link setup, ABA usage pattern
- `examples/basic_vectors.cpp` — Verbose output pattern, type aliases, Doxygen file comment
- `examples/inertia.cpp` — Inertia construction, COM offset pattern
- `examples/transforms.cpp` — Plücker transform usage, Rotation construction
- `examples/CMakeLists.txt` — Build target registration pattern
- `tests/TestDynamicsConsistency.cpp` — Gravity round-trip test pattern, precision expectations

### Secondary (MEDIUM confidence)
- [CITED: https://github.com/ros-industrial/universal_robot/blob/kinetic-devel/ur_description/urdf/ur5.urdf.xacro] — UR5 link masses, COMs, DH parameters
- [CITED: https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/] — Official UR5 DH parameters and dynamics data
- [CITED: Cambridge University Press 2024 — "Electromechanical modeling and identification of the UR5 e-series robot" (Tables II, III)] — UR5e link inertias (estimated geometrically)

### Tertiary (LOW confidence)
- None — all claims are either verified against source code or cited from published documentation.

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all API signatures verified from source code
- Architecture: HIGH — patterns verified from existing examples
- Pitfalls: HIGH — verified from test code and physics analysis
- UR5 parameters: MEDIUM — URDF data is authoritative, inertia values are geometric estimates from literature

**Research date:** 2026-06-06
**Valid until:** Project stable — no library changes expected in this area (Phase 18 is the last examples phase before RBDL comparison).
