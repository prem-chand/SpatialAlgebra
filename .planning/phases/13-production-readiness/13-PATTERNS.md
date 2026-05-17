# Phase 13: Production Readiness — Pattern Map

**Mapped:** 2026-05-17
**Files analyzed:** 33 (2 new, 28 modified, 3 removed)
**Analogs found:** 31 / 31

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `include/SpatialAlgebra.h` | config | request-response | no analog (new) | — |
| `.github/workflows/ci.yml` | config | event-driven | no analog (new) | — |
| `include/SpatialUtils.h` | utility | transform | itself (current file) | exact |
| `include/SpatialVector.h` | model | CRUD | itself (current file) | exact |
| `src/SpatialVector.cpp` | service | CRUD | itself (current file) | exact |
| `include/ForceVector.h` | model | CRUD | itself (current file) | exact |
| `src/ForceVector.cpp` | service | CRUD | itself (current file) | exact |
| `include/MotionVector.h` | model | CRUD | itself (current file) | exact |
| `src/MotionVector.cpp` | service | CRUD | itself (current file) | exact |
| `include/ArticulatedBodyInertia.h` | model | CRUD | itself (current file) | exact |
| `include/ForwardDynamics.h` | service | transform | itself (current file) | exact |
| `src/ForwardDynamics.cpp` | service | transform | itself (current file) | exact |
| `include/InverseDynamics.h` | service | transform | itself (current file) | exact |
| `src/InverseDynamics.cpp` | service | transform | itself (current file) | exact |
| `include/PluckerTransform.h` | model | transform | itself (current file) | exact |
| `src/PluckerTransform.cpp` | service | transform | itself (current file) | exact |
| `include/LowerTriangular.h` | model | CRUD | itself (current file) | exact |
| `src/LowerTriangular.cpp` | service | CRUD | itself (current file) | exact |
| `include/SpatialOperations.h` | utility | transform | itself (current file) | exact |
| `src/SpatialOperations.cpp` | utility | transform | itself (current file) | exact |
| `include/RigidBodyInertia.h` | model | CRUD | itself (current file) | exact |
| `CMakeLists.txt` | config | build | itself (current file) | exact |
| `tests/TestSpatialOperations.cpp` | test | CRUD | itself (current file) | exact |
| `tests/TestForwardDynamics.cpp` | test | transform | itself (current file) | exact |
| `tests/TestInverseDynamics.cpp` | test | transform | itself (current file) | exact |
| `tests/TestDynamicsConsistency.cpp` | test | transform | itself (current file) | exact |
| `tests/TestSpatialVector.cpp` | test | CRUD | itself (current file) | exact |
| `tests/TestRigidBodyInertia.cpp` | test | CRUD | itself (current file) | exact |
| `tests/TestArticulatedBodyInertia.cpp` | test | CRUD | itself (current file) | exact |
| `src/RigidBodyInertia.cpp` | — | — | — | REMOVE |
| `src/ArticulatedBodyInertia.cpp` | — | — | — | REMOVE |

## Pattern Assignments

### `include/SpatialAlgebra.h` (config, request-response) — NEW FILE

**Analog:** None found (first umbrella header). Use RESEARCH.md patterns.

**Pattern to follow (umbrella header style per D-28):**
```cpp
#ifndef SPATIAL_ALGEBRA_H
#define SPATIAL_ALGEBRA_H

/**
 * @file SpatialAlgebra.h
 * @brief Umbrella header for the SpatialAlgebra library
 * @details Includes all public headers in dependency order.
 *          Users can include this single header instead of individual headers.
 */

#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "Rotation.h"
#include "LowerTriangular.h"
#include "SpatialUtils.h"
#include "SpatialOperations.h"
#include "RigidBodyInertia.h"
#include "ArticulatedBodyInertia.h"
#include "PluckerTransform.h"
#include "ForwardDynamics.h"
#include "InverseDynamics.h"

#endif // SPATIAL_ALGEBRA_H
```

**Convention reference:** All existing headers use `#ifndef`/`#define` include guards matching the filename (`SPATIAL_VECTOR_H`, `FORWARD_DYNAMICS_H`, etc.). See `include/SpatialVector.h:1-2`.

---

### `.github/workflows/ci.yml` (config, event-driven) — NEW FILE

**Analog:** None found (no existing CI workflow in repo). Use RESEARCH.md patterns.

**Pattern to follow (from D-09 through D-12, D-11):**
```yaml
name: CI

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  build:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        compiler: [g++, clang++]
    
    runs-on: ${{ matrix.os }}
    
    steps:
    - uses: actions/checkout@v4
    
    - name: Install dependencies
      run: |
        brew install eigen googletest
    
    - name: Configure
      run: cmake -B build -DCMAKE_CXX_COMPILER=${{ matrix.compiler }}
    
    - name: Build
      run: cmake --build build
    
    - name: Test
      run: cd build && ctest --output-on-failure
```

---

### `include/SpatialUtils.h` (utility, transform) — MODIFIED

**Analog:** `include/SpatialUtils.h` (current file, lines 115-125)

**Current cross() force×force pattern (buggy — missing `-τ2×f1` term — already fixed in file on disk? Let me verify):**
Current file shows correct formula at lines 115-125. The research says the existing free function was buggy, but the file on disk already has the correct formula. This may have been partially fixed. The key change is D-01: make this the canonical implementation that all other methods delegate to.

**Current pattern (correct formula already in-place at lines 115-125):**
```cpp
inline ForceVector cross(const ForceVector& v1, const ForceVector& v2) noexcept {
    const Vector3d& t1 = v1.getAngular();
    const Vector3d& f1 = v1.getLinear();
    const Vector3d& t2 = v2.getAngular();
    const Vector3d& f2 = v2.getLinear();
    return ForceVector(
        t1.cross(t2) + f1.cross(f2),
        t1.cross(f2) - t2.cross(f1)
    );
}
```

**Existing free function overloads pattern (lines 77-106):**
```cpp
inline ForceVector cross(const MotionVector& v1, const ForceVector& v2) noexcept { ... }
inline MotionVector cross(const MotionVector& v1, const MotionVector& v2) noexcept { ... }
```

**Header structure (lines 1-17, 127-128):**
```cpp
#ifndef SPATIAL_UTILS_H
#define SPATIAL_UTILS_H
// ...
#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
namespace SpatialAlgebra {
// ... free functions ...
} // namespace SpatialAlgebra
#endif // SPATIAL_UTILS_H
```

---

### `src/SpatialVector.cpp` (service, CRUD) — MODIFIED

**Analog:** `src/SpatialVector.cpp` (current file, lines 48-54)

**Current crossForce pattern to MODIFY (D-01: delegate to canonical):**
```cpp
SpatialVector SpatialAlgebra::SpatialVector::crossForce(const SpatialVector &other) const
{
    return SpatialVector(
        angular.cross(other.angular) + linear.cross(other.linear),
        angular.cross(other.linear) - other.angular.cross(this->linear)
    );
}
```

**Target pattern (delegate to canonical free function):**
```cpp
SpatialVector SpatialAlgebra::SpatialVector::crossForce(const SpatialVector &other) const
{
    // Delegate to canonical free function in SpatialUtils.h
    return ForceVector(cross(ForceVector(*this), ForceVector(other)));
}
```

**Constructor pattern (lines 9-15):**
```cpp
SpatialVector::SpatialVector() : angular(Vector3d::Zero()), linear(Vector3d::Zero()) {}

SpatialVector::SpatialVector(const Vector3d &a, const Vector3d &l)
    : angular(a), linear(l) {}
```

---

### `src/ForceVector.cpp` (service, CRUD) — MODIFIED

**Analog:** `src/ForceVector.cpp` (current file)

**Current crossForce pattern to MODIFY (D-01: delegate to canonical):**
```cpp
ForceVector ForceVector::crossForce(const ForceVector &other) const   // lines 40-46
{
    return ForceVector(
        this->angular.cross(other.angular) + this->linear.cross(other.linear),
        this->angular.cross(other.linear) - other.angular.cross(this->linear)
    );
}
```

**Target pattern (delegate):**
```cpp
ForceVector ForceVector::crossForce(const ForceVector &other) const
{
    return cross(*this, other);  // Delegate to canonical free function
}
```

**Remove crossMotion (D-04): eliminate the entire `ForceVector::crossMotion` method at lines 31-38.**

**Constructor pattern (lines 6-12):**
```cpp
ForceVector::ForceVector() : SpatialVector() {}
ForceVector::ForceVector(const Vector3d &angular, const Vector3d &linear)
    : SpatialVector(angular, linear) {}
ForceVector::ForceVector(const SpatialVector &other)
    : SpatialVector(other) {}
```

---

### `include/ForceVector.h` (model, CRUD) — MODIFIED

**Analog:** `include/ForceVector.h` (current file, line 130)

**Remove declaration (D-04):** `ForceVector crossMotion(const ForceVector &other) const;` at line 130.

**Keep declaration:** `ForceVector crossForce(const ForceVector &other) const;` at line 140 (implementation now delegates).

---

### `include/MotionVector.h` (model, CRUD) — MODIFIED

**Analog:** `include/MotionVector.h` (current file, line 137)

**Remove declaration (D-04):** `MotionVector crossForce(const MotionVector &other) const;` at line 137.

**Keep declaration:** `MotionVector crossMotion(const MotionVector &other) const;` at line 129.

---

### `src/MotionVector.cpp` (service, CRUD) — MODIFIED

**Analog:** `src/MotionVector.cpp` (current file, lines 38-41)

**Remove implementation (D-04):** Eliminate entire `MotionVector::crossForce` method at lines 38-41.

---

### `include/ArticulatedBodyInertia.h` (model, CRUD) — MODIFIED

**Analog:** `include/ArticulatedBodyInertia.h` (current file, lines 150-155)

**Current operator+(RigidBodyInertia) pattern to FIX (CR-01 — swapped args, missing mass multiplier):**
```cpp
inline ArticulatedBodyInertia operator+(const RigidBodyInertia &other) const   // lines 150-155
{
    return ArticulatedBodyInertia(Inertia + other.getInertiaMatrixLT(),
                                H + other.getMass() * skew(other.getCom()),
                                M + lt::Identity(3) * other.getMass());
}
```

**Target pattern (fix: swap first/third args, add mass multiplier to coupling term):**
```cpp
inline ArticulatedBodyInertia operator+(const RigidBodyInertia &other) const
{
    double m = other.getMass();
    Vector3d c = other.getCom();
    return ArticulatedBodyInertia(
        Inertia + other.getInertiaMatrixLT(),
        H + skew(c) * m,
        M + lt::Identity(3) * m
    );
}
```

**Inline class pattern (entirely inline in header — lines 77-207):**
```cpp
class ArticulatedBodyInertia
{
private:
    lt Inertia;
    Eigen::Matrix3d H;
    lt M;
public:
    ArticulatedBodyInertia(const lt &inertia, const Eigen::Matrix3d &h, const lt &M)
        : Inertia(inertia), H(h), M(M) {}
    // ... all methods inline ...
};
```

---

### `include/ForwardDynamics.h` (service, transform) — MODIFIED

**Analog:** `include/ForwardDynamics.h` (current file, lines 139-155)

**Current computeAccelerations signature to MODIFY (D-07 — add gravity param):**
```cpp
void computeAccelerations(const Eigen::VectorXd& tau);   // line 155
```

**Target pattern:**
```cpp
void computeAccelerations(const Eigen::VectorXd& tau, 
                          const Vector3d& gravity = Vector3d::Zero());
```

**Link struct pattern (lines 79-111):**
```cpp
struct Link
{
    int parent;
    PluckerTransform X;
    RigidBodyInertia I;
    MotionVector S;
    double q, qdot, qddot;
    MotionVector v, c;
    ForceVector f;
    ArticulatedBodyInertia Ia;
    ForceVector pa;
    
    Link() : parent(-1),
             X(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero()),
             I(1.0, Vector3d::Zero(), lt::Identity(3)),
             S(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             q(0.0), qdot(0.0), qddot(0.0),
             v(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             c(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             f(ForceVector(Vector3d::Zero(), Vector3d::Zero())),
             Ia(lt::Identity(3), Eigen::Matrix3d::Zero(), lt::Identity(3)),
             pa(ForceVector(Vector3d::Zero(), Vector3d::Zero())) {}
};
```

---

### `src/ForwardDynamics.cpp` (service, transform) — MODIFIED

**Analog:** `src/ForwardDynamics.cpp` (current file)

**Current outwardPass pattern (lines 19-51):**
```cpp
void ForwardDynamics::outwardPass()
{
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        int parent = links[i].parent;
        if (parent == -1)
        {
            links[i].v = links[i].S * links[i].qdot;
            links[i].c = MotionVector(Vector3d::Zero(), Vector3d::Zero());
        }
        else
        {
            MotionVector vParent = links[parent].v;
            links[i].v = links[i].X.transformMotion(vParent) + links[i].S * links[i].qdot;
            MotionVector cParent = links[parent].c;
            links[i].c = links[i].X.transformMotion(cParent) + 
                         cross(links[i].v, links[i].S) * links[i].qdot;
        }
    }
}
```

**Current two-phase inwardPass (already implemented in file at lines 53-126 with the fix):**
The current file on disk has already been restructured into the clean two-phase approach (D-05). Phase 1 initializes Ia/pa (lines 56-73), Phase 2 accumulates (lines 77-93), Phase 3 solves (lines 95-125). This pattern should be preserved.

**Current computeAccelerations pattern to MODIFY (add gravity):**
```cpp
void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau)
{
    // Validate input
    // ...
    outwardPass();
    inwardPass(tau);
}
```

**Target computeAccelerations (with gravity per D-08):**
```cpp
void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau,
                                           const Vector3d& gravity)
{
    // Store gravity for use in outwardPass
    this->gravity = gravity;
    // ...
    outwardPass();
    inwardPass(tau);
}
```

**Validation pattern (lines 131-150):**
```cpp
if (tau.size() != static_cast<int>(links.size()))
{
    throw std::invalid_argument(
        "ForwardDynamics::computeAccelerations: tau size (" + 
        std::to_string(tau.size()) + ") does not match link count (" + 
        std::to_string(links.size()) + ")"
    );
}
```

---

### `include/InverseDynamics.h` (service, transform) — MODIFIED

**Analog:** `include/InverseDynamics.h` (current file, lines 127-143)

**Current computeTorques signature to MODIFY (D-07 — add gravity param):**
```cpp
Eigen::VectorXd computeTorques(const Eigen::VectorXd& qddot);   // line 143
```

**Target pattern:**
```cpp
Eigen::VectorXd computeTorques(const Eigen::VectorXd& qddot,
                               const Vector3d& gravity = Vector3d::Zero());
```

---

### `src/InverseDynamics.cpp` (service, transform) — MODIFIED

**Analog:** `src/InverseDynamics.cpp` (current file)

**Current outwardPass pattern to MODIFY (add gravity per D-08):**
```cpp
void InverseDynamics::outwardPass()
{
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        int parent = links[i].parent;
        if (parent == -1)
        {
            links[i].v = links[i].S * links[i].qdot;
            links[i].a = links[i].S * links[i].qddot;
        }
        else { /* ... */ }
    }
}
```

**Target outwardPass (with gravity as base acceleration):**
```cpp
void InverseDynamics::outwardPass()
{
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        int parent = links[i].parent;
        if (parent == -1)
        {
            links[i].v = links[i].S * links[i].qdot;
            // Gravity: base acceleration = S₀·q̈₀ - g
            links[i].a = links[i].S * links[i].qddot
                       - MotionVector(Vector3d::Zero(), gravity);
        }
        else { /* ... */ }
    }
}
```

**Current computeTorques pattern (lines 103-140) — add gravity member:**
```cpp
Eigen::VectorXd InverseDynamics::computeTorques(const Eigen::VectorXd& qddot)
{
    // Validate input
    // ...
    outwardPass();
    Eigen::VectorXd tau = inwardPass();
    return tau;
}
```

---

### `include/SpatialVector.h` (model, CRUD) — MODIFIED

**Analog:** `include/SpatialVector.h` (current file, lines 44, 90)

**Current global using Vector3d pattern to FIX (D-23):** Line 44: `using Vector3d = Eigen::Matrix<double, 3, 1>;` at global scope.

**Move into namespace SpatialAlgebra:**
```cpp
namespace SpatialAlgebra {
using Vector3d = Eigen::Matrix<double, 3, 1>;
// ... class definitions ...
}
```

**Constructor to add NaN assertions (D-13):**
```cpp
SpatialVector(const Vector3d &angular, const Vector3d &linear);
```
Add debug-mode assertions:
```cpp
SpatialVector(const Vector3d &angular, const Vector3d &linear)
    : angular(angular), linear(linear)
{
#ifndef NDEBUG
    if (angular.hasNaN() || linear.hasNaN() ||
        angular.array().isInf().any() || linear.array().isInf().any()) {
        std::cerr << "WARNING: NaN or Inf detected in SpatialVector constructor\n";
    }
#endif
}
```

**Existing NaN assertion pattern reference (LowerTriangular.h, lines 148-153):**
```cpp
#ifndef NDEBUG
if (i >= n || j >= n || i < 0 || j < 0)
    throw std::out_of_range("Index out of bounds");
#endif
```

---

### `include/RigidBodyInertia.h` (model, CRUD) — MODIFIED

**Analog:** `include/RigidBodyInertia.h` (current file)

**apply() method to add NaN assertions (D-13):** Line 88-101.

**Target pattern:**
```cpp
inline ForceVector apply(const MotionVector &mv) const
{
#ifndef NDEBUG
    if (mv.getAngular().hasNaN() || mv.getLinear().hasNaN())
        std::cerr << "WARNING: NaN detected in RigidBodyInertia::apply\n";
#endif
    // ... existing implementation ...
}
```

**Type alias location (D-23):** Line 18 `using Vector6d = ...` and line 19 `using lt = LowerTriangular;` already inside namespace. The `using Vector3d` from SpatialVector.h was leaking — now it moves inside namespace.

---

### `include/PluckerTransform.h` (model, transform) — MODIFIED

**Analog:** `include/PluckerTransform.h` (current file, lines 179-190)

**Current auto return type to FIX (D-21):**
```cpp
auto apply(const fv& v) const {
    return transformForce(v);
}

auto apply(const mv& v) const {
    return transformMotion(v);
}
```

**Target pattern (explicit return type — or move body after return type is visible):**
```cpp
ForceVector apply(const fv& v) const {
    return transformForce(v);
}

MotionVector apply(const mv& v) const {
    return transformMotion(v);
}
```

**multiply() return type issue (D-21):** The `multiply()` method at line 172 already returns `PluckerTransform` explicitly. The issue may be in the `.cpp` file's implementation. See `src/PluckerTransform.cpp:227-238` which already returns `PluckerTransform` correctly, but may need the signature aligned.

---

### `src/PluckerTransform.cpp` (service, transform) — MODIFIED

**Analog:** `src/PluckerTransform.cpp` (current file)

**multiply() implementation (lines 227-238) — already returns PluckerTransform explicitly:**
```cpp
PluckerTransform PluckerTransform::multiply(const PluckerTransform &X) const
{
    Rotation newRotation = rotation * X.rotation;
    Vector3d newTranslation = X.translation + static_cast<const Eigen::Matrix3d &>(X.rotation.transpose()) * translation;
    return PluckerTransform(newRotation, newTranslation);
}
```
No change needed; return type is already explicit.

**Apply methods (lines 240-243) — already delegates to multiply:**
```cpp
PluckerTransform PluckerTransform::apply(const PluckerTransform &X) const
{
    return multiply(X);
}
```

---

### `include/LowerTriangular.h` (model, CRUD) — MODIFIED

**Analog:** `include/LowerTriangular.h` (current file)

**Include guard: Already `#ifndef LOWER_TRIANGULAR_H` at line 1.** D-24 may already be done. Verify on disk.

**OpenMP pragma to remove (D-26):** Line 202: `#pragma omp parallel for collapse(2)`

**Target (remove the pragma line):**
```cpp
LowerTriangular operator*(const LowerTriangular &other) const
{
    // ... validation ...
    LowerTriangular result(n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            // ...
        }
    }
    return result;
}
```

**Global namespace using declaration to fix (line 569):** `using SpatialAlgebra::LowerTriangular;` — move inside namespace or remove.

---

### `src/LowerTriangular.cpp` (service, CRUD) — MODIFIED

**Analog:** `src/LowerTriangular.cpp` (current file, lines 28-53)

**Current pattern (no OMP pragma in .cpp — it's in the header):** No changes needed in .cpp for D-26. The OMP pragma is in the header at line 202.

---

### `include/SpatialOperations.h` (utility, transform) — MODIFIED

**Analog:** `include/SpatialOperations.h` (current file, lines 19-21)

**Current unsafe downcast signatures to FIX (D-22):**
```cpp
static SpatialVector crossProductMotion(const MotionVector& v1, const MotionVector& v2);
static SpatialVector crossProductForce(const MotionVector& v, const ForceVector& f);
```

These already accept the correct types (not using `SpatialVector` with static_cast). The implementation in `src/SpatialOperations.cpp` already delegates to the free functions correctly.

**No change needed** — the signatures already accept `MotionVector&` and `ForceVector&` directly.

---

### `src/SpatialOperations.cpp` (utility, transform) — MODIFIED

**Analog:** `src/SpatialOperations.cpp` (current file, lines 11-17)

**Current pattern (already delegates correctly):**
```cpp
SpatialVector SpatialOperations::crossProductMotion(const MotionVector& v1, const MotionVector& v2) {
    return cross(v1, v2);
}

SpatialVector SpatialOperations::crossProductForce(const MotionVector& v, const ForceVector& f) {
    return cross(v, f);
}
```

No change needed — matches D-22 pattern.

---

### `tests/TestSpatialOperations.cpp` (test, CRUD) — MODIFIED

**Analog:** `tests/TestSpatialOperations.cpp` (current file, lines 27-37)

**Helper function pattern to FIX (D-15 — needs to return actual matrices):**
```cpp
LowerTriangular createIdentityInertia() {           // lines 27-31
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    return LowerTriangular::fromFullMatrix(identity);
}

LowerTriangular createDiagonalInertia(double value) { // lines 33-37
    Eigen::Matrix3d diagonal = Eigen::Matrix3d::Identity() * value;
    return LowerTriangular::fromFullMatrix(diagonal);
}
```

**Test pattern (GTest with property-based tests, lines 56-114):**
```cpp
TEST(TestCrossProductMotion, SimpleRotationVectors) {
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);
    EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
    // ...
}
```

**GTest struct pattern (lines 46-50):**
```cpp
class TestCrossProductMotion : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};
```

**Main entry pattern (lines 316-320):**
```cpp
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

---

### `tests/TestForwardDynamics.cpp` (test, transform) — MODIFIED

**Analog:** `tests/TestForwardDynamics.cpp` (current file)

**Single-link test pattern (lines 19-43):**
```cpp
TEST(ForwardDynamicsTest, SingleLinkPendulum) {
    ForwardDynamics fd;
    Link link;
    link.parent = -1;
    link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
    link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    link.q = 0.0;
    link.qdot = 0.0;
    link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(link);
    
    Eigen::VectorXd tau(1);
    tau[0] = 1.0;
    fd.computeAccelerations(tau);
    
    EXPECT_NEAR(fd.links[0].qddot, 1.0, EPSILON);
}
```

**Multi-link test pattern (lines 51-92) — to be extended with known numerical values (D-17):**
```cpp
TEST(ForwardDynamicsTest, TwoLinkSerialChain) {
    ForwardDynamics fd;
    // Link 0 setup ...
    // Link 1 setup ...
    Eigen::VectorXd tau(2);
    fd.computeAccelerations(tau);
    EXPECT_GT(fd.links[0].qddot, 0.0);
    EXPECT_GT(fd.links[1].qddot, 0.0);
}
```

---

### `tests/TestInverseDynamics.cpp` (test, transform) — MODIFIED

**Analog:** `tests/TestInverseDynamics.cpp` (current file)

**Test pattern to extend for non-zero velocity (D-16):**
```cpp
TEST(InverseDynamicsTest, TwoLinkSerialChain) {
    InverseDynamics id;
    // Link 0 setup ...
    // Link 1 setup ...
    Eigen::VectorXd qddot(2);
    qddot[0] = 1.0;
    qddot[1] = 0.5;
    Eigen::VectorXd tau = id.computeTorques(qddot);
    EXPECT_GT(std::abs(tau[0]), 0.0);
}
```

**New tests needed: Set non-zero qdot and verify tau includes Coriolis terms.**

---

### `tests/TestSpatialVector.cpp` (test, CRUD) — MODIFIED

**Analog:** `tests/TestSpatialVector.cpp` (current file)

**Cross-force anti-commutativity test (lines 554-569) — verify with mixed inputs (D-03):**
```cpp
TEST(TestForceVector, CrossForceAntiCommutativity) {
    ForceVector a(Vector3d(1, 2, 3), Vector3d(4, 5, 6));
    ForceVector b(Vector3d(2, 3, 4), Vector3d(5, 6, 7));
    ForceVector a_cross_b = a.crossForce(b);
    ForceVector b_cross_a = b.crossForce(a);
    ForceVector neg_b_cross_a = b_cross_a * -1.0;
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[0], neg_b_cross_a.getAngular()[0]);
    // ...
}
```

**New test needed (mixed torque+force inputs — D-03):**
```cpp
TEST(TestForceVector, MixedTorqueForceCross) {
    // τ=(1,0,0), f=(0,1,0) — mixed torque-force input exercises all 3 terms
    ForceVector a(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    ForceVector b(Vector3d(0, 1, 0), Vector3d(0, 0, 1));
    ForceVector result = cross(a, b);
    // τ1×τ2 = (1,0,0)×(0,1,0) = (0,0,1)
    // f1×f2 = (0,1,0)×(0,0,1) = (1,0,0)
    // angular = (1,0,1)
    // τ1×f2 = (1,0,0)×(0,0,1) = (0,-1,0)
    // τ2×f1 = (0,1,0)×(0,1,0) = (0,0,0)
    // linear = (0,-1,0)
    EXPECT_NEAR(result.getAngular()[0], 1.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[1], 0.0, EPSILON);
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[0], 0.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[1], -1.0, EPSILON);
    EXPECT_NEAR(result.getLinear()[2], 0.0, EPSILON);
}
```

---

### `CMakeLists.txt` (config, build) — MODIFIED

**Analog:** `CMakeLists.txt` (current file)

**Current Eigen version pin to remove (D-18):** Line 12: `find_package(Eigen3 REQUIRED NO_MODULE)` — already unpinned.

**Current GTest pattern to MODIFY (D-25 — add FetchContent fallback):**
```cmake
find_package(GTest REQUIRED)      # line 15
```

**Target pattern with FetchContent fallback:**
```cmake
# Find GTest with FetchContent fallback
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/release-1.12.1.zip
    )
    FetchContent_MakeAvailable(googletest)
endif()
```

**Current stub source inclusion (D-27):** Line 18: `file(GLOB SOURCES "src/*.cpp")` — Remove stub files or change glob.

**Target:**
```cmake
# Source files (excluding empty stubs)
file(GLOB SOURCES "src/*.cpp")
# Remove empty stubs from the list
list(REMOVE_ITEM SOURCES 
    "${CMAKE_SOURCE_DIR}/src/RigidBodyInertia.cpp"
    "${CMAKE_SOURCE_DIR}/src/ArticulatedBodyInertia.cpp"
)
```

**Test target pattern (lines 30-34):**
```cmake
add_executable(TestForwardDynamics tests/TestForwardDynamics.cpp)
target_link_libraries(TestForwardDynamics SpatialAlgebra GTest::GTest GTest::Main)
add_test(NAME TestForwardDynamics COMMAND TestForwardDynamics)
```

---

## Shared Patterns

### Include Guard Convention
**Source:** All header files
**Apply to:** `include/SpatialAlgebra.h` (new), `include/LowerTriangular.h` (verify)
**Pattern:**
```cpp
#ifndef SPATIAL_ALGEBRA_H    // Upper-case, underscores, matches filename
#define SPATIAL_ALGEBRA_H
// ...
#endif // SPATIAL_ALGEBRA_H
```

### Doxygen Documentation Convention
**Source:** All header files
**Apply to:** All modified files — add/update doc blocks for changed signatures
**Pattern (from `include/ForwardDynamics.h:139-155`):**
```cpp
/**
 * @brief Compute joint accelerations from applied torques
 * @param tau Vector of joint torques
 * @param gravity Optional gravity vector (defaults to zero)
 * @details Main entry point for forward dynamics computation.
 * ...
 * @throws std::invalid_argument if tau.size() != links.size()
 */
void computeAccelerations(const Eigen::VectorXd& tau);
```

### Debug-mode NaN/Inf Assertion Pattern
**Source:** `include/LowerTriangular.h:148-153`
**Apply to:** `include/SpatialVector.h`, `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`, `include/SpatialUtils.h`
**Pattern:**
```cpp
#ifndef NDEBUG
if (angular.hasNaN() || linear.hasNaN() || 
    angular.array().isInf().any() || linear.array().isInf().any()) {
    std::cerr << "WARNING: NaN or Inf detected in SpatialVector constructor\n";
}
#endif
```

### Exception-based Error Handling
**Source:** All dynamics `.cpp` files
**Apply to:** All modified service files
**Pattern (from `src/ForwardDynamics.cpp:109-115`):**
```cpp
if (std::abs(denom) < EPSILON) {
    throw std::runtime_error(
        "ForwardDynamics::inwardPass: Near-zero inertia at joint " + 
        std::to_string(i) + " (denom=" + std::to_string(denom) + ")"
    );
}
```

### GTest Test Suite Pattern
**Source:** `tests/TestPluckerTransform.cpp` (most comprehensive example)
**Apply to:** All modified test files
**Pattern:**
```cpp
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
constexpr double EPSILON = 1e-10;

TEST(TestGroup, TestName) {
    // Arrange
    // Act
    // Assert
    EXPECT_NEAR(actual, expected, EPSILON);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

### Type Alias Convention
**Source:** `include/RigidBodyInertia.h:112`, `include/ForceVector.h:160`, `include/MotionVector.h:157`, `include/ArticulatedBodyInertia.h:210`, `include/PluckerTransform.h:207`
**Apply to:** All header files (existing pattern to preserve)
**Pattern:**
```cpp
// At end of namespace, after class definition:
using rbi = RigidBodyInertia;    // lowercase abbreviation
using fv = ForceVector;
using mv = MotionVector;
using abi = ArticulatedBodyInertia;
using plux = PluckerTransform;
using lt = LowerTriangular;
```

### Cross-Product Delegation Pattern
**Source:** RESEARCH.md (D-01 canonical approach)
**Apply to:** `src/SpatialVector.cpp:48-54`, `src/ForceVector.cpp:40-46`
**Pattern:**
```cpp
SpatialVector SpatialVector::crossForce(const SpatialVector &other) const {
    // Delegate to canonical free function in SpatialUtils.h
    return ForceVector(cross(ForceVector(*this), ForceVector(other)));
}
```

### Two-Phase ABA Inward Pass Pattern
**Source:** RESEARCH.md (D-05) and `src/ForwardDynamics.cpp:53-126` (already implemented)
**Apply to:** `src/ForwardDynamics.cpp` (preserve existing correct pattern)
**Pattern:**
```cpp
// Phase 1: Initialize Ia and pa from rigid body inertia (all links)
for (int i = 0; i < links.size(); i++) {
    links[i].Ia = ArticulatedBodyInertia(...);
    links[i].pa = ...;
}
// Phase 2: Accumulate child contributions (tip→base, no re-init)
for (int i = links.size()-1; i >= 0; i--) {
    // ... accumulate into parent without re-initializing
}
// Phase 3: Solve for joint accelerations
for (int i = links.size()-1; i >= 0; i--) {
    // ... compute qddot
}
```

## No Analog Found

Files with no close match in the codebase (planner should use RESEARCH.md patterns instead):

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| `.github/workflows/ci.yml` | config | event-driven | No existing CI infrastructure — first workflow in project |
| `include/SpatialAlgebra.h` | config | request-response | No existing umbrella header — all headers included individually |

## Metadata

**Analog search scope:** `include/`, `src/`, `tests/`, `CMakeLists.txt`
**Files scanned:** 35 (all source, header, and test files)
**Pattern extraction date:** 2026-05-17
