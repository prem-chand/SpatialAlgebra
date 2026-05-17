# SpatialAlgebra: Mathematical Conventions

## 1. Coordinate Frame Convention

The library uses a **right-handed Cartesian coordinate system** throughout.
All frame transforms follow the standard robotics convention:
the transform `X` from frame A to frame B satisfies `v_B = X · v_A`.

The gravity vector **g** is the gravitational acceleration vector expressed
in **world coordinates**. On Earth, the standard value is:

    g = [0, 0, -9.81]  (m/s²)

The magnitude `||g|| = 9.80665 m/s²` is used as the standard Earth surface
value in all test oracles.

## 2. Spatial Vector Representation

A spatial vector **s** is a 6D vector composed of two 3D vectors stacked
vertically:

    s = [angular; linear]

The interpretation of the two components depends on the vector type:

| Type            | Angular component   | Linear component |
|-----------------|---------------------|------------------|
| MotionVector    | ω (angular velocity) | v (linear velocity) |
| ForceVector     | τ (torque/moment)   | f (linear force) |

This ordering (angular first, linear second) follows Featherstone (2008)
§2.1.

## 3. Cross-Product Formulas

All cross products follow the Featherstone formulation (Chapter 2) with
the cross product operator × acting on spatial vectors. Four combinations
are defined via free functions in `SpatialUtils.h`:

### 3.1 Motion × Motion (crm) — Result: MotionVector

    [ω₁; v₁] × [ω₂; v₂] = [ ω₁ × ω₂ ; ω₁ × v₂ + v₁ × ω₂ ]

Used for velocity- and acceleration-related products in both ABA outward
pass (bias acceleration: `c = v × S · q̇`) and RNEA outward pass (Coriolis
acceleration: `v × S · q̇`).

### 3.2 Motion × Force (crf) — Result: ForceVector

    [ω₁; v₁] × [τ₂; f₂] = [ ω₁ × τ₂ + v₁ × f₂ ; ω₁ × f₂ ]

Used in ABA inward pass for Coriolis/centrifugal bias force computation:
`pₐ = Iₐ·c + v × Iₐ·v`, and in RNEA inward pass for Coriolis force:
`I·a + v × I·v`.

### 3.3 Force × Motion — Result: ForceVector

    [τ₁; f₁] × [ω₂; v₂] = [ τ₁ × ω₂ + f₁ × v₂ ; τ₁ × v₂ + f₁ × ω₂ ]

This is the dual of motion × force, differing in sign of the v×τ cross
term relative to motion × motion products.

### 3.4 Force × Force — Result: ForceVector

    [τ₁; f₁] × [τ₂; f₂] = [ τ₁ × τ₂ + f₁ × f₂ ; τ₁ × f₂ - τ₂ × f₁ ]

This cross product is anti-commutative: `f₁ × f₂ = -(f₂ × f₁)`.
It appears in force transformation operations and is the canonical
single-source formula for testing cross-product correctness.

## 4. Gravity Convention

Gravity enters both dynamics algorithms through the base link acceleration
term, following Featherstone §7.4 (D-07, D-08).

### 4.1 ABA (Forward Dynamics)

The base link bias acceleration is set to the **negative** of the gravity
vector:

    c₀ = -g

This is implemented in `ForwardDynamics::outwardPass()`:

```cpp
links[i].c = MotionVector(Vector3d::Zero(), -this->gravity);
```

For the default `gravity = Vector3d::Zero()`, the bias acceleration is zero
(corresponding to free-floating base).

### 4.2 RNEA (Inverse Dynamics)

The base link spatial acceleration is constructed as:

    a₀ = S₀·q̈₀ - [0; g]

where `[0; g]` is a MotionVector with zero angular component and the gravity
vector as the linear component. This is implemented in
`InverseDynamics::outwardPass()`:

```cpp
links[i].a = links[i].S * links[i].qddot
           - MotionVector(Vector3d::Zero(), gravity);
```

### 4.3 Forward Propagation

For non-base links, gravity propagates through the kinematic chain as part
of the recursive acceleration computation:

**RNEA outward pass:**
```
aᵢ = Xᵢ · a_parent + Sᵢ · q̈ᵢ + vᵢ × Sᵢ · q̇ᵢ
```

**ABA outward pass:**
```
cᵢ = Xᵢ · c_parent + vᵢ × Sᵢ · q̇ᵢ
```

This ensures that the gravity effect is consistently transmitted from base
to tip across both algorithms.

## 5. Plücker Transform Conventions

### 5.1 Transform Direction

The Plücker transform `X` stored in each link represents the transform
**from the parent frame to the child frame**. The motion transform
`X.transformMotion(v)` computes:

    v_child = X · v_parent

where `v_parent` is a motion vector expressed in the parent frame.

### 5.2 Force Transform

Force vectors transform contravariantly with respect to motion vectors:

    f_child = X.transformForce(f_parent) = X^{-T} · f_parent

This is computed by `PluckerTransform::transformForce()`.

### 5.3 Inverse Force Transform (Inward Pass)

In the inward pass of both algorithms, child forces are propagated back
to the parent frame using the inverse force transform:

    f_parent += X_child.inverseTransformForce(f_child)

where `X_child` is the Plücker transform from parent to child.

### 5.4 Articulated Body Inertia Transform

In the ABA inward pass, articulated body inertias are transformed using
`X.invtformABI(Ia)`, which computes:

    Ia_parent = X^{-1} · Ia_child · X^{-T}

This is the congruence transform that preserves the properties of the
articulated body inertia under coordinate changes.

## 6. Tolerance Strategy

Numerical tolerances are selected based on the complexity of the computation:

| Context         | EPSILON | Rationale                                  |
|-----------------|---------|--------------------------------------------|
| Unit tests      | 1e-10   | Direct algebraic operations (scalar ops)   |
| Round-trip dyn. | 1e-8    | Multi-step dynamics (accumulated FP error) |
| Multi-link dyn. | 1e-8    | Chain of 3+ transforms and force propag.   |

- **1e-10** is used in `TestInverseDynamics.cpp`, `TestForwardDynamics.cpp`,
  and `TestSpatialVector.cpp` for single-step verifications.
- **1e-8** is used in `TestDynamicsConsistency.cpp` for round-trip tests
  that involve both RNEA and ABA, where floating-point error accumulates
  across multiple O(n) passes.

The choice of 1e-10 for unit tests is justified because:
- All operations are double-precision (64-bit)
- Single-link scalar operations involve at most ~50 FLOPs
- Eigen3 uses fused multiply-add (FMA) instructions on Apple Silicon,
  reducing accumulation error

The choice of 1e-8 for round-trip tests is justified because:
- A full RNEA + ABA round-trip involves ~20n operations per link
- Double-precision machine epsilon is ~2.2e-16
- Accumulated error for 10-link chain is ~O(√n · ε) ≈ 7e-16 per operation,
  well within 1e-8 tolerance

## 7. API Change Log

This section documents significant API changes made during development.

| Date       | Change                                     | Rationale                             |
|------------|--------------------------------------------|---------------------------------------|
| 2026-05-14 | Removed `MotionVector::crossForce`         | Duplicated free-function `cross()`    |
| 2026-05-14 | Removed `ForceVector::crossMotion`         | Duplicated free-function `cross()`    |
| 2026-05-15 | Moved `Vector3d` typedef into `SpatialAlgebra` namespace | Avoid ODR violations with Eigen |
| 2026-05-15 | Created `include/SpatialAlgebra.h` umbrella header       | Single-include entry point     |
| 2026-05-16 | Changed `PluckerTransform::apply(const fv&)` return type to `ForceVector` | Fixed type-punning ODR hazard |
| 2026-05-16 | Changed `PluckerTransform::apply(const mv&)` return type to `MotionVector` | Fixed type-punning ODR hazard |
| 2026-05-16 | Removed unused `crossForce`/`crossMotion` overloads from `MotionVector.h`/`ForceVector.h` | Eliminated code duplication |
| 2026-05-17 | Restructured ABA inward pass into explicit Phase 1 (init) and Phase 2 (accumulate) | Fixed CR-02: child Ia overwriting parent |
| 2026-05-17 | Changed `RigidBodyInertia::apply()` to use `MotionVector` param and `ForceVector` return | Type-correct spatial inertia application |

---

*Maintained as part of the SpatialAlgebra project.*  
*Last updated: 2026-05-17*
