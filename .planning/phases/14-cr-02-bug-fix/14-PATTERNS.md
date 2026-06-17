# Phase 14: CR-02 Bug Fix - Pattern Map

**Mapped:** 2026-06-17
**Files analyzed:** 4 (modified) / 0 (new)
**Analogs found:** 4 / 4

## File Classification

| Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---------------|------|-----------|----------------|---------------|
| `src/ForwardDynamics.cpp` | service (algorithm) | transform-accumulate | `src/ForwardDynamics.cpp:82-179` (existing `inwardPass`) | exact (same file, restructured) |
| `include/ForwardDynamics.h` | header / config | — | `include/ForwardDynamics.h:5-16,117-191` (existing Doxygen) | exact (same file, docs-only) |
| `tests/TestForwardDynamics.cpp` | test | request-response | `tests/TestDynamicsConsistency.cpp:100-177` (`ThreeLinkSerialChain` round-trip) | role-match (cross-file test pattern) |
| `tests/TestDynamicsConsistency.cpp` | test | request-response | `tests/TestDynamicsConsistency.cpp:100-177` (itself — existing 3-link round-trip) | exact (same file, variant pattern) |

## Pattern Assignments

---

### `src/ForwardDynamics.cpp` — Restructure `inwardPass()` (service, transform-accumulate)

**Analog:** `src/ForwardDynamics.cpp` lines 82-179 (existing `inwardPass` — same file, algorithmic restructure)

**Imports pattern** (lines 1-16):
```cpp
#include "ForwardDynamics.h"
#include "LowerTriangular.h"
#include "SpatialUtils.h"
#include <stdexcept>
```

No import changes needed — all required includes already present. `LowerTriangular.h` provides `lt::fromFullMatrix()`; `SpatialUtils.h` provides `dot()`, `cross()`, `skew()`.

**Anonymous namespace helper — KEEP AS-IS** (lines 20-55):
```cpp
namespace {
    ArticulatedBodyInertia transformInertiaToParent(const PluckerTransform &X,
                                                    const ArticulatedBodyInertia &Ia) {
        // KEEP — computes X^T * Ia * X, correct for child-to-parent ABA propagation
        // [lines 22-53 unchanged]
    }
}
```

This helper is preserved verbatim. It converts ABI to 6×6, applies `X_6x6.transpose() * Ia_6x6 * X_6x6`, and symmetrizes `I_new`/`M_new`. Used in the condensation pass-to-parent step.

**Core pattern — single tip-to-base sweep** (replaces lines 82-179):

The existing 3-phase design:
1. Phase 1 (lines 84-102): Initialize Ia/pa for ALL links (base→tip)
2. Phase 2 (lines 104-157): Backward pass — condense Ia/pa, propagate to parent (tip→base)
3. Phase 3 (lines 161-179): Forward pass — correct qddot, compute spatial accelerations (base→tip)

Must be replaced with a single sweep (tip→base): accumulate → solve qddot → condense → pass to parent. Then compute spatial accelerations in a final forward pass.

**Phase 1 pattern (init all Ia/pa — KEEP, lines 84-102):**
```cpp
// Phase 1: Initialize Ia and pa from rigid body inertia
for (int i = 0; i < static_cast<int>(links.size()); i++)
{
    double mass = links[i].I.getMass();
    Vector3d com = links[i].I.getCom();
    links[i].Ia = ArticulatedBodyInertia(
        links[i].I.getInertiaMatrixLT(),
        skew(com) * mass,
        lt::Identity(3) * mass
    );

    ForceVector IaV = links[i].Ia.apply(links[i].v);
    links[i].pa = links[i].Ia.apply(links[i].c) + cross(links[i].v, IaV);

    links[i].pa = ForceVector(
        links[i].pa.getAngular() + links[i].f.getAngular(),
        links[i].pa.getLinear() + links[i].f.getLinear()
    );
}
```

This initialization is correct and must be kept. It sets up Ia from rigid body inertia and pa from bias forces + external forces for ALL links before the tip-to-base sweep begins.

**Phase 2 pattern — single tip-to-base sweep with condensation** (replaces lines 104-179):

For each link (tip→base), in order:
1. Children's condensed Ia/pa already accumulated from previous iterations
2. Compute `IaS = Ia.apply(S)`, `D = dot(S, IaS)`, `u = tau[i] - dot(S, pa)`
3. Solve `qddot = u / D`
4. If parent exists: condense Ia/pa, then pass to parent
5. Skip condensation for base link (parent == -1)

```cpp
// Phase 2: Single tip-to-base sweep — accumulate, solve, condense, pass
constexpr double EPSILON = 1e-10;

for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
{
    int parent = links[i].parent;

    // --- Children contributions already accumulated ---
    // (Children processed in previous iterations added their
    //  condensed Ia/pa to this link's Ia/pa)

    // --- Solve for qddot ---
    ForceVector IaS = links[i].Ia.apply(links[i].S);
    double D = dot(links[i].S, IaS);

    if (std::abs(D) < EPSILON)
    {
        throw std::runtime_error(
            "ForwardDynamics::inwardPass: Near-zero inertia at joint "
            + std::to_string(i) + " (D=" + std::to_string(D) + ")"
        );
    }

    double u = tau[i] - dot(links[i].S, links[i].pa);
    links[i].qddot = u / D;

    // --- Condense for parent ---
    if (parent != -1)
    {
        double invD = 1.0 / D;

        // Ia' = Ia - (Ia*S)*(1/D)*(Ia*S)^T
        Vector3d t = IaS.getAngular();
        Vector3d f = IaS.getLinear();

        lt inertiaCorr = lt::fromFullMatrix(t * t.transpose()) * invD;
        Eigen::Matrix3d HCorr = t * f.transpose() * invD;
        lt massCorr = lt::fromFullMatrix(f * f.transpose()) * invD;

        links[i].Ia = links[i].Ia
            + (ArticulatedBodyInertia(inertiaCorr, HCorr, massCorr) * (-1.0));

        // pa' = pa + Ia*S * qddot
        links[i].pa = links[i].pa + (IaS * links[i].qddot);

        // ---- Pass to parent ----
        links[parent].Ia = links[parent].Ia
            + transformInertiaToParent(links[i].X, links[i].Ia);
        links[parent].pa = ForceVector(
            links[parent].pa.getAngular()
                + links[i].X.inverseTransformForce(links[i].pa).getAngular(),
            links[parent].pa.getLinear()
                + links[i].X.inverseTransformForce(links[i].pa).getLinear()
        );
    }
}
```

**Phase 3 pattern — forward spatial acceleration pass** (replaces old Phase 3, lines 161-179):

After the tip-to-base sweep, compute spatial accelerations in a base→tip pass. No correction term needed — `qddot[i]` is already the final answer.

```cpp
// Phase 3: Forward pass — compute spatial accelerations
std::vector<MotionVector> a(links.size());
for (int i = 0; i < static_cast<int>(links.size()); i++)
{
    int parent = links[i].parent;

    MotionVector aParentInChild = (parent == -1)
        ? MotionVector(Vector3d::Zero(), Vector3d::Zero())
        : links[i].X.transformMotion(a[parent]);

    a[i] = aParentInChild + links[i].c + links[i].S * links[i].qddot;
}
```

**Error handling pattern** (from lines 118-124, 105):
```cpp
constexpr double EPSILON = 1e-10;

if (std::abs(D) < EPSILON)
{
    throw std::runtime_error(
        "ForwardDynamics::inwardPass: Near-zero inertia at joint "
        + std::to_string(i) + " (D=" + std::to_string(D) + ")"
    );
}
```

**Validation pattern** (from `computeAccelerations`, lines 182-209 — KEEP AS-IS):
```cpp
void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity)
{
    if (tau.size() != static_cast<int>(links.size())) { /* throw */ }
    for (int i = 0; i < tau.size(); i++) {
        if (std::isnan(tau[i]) || std::isinf(tau[i])) { /* throw */ }
    }
    this->gravity = gravity;
    outwardPass();
    inwardPass(tau);
}
```

**Key changes from existing code:**
- REMOVED: `std::vector<double> D_store(links.size())` — no longer needed (no correction phase)
- REMOVED: `std::vector<ArticulatedBodyInertia> Ia_unc(links.size())` — no longer needed
- REMOVED: Old Phase 3 correction loop (lines 161-179) — qddot is final after Phase 2 sweep
- CHANGED: `qddot` is now the FINAL answer (not "partial" requiring correction)
- CHANGED: Phase 3 becomes pure spatial acceleration computation (base→tip)
- KEPT: `outwardPass()` (lines 57-80) — no changes, gravity handling is correct
- KEPT: `transformInertiaToParent` helper (lines 22-53) — correct `X^T * Ia * X` transform
- KEPT: `computeAccelerations()` validation (lines 182-209) — no changes

**Critical ordering invariant:**
1. Accumulate children BEFORE solving qddot
2. Solve qddot BEFORE condensing
3. Condense BEFORE passing to parent
4. Skip condensation for base link (parent == -1)

---

### `include/ForwardDynamics.h` — Doxygen Updates (header, documentation)

**Analog:** `include/ForwardDynamics.h` lines 5-16, 117-191 (existing Doxygen on same file)

**File-level Doxygen pattern** (lines 5-16 — update description):
```cpp
/**
 * @file ForwardDynamics.h
 * @brief Articulated Body Algorithm (ABA) for forward dynamics computation
 * @details This file implements Featherstone's Articulated Body Algorithm (ABA)
 *          for computing joint accelerations from applied torques. The ABA is
 *          an O(n) recursive algorithm that efficiently solves the forward
 *          dynamics problem for serial and branching kinematic chains.
 * 
 *          Algorithm Overview:
 *          1. Outward pass (base to tip): Propagate velocities, compute bias accelerations
 *          2. Inward pass (tip to base): Accumulate articulated inertias,
 *             solve for qddot with condensation, pass to parent
 *          3. Forward pass (base to tip): Compute spatial accelerations from qddot
 * 
 *          Mathematical Foundation:
 *          TODO: update with condensation formula citation
```

Update the Algorithm Overview to reflect the single sweep: say "solve for qddot with condensation" instead of "solve for accelerations."

**Class-level Doxygen pattern** (lines 117-152 — update description):
```cpp
/**
 * @brief Forward dynamics solver using Articulated Body Algorithm
 * @details The ForwardDynamics class implements Featherstone's ABA (Algorithm 7.3)
 *          for efficient O(n) computation of joint accelerations. The algorithm
 *          handles both serial chains and branching kinematic trees.
 * 
 *          Algorithm Complexity:
 *          - Time: O(n) where n is the number of links
 *          - Space: O(n) for storing link states
 * 
 *          Usage Pattern:
 *          1. Setup kinematic tree by populating links vector
 *          2. Set joint states (q, qdot) and external forces (f)
 *          3. Call computeAccelerations(tau) to solve for qddot
 * 
 * @note Links must be ordered such that parents appear before children
 * @note For branching trees, all children of a link must have consecutive indices
 */
```

**`computeAccelerations` Doxygen pattern** (lines 139-156 — update description):
```cpp
/**
 * @brief Compute joint accelerations from applied torques
 * @param tau Vector of joint torques (must match links.size())
 * @details Main entry point for forward dynamics computation.
 *          Executes the three-pass ABA:
 *          1. Outward pass: propagate velocities, compute bias accelerations
 *          2. Inward pass (single sweep): accumulate articulated inertias,
 *             solve qddot with condensation, propagate to parent
 *          3. Forward pass: compute spatial accelerations from qddot
 * 
 *          Mathematical formulation:
 *          q̈_i = (τ_i - S_iᵀ·p_A_i) / (S_iᵀ·I_A_i·S_i)
 *          where I_A_i is the accumulated articulated inertia at joint i
 *          and p_A_i is the accumulated bias force.
 *          The condensation step removes the joint i DOF from I_A_i
 *          before passing to the parent: I_A' = I_A - I_A·S·(Sᵀ·I_A·S)⁻¹·Sᵀ·I_A
 * 
 * @param gravity Gravity vector (default zero) for gravity-aware dynamics
 * @throws std::invalid_argument if tau.size() != links.size()
 * @throws std::runtime_error if denominator is near zero (singular configuration)
 */
```

**`inwardPass` Doxygen pattern** (lines 175-190 — update description):
```cpp
/**
 * @brief Inward pass: accumulate articulated inertias, solve qddot with condensation
 * @param tau Vector of joint torques
 * @details Single tip-to-base sweep. For each link (tip to base):
 *          - Children's condensed Ia/pa already accumulated from prior iterations
 *          - Compute IaS = Ia·S, D = Sᵀ·IaS, u = τ - Sᵀ·pa
 *          - Solve: qddot = u / D
 *          - Condense: I_A' = I_A - (I_A·S)·D⁻¹·(I_A·S)ᵀ
 *          - Condense: p_A' = p_A + I_A·S·qddot
 *          - Pass condensed I_A'/p_A' to parent (if parent exists)
 * 
 *          Condensation (per Featherstone Algorithm 7.3) removes the joint
 *          motion subspace before propagating to the parent, ensuring the
 *          parent sees the correct reflected inertia of the subtree.
 * 
 *          After the sweep, a forward pass computes spatial accelerations
 *          from the final qddot values and bias accelerations.
 * 
 * @note Base link (parent == -1) skips condensation — no parent to pass to
 * @note Children must be processed before their parent (tip-to-base order)
 */
```

**Doxygen style conventions (from existing code):**
- Block comments `/** */` on every declaration
- `@brief` on first line, `@details` for extended description
- `@param`, `@return`, `@throws`, `@note`, `@see` tags
- Mathematical notation: Unicode subscripts (`q̈`, `q̇`, `Sᵀ`), Greek letters transliterated
- References: `@see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7`
- `@code{.cpp} ... @endcode` for examples

**No changes to:**
- `Link` struct (lines 79-111) — fields unchanged
- `ForwardDynamicsLink` alias (line 114) — unchanged
- `#include` block (lines 54-61) — unchanged
- Any other class members

---

### `tests/TestForwardDynamics.cpp` — Update ThreeLinkNumericalValidation (test, request-response)

**Analog:** `tests/TestDynamicsConsistency.cpp` lines 100-177 (`ThreeLinkSerialChain` — round-trip ID→FD check with non-zero COM)

**Imports pattern** (lines 1-7, same as existing):
```cpp
#include "ForwardDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-10;
```

**Existing test pattern to modify** (`ThreeLinkNumericalValidation`, lines 259-303):

The current test (lines 259-303) uses `Vector3d::Zero()` for COM on all links and checks finiteness + `qddot[0] < qddot[2]`. Replace with a round-trip ID→FD check with non-zero COM.

**New test body pattern** (copy from `TestDynamicsConsistency.cpp:100-177` — the `ThreeLinkSerialChain` round-trip):

```cpp
TEST(ForwardDynamicsTest, ThreeLinkNumericalValidation) {
    // Setup 3-link serial chain with non-zero COM
    InverseDynamics id;
    
    InverseDynamicsLink id_l0;
    id_l0.parent = -1;
    id_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));  // non-zero COM
    id_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l0.q = 0.0;
    id_l0.qdot = 0.0;
    id.links.push_back(id_l0);

    InverseDynamicsLink id_l1;
    id_l1.parent = 0;
    id_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));  // non-zero COM
    id_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l1.q = 0.0;
    id_l1.qdot = 0.0;
    id.links.push_back(id_l1);

    InverseDynamicsLink id_l2;
    id_l2.parent = 1;
    id_l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));  // non-zero COM
    id_l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l2.q = 0.0;
    id_l2.qdot = 0.0;
    id.links.push_back(id_l2);

    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.25;

    Eigen::VectorXd tau = id.computeTorques(qddot_input);

    // Forward dynamics
    ForwardDynamics fd;
    
    Link l0;
    l0.parent = -1;
    l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l0.q = 0.0; l0.qdot = 0.0;
    l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l0);

    Link l1;
    l1.parent = 0;
    l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l1.q = 0.0; l1.qdot = 0.0;
    l1.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l1);

    Link l2;
    l2.parent = 1;
    l2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    l2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    l2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    l2.q = 0.0; l2.qdot = 0.0;
    l2.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
    fd.links.push_back(l2);

    fd.computeAccelerations(tau);

    // Verify round-trip: ABA(RNEA(qddot_input)) ≈ qddot_input
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], EPSILON);

    // Structural invariants
    for (int i = 0; i < 3; i++) {
        EXPECT_TRUE(std::isfinite(fd.links[i].qddot));
    }
    EXPECT_LT(fd.links[0].qddot, fd.links[2].qddot);
}
```

**Key differences from existing test:**
- **COM values:** `Vector3d(0.1, 0, 0)` instead of `Vector3d::Zero()` — this is what exercises the bug
- **ID→FD round-trip:** Uses `InverseDynamics::computeTorques()` with non-zero COM → `fd.computeAccelerations()` → verify `qddot_output ≈ qddot_input`
- **EPSILON:** `1e-8` (from `TestDynamicsConsistency.cpp`) — for round-trip checks use `1e-8` not `1e-10`
- **Link naming:** Uses `InverseDynamicsLink` for ID side (from `InverseDynamics.h`), bare `Link` for FD side
- **Test file header comment:** Add `@brief` and `@details` block comment explaining the non-zero COM round-trip

**Important:** The existing `#include` block in `tests/TestForwardDynamics.cpp` is:
```cpp
#include "ForwardDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>
```

Add `#include "InverseDynamics.h"` for the round-trip check (same as `TestDynamicsConsistency.cpp` line 1).

**Test name convention:** Keep `ForwardDynamicsTest.ThreeLinkNumericalValidation` — do not rename.

**Other tests in this file that MUST NOT regress (keep as-is):**
- `SingleLinkPendulum` (lines 19-43)
- `TwoLinkSerialChain` (lines 51-92)
- `BranchingKinematicTree` (lines 100-150)
- `PluckerTransformUsage` (lines 158-197)
- `ZeroTorqueStaticEquilibrium` (lines 202-220)
- `LargeTorqueProportionalAcceleration` (lines 225-251)
- `CondensationReducesInertiaNorm` (lines 310-363)
- `ThreeLinkSingleTorque` (lines 371-413)
- `SingleLinkWithGravity` (lines 423-454)
- `TwoLinkWithGravity` (lines 459-491)
- `GravityEffectScalesWithMass` (lines 514-570)
- `GravityProportionalityInvariant` (lines 590-631)
- `ReleaseModeStability` (lines 642-667)
- `ZeroMassEdgeCase` (lines 676-704)

**main() pattern** (lines 706-709, unchanged):
```cpp
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

---

### `tests/TestDynamicsConsistency.cpp` — Add Non-Zero COM Variants (test, request-response)

**Analog:** `tests/TestDynamicsConsistency.cpp` lines 100-177 (`ThreeLinkSerialChain` — existing round-trip) and lines 184-261 (`BranchingYConfiguration`)

**Imports pattern** (lines 1-6, already present):
```cpp
#include "InverseDynamics.h"
#include "ForwardDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-8;
```

No import changes needed — `InverseDynamics.h` and `ForwardDynamics.h` already included.

**New test: Non-Zero COM ThreeLinkSerialChain** (insert after line 177, before line 179):

This is a copy-paste variant of `ConsistencyTest.ThreeLinkSerialChain` (lines 100-177) with `Vector3d::Zero()` replaced by `Vector3d(0.1, 0, 0)` for COM on all links.

```cpp
/**
 * @brief Three-link serial chain round-trip with non-zero COM
 * @details Same structure as ThreeLinkSerialChain but each link has
 *          COM offset [0.1, 0, 0] along the link X-axis. This exercises
 *          the condensation step in ABA (Phase 14 CR-02 fix).
 *          With zero COM, the skew(com)*mass coupling terms vanish,
 *          masking the missing articulation step. Non-zero COM creates
 *          non-block-diagonal Ia matrices, requiring correct condensation
 *          for multi-link round-trip consistency.
 * 
 *          Independent oracle: ID(FD(qddot)) == qddot by construction
 *          of Featherstone Algorithm 7.3.
 */
TEST(ConsistencyTest, ThreeLinkSerialChainNonZeroCOM) {
    // Inverse dynamics setup
    InverseDynamics id;
    
    InverseDynamicsLink id_link0;
    id_link0.parent = -1;
    id_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link0.q = 0.0;
    id_link0.qdot = 0.0;
    id.links.push_back(id_link0);
    
    InverseDynamicsLink id_link1;
    id_link1.parent = 0;
    id_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link1.q = 0.0;
    id_link1.qdot = 0.0;
    id.links.push_back(id_link1);
    
    InverseDynamicsLink id_link2;
    id_link2.parent = 1;
    id_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link2.q = 0.0;
    id_link2.qdot = 0.0;
    id.links.push_back(id_link2);
    
    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.25;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input);
    
    // Forward dynamics setup (same structure, non-zero COM)
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_link0;
    fd_link0.parent = -1;
    fd_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link0.q = 0.0;
    fd_link0.qdot = 0.0;
    fd.links.push_back(fd_link0);
    
    ForwardDynamicsLink fd_link1;
    fd_link1.parent = 0;
    fd_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link1.q = 0.0;
    fd_link1.qdot = 0.0;
    fd.links.push_back(fd_link1);
    
    ForwardDynamicsLink fd_link2;
    fd_link2.parent = 1;
    fd_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link2.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link2.q = 0.0;
    fd_link2.qdot = 0.0;
    fd.links.push_back(fd_link2);
    
    fd.computeAccelerations(tau);
    
    // Verify accelerations match
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], EPSILON);
}
```

**New test: Non-Zero COM BranchingYConfiguration** (insert after line 261, before line 263):

Copy-paste variant of `ConsistencyTest.BranchingYConfiguration` (lines 184-261) with non-zero COM.

```cpp
/**
 * @brief Branching Y-configuration round-trip with non-zero COM
 * @details Same structure as BranchingYConfiguration but each link has
 *          COM offset [0.1, 0, 0]. With non-zero COM, the articulated
 *          body inertia matrices are non-block-diagonal, requiring
 *          correct condensation for multi-link consistency.
 */
TEST(ConsistencyTest, BranchingYConfigurationNonZeroCOM) {
    // Inverse dynamics setup
    InverseDynamics id;
    
    InverseDynamicsLink id_link0;
    id_link0.parent = -1;
    id_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_link0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link0.q = 0.0;
    id_link0.qdot = 0.0;
    id.links.push_back(id_link0);
    
    InverseDynamicsLink id_link1;
    id_link1.parent = 0;
    id_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_link1.I = RigidBodyInertia(0.5, Vector3d(0.1, 0, 0), lt::Identity(3) * 0.5);
    id_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link1.q = 0.0;
    id_link1.qdot = 0.0;
    id.links.push_back(id_link1);
    
    InverseDynamicsLink id_link2;
    id_link2.parent = 0;
    id_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    id_link2.I = RigidBodyInertia(0.5, Vector3d(0.1, 0, 0), lt::Identity(3) * 0.5);
    id_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_link2.q = 0.0;
    id_link2.qdot = 0.0;
    id.links.push_back(id_link2);
    
    Eigen::VectorXd qddot_input(3);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    qddot_input[2] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input);
    
    // Forward dynamics setup (same structure, non-zero COM)
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_link0;
    fd_link0.parent = -1;
    fd_link0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_link0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_link0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link0.q = 0.0;
    fd_link0.qdot = 0.0;
    fd.links.push_back(fd_link0);
    
    ForwardDynamicsLink fd_link1;
    fd_link1.parent = 0;
    fd_link1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_link1.I = RigidBodyInertia(0.5, Vector3d(0.1, 0, 0), lt::Identity(3) * 0.5);
    fd_link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link1.q = 0.0;
    fd_link1.qdot = 0.0;
    fd.links.push_back(fd_link1);
    
    ForwardDynamicsLink fd_link2;
    fd_link2.parent = 0;
    fd_link2.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(-1, 0, 0));
    fd_link2.I = RigidBodyInertia(0.5, Vector3d(0.1, 0, 0), lt::Identity(3) * 0.5);
    fd_link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_link2.q = 0.0;
    fd_link2.qdot = 0.0;
    fd.links.push_back(fd_link2);
    
    fd.computeAccelerations(tau);
    
    // Verify accelerations match
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
    EXPECT_NEAR(fd.links[2].qddot, qddot_input[2], EPSILON);
    
    // Verify symmetric branches have same acceleration
    EXPECT_NEAR(fd.links[1].qddot, fd.links[2].qddot, EPSILON);
}
```

**New test: Non-Zero COM TwoLinkRoundTripWithGravity** (insert after line 364, before line 366):

Copy-paste variant of `ConsistencyTest.TwoLinkRoundTripWithGravity` (lines 309-364) with non-zero COM.

```cpp
/**
 * @brief Two-link round-trip consistency with gravity and non-zero COM
 * @details Same structure as TwoLinkRoundTripWithGravity but each link
 *          has COM offset [0.1, 0, 0]. Verifies that gravity coupling
 *          through non-zero COM terms doesn't break condensation.
 */
TEST(ConsistencyTest, TwoLinkRoundTripWithGravityNonZeroCOM) {
    Vector3d gravity(0, 0, -9.81);
    
    // RNEA setup
    InverseDynamics id;
    
    InverseDynamicsLink id_l0;
    id_l0.parent = -1;
    id_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    id_l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l0.q = 0.0;
    id_l0.qdot = 0.0;
    id.links.push_back(id_l0);
    
    InverseDynamicsLink id_l1;
    id_l1.parent = 0;
    id_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    id_l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    id_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    id_l1.q = 0.0;
    id_l1.qdot = 0.0;
    id.links.push_back(id_l1);
    
    Eigen::VectorXd qddot_input(2);
    qddot_input[0] = 1.0;
    qddot_input[1] = 0.5;
    
    Eigen::VectorXd tau = id.computeTorques(qddot_input, gravity);
    
    // ABA setup
    ForwardDynamics fd;
    
    ForwardDynamicsLink fd_l0;
    fd_l0.parent = -1;
    fd_l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
    fd_l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_l0.q = 0.0;
    fd_l0.qdot = 0.0;
    fd.links.push_back(fd_l0);
    
    ForwardDynamicsLink fd_l1;
    fd_l1.parent = 0;
    fd_l1.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d(1, 0, 0));
    fd_l1.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
    fd_l1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    fd_l1.q = 0.0;
    fd_l1.qdot = 0.0;
    fd.links.push_back(fd_l1);
    
    fd.computeAccelerations(tau, gravity);
    
    EXPECT_NEAR(fd.links[0].qddot, qddot_input[0], EPSILON);
    EXPECT_NEAR(fd.links[1].qddot, qddot_input[1], EPSILON);
}
```

**Test file structure — what to preserve:**
- Existing tests `RoundTripABARNEA` through `RoundTripABA_RNEA_DirectComparison` (lines 16-423) — NO changes
- New tests inserted after each corresponding zero-COM test
- `main()` function (lines 425-428) unchanged

---

## Shared Patterns

### Authentication / Guard
**Not applicable** — this is a pure computation C++ library. No auth, no middleware.

### Error Handling
**Source:** `src/ForwardDynamics.cpp` lines 105-124 and `include/ForwardDynamics.h` lines 148-155

**Apply to:** `inwardPass()` (all Joint inertia checks)

```cpp
constexpr double EPSILON = 1e-10;

if (std::abs(D) < EPSILON)
{
    throw std::runtime_error(
        "ForwardDynamics::inwardPass: Near-zero inertia at joint "
        + std::to_string(i) + " (D=" + std::to_string(D) + ")"
    );
}
```

And in `computeAccelerations()` (lines 182-202):
```cpp
if (tau.size() != static_cast<int>(links.size()))
{
    throw std::invalid_argument(
        "ForwardDynamics::computeAccelerations: tau size ("
        + std::to_string(tau.size()) + ") does not match link count ("
        + std::to_string(links.size()) + ")"
    );
}

for (int i = 0; i < tau.size(); i++)
{
    if (std::isnan(tau[i]) || std::isinf(tau[i]))
    {
        throw std::invalid_argument(
            "ForwardDynamics::computeAccelerations: Invalid torque at joint "
            + std::to_string(i)
        );
    }
}
```

### Testing
**Source:** `tests/TestDynamicsConsistency.cpp` lines 1-428 and `tests/TestForwardDynamics.cpp` lines 1-709

**Apply to:** All test files

**GTest framework pattern:**
```cpp
#include "ForwardDynamics.h"      // or "InverseDynamics.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;

constexpr double EPSILON = 1e-8;  // For round-trip checks
// OR
constexpr double EPSILON = 1e-10; // For exact-value checks

TEST(SuiteName, TestName) {
    // Setup
    // Execute
    // Verify with EXPECT_NEAR / EXPECT_DOUBLE_EQ / EXPECT_GT / EXPECT_TRUE
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

**Assertion patterns:**
- `EXPECT_NEAR(actual, expected, EPSILON)` — for computed value checks
- `EXPECT_DOUBLE_EQ(actual, expected)` — for exact value checks
- `EXPECT_GT(val, 0.0)` — positivity checks
- `EXPECT_TRUE(std::isfinite(val))` — finiteness checks
- `EXPECT_LT(val_a, val_b)` — ordering checks

**Round-trip test pattern (ID → FD):**
```cpp
// 1. Define qddot_input
Eigen::VectorXd qddot_input(n);
qddot_input[0] = 1.0; qddot_input[1] = 0.5; // etc.

// 2. RNEA: qddot → tau
InverseDynamics id;
// ... setup id.links with same structure ...
Eigen::VectorXd tau = id.computeTorques(qddot_input);

// 3. ABA: tau → qddot_output
ForwardDynamics fd;
// ... setup fd.links with same structure ...
fd.computeAccelerations(tau);

// 4. Verify
EXPECT_NEAR(fd.links[i].qddot, qddot_input[i], EPSILON);
```

**Link setup pattern (existing test style):**
```cpp
Link l0;
l0.parent = -1;
l0.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
l0.I = RigidBodyInertia(1.0, Vector3d(0.1, 0, 0), lt::Identity(3));
l0.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
l0.q = 0.0; l0.qdot = 0.0;
l0.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
fd.links.push_back(l0);
```

**Test naming convention:** PascalCase test names, `Doxygen` block comments with `@brief` and `@details` above each test.

### Doxygen Documentation
**Source:** `include/ForwardDynamics.h` lines 4-52 and `src/ForwardDynamics.cpp` lines 1-11

**Apply to:** All header modifications

```cpp
/**
 * @file FileName.h
 * @brief One-line summary
 * @details Extended description of algorithm, mathematical foundation,
 *          and usage patterns.
 * 
 *          Algorithm Overview:
 *          1. Step one
 *          2. Step two
 * 
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7
 */
```

### Code Style
**Source:** Conventions from AGENTS.md and existing code

| Convention | Pattern |
|-----------|---------|
| Indentation | 4 spaces |
| Braces (classes) | Allman style (newline before `{`) |
| Braces (functions) | K&R style (same line) |
| Namespace | Everything in `namespace SpatialAlgebra` |
| Type aliases | `using mv = MotionVector`, `using fv = ForceVector`, `using lt = LowerTriangular` |
| Include guards | `#ifndef FILE_NAME_H` / `#define FILE_NAME_H` |
| Eigen types | `Vector3d`, `Matrix3d`, `VectorXd` |
| Variable naming | camelCase: `invD`, `IaS`, `inertiaCorr`, `HCorr`, `massCorr` |
| Comments | Doxygen block comments on every function/class, inline comments for formulas |

### Build / CMake
**Source:** `CMakeLists.txt` lines 122-130 (existing `TestForwardDynamics`) and 152-160 (existing `TestDynamicsConsistency`)

**No changes needed** — test executables already registered. No new test files.

```cmake
add_executable(TestForwardDynamics tests/TestForwardDynamics.cpp)
target_link_libraries(TestForwardDynamics
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)
add_test(NAME TestForwardDynamics COMMAND TestForwardDynamics)

add_executable(TestDynamicsConsistency tests/TestDynamicsConsistency.cpp)
target_link_libraries(TestDynamicsConsistency
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)
add_test(NAME TestDynamicsConsistency COMMAND TestDynamicsConsistency)
```

---

## No Analog Found

All files have strong existing analogs. This phase modifies existing files following their own patterns; no new file types or patterns are introduced.

---

## Metadata

**Analog search scope:** `include/`, `src/`, `tests/`, `CMakeLists.txt`
**Files scanned:** 10 (4 target + 6 analog/support)
**Pattern extraction date:** 2026-06-17
**Phase reference:** Featherstone Algorithm 7.3 (condensation step in ABA inward pass)
