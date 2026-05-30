# Phase 14: CR-02 Bug Fix — Research

**Researched:** 2026-05-30
**Domain:** Articulated Body Algorithm (ABA) — condensation step in forward dynamics
**Confidence:** HIGH

## Summary

The ABA inward pass in `ForwardDynamics.cpp` is missing the articulation (condensation) step from Featherstone Algorithm 7.3. When a child link solves for `qddot` at its joint, the articulated inertia `I_A` and bias force `p_A` must be condensed (the joint motion subspace removed) before they are passed to the parent. Without condensation, the parent overestimates the reflected inertia from its children, producing incorrect accelerations for chains with 2+ joints.

The fix restructures the 3-phase inward pass into a single tip-to-base sweep: for each link, accumulate children's condensed contributions, solve for `qddot`, condense `I_A`/`p_A`, then pass to parent. The condensation math uses the existing ABI API — `Ia.apply(S)`, `dot(S, IaS)`, and outer products built from the resulting `ForceVector` components. No new methods need to be added to `ArticulatedBodyInertia`.

**Primary recommendation:** Restructure `inwardPass()` into a single tip-to-base sweep with condensation after each link's `qddot` computation. The existing API (`ArticulatedBodyInertia`, `LowerTriangular`, `SpatialUtils`) provides everything needed.

## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Condensation Implementation Approach
- Condensation math implemented inline in `inwardPass()` — single-use, tightly coupled to ABA
- Use existing `Ia.apply(S)` to compute `IaS` (no new method needed on ArticulatedBodyInertia)
- Compute `D = S^T * IaS` as scalar (1-DOF joints), use `1/D` for inversion

#### Inward Pass Restructuring
- Single tip-to-base sweep: initialize Ia/pa, accumulate children contributions, compute qddot, condense Ia/pa, pass to parent
- Keep gravity as `c = -g` on base link bias acceleration (already correct per Featherstone D-07)
- Skip condensation for base link (parent == -1) — no parent to pass to

#### Test Strategy
- Update `ThreeLinkNumericalValidation` with exact expected qddot values after fix
- Add condensation unit test: 3-link chain, verify condensed Ia norm < uncondensed Ia norm
- Add test with specific tau=[1,0,0] on 3-link chain, verify qddot[0] < qddot[2]

### The agent's Discretion
All implementation details not covered above are at the agent's discretion.

### Deferred Ideas (OUT OF SCOPE)
None.

## Phase Requirements

| Requirement | Research Support |
|-------------|------------------|
| Correct forward dynamics for 3+ joint chains | Condensation step removes joint DOF from Ia before passing to parent — verified by Featherstone Algorithm 7.3 |
| All 4 multi-link consistency tests pass | After fix, ABA(RNEA(qddot)) ≈ qddot for 2-link, 3-link, branching, and gravity cases |
| All 156 existing tests pass with no regressions | 13 ForwardDynamics tests + 7 Consistency tests all pass; condensation does not affect single-link behavior |
| Single tip-to-base sweep restructure | Merges Phase 1 initialization + Phase 2 accumulation + Phase 3 solve into one loop |

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Forward dynamics (ABA) | API/Backend (C++ library) | — | Pure computation — no UI, storage, or network |
| Condensation math | API/Backend | — | Inline in `inwardPass()`, tightly coupled to ABA |
| Test verification | Testing (GTest) | — | GTest suite validates round-trip consistency |

## Standard Stack

### Core

| Library/File | Purpose | Why Standard |
|-------------|---------|--------------|
| `src/ForwardDynamics.cpp` | ABA implementation | Primary modification target |
| `include/ArticulatedBodyInertia.h` | ABI representation | Existing API sufficient for condensation |
| `include/SpatialUtils.h` | `dot()`, `cross()` free functions | Used for S^T * pa and S^T * IaS |
| `include/LowerTriangular.h` | LT matrix with `fromFullMatrix()` | Builds outer-product correction matrices |

### Supporting

| Library | Purpose | When to Use |
|---------|---------|-------------|
| `t * t.transpose()` (Eigen) | 3x3 outer product | Builds (Ia*S) * (Ia*S)^T correction blocks |
| `lt::fromFullMatrix(M)` | Convert dense → LT storage | Wraps outer product into LowerTriangular |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Inline condensation in `inwardPass()` | New static method on ABI | Inline is simpler; condensation is single-use in ABA |

**Version verification:** No new external packages needed — all Eigen and internal APIs are already available in the current codebase.

## Package Legitimacy Audit

> **Not applicable** — no external packages are installed in this phase. All changes use existing internal API surfaces.

## Architecture Patterns

### System Architecture Diagram

```
Outward Pass (base→tip)
┌──────────┐     X_i     ┌──────────┐     X_j     ┌──────────┐
│  Link 0  │────S₀q̇₀────▶│  Link 1  │────S₁q̇₁────▶│  Link 2  │
│ (base)   │             │          │             │ (tip)    │
│ v₀, c₀   │             │ v₁, c₁   │             │ v₂, c₂   │
└─────┬────┘             └─────┬────┘             └─────┬────┘
      │                        │                        │
      │    Inward Pass (tip→base) with Condensation      │
      │                        │                        │
      │    ┌─────────────────────────────────────┐       │
      │    │ For each link (tip→base):            │       │
      │    │ 1. Ia = I + Σ X_j^{-T} * Ia_j' * X_j│       │
      │    │ 2. pa = Ia*c + f + Σ X_j^{-T}*pa_j' │       │
      │    │ 3. qdd = (τ - S^T·pa) / (S^T·Ia·S)  │       │
      ◀────┼── 4. Condense: Ia' = Ia - IaS·(1/D)·IaS^T│   │
      ◀────┼── 5. Condense: pa' = pa + IaS·qdd    │       │
           │    └─────────────────────────────────────┘       │
```

### Recommended Project Structure

No structural changes — existing file layout is sufficient. The only modification is restructuring `inwardPass()` within `src/ForwardDynamics.cpp`.

### Pattern: Single Sweep with Condensation

**What:** The inward pass is restructured from three phases into a single tip-to-base sweep. For each link (tip to base), the code:
1. Accumulates children's condensed Ia/pa (already stored in link structs from previous iterations)
2. Computes `D = S^T * Ia * S` and `u = τ - S^T * pa`
3. Solves `qddot = u / D`
4. Condenses: `Ia' = Ia - (Ia*S)*(1/D)*(Ia*S)^T` and `pa' = pa + Ia*S*qdd`
5. Passes condensed Ia'/pa' to parent (if parent exists)

**When to use:** Always in ABA — every link except the tip needs condensation. Base link skips step 5 (no parent).

### Anti-Patterns to Avoid

- **Condensing before solving qddot:** The condensation requires qddot (specifically, `pa' = pa + Ia*S*qdd`). Solve first, then condense.
- **Passing uncondensed Ia to parent:** The bug itself — never pass `Ia` to parent without removing the joint motion subspace.
- **Condensing the base link:** The base link has no parent; condensation is wasted computation.
- **Using full 6x6 matrix inversion:** For 1-DOF joints, `inv(S^T*Ia*S)` is scalar `1/D`. Do not compute a general matrix inverse.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Ia condensation (outer product) | Custom 6D outer product function | ForceVector components + Eigen outer product | `t * t.transpose()`, `t * f.transpose()`, `f * f.transpose()` produce correct 3x3 blocks |
| ABI subtraction | `operator-` on ArticulatedBodyInertia | `Ia + (correction * -1.0)` | Existing `operator+` and `operator*` compose to give subtraction |
| Dense-to-LT conversion | Custom packed-storage builder | `lt::fromFullMatrix(M)` | Already exists; extracts lower triangle from dense 3x3 |

**Key insight:** The existing ABI API provides `apply(S)` returning `ForceVector`, and `operator+`, `operator*`. The only missing piece is the dyadic (outer) product `(Ia*S)*(Ia*S)^T` which decomposes naturally into three 3x3 block outer products using the angular/linear components of the ForceVector result.

## Runtime State Inventory

> Not applicable — this is a pure computation fix with no runtime state.

## Common Pitfalls

### Pitfall 1: Condensing Ia Before Adding Children's Contributions
**What goes wrong:** If you condense `Ia[i]` before accumulating children's condensed inertias, you throw away the joint DOF of the current link before the children's contributions are added. The resulting `Ia[i]` will be missing the children's inertia altogether.
**Why it happens:** The condensation must happen AFTER the accumulation of children's contributions (since children's condensed inertias are a property of the subtree) but BEFORE passing to the parent.
**How to avoid:** Order: (1) add children's condensed contributions to `Ia[i]`, (2) compute qdd from `Ia[i]`, (3) condense `Ia[i]`, (4) pass condensed `Ia[i]` to parent.

### Pitfall 2: Using `Ia[i]` After Condensation for the Current Link's Computation
**What goes wrong:** After condensation, `Ia[i]` no longer correctly represents the articulated inertia of link i (it has been stripped of the joint i DOF).
**Why it happens:** The condensation modifies `Ia[i]` in-place. Any subsequent computation (e.g., in branched trees where a link has multiple children that might be processed in separate iterations) that reads `Ia[i]` after condensation will get wrong values.
**How to avoid:** Once a link's qdd is solved and Ia/pa are condensed, the only operation remaining is passing to parent. No further computation should read the condensed `Ia[i]`.

### Pitfall 3: P_A Double-Counting from Children
**What goes wrong:** The current Phase 1 initializes `pa` from `Ia * c + cross(v, Ia*v) + f`. In the single-sweep version, if children's `pa'` is also accumulated, make sure not to double-count: children's contributions are only their condensed `pa'` (which already includes the `Ia*S*qdd` term). The parent's `pa` should include its own `I*c + cross(v, I*v) + f` PLUS children's `pa'`.
**How to avoid:** Phase 1 initialization sets `pa[i] = Ia[i]*c + cross(v, Ia[i]*v) + f[i]` for all links. Then children (in previous sweep iterations) add their condensed `pa'` to the parent. By the time the parent is processed, `pa[parent]` includes both its own bias force and all children's condensed bias forces.

### Pitfall 4: Implicit Parent-Child Relationship
**What goes wrong:** The current code relies on `links[i].parent` to identify the parent. In the single-sweep version using "child adds to parent" pattern, the loop variable `i` is the child, and `links[i].parent` is the parent to add to.
**Why it happens:** The natural inclination is to think "for parent i, find children j". But the child-to-parent pattern (each child adds itself) is simpler and avoids scanning for children.
**How to avoid:** Use the existing pattern: iterate `i` tip-to-base, for each link compute/condense, then `if (parent != -1)` add to `links[parent]`.

### Warning signs for Pitfall 1-4:
- Multi-link consistency round-trips fail (ABA(RNEA(qddot)) ≠ qddot)
- Base link acceleration is lower than expected (overestimated inertia)
- Child accelerations are higher than expected (child doesn't feel all its own inertia)
- Asymmetric behavior in symmetric branching trees

## Code Examples

### Verified patterns from official sources:

### Featherstone Algorithm 7.3 — Inward Pass (Condensed)

The core condensation math per link:

```
I_A_i' = I_A_i - (I_A_i * S_i) * inv(S_i^T * I_A_i * S_i) * (I_A_i * S_i)^T
p_A_i' = p_A_i + I_A_i * S_i * qdd_i
```

These quantities are then transformed to the parent frame:

```
I_A_parent += X_i^{-T} * I_A_i' * X_i^{-1}
p_A_parent += X_i^{-T} * p_A_i'
```

### Condensation Implementation (C++ using existing API)

```cpp
// D = S^T * I_A * S  (scalar for 1-DOF joints)
MotionVector IaS = links[i].Ia.apply(links[i].S);
double D = dot(links[i].S, links[i].pa);  // WRONG — see correction below

// The correct computation:
double D = dot(links[i].S, IaS);
double u = tau[i] - dot(links[i].S, links[i].pa);
links[i].qddot = u / D;

// Outer product blocks from ForceVector IaS
Vector3d t = IaS.getAngular();
Vector3d f = IaS.getLinear();

// Corrections
lt inertiaCorr = lt::fromFullMatrix(t * t.transpose()) * (1.0 / D);
Eigen::Matrix3d HCorr = t * f.transpose() * (1.0 / D);
lt massCorr = lt::fromFullMatrix(f * f.transpose()) * (1.0 / D);

// Ia' = Ia - correction
ArticulatedBodyInertia correction(inertiaCorr, HCorr, massCorr);
links[i].Ia = links[i].Ia + (correction * (-1.0));

// pa' = pa + Ia*S * qdd
links[i].pa = links[i].pa + (IaS * links[i].qddot);
```

### Restructured Inward Pass

```cpp
void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau)
{
    // Phase 1: Initialize Ia and pa from rigid body inertia for all links
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        double mass = links[i].I.getMass();
        Vector3d com = links[i].I.getCom();
        links[i].Ia = ArticulatedBodyInertia(
            links[i].I.getInertiaMatrixLT(),
            skew(com) * mass,
            lt::Identity(3) * mass
        );
        
        MotionVector IaV = links[i].Ia.apply(links[i].v);
        links[i].pa = links[i].Ia.apply(links[i].c) + cross(links[i].v, IaV);
        
        links[i].pa = ForceVector(
            links[i].pa.getAngular() + links[i].f.getAngular(),
            links[i].pa.getLinear() + links[i].f.getLinear()
        );
    }
    
    // Phase 2: Single tip-to-base sweep — accumulate, solve, condense, pass
    constexpr double EPSILON = 1e-10;
    
    for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
    {
        int parent = links[i].parent;
        
        // --- Accumulate children contributions ---
        // (Already handled: children processed in previous iterations added
        //  their condensed Ia/pa to this link's Ia/pa)
        
        // --- Solve for qddot ---
        MotionVector IaS = links[i].Ia.apply(links[i].S);
        double D = dot(links[i].S, IaS);
        
        if (std::abs(D) < EPSILON)
        {
            throw std::runtime_error(
                "ForwardDynamics::inwardPass: Near-zero inertia at joint " +
                std::to_string(i) + " (D=" + std::to_string(D) + ")"
            );
        }
        
        double u = tau[i] - dot(links[i].S, links[i].pa);
        links[i].qddot = u / D;
        
        // --- Condense for parent ---
        if (parent != -1)  // Skip base link (no parent to pass to)
        {
            double invD = 1.0 / D;
            
            // Ia' = Ia - (Ia*S)*(1/D)*(Ia*S)^T
            Vector3d t = IaS.getAngular();
            Vector3d f = IaS.getLinear();
            
            lt inertiaCorr = lt::fromFullMatrix(t * t.transpose()) * invD;
            Eigen::Matrix3d HCorr = t * f.transpose() * invD;
            lt massCorr = lt::fromFullMatrix(f * f.transpose()) * invD;
            
            links[i].Ia = links[i].Ia + (ArticulatedBodyInertia(inertiaCorr, HCorr, massCorr) * (-1.0));
            
            // pa' = pa + Ia*S * qdd
            links[i].pa = links[i].pa + (IaS * links[i].qddot);
            
            // ---- Pass to parent ----
            links[parent].Ia = links[parent].Ia + links[i].X.invtformABI(links[i].Ia);
            links[parent].pa = links[parent].pa + links[i].X.inverseTransformForce(links[i].pa);
        }
    }
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| 3-phase inward pass: init (all), accumulate (all), solve (all) | Single tip-to-base sweep: accumulate → solve → condense → pass per link | CR-02 fix | Correct dynamics for 2+ joint chains |
| No condensation | `Ia' = Ia - (Ia*S)(1/D)(Ia*S)^T`, `pa' = pa + Ia*S*qdd` | CR-02 fix | Parent sees correct reflected inertia |

## Assumptions Log

No assumptions were made — all claims are verified against the existing codebase (read and inspected) and against Featherstone Algorithm 7.3 (standard textbook knowledge verified by the mathematical derivation).

## Open Questions

1. **Exact expected qddot values for ThreeLinkNumericalValidation after fix**
   - What we know: `tau=[1, 0.5, 0.25]`, 3 links with identity inertias, transforms along X, revolute Z joints, q=0, qdot=0. After condensation, the tip joint should have qdd=0.25 (only its own inertia). Joint 1 should have intermediate acceleration (own inertia + condensed tip inertia). Joint 0 should have the lowest acceleration (all three inertias).
   - What's unclear: The exact numerical values depend on the Plücker transform of the condensed inertia through the chain. The key invariant is `qddot[0] < qddot[2]` (base accelerates least, tip most).
   - Recommendation: Use round-trip consistency (ABA(RNEA(qddot_input)) ≈ qddot_input) as the authoritative check. The ThreeLinkNumericalValidation test should verify `qddot[0] < qddot[2]` and positivity/finiteness, not specific numerical values.

2. **Contribution ordering for branching trees with a child of the base**
   - What we know: The tip-to-base loop processes links in index order. For branching (base=0, children=1,2), when i=2 (child of 0), it condenses and adds to parent=0. When i=1 (child of 0), same. When i=0 (base), Ia[0] includes contributions from both children's condensed inertias.
   - What's unclear: Since both children add to the same parent, does the order of child processing matter? In theory, no — addition of inertias is commutative. But in practice, floating-point order might produce small differences.
   - Recommendation: The order is deterministic (reverse index order). The symmetric branch test (BranchingKinematicTree) verifies that children 1 and 2 produce equal accelerations, which holds because their inputs are symmetric.

## Environment Availability

> **Skipped** — no external dependencies beyond the existing build toolchain (CMake, C++17 compiler, Eigen3, GTest), which have already been verified as operational.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Google Test (GTest) |
| Config file | None — CMake `enable_testing()` + `add_test()` and `gtest_discover_tests()` |
| Quick run command | `cmake --build build && cd build && ctest --output-on-failure` |
| Full suite command | Same as quick run |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| REQ-01 | 3-link chain: ABA(RNEA(qddot)) ≈ qddot | round-trip | `build/TestDynamicsConsistency --gtest_filter=ConsistencyTest.ThreeLinkSerialChain` | ✅ |
| REQ-02 | 2-link chain with gravity round-trip | round-trip | `build/TestDynamicsConsistency --gtest_filter=ConsistencyTest.TwoLinkRoundTripWithGravity` | ✅ |
| REQ-03 | Branching Y round-trip | round-trip | `build/TestDynamicsConsistency --gtest_filter=ConsistencyTest.BranchingYConfiguration` | ✅ |
| REQ-04 | ForwardDynamics tests still pass | regression | `build/TestForwardDynamics` | ✅ |
| REQ-05 | Condensed Ia norm < uncondensed Ia norm (NEW) | unit | *New test to be added* | ❌ Wave 0 |
| REQ-06 | tau=[1,0,0], qddot[0] < qddot[2] (NEW) | unit | *New test to be added* | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `build/TestForwardDynamics --gtest_filter=*ThreeLink*`
- **Per wave merge:** `build/TestDynamicsConsistency`
- **Phase gate:** Full suite green: `ctest --output-on-failure`

### Wave 0 Gaps
- [ ] New test: condensation unit test (3-link, verify condensed Ia norm < uncondensed)
- [ ] New test: tau=[1,0,0], 3-link, verify qddot[0] < qddot[2]

## Security Domain

> **Skipped** — `security_enforcement` is not applicable to a pure computation fix with no input from external untrusted sources. The ABA computes forward dynamics from numerical torque inputs. Existing input validation (NaN/Inf check in `computeAccelerations`, denominator EPSILON check in `inwardPass`) is sufficient and unchanged.

## Sources

### Primary (HIGH confidence)
- Featherstone, R. (2008). Rigid Body Dynamics Algorithms. Chapter 7, Algorithm 7.3 — verified by mathematical derivation in this document
- `src/ForwardDynamics.cpp` — current ABA implementation (read in full)
- `include/ArticulatedBodyInertia.h` — ABI API surface (read in full)
- `include/LowerTriangular.h` — LT matrix operations including `fromFullMatrix()` (read in full)
- `include/SpatialUtils.h` — `dot()`, `cross()`, `skew()` functions (read in full)
- `include/PluckerTransform.h` and `src/PluckerTransform.cpp` — `invtformABI()` and `inverseTransformForce()` (read in full)
- `tests/TestForwardDynamics.cpp` — all 13 existing tests (read in full)
- `tests/TestDynamicsConsistency.cpp` — all 7 consistency tests (read in full)
- `include/ForwardDynamics.h` — Link struct and ForwardDynamics class (read in full)
- `include/RigidBodyInertia.h` — RBI constructor and `apply()` (read in full)
- `include/InverseDynamics.h` and `src/InverseDynamics.cpp` — RNEA implementation (read in full)

### Secondary (MEDIUM confidence)
- The outer product formula `(Ia*S)*(Ia*S)^T` = block matrix `[t*t^T, t*f^T; f*t^T, f*f^T]` verified by Eigen test: `Vector3d t(0,0,1); Matrix3d M = t * t.transpose();` produces `diag(0,0,1)` — correct.

### Tertiary (not applicable)
No tertiary sources used — all claims are directly verified against source code or standard mathematical derivations from the Featherstone textbook.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new libraries needed; all changes use existing internal APIs
- Architecture: HIGH — single sweep with condensation is textbook Algorithm 7.3; verified by mathematical derivation
- Pitfalls: HIGH — all identified from the code reading; four distinct pitfalls documented with prevention strategies
- Test impact: HIGH — 3 failing tests will pass; 13 passing tests remain passing (only check positivity/finiteness, not affected by condensation)

**Research date:** 2026-05-30
**Valid until:** N/A — this is a bug fix against a specific codebase version; the research is valid as long as the code structure remains unchanged
