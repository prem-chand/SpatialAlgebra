# SpatialAlgebra × Pinocchio — Cross-Validation Report

**Generated**: 2026-06-18T12:38:20.756996Z
**Solver A**: SpatialAlgebra (C++)
**Solver B**: Pinocchio (C++)
**Benchmark**: `/Users/premchand/Documents/GitHub/SpatialAlgebra/build/tests/test-models/bench_pinocchio_comparison`
**Tolerance**: kinematics ≤ 1e-8, dynamics ≤ 1e-6

**Summary**: 4/8 models pass, 4 fail

---
## BranchingYConfiguration

**DOF**: 3  |  **Layout**: 3-DOF tree: base → left (Z-rev) + right (Z-rev), each offset  |  **Status**: ✅ PASS

*Branching topology test*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| ABA (qddot = FD(τ, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Mass matrix (CRBA / RNEA-column) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## SingleLinkChain

**DOF**: 1  |  **Layout**: 1-DOF Z-revolute at origin  |  **Status**: ✅ PASS

*Featherstone (2008) §5.1 — simplest non-trivial RNEA*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| ABA (qddot = FD(τ, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Mass matrix (CRBA / RNEA-column) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## SingleLinkXAxisCOMOffset

**DOF**: 1  |  **Layout**: 1-DOF X-revolute with COM offset at (0, 0.5, 0)  |  **Status**: ❌ FAIL

*Parallel axis: I_joint(XX) = 1 + 1*0.25 = 1.25*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 3.74e-01 | 3.74e-01 | 1e-06 | ✗ |
| ABA (qddot = FD(τ, g)) | 1.52e-01 | 1.52e-01 | 1e-06 | ✗ |
| Mass matrix (CRBA / RNEA-column) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Gravity compensation torques | 1.29e+00 | 1.29e+00 | 1e-06 | ✗ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## ThreeLinkSerialChain

**DOF**: 3  |  **Layout**: 3-DOF serial, Z-revolute, 1m X spacing  |  **Status**: ❌ FAIL

*3-link; M[0,0]=8 at q=0 (parallel axis accumulation)*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 1.24e+00 | 5.56e-01 | 1e-06 | ✗ |
| ABA (qddot = FD(τ, g)) | 9.53e-01 | 6.30e-01 | 1e-06 | ✗ |
| Mass matrix (CRBA / RNEA-column) | 3.52e-01 | 1.17e-01 | 1e-06 | ✗ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## ThreeLinkSerialNonZeroCOM

**DOF**: 3  |  **Layout**: 3-DOF Z-revolute, each link COM=(0.1, 0, 0), 1m spacing  |  **Status**: ❌ FAIL

*3-link with parallel axis corrections*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 1.26e+00 | 6.59e-01 | 1e-06 | ✗ |
| ABA (qddot = FD(τ, g)) | 9.62e-01 | 6.86e-01 | 1e-06 | ✗ |
| Mass matrix (CRBA / RNEA-column) | 4.04e-01 | 1.87e-01 | 1e-06 | ✗ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## TwoLinkSerialChain

**DOF**: 2  |  **Layout**: 2-DOF serial, Z-revolute, 1m X offset  |  **Status**: ✅ PASS

*Two-link serial; M[0,0]=3 at q=0, M[1,1]=1*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| ABA (qddot = FD(τ, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Mass matrix (CRBA / RNEA-column) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## TwoLinkWithCOMOffset

**DOF**: 2  |  **Layout**: 2-DOF Z-revolute, each link COM=(0.1, 0, 0)  |  **Status**: ❌ FAIL

*2-link with parallel axis from COM drift*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 7.71e-02 | 5.57e-02 | 1e-06 | ✗ |
| ABA (qddot = FD(τ, g)) | 3.61e-01 | 2.00e-01 | 1e-06 | ✗ |
| Mass matrix (CRBA / RNEA-column) | 9.05e-02 | 6.04e-02 | 1e-06 | ✗ |
| Gravity compensation torques | 4.53e-01 | 3.53e-01 | 1e-06 | ✗ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## TwoLinkWithGravity

**DOF**: 2  |  **Layout**: 2-DOF Z-revolute, gravity compensation test  |  **Status**: ✅ PASS

*Gravity compensation: g=(0,0,-9.81), tau_g = J^T * F_g*

| Check | Max Rel Err | Mean Rel Err | Threshold | Status |
|-------|-------------|-------------|-----------|--------|
| RNEA (τ = ID(qddot, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| ABA (qddot = FD(τ, g)) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Mass matrix (CRBA / RNEA-column) | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Gravity compensation torques | 0.00e+00 | 0.00e+00 | 1e-06 | ✓ |
| Forward kinematics (joint transforms) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |
| Center of mass (world-frame COM) | 0.00e+00 | 0.00e+00 | 1e-08 | ✓ |

---
## Cross-Model Summary

| Model | Overall | RNEA | ABA | Mass Matrix | Gravity | FK | COM |
|-------|---------|------|-----|-------------|---------|-----|-------|
| BranchingYConfiguration | ✅ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| SingleLinkChain | ✅ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| SingleLinkXAxisCOMOffset | ❌ | ✗ | ✗ | ✓ | ✗ | ✓ | ✓ |
| ThreeLinkSerialChain | ❌ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ |
| ThreeLinkSerialNonZeroCOM | ❌ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ |
| TwoLinkSerialChain | ✅ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| TwoLinkWithCOMOffset | ❌ | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ |
| TwoLinkWithGravity | ✅ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

---
*Report generated by Phase 24 — Result Reporting*
