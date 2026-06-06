# Phase 18 Discussion Log

**Date:** 2026-06-06
**Phase:** 18 — Robot Examples

## Areas Discussed

### 1. Link Parameter Realism
- **Options presented:** (a) Thin-rod model (mass = 1 kg/m, COM at midpoint), (b) Specific robot params
- **User selection:** Specific robot params
- **Follow-up on which robot:** UR5-like parameters
- **Decision:** D-01 (UR5-like mass/inertia/length, Z-down gravity), D-02 (gravity = Vector3d(0,0,-9.81))

### 2. Gravity Verification Strategy
- **Options presented:** (a) FD+ID cross-validation, (b) Gravity torques only
- **User selection:** FD+ID cross-validation
- **Decision:** D-03, D-04, D-12 (ID computes gravity torques → FD verifies zero acceleration)

### 3. 2-link vs 3-link Code Reuse
- **Options presented:** (a) Fully self-contained, (b) Shared robot_utils.h
- **User selection:** Fully self-contained
- **Decision:** D-05, D-06, D-07 (standalone files, explicit build targets)

### 4. Output Formatting
- **Options presented:** (a) Verbose, (b) Concise/table-oriented
- **User selection:** Verbose (matching existing examples)
- **Decision:** D-08 (section headers, explanatory text, physical interpretation)

## Deferred Ideas

None.
