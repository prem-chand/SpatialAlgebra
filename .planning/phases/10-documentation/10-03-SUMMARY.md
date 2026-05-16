---
phase: 10-documentation
plan: 03
subsystem: documentation
tags: [examples, tutorials, cmake]
dependency_graph:
  requires: [10-01]
  provides: [DOC-03]
  affects: [README.md, CMakeLists.txt]
tech_stack:
  added: []
  patterns: []
key_files:
  created: [examples/CMakeLists.txt, examples/basic_vectors.cpp, examples/transforms.cpp, examples/inertia.cpp, examples/dynamics.cpp]
  modified: [CMakeLists.txt, README.md]
decisions:
  - "Used direct #include without SpatialAlgebra/ prefix (headers in include/ directory)"
  - "Simplified dynamics example to match actual ForwardDynamics API (uses links vector)"
  - "Used LowerTriangular packed format for inertia tensors"
metrics:
  duration_minutes: 25
  completed: "2026-05-16"
---

# Phase 10 Plan 03: Compilable Examples Summary

**One-liner:** Created 4 working example programs demonstrating vectors, transforms, inertia, and forward dynamics with CMake integration.

## Summary

Created comprehensive usage examples in `examples/` directory:

1. **basic_vectors.cpp** (108 lines) - Demonstrates:
   - MotionVector and ForceVector creation
   - Vector operations (add, subtract, scale)
   - Cross products (crossMotion, crossForce)
   - Dot product for power computation
   - Type aliases (mv, fv)

2. **transforms.cpp** (120 lines) - Demonstrates:
   - Rotation from angle-axis (Eigen::AngleAxisd)
   - PluckerTransform creation and usage
   - transformMotion() and transformForce()
   - Inverse transform computation
   - Physical interpretation of coordinate changes

3. **inertia.cpp** (130 lines) - Demonstrates:
   - RigidBodyInertia with LowerTriangular packed format
   - apply() to compute force from motion
   - Offset center of mass effects
   - tformRBI() to transform inertia
   - Property getters (getMass, getCom, getInertiaMatrixLT)

4. **dynamics.cpp** (156 lines) - Demonstrates:
   - ForwardDynamics with Link structure
   - Kinematic chain setup (parent indices, transforms, joint axes)
   - Articulated Body Algorithm (ABA) for forward dynamics
   - Multiple test cases with different torque configurations

## Verification

All examples compile and run successfully:

```
✓ example_vectors - Shows vector operations with correct output
✓ example_transforms - Demonstrates 90° Z rotation correctly
✓ example_inertia - Shows F=ma and torque from offset COM
✓ example_dynamics - ABA produces physically sensible accelerations
```

## Deviations from Plan

### API Adjustments (Rule 3 - Auto-fix blocking issues)

1. **Rotation API**: `setFromAngleAxis()` takes `Eigen::AngleAxisd`, not separate angle/axis
2. **LowerTriangular API**: Uses `operator()` for element access, not `set()` method
3. **RigidBodyInertia**: Requires `LowerTriangular` for inertia tensor, not `Matrix3d`
4. **ForwardDynamics**: Uses `links` vector with Link structs, no `setLink()` method
5. **ArticulatedBodyInertia**: No constructor from RigidBodyInertia - removed from example
6. **PluckerTransform**: No `operator*` for chaining - noted in example comments

### Documentation Updates

- README.md paths updated: `./build/examples/` not `./build/`
- Removed transform chaining from transforms.cpp (not implemented)
- Simplified dynamics example to match actual API

## Files Created/Modified

| File | Action | Lines | Purpose |
|------|--------|-------|---------|
| examples/CMakeLists.txt | Created | 20 | Build configuration |
| examples/basic_vectors.cpp | Created | 108 | Vector operations |
| examples/transforms.cpp | Created | 120 | Plücker transforms |
| examples/inertia.cpp | Created | 130 | Rigid body inertia |
| examples/dynamics.cpp | Created | 156 | Forward dynamics ABA |
| CMakeLists.txt | Modified | +2 | Add examples subdirectory |
| README.md | Modified | ~5 | Update example paths |

## Commits

- `8a0b11b`: docs(10-03): add compilable usage examples

## Example Output Samples

**example_vectors**: Correctly shows cross products, dot products, vector arithmetic
**example_transforms**: 90° Z rotation transforms X-axis to Y-axis as expected
**example_inertia**: F=ma=2N for 2kg mass, torque from offset COM visible
**example_dynamics**: ABA produces 10 rad/s² and 5 rad/s² for 1.0 and 0.5 N·m torques
