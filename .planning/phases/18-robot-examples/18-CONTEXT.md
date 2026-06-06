# Phase 18: Robot Examples - Context

**Gathered:** 2026-06-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Create two example executables demonstrating real-world robot physics: a 2-link Z-Z planar arm and a 3-link Z-Y-Z spatial RRR arm. Both examples must demonstrate forward dynamics (ABA) and inverse dynamics (RNEA) with gravity compensation, using FD+ID cross-validation to verify static equilibrium.

Requirements EX-01 and EX-02 define WHAT — this context captures HOW.

</domain>

<decisions>
## Implementation Decisions

### Link Parameter Model
- **D-01:** Use UR5-like parameters (mass, inertia, link length) adapted for Z-Z planar and Z-Y-Z spatial configurations. Parameter values from published UR5 datasheet (link masses, approximate link lengths).
- **D-02:** Gravity vector = `Vector3d(0, 0, -9.81)` — standard Z-down robotics convention. Both solvers accept gravity as a `Vector3d` parameter (ForwardDynamics::computeAccelerations, InverseDynamics::computeTorques).

### Gravity Verification Strategy
- **D-03:** FD+ID cross-validation pattern in both examples:
  1. Set up robot model with UR5-like parameters
  2. Set static pose (q arbitrary, qdot=0, qddot=0)
  3. Call ID with gravity = [0,0,-9.81] → produces gravity compensation torques
  4. Feed gravity torques into FD with same gravity → should produce zero acceleration
- **D-04:** This demonstrates static equilibrium: gravity torques exactly cancel gravitational forces, so no acceleration occurs.

### Code Organization
- **D-05:** Fully self-contained .cpp files — each example is standalone (no shared robot_utils.h). Matches existing example pattern (`examples/{basic_vectors,transforms,inertia,dynamics}.cpp`).
- **D-06:** File names: `examples/example_robot_2link.cpp` and `examples/example_robot_3link.cpp`.
- **D-07:** Executable names registered in `examples/CMakeLists.txt`: `example_robot_2link` and `example_robot_3link`, both linking `SpatialAlgebra Eigen3::Eigen`.

### Output Format
- **D-08:** Verbose output matching existing examples: section headers (`=== Title ===`), explanatory text before each computation, print intermediate results with physical interpretation, summary at end.

### Example Structure (both examples)
- **D-09:** Each example demonstrates:
  1. Robot model setup with realistic parameters (mass, length, inertia, joint axes)
  2. Forward dynamics (ABA): apply known torques, compute accelerations
  3. Inverse dynamics (RNEA): compute gravity torques for static pose
  4. Cross-validation: gravity torques → FD → verify zero acceleration
  5. Print physical interpretation of results

### Joint Configuration
- **D-10:** 2-link Z-Z planar arm: both joints rotate about local Z axis, arm operates in XY plane. Gravity [0,0,-9.81] creates non-zero gravity torques via COM moment arms about joint axes.
- **D-11:** 3-link Z-Y-Z spatial RRR arm: joint 1 about Z (waist yaw), joint 2 about Y (shoulder pitch), joint 3 about Z (elbow roll). Standard anthropomorphic arm configuration.

### Gravity Torque Sign Convention (for static equilibrium check)
- **D-12:** The cross-validation is: `ID(q, 0, 0, gravity) = tau_gravity` → `FD(tau_gravity, gravity) → qddot ≈ 0` (within numerical precision ~1e-12). This confirms the solvers are consistent: gravity torques from ID exactly cancel gravity forces in FD.

### Claude's Discretion
- Exact UR5 parameter values (masses, inertias, link lengths) to use
- Specific static pose(s) used for gravity demonstration
- Number of test cases per example (1-3 different poses/torques)
- Exact cout output format within the verbose style
- Include guards and Doxygen for example files (following existing pattern)
- CMake `add_executable` details in examples/CMakeLists.txt

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Example Patterns
- `examples/dynamics.cpp` — Existing 2-link ABA example (reference pattern for Phase 18)
- `examples/basic_vectors.cpp` — Existing verbose output pattern
- `examples/transforms.cpp` — Existing verbose output pattern
- `examples/inertia.cpp` — Existing verbose output pattern
- `examples/CMakeLists.txt` — Example build target registration pattern

### Solver APIs (Used in Examples)
- `include/ForwardDynamics.h:156` — `computeAccelerations(tau, gravity)` API
- `include/InverseDynamics.h:144` — `computeTorques(qddot, gravity)` API
- `include/RigidBodyInertia.h` — Inertia construction (mass, COM, inertia tensor)
- `include/PluckerTransform.h` — Link frame transform construction
- `include/MotionVector.h` — Joint motion axis (screw axis) construction

### Requirements
- `.planning/REQUIREMENTS.md` §EX-01, EX-02 — Requirement definitions
- `.planning/ROADMAP.md` §Phase 18 — Success criteria

### Prior Phase Context
- `.planning/phases/14-cr-02-bug-fix/14-CONTEXT.md` — ABA gravity parameter fix
- `.planning/phases/16-benchmark-infrastructure/16-CONTEXT.md` §D-01–D-04 — ModelFactory pattern (not used directly, but informs link construction)

### Build System
- `CMakeLists.txt` — Root build config
- `examples/CMakeLists.txt` — Target registration pattern

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `ForwardDynamics.h:156` — `computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity)` — FD solver with gravity support
- `InverseDynamics.h:144` — `computeTorques(const Eigen::VectorXd& qddot, const Vector3d& gravity)` — ID solver with gravity support
- `examples/dynamics.cpp` — Existing 2-link ABA example (demonstrates Link setup, solver usage, output pattern)
- `examples/CMakeLists.txt` — Known-good add_executable pattern with SpatialAlgebra link

### Established Patterns
- Self-contained `main()` with verbose `std::cout` output
- Doxygen block comments at file top explaining what the example demonstrates
- Using `std::cout << "=== Section ===" << std::endl;` for section headers
- `using namespace SpatialAlgebra;` + type aliases at file scope

### Integration Points
- `examples/CMakeLists.txt` — Add two new `add_executable` entries
- `examples/example_robot_2link.cpp` — New file (Z-Z planar arm)
- `examples/example_robot_3link.cpp` — New file (Z-Y-Z spatial arm)

</code_context>

<specifics>
## Specific Ideas

- Gravity cross-validation pattern: `tau_gravity = id.computeTorques(VectorXd::Zero(n), gravity)` → `fd.computeAccelerations(tau_gravity, gravity)` → verify `fd.links[i].qddot ≈ 0`
- For the 2-link example, demonstrate at least one pose with arm extended (horizontal) to show maximal gravity torque
- For the 3-link example, demonstrate at least one pose with arm extended horizontally to show joint 2 (shoulder) bearing the gravity load
- Use `Rotation::Identity()` for transforms where no rotation occurs between consecutive link frames (pure translation only)
- UR5 parameter sources: UR5 technical specification (mass: link 1 ~3.7 kg, link 2 ~8.4 kg, link 3 ~2.3 kg; approximate link lengths: 0.089m, 0.425m, 0.392m)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 18-robot-examples*
*Context gathered: 2026-06-06*
