# Phase 20: Test Model Library - Context

**Gathered:** 2026-06-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Extract all 11 existing test domains into a standalone, header-only test model library that depends only on Eigen3 (zero SpatialAlgebra dependency). Create an abstract adapter interface (`RobotSolver`) with a full robotics solver API, a `RobotModel` POD struct for defining kinematic chains in pure data form, and a SpatialAlgebra-to-adapter bridge. The library lives in `tests/test-models/` as an INTERFACE CMake target that existing tests consume.

Requirements TML-01 through TML-05 define WHAT — this context captures HOW.

</domain>

<decisions>
## Implementation Decisions

### Adapter Interface — RobotSolver (D-01 through D-06)
- **D-01:** Abstract base class with pure virtual methods (not CRTP, not C++20 concepts). Standard polymorphism — any solver that inherits and implements the interface can be swapped in. GTest-friendly, easy to mock.
- **D-02:** Full robotics solver API surface. Methods:
  - `setState(q, qdot)` — load joint configuration and velocity
  - `computeTorques(qddot, gravity)` → tau — inverse dynamics (RNEA)
  - `computeAccelerations(tau, gravity)` → qddot — forward dynamics (ABA)
  - `forwardKinematics()` — update internal joint transforms
  - `getJointTransform(idx)` → world-frame transform (POD, not PluckerTransform)
  - `computeMassMatrix()` → MatrixXd — joint-space inertia matrix H(q)
  - `computeGravityTorques(gravity)` → tau_gravity — gravity compensation torques
  - `computeJointSpaceJacobian(idx)` → MatrixXd — J for a given link
  - `getLinkCOM(idx)` → Vector3d — center of mass in world frame
  - `getDOF()` → int — number of degrees of freedom
- **D-03:** All public methods use dynamic Eigen types (`Eigen::VectorXd`, `Eigen::MatrixXd`). Templates on scalar or DOF count are NOT used — keeps the interface simple and avoids template proliferation across solver implementations.
- **D-04:** Hybrid state model. The `RobotModel` struct holds the reference joint configuration (q, qdot, qddot). The adapter maintains cached internal state (like Pinocchio's `Data` structure). When `setState()` is called, the adapter copies values and rebuilds internal cache as needed. The model reference state is accessible but the adapter's internal state is the computation source of truth.
- **D-05:** Builder pattern for construction. `RobotModel` is a pure data description of the kinematic chain. Each adapter provides a `build(const RobotModel&)` method that constructs the solver's internal representation and returns a ready-to-use solver. This separates model definition from solver instantiation — the same `RobotModel` can be passed to different adapters (SpatialAlgebra, Pinocchio, RBDL).
- **D-06:** Naming convention: `RobotSolver` for the abstract base class, `RobotModel` for the chain description POD, `JointSpec` for per-link parameter struct.

### Model Definition Format — RobotModel (D-07 through D-09)
- **D-07:** Pure data POD structs. Zero methods beyond simple accessors — no computation, no virtual dispatch, no inheritance. Fields use only Eigen3 types (`Eigen::Vector3d`, `Eigen::Matrix3d`, `Eigen::Matrix4d` for transforms). No SpatialAlgebra headers included.
- **D-08:** `RobotModel` contains an array of `JointSpec` structures. Each `JointSpec` describes one link:
  ```
  struct JointSpec {
      int parent;                          // Parent link index (-1 for base)
      Eigen::Matrix4d parentToJoint;       // Homogeneous transform (R|t)
      Eigen::Vector3d jointAxis;           // Revolute axis (ω), or prismatic direction (v)
      JointType type;                      // REVOLUTE, PRISMATIC, FIXED
      double mass;                         // Link mass
      Eigen::Vector3d com;                 // Center of mass in link frame
      Eigen::Matrix3d inertia;             // 3×3 rotational inertia at COM (dense, not LT)
      std::string name;                    // Optional link name for debugging
  };
  ```
  Transform stored as a 4×4 homogeneous matrix (not PluckerTransform) — zero SA dependency. Adapters convert to solver-specific transform types internally.
- **D-09:** Joint types: `enum class JointType { REVOLUTE, PRISMATIC, FIXED }`. Start with these three — extend in Phase 21 (Test Refinement) if needed for URDF-like chains. Revolute joints use the `jointAxis` as the rotation axis (ω). Prismatic joints use `jointAxis` as the translation direction (v).

### Directory Structure and Organization (D-10 through D-12)
- **D-10:** Library lives at `tests/test-models/`. Directory layout:
  ```
  tests/test-models/
  ├── CMakeLists.txt              # INTERFACE library: test_models → Eigen3::Eigen
  ├── robot_model.h               # RobotModel + JointSpec structs, JointType enum
  ├── robot_solver.h              # RobotSolver abstract base class
  ├── sa_adapter.h                # SpatialAlgebra adapter (links SpatialAlgebra)
  ├── sa_adapter.cpp              # Adapter implementation
  └── chains/                     # Model definitions, one header per test domain
      ├── spatial_vectors.h       # Single/multi-link vector test models
      ├── plucker_transforms.h    # Transform test models
      ├── rotation.h              # Rotation test models
      ├── lower_triangular.h      # Matrix operation test models
      ├── rigid_body_inertia.h    # RBI test models
      ├── articulated_body.h      # ABI test models
      ├── spatial_utils.h         # Utility test models
      ├── inverse_dynamics.h      # RNEA test models
      ├── forward_dynamics.h      # ABA test models
      ├── consistency.h           # Round-trip test models
      └── spatial_operations.h    # SpatialOps test models
  ```
- **D-11:** One header per test domain in `chains/`. Each header declares model-building functions that return `RobotModel` instances. Function naming: `makeThreeLinkSerialChain()`, `makeBranchingYConfiguration()`, etc. Functions are declared in headers, defined inline or in a `.ipp` file if templated.
- **D-12:** The SpatialAlgebra adapter (`sa_adapter.h/.cpp`) lives in the same directory but is the ONLY file that includes SpatialAlgebra headers. The INTERFACE library target does NOT link SpatialAlgebra — test executables that use the SA adapter add the `SpatialAlgebra` link separately. This creates a hard compilation boundary: including `robot_model.h` or `robot_solver.h` cannot accidentally pull in SA symbols.

### CMake Strategy (D-13 through D-15)
- **D-13:** INTERFACE library target named `test_models`. CMakeLists.txt:
  ```cmake
  add_library(test_models INTERFACE)
  target_include_directories(test_models INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_models INTERFACE Eigen3::Eigen)
  ```
  This is header-only — no compiled source. Test executables do `target_link_libraries(TestX PRIVATE test_models SpatialAlgebra)`.
- **D-14:** The SpatialAlgebra adapter is compiled separately via `add_library(sa_test_adapter sa_adapter.cpp)` in the same CMakeLists.txt, linking `SpatialAlgebra`. Tests that use it link both `test_models` and `sa_test_adapter`. Tests that only define models (no solver comparison) link just `test_models`.
- **D-15:** Existing tests are NOT immediately converted to use the new library — that's a Phase 21+ concern. Phase 20 creates the library infrastructure and adapter. A smoke test verifies the library compiles without SA dependency and the adapter correctly wraps existing solvers.

### the agent's Discretion
- Exact `RobotSolver` method signatures (const-correctness, return-by-value vs return-by-reference for large matrices)
- Internal caching strategy for the SA adapter (what gets cached, when to invalidate)
- Whether `forwardKinematics()` is called automatically before `getJointTransform()` or must be explicit
- Error handling convention (exceptions vs error codes for invalid state)
- Whether `RobotModel` stores DOF count explicitly or derives it from `JointSpec` array size

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase Requirements
- `.planning/REQUIREMENTS.md` — TML-01 through TML-05 (Test Model Library requirements)
- `.planning/ROADMAP.md` — Phase 20 goal and success criteria
- `.planning/PROJECT.md` — Project constraints: Eigen3 compatibility, C++17, existing code conventions

### Existing Model Patterns
- `benchmarks/common/model_factory.h` — Existing `JointConfig` struct and `ModelFactory` class (closest existing model-building pattern). The new `RobotModel`/`JointSpec` is a generalization of this to be solver-agnostic.
- `benchmarks/common/model_factory.cpp` — Factory implementation showing how `JointConfig` gets converted to `ForwardDynamicsLink`/`InverseDynamicsLink`. The SA adapter follows this pattern internally.

### Test Patterns
- `tests/TestDynamicsConsistency.cpp` — Round-trip test pattern (ID→tau→FD→qddot_check) used across all consistency tests. The adapter interface must support this two-solver pattern.
- `tests/TestForwardDynamics.cpp` — Direct dynamics tests showing model construction patterns (Link struct, PluckerTransform, RigidBodyInertia constructors).
- `tests/TestInterfaceContracts.md` — Existing test interface documentation — reference for understanding current test organization.

### Prior Phase Context
- `.planning/phases/14-cr-02-bug-fix/14-CONTEXT.md` — Test structure decisions from CR-02 fix (TDD pattern, test file organization, non-zero COM model construction)
- `.planning/phases/18-robot-examples/18-CONTEXT.md` — Self-contained example pattern (standalone .cpp files with model definitions)

### Codebase Structure
- `CMakeLists.txt` (root) — How tests are currently registered via `add_test` and linked against `SpatialAlgebra`
- `tests/` directory — 11 test executables, each with its own test domain

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `benchmarks/common/model_factory.h:JointConfig` — Existing per-link parameter struct with `axis`, `translation`, `mass`, `com` fields. The new `JointSpec` generalizes this by adding `parent`, `parentToJoint` (4×4 transform), `inertia` (3×3 dense), and `JointType` enum.
- `benchmarks/common/model_factory.h:ModelFactory` — Creates `ForwardDynamics` and `InverseDynamics` solver objects from `JointConfig` arrays. The SA adapter inherits from `RobotSolver` and uses this factory pattern internally.
- `include/ForwardDynamics.h:Link` (aka `ForwardDynamicsLink`) — Existing link struct with `parent`, `X` (PluckerTransform), `I` (RigidBodyInertia), `S` (MotionVector), `q`, `qdot`, `qddot`, `v`, `c`, `f`, `Ia`, `pa`. The SA adapter converts `JointSpec` arrays into these link structs.
- `include/InverseDynamics.h:InverseDynamicsLink` — Mirror struct for RNEA solver (different `Link` to avoid POSIX conflict). Same conversion needed in SA adapter.

### Established Patterns
- **Test model construction pattern** (all 11 test files): Build ID chain → computeTorques(qddot) → Build FD chain with same params → computeAccelerations(tau) → EXPECT_NEAR assertions. The adapter interface must support constructing two separate solver instances from the same `RobotModel`.
- **SpatialTransform convention**: PluckerTransform uses `[R, 0; -R*skew(r), R]` for motion and `[R^T, skew(r)*R^T; 0, R^T]` for force. The `JointSpec` uses a 4×4 homogeneous `[R|t; 0|1]` matrix — SA adapter converts to PluckerTransform via `Rotation(R.topLeftCorner<3,3>())` and `Vector3d(T.topRightCorner<3,1>())`.
- **Joint axis representation**: SpatialAlgebra uses `MotionVector(angular, linear)` for screw axes. The `JointSpec` uses `Vector3d` for the axis direction — SA adapter converts: REVOLUTE → `MotionVector(axis, Zero)`, PRISMATIC → `MotionVector(Zero, axis)`.
- **Doxygen on all declarations**: Existing codebase convention. All new types in `tests/test-models/` must follow this pattern.

### Integration Points
- `tests/CMakeLists.txt` — Add `add_subdirectory(test-models)` and update test targets to link `test_models`
- `tests/compile_smoke_test.cpp` — Add a smoke test that `#include "robot_model.h"` and `#include "robot_solver.h"` compile WITHOUT including any SpatialAlgebra headers (verifying D-10 zero-dependency guarantee)
- Existing test executables (`TestForwardDynamics`, `TestDynamicsConsistency`, etc.) — NOT modified in Phase 20; they continue using direct SA types

</code_context>

<specifics>
## Specific Ideas

- The `RobotModel` struct should support both serial chains (parent[i] = i-1 for i>0) and branching trees (multiple children per parent). Branching is needed for `BranchingYConfiguration` and the 3-link branching examples.
- For transforms stored in `JointSpec`, use the bottom row `[0 0 0 1]` convention (standard homogeneous matrix). This is trivial to convert to any solver's internal representation.
- The SA adapter should cache the `computeMassMatrix()` result until `setState()` is called again (like Pinocchio's `Data` pattern). Avoids recomputing H(q) when multiple queries are made for the same configuration.
- Consider a `makeDefault*()` naming convention for model factory functions in `chains/`: `makeDefaultThreeLinkSerial()` for the standard version, with overloads accepting custom mass/COM/inertia. Keeps the common case simple.
- `JointSpec::inertia` stores the 3×3 rotational inertia at the COM (dense `Eigen::Matrix3d`, not `LowerTriangular`). This is the most portable format — all solvers (SA, Pinocchio, RBDL) can consume 3×3 dense inertia matrices.

</specifics>

<deferred>
## Deferred Ideas

- **Converting existing tests to use the new library** — Phase 21 (Test Refinement) will update existing test executables to consume `test_models` and `sa_test_adapter`. Phase 20 only creates the infrastructure.
- **Pinocchio and RBDL adapters** — Phase 22 (PCC) and future phases implement non-SA adapters. The `RobotSolver` interface is designed to accommodate them.
- **Python bindings for model definitions** — Out of scope for v1.3 C++ milestone. Python comparison harness (Phase 23) duplicates models in Python, not via bindings.
- **URDF/SDF/MJCF parser integration** — Explicitly excluded. Models are programmatically constructed. URDF loading would bring a massive dependency chain.

</deferred>

---

*Phase: 20-test-model-library*
*Context gathered: 2026-06-17*
