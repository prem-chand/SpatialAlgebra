# Codebase Structure

**Analysis Date:** 2026-06-05

## Directory Layout

```
SpatialAlgebra/
├── include/               # Public headers (13 files)
│   ├── SpatialAlgebra.h   # Umbrella header (includes all, in order)
│   ├── SpatialVector.h    # 6D spatial vector base class
│   ├── MotionVector.h     # Twist (motion) vector
│   ├── ForceVector.h       # Wrench (force) vector
│   ├── Rotation.h         # 3D rotation matrix extending Eigen::Matrix3d
│   ├── LowerTriangular.h  # Packed-storage lower triangular matrix
│   ├── SpatialUtils.h     # Free functions: skew(), dot(), cross()
│   ├── SpatialOperations.h # Static utility wrappers
│   ├── PluckerTransform.h # 6x6 Plücker coordinate transform
│   ├── RigidBodyInertia.h # Body mass properties
│   ├── ArticulatedBodyInertia.h # Articulated body inertia [I,H;H^T,M]
│   ├── ForwardDynamics.h  # Articulated Body Algorithm (ABA)
│   └── InverseDynamics.h  # Recursive Newton-Euler Algorithm (RNEA)
├── src/                   # Implementation files (10 files)
│   ├── SpatialVector.cpp
│   ├── MotionVector.cpp
│   ├── ForceVector.cpp
│   ├── Rotation.cpp
│   ├── LowerTriangular.cpp
│   ├── PluckerTransform.cpp
│   ├── SpatialOperations.cpp
│   ├── ForwardDynamics.cpp
│   ├── InverseDynamics.cpp
│   └── main.cpp           # Demo usage
├── tests/                 # Test files (13 files)
│   ├── TestSpatialVector.cpp
│   ├── TestPluckerTransform.cpp
│   ├── TestRotation.cpp
│   ├── TestLowerTriangular.cpp
│   ├── TestSpatialUtils.cpp
│   ├── TestRigidBodyInertia.cpp
│   ├── TestArticulatedBodyInertia.cpp
│   ├── TestForwardDynamics.cpp
│   ├── TestInverseDynamics.cpp
│   ├── TestDynamicsConsistency.cpp
│   ├── TestSpatialOperations.cpp
│   ├── compile_smoke_test.cpp
│   └── TestInterfaceContracts.md  # Interface contract documentation
├── examples/              # Usage examples (5 files)
│   ├── basic_vectors.cpp  # MotionVector / ForceVector operations
│   ├── transforms.cpp     # Plücker transform demonstration
│   ├── inertia.cpp        # Rigid body inertia operations
│   ├── dynamics.cpp       # Forward dynamics (ABA) demonstration
│   └── CMakeLists.txt     # Example build targets
├── robot_dynamics/        # Python implementation
│   └── rnea.py            # Standalone RNEA using NumPy (not integrated)
├── eigen-5.0.1/           # Vendored Eigen5 (alternative build)
├── docs/                  # Generated documentation
│   ├── html/              # Doxygen HTML output
│   └── latex/             # Doxygen LaTeX output
├── .planning/             # GSD planning artifacts
│   ├── codebase/          # Codebase maps (this file)
│   ├── research/          # Research documents
│   ├── phases/            # Phase plans and summaries
│   ├── milestones/        # Milestone artifacts
│   ├── ROADMAP.md
│   ├── PROJECT.md
│   ├── REQUIREMENTS.md
│   └── STATE.md
├── .vscode/               # VS Code configuration
│   ├── c_cpp_properties.json  # IntelliSense config (Eigen3 path)
│   └── settings.json
├── build/                 # Default build directory
├── build-eigen5/          # Eigen5 alternative build
├── CMakeLists.txt         # Root CMake build configuration
├── Doxyfile               # Doxygen configuration (1.12.0)
├── MATERIAL_CONVENTIONS.md # Mathematical conventions document
├── CLAUDE.md              # Claude agent instructions
├── AGENTS.md              # Agent instructions (build/test/arch)
├── README.md              # Project README
├── LICENSE                # License file
└── .gitignore
```

## Directory Purposes

**`include/`:**
- Purpose: All public API headers. Single-file-per-class, with `SpatialAlgebra.h` umbrella header.
- Contains: Header declarations with full Doxygen documentation. All classes in `namespace SpatialAlgebra`.
- Key files: `SpatialAlgebra.h` (umbrella), `SpatialVector.h` (foundation), `ForwardDynamics.h` (ABA), `InverseDynamics.h` (RNEA)

**`src/`:**
- Purpose: Implementation files for headers. One `.cpp` per header.
- Contains: Method definitions, algorithm implementations, demo `main.cpp`.
- Key files: `PluckerTransform.cpp` (most complex ~252 lines), `ForwardDynamics.cpp` (ABA ~210 lines), `InverseDynamics.cpp` (RNEA ~137 lines)

**`tests/`:**
- Purpose: GTest test suites. One test file per class/component, plus cross-cutting tests.
- Contains: GTest `TEST()` / `TEST_F()` test cases for all components.
- Key files: `TestPluckerTransform.cpp` (most comprehensive, ~940 lines), `TestSpatialVector.cpp` (~609 lines), `TestDynamicsConsistency.cpp` (cross-check)

**`examples/`:**
- Purpose: Runnable usage demonstrations, built as separate executables.
- Contains: Self-contained `main()` functions demonstrating library usage with descriptive output.
- Key files: `dynamics.cpp` (ABA workflow), `transforms.cpp` (Plücker transform chain)

**`eigen-5.0.1/`:**
- Purpose: Vendored copy of Eigen 5.0.1 for alternative build configuration (`build-eigen5/`).
- Contains: Full Eigen 5.0.1 source distribution.
- Committed: Yes (vendored dependency).

**`robot_dynamics/`:**
- Purpose: Python reference implementation for algorithm validation.
- Contains: NumPy-based RNEA (`rnea.py`). Not integrated with C++ library.
- Committed: Yes.

**`docs/`:**
- Purpose: Generated documentation artifacts (not hand-edited).
- Contains: Doxygen HTML and LaTeX output. Regenerated via `doxygen Doxyfile`.
- Generated: Yes. Committed: Yes.

**`.planning/`:**
- Purpose: GSD planning and project management artifacts.
- Contains: Codebase maps, phase plans, milestones, research notes, state tracking.
- Generated: Planning artifacts. Committed: Yes (by GSD conventions).

**`build/`, `build-eigen5/`:**
- Purpose: CMake build output directories.
- Contains: Compiled static library, test executables, example executables, CMake cache.
- Generated: Yes. Not committed (in `.gitignore`).

## Key File Locations

**Entry Points:**
- `src/main.cpp`: Demo entry point (shows vectors, transforms, inertia)
- `examples/basic_vectors.cpp`: Vector operations example
- `examples/transforms.cpp`: Plücker transform example
- `examples/inertia.cpp`: Inertia operations example
- `examples/dynamics.cpp`: Forward dynamics (ABA) example

**Configuration:**
- `CMakeLists.txt`: Root build configuration (C++17, Eigen3, GTest, coverage options)
- `Doxyfile`: Doxygen 1.12.0 documentation generation
- `.vscode/c_cpp_properties.json`: VS Code IntelliSense (Eigen 3.4.0 include path)
- `.vscode/settings.json`: VS Code editor settings

**Core Logic:**
- `include/SpatialVector.h` + `src/SpatialVector.cpp`: Spatial vector base (75 lines impl)
- `include/PluckerTransform.h` + `src/PluckerTransform.cpp`: Transform logic (252 lines)
- `include/ForwardDynamics.h` + `src/ForwardDynamics.cpp`: ABA solver (210 lines)
- `include/InverseDynamics.h` + `src/InverseDynamics.cpp`: RNEA solver (137 lines)
- `include/RigidBodyInertia.h`: Fully inline inertia (122 lines)
- `include/ArticulatedBodyInertia.h`: Fully inline articulated inertia (220 lines)
- `include/LowerTriangular.h` + `src/LowerTriangular.cpp`: Packed storage matrix (55 lines impl)
- `include/SpatialUtils.h`: Free function utilities, all inline (129 lines)

**Testing:**
- `tests/TestSpatialVector.cpp`: 11 tests for SpatialVector/MotionVector/ForceVector (609 lines)
- `tests/TestPluckerTransform.cpp`: ~33 tests across 7 test suites (940 lines)
- `tests/TestRotation.cpp`: Rotation tests
- `tests/TestLowerTriangular.cpp`: LowerTriangular tests
- `tests/TestSpatialUtils.cpp`: SpatialUtils free function tests
- `tests/TestRigidBodyInertia.cpp`: RigidBodyInertia tests
- `tests/TestArticulatedBodyInertia.cpp`: ArticulatedBodyInertia tests
- `tests/TestForwardDynamics.cpp`: ABA tests
- `tests/TestInverseDynamics.cpp`: RNEA tests
- `tests/TestDynamicsConsistency.cpp`: Cross-validation between inverse and forward dynamics
- `tests/TestSpatialOperations.cpp`: SpatialOperations tests
- `tests/compile_smoke_test.cpp`: Simple compilation verification

**Documentation:**
- `AGENTS.md`: Agent instructions (build, test, arch, conventions)
- `CLAUDE.md`: Claude-specific configuration
- `MATHEMATICAL_CONVENTIONS.md`: Mathematical notation reference
- `README.md`: Project overview and setup instructions
- `REVIEW.md`: Review documentation
- `v1.0-VERIFICATION.md`: v1.0 verification checklist

**Dependencies:**
- Eigen 3.3+ (system): `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` in `CMakeLists.txt:12`
- Google Test (system or FetchContent fallback): `find_package(GTest QUIET)` in `CMakeLists.txt:15-24`
- Eigen 5.0.1 (vendored): `eigen-5.0.1/` for alternative build

## Naming Conventions

**Files:**
- PascalCase for classes: `SpatialVector.h`, `PluckerTransform.h`, `ForwardDynamics.h`
- Test files prefixed with `Test`: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`
- Implementation files match header names: `SpatialVector.cpp` implements `SpatialVector.h`
- Purpose-prefixed in `src/`: `main.cpp` (demo), project classes match headers
- Examples use lowercase descriptive names: `basic_vectors.cpp`, `transforms.cpp`

**Directories:**
- Lowercase singular: `include/`, `src/`, `tests/`, `examples/`, `docs/`, `build/`
- `robot_dynamics/`: snake_case for Python package compatibility
- `eigen-5.0.1/`: versioned vendor directory

## Where to Add New Code

**New Feature / Class:**
- Header declaration: `include/<ClassName>.h`
- Implementation: `src/<ClassName>.cpp`
- Tests: `tests/Test<ClassName>.cpp`
- Example usage: `examples/<feature_name>.cpp`
- Register in `CMakeLists.txt`:
  - Source files are auto-globbed from `src/*.cpp` (line 27), but stubs are explicitly removed (lines 30-33)
  - Test executables must be added manually with `add_executable()` + `target_link_libraries()` + `add_test()`
  - Examples registered in `examples/CMakeLists.txt`

**New Utility / Helper:**
- If a free function, add to `include/SpatialUtils.h` (all inline) or create a new header if group is significant
- If a static method class, add to existing or new header

**New dynamics algorithm:**
- Follow the pattern of `ForwardDynamics` / `InverseDynamics`:
  - Header: define a `struct Link` for the kinematic chain
  - Header: define a solver class with `links` vector and compute method
  - Source: implement outward/inward pass methods
  - Tests: validate against known results or cross-validate with existing algorithm

**Type aliases:**
- Add namespace-level aliases in the class header: `using mv = MotionVector`, `using plux = PluckerTransform`
- Keep aliases short (2-4 characters), lowercase, abbreviations of full names

## Special Directories

**`eigen-5.0.1/`:**
- Purpose: Vendored Eigen 5.0.1 for independent build verification
- Generated: No (downloaded/distributed dependency)
- Committed: Yes

**`docs/`:**
- Purpose: Generated Doxygen HTML and LaTeX documentation
- Generated: Yes (by `doxygen Doxyfile`)
- Committed: Yes

**`build/`, `build-eigen5/`:**
- Purpose: CMake build artifacts (static library, test executables, examples)
- Generated: Yes
- Committed: No (in `.gitignore`)

---

*Structure analysis: 2026-06-05*
