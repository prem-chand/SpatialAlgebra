# Codebase Structure

**Analysis Date:** 2026-05-17

## Directory Layout

```
SpatialAlgebra/
├── include/               # Public header files (library interface)
│   ├── SpatialVector.h    # 6D spatial vector base class
│   ├── MotionVector.h     # Twist type (covariant transforms)
│   ├── ForceVector.h      # Wrench type (contravariant transforms)
│   ├── Rotation.h         # 3D rotation (extends Eigen::Matrix3d)
│   ├── PluckerTransform.h # 6×6 Plücker coordinate transform
│   ├── RigidBodyInertia.h # Rigid body mass properties (all inline)
│   ├── ArticulatedBodyInertia.h # Composite inertia for ABA (all inline)
│   ├── LowerTriangular.h  # Packed lower-triangular matrix storage
│   ├── SpatialUtils.h     # Free functions: skew(), dot(), cross()
│   ├── SpatialOperations.h # Static utility class wrapper
│   ├── ForwardDynamics.h  # ABA forward dynamics solver + Link struct
│   └── InverseDynamics.h  # RNEA inverse dynamics solver + InverseDynamicsLink struct
│
├── src/                   # Library source files
│   ├── SpatialVector.cpp
│   ├── MotionVector.cpp
│   ├── ForceVector.cpp
│   ├── Rotation.cpp
│   ├── PluckerTransform.cpp
│   ├── RigidBodyInertia.cpp      # [Stub] All inline in header
│   ├── ArticulatedBodyInertia.cpp # [Stub] All inline in header
│   ├── LowerTriangular.cpp
│   ├── SpatialOperations.cpp
│   ├── ForwardDynamics.cpp
│   ├── InverseDynamics.cpp
│   └── main.cpp                  # Minimal library demo
│
├── tests/                 # Test files (GTest + bare assert)
│   ├── TestSpatialVector.cpp          # 609 lines — SpatialVector, MotionVector, ForceVector
│   ├── TestPluckerTransform.cpp       # 946 lines — PluckerTransform + RBI + ABI
│   ├── TestRotation.cpp               # Rotation class
│   ├── TestLowerTriangular.cpp        # LowerTriangular matrix
│   ├── TestSpatialUtils.cpp           # Spatial utility functions
│   ├── TestRigidBodyInertia.cpp       # Rigid body inertia
│   ├── TestArticulatedBodyInertia.cpp # Articulated body inertia
│   ├── TestForwardDynamics.cpp        # ABA forward dynamics
│   ├── TestInverseDynamics.cpp        # RNEA inverse dynamics
│   ├── TestDynamicsConsistency.cpp    # RNEA↔ABA cross-validation
│   └── TestInterfaceContracts.md      # Interface documentation (markdown, not compiled)
│
├── examples/             # Usage demonstration executables
│   ├── CMakeLists.txt    # Builds 4 example executables
│   ├── basic_vectors.cpp # Motion/force vector operations demo
│   ├── transforms.cpp    # Plücker transform demo
│   ├── inertia.cpp       # Rigid body inertia demo
│   └── dynamics.cpp      # Forward dynamics (ABA) demo
│
├── robot_dynamics/       # Standalone Python implementation
│   └── rnea.py           # RNEA using NumPy (not integrated with C++)
│
├── build/                # CMake build output (gitignored)
│   ├── libSpatialAlgebra.a   # Static library
│   ├── TestSpatialVector     # Test executables
│   ├── TestPluckerTransform
│   ├── ...                   # Other test executables
│   └── example_*             # Example executables
│
├── docs/                 # Generated Doxygen documentation
│   ├── html/             # HTML documentation
│   └── latex/            # LaTeX documentation
│
├── .planning/            # GSD project planning artifacts
│   ├── codebase/         # Codebase analysis documents (this directory)
│   ├── phases/           # Phase artifacts
│   ├── research/         # Research artifacts
│   ├── PROJECT.md        # Project definition
│   ├── ROADMAP.md        # Milestone plans
│   ├── REQUIREMENTS.md   # Requirements tracking
│   ├── STATE.md          # Current milestone state
│   └── config.json       # GSD configuration
│
├── .vscode/              # VSCode configuration
│   └── c_cpp_properties.json  # IntelliSense includes for Eigen3
│
├── .github/              # GitHub configuration
│
├── CMakeLists.txt        # Root build configuration
├── Doxyfile              # Doxygen documentation config
├── Doxyfile.bak          # Doxygen config backup
├── AGENTS.md             # Agent instructions for AI tools
├── REVIEW.md             # Code review checklist/documentation
├── v1.0-VERIFICATION.md  # v1.0 milestone verification
├── LICENSE               # License file
└── README.md             # Project readme
```

## Directory Purposes

**`include/`:**
- Purpose: Public library API — all headers users must include
- Contains: 12 header files, one per class/module
- Key files: `SpatialVector.h:68` (base class), `PluckerTransform.h:77` (transform), `LowerTriangular.h:75` (packed matrix), `ForwardDynamics.h:134` (ABA solver)
- No umbrella header: users include individual headers as needed

**`src/`:**
- Purpose: Library implementation — compiled into `libSpatialAlgebra.a`
- Contains: 12 `.cpp` files matching headers + `main.cpp` demo
- Key files: `PluckerTransform.cpp:250` lines (most complex implementation), `ForwardDynamics.cpp:158` lines (ABA algorithm), `InverseDynamics.cpp:141` lines (RNEA algorithm)
- Note: `RigidBodyInertia.cpp` and `ArticulatedBodyInertia.cpp` are empty stubs

**`tests/`:**
- Purpose: Unit tests and property-based tests
- Contains: 9 GTest test executables + 1 markdown contract doc
- Key files: `TestPluckerTransform.cpp:946` lines (most comprehensive tests), `TestSpatialVector.cpp:609` lines, `TestForwardDynamics.cpp`, `TestInverseDynamics.cpp`
- 2 test executables registered in `CMakeLists.txt` (`TestSpatialVector`, `TestPluckerTransform` at lines 30-31); 7 more added later (lines 67-132)
- All test files have their own `main()` calling `RUN_ALL_TESTS()`

**`examples/`:**
- Purpose: Usage demonstrations for all major library features
- Contains: 4 example source files, each produces its own executable
- Key files: `basic_vectors.cpp`, `transforms.cpp`, `inertia.cpp`, `dynamics.cpp`
- Built by `examples/CMakeLists.txt` which links each executable to `SpatialAlgebra` library

**`build/`:**
- Purpose: CMake build output directory
- Contains: Static library + test/example executables (gitignored)
- Build command: `cmake -B build && cmake --build build`
- Test command: `cmake --build build && cd build && ctest --output-on-failure`

**`robot_dynamics/`:**
- Purpose: Standalone Python implementation
- Contains: `rnea.py` — RNEA using NumPy with `RigidBodyParams` dataclass and `MultiBodySystem` class
- Not integrated with C++ library

**`.planning/`:**
- Purpose: GSD project management artifacts
- Contains: Codebase analysis, milestone context, phase plans, requirements, state tracking
- Key files: `STATE.md` (current milestone status), `PROJECT.md` (project definition), `ROADMAP.md` (planned phases/features)

## Key File Locations

**Entry Points:**
- `src/main.cpp`: Minimal library demo (creates vectors, transforms, inertias). Built as part of `SpatialAlgebra` library, not a separate executable.
- `tests/TestSpatialVector.cpp:605-608`: `main()` running GTest for vector tests
- `tests/TestPluckerTransform.cpp:942-945`: `main()` running GTest for transform tests
- Every test file (`tests/Test*.cpp`) has its own `main()` with `RUN_ALL_TESTS()`
- `examples/basic_vectors.cpp:24`: Example entry point
- `examples/transforms.cpp:30`: Example entry point
- `examples/inertia.cpp:30`: Example entry point
- `examples/dynamics.cpp:30`: Example entry point

**Configuration:**
- `CMakeLists.txt`: Root build configuration — C++17, Eigen3, GTest, library, tests, examples (135 lines)
- `examples/CMakeLists.txt`: Example build configuration (23 lines)
- `Doxyfile`: Doxygen documentation generation configuration
- `.vscode/c_cpp_properties.json`: VSCode IntelliSense includes path for Eigen3
- `AGENTS.md`: Agent instructions for AI tooling (C++17 build, test commands, architecture overview, coding conventions)

**Core Logic:**
- `include/SpatialVector.h`: Base class definition (177 lines of header, 65 lines impl)
- `include/MotionVector.h`: Motion vector specialization (161 lines header, 55 lines impl)
- `include/ForceVector.h`: Force vector specialization (164 lines header, 60 lines impl)
- `include/PluckerTransform.h`: Plücker transform (211 lines header, 250 lines impl) — most complex class
- `include/LowerTriangular.h`: Packed matrix (571 lines header, 55 lines impl) — largest header
- `include/Rotation.h`: 3D rotation (177 lines header, 71 lines impl)
- `include/ForwardDynamics.h`: ABA solver (191 lines header, 158 lines impl)
- `include/InverseDynamics.h`: RNEA solver (175 lines header, 141 lines impl)

**Testing:**
- `tests/TestSpatialVector.cpp`: 609 lines — comprehensive vector/motion/force tests
- `tests/TestPluckerTransform.cpp`: 946 lines — transform/RBI/ABI tests with property tests
- Other test files exist but vary in implementation status

## Naming Conventions

**Files:**
- PascalCase with descriptive names: `SpatialVector.h`, `PluckerTransform.cpp`, `ForwardDynamics.h`
- Test files prefixed with `Test`: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`
- Example files lowercase with underscores: `basic_vectors.cpp`, `transforms.cpp`

**Classes:**
- PascalCase: `SpatialVector`, `MotionVector`, `ForceVector`, `PluckerTransform`, `RigidBodyInertia`, `ArticulatedBodyInertia`, `LowerTriangular`, `ForwardDynamics`, `InverseDynamics`, `SpatialOperations`, `Rotation`
- Structs also PascalCase: `Link`, `InverseDynamicsLink`, `RigidBodyParams` (Python)

**Type Aliases (at namespace scope):**
- Lowercase abbreviations for concise notation: `mv` = MotionVector, `fv` = ForceVector, `plux` = PluckerTransform, `rbi` = RigidBodyInertia, `abi` = ArticulatedBodyInertia, `lt` = LowerTriangular

**Methods:**
- camelCase verbs: `transformMotion`, `crossProductMotion`, `getAngular`, `setFromAngleAxis`, `computeAccelerations`, `inverseTransformMotion`
- Getters prefixed with `get`: `getAngular()`, `getLinear()`, `getMass()`, `getCom()`, `getInertiaMatrixLT()`
- Operators as standard C++: `operator+`, `operator-`, `operator*`, `operator()`

**Members:**
- Private: lowercase (e.g., `angular`, `linear`, `rotation`, `translation`, `mass`, `com`, `Inertia`, `H`, `M`, `data`, `n`)
- Note: `ArticulatedBodyInertia` uses PascalCase members (`Inertia`, `H`, `M`) — inconsistent with other classes

**Namespaces:**
- `SpatialAlgebra` — single namespace for entire library

**Include Guards:**
- `#ifndef`/`#define`/`#endif` matching filename: `SPATIAL_VECTOR_H`, `PLUCKER_TRANSFORM_H`, `FORWARD_DYNAMICS_H`
- Exception: `LowerTriangular.h` uses `#pragma once`

**Source directory names:**
- All lowercase: `include/`, `src/`, `tests/`, `examples/`, `docs/`, `build/`, `robot_dynamics/`

## Where to Add New Code

**New Feature (new class):**
- Header: `include/NewClass.h`
- Source: `src/NewClass.cpp`
- Tests: `tests/TestNewClass.cpp`
- CMake: Add test executable in `CMakeLists.txt` following pattern (lines 30-31 or 67-132)
- Type alias: Add `using nc = NewClass;` at namespace scope in header

**New Method (existing class):**
- Declare in the matching header file (`include/ClassName.h`) inside `namespace SpatialAlgebra`
- Define in the matching source file (`src/ClassName.cpp`) or inline in header for trivial implementations
- Add Doxygen `@brief`, `@details`, `@param`, `@return` comments

**New Algorithm (dynamics variant):**
- Header: `include/NewAlgorithm.h` (model after `include/ForwardDynamics.h`)
- Source: `src/NewAlgorithm.cpp` (model after `src/ForwardDynamics.cpp`)
- Tests: `tests/TestNewAlgorithm.cpp` with GTest
- Register in `CMakeLists.txt` following the pattern at lines 94-102 (or any line range following the established add_executable → target_link_libraries → add_test sequence)

**New Example:**
- Source: `examples/new_example.cpp`
- Build: Add to `examples/CMakeLists.txt` following existing pattern (line 10-12)
- Link against `SpatialAlgebra` and `Eigen3::Eigen`

**New Test:**
- Source: `tests/TestNewFeature.cpp` with `main()` calling `RUN_ALL_TESTS()`
- CMake: Add in root `CMakeLists.txt` following pattern:
  ```
  add_executable(TestNewFeature tests/TestNewFeature.cpp)
  target_link_libraries(TestNewFeature SpatialAlgebra GTest::GTest GTest::Main)
  add_test(NAME TestNewFeature COMMAND TestNewFeature)
  ```
- Use GTest `TEST()` macros with `EXPECT_DOUBLE_EQ` / `EXPECT_NEAR`

**New Utility Function:**
- If closely related to existing class: add as method in that class's header/source
- If general spatial algebra utility: add as free function in `SpatialUtils.h` (inline in header)
- If static utility wrapper: add as static method in `SpatialOperations.h/cpp`

## Special Directories

**`build/`:**
- Purpose: Build output directory
- Generated: Yes (by CMake)
- Committed: No (gitignored)

**`docs/`:**
- Purpose: Generated Doxygen documentation
- Generated: Yes (`doxygen Doxyfile`)
- Committed: Yes (HTML and LaTeX output)

**`.planning/`:**
- Purpose: GSD project planning artifacts
- Generated: No (manually created by developers/AI agents)
- Committed: Yes

**`robot_dynamics/`:**
- Purpose: Standalone Python implementation, not part of C++ build
- Generated: No
- Committed: Yes

---

*Structure analysis: 2026-05-17*
