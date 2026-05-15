# Codebase Structure

**Analysis Date:** 2026-05-15

## Directory Layout

```
SpatialAlgebra/
├── include/              # Header files (public API)
├── src/                  # Implementation files
├── tests/                # Test files
├── robot_dynamics/       # Python RNEA implementation
├── docs/                 # Generated documentation
│   ├── html/             # Doxygen HTML output
│   └── latex/            # Doxygen LaTeX output
├── build/                # CMake build artifacts (gitignored)
├── .vscode/              # VSCode configuration
├── .github/              # GitHub configuration
├── .planning/            # GSD planning artifacts
├── CMakeLists.txt        # Build configuration
├── Doxyfile              # Doxygen documentation config
├── AGENTS.md             # Agent instructions
└── README.md             # Project readme
```

## Directory Purposes

**include/:**
- Purpose: Public header files defining the library API
- Contains: 12 header files (.h)
- Key files:
  - `SpatialVector.h` - Base spatial vector class
  - `MotionVector.h`, `ForceVector.h` - Derived spatial vectors
  - `Rotation.h` - 3D rotation class
  - `PluckerTransform.h` - 6x6 spatial transform
  - `RigidBodyInertia.h`, `ArticulatedBodyInertia.h` - Inertia representations
  - `LowerTriangular.h` - Packed matrix storage
  - `SpatialOperations.h`, `SpatialUtils.h` - Utility functions

**src/:**
- Purpose: Implementation files for header declarations
- Contains: 10 source files (.cpp)
- Key files:
  - `SpatialVector.cpp` - Base class implementation
  - `PluckerTransform.cpp` - Transform operations (154 lines)
  - `Rotation.cpp` - Rotation conversions
  - `MotionVector.cpp`, `ForceVector.cpp` - Derived class implementations
  - `main.cpp` - Demo executable

**tests/:**
- Purpose: Unit test files
- Contains: 5 test files (.cpp)
- Key files:
  - `TestSpatialVector.cpp` - Basic assert-based tests (17 lines)
  - `TestPluckerTransform.cpp` - GTest-based tests (63 lines)
  - `TestArticulatedBodyInertia.cpp` - Empty stub
  - `TestRigidBodyInertia.cpp` - Empty stub
  - `TestSpatialOperations.cpp` - Empty stub

**robot_dynamics/:**
- Purpose: Python implementation of RNEA algorithm
- Contains: `rnea.py` (120 lines)
- Note: Standalone, NOT integrated with C++ library

**docs/:**
- Purpose: Generated Doxygen documentation
- Contains: html/, latex/ subdirectories
- Generated via: `doxygen Doxyfile`

## Key File Locations

**Entry Points:**
- `src/main.cpp`: Demo application showing library usage
- `include/SpatialVector.h`: Primary entry point for library users

**Configuration:**
- `CMakeLists.txt`: Build system configuration
- `Doxyfile`: Documentation generation settings
- `.vscode/c_cpp_properties.json`: VSCode IntelliSense config

**Core Logic:**
- `include/` directory: All class definitions
- `src/` directory: All implementations

**Testing:**
- `tests/` directory: All test files
- Only 2 of 5 test files have implementations

## Naming Conventions

**Files:**
- Headers: PascalCase.h (e.g., `SpatialVector.h`, `PluckerTransform.h`)
- Sources: PascalCase.cpp (e.g., `SpatialVector.cpp`, `PluckerTransform.cpp`)
- Tests: PascalCase with Test prefix (e.g., `TestSpatialVector.cpp`)

**Classes:**
- PascalCase (e.g., `SpatialVector`, `MotionVector`, `PluckerTransform`)

**Functions:**
- camelCase (e.g., `transformMotion`, `crossProductMotion`, `getAngular`)

**Variables:**
- camelCase for locals (e.g., `transformedAngular`, `newRotation`)
- Snake case for private members (e.g., `angular`, `linear`, `rotation`, `translation`)

**Type Aliases:**
- Lowercase abbreviations: `mv`, `fv`, `plux`, `rbi`, `abi`, `lt`

## Where to Add New Code

**New Spatial Vector Type:**
- Header: `include/NewVectorType.h`
- Implementation: `src/NewVectorType.cpp`
- Tests: `tests/TestNewVectorType.cpp`

**New Transform/Operation:**
- Add to existing class in `include/` and `src/`
- Or create new utility in `include/SpatialUtils.h`

**New Test:**
- Tests: `tests/Test<ClassName>.cpp`
- Register in `CMakeLists.txt`:
  ```cmake
  add_executable(TestNewFeature tests/TestNewFeature.cpp)
  target_link_libraries(TestNewFeature SpatialAlgebra GTest::GTest GTest::Main)
  add_test(NAME TestNewFeature COMMAND TestNewFeature)
  ```

**Utilities:**
- Free functions: `include/SpatialUtils.h`
- Static class methods: `include/SpatialOperations.h`

## Special Directories

**build/:**
- Purpose: CMake build artifacts
- Generated: Yes (by CMake)
- Committed: No (in .gitignore)

**docs/html/, docs/latex/:**
- Purpose: Generated documentation
- Generated: Yes (by Doxygen)
- Committed: No (in .gitignore)

**.vscode/:**
- Purpose: Editor configuration
- Generated: Manually
- Committed: No (in .gitignore)

**.github/:**
- Purpose: GitHub configuration
- Contains: `copilot-instructions.md`
- Committed: No (in .gitignore)

## File Statistics

**Headers:** 12 files, ~1,780 lines total
**Sources:** 10 files, ~467 lines total
**Tests:** 5 files, ~183 lines total (only 2 non-empty)
**Python:** 1 file, 120 lines

---

*Structure analysis: 2026-05-15*
