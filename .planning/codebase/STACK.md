# Technology Stack

**Analysis Date:** 2026-05-17

## Languages

**Primary:**
- C++17 - Core library implementation in `include/` (12 headers) and `src/` (12 source files + `main.cpp`)
- C++17 - Test code in `tests/` (10 test files)

**Secondary:**
- Python 3 - Standalone RNEA implementation in `robot_dynamics/rnea.py`
- CMake - Build system language in `CMakeLists.txt` and `examples/CMakeLists.txt`

## Runtime

**Environment:**
- Native compiled C++ library (no interpreter runtime)
- macOS (Darwin) development environment (Apple clang via `/usr/bin/clang`)

**Package Manager:**
- None (CMake-based build with `find_package` for system dependencies)
- Homebrew for system dependencies (`brew install eigen`, `brew install googletest`)
- Lockfile: Not applicable

## Frameworks

**Core:**
- Eigen3 3.3+ (compatible with Eigen 3.4.0_1 currently installed) - Linear algebra backend for all matrix/vector operations, found via `find_package(Eigen3 REQUIRED NO_MODULE)` in `CMakeLists.txt:12`

**Testing:**
- Google Test (GTest, no version pin) - Test framework used by all 10 test executables, found via `find_package(GTest REQUIRED)` in `CMakeLists.txt:15`

**Build/Dev:**
- CMake 3.10+ - Build system configuration (`CMakeLists.txt:1`)
- Doxygen 1.12.0 - Documentation generation (`Doxyfile:1`)
- VSCode - IDE configuration in `.vscode/c_cpp_properties.json` and `.vscode/settings.json`

## Key Dependencies

**Critical:**
- Eigen3 (3.3+ / 3.4.0_1 installed) - Required for all linear algebra. Provides `Eigen::Matrix3d`, `Eigen::Vector3d`, `Eigen::Quaterniond`, `Eigen::AngleAxisd`. Installed via Homebrew at `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`.

**Infrastructure:**
- Google Test (no version pin) - Unit testing framework for all test executables. Installed via Homebrew (`brew install googletest`).
- Standard C++ Library - No additional runtime libraries beyond C++17 standard library and pthreads (OpenMP available in `LowerTriangular.h`).

## Configuration

**Environment:**
- No runtime environment variables required
- Pure library with compile-time configuration only
- Eigen 5.x compatibility issue documented: CMakeLists.txt version pin `Eigen3 3.3` may fail with Eigen 5.x from Homebrew; workaround is to remove version pin or set `-DEigen3_DIR`

**Build:**
- `CMakeLists.txt` - Main build configuration (CMake 3.10+, C++17, Eigen3+GTest linking, static library)
- `examples/CMakeLists.txt` - Example executables build config
- `Doxyfile` - Doxygen documentation generation config (2884 lines)
- `.vscode/c_cpp_properties.json` - VSCode IntelliSense configuration (cppStandard: c++17, compiler: clang)
- `.gitignore` - 41 entries ignoring build artifacts, docs output, and IDE files

## Platform Requirements

**Development:**
- CMake 3.10+
- C++17 compatible compiler (g++ or clang++, Apple clang via Xcode CLT)
- Eigen3 3.3+ (`brew install eigen`)
- Google Test (`brew install googletest`)

**Production:**
- No deployment target specified
- Builds to static library: `build/libSpatialAlgebra.a`
- Header-only style with separate implementation files
- No runtime dependencies beyond standard C++ library

## Build Artifacts

**Library:**
- `build/libSpatialAlgebra.a` - Static library from all `src/*.cpp` sources

**Test Executables:**
- `build/TestSpatialVector` - Spatial vector tests
- `build/TestPluckerTransform` - Plücker transform tests
- `build/TestRotation` - Rotation tests
- `build/TestLowerTriangular` - Lower triangular matrix tests
- `build/TestSpatialUtils` - Spatial utilities tests
- `build/TestRigidBodyInertia` - Rigid body inertia tests
- `build/TestArticulatedBodyInertia` - Articulated body inertia tests
- `build/TestForwardDynamics` - Forward dynamics (ABA) tests
- `build/TestInverseDynamics` - Inverse dynamics (RNEA) tests
- `build/TestDynamicsConsistency` - Dynamics consistency tests
- `build/TestSpatialOperations` - Spatial operations tests

**Examples:**
- `build/examples/example_vectors` - Basic vectors example
- `build/examples/example_transforms` - Transforms example
- `build/examples/example_inertia` - Inertia example
- `build/examples/example_dynamics` - Dynamics example (ABA)

**Documentation:**
- `docs/html/` - Generated Doxygen HTML documentation
- `docs/latex/` - Generated Doxygen LaTeX documentation

## Compiler Configuration

**Flags (default from CMake):**
- C++17 standard via `set(CMAKE_CXX_STANDARD 17)` in `CMakeLists.txt:4`
- No custom compiler flags (no `-O2`, `-Wall`, etc. explicitly set)
- Apple Clang via `/usr/bin/clang` (`compilerPath` in `.vscode/c_cpp_properties.json:13`)

---

*Stack analysis: 2026-05-17*
