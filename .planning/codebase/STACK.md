# Technology Stack

**Analysis Date:** 2026-06-05

## Languages

**Primary:**
- C++17 - Core library implementation in `include/` (13 headers) and `src/` (10 `.cpp` files)
- All library types, algorithms, and tests are written in C++17

**Secondary:**
- Python 3 - Standalone RNEA implementation in `robot_dynamics/rnea.py` (not integrated with C++ library)
- CMake - Build system configuration (`CMakeLists.txt`)
- YAML - CI pipeline definition (`.github/workflows/ci.yml`)

## Runtime

**Environment:**
- Native compiled C++ static library (`build/libSpatialAlgebra.a`)
- No runtime interpreter or VM
- Build artifacts are compiled executables and a static library

**Package Manager:**
- CMake 3.19+ (as build system, not a package manager)
- Homebrew (development dependency installation: `brew install eigen`, `brew install googletest`)
- FetchContent (CMake module for GTest fallback when not installed system-wide)
- Lockfile: Not applicable (C++ library with no language-level package manager)

## Frameworks

**Core:**
- Eigen3 3.4.x through 5.x - Linear algebra backend for all matrix/vector operations
  - Used via `#include <Eigen/Dense>` and `#include <Eigen/Geometry>`
  - `Rotation` class inherits from `Eigen::Matrix3d`
  - All `Vector3d`, `Matrix3d`, `VectorXd` types are Eigen typedefs
  - CMake: `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` at `CMakeLists.txt:12`

**Testing:**
- Google Test (GTest) 1.12.1+ - Unit testing framework for all test executables
  - Used via `#include <gtest/gtest.h>`
  - `TEST()` / `EXPECT_DOUBLE_EQ()` / `EXPECT_NEAR()` macros throughout
  - System install preferred; FetchContent fallback to `release-1.12.1.zip`
  - 10 test executables registered in CMakeLists.txt
  - Invoked via `ctest --output-on-failure`

**Build/Dev:**
- CMake 3.19+ - Build system configuration (`CMakeLists.txt`)
- Doxygen 1.12.0 - API documentation generation (`Doxyfile`)
- g++ / clang++ - C++17 compilers (CI tests both on Ubuntu and macOS)
- lcov - Code coverage generation (CI, Ubuntu + g++ + Eigen 3.4 only)

**CI:**
- GitHub Actions - CI pipeline defined in `.github/workflows/ci.yml`
- Testing matrix: 8 configurations (2 OS × 2 compilers × 2 Eigen versions)
- Codecov - Coverage upload via `codecov/codecov-action@v4`

## Key Dependencies

**Critical:**
- Eigen3 3.4+ / 5.x - Required at compile time. All linear algebra depends on it. No runtime dependency after compilation.
  - System install via Homebrew on macOS or `apt-get` on Ubuntu
  - Bundled copy at `eigen-5.0.1/` (vendored, not used by the build by default)
  - VSCode config points to `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3` for IntelliSense

**Development:**
- Google Test - Unit test framework, linked into all 10 test executables
  - Falls back to FetchContent download if not installed system-wide
- Doxygen - Documentation generation; not a build dependency, only for docs
- lcov - Code coverage instrumentation (optional, via `-DENABLE_COVERAGE=ON`)

## Configuration

**Build:**
- `CMakeLists.txt` (157 lines) - Main build configuration
  - `CMAKE_CXX_STANDARD 17` at line 4
  - Option: `ENABLE_COVERAGE` (OFF by default) at line 42
  - Examples built via `add_subdirectory(examples)` at line 157
- `examples/CMakeLists.txt` - Example executables configuration
- `.gitignore` - Ignores `build/`, `docs/` (generated), standard C++ artifacts

**Environment:**
- No runtime environment variables required
- No `.env` files present
- VSCode config in `.vscode/c_cpp_properties.json` specifies compiler and include paths

**Documentation:**
- `Doxyfile` (2884 lines, Doxygen 1.12.0) - Output to `docs/`
- Generated HTML at `docs/html/`, LaTeX at `docs/latex/`
- Run via `doxygen Doxyfile`

**CI:**
- `.github/workflows/ci.yml` (80 lines) - GitHub Actions workflow
- `.github/copilot-instructions.md` - GitHub Copilot code generation guidelines

## Platform Requirements

**Development:**
- macOS (Darwin) or Linux
- CMake 3.19+
- C++17 compatible compiler (g++ or clang++)
- Eigen3 3.4.x or 5.x (`brew install eigen` or `apt-get install libeigen3-dev`)
- Google Test (`brew install googletest` or `apt-get install libgtest-dev`)
- Doxygen (optional, for docs)

**Production:**
- Static library: `build/libSpatialAlgebra.a`
- No runtime dependencies beyond standard C++ library
- Link against any C++17 project with Eigen3 available at build time

## Build Artifacts

| Artifact | Source | Description |
|----------|--------|-------------|
| `build/libSpatialAlgebra.a` | All `src/*.cpp` | Core static library |
| `build/TestSpatialVector` | `tests/TestSpatialVector.cpp` | Spatial vector tests |
| `build/TestPluckerTransform` | `tests/TestPluckerTransform.cpp` | Plücker transform tests |
| `build/TestRotation` | `tests/TestRotation.cpp` | Rotation tests |
| `build/TestLowerTriangular` | `tests/TestLowerTriangular.cpp` | Lower triangular tests |
| `build/TestSpatialUtils` | `tests/TestSpatialUtils.cpp` | Spatial utils tests |
| `build/TestRigidBodyInertia` | `tests/TestRigidBodyInertia.cpp` | Rigid body inertia tests |
| `build/TestArticulatedBodyInertia` | `tests/TestArticulatedBodyInertia.cpp` | Articulated body inertia tests |
| `build/TestForwardDynamics` | `tests/TestForwardDynamics.cpp` | Forward dynamics tests |
| `build/TestInverseDynamics` | `tests/TestInverseDynamics.cpp` | Inverse dynamics tests |
| `build/TestDynamicsConsistency` | `tests/TestDynamicsConsistency.cpp` | Dynamics consistency tests |
| `build/TestSpatialOperations` | `tests/TestSpatialOperations.cpp` | Spatial operations tests |
| `build/examples/example_*` | `examples/*.cpp` | Usage examples |
| `docs/html/` | Doxygen output | API documentation (HTML) |
| `docs/latex/` | Doxygen output | API documentation (LaTeX) |

## Compiler & Standards

- **C++ Standard:** C++17 (`set(CMAKE_CXX_STANDARD 17)` in `CMakeLists.txt:4`)
- **Compiler:** System default (no hardcoded `g++` in CMake -- changed from earlier version at `CMakeLists.txt:6`)
- **VSCode IntelliSense:** clang++ on macOS (`/usr/bin/clang`), C++17 standard
- **No linter or formatter configured** (no `.clang-format`, `.clang-tidy` in project root; Eigen bundle has its own `.clang-format`)
- **OpenMP:** Not required (removed in v1.1, per README)

## License

- **Library:** GNU General Public License v3.0 (`LICENSE`)
- **Eigen (bundled):** MPL 2.0 with GPL exceptions (per `COPYING.MPL2`)

---

*Stack analysis: 2026-06-05*
