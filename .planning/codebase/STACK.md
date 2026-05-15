# Technology Stack

**Analysis Date:** 2026-05-15

## Languages

**Primary:**
- C++17 - Core library implementation in `include/` and `src/`

**Secondary:**
- Python 3 - Standalone RNEA implementation in `robot_dynamics/rnea.py`

## Runtime

**Environment:**
- Native compiled C++ library (no runtime)
- macOS (Darwin) development environment

**Package Manager:**
- None (CMake-based build)
- Homebrew for system dependencies (Eigen3, Google Test)

## Frameworks

**Core:**
- Eigen3 3.3+ - Linear algebra backend for all matrix/vector operations
  - `Eigen::Matrix3d` - 3x3 rotation matrices
  - `Eigen::Vector3d` - 3D vectors
  - `Eigen::Matrix<double, 6, 1>` - 6D spatial vectors
  - `Eigen::AngleAxisd` - Angle-axis rotation representation
  - `Eigen::Quaterniond` - Quaternion rotation representation

**Testing:**
- Google Test (GTest) - Test framework for `TestPluckerTransform.cpp`
- Basic `assert()` - Minimal testing in `TestSpatialVector.cpp`

**Build:**
- CMake 3.10+ - Build system configuration
- g++ (GCC) - C++ compiler (configured in `CMakeLists.txt:7`)

## Key Dependencies

**Critical:**
- Eigen3 3.3+ - Required for all linear algebra operations
  - Location: System-wide (Homebrew: `/usr/local/Cellar/eigen/`)
  - VSCode config references: `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`

**Infrastructure:**
- Google Test - Unit testing framework
  - Linked via `GTest::GTest` and `GTest::Main` targets

## Configuration

**Environment:**
- No runtime environment variables required
- Pure library with compile-time configuration

**Build:**
- `CMakeLists.txt` - Main build configuration
- `.vscode/c_cpp_properties.json` - VSCode IntelliSense configuration
- `Doxyfile` - Doxygen documentation generation

## Platform Requirements

**Development:**
- CMake 3.10+
- C++17 compatible compiler (g++, clang++)
- Eigen3 3.3+ (`brew install eigen`)
- Google Test (`brew install googletest`)

**Production:**
- Compiled static library: `build/libSpatialAlgebra.a`
- Header-only style with separate implementation files
- No runtime dependencies beyond standard C++ library

## Build Artifacts

**Library:**
- `build/libSpatialAlgebra.a` - Static library

**Executables:**
- `build/TestSpatialVector` - Spatial vector tests
- `build/TestPluckerTransform` - Plücker transform tests

**Documentation:**
- `docs/html/` - Generated Doxygen HTML documentation
- `docs/latex/` - Generated Doxygen LaTeX documentation

---

*Stack analysis: 2026-05-15*
