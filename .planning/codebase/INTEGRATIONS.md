# External Integrations

**Analysis Date:** 2026-05-15

## Libraries & External Services

**Linear Algebra:**
- Eigen3 - Core dependency for all matrix/vector operations
  - SDK: Eigen header-only library
  - Integration: `#include <Eigen/Dense>`, `#include <Eigen/Geometry>`
  - Used in: All header files for vector/matrix types

## Data Storage

**Databases:**
- None - Pure computational library with no persistence

**File Storage:**
- None - In-memory computations only

**Caching:**
- None

## Authentication & Identity

**Auth Provider:**
- Not applicable - Local computational library

## Monitoring & Observability

**Error Handling:**
- C++ exceptions for invalid operations (`std::invalid_argument`, `std::out_of_range`)
- `assert()` for debug-mode checks
- No structured logging framework

**Logs:**
- `std::cout` for `print()` methods in all classes
- No formal logging infrastructure

## CI/CD & Deployment

**Hosting:**
- GitHub repository (local development only)

**CI Pipeline:**
- None configured
- No GitHub Actions workflows

**Build System:**
- CMake with manual build process:
  ```sh
  cmake -B build
  cmake --build build
  ```

## Environment Configuration

**Required env vars:**
- None

**Secrets location:**
- Not applicable - no external services

## Python Integration

**Standalone Implementation:**
- `robot_dynamics/rnea.py` - Recursive Newton-Euler Algorithm
  - Uses NumPy for array operations
  - NOT integrated with C++ library
  - Separate implementation for Python ecosystem

**Type Hints:**
- Python file uses `typing.List`, `typing.Optional`, `dataclasses`

## Webhooks & Callbacks

**Incoming:**
- None

**Outgoing:**
- None

## Type Aliases (Internal Integration)

**Convenience aliases in `namespace SpatialAlgebra`:**
- `mv` = `MotionVector`
- `fv` = `ForceVector`
- `plux` = `PluckerTransform`
- `rbi` = `RigidBodyInertia`
- `abi` = `ArticulatedBodyInertia`
- `lt` = `LowerTriangular`

## Documentation Generation

**Doxygen:**
- Config: `Doxyfile`
- Command: `doxygen Doxyfile`
- Output: `docs/html/`, `docs/latex/`

---

*Integration audit: 2026-05-15*
