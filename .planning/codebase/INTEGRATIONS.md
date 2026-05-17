# External Integrations

**Analysis Date:** 2026-05-17

## APIs & External Services

**None detected.** The library has zero runtime HTTP, gRPC, or IPC integrations. All computation is local, single-process C++.

## Data Storage

**Databases:**
- None detected. No database client libraries, ORM, or embedded databases (SQLite, LevelDB, etc.) are used.

**File Storage:**
- Local filesystem only. No cloud storage SDKs (S3, GCS, etc.) detected.

**Caching:**
- None detected. No Redis, Memcached, or in-memory cache libraries.

## Authentication & Identity

**Auth Provider:**
- None. No authentication, authorization, or identity management libraries detected. No user model exists.

## Monitoring & Observability

**Error Tracking:**
- None. No Sentry, Datadog, or similar libraries detected.

**Logging:**
- Console output only via `print()` methods on each class (e.g., `SpatialVector::print()` at `include/SpatialVector.h:172`)
- No structured logging framework (spdlog, glog, etc.)
- No log levels
- No log file output

**Metrics:**
- None detected.

**Tracing:**
- None detected.

## CI/CD & Deployment

**Hosting:**
- Not applicable. Static library with no deployment target.

**CI Pipeline:**
- None detected. No `.github/` workflow files exist.

**Package Publishing:**
- No package registry detected (no vcpkg, Conan, or Homebrew formula).

## Environment Configuration

**Required env vars:**
- None. The library has no runtime environment variable dependencies.

**Secrets location:**
- Not applicable. No secrets are used.

## Webhooks & Callbacks

**Incoming:**
- None.

**Outgoing:**
- None.

## External SDK / System Library Dependencies

**System Libraries:**
- Eigen3 (`find_package(Eigen3 REQUIRED NO_MODULE)`) — the only external library dependency
- Google Test (`find_package(GTest REQUIRED)`) — test-only dependency, not shipped with library

**Python Integration:**
- `robot_dynamics/rnea.py` contains a standalone Python 3 NumPy-based RNEA implementation
- **Not integrated** with the C++ library — no pybind11, Cython, or any Python binding mechanism
- Python dependency: `numpy` (used in `robot_dynamics/rnea.py:2`)
- Python stdlib: `dataclasses`, `typing` (stdlib only)

## Build-Time Integrations

**CMake External Dependencies:**
- `find_package(Eigen3)` — system-installed via Homebrew at `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`
- `find_package(GTest)` — system-installed via Homebrew
- No FetchContent, ExternalProject, or vcpkg/Conan package management

**Compiler/Platform:**
- Apple Clang via Xcode Command Line Tools (`/usr/bin/clang`)
- VSCode with `clangd` or C++ IntelliSense configured in `.vscode/c_cpp_properties.json`

## Document Generation

**Doxygen:**
- Configuration: `Doxyfile`
- Generate: `doxygen Doxyfile` produces `docs/html/` and `docs/latex/`
- No integration with CI or automated publishing

---

*Integration audit: 2026-05-17*
