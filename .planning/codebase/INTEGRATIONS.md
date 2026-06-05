# External Integrations

**Analysis Date:** 2026-06-05

## APIs & External Services

**None.** This is a compiled C++17 static library with no runtime HTTP or network dependencies. The library has zero external integrations at runtime.

All external interactions are **build-time or development-time only**:
- **Eigen3** — Linear algebra library (compile-time dependency)
- **Google Test** — Unit testing framework (test-time dependency)
- **Doxygen** — Documentation generation (dev-time dependency)
- **Homebrew / apt** — System package managers for installing Eigen3 and GTest
- **FetchContent** — CMake module that downloads GTest from GitHub as fallback

## Data Storage

**Databases:**
- None. The library has no persistent data storage.

**File Storage:**
- Local filesystem only. Build artifacts (`build/`), documentation (`docs/`), and source code are stored on the local filesystem.

**Caching:**
- None. No caching layer used.

## Authentication & Identity

**Auth Provider:**
- None. The library has no user authentication, API keys, or identity management.

## Monitoring & Observability

**Error Tracking:**
- None. No Sentry, Datadog, or similar services.

**Logs:**
- No structured logging framework. All output goes to `std::cout` via `.print()` methods on each class (for debugging/display).

## CI/CD & Deployment

**Hosting:**
- Not applicable. This is a C++ static library, not a deployed service.
- Distribution is via source code (GitHub repository).

**CI Pipeline:**
- **GitHub Actions** — Defined in `.github/workflows/ci.yml`
  - Triggers: `push` and `pull_request` on `main` branch
  - **Matrix:**
    | OS | Compiler | Eigen Version |
    |----|----------|---------------|
    | ubuntu-latest | g++ | 3.4 |
    | ubuntu-latest | g++ | 5.0 |
    | ubuntu-latest | clang++ | 3.4 |
    | ubuntu-latest | clang++ | 5.0 |
    | macos-latest | g++ | 3.4 |
    | macos-latest | g++ | 5.0 |
    | macos-latest | clang++ | 3.4 |
    | macos-latest | clang++ | 5.0 |
  - **Coverage:** Ubuntu + g++ + Eigen 3.4 only, uploaded to Codecov
  - **Eigen install:** Built from source for non-standard versions (macOS 5.0, Linux 5.0)
  - **GTest install:** `brew install googletest` (macOS) or `apt-get install libgtest-dev` (Linux)
  - **Steps:** Checkout → Install dependencies → Configure → Build → Test → Upload coverage

**Coverage Reporting:**
- **Codecov** — Coverage upload via `codecov/codecov-action@v4` at `.github/workflows/ci.yml:77`
  - Only on ubuntu-latest + g++ + Eigen 3.4
  - `fail_ci_if_error: false`

## Environment Configuration

**Required env vars:**
- None. The library has zero runtime environment variable dependencies.

**Optional env vars:**
- `Eigen3_DIR` — CMake variable for finding Eigen3 installation directory (used in CI)
- `COVERAGE_FLAG` — CI-internal flag for enabling coverage build

**Secrets location:**
- None. No secrets required.

## Webhooks & Callbacks

**Incoming:**
- None. No webhook endpoints.

**Outgoing:**
- None. No webhook callbacks.

## Network Dependencies

**Build-time:**
- `https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz` — Eigen 3.4 source (CI fallback)
- `https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz` — Eigen 5.0.1 source (CI fallback)
- `https://github.com/google/googletest/archive/release-1.12.1.zip` — GTest FetchContent fallback
- `codecov/codecov-action@v4` — Coverage upload (CI only)

**Runtime:**
- None. No network calls at runtime.

## Vendored Dependencies

- **Eigen 5.0.1** — Full copy vendored at `eigen-5.0.1/`. Not used by the default build (which uses system-installed Eigen), but present for convenience/testing. Contains its own CI, tests, benchmarks, and build system.

## Development Tools

**VSCode Configuration:**
- `.vscode/c_cpp_properties.json` — IntelliSense configuration pointing to Eigen 3.4.0 headers at `/usr/local/Cellar/eigen/3.4.0_1/include/eigen3`
- `.vscode/settings.json` — Minimal: `files.associations` for Makefile and array

**GitHub Copilot:**
- `.github/copilot-instructions.md` — Instructions for Copilot code generation: detailed comments, type hints, modular code, security considerations

---

*Integration audit: 2026-06-05*
