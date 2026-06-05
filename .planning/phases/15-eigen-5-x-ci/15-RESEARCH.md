# Phase 15: Eigen 5.x CI — Research

**Researched:** 2026-06-05
**Domain:** Build system (CMake), CI/CD (GitHub Actions), C++ library version compatibility (Eigen 3.4.x ↔ 5.x)
**Confidence:** HIGH

## Summary

This phase adds Eigen 5.x to the CI build matrix while maintaining backward compatibility with Eigen 3.4.x. The work spans three areas: (1) CMakeLists.txt version range syntax change to accept both 3.4.x and 5.x, (2) CI matrix expansion from 4 jobs (2 OS × 2 compilers) to 8 jobs (2 OS × 2 compilers × 2 Eigen versions), and (3) optional preprocessor guards for any Eigen 5.x API differences discovered during compilation.

**Primary recommendation:** Use `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` version-range syntax (requires bumping CMake minimum to 3.19), expand the CI matrix with an `eigen-version` dimension, install Eigen 3.4 via system packages and Eigen 5.0.1 from source tarball, and let compilation errors (if any) determine the need for `#if EIGEN_MAJOR_VERSION >= 5` guards. Codebase analysis suggests zero source code changes are likely — this library uses only stable Eigen APIs (`Matrix3d::cross()`, `transpose()`, `block<>()`, `AngleAxisd`, `Quaterniond`, `Identity`, `eval()`).

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| CI-01 | Eigen 5.x compatible via version range syntax `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` | Verified syntax from official Eigen documentation. Requires CMake 3.19+. CMakeLists.txt currently requires 3.10 — must bump. |
| CI-02 | Eigen 5.x added to CI build matrix | Verified standard GitHub Actions matrix expansion pattern. 8 jobs: 2 OS × 2 compilers × 2 Eigen versions. Install strategy documented per OS. |
</phase_requirements>

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Use `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` version range syntax. Requires CMake 3.19+ (already satisfied by CI runners). Declares intent clearly while allowing 3.4.x through 5.x.
- **D-02:** Add `eigen-version` dimension to the CI matrix, producing 8 jobs (2 OS × 2 compilers × 2 Eigen versions: `3.4` and `5.0`). Explicit testing of both versions on all platforms.
- **D-03:** Install Eigen versions explicitly per job:
  - macOS: `brew install eigen@3.4` (or from source if tap unavailable) and `brew install eigen` (latest = 5.x)
  - Ubuntu: `apt install libeigen3-dev` (3.4.x) and build Eigen 5.0.1 from source tarball
  - Use `-DEigen3_DIR` to point CMake at the correct installation
- **D-04:** Fix forward-compatible code using `#if EIGEN_VERSION` preprocessor guards where 3.4.x and 5.x need different code paths. Clean builds under both versions — no warning suppression.
- **D-05:** Keep single coverage run on `ubuntu-latest` + `g++` with Eigen 3.4 (the existing config). Coverage reporting does not extend to Eigen 5.x jobs. Add CodeCov `flags` to distinguish coverage source if desired.

### The Agent's Discretion

- Specific Eigen 5.x API changes that need version-guarded fixes (to be discovered during compilation)
- Download/checksum details for Eigen 5.0.1 source tarball in CI
- Exact `apt` and `brew` version command details

### Deferred Ideas (OUT OF SCOPE)

None — discussion stayed within phase scope.
</user_constraints>

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| CMake version range syntax | Build (CMakeLists.txt) | — | Single-line change to `find_package()` call, owned by build system |
| CI matrix expansion | CI/CD (.github/workflows/) | — | Workflow YAML change — no build system or code involvement |
| Eigen source install in CI | CI/CD | Build | Download/configure/install tarball in CI step; build system only sees installed path via `-DEigen3_DIR` |
| Preprocessor guards for API diffs | Source code (headers) | — | Only needed if compilation fails; guards in header files using Eigen's version macros |
| Coverage preservation | CI/CD | — | Conditional step gated on `matrix.eigen-version == '3.4'` in existing coverage upload |

## Standard Stack

### Core

| Component | Version | Purpose | Why Standard |
|-----------|---------|---------|--------------|
| GitHub Actions | — | CI/CD orchestration | Already used; no migration needed |
| CMake version range | 3.19+ syntax | Accept Eigen 3.4.x through 5.x | Eigen 3.4.1+ CMake guide documents `3.4...5` syntax. Requires CMake 3.19+. |

### Eigen Version Detection (Preprocessor)

| Macro | Eigen 3.4 Value | Eigen 5.x Value | Notes |
|-------|-----------------|------------------|-------|
| `EIGEN_WORLD_VERSION` | `3` | `3` | Always `3` — for posterity per Eigen 5.0 release notes |
| `EIGEN_MAJOR_VERSION` | `4` | `5` | **The key discriminator** |
| `EIGEN_MINOR_VERSION` | `0` | `0` | Minor version |
| `EIGEN_PATCH_VERSION` | — | patch level | Set for patch releases (Eigen adopted semver with 5.0) |

**Guard pattern for code that needs different paths:**
```cpp
#include <Eigen/Core>  // defines EIGEN_MAJOR_VERSION
#if EIGEN_MAJOR_VERSION >= 5
    // Eigen 5.x path
#else
    // Eigen 3.4.x path
#endif
```

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Source build Eigen 5.x in CI | Homebrew `brew install eigen` on macOS, `apt` on Ubuntu | CI runners have Eigen 3.4 from apt. Homebrew latest is 5.x. Source build gives deterministic version control. |
| `FetchContent` for Eigen 5.x | Source tarball download + `cmake --install` | `FetchContent` integrates into your own build tree; separate install prefix is cleaner for `find_package` with `-DEigen3_DIR` |

## Package Legitimacy Audit

> No external packages are installed in this phase — only Eigen itself, which is a well-known reputable library (10+ years, 500+ executables test suite, MPL 2.0). No slopcheck needed.

| Package | Registry | Age | Downloads | Source Repo | Disposition |
|---------|----------|-----|-----------|-------------|-------------|
| Eigen 5.0.1 | Source tarball, GitLab | 8 mos | N/A (headers-only) | gitlab.com/libeigen/eigen | Approved — no registry risk |

## Architecture Patterns

### CI Matrix Structure
```
Current (4 jobs):          New (8 jobs):
os: [ubuntu, macos]        os: [ubuntu, macos]
compiler: [g++, clang++]   compiler: [g++, clang++]
                           eigen-version: [3.4, 5.0]
```

### Pattern 1: Matrix Expansion with `eigen-version`
**What:** Add a third dimension to the existing 2D build matrix.

**When to use:** Any time a dependency version needs cross-product testing with OS and compiler.

**CI workflow structure:**
```yaml
jobs:
  build:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        compiler: [g++, clang++]
        eigen-version: [3.4, 5.0]
    
    runs-on: ${{ matrix.os }}
    
    steps:
    - uses: actions/checkout@v4

    - name: Install dependencies (macOS) — Eigen 3.4
      if: runner.os == 'macOS' && matrix.eigen-version == '3.4'
      run: |
        brew update
        brew install eigen  # This is 5.x — see note below

    - name: Install Eigen 5.0.1 from source (macOS)
      if: runner.os == 'macOS' && matrix.eigen-version == '5.0'
      run: |
        brew update
        brew install eigen  # Homebrew now ships 5.x; verify version
        # OR build from source

    - name: Install dependencies (Ubuntu) — Eigen 3.4
      if: runner.os == 'Linux' && matrix.eigen-version == '3.4'
      run: |
        sudo apt-get update
        sudo apt-get install -y libeigen3-dev ...

    - name: Install Eigen 5.0.1 from source (Ubuntu)
      if: runner.os == 'Linux' && matrix.eigen-version == '5.0'
      run: |
        curl -L https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz | tar xz
        cd eigen-5.0.1 && cmake -B build -DCMAKE_INSTALL_PREFIX=$HOME/eigen-5.0.1
        cmake --build build && cmake --install build

    - name: Configure
      run: |
        EIGEN_FLAG=""
        if [ "${{ matrix.eigen-version }}" = "5.0" ]; then
          EIGEN_FLAG="-DEigen3_DIR=$HOME/eigen-5.0.1/share/eigen3/cmake"
        fi
        cmake -B build -DCMAKE_CXX_COMPILER=${{ matrix.compiler }} $EIGEN_FLAG ...
```

### Coverate Preservation Pattern
```yaml
    - name: Upload coverage
      if: matrix.os == 'ubuntu-latest' && matrix.compiler == 'g++' && matrix.eigen-version == '3.4'
      uses: codecov/codecov-action@v4
      with:
        directory: ./build
        fail_ci_if_error: false
```

### Anti-Patterns to Avoid
- **Using `exclude:` to prune Eigen 5.0 from certain OS/compiler combos:** All 8 combinations are valid and informative. Pruning risks missing a platform-specific Eigen 5.x issue.
- **Setting `fail-fast: true`:** One Eigen 5.0 failure would cancel all Eigen 3.4 jobs. Keep the existing default or explicitly set `fail-fast: false`.
- **Hardcoding version pin in CMakeLists.txt:** The whole point of the version range is to accept both. Don't pin to `3.4` or `5.0` in `find_package`.
- **Mixing Eigen versions in the same build directory:** Always use separate build dirs or clean when switching. The `-DEigen3_DIR` approach ensures deterministic selection.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Eigen version detection in C++ | Custom version-detection logic | `#include <Eigen/Core>` provides `EIGEN_MAJOR_VERSION` / `EIGEN_WORLD_VERSION` | Eigen ships these macros; they are tested, stable, and portable. |
| CMake version range upper bound | `if(VERSION ...)` logic | `find_package(Eigen3 3.4...5 ...)` | CMake's `...` syntax handles the upper-bound exclusion semantics correctly. |
| Eigen tarball download and extraction | Custom download script | `curl -L ... | tar xz` in CI step, or use `actions/cache` | Simple, auditable, and CI-native. No extra actions required. |

## Common Pitfalls

### Pitfall 1: `find_package(Eigen3 3.3 REQUIRED)` rejects Eigen 5.x
**What goes wrong:** Eigen 5.x installed but `find_package(Eigen3 3.3 REQUIRED)` fails because CMake interprets `3.3` as ">=3.3.0, <3.4.0" (COMpatible semantic), and Eigen 5.0 doesn't match.
**Why it happens:** CMake's version-matching treats a bare version as an exact upper-bound when the package uses proper versioning. Eigen 5.0 ships `Eigen3ConfigVersion.cmake` that conforms to semver.
**How to avoid:** Use version range syntax `find_package(Eigen3 3.4...5 REQUIRED)` instead of bare version.
**Warning signs:** CI step fails with `Could not find a configuration file for package "Eigen3" that is compatible with requested version "3.3"`.

### Pitfall 2: Eigen 5.0.1 vs 5.0.0 version incompatibility in CMake
**What goes wrong:** `find_package(Eigen3 5.0.0 REQUIRED)` fails with Eigen 5.0.1 installed.
**Why it happens:** Eigen 5.0.1's `Eigen3ConfigVersion.cmake` has a version-compatibility check that makes `5.0.1` incompatible with `5.0.0` as an exact request.
**How to avoid:** Don't request exact Eigen versions. Use the version range `3.4...5` which avoids pinning to a specific patch.
**Warning signs:** CI step fails when Eigen 5.0.1 is installed but `find_package` requests `5.0.0`.

### Pitfall 3: macOS `brew install eigen` installs 5.x, not 3.4.x
**What goes wrong:** Both the "3.4" and "5.0" macOS matrix jobs install the same Eigen because Homebrew's `eigen` formula now points to 5.x.
**Why it happens:** Homebrew updated `eigen` formula to Eigen 5.0.x. The old `eigen@3.4` tap may not exist or be outdated.
**How to avoid:** For the Eigen 3.4 job on macOS, build from source tarball the same way as Ubuntu, or use `brew extract` to pin 3.4. Simplest: build both versions from source on macOS for deterministic results.
**Workaround in CONTEXT.md:** D-03 already says "or from source if tap unavailable," which is the correct fallback. Source build is recommended for all Eigen 3.4 CI jobs to avoid this ambiguity.

### Pitfall 4: Eigen 5.x requires C++14, but project uses C++17
**What goes wrong:** Nothing — already on C++17. Eigen 5.0 says it's the last release to support C++14. This project already sets `CMAKE_CXX_STANDARD 17`.
**How to avoid:** No action needed.
**Warning signs:** Eigen 5.x compiler errors about C++14 features. (Not expected here.)

### Pitfall 5: `hasNaN()` / `array().isInf().any()` patterns change in Eigen 5
**What goes wrong:** These are used in `SpatialVector.cpp`, `RigidBodyInertia.h`, `ArticulatedBodyInertia.h` for debug assertions. If Eigen changed internal NaN/Inf handling, these could produce different results.
**Why it happens:** Eigen 5 overhauled vectorization and packet math.
**How to avoid:** Compile-test only. These are `#ifndef NDEBUG` guarded checks. If Eigen 5's scalar comparison behavior change (mentioned in release notes — scalar comparisons now return `Scalar(1)` instead of all-bits-set) affects expression templates, these should still work because they use `any()` and `isInf()` methods which return bool.
**Warning signs:** Compilation errors in `SpatialVector.cpp` or `RigidBodyInertia.h` around `hasNaN()`, `isInf()`, `any()`.

## Code Examples

### CMakeLists.txt — Version range change
```cmake
# Current (line 12):
find_package(Eigen3 REQUIRED NO_MODULE)

# New:
find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)
```

Also bump `cmake_minimum_required`:
```cmake
# Current (line 1):
cmake_minimum_required(VERSION 3.10)

# New — required for version range syntax:
cmake_minimum_required(VERSION 3.19)
```
[VERIFIED: Official Eigen CMake Guide at libeigen.gitlab.io/eigen/docs-5.0/TopicCMakeGuide.html]

### CI Workflow — Full expanded matrix
```yaml
jobs:
  build:
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
        compiler: [g++, clang++]
        eigen-version: [3.4, 5.0]
    
    runs-on: ${{ matrix.os }}

    steps:
    - uses: actions/checkout@v4

    - name: Install Eigen 3.4 (Ubuntu)
      if: runner.os == 'Linux' && matrix.eigen-version == '3.4'
      run: |
        sudo apt-get update
        sudo apt-get install -y libeigen3-dev libgtest-dev lcov

    - name: Build and install Eigen 5.0.1 (Linux)
      if: runner.os == 'Linux' && matrix.eigen-version == '5.0'
      run: |
        curl -sL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz | tar xz
        cmake -B eigen-5.0.1/build -S eigen-5.0.1 -DCMAKE_INSTALL_PREFIX=$HOME/eigen-5.0.1
        cmake --build eigen-5.0.1/build
        cmake --install eigen-5.0.1/build
        sudo apt-get update
        sudo apt-get install -y libgtest-dev lcov

    - name: Install Eigen 3.4 (macOS) — from source
      if: runner.os == 'macOS' && matrix.eigen-version == '3.4'
      run: |
        # Build 3.4 from source since Homebrew may not have it
        curl -sL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz | tar xz
        cmake -B eigen-3.4.0/build -S eigen-3.4.0 -DCMAKE_INSTALL_PREFIX=$HOME/eigen-3.4.0
        cmake --build eigen-3.4.0/build
        cmake --install eigen-3.4.0/build
        brew install googletest

    - name: Install Eigen 5.0.1 (macOS) — from source
      if: runner.os == 'macOS' && matrix.eigen-version == '5.0'
      run: |
        curl -sL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz | tar xz
        cmake -B eigen-5.0.1/build -S eigen-5.0.1 -DCMAKE_INSTALL_PREFIX=$HOME/eigen-5.0.1
        cmake --build eigen-5.0.1/build
        cmake --install eigen-5.0.1/build
        brew install googletest

    - name: Configure
      run: |
        EIGEN_DIR_FLAG=""
        if [ "${{ matrix.eigen-version }}" = "5.0" ]; then
          # Ubuntu and macOS both use $HOME prefix for source build
          EIGEN_DIR_FLAG="-DEigen3_DIR=$HOME/eigen-5.0.1/share/eigen3/cmake"
        else
          EIGEN_DIR_FLAG="-DEigen3_DIR=$HOME/eigen-3.4.0/share/eigen3/cmake"
        fi
        # Coverage only on ubuntu + g++ + eigen 3.4
        COVERAGE_FLAG="${{ (matrix.os == 'ubuntu-latest' && matrix.compiler == 'g++' && matrix.eigen-version == '3.4') && '-DENABLE_COVERAGE=ON' || '' }}"
        cmake -B build -DCMAKE_CXX_COMPILER=${{ matrix.compiler }} $EIGEN_DIR_FLAG $COVERAGE_FLAG

    - name: Build
      run: cmake --build build

    - name: Test
      run: cd build && ctest --output-on-failure

    - name: Upload coverage
      if: matrix.os == 'ubuntu-latest' && matrix.compiler == 'g++' && matrix.eigen-version == '3.4'
      uses: codecov/codecov-action@v4
      with:
        directory: ./build
        fail_ci_if_error: false
```

### Preprocessor guard pattern (if API differences found)
```cpp
#include <Eigen/Core>

class Rotation : public Eigen::Matrix3d {
public:
    // Existing constructors — no change needed for Eigen 5.x
    
    // If Eigen 5.x deprecates some conversion, guard it:
    Vector3d operator*(const Vector3d &vector) const
    {
#if EIGEN_MAJOR_VERSION >= 5
        // Eigen 5.x path (if different)
        return Eigen::Matrix3d(*this) * vector;
#else
        return static_cast<const Eigen::Matrix3d&>(*this) * vector;
#endif
    }
};
```

**Key insight:** The codebase already uses explicit `static_cast<const Eigen::Matrix3d&>(*this)` for all Eigen conversions, which is the safe pattern for both versions. No changes are expected.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` (version pin) | `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` (version range) | Eigen 3.4.1+ (2025) | Accepts both 3.4.x and 5.x. Requires CMake 3.19+. |
| Eigen `WORLD.MAJOR.MINOR` versioning | Eigen MAJOR.MINOR.PATCH (semver) | Eigen 5.0 (2025-09-30) | CMake version range syntax `3.4...5` means `>=3.4, <6` — fine regardless. |
| Source install via `make install` | `cmake --install` | CMake 3.15+ | Cleaner CI: `cmake --build build && cmake --install build` |

**Deprecated/outdated:**
- `EIGEN_HAS_CXX11` macro: removed in Eigen 5.x. This codebase does not use it.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Ubuntu `apt install libeigen3-dev` provides Eigen 3.4.x (not 5.x) | CI Matrix | Low — `apt-cache show libeigen3-dev` in a CI run will confirm. If wrong, switch to source install for ubuntu 3.4 too. |
| A2 | `brew install eigen` on macOS-latest will install Eigen 5.x | CI Matrix | Medium — `brew info eigen` at research time showed latest is 5.0.1 for Arch Linux Homebrew; macOS may differ. Build both from source for determinism. |
| A3 | No source code changes needed — Eigen 5.x API is backward-compatible for the APIs used | Preprocessor Guards | Medium — the only way to confirm is to compile. If compilation errors occur, add `#if EIGEN_MAJOR_VERSION >= 5` guards. The codebase uses very basic Eigen APIs (Vector3d, Matrix3d, block<>, cross, transpose, Identity, AngleAxisd, Quaterniond, eval()) which are all core and stable. |

## Open Questions

1. **What does `apt show libeigen3-dev` return on ubuntu-latest?**
   - What we know: Previous CI logs show it installs Eigen 3.4.0 on ubuntu-latest.
   - Recommendation: Confirm by checking CI log or running `apt-cache policy libeigen3-dev` in a test CI run.

2. **Does Homebrew's `eigen` formula require any special handling to work with Eigen 5.x?**
   - What we know: AGENTS.md documents that Eigen 5.0.1 installed via Homebrew may fail `find_package(Eigen3 3.3 REQUIRED)` because of CMake version range changes. This is the exact issue D-01 fixes.
   - Recommendation: The version range syntax fix also resolves Homebrew installations.

## Environment Availability

> Skipping this section — this phase is purely about CI/CD and build configuration changes. No local development tools need verification. CI runners provide their own environments.

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | CTest (wrapping GTest executables) |
| Config file | Inline in CMakeLists.txt via `enable_testing()` + `add_test()` |
| Quick run command | `cmake --build build && cd build && ctest --output-on-failure` |
| Full suite command | Same — all tests run in < 30s |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CI-01 | Library compiles with Eigen 5.x headers | smoke | `cmake --build build` | ✅ Phase execution |
| CI-01 | `find_package(Eigen3 3.4...5 REQUIRED)` accepts both versions | smoke | `cmake -B build -DEigen3_DIR=...` | ✅ CI step |
| CI-02 | Full test suite passes with Eigen 5.x | integration | `cd build && ctest --output-on-failure` | ✅ in CI workflow |
| CI-02 | 8 CI matrix jobs all green | e2e | GitHub Actions UI | ✅ 8 jobs in workflow |

### Sampling Rate
- **Per task commit:** N/A — build & test only in CI
- **Per wave merge:** Ensure all 8 CI jobs green
- **Phase gate:** CI-01 and CI-02 requirements satisfied (verify in GitHub Actions)

### Wave 0 Gaps
- None — existing test infrastructure (`tests/*.cpp`, CTest, CI workflow) covers all phase requirements. No new test files needed.

## Security Domain

> Not applicable — this phase touches only build configuration and CI/CD. No input validation, authentication, session management, access control, cryptography, or data handling. Skipping ASVS analysis.

## Sources

### Primary (HIGH confidence)
- [Eigen 5.0 CMake Guide](https://libeigen.gitlab.io/eigen/docs-5.0/TopicCMakeGuide.html) — Version range syntax `3.4...5` documented with example. Source: fetched via WebFetch.
- [Eigen 5.0 Release Notes](https://libeigen.gitlab.io/releases/5.0/) — Breaking changes list, C++14 requirement, versioning scheme. Source: fetched via WebFetch.
- [Eigen 5.0.1 tarball](https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz) — Installation via `cmake -B build -DCMAKE_INSTALL_PREFIX=...` pattern. Source: GitLab releases.
- [Eigen Macros.h](https://gitlab.com/libeigen/eigen/-/raw/master/Eigen/src/Core/util/Macros.h) — `EIGEN_MAJOR_VERSION`, `EIGEN_MINOR_VERSION`, `EIGEN_PATCH_VERSION` and `EIGEN_VERSION_AT_LEAST` macro definitions. Source: raw GitLab fetch.

### Secondary (MEDIUM confidence)
- [Eigen 5.0 issue #3004](https://gitlab.com/libeigen/eigen/-/issues/3004) — Documents 5.0.0 vs 5.0.1 CMake version incompatibility. Confirms the version range syntax is the correct fix. Source: WebSearch, GitLab issue.
- [Ceres Solver issue #1196](https://github.com/ceres-solver/ceres-solver/issues/1196) — Real-world example of `find_package(Eigen3 3.3)` failing with Eigen 5.0.1. Source: WebSearch, GitHub.
- Existing CI workflow and CMakeLists.txt — Verified by reading actual files. Source: codebase.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — CMake version range syntax directly from Eigen docs; CI matrix pattern is standard GitHub Actions technique.
- Architecture: HIGH — Matrix expansion, conditional install steps, coverage preservation all follow established patterns from existing CI workflow.
- Pitfalls: HIGH — Verified via Eigen GitLab issues, Ceres Solver issue, and release notes.
- API compatibility: MEDIUM — The primary risk (Eigen 5.x breaking this codebase's specific API usage) cannot be confirmed without compilation. However, the APIs used are the most stable core Eigen features.

**Research date:** 2026-06-05
**Valid until:** 2026-09-05 (Eigen 5.x release notes are stable; URLs may change)
