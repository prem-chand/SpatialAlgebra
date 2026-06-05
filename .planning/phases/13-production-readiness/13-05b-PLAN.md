---
phase: 13-production-readiness
plan: 05b
type: execute
wave: 3
depends_on:
  - 13-01
  - 13-04
files_modified:
  - CMakeLists.txt
  - .github/workflows/ci.yml
  - .gitignore
  - tests/compile_smoke_test.cpp
autonomous: true
requirements:
  - VEC-01
must_haves:
  truths:
    - "Library builds without system-installed GTest via FetchContent fallback (GTest QUIET + FetchContent)"
    - "CI runs all tests on push/PR to main branch with 4-matrix build (ubuntu/macos × g++/clang++)"
    - "Umbrella header compile smoke test exists (tiny target includes only SpatialAlgebra.h)"
    - "CI matrix has explicit Eigen 5.x compatibility job"
    - "Coverage build uses CMake ENABLE_COVERAGE option (default OFF, single CI job)"
    - ".gitignore does not exclude .github/ directory from version control"
  artifacts:
    - path: "CMakeLists.txt"
      provides: "GTest FetchContent fallback, ENABLE_COVERAGE option, stub exclusion, umbrella header smoke test target"
      contains: "FetchContent_MakeAvailable"
    - path: ".github/workflows/ci.yml"
      provides: "CI pipeline with 4-matrix build + Eigen 5.x compat job + CodeCov upload"
      contains: "eigen5-compat"
    - path: "tests/compile_smoke_test.cpp"
      provides: "Compile-only test verifying SpatialAlgebra.h is self-contained"
      contains: "#include \"SpatialAlgebra.h\""
    - path: ".gitignore"
      provides: "CI workflow file un-ignored"
      does_not_contain: "^\.github"
  key_links:
    - from: "CMakeLists.txt"
      to: "tests/compile_smoke_test.cpp"
      via: "add_executable(SpatialAlgebraCompileSmoke ...) — verifies SpatialAlgebra.h compiles alone"
      pattern: "SpatialAlgebraCompileSmoke"
    - from: ".github/workflows/ci.yml"
      to: "CMakeLists.txt"
      via: "ci workflow runs cmake configure with various compilers"
      pattern: "cmake -B build"
user_setup:
  - service: github-actions
    why: "CI pipeline runs on GitHub"
    env_vars:
      - name: GITHUB_TOKEN
        source: "GitHub automatic token (provided by Actions runtime)"
---
<objective>
Update CMake (GTest FetchContent fallback, stub removal, coverage option, compile smoke test), add CI workflow with explicit Eigen 5.x compatibility job, create compile smoke test, and update .gitignore.

**Purpose:** Make the library buildable without system-installed GTest and automatically verified on every push/PR.
**Output:** Buildable library with CI, compile smoke test, clean CMake, `.gitignore` updated for CI workflow.
</objective>

<execution_context>
@/Users/premchand/.config/opencode/get-shit-done/workflows/execute-plan.md
@/Users/premchand/.config/opencode/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/13-production-readiness/13-CONTEXT.md (D-25, D-27, D-11, D-28 smoke test, D-09 through D-12)
@.planning/phases/13-production-readiness/13-RESEARCH.md (CI patterns lines 85-120)
@.planning/phases/13-production-readiness/13-PATTERNS.md (All patterns for CI, CMake)
@.planning/phases/13-production-readiness/13-REVIEWS.md (13-05 concerns — FetchContent optional, Eigen 5 CI job, umbrella smoke test)

<interfaces>
<!-- Current CMakeLists.txt:12-18 -->
Line 12: `find_package(Eigen3 REQUIRED NO_MODULE)` — no version pin (D-18 already done)
Line 15: `find_package(GTest REQUIRED)` — change to QUIET + FetchContent fallback (D-25)
Line 18: `file(GLOB SOURCES "src/*.cpp")` — add stub exclusion (D-27)

<!-- Umbrella smoke test pattern -->
NEW FILE: tests/compile_smoke_test.cpp:
```cpp
// Compile smoke test: verifies SpatialAlgebra.h is self-contained
// and compiles standalone (no prior includes needed)
#include "SpatialAlgebra.h"
int main() { return 0; }
```

<!-- Current .gitignore -->
May contain `.github/` which blocks CI workflow from being tracked. Remove that line.
</interfaces>
</context>

<tasks>

<task type="auto">
  <name>Task 1 (Group B — Build Infrastructure): CMake updates — GTest fallback, stub removal, coverage, umbrella smoke test, .gitignore</name>
  <files>CMakeLists.txt, .gitignore</files>
  <read_first>
    CMakeLists.txt:1-35      (full CMakeLists.txt — project setup, find_package, sources, targets, testing)
    .gitignore:1-10           (current gitignore — check if .github/ is ignored)
  </read_first>
  <action>
    Per D-25, D-27, D-11, D-28 (smoke test), and review concern about FetchContent Homebrew conflicts and `.gitignore` blocking CI workflow.

    **CMakeLists.txt changes:**

    1. **GTest FetchContent fallback with Homebrew-friendly gating (D-25):**
       Replace line 15: `find_package(GTest REQUIRED)`
       with:
       ```cmake
       # Find GTest optionally; fall back to FetchContent if not installed
       find_package(GTest QUIET)
       if(NOT GTest_FOUND)
           include(FetchContent)
           FetchContent_Declare(
               googletest
               URL https://github.com/google/googletest/archive/release-1.12.1.zip
           )
           FetchContent_MakeAvailable(googletest)
           include_directories(${googletest_SOURCE_DIR}/googletest/include)
           message(STATUS "GTest not found — using FetchContent fallback")
       else()
           message(STATUS "GTest found at ${GTEST_LIBRARIES}")
       endif()
       ```
       The `QUIET` mode means system-installed GTest (e.g., `brew install googletest`) takes priority. Only when absent does FetchContent download.

    2. **Empty stub removal from source glob (D-27):**
       After the `file(GLOB SOURCES "src/*.cpp")` line, add:
       ```cmake
       # Remove empty stub files from compilation
       list(REMOVE_ITEM SOURCES
           "${CMAKE_SOURCE_DIR}/src/RigidBodyInertia.cpp"
           "${CMAKE_SOURCE_DIR}/src/ArticulatedBodyInertia.cpp"
       )
       ```
       (Stub files will be deleted by the agent during execution — the CMake exclusion ensures the build succeeds whether or not the physical files exist.)

    3. **Coverage option (D-11) — single CI job only:**
       After `enable_testing()`, add:
       ```cmake
       # Code coverage build (CI only — default OFF for local builds)
       option(ENABLE_COVERAGE "Enable coverage flags for gcov" OFF)
       if(ENABLE_COVERAGE)
           set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
           set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
       endif()
       ```
       The option defaults to OFF. Only the ubuntu+g++ CI job passes `-DENABLE_COVERAGE=ON`.

    4. **Umbrella header compile smoke test (D-28):**
       After the existing test registration (after `add_test` lines), add:
       ```cmake
       # Compile smoke test: verify SpatialAlgebra.h is self-contained
       add_executable(SpatialAlgebraCompileSmoke
           "${CMAKE_SOURCE_DIR}/tests/compile_smoke_test.cpp"
       )
       target_include_directories(SpatialAlgebraCompileSmoke PRIVATE
           "${CMAKE_SOURCE_DIR}/include"
       )
       target_link_libraries(SpatialAlgebraCompileSmoke PRIVATE SpatialAlgebra)
       add_test(NAME CompileSmokeTest COMMAND SpatialAlgebraCompileSmoke)
       ```
       The test source file (`tests/compile_smoke_test.cpp`) will be created during execution with:
       ```cpp
       // Compile smoke test: verifies SpatialAlgebra.h is self-contained
       #include "SpatialAlgebra.h"
       int main() { return 0; }
       ```

    5. **.gitignore:**
       Read `.gitignore` and remove any line that ignores `.github/`. Specifically, if `.github` or `.github/` appears in `.gitignore`, remove it. This ensures the CI workflow file is tracked by git.

    **Git commit boundaries:**
    This single commit covers all CMake and .gitignore changes.
  </action>
  <verify>
    <automated>cmake -B build 2>&1 | tail -10 && cmake --build build 2>&1 | tail -10 && cd build && ctest --output-on-failure 2>&1 | tail -20</automated>
    <human-check>Confirm compile smoke test executable exists and runs: `ls build/compile* 2>/dev/null || ls build/SpatialAlgebraCompile* 2>/dev/null`</human-check>
  </verify>
  <acceptance_criteria>
    - CMakeLists.txt: `find_package(GTest QUIET)` + FetchContent fallback block (system GTest preferred)
    - CMakeLists.txt: `list(REMOVE_ITEM SOURCES ...)` for both stub files
    - CMakeLists.txt: `option(ENABLE_COVERAGE ... OFF)` block — defaults OFF
    - CMakeLists.txt: `add_executable(SpatialAlgebraCompileSmoke ...)` block
    - `tests/compile_smoke_test.cpp` exists and contains `#include "SpatialAlgebra.h"`
    - `cmake -B build && cmake --build build` succeeds on a system with and without system-installed GTest
    - `cd build && ctest -R CompileSmokeTest` passes
    - `.gitignore` does NOT contain `.github/` or `.github` exclusion
    - **Verification command:** `grep -c "github" .gitignore` returns 0
  </acceptance_criteria>
</task>

<task type="auto">
  <name>Task 2 (Group B — CI Infrastructure): Create GitHub Actions CI workflow with Eigen 5.x compat job and Ubuntu+clang++ coverage job</name>
  <files>.github/workflows/ci.yml</files>
  <read_first>
    (No existing CI file — first CI workflow for this repository)
  </read_first>
  <action>
    Per D-09 through D-12, review concern about Eigen 5.x compatibility (single job), and concern about coverage in a single job only.

    Create NEW file `.github/workflows/ci.yml`:

    ```yaml
    name: CI

    on:
      push:
        branches: [ main ]
      pull_request:
        branches: [ main ]

    jobs:
      # Main build matrix: 4 configurations covering OS × compiler
      build:
        strategy:
          matrix:
            os: [ubuntu-latest, macos-latest]
            compiler: [g++, clang++]

        runs-on: ${{ matrix.os }}

        steps:
        - uses: actions/checkout@v4

        - name: Install dependencies (macOS)
          if: runner.os == 'macOS'
          run: |
            brew update
            brew install eigen googletest

        - name: Install dependencies (Ubuntu)
          if: runner.os == 'Linux'
          run: |
            sudo apt-get update
            sudo apt-get install -y libeigen3-dev libgtest-dev lcov

        - name: Configure
          run: |
            cmake -B build \
              -DCMAKE_CXX_COMPILER=${{ matrix.compiler }} \
              ${{ (matrix.os == 'ubuntu-latest' && matrix.compiler == 'clang++') && '-DENABLE_COVERAGE=ON' || '' }}

        - name: Build
          run: cmake --build build

        - name: Test
          run: cd build && ctest --output-on-failure

        - name: Upload coverage (ubuntu+clang++ only)
          if: matrix.os == 'ubuntu-latest' && matrix.compiler == 'clang++'
          uses: codecov/codecov-action@v4
          with:
            directory: ./build
            fail_ci_if_error: false

      # Eigen 5.x compatibility: single explicit job using Homebrew-installed Eigen 5
      eigen5-compat:
        runs-on: macos-latest
        steps:
        - uses: actions/checkout@v4

        - name: Install Eigen 5.x (Homebrew edge)
          run: |
            brew update
            brew install eigen

        - name: Configure with Eigen 5.x
          run: cmake -B build-eigen5 \
            -DEigen3_DIR="$(brew --prefix eigen)/share/eigen3/cmake"

        - name: Build
          run: cmake --build build-eigen5

        - name: Test
          run: cd build-eigen5 && ctest --output-on-failure
    ```

    Notes:
    - **Coverage:** Single job (ubuntu+clang++ only)
    - **Eigen 5.x job:** Separate job on macOS with explicit `-DEigen3_DIR=` pointing to Homebrew's Eigen 5 cmake config.
    - **4-matrix build:** ubuntu+g++, ubuntu+clang++, macos+g++, macos+clang++.
    - **FetchContent not needed in CI:** Both Ubuntu and macOS install system GTest via apt/brew.
    - **CodeCov fails silently:** `fail_ci_if_error: false` ensures coverage upload issues don't fail the CI.
  </action>
  <verify>
    <automated>python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml')); print('Valid YAML')" 2>&1</automated>
    <human-check>Verify YAML structure is correct by inspection.</human-check>
  </verify>
  <acceptance_criteria>
    - `.github/workflows/ci.yml` exists and is valid YAML
    - CI triggers on push and PR to main branch
    - **build job:** 4-matrix (ubuntu/macos × g++/clang++) with cmake configure/build/test
    - **Coverage upload step:** single job only (ubuntu+clang++), not all matrix jobs
    - **eigen5-compat job:** separate job on macOS with `-DEigen3_DIR` pointing to Homebrew Eigen 5
    - Dependencies: eigen3 + googletest installed via brew (macOS) or apt (Ubuntu)
    - No coverage flags applied to non-ubuntu+clang++ jobs
    - **Verification command:** `python3 -c "import yaml; d=yaml.safe_load(open('.github/workflows/ci.yml')); assert 'eigen5-compat' in d.get('jobs',{}), 'Missing eigen5-compat job'"`
  </acceptance_criteria>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| CI workflow → GitHub Actions | YAML config executes on GitHub runners; external dependency downloads from brew/apt |
| FetchContent → internet | CMake downloads googletest release zip from github.com — pinned to release-1.12.1 tag |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-13-06 | Tampering | CMake GTest FetchContent URL | mitigate | URL pinned to GitHub release-1.12.1.zip (not mutable branch reference) |
| T-13-07 | Tampering | CI workflow actions/checkout | accept | `actions/checkout@v4` pinned to major version v4 |
| T-13-09 | Tampering | Compile smoke test | mitigate | New test target verifies umbrella header is self-contained — catches missing includes |
| T-13-SC | Tampering | brew/apt package downloads | accept | Packages from official OS package managers; Eigen 5 from Homebrew |
</threat_model>

<verification>
- `cmake -B build && cmake --build build && cd build && ctest --output-on-failure` — full suite green (excluding pre-existing CR-02 failures)
- `tests/compile_smoke_test.cpp` exists and standalone compile test passes
- `.github/workflows/ci.yml` is valid YAML with eigen5-compat job
- `.gitignore` does not exclude `.github/`
</verification>

<success_criteria>
1. **CMake:** GTest QUIET+FetchContent fallback; stub exclusion; coverage option default OFF; compile smoke test target
2. **CI workflow:** 4-matrix build + separate eigen5-compat job + coverage on single ubuntu+clang++ job
3. `.gitignore` updated to track `.github/workflows/ci.yml`
4. Full build succeeds, compile smoke test passes, all existing tests pass
</success_criteria>

<output>
Create `.planning/phases/13-production-readiness/13-05b-SUMMARY.md` when done
</output>
