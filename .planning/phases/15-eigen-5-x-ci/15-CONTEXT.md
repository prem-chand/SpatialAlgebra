# Phase 15: Eigen 5.x CI - Context

**Gathered:** 2026-06-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Add Eigen 5.x to the CI build matrix while maintaining backward compatibility with Eigen 3.4.x. This is an infrastructure/CI phase: update CMakeLists.txt version range, expand CI matrix to include explicit Eigen 5.x jobs, fix any forward-compatible API deprecations with version guards, and verify all tests pass under both Eigen versions.

Requirements CI-01 (version range syntax) and CI-02 (Eigen 5.x in CI matrix) define WHAT — this context captures HOW.

</domain>

<decisions>
## Implementation Decisions

### CMake Version Range
- **D-01:** Use `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` version range syntax. Requires CMake 3.19+ (already satisfied by CI runners). Declares intent clearly while allowing 3.4.x through 5.x.

### CI Matrix Structure
- **D-02:** Add `eigen-version` dimension to the CI matrix, producing 8 jobs (2 OS × 2 compilers × 2 Eigen versions: `3.4` and `5.0`). Explicit testing of both versions on all platforms.
- **D-03:** Install Eigen versions explicitly per job:
  - macOS: `brew install eigen@3.4` (or from source if tap unavailable) and `brew install eigen` (latest = 5.x)
  - Ubuntu: `apt install libeigen3-dev` (3.4.x) and build Eigen 5.0.1 from source tarball
  - Use `-DEigen3_DIR` to point CMake at the correct installation

### Eigen 5.x API Warnings
- **D-04:** Fix forward-compatible code using `#if EIGEN_VERSION` preprocessor guards where 3.4.x and 5.x need different code paths. Clean builds under both versions — no warning suppression.

### Coverage
- **D-05:** Keep single coverage run on `ubuntu-latest` + `g++` with Eigen 3.4 (the existing config). Coverage reporting does not extend to Eigen 5.x jobs. Add CodeCov `flags` to distinguish coverage source if desired.

### the agent's Discretion
- Specific Eigen 5.x API changes that need version-guarded fixes (to be discovered during compilation)
- Download/checksum details for Eigen 5.0.1 source tarball in CI
- Exact `apt` and `brew` version command details

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Build & CI
- `CMakeLists.txt:12` — Current `find_package(Eigen3 REQUIRED NO_MODULE)` (target for D-01)
- `.github/workflows/ci.yml` — Current 4-matrix CI workflow (target for D-02/D-03)
- `.planning/REQUIREMENTS.md` §CI-01, CI-02 — Requirement definitions

### Prior Context
- `.planning/phases/13-production-readiness/13-CONTEXT.md` §D-18 — Version pin removal decision
- `.planning/phases/13-production-readiness/13-CONTEXT.md` §D-09–D-12 — CI pipeline pattern established

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Existing CI workflow** (`.github/workflows/ci.yml`) — 4-matrix pattern ready to extend with `eigen-version` dimension
- **CMakeLists.txt** — No version pin currently, clean target for version range syntax
- **FetchContent pattern** — Already used for GTest fallback; can be adapted for Eigen 5.x from-source build in CI

### Established Patterns
- **Matrix strategy** — Existing `os` + `compiler` matrix in CI; adding `eigen-version` follows same pattern
- **Conditional installation** — `if: runner.os == 'macOS'` / `if: runner.os == 'Linux'` blocks for OS-specific package install

### Integration Points
- `CMakeLists.txt:12` — Single line change: version range syntax
- `.github/workflows/ci.yml` — Add `eigen-version` to matrix, add install/build steps per version
- No source code changes expected beyond potential `#if EIGEN_VERSION` guards

</code_context>

<specifics>
## Specific Ideas

No specific external examples cited — open to standard approaches for CI matrix expansion and Eigen version management.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 15-eigen-5-x-ci*
*Context gathered: 2026-06-05*
