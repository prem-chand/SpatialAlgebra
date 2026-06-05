# Phase 15: Eigen 5.x CI - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-05
**Phase:** 15-eigen-5-x-ci
**Areas discussed:** Version range syntax, CI matrix structure, Eigen 5.x API warnings, Coverage strategy

---

## Version Range Syntax

| Option | Description | Selected |
|--------|-------------|----------|
| Version range syntax | `find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)` — CMake 3.19+ version range | ✓ |
| Keep no-version-pin | Stay with current `find_package(Eigen3 REQUIRED NO_MODULE)` | |
| Dual find_package with fallback | Try 5.x first, fallback to 3.4 | |

**User's choice:** Version range syntax
**Notes:** Matches success criteria, declares intent clearly.

---

## CI Matrix Structure

| Option | Description | Selected |
|--------|-------------|----------|
| Add eigen-version to matrix | 8 jobs (2 OS × 2 compilers × 2 Eigen versions) | ✓ |
| OS-asymmetry hybrid | macOS → 5.x, Ubuntu → 3.4.x, 4 jobs | |
| Separate eigen5 job | 4-matrix + one explicit Eigen 5.x job | |

**User's choice:** Add eigen-version to matrix
**Notes:** Explicit testing of both versions on all platforms.

| Option | Description | Selected |
|--------|-------------|----------|
| brew install + install from source | macOS brew for 5.x, build 3.4 from source; Ubuntu apt for 3.4.x, build 5.x from source | ✓ |
| brew switch + apt versions | macOS brew switch, Ubuntu PPA 5.x | |
| git checkout of Eigen tag | Clone Eigen, checkout tags, set -DEigen3_DIR | |

**User's choice:** brew install + install from source
**Notes:** Most explicit version control at the cost of build-from-source time.

---

## Eigen 5.x API Warnings

| Option | Description | Selected |
|--------|-------------|----------|
| Fix forward-compatible code | Replace deprecated API with `#if EIGEN_VERSION` guards | ✓ |
| Suppress with -Wno-deprecated | Add flag to Eigen 5.x CI jobs | |
| Best effort | Fix easy ones, suppress the rest | |

**User's choice:** Fix forward-compatible code
**Notes:** Clean builds under both versions — no warning suppression.

---

## Coverage Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Single coverage config | Keep coverage on ubuntu+g++ with Eigen 3.4.x only | ✓ |
| Coverage on both | Run coverage on both Eigen versions, upload with separate flags | |
| Coverage on Eigen 5.x only | Move coverage from 3.4 to 5.x | |

**User's choice:** Single coverage config
**Notes:** Simpler reporting, avoid CodeCov merge issues.

---

## the agent's Discretion

- Specific Eigen 5.x API changes that need version-guarded fixes (discovered during compilation)
- Download/checksum details for Eigen 5.0.1 source tarball in CI
- Exact `apt` and `brew` version command details

## Deferred Ideas

None — discussion stayed within phase scope.
