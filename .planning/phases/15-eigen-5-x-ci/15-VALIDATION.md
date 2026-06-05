# Validation Strategy — Phase 15

**Phase:** 15 — eigen-5-x-ci
**Date:** 2026-06-05

## Test Framework

| Property | Value |
|----------|-------|
| Framework | CTest (wrapping GTest executables) |
| Config file | Inline in CMakeLists.txt via `enable_testing()` + `add_test()` |
| Quick run command | `cmake --build build && cd build && ctest --output-on-failure` |
| Full suite command | Same — all tests run in < 30s |

## Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | Verification Method |
|--------|----------|-----------|-------------------|---------------------|
| CI-01 | Library compiles with Eigen 5.x headers | smoke | `cmake --build build` | CI step succeeds |
| CI-01 | `find_package(Eigen3 3.4...5 REQUIRED)` accepts both versions | smoke | `cmake -B build -DEigen3_DIR=...` | CI configure step succeeds |
| CI-02 | Full test suite passes with Eigen 5.x | integration | `cd build && ctest --output-on-failure` | CTest exit code 0 |
| CI-02 | 8 CI matrix jobs all green | e2e | GitHub Actions UI | All 8 jobs pass |

## Sampling Rate

- **Per wave merge:** Ensure all 8 CI jobs green
- **Phase gate:** CI-01 and CI-02 requirements satisfied (verify in GitHub Actions)

## Wave 0 Gaps

None — existing test infrastructure (`tests/*.cpp`, CTest, CI workflow) covers all phase requirements. No new test files needed.

---

*Derived from RESEARCH.md Validation Architecture section*
