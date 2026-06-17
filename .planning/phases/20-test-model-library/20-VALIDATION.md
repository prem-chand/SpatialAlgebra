---
phase: 20
slug: test-model-library
status: draft
nyquist_compliant: true
wave_0_complete: true
created: 2026-06-17
---

# Phase 20 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Test 1.12.1 (GTest, via CMake FetchContent) |
| **Config file** | none — GTest auto-detected by CMake `find_package(GTest)` |
| **Quick run command** | `cmake --build build && cd build && ./compile_smoke_test` |
| **Full suite command** | `cmake --build build && cd build && ctest --output-on-failure` |
| **Estimated runtime** | ~5 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cmake --build build && cd build && ctest -R smoke`
- **After every plan wave:** Run `cmake --build build && cd build && ctest --output-on-failure`
- **Before `/gsd-verify-work`:** All 11 existing test executables pass (no regression) + smoke test passes
- **Max feedback latency:** 10 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 20-01-01 | 01 | 1 | TML-01 | T-20-01 | Compilation firewall: no SA types in test_models headers | compile | `cmake --build build --target compile_smoke_test` | ❌ W0 | ⬜ pending |
| 20-01-02 | 01 | 1 | TML-05 | T-20-04 | INTERFACE target links only Eigen3 | build | `cmake --build build --target test_models` | ❌ W0 | ⬜ pending |
| 20-02-01 | 02 | 1 | TML-03 | — | RobotSolver pure virtual compiles and is mockable | compile | `cmake --build build --target compile_smoke_test` | ❌ W0 | ⬜ pending |
| 20-02-02 | 02 | 2 | TML-02 | — | 11 chains/ headers produce valid RobotModel | unit | `cd build && ./compile_smoke_test` | ❌ W0 | ⬜ pending |
| 20-03-01 | 03 | 2 | TML-04 | T-20-01, T-20-02, T-20-03 | Adapter wraps ID/FD; size checks throw std::invalid_argument; idx bounds throw std::out_of_range | integration | `cd build && ctest -R smoke` | ❌ W0 | ⬜ pending |
| 20-03-02 | 03 | 2 | TML-04 | T-20-05 | PIMPL pattern hides SA types; unique_ptr prevents use-after-move | unit | `cmake --build build --target sa_test_adapter` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test-models/CMakeLists.txt` — new file: INTERFACE target `test_models`, `sa_test_adapter`, `compile_smoke_test`
- [ ] `tests/test-models/robot_model.h` — new file: JointSpec, RobotModel, JointType enum
- [ ] `tests/test-models/robot_solver.h` — new file: RobotSolver abstract base class
- [ ] `tests/test-models/sa_adapter.h` — new file: SpatialAlgebraAdapter (PIMPL) declaration
- [ ] `tests/test-models/sa_adapter.cpp` — new file: SA adapter implementation
- [ ] `tests/test-models/chains/*.h` — 11 new files: model factory functions per test domain
- [ ] `tests/compile_smoke_test.cpp` — update: include `robot_model.h`, `robot_solver.h`, instantiate types
- [ ] Root `CMakeLists.txt` — add `add_subdirectory(tests/test-models)` after `enable_testing()`

*Major gaps: the entire `tests/test-models/` directory is new. All files are Wave 0.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| All 11 test domains covered (TML-02) | TML-02 | Manual audit of chains/ directory — 11 headers with 2-4 factory functions each | Count files in `tests/test-models/chains/`; verify each matches a test domain from CONTEXT.md D-10 |
| Zero SA dependency in headers (TML-01) | TML-01 | Verified by compile smoke test + manual grep | `grep -r "SpatialAlgebra" tests/test-models/robot_model.h tests/test-models/robot_solver.h` — must return zero matches |
| Existing tests not modified (TML-04) | TML-04 | Verified by git diff scope | `git diff --name-only` — no changes under `tests/Test*.cpp` except `compile_smoke_test.cpp` |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 10s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
