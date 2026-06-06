# Validation Strategy — Phase 18

**Phase:** 18 — robot-examples
**Date:** 2026-06-06

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | CMake build + runtime execution |
| **Config file** | `examples/CMakeLists.txt` |
| **Quick build command** | `cmake --build build --target example_robot_2link --target example_robot_3link` |
| **Full build command** | `cmake --build build` |
| **Estimated build time** | ~15 seconds (incremental) |

## Sampling Rate

- **After every task commit:** `cmake --build build` compilation check
- **After every plan wave:** Full runtime test of both executables
- **Before `/gsd-verify-work`:** Both examples compile, run, and produce physically correct output
- **Max feedback latency:** 60 seconds

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 18-01-01 | 01 | 1 | EX-01 | T-18-01 / — | N/A | smoke | `cmake --build build && build/example_robot_2link` | ❌ W0 | ⬜ pending |
| 18-01-02 | 01 | 1 | EX-02 | T-18-01 / — | N/A | smoke | `cmake --build build && build/example_robot_3link` | ❌ W0 | ⬜ pending |
| 18-01-03 | 01 | 1 | EX-01, EX-02 | T-18-01 / — | N/A | build | `cmake -B build && cmake --build build` | example_robot_2link: ❌ W0, example_robot_3link: ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

## Wave 0 Requirements

- [ ] `examples/example_robot_2link.cpp` — 2-link Z-Z planar arm example
- [ ] `examples/example_robot_3link.cpp` — 3-link Z-Y-Z spatial arm example
- [ ] `examples/CMakeLists.txt` — edit to add two new targets

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Physical interpretation text | EX-01, EX-02 | Examples output to stdout with no assertions; correctness by inspection | Run both executables and verify output includes explanatory text with correct physics (Z-Z arm zero gravity torque, Z-Y-Z arm ~30 N·m on joint 2) |
| Cross-validation residual < 1e-10 | EX-01, EX-02 | Numerical tolerance checked during manual inspection | Verify printed qddot values from FD+ID round-trip are ≈ 0 (1e-12 to 1e-10) |

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
