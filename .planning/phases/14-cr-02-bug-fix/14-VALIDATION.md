---
phase: 14
slug: cr-02-bug-fix
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-04
---

# Phase 14 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Test (GTest) |
| **Config file** | CMakeLists.txt — `enable_testing()` + `add_test()` + `gtest_discover_tests()` |
| **Quick run command** | `cmake --build build && cd build && ctest --output-on-failure` |
| **Full suite command** | `cmake --build build && cd build && ctest --output-on-failure` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cmake --build build && cd build && ctest --output-on-failure`
- **After every plan wave:** Run `cmake --build build && cd build && ctest --output-on-failure`
- **Before `/gsd-verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 14-01 | 01 | 1 | BFIX-01 | N/A | N/A — pure computation, no external input | round-trip | `build/TestDynamicsConsistency --gtest_filter=*ThreeLink*` | ✅ | ⬜ pending |
| 14-01 | 01 | 1 | BFIX-01 | N/A | N/A | unit | `build/TestForwardDynamics --gtest_filter=*Condensation*` | ❌ W0 | ⬜ pending |
| 14-01 | 01 | 1 | BFIX-01 | N/A | N/A | unit | `build/TestForwardDynamics --gtest_filter=*TauOrdering*` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/TestForwardDynamics.cpp` — condensation unit test: 3-link chain, condensed Ia norm < uncondensed Ia norm
- [ ] `tests/TestForwardDynamics.cpp` — tau=[1,0,0] on 3-link chain, verify qddot[0] < qddot[2]

*If none: "Existing infrastructure covers all phase requirements."*

---

## Manual-Only Verifications

All phase behaviors have automated verification.

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
