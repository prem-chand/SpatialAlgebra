---
phase: 13
slug: production-readiness
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-05-17
---

# Phase 13 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Test (GTest) |
| **Config file** | `CMakeLists.txt` (test executables registered via `add_test()`) |
| **Quick run command** | `cmake --build build && cd build && ctest --output-on-failure -R Test` |
| **Full suite command** | `cmake --build build && cd build && ctest --output-on-failure` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run quick command (affected test executable)
- **After every plan wave:** Run full suite
- **Before `/gsd-verify-work`:** Full suite must be green (all tests passing)
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| TBD | TBD | 1 | VEC-01, UTL-03 | — | unit | `build/TestSpatialUtils` | ❌ W0 | ⬜ pending |
| TBD | TBD | 1 | INR-01 | — | unit | `build/TestArticulatedBodyInertia` | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | ABA-01, ABA-02 | — | unit | `build/TestForwardDynamics` | ✅ | ⬜ pending |
| TBD | TBD | 2 | TST-07 | — | integration | `build/TestDynamicsConsistency` | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements. GTest is already configured in CMakeLists.txt. Test executables already compiled via `cmake --build build`.

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
