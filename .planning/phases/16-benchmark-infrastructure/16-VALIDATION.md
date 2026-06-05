---
phase: 16
slug: benchmark-infrastructure
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-05
---

# Phase 16 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Build Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Benchmark v1.9.5 (via FetchContent) |
| **Config file** | `benchmarks/CMakeLists.txt` |
| **Quick build command** | `cmake --build build --target bench_all` |
| **Full build command** | `cmake --build build` |
| **Estimated build time** | ~30 seconds (incremental) / ~120 seconds (clean) |

---

## Sampling Rate

- **After every task commit:** Run `cmake --build build --target bench_all 2>&1 | tail -20`
- **After every plan wave:** Run full `cmake --build build`, then `build/bench_all --benchmark_list_tests`
- **Before `/gsd-verify-work`:** Full build must succeed, `--benchmark_list_tests` must enumerate all expected stubs
- **Max feedback latency:** 120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 16-01-01 | 01 | 1 | BINF-01 | T-16-01 / — | N/A | compile | `cmake --build build --target bench_all 2>&1` | ❌ W0 | ⬜ pending |
| 16-01-02 | 01 | 1 | BINF-02 | T-16-01 / — | N/A | compile | `test -d benchmarks/aba && test -d benchmarks/rnea && test -d benchmarks/core && test -d benchmarks/common` | ❌ W0 | ⬜ pending |
| 16-01-03 | 01 | 1 | BINF-03 | T-16-01 / — | N/A | compile | `cmake --build build --target bench_all 2>&1` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `benchmarks/CMakeLists.txt` — Google Benchmark FetchContent + build targets
- [ ] `benchmarks/common/model_factory.h` — shared model factory header
- [ ] `benchmarks/common/model_factory.cpp` — shared model factory implementation

*If none: "Existing infrastructure covers all phase requirements."*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Google Benchmark v1.9.5 version check | BINF-01 | Version pin verified at FetchContent declaration time | Check `benchmarks/CMakeLists.txt` for `GIT_TAG v1.9.5` |

*If none: "All phase behaviors have automated verification."*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
