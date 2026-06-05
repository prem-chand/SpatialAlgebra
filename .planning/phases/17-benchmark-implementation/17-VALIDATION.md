---
phase: 17
slug: benchmark-implementation
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-05
---

# Phase 17 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Benchmark v1.9.5 + Google Test (existing) |
| **Config file** | `benchmarks/CMakeLists.txt` — FetchContent, bench_all target |
| **Quick run command** | `cmake --build build_bench 2>&1 \| tail -3` |
| **Full suite command** | `cmake -B build_bench -DSA_BUILD_BENCHMARKS=ON && cmake --build build_bench && build_bench/benchmarks/bench_all --benchmark_min_time=0.1` |
| **Estimated runtime** | ~30 seconds (short min_time for fast feedback) |

---

## Sampling Rate

- **After every task commit:** `cmake --build build_bench 2>&1 | tail -5`
- **After every plan wave:** Full build + `bench_all --benchmark_min_time=0.1`
- **Before `/gsd:verify-work`:** Full build + full bench_all run
- **Max feedback latency:** ~30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 17-01-01 | 01 | 1 | BENCH-01 | — | N/A | build | `cmake --build build_bench` | ⬜ W0 | ⬜ pending |
| 17-01-02 | 01 | 1 | BENCH-01 | — | N/A | runtime | `build_bench/benchmarks/bench_all --benchmark_filter="BM_ABA" --benchmark_min_time=0.05` | ⬜ W0 | ⬜ pending |
| 17-02-01 | 02 | 2 | BENCH-02 | — | N/A | build | `cmake --build build_bench` | ⬜ W0 | ⬜ pending |
| 17-02-02 | 02 | 2 | BENCH-02 | — | N/A | runtime | `build_bench/benchmarks/bench_all --benchmark_filter="BM_RNEA" --benchmark_min_time=0.05` | ⬜ W0 | ⬜ pending |
| 17-03-01 | 03 | 1 | BENCH-03 | — | N/A | build | `cmake --build build_bench` | ⬜ W0 | ⬜ pending |
| 17-03-02 | 03 | 1 | BENCH-03 | — | N/A | runtime | `build_bench/benchmarks/bench_all --benchmark_filter="BM_Plucker\|BM_Cross\|BM_Inertia" --benchmark_min_time=0.05` | ⬜ W0 | ⬜ pending |

---

## Wave 0 Requirements

- `benchmarks/` subdirectory builds with `-DSA_BUILD_BENCHMARKS=ON`
- Google Benchmark v1.9.5 fetched via FetchContent
- Existing stub files compile and link

*Wave 0 already satisfied by Phase 16.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Timing results are physically plausible | BENCH-01, BENCH-02, BENCH-03 | Timing is stochastic — no deterministic pass/fail threshold | Run `bench_all --benchmark_format=csv` and verify: ABA/RNEA times increase with DOF, microbenchmarks are stable within 20% across runs |

*All phase behaviors have automated build/run verification. Timing plausibility is human-reviewed.*

---

## Validation Sign-Off

- [ ] All tasks have build or runtime automated verification
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
