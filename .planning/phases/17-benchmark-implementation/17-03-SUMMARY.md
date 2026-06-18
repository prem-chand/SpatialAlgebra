# Plan 17-03 Summary: Solver benchmarks

## Status: Complete

### Changes
- **benchmarks/aba/bench_aba_stub.cpp**: Replaced stub with real ABA benchmark — ModelFactory::createFD, q/qdot setup, fillTorques inside loop, computeAccelerations, DoNotOptimize
- **benchmarks/rnea/bench_rnea_stub.cpp**: Replaced stub with real RNEA benchmark — ModelFactory::createID, q/qdot setup, fillAccelerations inside loop, computeTorques, DoNotOptimize

### Key design decisions
- Zero-alloc fill methods used inside timed loops (D-19)
- Joint state (q, qdot) set once outside timed loop
- RNEA's computeTorques allocates internally (inherent to current API — RESEARCH.md §Pitfall 1)

### Verification
- 0 `(void)state` or `(void)nDOF` casts in either file
- Key function calls present in both files
- `cmake --build build_bench` succeeds
- ABA 1DOF: ~84.5 µs timing output
- RNEA 1DOF: ~29.4 µs timing output
