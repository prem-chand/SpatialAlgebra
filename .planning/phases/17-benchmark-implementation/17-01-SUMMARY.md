# Plan 17-01 Summary: Shared benchmark infrastructure

## Status: Complete

### Changes
- **benchmarks/common/random_state.h**: Added 6 new method declarations with Doxygen comments (randomTorques, randomAccelerations, fillPositions, fillVelocities, fillTorques, fillAccelerations)
- **benchmarks/common/random_state.cpp**: Added 6 implementations with correct distribution ranges (D-16: [-10,10] Nm, D-17: [-5,5] rad/s²)
- **benchmarks/bench_all.cpp**: Replaced 2 placeholder forward declarations/registrations with 18 individual core microbenchmark registrations (adjusted from 19 because cross(FV,MV) doesn't exist); removed `#include "model_factory.h"` and `#include "random_state.h"`; kept ABA/RNEA DOF sweep loop unchanged

### Verification
- `cmake --build build_bench --target bench_common` succeeds
- `cmake --build build_bench --target bench_all` succeeds
- Smoke test produces valid timing output
