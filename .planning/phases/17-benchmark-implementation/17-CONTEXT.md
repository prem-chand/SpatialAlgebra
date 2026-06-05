# Phase 17: Benchmark Implementation - Context

**Gathered:** 2026-06-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the stub benchmark functions created in Phase 16 with real timing loops for ABA forward dynamics (BENCH-01), RNEA inverse dynamics (BENCH-02), and core operations microbenchmarks (Plücker transforms, cross products — BENCH-03). All benchmark source files already exist as stubs in `benchmarks/{aba,rnea,core}/`; this phase fills the implementation bodies and extends `RandomState` with new generation methods.

</domain>

<decisions>
## Implementation Decisions

### ABA Benchmark Torque Strategy
- **D-16:** Random torques via `RandomState::randomTorques(nDOF)` with uniform range [-10, 10] Nm — user chose random over all-ones for realistic torque profiles, despite the researcher's recommendation for all-ones.

### RNEA Acceleration Strategy
- **D-17:** Random qddot via `RandomState::randomAccelerations(nDOF)` with uniform range [-5, 5] rad/s² — exercises all RNEA code paths including the S·q̈ term in the outward pass (InverseDynamics.cpp:47-50).

### Microbenchmark Scope
- **D-18:** Full API surface — all 20 operations benchmarked:
  - Plücker transforms (8): `transformMotion`, `transformForce`, `inverseTransformMotion`, `inverseTransformForce`, `tformRBI`, `invtformRBI`, `tformABI`, `invtformABI`
  - Cross products (4): `cross(mv,mv)`, `cross(mv,fv)`, `cross(fv,mv)`, `cross(fv,fv)`
  - Inertia operations (8): All `apply` variants on RigidBodyInertia and ArticulatedBodyInertia
  User chose full API surface over the researcher's recommended algorithm-critical subset for complete performance profiling.

### State Regeneration Strategy
- **D-19:** Per-iteration random state generation using pre-allocated `Eigen::VectorXd` buffers — add non-allocating `fillPositions(double*, int)` and `fillVelocities(double*, int)` methods to `RandomState` (and corresponding `fillTorques`, `fillAccelerations`) to avoid heap allocation inside the timed loop.

### RandomState Extensions (new methods)
- `randomTorques(int nDOF)` returning `Eigen::VectorXd` — uniform [-10, 10] Nm (D-16)
- `randomAccelerations(int nDOF)` returning `Eigen::VectorXd` — uniform [-5, 5] rad/s² (D-17)
- `fillPositions(double* buf, int n)` / `fillVelocities(double* buf, int n)` — zero-alloc variants for per-iteration use (D-19)
- `fillTorques(double* buf, int n)` / `fillAccelerations(double* buf, int n)` — zero-alloc variants (D-19)

### Inherited Decisions (from Phase 16)
- **D-05:** Single `bench_all` executable with `RegisterBenchmark` — sub-benchmarks filterable via `--benchmark_filter`
- **D-06:** Subdirectory layout: `benchmarks/{aba,rnea,core,common}/`
- **D-08:** All source files already exist as stubs
- **D-09:** Joint positions: uniform random in [-π/2, π/2]
- **D-10:** Joint velocities: uniform random in [-5, 5] rad/s
- **D-11:** Fixed random seed 42

### Claude's Discretion
- Exact benchmark function signatures (within `void BM_*(benchmark::State&, int nDOF)` pattern)
- DOF iteration step size (1 or configurable) within the 1..20 range
- Google Benchmark `MinTime` / iterations configuration
- Specific `DoNotOptimize` / `ClobberMemory` placement in microbenchmark loops

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Build System & Infrastructure
- `CMakeLists.txt` — Root build config (FetchContent pattern, option pattern)
- `benchmarks/CMakeLists.txt` — Google Benchmark FetchContent, bench_all target
- `.planning/REQUIREMENTS.md` §BENCH-01, BENCH-02, BENCH-03 — Requirement definitions

### Existing Code (Stubs to Fill)
- `benchmarks/aba/bench_aba_stub.cpp` — ABA benchmark stub (TODO for Phase 17)
- `benchmarks/rnea/bench_rnea_stub.cpp` — RNEA benchmark stub (TODO for Phase 17)
- `benchmarks/core/bench_core_stub.cpp` — Core microbenchmark stub (TODO for Phase 17)
- `benchmarks/bench_all.cpp` — Main benchmark registration (DOF sweep loop)
- `benchmarks/common/model_factory.h` / `model_factory.cpp` — Factory API for solver construction
- `benchmarks/common/random_state.h` / `random_state.cpp` — Random state generator (needs extension)

### Solver APIs (Benchmarked Code)
- `include/ForwardDynamics.h` — `ForwardDynamics::computeAccelerations(tau)` interface
- `include/InverseDynamics.h` — `InverseDynamics::computeTorques(q, qdot, qddot)` interface
- `include/PluckerTransform.h` — All transform methods (motion, force, inertia variants)
- `include/RigidBodyInertia.h` — Inertia operations
- `include/ArticulatedBodyInertia.h` — Articulated inertia operations
- `include/SpatialUtils.h` — `cross()` free functions

### Prior Phase Context
- `.planning/phases/16-benchmark-infrastructure/16-CONTEXT.md` — Locked decisions D-01 through D-15
- `.planning/phases/14-cr-02-bug-fix/14-CONTEXT.md` — ABA condensation fix context

### Project Context
- `.planning/PROJECT.md` §Out of Scope — Confirms Python benchmark harness out of scope

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Benchmark stubs** (`benchmarks/{aba,rnea,core}/*_stub.cpp`) — Already include correct headers and function signatures, just need body implementation
- **RandomState** (`benchmarks/common/random_state.h`) — Existing class with `randomPositions()`, `randomVelocities()` returning `std::vector<double>`. Needs new methods for torques, accelerations, and zero-alloc fill variants.
- **ModelFactory** (`benchmarks/common/model_factory.h`) — `createFD(nDOF)` returns populated `ForwardDynamics` solver, `createID(nDOF)` returns `InverseDynamics` solver
- **Google Benchmark** — `DoNotOptimize()`, `ClobberMemory()`, `State` iterator for timing loops

### Established Patterns
- **Programmatic registration**: `benchmark::RegisterBenchmark(name, fn, args...)` in main() for DOF sweep n=1..20
- **Per-iteration state pattern**: `for (auto _ : state) { operation(); DoNotOptimize(result); }`
- **Setup outside timing loop**: ModelFactory + buffer pre-allocation before `for (auto _ : state)`

### Integration Points
- `benchmarks/common/random_state.h` — Add `randomTorques()`, `randomAccelerations()`, `fill*()` methods
- `benchmarks/aba/bench_aba_stub.cpp` — Replace `(void)state; (void)nDOF;` stub with real implementation
- `benchmarks/rnea/bench_rnea_stub.cpp` — Same for RNEA
- `benchmarks/core/bench_core_stub.cpp` — Same for microbenchmarks

</code_context>

<specifics>
## Specific Ideas

- ABA benchmark flow: ModelFactory::createFD(nDOF) → RandomState → applyState(fd, q, qdot) → computeAccelerations(randomTorques) → DoNotOptimize(result)
- RNEA benchmark flow: ModelFactory::createID(nDOF) → RandomState → applyState(id, q, qdot, qddot) → computeTorques() → DoNotOptimize(result)
- Microbenchmarks: Instantiate single PlückerTransform / MotionVector / ForceVector, loop operation with DoNotOptimize output
- Random qddot range [-5, 5] rad/s² matches D-10 velocity range (consistency)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 17-benchmark-implementation*
*Context gathered: 2026-06-05*
