# Phase 16: Benchmark Infrastructure - Context

**Gathered:** 2026-06-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Build system integration for Google Benchmark v1.9.5 (via FetchContent), benchmarks/ directory structure, and shared benchmark utilities (model factory, random joint state generator). This is infrastructure for Phases 17-18 — creates the skeleton that benchmark authors will fill.

Requirements BINF-01 (FetchContent), BINF-02 (benchmarks/ directory with SA_BUILD_BENCHMARKS guard), and BINF-03 (shared utilities) define WHAT — this context captures HOW.

</domain>

<decisions>
## Implementation Decisions

### Model Factory Design
- **D-01:** Unified API returning both ForwardDynamics and InverseDynamics solver objects — caller specifies which solver type to construct
- **D-02:** Configurable joint axis per joint (Z, X, Y, or arbitrary screw axis per link)
- **D-03:** Support both serial chains and branching (Y-shaped) configurations
- **D-04:** Link inertia configurable via parameters (mass range, COM range, inertia tensor range) — allows uniform, random, or explicit per-link specification

### Benchmark Structure
- **D-05:** Single unified executable (`bench_all`) with Google Benchmark sub-benchmark registration — benchmarks filterable at runtime via `--benchmark_filter`
- **D-06:** Subdirectories per domain: `benchmarks/aba/`, `benchmarks/rnea/`, `benchmarks/core/` (Plücker, cross-product), `benchmarks/common/` (shared utilities)
- **D-07:** Google Benchmark FetchContent lives in `benchmarks/CMakeLists.txt` — main CMakeLists.txt only adds `add_subdirectory(benchmarks)` guarded by `SA_BUILD_BENCHMARKS`
- **D-08:** Phase 16 creates all benchmark source files (stubs for Phases 17-18 to fill with benchmark logic)

### Random State Generation
- **D-09:** Joint positions: uniform random in [-π/2, π/2]
- **D-10:** Joint velocities: uniform random in [-5, 5] rad/s
- **D-11:** Fixed random seed 42 for reproducibility across runs

### Build Configuration
- **D-12:** Both Debug and Release build profiles supported
- **D-13:** Default build type inherits from parent CMake project (no forced Release)
- **D-14:** LTO (`-flto`) enabled for benchmark targets
- **D-15:** Guard option `SA_BUILD_BENCHMARKS` (default OFF) in top-level CMakeLists.txt

### the agent's Discretion
- Specific benchmark executable name (within `bench_all` convention)
- Exact subdirectory structure within each domain directory
- CMake minimum version for `benchmarks/CMakeLists.txt` (inherit from parent)
- Google Benchmark version pin (v1.9.5 per requirements) and FetchContent URL details
- Shared utility function signatures and header file organization
</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Build System
- `CMakeLists.txt:17-23` — Existing FetchContent pattern for GTest (reuse for Google Benchmark)
- `examples/CMakeLists.txt` — Existing subdirectory CMakeLists.txt pattern (reuse for benchmarks/)
- `.planning/REQUIREMENTS.md` §BINF-01, BINF-02, BINF-03 — Requirement definitions

### Existing Code
- `include/ForwardDynamics.h:79-111` — ForwardDynamicsLink struct (determines factory output API)
- `include/InverseDynamics.h` — InverseDynamicsLink struct (for factory output API)
- `include/ForwardDynamics.h:134-156` — ForwardDynamics class interface
- `include/InverseDynamics.h` — InverseDynamics class interface

### Prior Context
- `.planning/phases/14-cr-02-bug-fix/14-CONTEXT.md` — ABA inward pass condensation (factory must produce correct multi-link chains)
- `.planning/phases/13-production-readiness/13-CONTEXT.md` §D-20 — CI pipeline established (benchmarks need SA_BUILD_BENCHMARKS guard)

### Project Context
- `.planning/PROJECT.md` §Out of Scope — Confirms Python benchmark harness out of scope, Google Benchmark is the standard
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **FetchContent pattern** (`CMakeLists.txt:17-23`) — Established pattern for external dependencies (used for GTest). Reuse for Google Benchmark in `benchmarks/CMakeLists.txt`.
- **Subdirectory CMakeLists.txt** (`examples/CMakeLists.txt`) — Established pattern for building sub-projects with correct include/link paths.
- **ForwardDynamicsLink struct** (`include/ForwardDynamics.h:79-111`) — Contains all fields needed for factory output (parent, X, I, S, q, qdot, etc.)
- **Link construction** — Existing test setup code in tests/ is a pattern reference for factory construction

### Established Patterns
- **add_subdirectory** for modular builds (examples/CMakeLists.txt included from root)
- **option()** for build guards (ENABLE_COVERAGE in CMakeLists.txt:42)
- **FetchContent + QUIET** for optional external deps (GTest pattern at line 15)
- **Link ordering** — parents before children in the links vector (established in ForwardDynamics.h)

### Integration Points
- `CMakeLists.txt` — Add `option(SA_BUILD_BENCHMARKS OFF)` + `add_subdirectory(benchmarks)` guard
- `benchmarks/CMakeLists.txt` — New file: FetchContent for Google Benchmark, build shared utils + all benchmark stubs
- `benchmarks/common/` — New directory: model_factory.h and model_factory.cpp with unified API (D-01)
- `benchmarks/aba/`, `benchmarks/rnea/`, `benchmarks/core/` — New directories: skeleton benchmark files
</code_context>

<specifics>
## Specific Ideas

No specific external examples cited — open to standard Google Benchmark practices.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 16-benchmark-infrastructure*
*Context gathered: 2026-06-05*
