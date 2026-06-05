# Phase 16: Benchmark Infrastructure - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-05
**Phase:** 16-benchmark-infrastructure
**Areas discussed:** Model factory design, Benchmark structure, Random state generation, Benchmark compilation profile

---

## Model Factory Design

| Option | Description | Selected |
|--------|-------------|----------|
| ForwardDynamics object | Factory returns ready-to-use ForwardDynamics | |
| Standalone link vectors | Factory returns std::vector\<ForwardDynamicsLink\> | |
| Unified API returning both | Single factory producing either solver type | ✓ |

**User's choice:** Unified API returning both
**Notes:** Factory should be shared across ABA and RNEA benchmarks.

| Option | Description | Selected |
|--------|-------------|----------|
| Revolute Z only | Z-axis revolute joints only | |
| Configurable axis per joint | XYZ axis per joint as parameter | ✓ |

**User's choice:** Configurable axis per joint

| Option | Description | Selected |
|--------|-------------|----------|
| Serial chains only | Linear chain topology | |
| Serial + branching | Support branching Y-shaped chains | ✓ |

**User's choice:** Serial + branching

| Option | Description | Selected |
|--------|-------------|----------|
| Uniform per-link | All links get same mass and inertia | |
| Random per-link | Links get random mass/COM | |
| Configurable parameter | Factory takes optional mass range, COM range | ✓ |

**User's choice:** Configurable parameter

---

## Benchmark Structure

| Option | Description | Selected |
|--------|-------------|----------|
| One executable per benchmark | Separate executables for each benchmark | |
| One unified executable | Single binary with sub-benchmark registration | ✓ |

**User's choice:** One unified executable

| Option | Description | Selected |
|--------|-------------|----------|
| Single CMakeLists.txt | One benchmarks/CMakeLists.txt for everything | |
| Subdirectories per domain | aba/, rnea/, core/, common/ subdirs | ✓ |

**User's choice:** Subdirectories per domain

| Option | Description | Selected |
|--------|-------------|----------|
| In benchmarks/CMakeLists.txt | FetchContent in subdirectory | ✓ |
| In top-level CMakeLists.txt | Central FetchContent | |

**User's choice:** In benchmarks/CMakeLists.txt

| Option | Description | Selected |
|--------|-------------|----------|
| Phase 16: common + skeleton | Shared utilities only now | |
| Phase 16: all source files | Create stubs for all executables now | ✓ |

**User's choice:** Phase 16: all source files

---

## Random State Generation

| Option | Description | Selected |
|--------|-------------|----------|
| Full rotation [0, 2π) | Uniform random in [0, 2π) | |
| Symmetrical [-π, π] | Uniform random in [-π, π] | |
| Constrained [-π/2, π/2] | Avoids kinematic singularities | ✓ |

**User's choice:** Constrained [-π/2, π/2]

| Option | Description | Selected |
|--------|-------------|----------|
| Small [-1, 1] rad/s | Low-velocity regime | |
| Moderate [-5, 5] rad/s | Realistic operating velocities | ✓ |
| Wide [-50, 50] rad/s | High-speed regime | |

**User's choice:** Moderate [-5, 5] rad/s

| Option | Description | Selected |
|--------|-------------|----------|
| Fixed seed (42) | Fully deterministic | ✓ |
| Time-based seed | Different chain each run | |

**User's choice:** Fixed seed (42)

---

## Benchmark Compilation Profile

| Option | Description | Selected |
|--------|-------------|----------|
| Release-only (-O3 -DNDEBUG) | Single profile | |
| Debug + Release | Both profiles | ✓ |

**User's choice:** Debug + Release

| Option | Description | Selected |
|--------|-------------|----------|
| Yes — Release by default | SA_BUILD_BENCHMARKS implies Release | |
| No — inherit from parent | Standard CMake behavior | ✓ |

**User's choice:** No — inherit from parent

| Option | Description | Selected |
|--------|-------------|----------|
| Yes | Enable -flto for benchmark targets | ✓ |
| No | Keep LTO off | |

**User's choice:** Yes (LTO enabled)

---

## the agent's Discretion

- Specific benchmark executable name (within `bench_all` convention)
- Exact subdirectory structure within each domain directory
- CMake minimum version for benchmarks/CMakeLists.txt
- Google Benchmark version pin details
- Shared utility function signatures and header file organization

## Deferred Ideas

None — discussion stayed within phase scope.
