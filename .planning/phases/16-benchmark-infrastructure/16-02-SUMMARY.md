---
phase: 16-benchmark-infrastructure
plan: 02
subsystem: testing
tags: benchmark, model-factory, random-state, mt19937, deterministic-rng

# Dependency graph
requires:
  - phase: 16-01
    provides: benchmarks/ CMakeLists.txt with FetchContent, bench_common target
provides:
  - ModelFactory class for constructing arbitrary n-DOF ForwardDynamics/InverseDynamics solver objects
  - RandomState class for deterministic random joint state generation with seed 42
affects:
  - Phase 17 (ABA benchmarks)
  - Phase 18 (RNEA benchmarks)

# Tech tracking
tech-stack:
  added:
    - C++ <random> (std::mt19937 with seed 42, std::uniform_real_distribution)
  patterns:
    - Factory pattern for solver object construction (ModelFactory with 6 method overloads)
    - Deterministic RNG with fixed seed for reproducible benchmark inputs
    - #pragma once for header guards (LowerTriangular.h convention)
    - namespace SpatialAlgebra::Bench for benchmark utilities

key-files:
  created:
    - benchmarks/common/model_factory.h
    - benchmarks/common/model_factory.cpp
    - benchmarks/common/random_state.h
    - benchmarks/common/random_state.cpp
  modified: []

key-decisions:
  - "Followed plan exactly — no deviations. All D-01 through D-11 design decisions implemented as specified."
  - "JointConfig defaults: Z revolute axis, half-meter X translation, 1kg mass, COM at origin — per D-02/D-04."
  - "Branching implementation: link at branchPoint+2 has parent overwritten to branchPoint (both children of branchPoint have consecutive indices for the ABA tree invariant)."
  - "Removed _init.cpp placeholder from benchmarks/common/ — no longer needed now that model_factory.cpp and random_state.cpp exist."

patterns-established:
  - "Benchmark utility code lives in benchmarks/common/, in namespace SpatialAlgebra::Bench."
  - "Use std::uniform_real_distribution for constrained-range random generation with std::mt19937."
  - "Use std::move(link) on push_back for efficiency when constructing solver link vectors."

requirements-completed:
  - BINF-03

# Metrics
duration: 2min
completed: 2026-06-05
---

# Phase 16: Benchmark Infrastructure — Plan 02 Summary

**Shared benchmark utilities: ModelFactory for n-DOF solver construction and RandomState for deterministic random joint state generation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-06-05T13:16:15Z
- **Completed:** 2026-06-05T13:18:09Z
- **Tasks:** 3 (all auto)
- **Files modified:** 5 (4 created, 1 deleted)

## Accomplishments

- **ModelFactory** with 6 method overloads: createFD (serial), createID (serial), createFD (per-vector), createID (per-vector), createFDBranching, createIDBranching — covers serial chains, branching (Y-shaped) chains, and per-link custom configuration
- **JointConfig** struct with configurable joint axis (default: Z revolute), translation (default: 0.5m along X), mass (default: 1.0 kg), and COM (default: origin)
- **RandomState** class with fixed seed 42 (std::mt19937), randomPositions() in [-π/2, π/2], randomVelocities() in [-5, 5] rad/s
- Removed `_init.cpp` placeholder — real `.cpp` files now satisfy the CMake GLOB requirement

## Task Commits

Each task was committed atomically:

1. **Task 1: Create model_factory.h** — `dbfd607` (feat)
2. **Task 2: Create model_factory.cpp** — `82fa8cd` (feat)
3. **Task 3: Create random_state.h/.cpp** — `998c76a` (feat)

**Plan metadata:** (committed with SUMMARY.md)

## Files Created/Modified

- `benchmarks/common/model_factory.h` — JointConfig struct and ModelFactory class declaration with 6 methods
- `benchmarks/common/model_factory.cpp` — Full implementations for all 6 factory methods
- `benchmarks/common/random_state.h` — RandomState class declaration with seed 42
- `benchmarks/common/random_state.cpp` — randomPositions and randomVelocities implementations
- `benchmarks/common/_init.cpp` — **Deleted:** Placeholder no longer needed

## Decisions Made

- **Followed plan exactly** — no architectural decisions needed. All D-01 through D-11 design decisions from CONTEXT.md were implemented as specified.
- **Branching implementation:** Link at branchPoint+2 has parent overwritten to branchPoint, creating a Y-shaped chain where both children of branchPoint (indices branchPoint+1 and branchPoint+2) have consecutive indices, satisfying the ABA tree invariant.
- **_init.cpp removed:** The CMake GLOB now finds model_factory.cpp and random_state.cpp, making the placeholder redundant.

## Deviations from Plan

None - plan executed exactly as written.

**Total deviations:** 0
**Impact on plan:** N/A

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required. All code is C++17 with standard library dependencies.

## Next Phase Readiness

- Shared utilities are complete — ModelFactory and RandomState are ready for consumption by Phase 17 (ABA benchmarks) and Phase 18 (RNEA benchmarks).
- Benchmarks/common/CMakeLists.txt auto-discovers sources via GLOB — no CMake changes needed when adding new common utilities.
- Next: Plan 03 (bench_all executable with stub benchmark registrations and LTO profile).

---

*Phase: 16-benchmark-infrastructure*
*Completed: 2026-06-05*
