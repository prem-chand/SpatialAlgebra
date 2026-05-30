# Pitfalls Research

**Domain:** C++ Robotics Dynamics — Benchmarks, Stability & CI Compatibility
**Researched:** 2026-05-30
**Confidence:** HIGH (verified against codebase analysis, Eigen 5.x official docs, PCL issue #6351, RBDL/Pinocchio API docs, Google Benchmark documentation, Featherstone textbook analysis)

---

## Critical Pitfalls

### Pitfall 1: CR-02 Inward Pass — Bias Force (pa) Initialized Before Child Inertia Accumulation

**What goes wrong:**
The ABA inward pass in `src/ForwardDynamics.cpp` uses a two-phase approach that produces incorrect accelerations for 3+ link chains. Tested through `TestDynamicsConsistency.cpp` — the `ThreeLinkSerialChain` and `BranchingYConfiguration` tests fail because the round-trip RNEA(ABA(tau)) ≠ tau for multi-link chains. Single-link and two-link tests pass because there's only one accumulation step.

**Root cause:**
In `inwardPass()`:

```cpp
// Phase 1 (lines 56-73): Initialize ALL links
for (int i = 0; i < links.size(); i++) {
    Ia[i] = I_i;                        // Rigid body inertia only
    pa[i] = I_i * c_i + v_i × (I_i * v_i) + f_i;  // Uses I_i, not full I_A
}

// Phase 2 (lines 77-93): Accumulate children
for (int i = n-1; i >= 0; i--) {
    Ia[parent] += transform(Ia[i]);     // Correct: I_A accumulates child Ia
    pa[parent] += transform(pa[i]);     // WRONG: doesn't account for I_A × c
}
```

Per Featherstone Algorithm 7.3, each link's bias force **must** be computed using the **final** articulated inertia I_A (which includes all children):
```
p_A_i = I_A_i × c_i + v_i × (I_A_i × v_i) - f_i_ext
```

The current Phase 1 computes `pa[i]` using only the rigid body inertia `I_i`. When children are accumulated in Phase 2, the `I_A_i × c_i` term never gets recomputed with the full inertia. The transform accumulation of child pa only partially compensates — the missing cross-coupling term `(child_I_A_transformed) × c_i` is silently dropped.

**Why it's subtle:**
- Single and 2-link chains work because there's only one accumulation hop — the missing term is numerically small or zero for identity inertias at zero configuration
- 3+ link chains expose the cumulative error because each intermediate link must compute p_A using its own I_A, not just additive child contributions
- The bug only manifests in dynamic scenarios (Coriolis, coupled accelerations) — static/gravity-only tests may appear correct
- All 156 passing tests don't catch it because they're single-link or 2-link only

**Consequences:**
- Multi-link round-trip consistency fails: RNEA(ABA(tau)) ≠ tau for chains with 3+ links
- Branching trees (Y-configuration) get wrong acceleration distribution
- Forward dynamics is unreliable for any robot with more than 2 joints
- Any downstream use (simulation, MPC, control) would produce incorrect physics

**Prevention:**
Restructure the inward pass as a single pass per link (tip to base):

```cpp
for (int i = n-1; i >= 0; i--) {
    // Start with rigid body inertia
    Ia_i = I_i;
    pa_i = -cross(v_i, I_i * v_i);   // Velocity product term only (no I_A*c yet)

    // Add children (already fully computed since we go tip→base)
    for (int child : children_of(i)) {
        Ia_i += transform(Ia_child);
        pa_i += transform(pa_child);
    }

    // NOW compute bias force with FULL I_A (includes children)
    pa_i += Ia_i * c_i;
    pa_i -= f_i_ext;                  // External forces

    // If not base, transform to parent frame (handled in next iteration)
}
```

**Detection:**
- `TestDynamicsConsistency::ThreeLinkSerialChain` fails (test exists)
- `TestDynamicsConsistency::BranchingYConfiguration` fails (test exists)
- `ForwardDynamicsTest::ThreeLinkNumericalValidation` has a `TODO(CR-02)` marker at line 309 expecting `qddot[0] < qddot[2]`

**Which phase to address:** Phase: "Fix CR-02 ABA Inward Pass" — must come BEFORE benchmarks or examples, since benchmarking a buggy solver produces meaningless numbers.

---

### Pitfall 2: Eigen 5.x CMake `find_package` Version Pin Breakage

**What goes wrong:**
The current `CMakeLists.txt:12` uses:
```cmake
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
```

The `3.3` version argument restricts `find_package` to Eigen versions ≥3.3.0 but <4.0.0 in CMake's default version-compatibility logic. Eigen 5.0.0 was released with new semantic versioning (jumping from 3.4 to 5.0), which means CMake's range-check rejects it because:
- 5.0.0 ≥ 3.3.0 ✓ (lower bound satisfied)
- 5.0.0 < 4.0.0 ✗ (upper bound implied by `3.3` format fails)

This is **exactly** the issue documented in `PointCloudLibrary/pcl#6351` — a known breakage pattern affecting projects that pin Eigen 3.x.

**Why it's subtle:**
- `brew install eigen` on macOS typically installs the latest (5.x as of 2026)
- Local builds work because the user has `Eigen3_DIR` set or uses the system-installed 3.4.x
- CI builds work because `apt-get install libeigen3-dev` on Ubuntu 24.04 still installs 3.4.x
- The breakage is silent until Homebrew updates, then `cmake -B build` fails with:
  ```
  Could NOT find Eigen3 (missing: Eigen3_DIR)
  ```
  — which looks like an installation issue, not a version issue

**Consequences:**
- Build breakage on any machine with Eigen 5.x installed
- Devs waste time debugging "Eigen not found" when it IS installed — it's just version-incompatible
- Inconsistent behavior between CI (Ubuntu 3.4.x) and local dev (macOS 5.x)

**Prevention:**
Replace the version-pinned `find_package` with a range-aware approach:

```cmake
# Option A: Remove version pin entirely (if no API-breaking changes)
find_package(Eigen3 REQUIRED NO_MODULE)

# Option B: Use version range (supported since Eigen 3.4.1, CMake 3.19+)
find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)
# This accepts: >=3.4.1 AND <6.0.0
# Which covers both 3.4.x and 5.x

# Option C: Test version at compile time (most robust for API-sensitive code)
find_package(Eigen3 REQUIRED NO_MODULE)
target_compile_definitions(SpatialAlgebra PRIVATE
    EIGEN_VERSION_CHECKED=1
)
# Guard Eigen 5.x- specific workarounds with:
# #if EIGEN_VERSION_AT_LEAST(5,0,0)
```

**Detection:**
- CI should have at least one matrix entry explicitly installing Eigen 5.x (e.g., `brew install eigen@5` or fetching from GitLab)
- The `cmake -B build` step should fail loudly if there's a version issue, not silently find the wrong version
- A compile-time check in a header: `#if EIGEN_VERSION_AT_LEAST(3, 3, 0)` guard for required features

**Which phase to address:** Phase: "Eigen 5.x CI Compatibility" — can be done in parallel with CR-02 fix.

---

### Pitfall 3: Frame Convention Mismatch in RBDL/Pinocchio Comparison Benchmarks

**What goes wrong:**
When comparing SpatialAlgebra against RBDL or Pinocchio, naive cross-validation produces different numerical results even when both libraries implement the same algorithm (Featherstone ABA/RNEA). The comparison code silently measures "difference due to convention" rather than "difference due to performance/correctness."

**Root cause — three independent convention axes:**

1. **Motion/Force vector ordering:**
   - SpatialAlgebra: `[angular; linear]` — angular components first (Featherstone standard)
   - RBDL: Also `[angular; linear]` per Featherstone
   - Pinocchio: Also `[angular; linear]` in SE3 notation
   - *Risk: LOW* — all three use Featherstone ordering, but MUST verify when constructing from raw arrays

2. **Transform direction:**
   - SpatialAlgebra: `X` transforms parent → child (motion), `X^{-T}` transforms child → parent (force)
   - RBDL: `SpatialTransform` with `E` and `r`. `apply()` transforms child → parent for motion vectors
   - Pinocchio: `SE3` placement. `act()` transforms from joint frame to parent frame
   - *Risk: HIGH* — each library has a different API for "which direction does this transform go?"
   - Example: RBDL's `X.apply(v)` = parent-to-child motion in SpatialAlgebra terms

3. **Joint screw axis (S) convention:**
   - SpatialAlgebra: MotionVector with angular and linear parts. Revolute Z = `(0,0,1, 0,0,0)`
   - RBDL: Joint axis defined by `JointType` (RevoluteZ, RevoluteX, etc.) with implicit screw
   - Pinocchio: Joint model defines S via `JointModel` derived classes (RX, RY, RZ)
   - *Risk: MEDIUM* — all are Featherstone-based, so cross-axis mapping is 1:1, but floating-base and multi-DOF joints add complexity

4. **Floating base vs. fixed base:**
   - SpatialAlgebra: Fixed-base only (parent=-1 for base, no 6-DOF root joint)
   - RBDL: Supports both. Floating base uses a 6-DOF joint at root
   - Pinocchio: Same as RBDL — `JointModelFreeFlyer` for floating base
   - *Risk: HIGH* — comparing a fixed-base solver against a library's floating-base output will produce different torque/acceleration values

**Consequences:**
- Comparison benchmark "fails" or produces misleading ratios
- Cross-validation test "discovers" a bug that's actually a convention difference
- Developer wastes time debugging non-existent issues
- Published benchmark numbers are not reproducible by others

**Prevention:**

1. **Establish exact equivalence before running benchmarks.** Create a test that:
   - Builds identical kinematic chain (same masses, inertias, transforms, joint axes)
   - Runs RNEA with identical inputs in both libraries
   - Verifies torque output matches to machine epsilon (1e-12)
   This test validates convention alignment, not dynamics correctness.

2. **Document the convention mapping:**
   ```
   SpatialAlgebra → RBDL
   ===========================
   X.transformMotion(v) → X.apply(v)  // verify direction!
   inverseTransformForce(f) → X.applyTranspose(f)  // may need transpose
   S = (axis, 0) → JointType::RevoluteZ  // check axis mapping
   ```

3. **Use identical model files.** If both libraries support URDF, load the same URDF rather than manually constructing links in each API.

4. **Compare ratios, not absolute values.** If conventions differ, compare `ratio_of_accelerations` rather than `acceleration[0]` directly.

**Detection:**
- Comparison benchmark produces identical results for single-link but diverges for multi-link
- Sign flips between libraries (negative vs positive torque for same input)
- Gravity compensation torques differ systematically

**Which phase to address:** Phase: "RBDL/Pinocchio Comparison Benchmarks" — after CR-02 fix (otherwise buggy solver corrupts comparison).

---

### Pitfall 4: Benchmarking an Unfixed Bug

**What goes wrong:**
Performance benchmarks are added while CR-02 (multi-link ABA bug) is still present. The benchmark numbers measure "how fast the wrong answer is computed" — producing meaningless metrics that waste everyone's time.

**Why it happens:**
- Benchmarking is a "safe" task that can be parallelized
- The bug fix seems complex, benchmarks seem simple
- "We can replace numbers later" — they never get replaced
- Published benchmarks with wrong physics damage library credibility

**Consequences:**
- All multi-link benchmark results invalidated when CR-02 is fixed
- Performance regression from fix appears as "we got slower" — but actually the old number was computing the wrong thing
- Comparison against RBDL/Pinocchio shows SpatialAlgebra "faster but wrong" — which is worse than useless

**Prevention (hard rule):**
```
BENCHMARKS MUST NOT BE ADDED UNTIL CR-02 IS FIXED AND ALL 4 CONSISTENCY TESTS PASS.
```

The acceptance gate is: `TestDynamicsConsistency.cpp` must show 0 failures (4/4 tests passing).

**Detection:**
- If `CTest` output shows `3 tests FAILED` in any consistency test, benchmarks are premature
- The `ThreeLinkNumericalValidation` test's `TODO(CR-02)` comment at line 309 is the sentinel

**Which phase to address:** Phase ordering enforcement — the CR-02 fix phase must complete verification before the benchmarks phase begins.

---

## Performance Traps

| Trap | Symptoms | Prevention |
|------|----------|------------|
| **Debug mode benchmarks** | Unrealistically slow numbers (10-100x). Algorithmic O(n) appears O(n²) due to assertions. | Always benchmark Release builds (`-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG"`). |
| **Compiler optimizing away benchmarked code** | Zero-time reported. Loop body eliminated as dead code. | Use `benchmark::DoNotOptimize(result)` and `benchmark::ClobberMemory()` for every output value. Never return unused values. |
| **Cold cache vs warm cache** | First iteration 10x slower than median. Benchmark reports "average" masking startup cost. | Run 5+ warmup iterations before timing. Google Benchmark does this automatically. For manual benchmarks, discard first N measurements. |
| **Including setup in timing** | Benchmark measures setup+computation, overstating real cost. Vector/matrix allocation dominates. | Move allocation outside timed loop. Reserve capacity. Use `PauseTiming()`/`ResumeTiming()` for setup. |
| **CPU frequency scaling** | Noisy results (±20%). Turbo Boost kicks in/out during run. | Pin CPU frequency if possible (`cpupower frequency-set -g performance`). Run multiple repetitions and check min/max spread. |
| **Dynamic memory in hot path** | Consistent 2-3x slowdown from `malloc` in ABA loop. Benchmarks hide allocation in "one-time setup." | Pre-allocate all link vectors. Profile with `perf` to check for `malloc` calls in hot path. SpatialAlgebra already stores links in `std::vector` — ensure no reallocation. |
| **Small-N benchmark misleading extrapolation** | O(n) looks O(1) for n=3. Benchmark only tests 3-link chains. | Test across N = {2, 3, 6, 10, 20, 50} links to establish scaling curve. |
| **Compiler autovectorization differences** | Clang 2x faster than GCC due to different SLP vectorizer behavior on Eigen expressions. | Test with both compilers. Report per-compiler results separately. Don't average across compilers. |
| **Google Benchmark vs manual `std::chrono`** | Manual timing has ±5-15% noise floor. Google Benchmark reduces to ±1-3% with statistical iteration control. | Always use a proper framework (Google Benchmark) for microbenchmarks. Manual timing is only acceptable for coarse (algorithmic) measurements. |

---

## "Looks Done But Isn't" Checklist

### CR-02 Bug Fix
- [ ] **All 4 consistency tests pass** — `TestDynamicsConsistency` shows 0 failures, not 3
- [ ] **`ThreeLinkNumericalValidation`** — `EXPECT_NE(fd.links[0].qddot, fd.links[2].qddot)` passes (base ≠ tip)
- [ ] **Branching tree symmetry** — `fd.links[1].qddot ≈ fd.links[2].qddot` in Y-configuration
- [ ] **Round-trip with random states** — not just zero-initialized: test with non-zero `q`, `qdot`, random inertias
- [ ] **Existing 156 tests still pass** — fix doesn't break single/two-link behavior
- [ ] **NaN/Inf guards still trigger** — debug-mode assertions not accidentally removed

### Performance Benchmarks
- [ ] **Built with `-O3 -DNDEBUG`** — verified via `CMAKE_BUILD_TYPE=Release`, not just assumed
- [ ] **Google Benchmark (or equivalent)** — `DoNotOptimize`/`ClobberMemory` on all outputs
- [ ] **Warm cache results** — at least 5 warmup iterations discarded
- [ ] **Multiple link counts tested** — {2, 3, 6, 10, 20} minimum to establish O(n) scaling
- [ ] **Per-compiler results separated** — not averaged across g++/clang++
- [ ] **Statistical significance** — coefficient of variation < 5% across runs
- [ ] **Benchmark code in separate file** — not mixed with test code in `tests/`
- [ ] **No `std::cout` or logging in hot path** — I/O operations skew timing

### Eigen 5.x CI Compatibility
- [ ] **CI matrix includes Eigen 5.x** — at least one job installs Eigen 5.x (not just whatever `apt` provides)
- [ ] **`find_package` without version pin** — or uses range syntax `3.4...5`
- [ ] **Compile-time version check** — `EIGEN_VERSION_AT_LEAST(5,0,0)` guards for any API differences
- [ ] **Both macOS and Linux tested** — Homebrew tends to be ahead of apt in Eigen version
- [ ] **No `#pragma GCC` or `#pragma clang` warnings** from Eigen 5.x headers (newer Eigen adds `[[nodiscard]]` etc.)

### RBDL/Pinocchio Comparison
- [ ] **Frame convention verified** — match on single-link exact numerical equivalence before multi-link comparison
- [ ] **Identical model used** — same masses, inertias, transforms, axes in both libraries
- [ ] **RNEA comparison matches to 1e-12** — not just "close enough"
- [ ] **Floating base handled explicitly** — SpatialAlgebra is fixed-base only; comparison must use fixed-base models
- [ ] **Gravity vector identical** — both libraries get same `g = (0, 0, -9.81)` in same frame
- [ ] **Performance comparison fair** — same compiler flags, same CPU, same warmup strategy

---

## Phase-Specific Warnings

| Phase | Likely Pitfall | Severity | Mitigation |
|-------|---------------|----------|------------|
| CR-02 Bug Fix | Fixing only symptoms (swapping sign) instead of root cause (pa initialization order) | CRITICAL | Verify against Featherstone Algorithm 7.3 pseudocode line-by-line. Single-pass restructure, not band-aid. |
| CR-02 Bug Fix | Breaking existing 156 passing tests | HIGH | Run full test suite after fix. Add `ThreeLinkRandomState` test. |
| CR-02 Bug Fix | Introducing new NaN/Inf edge case in restructured inward pass | MEDIUM | Keep debug-mode assertions on all intermediate values in the new single-pass loop. |
| Eigen 5.x CI | Adding Eigen 5.x to CI matrix but using same version-pinned `find_package` — CI fails silently | HIGH | Remove/patch version pin first, then add CI job. |
| Eigen 5.x CI | Eigen 5.x API changes ([[nodiscard]], deleted copy ops) causing compilation errors in SpatialAlgebra headers | MEDIUM | Build SpatialAlgebra against Eigen 5.x HEAD before CI configuration. Check `Eigen/src/Core/util/Macros.h` for version. |
| Benchmarks | Mixing benchmark and test executables in same CMake target | LOW | Add separate `benchmarks/` directory with its own CMakeLists.txt. Don't link benchmark library into test executables. |
| Benchmarks | Publishing benchmark numbers with CR-02 unfixed | CRITICAL | GATE: benchmarks phase must be AFTER CR-02 fix phase in ROADMAP.md ordering. |
| Benchmarks | Benchmark only forward dynamics, not inverse dynamics | MEDIUM | Both ABA and RNEA need benchmarks (users need both). |
| RBDL Comparison | Using incompatible RBDL version (RBDL 1.0 vs 2.0 have different APIs) | MEDIUM | Pin specific RBDL version in FetchContent. Document API version tested. |
| RBDL/Pinocchio Comparison | Not accounting for missing floating base support in SpatialAlgebra | HIGH | Only compare serial chains with fixed-base models. Document "fixed-base only" in benchmark methodology. |
| Real-world Examples | Examples that look like they work but use buggy CR-02 behavior | HIGH | Only add examples AFTER CR-02 fix. Validate examples produce physically correct output (energy conservation, gravity compensation). |

---

## Cross-Cutting Concerns

### Backward Compatibility During Fixes
- CR-02 fix changes `Link::Ia` and `Link::pa` intermediate values — any external code that reads these mid-solver will see different values
- Public API (`computeAccelerations`) unchanged — input/output contract preserved
- `ForwardDynamicsLink` type alias unchanged — user code that constructs links is unaffected

### Featherstone Reference Alignment
- Current implementation deviates from Algorithm 7.3 in the inward pass structure (two-phase vs single-pass)
- After fix, the implementation should match Algorithm 7.3 line-by-line, making it easier to verify against the textbook
- Keep the textbook pseudocode as a comment block above `inwardPass()` for future maintainers

### Eigen Version Detection Strategy
During the transition period where both Eigen 3.4.x and 5.x are in common use:
```cpp
#include <Eigen/src/Core/util/Macros.h>  // Defines EIGEN_WORLD_VERSION etc.
#if EIGEN_WORLD_VERSION >= 5
// Eigen 5.x code paths
#else
// Eigen 3.4 code paths
#endif
```
This is preferable to CMake version detection because it's compile-time and can't get out of sync with what's actually included.

---

## Sources

| Finding | Source | Confidence |
|---------|--------|------------|
| CR-02 root cause (pa before Ia accumulation) | Code analysis: `src/ForwardDynamics.cpp:56-93` | HIGH |
| Eigen 5.x find_package incompatibility | [PCL Issue #6351](https://github.com/PointCloudLibrary/pcl/issues/6351) + Eigen docs: [TopicCMakeGuide](https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html) | HIGH |
| Eigen version range syntax 3.4...5 | [Eigen Nightly Docs](https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html) | HIGH |
| RBDL spatial algebra implementation | [RBDL DeepWiki: Mathematical Foundation](https://deepwiki.com/rbdl/rbdl/3.6-mathematical-foundation) | MEDIUM |
| Pinocchio Featherstone foundation | [Pinocchio docs](https://stack-of-tasks.github.io/pinocchio/) | MEDIUM |
| RBDL/Pinocchio convention mismatch reports | [Pinocchio Issue #1721](https://github.com/stack-of-tasks/pinocchio/issues/1721) | MEDIUM |
| Google Benchmark optimization barriers | [Google Benchmark docs](https://github.com/google/benchmark) | HIGH |
| C++ microbenchmark pitfalls | Multiple sources (Google Benchmark docs, SO, Ash Vardanian) | HIGH |
| Featherstone Algorithm 7.3 structure | Featherstone (2008) Rigid Body Dynamics Algorithms, §7.2.1 | HIGH |
| Current CI matrix (4-matrix, no Eigen 5.x) | `.github/workflows/ci.yml` | HIGH |
| Current test gaps (156/158 passing, 3 CR-02 failures) | `TestDynamicsConsistency.cpp` + codebase analysis | HIGH |
