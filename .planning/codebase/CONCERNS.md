# Codebase Concerns

**Analysis Date:** 2026-05-17

## Tech Debt

### CR-01: Swapped Constructor Arguments in `ArticulatedBodyInertia::operator+(RigidBodyInertia)`

- **Issue:** The `operator+(const RigidBodyInertia &)` method passes arguments to the `ArticulatedBodyInertia` constructor in the wrong order. The first argument receives the `M` (mass matrix) field instead of `Inertia` (rotational inertia), and the third argument receives `Inertia` instead of `M`. Additionally, the coupling term is missing the `mass` multiplier — it uses `skew(other.getCom())` instead of `other.getMass() * skew(other.getCom())`.
- **Files:**
  - `include/ArticulatedBodyInertia.h:150-155`
- **Impact:** Combining articulated body inertias with rigid body inertias (e.g., during ABA) produces numerically incorrect results in all components. This affects any code path that accumulates inertias.
- **Fix approach:** Swap the first and third constructor arguments and add the mass multiplier to the coupling term.

### CR-02: ABA Forward Dynamics — Inward Pass Overwrites Accumulated Child Inertias

- **Issue:** The Articulated Body Algorithm (ABA) inward pass (`ForwardDynamics::inwardPass` at `src/ForwardDynamics.cpp:54-107`) initializes `links[i].Ia` and `links[i].pa` for each link while iterating tip-to-base, but this **overwrites** the contributions that child links already added to their parent. For a 3-link chain (0→1→2): link 2's contribution is added to link 1, then link 1 re-initializes itself (destroying link 2's contribution), then propagates only its own (wrong) inertia to link 0. The result: the ABA computes accelerations as if each link is isolated, ignoring all downstream inertias for non-terminal links.
- **Files:** `src/ForwardDynamics.cpp:54-107`
- **Impact:** Multi-link forward dynamics are incorrect. The single-link pendulum test passes (no children to accumulate). The two-link test only checks "positive and finite" which masks the bug. The round-trip consistency tests compound errors. **This makes ABA effectively broken for all multi-link systems.**
- **Fix approach:** Restructure the inward pass with a two-phase approach: initialize all Ia/pa from rigid body inertia first, then accumulate from children to parents without re-initialization.

### CR-03: Three Mutually-Inconsistent Implementations of Force×Force Cross Product

- **Issue:** The mathematically correct force-force cross product formula is `[τ1×τ2 + f1×f2; τ1×f2 - τ2×f1]` (anti-commutative). It has three separate implementations, **none of which match the correct formula**:
  - `SpatialVector::crossForce` (`src/SpatialVector.cpp:48-55`): Angular = ω1×ω2 + v1×v2 (correct), Linear = ω1×v2 — **missing `- v2×ω1`**
  - `ForceVector::crossForce` (`src/ForceVector.cpp:40-46`): Angular = τ1×τ2 — **missing `f1×f2`**, Linear = τ1×f2 - τ2×f1 (correct)
  - `cross(ForceVector,ForceVector)` (`include/SpatialUtils.h:115-125`): Angular = τ1×τ2 + f1×f2 (correct), Linear = τ1×f2 — **missing `- τ2×f1`**
- **Files:**
  - `src/SpatialVector.cpp:48-55`
  - `src/ForceVector.cpp:40-46`
  - `include/SpatialUtils.h:115-125`
- **Impact:** Any code using cross-product operations with force vectors will produce incorrect results that depend on which overload happens to be called. Violates anti-commutativity, a fundamental physical property.
- **Fix approach:** Unify all three to use the single correct formula: `angular = t1.cross(t2) + f1.cross(f2)`, `linear = t1.cross(f2) - t2.cross(f1)`.

### Swapped Constructor Arguments in ArticulatedBodyInertia (Documented but Unfixed)

- **Files:** `include/ArticulatedBodyInertia.h:95-96`
- **Issue:** Same bug as CR-01 but already identified in a prior review (`REVIEW.md:73-102`). The code is in a verified-known-broken state.
- **Fix approach:** See CR-01 fix.

### `MotionVector::crossForce` Has No Physical Meaning

- **Issue:** `MotionVector::crossForce` computes `[ω1×v2; v1×ω2]` which does not correspond to any standard spatial algebra operation. The force cross product (`crf`) is defined as `motion × force → force`. Applying `crossForce` to two motion vectors has no physical meaning. Additionally, `ForceVector::crossMotion` applies the motion-cross formula to two force vectors, which is a category error.
- **Files:**
  - `src/MotionVector.cpp:38-41`
  - `src/ForceVector.cpp:31-38`
  - `include/MotionVector.h:131-137` (Doxygen comment says "crf operation" but returns `MotionVector`)
- **Impact:** These methods compile and run, but produce mathematically and physically meaningless results. They mislead users into thinking there's a valid operation here.
- **Fix approach:** Either remove `MotionVector::crossForce` and `ForceVector::crossMotion`, or mark them as deprecated with clear documentation that they have no physical interpretation.

### `PluckerTransform::multiply` Has `auto` Return Type — Linkage Hazard

- **Issue:** `multiply` is declared with `auto` return type in the header (`include/PluckerTransform.h:172`) and defined in the `.cpp` file (`src/PluckerTransform.cpp:228-239`). In C++17, `auto` return type deduction requires the definition to be visible at the call site. Calling `multiply` from any other translation unit will fail to compile/link. Currently it's only called from `apply()` within the same `.cpp`, so it works — but as a public API this is a trap.
- **Files:** `include/PluckerTransform.h:172`, `src/PluckerTransform.cpp:228-239`
- **Impact:** Any user including `PluckerTransform.h` and calling `multiply()` will get a compile/link error.
- **Fix approach:** Change the return type from `auto` to `PluckerTransform` explicitly.

### `Rotation::transpose()` Invokes Method Without Explicit Object

- **Issue:** `Rotation::transpose()` uses `Eigen::Matrix3d::transpose()` which relies on implicit `*this` to call a non-static member function. This works only because we're inside a member function of the derived class, but it's fragile — minor changes to Eigen's API or context could silently break. The `inverse()` method correctly uses `this->transpose()`.
- **Files:** `src/Rotation.cpp:54-57`
- **Impact:** Fragile code that could silently produce wrong results or fail to compile with Eigen version updates.
- **Fix approach:** Use `Eigen::Matrix3d::transpose()` instead of `Eigen::Matrix3d::transpose()` — wait, that's the issue. Use `this->transpose()` consistently.

### Test Helper Functions Return Zero Instead of Identity/Diagonal

- **Issue:** Helper functions `createIdentityInertia()` and `createDiagonalInertia(double)` in `tests/TestSpatialOperations.cpp` create a `LowerTriangular(3)` (zero-initialized) and return it without setting any values. The data vector they construct is dead code — the return value is always the zero matrix. Five tests use these helpers and only verify `getMass()` and `getCom()`, not the inertia matrix values, so they pass despite using wrong input inertias.
- **Files:** `tests/TestSpatialOperations.cpp:27-37, 206-218`
- **Impact:** Test coverage is illusory — these tests verify minimal properties of inertia transforms while using zero inertia inputs. A zero inertia matrix is physically invalid for a body with mass.
- **Fix approach:** Actually populate the LowerTriangular matrix with identity/diagonal values before returning.

### Global Namespace Pollution with `using Vector3d`

- **Issue:** Six header files each declare `using Vector3d = Eigen::Matrix<double, 3, 1>` at **global namespace scope** (before the `namespace SpatialAlgebra` block). This leaks the alias into the global namespace of every translation unit that includes any SpatialAlgebra header. Only one declaration is needed, and it should be inside `namespace SpatialAlgebra`.
- **Files:**
  - `include/SpatialVector.h:44`
  - `include/LowerTriangular.h:54`
  - `include/PluckerTransform.h:54`
  - `include/RigidBodyInertia.h:17`
  - `include/InverseDynamics.h:59`
  - `include/ForwardDynamics.h:63`
- **Impact:** ODR violations or ambiguity if users have their own `Vector3d` type, or if Eigen ever adds a `Vector3d` typedef.
- **Fix approach:** Move the `using` declaration inside `namespace SpatialAlgebra` in a single header (e.g., `SpatialVector.h` or a common base header).

### `MotionVector::getAngular/getLinear` Shadow Parent Methods Inconsistently

- **Issue:** `MotionVector` and `ForceVector` declare their own `getAngular()` and `getLinear()` methods that return by value (hiding `SpatialVector::getAngular()/getLinear()` which return by const reference). There is no additional behavior or type narrowing — these are redundant and create confusion about which method is called through a base-class reference.
- **Files:**
  - `include/MotionVector.h:101-107` (getters not explicitly shown but inherited from SpatialVector — MotionVector.cpp has no getter overrides)
  - `include/ForceVector.h` (same)
- **Impact:** Hybrid dispatch behavior: calling through a `SpatialVector&` returns const ref, calling through `MotionVector&` returns by value. This can cause subtle lifetime issues.
- **Fix approach:** Remove the redundant overrides and let the base class methods be inherited.

### Empty `src/ArticulatedBodyInertia.cpp` and `src/RigidBodyInertia.cpp`

- **Issue:** Both files contain only a comment saying "No implementation needed — all methods are inline in header". These are compiled into the library, producing empty object files. They serve no purpose and are included by CMake's `file(GLOB SOURCES "src/*.cpp")`.
- **Files:**
  - `src/ArticulatedBodyInertia.cpp` (8 lines, no code)
  - `src/RigidBodyInertia.cpp` (8 lines, no code)
- **Impact:** Trivial — small waste of compilation time. Indicates incomplete migration from header-only to separate compilation.
- **Fix approach:** Remove the empty files or add real implementations.

## Known Bugs

### ABI Transform (tformABI/invtformABI) Formula Errors — 5 Failing Tests

- **Symptoms:** Five tests fail in `TestPluckerTransform.cpp`:
  - `TransformABITest.Property_Symmetric`
  - `InverseTransformABITest.InverseIsIdentity`
  - `InverseTransformABITest.Property_Symmetric`
  - `InverseTransformABITest.RoundTrip`
  - `TestInverse.MultiplyWithInverseIsIdentity`
- **Files:** `src/PluckerTransform.cpp:104-161, 163-213`
- **Trigger:** Transforming articulated body inertias with non-zero coupling matrices (H ≠ 0).
- **Root Cause:** Formula derivation for `tformABI()` and `invtformABI()` does not correctly implement Featherstone Eq 7.16. The current implementation has incorrect handling of the coupling matrix H transformation and sign errors in the inverse transform formula.
- **Workaround:** For simple chains with H = 0 (no rotational-linear coupling), the algorithm works correctly. Avoid transforming ABI with non-zero H.
- **Fix approach:** Re-derive the block formula based on Featherstone Eq 7.16 and reimplement. Verification required against known test cases.
- **Verification state:** Identified in `v1.0-VERIFICATION.md` as known issue. Unfixed as of 2026-05-17.

### Multi-Link Dynamics Consistency Between RNEA and ABA — 2 Failing Tests

- **Symptoms:** `TestDynamicsConsistency.cpp` tests fail for multi-link systems:
  - `ConsistencyTest.ThreeLinkSerialChain`
  - `ConsistencyTest.BranchingYConfiguration`
- **Files:** `tests/TestDynamicsConsistency.cpp`
- **Trigger:** Cross-validation between inverse dynamics (RNEA) and forward dynamics (ABA) for multi-link chains.
- **Root Cause:** RNEA and ABA use different conventions for joint acceleration computation, different bias acceleration handling, and may have different joint axis conventions. The CR-02 bug (ABA overwriting child inertias) compounds this.
- **Workaround:** Single-link consistency verified and working. Do not rely on RNEA↔ABA consistency for multi-link systems.
- **Fix approach:** Fix CR-02 first, then align RNEA/ABA conventions. Verify with the consistency tests.

## Security Considerations

### No Input Validation on SpatialVector Operations

- **Risk:** Spatial vector constructors and arithmetic operations do not validate their inputs. NaN or Inf values propagate silently through the system. While this is a mathematical library (not a network-facing service), NaN propagation can cause silent failures in dependent applications (e.g., robot controllers making unsafe movements).
- **Files:** All `SpatialVector` operations in `src/SpatialVector.cpp`, `src/MotionVector.cpp`, `src/ForceVector.cpp`
- **Current mitigation:** Input validation (`std::isnan`/`std::isinf` checks) exists only in `ForwardDynamics::computeAccelerations` and `InverseDynamics::computeTorques`. No validation in core vector/inertia operations.
- **Recommendation:** Add debug-mode assertions for NaN/Inf in critical arithmetic operations. Consider a validation layer for production use.

### No Bounds Checking in Release Mode

- **Risk:** `LowerTriangular` has bounds checking guarded by `#ifndef NDEBUG`, meaning it's disabled in release builds. If a user accesses an out-of-bounds element in release mode, the behavior is undefined (silent memory corruption via Eigen's `VectorXd::operator[]`).
- **Files:** `include/LowerTriangular.h:146-173`
- **Current mitigation:** Debug-mode bounds checking only.
- **Recommendation:** Consider `std::out_of_range` throws unconditionally, or document the performance/security tradeoff.

### Potential Forthcoming Eigen 5.x Incompatibility

- **Risk:** CMakeLists.txt uses `find_package(Eigen3 3.3 REQUIRED NO_MODULE)` which may fail with Eigen 5.x because Eigen 5 changed its CMake version-compatibility range. Users on Homebrew (which now ships Eigen 5.0.1) will get a build error.
- **Files:** `CMakeLists.txt:12`
- **Current mitigation:** Either remove the version pin or set `-DEigen3_DIR=$(brew --prefix eigen)/share/eigen3/cmake`. Documented in `AGENTS.md`.
- **Recommendation:** Remove the version pin or bump to support both 3.3+ and 5.x.

## Performance Bottlenecks

### ABI Transform Uses Full 6×6 Matrix Construction Instead of Block Operations

- **Problem:** `PluckerTransform::tformABI` and `invtformABI` build full 6×6 matrices using `MatrixXd::Zero(6,6)` then extract 3×3 blocks. This involves heap allocation for every inertia transform and O(6³) operations when O(3³) block operations would suffice.
- **Files:** `src/PluckerTransform.cpp:104-161, 163-213`
- **Cause:** Implementation shortcuts — constructing the full 6×6 matrix and using dense multiplication is simpler to code but far from optimal.
- **Improvement path:** Implement direct block-wise formulas for ABI transformation. The three components (I, H, M) can each be transformed using 3×3 matrix multiplications and skew operations, avoiding the 6×6 overhead.

### OpenMP Dependency with Unclear Availability

- **Problem:** `LowerTriangular::operator*(const LowerTriangular&)` uses `#pragma omp parallel for collapse(2)`. OpenMP support on macOS (Apple Clang) is not available by default — users must install `libomp` via Homebrew. The CMakeLists.txt does not check for or enable OpenMP, so this pragma is silently ignored on most macOS builds. On platforms where it does work, it adds threading overhead for small matrices (n=3 typically) where it's unnecessary.
- **Files:** `include/LowerTriangular.h:202-203`
- **Cause:** The `LowerTriangular` class was designed for general n×n matrices, but in practice it's always 3×3 (inertia tensors). The OpenMP parallelization is overkill for 3×3 matrices and may actually be slower due to thread spawning overhead.
- **Improvement path:** Remove OpenMP dependency entirely (it's unused in practice for n=3). Specialize operations for 3×3 case with explicit unrolled loops.

### Redundant Cross Product Implementations

- **Problem:** The force×force cross product is implemented (incorrectly) in three separate places. Even after fixing the correctness issue, having three implementations means three code paths to maintain and optimize. Each does essentially the same work.
- **Files:**
  - `src/SpatialVector.cpp:48-55`
  - `src/ForceVector.cpp:40-46`
  - `include/SpatialUtils.h:115-125`
- **Cause:** Lack of a single canonical implementation. The base class `SpatialVector` provides a version, derived classes override, and free functions provide yet another.
- **Improvement path:** Pick one canonical implementation (recommended: the free function in `SpatialUtils.h`) and delegate all other implementations to it.

## Fragile Areas

### ABA Forward Dynamics Implementation

- **Files:** `src/ForwardDynamics.cpp` (entire file, 158 lines)
- **Why fragile:** The ABA implementation has a known critical bug (CR-02) that makes multi-link results incorrect. The bias acceleration handling (WR-05 from REVIEW.md) is coupled with the Ia initialization issue — fixing one without the other will not produce correct results. The algorithm has no gravity term (`gravityMode` is not implemented), and this omission is undocumented. The `inwardPass` method does both accumulation and solve in one phase, making the logic hard to verify.
- **Test coverage:** 6 tests, all single-link or verify only "positive and finite". The multi-link tests pass despite the bug. Tests do not verify correct numerical values for multi-link scenarios.
- **Safe modification:** Refactor the inward pass to clearly separate initialization, accumulation, and solve phases. Add gravity support as a separate concern.

### RNEA Inverse Dynamics Implementation

- **Files:** `src/InverseDynamics.cpp` (entire file, 141 lines)
- **Why fragile:** The inward pass (`inwardPass`) iterates through all links for each link to find children (`for (int j = 0; j < links.size(); j++)`). This is O(n²) instead of the expected O(n) for RNEA. The algorithm assumes links are ordered with children after parents but does not enforce or verify this. There is no gravity term (documented as "add external forces").
- **Test coverage:** 5 tests, all single-link or two-link. No tests with non-zero velocity (Coriolis/centrifugal terms are computed but not verified).
- **Safe modification:** Restructure the inward pass to use a precomputed child list for O(n) performance. Add gravity support. Validate link ordering.

### LowerTriangular `inverse()` Implementation

- **Files:** `include/LowerTriangular.h:424`, `src/LowerTriangular.cpp:25-53`
- **Why fragile:** The inverse uses forward substitution with a hardcoded singularity threshold of `1e-15`. This threshold is arbitrary and may be too strict for some applications or not strict enough for others. There's no documentation about the expected numerical accuracy. The method throws `std::runtime_error` for singular matrices but doesn't provide diagnostic information about the condition number.
- **Test coverage:** 5 tests covering inverse correctness, diagonal elements, identity, 2×2 case. No tests for near-singular matrices or numerical stability.
- **Safe modification:** Document the singularity threshold. Consider using a relative tolerance based on matrix norm. Add a `conditionNumber()` estimate method.

### `SpatialOperations` Static Methods with Unsafe Downcasts

- **Files:** `src/SpatialOperations.cpp:11-25`
- **Why fragile:** `crossProductMotion` and `crossProductForce` use `static_cast` from `const SpatialVector&` to `const MotionVector&` / `const ForceVector&`. If a user passes the wrong concrete type, the cast succeeds silently and produces garbage results. There is no runtime type checking. This bypasses the type safety that `MotionVector` and `ForceVector` are designed to provide.
- **Test coverage:** 11 tests. Tests pass correct types, so the issue is not caught.
- **Safe modification:** Change the function signatures to accept the correct concrete types directly (`const MotionVector&`, `const ForceVector&`), eliminating the need for casts.

## Scaling Limits

### No Branching Tree Support in ABA (Effectively)

- **Current capacity:** Single serial chain forward dynamics works (after fixing CR-02).
- **Limit:** Branching trees produce incorrect results due to CR-02 and the inward pass structure. The O(n²) child lookup in RNEA also limits practical tree sizes.
- **Scaling path:** Fix CR-02, then add explicit child-link index lists for O(n) branching tree traversal.

### No Joint Limit or Singularity Handling

- **Current capacity:** The singularity check in ABA (`std::abs(denom) < EPSILON` with `EPSILON = 1e-10`) detects exactly-zero inertia projections but does not handle near-singular configurations gracefully.
- **Limit:** At kinematic singularities, the ABA throws a `std::runtime_error`. There's no damping, regularization, or fallback behavior. For real robot control, this would cause abrupt failure.
- **Scaling path:** Add damped least-squares (DLS) / Levenberg-Marquardt regularization for near-singular configurations.

## Dependencies at Risk

### Eigen3 Version Compatibility (3.3+ Pinned)

- **Risk:** CMakeLists.txt pins `find_package(Eigen3 3.3 REQUIRED NO_MODULE)`. Eigen 5.x changed its CMake version range, causing `find_package` to fail even though the actual API used is compatible. Users on macOS Homebrew (which now ships Eigen 5.x) cannot build without manual workarounds.
- **Impact:** Build failure on systems with Eigen 5.x.
- **Migration plan:** Remove the version pin entirely, or add a version compatibility fallback. Test against both Eigen 3.4 and 5.x.

### Google Test Dependency

- **Risk:** CMakeLists.txt uses `find_package(GTest REQUIRED)` which expects GTest to be system-installed. There's no fallback to FetchContent or bundled GTest. This is fragile for CI or cross-platform builds.
- **Impact:** Build failure if GTest is not installed.
- **Migration plan:** Add a FetchContent fallback in CMakeLists.txt for GTest, or document the system dependency clearly.

## Missing Critical Features

### No Gravity Term in Dynamics Algorithms

- **Problem:** Both RNEA and ABA implementations explicitly state "Gravity forces (not implemented, add external forces)". There is no `gravityMode`, no gravity vector parameter, and no documentation on how to correctly add gravity via external forces. This is a significant omission for robotics users — gravity compensation is essential for all practical applications.
- **Files:**
  - `include/ForwardDynamics.h:21` (comment)
  - `include/InverseDynamics.h:21` (comment)
- **Blocks:** Any practical robotics simulation or control application.
- **Fix approach:** Add an optional gravity vector parameter to both `computeAccelerations` and `computeTorques`. The RNEA gravity acceleration propagates through the outward pass as an additional acceleration term in the base link.

### No Python Bindings or C++/Python Integration

- **Problem:** There's a standalone Python RNEA implementation (`robot_dynamics/rnea.py`) that is completely independent from the C++ library. It implements fewer features (no ABA, no Plücker transforms, no articulated body inertia, uses plain 3D vectors instead of spatial vectors). There is no Pybind11, Cython, or any C++/Python bridge.
- **Files:** `robot_dynamics/rnea.py` (does not import from C++)
- **Blocks:** Users who want to use the library from Python (common in robotics) cannot do so.
- **Fix approach:** Add pybind11 bindings for core classes, or at minimum document the gap.

### No Umbrella Header for Library Inclusion

- **Problem:** Users must include individual headers (`SpatialVector.h`, `PluckerTransform.h`, etc.) explicitly. There is no `include/SpatialAlgebra.h` that includes all public headers. This is a discoverability and convenience issue.
- **Files:** (missing) `include/SpatialAlgebra.h`
- **Fix approach:** Create a single umbrella header that includes all public headers, with documented include order.

## Test Coverage Gaps

### Untested: ABA Multi-Link Correctness

- **What's not tested:** Correct numerical acceleration values for multi-link serial chains and branching trees. The existing tests only verify "positive and finite" or single-link correctness.
- **Files:** `tests/TestForwardDynamics.cpp`
- **Risk:** The CR-02 bug exists undetected in production. Any user computing forward dynamics for a 2+ link chain gets wrong results.
- **Priority:** HIGH

### Untested: RNEA with Non-Zero Velocity

- **What's not tested:** All 5 RNEA tests use `qdot = 0.0`. The Coriolis/centrifugal computation is never exercised.
- **Files:** `tests/TestInverseDynamics.cpp`
- **Risk:** The Coriolis bias acceleration calculation in RNEA may be incorrect without detection.
- **Priority:** HIGH

### Untested: Force×Force Cross Product Correctness

- **What's not tested:** The cross product tests use pure torque or pure force inputs where the buggy formulas happen to produce correct results. Tests with mixed torque+force inputs (which trigger the inconsistencies) are absent.
- **Files:** `tests/TestSpatialVector.cpp`, `tests/TestSpatialUtils.cpp`
- **Risk:** Users relying on force cross product for any non-trivial wrench get silently incorrect results.
- **Priority:** HIGH

### Untested: Edge Cases (NaN, Inf, Extreme Values)

- **What's not tested:** NaN and Inf propagation through spatial vector operations. Extreme values (very large/small numbers), degenerate geometries (zero-length links, zero mass).
- **Files:** All test files.
- **Risk:** Numerical edge cases in downstream applications may produce silently incorrect physical results.
- **Priority:** MEDIUM

### Untested: Inertia Addition Edge Cases (Zero Combined Mass)

- **What's not tested:** `RigidBodyInertia::operator+` with zero combined mass (both bodies massless, or complementary masses that cancel). The code handles this case (`if (newMass > 0.0)`) but it's untested.
- **Files:** `include/RigidBodyInertia.h:64-68`
- **Risk:** Division by zero scenario exists but is guarded; the code path for `newMass == 0` returns uninitialized COM.
- **Priority:** LOW

### Untested: PluckerTransform Inertia Transform (Rotation Only)

- **What's not tested:** In `TestSpatialOperations.cpp`, the `createIdentityInertia()` helper returns a zero matrix instead of identity. The tests using this helper only verify `getMass()` and `getCom()`, not inertia matrix values. The actual inertia tensor transformation is never verified for numerical correctness.
- **Files:** `tests/TestSpatialOperations.cpp:27-37, 218-293`
- **Risk:** The inertia transformation via `SpatialOperations::transformInertia` / `PluckerTransform::tformRBI` may produce incorrect inertia matrices with undetected bugs.
- **Priority:** MEDIUM

---

*Concerns audit: 2026-05-17*
