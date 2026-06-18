# Phase 17: Benchmark Implementation - Pattern Map

**Mapped:** 2026-06-05
**Files analyzed:** 6 modified
**Analogs found:** 6 / 6

## File Classification

| Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---|---|---|---|---|
| `benchmarks/common/random_state.h` | utility | request-response | `benchmarks/common/random_state.h` (extend existing) | exact |
| `benchmarks/common/random_state.cpp` | utility | request-response | `benchmarks/common/random_state.cpp` (extend existing) | exact |
| `benchmarks/aba/bench_aba_stub.cpp` | test (benchmark) | CRUD | `benchmarks/aba/bench_aba_stub.cpp` (replace stub body) | exact |
| `benchmarks/rnea/bench_rnea_stub.cpp` | test (benchmark) | CRUD | `benchmarks/rnea/bench_rnea_stub.cpp` (replace stub body) | exact |
| `benchmarks/core/bench_core_stub.cpp` | test (benchmark) | transform | `benchmarks/core/bench_core_stub.cpp` (replace stubs with 20 functions) | exact |
| `benchmarks/bench_all.cpp` | test (benchmark entry) | request-response | `benchmarks/bench_all.cpp` (replace core registrations) | exact |

## Pattern Assignments

### `benchmarks/common/random_state.h` (utility, request-response)

**Analog:** Existing `benchmarks/common/random_state.h` (lines 1-60) extended with RESEARCH.md §Pattern 4 additions.

**Existing structure** (lines 16-60):
```cpp
#pragma once

#include <random>
#include <vector>

namespace SpatialAlgebra::Bench {

class RandomState {
public:
    RandomState() : rng_(42) {}

    std::vector<double> randomPositions(int nDOF);

    std::vector<double> randomVelocities(int nDOF);

private:
    std::mt19937 rng_;
};

}
```

**New methods to add** (after `randomVelocities` declaration, before `private:`):
```cpp
    // Phase 17: New methods
    std::vector<double> randomTorques(int nDOF);          // D-16: [-10, 10] Nm
    std::vector<double> randomAccelerations(int nDOF);    // D-17: [-5, 5] rad/s²

    // Zero-alloc fill variants (D-19)
    void fillPositions(double* buf, int n);
    void fillVelocities(double* buf, int n);
    void fillTorques(double* buf, int n);
    void fillAccelerations(double* buf, int n);
```

**Also add include:** Add `#include <Eigen/Dense>` after `#include <vector>` on line 17 (needed by research pattern, though not strictly used in random_state.h itself — included for downstream callers that use `Eigen::VectorXd` buffers).

---

### `benchmarks/common/random_state.cpp` (utility, request-response)

**Analog:** Existing `benchmarks/common/random_state.cpp` (lines 1-32) extended with RESEARCH.md §Pattern 4 additions.

**Existing implementation** (lines 9-32):
```cpp
#include "random_state.h"
#include <cmath>

namespace SpatialAlgebra::Bench {

std::vector<double> RandomState::randomPositions(int nDOF) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    std::vector<double> positions(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        positions[i] = dist(rng_);
    }
    return positions;
}

std::vector<double> RandomState::randomVelocities(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> velocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        velocities[i] = dist(rng_);
    }
    return velocities;
}

}
```

**New methods to add after `randomVelocities`** (before closing `}` of namespace):
```cpp
// ── Phase 17: New returning variants ──────────────────────────────
std::vector<double> RandomState::randomTorques(int nDOF) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);  // D-16
    std::vector<double> torques(nDOF);
    for (int i = 0; i < nDOF; ++i)
        torques[i] = dist(rng_);
    return torques;
}

std::vector<double> RandomState::randomAccelerations(int nDOF) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);  // D-17
    std::vector<double> accelerations(nDOF);
    for (int i = 0; i < nDOF; ++i)
        accelerations[i] = dist(rng_);
    return accelerations;
}

// ── Phase 17: Zero-alloc fill variants (D-19) ─────────────────────
void RandomState::fillPositions(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-M_PI_2, M_PI_2);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillVelocities(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillTorques(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}

void RandomState::fillAccelerations(double* buf, int n) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}
```

**Doxygen style:** Follow existing file pattern — add `@brief` / `@details` to each new method declaration in the header.

---

### `benchmarks/aba/bench_aba_stub.cpp` (test/benchmark, CRUD)

**Analog:** Existing `benchmarks/aba/bench_aba_stub.cpp` (lines 1-32) — replaces stub body with real benchmark logic.
**Solver API source:** `include/ForwardDynamics.h:156` — `computeAccelerations(const Eigen::VectorXd&, const Vector3d&)`.

**Existing file structure** (lines 1-32) — keep file header + includes + namespace, replace function body:
```cpp
/**
 * @file bench_aba_stub.cpp
 * @brief ABA forward dynamics benchmark stub
 * @details ... (keep existing header, update @brief to remove "stub")
 */
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // (void)state; (void)nDOF;  ← REMOVE these
```

**Replacement body** (from RESEARCH.md §Pattern 1, lines 236-266):
```cpp
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // ── Setup (outside timed loop) ──────────────────────────────────
    ModelFactory factory;
    RandomState rng;

    // Create n-DOF serial chain (D-01/D-02)
    ForwardDynamics fd = factory.createFD(nDOF);

    // Apply random joint state (D-09/D-10) — set once, reused
    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        fd.links[i].q = q[i];
        fd.links[i].qdot = qdot[i];
    }

    // Pre-allocated torque buffer (D-19: zero-alloc inside timed loop)
    Eigen::VectorXd tau(nDOF);

    // ── Timed loop ──────────────────────────────────────────────────
    for (auto _ : state) {
        // D-19: regenerated torques each iteration (zero-alloc via fill)
        rng.fillTorques(tau.data(), nDOF);
        // D-16: random torques; computeAccelerations is the measured operation
        fd.computeAccelerations(tau);
        // Prevent dead-code elimination
        benchmark::DoNotOptimize(fd.links[0].qddot);
    }
}
```

**Key solver API pattern** (`include/ForwardDynamics.h:156`):
```cpp
void computeAccelerations(const Eigen::VectorXd& tau,
                          const Vector3d& gravity = Vector3d::Zero());
```
- Takes `const Eigen::VectorXd& tau` — pre-allocate buffer, write through `tau.data()` via `fillTorques`.
- Also accepts optional `gravity` parameter — not used in benchmark (defaults to zero).
- Writes output to `fd.links[i].qddot` — use `DoNotOptimize(fd.links[0].qddot)` to consume.

**ModelFactory API** (`benchmarks/common/model_factory.h:76`):
```cpp
ForwardDynamics createFD(int nDOF, const JointConfig& cfg = JointConfig());
```
- Returns fully-populated `ForwardDynamics` with `nDOF` links.
- `JointConfig` default: Z revolute axis, half-meter X translation, 1.0 kg mass.
- Links initialized with zero state — benchmark sets `q` and `qdot` manually.

**`ForwardDynamics::Link` struct** (`include/ForwardDynamics.h:79-111`):
```cpp
struct Link {
    int parent;
    PluckerTransform X;
    RigidBodyInertia I;
    MotionVector S;
    double q, qdot, qddot;      // qddot is output from computeAccelerations
    MotionVector v, c;
    ForceVector f;
    ArticulatedBodyInertia Ia;
    ForceVector pa;
};
```
- Link state is accessed via `fd.links[i]` — benchmark writes `q`, `qdot`, reads `qddot`.
- `computeAccelerations` fully recomputes `v`, `c`, `Ia`, `pa`, `qddot` each call — no stale state carryover.

---

### `benchmarks/rnea/bench_rnea_stub.cpp` (test/benchmark, CRUD)

**Analog:** Existing `benchmarks/rnea/bench_rnea_stub.cpp` (lines 1-32) — replaces stub body with real benchmark logic.
**Solver API source:** `include/InverseDynamics.h:144` — `computeTorques(const Eigen::VectorXd&, const Vector3d&)`.

**Existing file structure** (lines 1-32) — keep file header + includes + namespace, replace function body:
```cpp
/**
 * @file bench_rnea_stub.cpp
 * @brief RNEA inverse dynamics benchmark stub
 * @details ... (keep existing header, update @brief to remove "stub")
 */
#include <benchmark/benchmark.h>
#include "InverseDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF) {
    // (void)state; (void)nDOF;  ← REMOVE these
```

**Replacement body** (from RESEARCH.md §Pattern 2, lines 291-323):
```cpp
void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF) {
    // ── Setup (outside timed loop) ──────────────────────────────────
    ModelFactory factory;
    RandomState rng;

    // Create n-DOF serial chain (D-01/D-02)
    InverseDynamics id = factory.createID(nDOF);

    // Apply random joint state (D-09/D-10) — set once, reused
    std::vector<double> q = rng.randomPositions(nDOF);
    std::vector<double> qdot = rng.randomVelocities(nDOF);
    for (int i = 0; i < nDOF; ++i) {
        id.links[i].q = q[i];
        id.links[i].qdot = qdot[i];
    }

    // Pre-allocated acceleration buffer (D-19: zero-alloc inside timed loop)
    // Range: [-5, 5] rad/s² (D-17)
    Eigen::VectorXd qddot(nDOF);

    // ── Timed loop ──────────────────────────────────────────────────
    for (auto _ : state) {
        // D-19: regenerated accelerations each iteration
        rng.fillAccelerations(qddot.data(), nDOF);

        // D-17: random qddot; computeTorques is the measured operation
        // Returns Eigen::VectorXd (heap-allocated each call — inherent to API)
        Eigen::VectorXd tau = id.computeTorques(qddot);

        // Prevent dead-code elimination
        benchmark::DoNotOptimize(tau[0]);
    }
}
```

**Key solver API pattern** (`include/InverseDynamics.h:144`):
```cpp
Eigen::VectorXd computeTorques(const Eigen::VectorXd& qddot,
                                const Vector3d& gravity = Vector3d::Zero());
```
- Takes `const Eigen::VectorXd& qddot` — pre-allocate buffer, write through `qddot.data()` via `fillAccelerations`.
- **Returns** a newly heap-allocated `Eigen::VectorXd` each call — this allocation IS inside the timed loop and IS part of the measured time. Do not try to eliminate it.
- Accepts optional `gravity` parameter — not used in benchmark (defaults to zero).

**ModelFactory API** (`benchmarks/common/model_factory.h:84`):
```cpp
InverseDynamics createID(int nDOF, const JointConfig& cfg = JointConfig());
```
- Returns fully-populated `InverseDynamics` with `nDOF` `InverseDynamicsLink` objects.

**`InverseDynamicsLink` struct** (`include/InverseDynamics.h:75-100`):
```cpp
struct InverseDynamicsLink {
    int parent;
    PluckerTransform X;
    RigidBodyInertia I;
    MotionVector S;
    double q, qdot, qddot;      // qddot is INPUT to computeTorques
    MotionVector v, a;
};
```
- Note: `qddot` is an **input** field in RNEA (vs output in ABA). The `computeTorques` call copies `qddot[i]` from the input vector into `links[i].qddot` internally.

---

### `benchmarks/core/bench_core_stub.cpp` (test/benchmark, transform)

**Analog:** Existing `benchmarks/core/bench_core_stub.cpp` (lines 1-43) — replaces 2 placeholder stubs with 20 individual benchmark functions.
**API sources:** `include/PluckerTransform.h:77-216`, `include/SpatialUtils.h:78-126`, `include/RigidBodyInertia.h:27-122`, `include/ArticulatedBodyInertia.h:77-220`.

**Existing includes** (lines 21-26) — these chain-include all needed headers via PluckerTransform.h, but explicit includes are safer:
```cpp
#include <benchmark/benchmark.h>
#include "PluckerTransform.h"
#include "SpatialUtils.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
```
**Updated includes** (remove unused `model_factory.h`/`random_state.h`, add inertia headers):
```cpp
#include <benchmark/benchmark.h>
#include "PluckerTransform.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "SpatialUtils.h"
#include "RigidBodyInertia.h"
#include "ArticulatedBodyInertia.h"
#include "LowerTriangular.h"

using namespace SpatialAlgebra;
```

**Test object factory pattern** — static helper functions creating pre-constructed objects (from RESEARCH.md §Pattern 3, lines 352-378):
```cpp
static PluckerTransform make_test_transform() {
    return PluckerTransform(
        Rotation(Eigen::AngleAxisd(0.5, Vector3d::UnitZ())),
        Vector3d(0.3, 0.1, 0.7)
    );
}

static MotionVector make_test_mv() {
    return MotionVector(Vector3d(0.2, 0.5, 0.8), Vector3d(0.1, 0.3, 0.6));
}

static ForceVector make_test_fv() {
    return ForceVector(Vector3d(0.4, 0.9, 0.2), Vector3d(0.7, 0.1, 0.5));
}

static RigidBodyInertia make_test_rbi() {
    return RigidBodyInertia(1.5, Vector3d(0.1, 0.2, 0.3),
                            LowerTriangular::Identity(3));
}

static ArticulatedBodyInertia make_test_abi() {
    return ArticulatedBodyInertia(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );
}
```

**Plücker transform benchmark pattern** (example — repeat for all 8):
```cpp
void BM_TransformMotion(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    MotionVector mv_in = make_test_mv();

    for (auto _ : state) {
        MotionVector mv_out = X.transformMotion(mv_in);
        benchmark::DoNotOptimize(mv_out);
    }
}
```

The 8 Plücker transform benchmarks. The only thing that changes is the method call and return type:
| Function | Method | Return Type |
|---|---|---|
| `BM_TransformMotion` | `X.transformMotion(mv_in)` | `MotionVector` |
| `BM_TransformForce` | `X.transformForce(fv_in)` | `ForceVector` |
| `BM_InverseTransformMotion` | `X.inverseTransformMotion(mv_in)` | `MotionVector` |
| `BM_InverseTransformForce` | `X.inverseTransformForce(fv_in)` | `ForceVector` |
| `BM_TformRBI` | `X.tformRBI(rbi_in)` | `RigidBodyInertia` |
| `BM_InvtformRBI` | `X.invtformRBI(rbi_in)` | `RigidBodyInertia` |
| `BM_TformABI` | `X.tformABI(abi_in)` | `ArticulatedBodyInertia` |
| `BM_InvtformABI` | `X.invtformABI(abi_in)` | `ArticulatedBodyInertia` |

**Cross product benchmark pattern** (from RESEARCH.md §Pattern 3, lines 462-500):
```cpp
void BM_CrossMvMv(benchmark::State& state) {
    MotionVector mv1 = make_test_mv();
    MotionVector mv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));

    for (auto _ : state) {
        MotionVector result = cross(mv1, mv2);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossMvFv(benchmark::State& state) {
    MotionVector mv1 = make_test_mv();
    ForceVector fv1 = make_test_fv();

    for (auto _ : state) {
        ForceVector result = cross(mv1, fv1);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossFvMv(benchmark::State& state) {
    ForceVector fv1 = make_test_fv();
    MotionVector mv1 = make_test_mv();

    for (auto _ : state) {
        ForceVector result = cross(fv1, mv1);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossFvFv(benchmark::State& state) {
    ForceVector fv1 = make_test_fv();
    ForceVector fv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));

    for (auto _ : state) {
        ForceVector result = cross(fv1, fv2);
        benchmark::DoNotOptimize(result);
    }
}
```

**Inertia operation benchmark pattern** (from RESEARCH.md §Pattern 3, lines 503-575):
RBI operations (3 benchmarks):
```cpp
void BM_RBIApply(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        ForceVector fv = rbi.apply(mv);
        benchmark::DoNotOptimize(fv);
    }
}

void BM_RBIAdd(benchmark::State& state) {
    RigidBodyInertia rbi1 = make_test_rbi();
    RigidBodyInertia rbi2(2.0, Vector3d(0.05, 0.1, 0.15),
                          LowerTriangular::Identity(3));
    for (auto _ : state) {
        RigidBodyInertia result = rbi1 + rbi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_RBIScale(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        RigidBodyInertia result = rbi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}
```

ABI operations (4 benchmarks):
```cpp
void BM_ABIApply(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        ForceVector fv = abi.apply(mv);
        benchmark::DoNotOptimize(fv);
    }
}

void BM_ABIAdd(benchmark::State& state) {
    ArticulatedBodyInertia abi1 = make_test_abi();
    ArticulatedBodyInertia abi2(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi1 + abi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIAddRBI(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi + rbi;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIScale(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}
```

**Note on signature:** All 20 core benchmark functions use `(benchmark::State& state)` — no `int nDOF` parameter (the O(1) operations have no DOF dependency).

---

### `benchmarks/bench_all.cpp` (test/benchmark entry, request-response)

**Analog:** Existing `benchmarks/bench_all.cpp` (lines 1-52) — replace 2 core placeholder forward declarations + registrations with 20 individual ones.

**Existing forward declarations** (lines 20-24) — replace:
```cpp
// Old (lines 20-24):
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF);
void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF);
void BM_PluckerTransform(benchmark::State& state, int nDOF);
void BM_CrossProduct(benchmark::State& state, int nDOF);

// New:
void BM_ABA_ForwardDynamics(benchmark::State&, int);
void BM_RNEA_InverseDynamics(benchmark::State&, int);

void BM_TransformMotion(benchmark::State&);
void BM_TransformForce(benchmark::State&);
void BM_InverseTransformMotion(benchmark::State&);
void BM_InverseTransformForce(benchmark::State&);
void BM_TformRBI(benchmark::State&);
void BM_InvtformRBI(benchmark::State&);
void BM_TformABI(benchmark::State&);
void BM_InvtformABI(benchmark::State&);
void BM_CrossMvMv(benchmark::State&);
void BM_CrossMvFv(benchmark::State&);
void BM_CrossFvMv(benchmark::State&);
void BM_CrossFvFv(benchmark::State&);
void BM_RBIApply(benchmark::State&);
void BM_RBIAdd(benchmark::State&);
void BM_RBIScale(benchmark::State&);
void BM_ABIApply(benchmark::State&);
void BM_ABIAdd(benchmark::State&);
void BM_ABIAddRBI(benchmark::State&);
void BM_ABIScale(benchmark::State&);
```

**Existing DOF sweep loop** (lines 29-46) — keep ABA/RNEA, remove old core stubs, add individual core registrations after the loop:
```cpp
int main(int argc, char** argv) {
    // ABA/RNEA: DOF sweep n=1..20
    for (int n = 1; n <= 20; ++n) {
        benchmark::RegisterBenchmark(
            ("BM_ABA_ForwardDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_ABA_ForwardDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_RNEA_InverseDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_RNEA_InverseDynamics, n
        );
    }

    // Core microbenchmarks: 20 individual operations (no DOF sweep)
    benchmark::RegisterBenchmark("BM_TransformMotion", BM_TransformMotion);
    benchmark::RegisterBenchmark("BM_TransformForce", BM_TransformForce);
    benchmark::RegisterBenchmark("BM_InverseTransformMotion", BM_InverseTransformMotion);
    benchmark::RegisterBenchmark("BM_InverseTransformForce", BM_InverseTransformForce);
    benchmark::RegisterBenchmark("BM_TformRBI", BM_TformRBI);
    benchmark::RegisterBenchmark("BM_InvtformRBI", BM_InvtformRBI);
    benchmark::RegisterBenchmark("BM_TformABI", BM_TformABI);
    benchmark::RegisterBenchmark("BM_InvtformABI", BM_InvtformABI);
    benchmark::RegisterBenchmark("BM_CrossMvMv", BM_CrossMvMv);
    benchmark::RegisterBenchmark("BM_CrossMvFv", BM_CrossMvFv);
    benchmark::RegisterBenchmark("BM_CrossFvMv", BM_CrossFvMv);
    benchmark::RegisterBenchmark("BM_CrossFvFv", BM_CrossFvFv);
    benchmark::RegisterBenchmark("BM_RBIApply", BM_RBIApply);
    benchmark::RegisterBenchmark("BM_RBIAdd", BM_RBIAdd);
    benchmark::RegisterBenchmark("BM_RBIScale", BM_RBIScale);
    benchmark::RegisterBenchmark("BM_ABIApply", BM_ABIApply);
    benchmark::RegisterBenchmark("BM_ABIAdd", BM_ABIAdd);
    benchmark::RegisterBenchmark("BM_ABIAddRBI", BM_ABIAddRBI);
    benchmark::RegisterBenchmark("BM_ABIScale", BM_ABIScale);

    // Unchanged boilerplate
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
```

**Key difference:** Core benchmarks are registered WITHOUT `n` argument and WITHOUT the `"/{N}DOF"` name suffix, because the functions take `(benchmark::State&)` only (no DOF parameter).

## Shared Patterns

### Per-Iteration State Regeneration (Zero-Alloc Fill)
**Source:** `benchmarks/aba/bench_aba_stub.cpp` (ABA pattern), `benchmarks/rnea/bench_rnea_stub.cpp` (RNEA pattern)
**Apply to:** All solver benchmarks with regenerated inputs
```cpp
// Pre-allocate buffer outside timed loop:
Eigen::VectorXd tau(nDOF);

for (auto _ : state) {
    // Fill through Eigen's contiguous storage (zero-alloc):
    rng.fillTorques(tau.data(), nDOF);
    // Use the buffer in the measured operation:
    fd.computeAccelerations(tau);
    // Consume output to prevent dead-code elimination:
    benchmark::DoNotOptimize(fd.links[0].qddot);
}
```

### Solver State Application Pattern
**Source:** RESEARCH.md §Common Operation 4
**Apply to:** ABA and RNEA benchmark setup
```cpp
// After factory.createFD(nDOF) / createID(nDOF):
std::vector<double> q = rng.randomPositions(nDOF);
std::vector<double> qdot = rng.randomVelocities(nDOF);
for (int i = 0; i < nDOF; ++i) {
    fd.links[i].q = q[i];
    fd.links[i].qdot = qdot[i];
}
```

### DoNotOptimize Placement
**Source:** Google Benchmark best practices (verified for SpatialAlgebra types)
**Apply to:** All benchmark loops
- For `SpatialVector`/`MotionVector`/`ForceVector`: `benchmark::DoNotOptimize(result)` on the whole object works because Eigen types write to memory.
- For scalar output (`qddot`, `tau[i]`): `benchmark::DoNotOptimize(result)` directly.
- For `PluckerTransform`/`RigidBodyInertia`/`ArticulatedBodyInertia`: `DoNotOptimize(result)` on the returned object.

### Google Benchmark Registration Pattern
**Source:** `benchmarks/bench_all.cpp` existing pattern
**Apply to:** All benchmark registrations in `bench_all.cpp`

With DOF parameter:
```cpp
benchmark::RegisterBenchmark("BM_Name/{N}DOF", BM_Function, nDOF);
```
Without DOF parameter (core microbenchmarks):
```cpp
benchmark::RegisterBenchmark("BM_Name", BM_Function);
// Note: no n argument passed, no DOF suffix in name
```

### Doxygen Style for Benchmark Files
**Source:** Existing stub files (e.g., `bench_aba_stub.cpp:1-14`)
**Apply to:** All modified benchmark files
```cpp
/**
 * @file bench_aba.cpp
 * @brief ABA forward dynamics benchmark
 * @details Benchmarks ForwardDynamics::computeAccelerations(tau)
 *          with a DOF sweep n=1..20. ...
 */
```

### RandomState Fill Variant Pattern
**Source:** RESEARCH.md §Pattern 4
**Apply to:** All fill methods in `random_state.cpp`
```cpp
void RandomState::fillTorques(double* buf, int n) {
    // Create distribution each call (very cheap — two doubles)
    // to guarantee identical sequence to returning variant
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < n; ++i)
        buf[i] = dist(rng_);
}
```
**Key constraint:** Each `fill*` method must create a new `std::uniform_real_distribution` each call (not reuse a stored one) to produce the same sequence as the corresponding returning method for the same `n`.

## No Analog Found

All 6 files have exact matches in the existing codebase (each file extends/replaces its own stub). No files require RESEARCH.md patterns as the sole reference — RESEARCH.md provides the concrete implementation pattern but the file structure, include style, and namespace conventions come from the existing stubs.

## Metadata

**Analog search scope:** `benchmarks/*.cpp`, `benchmarks/common/*.h/.cpp`, `include/ForwardDynamics.h`, `include/InverseDynamics.h`, `include/PluckerTransform.h`, `include/RigidBodyInertia.h`, `include/ArticulatedBodyInertia.h`, `include/SpatialUtils.h`, `include/LowerTriangular.h`, `include/MotionVector.h`, `include/ForceVector.h`, `include/Rotation.h`
**Files scanned:** ~18
**Pattern extraction date:** 2026-06-05
