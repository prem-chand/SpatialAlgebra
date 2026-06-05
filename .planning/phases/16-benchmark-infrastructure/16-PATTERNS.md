# Phase 16: Benchmark Infrastructure - Pattern Map

**Mapped:** 2026-06-05
**Files analyzed:** 10 new/modified
**Analogs found:** 10 / 10

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `CMakeLists.txt` (modify) | config | build-time | `CMakeLists.txt:42-46` (ENABLE_COVERAGE option) | exact |
| `benchmarks/CMakeLists.txt` | config | build-time | `CMakeLists.txt:17-23` (FetchContent GTest) | exact |
| `benchmarks/common/model_factory.h` | utility | CRUD | `include/SpatialOperations.h` (static utility class) | role-match |
| `benchmarks/common/model_factory.cpp` | utility | CRUD | `src/SpatialOperations.cpp` (utility impl) | role-match |
| `benchmarks/common/random_state.h` | utility | — | `include/LowerTriangular.h` (class with priv data + methods) | partial |
| `benchmarks/common/random_state.cpp` | utility | — | `src/SpatialOperations.cpp` (utility impl) | role-match |
| `benchmarks/aba/bench_aba_stub.cpp` | test | request-response | `tests/TestPluckerTransform.cpp` (GTest file) | exact |
| `benchmarks/rnea/bench_rnea_stub.cpp` | test | request-response | `tests/TestPluckerTransform.cpp` (GTest file) | exact |
| `benchmarks/core/bench_core_stub.cpp` | test | request-response | `tests/TestPluckerTransform.cpp` (GTest file) | exact |
| `benchmarks/bench_all.cpp` | test | request-response | `tests/TestSpatialVector.cpp:605-609` (main()) | role-match |

## Pattern Assignments

### `CMakeLists.txt` (root, modify) — config, build-time

**Analog:** `CMakeLists.txt:42-46` (option + guard pattern)

**Option + guard pattern** (lines 1-4 of the modification):
```cmake
# Analog: lines 42-46 of CMakeLists.txt — ENABLE_COVERAGE pattern
# option(ENABLE_COVERAGE "Enable coverage flags for CI" OFF)
# if(ENABLE_COVERAGE)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
# endif()

# New code to add (after line 46, before "Build tests" line 48):
option(SA_BUILD_BENCHMARKS "Build performance benchmarks" OFF)
if(SA_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
```

**Insertion point** (line 48 marker):
- Insert after the `ENABLE_COVERAGE` block (after line 46 `endif()`), before `# Build tests` on line 48.
- The project also has `add_subdirectory(examples)` at line 157 — `add_subdirectory(benchmarks)` follows the same pattern.

---

### `benchmarks/CMakeLists.txt` — config, build-time

**Analog:** `CMakeLists.txt:17-23` (FetchContent pattern for GTest) + `examples/CMakeLists.txt` (subdirectory CMakeLists.txt pattern)

**FetchContent pattern** (lines 17-23):
```cmake
# Analog: CMakeLists.txt:17-23 — GTest FetchContent pattern
# find_package(GTest QUIET)
# if(NOT GTest_FOUND)
#     include(FetchContent)
#     FetchContent_Declare(
#         googletest
#         URL https://github.com/google/googletest/archive/release-1.12.1.zip
#     )
#     FetchContent_MakeAvailable(googletest)
#     include_directories(${googletest_SOURCE_DIR}/googletest/include)
# endif()

# Benchmarks CMakeLists.txt should use:
cmake_minimum_required(VERSION 3.19)
project(SpatialAlgebraBenchmarks LANGUAGES CXX)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
)
FetchContent_MakeAvailable(googlebenchmark)

set(CMAKE_CXX_STANDARD 17)
```

**Subdirectory CMakeLists.txt pattern** (from `examples/CMakeLists.txt` lines 1-3):
```cmake
# Analog: examples/CMakeLists.txt:1-3
# # Examples CMakeLists.txt
# # Build usage examples for SpatialAlgebra library
# set(CMAKE_CXX_STANDARD 17)
```

**Executable + link pattern** (from `examples/CMakeLists.txt:10-11`):
```cmake
# Analog: examples/CMakeLists.txt:10-11 — executable + link
# add_executable(example_vectors basic_vectors.cpp)
# target_link_libraries(example_vectors SpatialAlgebra Eigen3::Eigen)
```

**LTO pattern** (from RESEARCH.md, no direct analog — compiler-specific):
```cmake
# D-14: LTO enabled for benchmark targets
target_compile_options(bench_all PRIVATE -O3 -DNDEBUG -flto)
target_link_options(bench_all PRIVATE -flto)
```

**`bench_all` executable linking pattern** (references `benchmark::benchmark`, Eigen3, and the library):
```cmake
add_executable(bench_all
    bench_all.cpp
    ${ABA_SOURCES}
    ${RNEA_SOURCES}
    ${CORE_SOURCES}
)

target_link_libraries(bench_all
    PRIVATE
    SpatialAlgebra
    benchmark::benchmark
    Eigen3::Eigen
    bench_common
)
```

---

### `benchmarks/common/model_factory.h` — utility, CRUD

**Analog:** `include/SpatialOperations.h` (static utility class pattern, lines 1-33)

**Imports + namespace pattern** (SpatialOperations.h:1-16):
```cpp
// Analog: SpatialOperations.h:1-16 — utility class with namespace
#ifndef SPATIAL_OPERATIONS_H
#define SPATIAL_OPERATIONS_H

#include "MotionVector.h"
#include "ForceVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include <Eigen/Geometry>

namespace SpatialAlgebra {
// ... class declaration ...
}  // namespace SpatialAlgebra

#endif
```

**For model_factory.h, use `#pragma once`** (LowerTriangular.h exception at line 1):
```cpp
// Following LowerTriangular.h's #pragma once pattern (belongs to the project conventions):
#pragma once

#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "MotionVector.h"
#include <Eigen/Dense>
#include <vector>

namespace SpatialAlgebra::Bench {

struct JointConfig {
    MotionVector axis;      // Screw axis (default: Z revolute)
    Eigen::Vector3d translation;  // Transform from parent
    double mass = 1.0;
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
};

class ModelFactory {
public:
    ForwardDynamics createFD(int nDOF, const JointConfig& cfg = JointConfig());
    InverseDynamics createID(int nDOF, const JointConfig& cfg = JointConfig());
    // ... branching overloads ...
};

}  // namespace SpatialAlgebra::Bench
```

**ForwardDynamicsLink struct pattern** (ForwardDynamics.h:79-111 — for matching factory output):
```cpp
// Analog: ForwardDynamics.h:79-111 — Link struct that factory must produce
struct Link
{
    int parent;                     ///< Parent link index (-1 for base)
    PluckerTransform X;             ///< Transform from parent to this link
    RigidBodyInertia I;             ///< Rigid body inertia
    MotionVector S;                 ///< Joint motion axis (screw axis)
    double q;                       ///< Joint position
    double qdot;                    ///< Joint velocity
    double qddot;                   ///< Joint acceleration (output)
    MotionVector v;                 ///< Spatial velocity
    MotionVector c;                 ///< Bias acceleration
    ForceVector f;                  ///< Spatial force (external forces)
    ArticulatedBodyInertia Ia;      ///< Articulated body inertia
    ForceVector pa;                 ///< Bias force

    Link() : parent(-1),
             X(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero()),
             I(1.0, Vector3d::Zero(), lt::Identity(3)),
             S(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             q(0.0), qdot(0.0), qddot(0.0),
             v(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             c(MotionVector(Vector3d::Zero(), Vector3d::Zero())),
             f(ForceVector(Vector3d::Zero(), Vector3d::Zero())),
             Ia(lt::Identity(3), Eigen::Matrix3d::Zero(), lt::Identity(3)),
             pa(ForceVector(Vector3d::Zero(), Vector3d::Zero())) {}
};
```

---

### `benchmarks/common/model_factory.cpp` — utility, CRUD

**Analog:** `src/SpatialOperations.cpp` (utility implementation, lines 1-26)

**Implementation pattern** (SpatialOperations.cpp:1-26):
```cpp
// Analog: SpatialOperations.cpp:1-26 — utility class implementation
/**
 * @file SpatialOperations.cpp
 * @brief Implementation of static utility class for spatial algebra operations
 */

#include "SpatialOperations.h"
#include "SpatialUtils.h"

namespace SpatialAlgebra {

SpatialVector SpatialOperations::crossProductMotion(const MotionVector& v1, const MotionVector& v2) {
    return cross(v1, v2);
}
// ... more methods ...

} // namespace SpatialAlgebra
```

**For model_factory.cpp, the link construction pattern** (from examples/dynamics.cpp:43-70):
```cpp
// Analog: examples/dynamics.cpp:43-70 — Link construction pattern
// Link 1: Base link (connected to world)
Link link1;
link1.parent = -1;  // Base link (no parent)

// Set inertia: mass=1kg, COM at origin, diagonal inertia
lt I1_tensor(3);
I1_tensor(0, 0) = 0.1; I1_tensor(1, 0) = 0.0; I1_tensor(1, 1) = 0.1;
I1_tensor(2, 0) = 0.0; I1_tensor(2, 1) = 0.0; I1_tensor(2, 2) = 0.1;
link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), I1_tensor);

// Transform from parent: identity
Rotation R1;
R1.setIdentity();
link1.X = PluckerTransform(R1, Vector3d::Zero());

// Joint axis: rotation about Z axis (revolute joint)
link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());

link1.q = 0.0;
link1.qdot = 0.0;
```

**Factory implementation pattern** (from examples/dynamics.cpp, model_factory creates multiple links in a loop):
```cpp
// ModelFactory::createFD implementation follows the link construction pattern
// but in a loop with configurable parameters:
ForwardDynamics ModelFactory::createFD(int nDOF, const JointConfig& cfg) {
    ForwardDynamics fd;
    fd.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        ForwardDynamics::Link link;
        link.parent = i - 1;  // -1 for base, i-1 for children

        // Translation along X (from cfg)
        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        // Inertia from cfg
        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));

        // Joint axis from cfg
        link.S = cfg.axis;

        link.q = 0.0;
        link.qdot = 0.0;
        fd.links.push_back(std::move(link));
    }

    return fd;
}
```

---

### `benchmarks/common/random_state.h` — utility, helper

**Analog:** `include/LowerTriangular.h` (class with private data + public methods, lines 75-125)

**Class pattern** (LowerTriangular.h:75-120):
```cpp
// Analog: LowerTriangular.h:75-120 — class with private data + public interface
class LowerTriangular
{
private:
    Eigen::VectorXd data;
    int n;

public:
    explicit LowerTriangular(int size) : n(size)
    {
        data.resize(size * (size + 1) / 2);
        data.setZero();
    }
    // ... methods ...
};
```

**For random_state.h:**
```cpp
#pragma once
#include <random>
#include <vector>

namespace SpatialAlgebra::Bench {

class RandomState {
public:
    RandomState() : rng_(42) {}  // D-11: Fixed seed

    std::vector<double> randomPositions(int nDOF);
    std::vector<double> randomVelocities(int nDOF);

private:
    std::mt19937 rng_;
};

}
```

---

### `benchmarks/common/random_state.cpp` — utility, helper

**Analog:** `src/SpatialOperations.cpp` (utility implementation pattern)

**Implementation pattern** (SpatialOperations.cpp:1-26 — adapted):
```cpp
// Pattern: utility implementation using standard C++ <random>
// D-09: positions uniform in [-π/2, π/2]
// D-10: velocities uniform in [-5, 5]
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

---

### `benchmarks/aba/bench_aba_stub.cpp` — test (benchmark), request-response

**Analog:** `tests/TestSpatialVector.cpp` (GTest includes pattern, lines 1-7)

**Includes pattern** (TestSpatialVector.cpp:1-7):
```cpp
// Analog: TestSpatialVector.cpp:1-7 — includes pattern for test/benchmark files
#include "SpatialVector.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include <gtest/gtest.h>

using namespace SpatialAlgebra;
```

**Benchmark stub pattern** (from RESEARCH.md Pattern 6):
```cpp
// benchmarks/aba/bench_aba_stub.cpp (STUB — Phase 17 fills implementation)
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement ABA benchmark with DOF sweep
    (void)state;
    (void)nDOF;
}
```

---

### `benchmarks/rnea/bench_rnea_stub.cpp` — test (benchmark), request-response

**Analog:** Same includes pattern as `bench_aba_stub.cpp`

```cpp
// benchmarks/rnea/bench_rnea_stub.cpp (STUB — Phase 17 fills implementation)
#include <benchmark/benchmark.h>
#include "InverseDynamics.h"
#include "common/model_factory.h"
#include "common/random_state.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement RNEA benchmark with DOF sweep
    (void)state;
    (void)nDOF;
}
```

---

### `benchmarks/core/bench_core_stub.cpp` — test (benchmark), request-response

**Analog:** Same includes pattern but for core ops:

```cpp
// benchmarks/core/bench_core_stub.cpp (STUB — Phase 17 fills implementation)
#include <benchmark/benchmark.h>
#include "PluckerTransform.h"
#include "SpatialUtils.h"

using namespace SpatialAlgebra;

void BM_PluckerTransform(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement Plücker transform benchmark
    (void)state;
    (void)nDOF;
}

void BM_CrossProduct(benchmark::State& state, int nDOF) {
    // TODO: Phase 17 — implement cross product benchmark
    (void)state;
    (void)nDOF;
}
```

---

### `benchmarks/bench_all.cpp` — test (benchmark entry point), request-response

**Analog:** `tests/TestSpatialVector.cpp:605-609` (main function pattern)

**Test main() pattern** (TestSpatialVector.cpp:605-609):
```cpp
// Analog: TestSpatialVector.cpp:605-609 — main() pattern
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

**Benchmark main() pattern** (from RESEARCH.md Pattern 2):
```cpp
// Analog adapts GTest main → Google Benchmark custom main with RegisterBenchmark
#include <benchmark/benchmark.h>
#include "common/model_factory.h"
#include "common/random_state.h"

// Forward declarations for benchmark functions (defined in domain stubs)
void BM_ABA_ForwardDynamics(benchmark::State& state, int nDOF);
void BM_RNEA_InverseDynamics(benchmark::State& state, int nDOF);
void BM_PluckerTransform(benchmark::State& state, int nDOF);
void BM_CrossProduct(benchmark::State& state, int nDOF);

int main(int argc, char** argv) {
    // Register ABA benchmarks for DOF sweep n=1..20
    for (int n = 1; n <= 20; ++n) {
        benchmark::RegisterBenchmark(
            ("BM_ABA_ForwardDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_ABA_ForwardDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_RNEA_InverseDynamics/" + std::to_string(n) + "DOF").c_str(),
            BM_RNEA_InverseDynamics, n
        );
        benchmark::RegisterBenchmark(
            ("BM_PluckerTransform/" + std::to_string(n) + "DOF").c_str(),
            BM_PluckerTransform, n
        );
        benchmark::RegisterBenchmark(
            ("BM_CrossProduct/" + std::to_string(n) + "DOF").c_str(),
            BM_CrossProduct, n
        );
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
```

---

## Shared Patterns

### FetchContent for External Dependencies
**Source:** `CMakeLists.txt:17-23`
**Apply to:** `benchmarks/CMakeLists.txt`
```cmake
# Pattern for adding external dependencies via FetchContent
include(FetchContent)
FetchContent_Declare(
    <target_name>
    GIT_REPOSITORY <url>     # or URL <zip>
    GIT_TAG <version>
)
FetchContent_MakeAvailable(<target_name>)
```
Key adaptation for benchmark: Must set `BENCHMARK_ENABLE_TESTING OFF` and `BENCHMARK_ENABLE_INSTALL OFF` before `FetchContent_MakeAvailable` to avoid pulling in GTest as a transitive dependency.

### CMake `option()` Guard Pattern
**Source:** `CMakeLists.txt:42-46`
**Apply to:** Root `CMakeLists.txt` modification
```cmake
option(SA_BUILD_BENCHMARKS "Build performance benchmarks" OFF)
if(SA_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
```
Insertion point: after the ENABLE_COVERAGE block (line 46), before `# Build tests` (line 48).

### Executable + Link Libraries Pattern
**Source:** `examples/CMakeLists.txt:10-11`, `CMakeLists.txt:59-63`
**Apply to:** `benchmarks/CMakeLists.txt` for `bench_all`
```cmake
# Pattern: add_executable + target_link_libraries
add_executable(<target> <sources>)
target_link_libraries(<target>
    PRIVATE
    SpatialAlgebra
    benchmark::benchmark
    Eigen3::Eigen
    bench_common
)
```

### Includes + Using Namespace Pattern
**Source:** All test files (e.g., `TestSpatialVector.cpp:1-13`, `TestPluckerTransform.cpp:1-14`)
**Apply to:** All benchmark stub files
```cpp
#include <benchmark/benchmark.h>
#include "ForwardDynamics.h"
// ... other SpatialAlgebra headers ...

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Bench;

// Benchmark functions declared with:
void BM_*(benchmark::State& state, int nDOF);
```

### Link Struct Construction Pattern
**Source:** `examples/dynamics.cpp:43-70`, `ForwardDynamics.h:79-111`
**Apply to:** `benchmarks/common/model_factory.cpp`
```cpp
Link link;
link.parent = -1;  // Base link
link.X = PluckerTransform(Rotation::Identity(), Vector3d::Zero());
link.I = RigidBodyInertia(mass, com, lt::Identity(3));
link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());  // Z revolute
link.q = 0.0;
link.qdot = 0.0;
fd.links.push_back(link);
```

### Main Function Pattern
**Source:** `tests/TestSpatialVector.cpp:605-609`, `tests/TestPluckerTransform.cpp:936-940`
**Apply to:** `benchmarks/bench_all.cpp`

GTest pattern:
```cpp
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

Google Benchmark adaptation:
```cpp
int main(int argc, char** argv) {
    benchmark::RegisterBenchmark("BM_*", BM_Fn, args...);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
```

### Low-Level Include + Namespace Aliases Pattern
**Source:** Many files use `using lt = LowerTriangular` — e.g., `src/main.cpp:8`, `examples/dynamics.cpp:28`
**Apply to:** `benchmarks/common/model_factory.h/.cpp`
```cpp
#include "SpatialVector.h"
#include "LowerTriangular.h"

namespace SpatialAlgebra {
using lt = LowerTriangular;
using Vector3d = Eigen::Vector3d;
// Type alias pattern: using mv = MotionVector, using fv = ForceVector
}  // namespace SpatialAlgebra
```

## Shared Utilities Placement

The shared utilities (`model_factory.h/.cpp`, `random_state.h/.cpp`) live in `benchmarks/common/` and are built as a static library `bench_common` by a `benchmarks/common/CMakeLists.txt`:

```cmake
# benchmarks/common/CMakeLists.txt
add_library(bench_common STATIC
    model_factory.cpp
    random_state.cpp
)
target_include_directories(bench_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(bench_common PUBLIC SpatialAlgebra Eigen3::Eigen)
```

## No Analog Found

All files have close matches in the existing codebase. No files require RESEARCH.md patterns as the sole reference.

## Metadata

**Analog search scope:** `CMakeLists.txt`, `examples/CMakeLists.txt`, `include/*.h`, `src/*.cpp`, `tests/*.cpp`, `examples/*.cpp`
**Files scanned:** ~20
**Pattern extraction date:** 2026-06-05
