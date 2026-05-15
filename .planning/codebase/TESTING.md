# Testing Patterns

**Analysis Date:** 2026-05-15

## Test Framework

**Runner:**
- Google Test (GTest) - Used in `TestPluckerTransform.cpp`
- Basic `assert()` - Used in `TestSpatialVector.cpp`
- Config: Linked via CMake (`GTest::GTest`, `GTest::Main`)

**Assertion Library:**
- Google Test assertions: `TEST()`, `EXPECT_DOUBLE_EQ()`
- Standard assert: `assert()` from `<cassert>`

**Run Commands:**
```sh
cmake -B build && cmake --build build
cd build && ctest --output-on-failure
# or run directly:
./build/TestSpatialVector
./build/TestPluckerTransform
```

## Test File Organization

**Location:**
- Separate `tests/` directory
- Co-located with source would be preferred but not current pattern

**Naming:**
- `Test<ClassName>.cpp` pattern
- Examples: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`

**Structure:**
```
tests/
├── TestSpatialVector.cpp      # Basic assert tests (17 lines)
├── TestPluckerTransform.cpp   # GTest tests (63 lines)
├── TestArticulatedBodyInertia.cpp  # Empty stub
├── TestRigidBodyInertia.cpp   # Empty stub
└── TestSpatialOperations.cpp  # Empty stub
```

## Test Structure

**GTest Suite Organization:**
```cpp
#include <gtest/gtest.h>
#include "PluckerTransform.h"

using namespace SpatialAlgebra;

TEST(PluckerTransform, TransformMotion)
{
    // Arrange
    Rotation E{Eigen::Matrix3d::Identity()};
    PluckerTransform transform(E, Vector3d::Zero());
    SpatialVector motion({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});

    // Act
    SpatialVector transformed = transform.transformMotion(motion);

    // Assert
    EXPECT_DOUBLE_EQ(transformed.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(transformed.getAngular()[1], 2.0);
    // ...
}
```

**Assert-Based Tests:**
```cpp
#include "SpatialVector.h"
#include <cassert>

using namespace SpatialAlgebra;

int main()
{
    SpatialVector v1({1, 2, 3}, {4, 5, 6});
    SpatialVector v2 = v1 * 2.0;

    assert(v2.getAngular()[0] == 2);
    assert(v2.getLinear()[1] == 10);

    v1.print();
    return 0;
}
```

**Patterns:**
- No setup/teardown fixtures currently used
- No parameterized tests
- No test helpers or utilities

## Mocking

**Framework:** None

**Patterns:**
- No mocking infrastructure
- All tests use real implementations
- Eigen types not mocked

**What to Mock:**
- Not applicable - library has no external dependencies to mock

**What NOT to Mock:**
- Core algebra operations (should test real behavior)

## Fixtures and Factories

**Test Data:**
```cpp
// Inline object creation
SpatialVector motion({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});
Rotation E{Eigen::Matrix3d::Identity()};
PluckerTransform transform(E, Vector3d::Zero());
```

**Location:**
- No shared fixtures
- Each test creates its own objects

## Coverage

**Requirements:** None enforced

**Current State:**
- Only 2 of 5 test files have implementations
- `TestArticulatedBodyInertia.cpp` - Empty (0 lines)
- `TestRigidBodyInertia.cpp` - Empty (0 lines)
- `TestSpatialOperations.cpp` - Empty (0 lines)
- `AxialScrewTransform.h/.cpp` - Empty stubs (0 lines)

**View Coverage:**
- No coverage tool configured
- Would require: `cmake --build build -DCMAKE_BUILD_TYPE=Debug --target coverage`

## Test Types

**Unit Tests:**
- Scope: Individual class methods
- Approach: Direct instantiation and method calls
- Examples: Testing `transformMotion()`, operator overloads

**Integration Tests:**
- None currently
- Would test interactions between multiple classes

**E2E Tests:**
- Not applicable (library, not application)

## Common Patterns

**Identity Tests:**
```cpp
// Test with identity transform
Rotation E{Eigen::Matrix3d::Identity()};
PluckerTransform transform(E, Vector3d::Zero());
// Expect input == output
```

**Type Alias Usage:**
```cpp
// Using concise aliases in tests
mv m1({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});
fv f1({2.0, 4.0, 6.0}, {4.0, 5.0, 6.0});
plux p1(E1, Vector3d::Zero());
```

**Print for Debugging:**
```cpp
m2.print();
f2.print();
p3.print();
```

## Test Registration (CMake)

```cmake
# CMakeLists.txt
add_executable(TestSpatialVector tests/TestSpatialVector.cpp)
add_executable(TestPluckerTransform tests/TestPluckerTransform.cpp)

target_link_libraries(TestSpatialVector 
    SpatialAlgebra 
    GTest::GTest 
    GTest::Main
)

target_link_libraries(TestPluckerTransform
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestSpatialVector COMMAND TestSpatialVector)
add_test(NAME TestPluckerTransform COMMAND TestPluckerTransform)
```

## Test Gaps

**Missing Tests:**
- `Rotation` class - No dedicated tests
- `LowerTriangular` class - No dedicated tests
- `RigidBodyInertia` class - Stub file only
- `ArticulatedBodyInertia` class - Stub file only
- `SpatialOperations` class - Stub file only
- `SpatialUtils` functions - No dedicated tests

**Recommendations:**
1. Implement tests in stub files
2. Add property-based tests for algebraic laws
3. Add numerical precision tests
4. Test edge cases (zero vectors, identity transforms)

---

*Testing analysis: 2026-05-15*
