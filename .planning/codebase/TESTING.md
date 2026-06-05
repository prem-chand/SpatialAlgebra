# Testing Patterns

**Analysis Date:** 2026-06-05

## Test Framework

**Runner:**
- Google Test (GTest) — version release-1.12.1 (via FetchContent fallback)
- Config: Built-in CMake `enable_testing()` + `add_test()` commands in `CMakeLists.txt`
- System or Homebrew GTest is tried first via `find_package(GTest QUIET)`, fallback to FetchContent

**Assertion Library:**
- GTest built-in macros: `EXPECT_DOUBLE_EQ`, `EXPECT_NEAR`, `EXPECT_EQ`, `EXPECT_TRUE`, `EXPECT_FALSE`, `ASSERT_NEAR`, `EXPECT_THROW`
- No additional assertion libraries (no Catch2, no Boost.Test)

**Run Commands:**
```bash
cmake --build build && cd build && ctest --output-on-failure    # Run all tests
cmake --build build && cd build && ctest -R TestRotation        # Run a single test suite
cd build && ./TestSpatialVector                                  # Run a single test executable directly
cd build && ./TestForwardDynamics --gtest_filter=*TwoLink*      # GTest filter pattern
cd build && ./TestArticulatedBodyInertia --gtest_list_tests     # List test names
```

**Coverage:**
```bash
cmake -B build -DENABLE_COVERAGE=ON
cmake --build build
cd build && ctest --output-on-failure
# then use gcov or lcov to generate coverage reports
```

## Test File Organization

**Location:** All test files in `tests/` directory — co-located at project level, not with source.

**Naming:**
- Files prefixed with `Test` + class name: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`, `TestRotation.cpp`, `TestLowerTriangular.cpp`, `TestRigidBodyInertia.cpp`, `TestArticulatedBodyInertia.cpp`, `TestSpatialUtils.cpp`, `TestSpatialOperations.cpp`, `TestInverseDynamics.cpp`, `TestForwardDynamics.cpp`, `TestDynamicsConsistency.cpp`
- Each file corresponds to one class or one module

**Executable Registration:**
- Each test file maps to a separate CMake executable target defined in `CMakeLists.txt`
- Each target links against `SpatialAlgebra`, `GTest::GTest`, and `GTest::Main`
- Each target has a corresponding `add_test()` registration for `ctest`

**Registered Test Executables (from `CMakeLists.txt`):**
| Executable | Lines in CMake |
|---|---|
| `TestSpatialVector` | `add_executable` at line 52 |
| `TestPluckerTransform` | line 53 |
| `TestRotation` | line 54 |
| `TestLowerTriangular` | line 55 |
| `TestSpatialUtils` | line 56 |
| `TestRigidBodyInertia` | line 89 |
| `TestArticulatedBodyInertia` | line 106 |
| `TestForwardDynamics` | line 116 |
| `TestSpatialOperations` | line 126 |
| `TestInverseDynamics` | line 136 |
| `TestDynamicsConsistency` | line 146 |

**Build Artifacts (from `build/`):**
```
TestArticulatedBodyInertia      (1.2 MB)
TestDynamicsConsistency         (3.2 MB)
TestForwardDynamics             (3.2 MB)
TestInverseDynamics             (2.7 MB)
TestLowerTriangular             (790 KB)
TestPluckerTransform            (6.2 MB)
TestRigidBodyInertia            (919 KB)
TestRotation                    (1.2 MB)
TestSpatialOperations           (2.7 MB)
TestSpatialUtils                (2.6 MB)
TestSpatialVector               (843 KB)
```

## Test Structure

**Suite Organization:**

```cpp
// ============================================================================
// [Category Name] Tests
// ============================================================================

/**
 * @brief Test [what is being tested]
 * @details [Scenario description]
 */
TEST(TestSuiteName, TestCaseName)
{
    // [optional: Arrange]
    Vector3d angular(1.0, 2.0, 3.0);
    Vector3d linear(4.0, 5.0, 6.0);

    // [Act]
    SpatialVector v(angular, linear);

    // [Assert]
    EXPECT_DOUBLE_EQ(v.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(v.getAngular()[1], 2.0);
    EXPECT_DOUBLE_EQ(v.getAngular()[2], 3.0);
}
```

**Suite Name Patterns (observed):**
- `TestSpatialVector` — tests for `SpatialVector` class (no fixture)
- `TestMotionVector` — tests for `MotionVector` class (no fixture)
- `TestForceVector` — tests for `ForceVector` class (no fixture)
- `TransformMotionTest` — tests for `transformMotion` method (`TestPluckerTransform.cpp`)
- `TransformForceTest` — tests for `transformForce` method (`TestPluckerTransform.cpp`)
- `PluckerPropertyTest` — property-based tests for Plücker transforms (`TestPluckerTransform.cpp`)
- `RotationTest` — tests for `Rotation` class (`TestRotation.cpp`)
- `LowerTriangularTest` — tests for `LowerTriangular` class (`TestLowerTriangular.cpp`)
- `RigidBodyInertiaTest` — tests for `RigidBodyInertia` class (`TestRigidBodyInertia.cpp`)
- `ArticulatedBodyInertiaTest` — tests for `ArticulatedBodyInertia` class (`TestArticulatedBodyInertia.cpp`)
- `InverseDynamicsTest` — tests for `InverseDynamics` (`TestInverseDynamics.cpp`)
- `ForwardDynamicsTest` — tests for `ForwardDynamics` (`TestForwardDynamics.cpp`)
- `ConsistencyTest` — round-trip tests (`TestDynamicsConsistency.cpp`)
- `TestSkew`, `TestDot`, `TestCross` — fixture-based test suites (`TestSpatialUtils.cpp`)
- `TestCrossProductMotion`, `TestCrossProductForce`, `TestTransformInertia` — fixture-based test suites (`TestSpatialOperations.cpp`)

**Test Case Name Patterns:**
- VerbNoun: `Constructor`, `Getters`, `Addition`, `Subtraction`, `DotProduct`, `CrossProductAntiCommutativity`
- Scenario-based: `CrossMotionWithLinearComponents`, `SingleLinkPendulum`, `TwoLinkSerialChain`
- Property-based: `Property_SkewSymmetric`, `Property_UnitDeterminant`, `Property_Orthogonality`
- Requirement-referencing: `PackedStorageSize` (referencing LTR-01), `PackedStorageIndexing` (LTR-01), `UpperTriangularReturnsZero` (LTR-01)

**Main Function:**
- `TestSpatialVector.cpp` defines its own `main()`:
  ```cpp
  int main(int argc, char **argv)
  {
      ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
  }
  ```
- All other test files rely on `GTest::Main` linkage (no explicit `main()`)

## Assertion Styles

**Patterns used:**

- `EXPECT_DOUBLE_EQ(a, b)` — preferred for exact floating-point equality (zero-error cases)
- `EXPECT_NEAR(a, b, tolerance)` — used with `EPSILON` or `TOLERANCE` constants for computed values
- `ASSERT_NEAR` — used where subsequent assertions depend on the result (less common)
- `EXPECT_TRUE(cond)` / `EXPECT_FALSE(cond)` — boolean conditions
- `EXPECT_EQ(a, b)` — integer equality
- `EXPECT_THROW(expr, exception_type)` — exception testing (`TestLowerTriangular.cpp`)

**Floating-Point Tolerance Constants:**
- `constexpr double EPSILON = 1e-10;` — used in `TestPluckerTransform.cpp`, `TestSpatialOperations.cpp`, `TestInverseDynamics.cpp`, `TestForwardDynamics.cpp`, `TestDynamicsConsistency.cpp`
- `const double TOLERANCE = 1e-10;` — used in `TestRotation.cpp`, `TestLowerTriangular.cpp`, `TestRigidBodyInertia.cpp`, `TestArticulatedBodyInertia.cpp`
- `constexpr double EPSILON = 1e-8;` — slightly looser tolerance in `TestDynamicsConsistency.cpp` (line 8)
- Both `double` and `constexpr double` used inconsistently — no project-wide standard

## Test Data / Fixtures

**No Fixture Tests (most common pattern):**
- Tests use `TEST()` (not `TEST_F()`) and define data inline within each test case
- Example: `TestSpatialVector.cpp`, `TestPluckerTransform.cpp`, `TestRotation.cpp`, `TestLowerTriangular.cpp`, `TestRigidBodyInertia.cpp`, `TestArticulatedBodyInertia.cpp`

**Fixture Tests (used in `TestSpatialUtils.cpp` and `TestSpatialOperations.cpp`):**
```cpp
class TestSkew : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST(TestSkew, CreatesSkewSymmetricMatrix) {
    // Arrange
    Vector3d v(1.0, 2.0, 3.0);
    
    // Act
    Eigen::Matrix3d S = skew(v);
    
    // Assert
    EXPECT_DOUBLE_EQ(S(0, 0), 0.0);
    // ...
}
```

**Note:** Most fixture classes have empty `SetUp()` / `TearDown()`. Fixtures are used as organizational grouping rather than setup/reuse.

**Test Data:**
- All test data defined inline within test functions — no external test data files, no JSON fixtures, no YAML configs
- `TestSpatialOperations.cpp` defines helper functions for creating test matrices:
  ```cpp
  LowerTriangular createIdentityInertia() {
      Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
      return LowerTriangular::fromFullMatrix(identity);
  }
  
  LowerTriangular createDiagonalInertia(double value) {
      Eigen::Matrix3d diagonal = Eigen::Matrix3d::Identity() * value;
      return LowerTriangular::fromFullMatrix(diagonal);
  }
  ```

## Mocking

**No mocking framework used.** No mock objects, no test doubles, no dependency injection patterns.

**Why:** The library is a pure math library with no external dependencies. Classes depend directly on concrete types (Eigen matrices, other SpatialAlgebra classes). There are no interfaces to mock.

**What would benefit from mocking:** If network/database IO were added, or if a plugin architecture were introduced. Currently not applicable.

## Property-Based Tests

Several test files include property-based tests that verify mathematical invariants:

```cpp
// Anti-commutativity: a×b = -(b×a)
TEST(TestMotionVector, CrossProductAntiCommutativity) {
    MotionVector a(...), b(...);
    MotionVector a_cross_b = a.crossMotion(b);
    MotionVector b_cross_a = b.crossMotion(a);
    MotionVector neg_b_cross_a = b_cross_a * -1.0;
    EXPECT_DOUBLE_EQ(a_cross_b.getAngular()[0], neg_b_cross_a.getAngular()[0]);
    // ... all 6 components
}

// Distributivity: a×(b+c) = a×b + a×c
TEST(TestMotionVector, CrossProductDistributivity) { ... }

// Scalar multiplication property: (k*a)×b = k*(a×b)
TEST(TestMotionVector, ScalarMultiplicationProperty) { ... }

// Skew-symmetry: S + S^T = 0
TEST(TestSkew, Property_SkewSymmetric) {
    Eigen::Matrix3d sum = S + S.transpose();
    EXPECT_NEAR(sum.norm(), 0.0, 1e-10);
}

// Determinant of rotation = 1
TEST(RotationTest, Property_UnitDeterminant) {
    Rotation rot(angleAxis);
    EXPECT_NEAR(rot.determinant(), 1.0, TOLERANCE);
}
```

## Round-Trip / Consistency Tests

`TestDynamicsConsistency.cpp` verifies the round-trip property:
- Start with random torque → compute accelerations via ABA → compute torques via RNEA → should get original torque
- Start with random acceleration → compute torques via RNEA → compute accelerations via ABA → should get original acceleration

```cpp
TEST(ConsistencyTest, RoundTripABARNEA) {
    Eigen::VectorXd tau_input(1);
    tau_input[0] = 1.0;
    // ... setup ...
    fd.computeAccelerations(tau_input);
    double qddot_result = fd.links[0].qddot;
    // ... inverse dynamics ...
    Eigen::VectorXd tau_output = id.computeTorques(qddot_vec);
    EXPECT_NEAR(tau_output[0], tau_input[0], EPSILON);
}
```

## Coverage

**Requirements:**
- Not enforced (no CI, no coverage gate)
- Optional CMake flag `ENABLE_COVERAGE` in `CMakeLists.txt:42-46`:
  ```cmake
  option(ENABLE_COVERAGE "Enable coverage flags for CI" OFF)
  if(ENABLE_COVERAGE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  endif()
  ```
- Uses `gcov` for profiling (no `gcovr` or `lcov` configured)

**Current state:** No coverage tracking in regular builds. Must enable via `-DENABLE_COVERAGE=ON`.

## Test Types

**Unit Tests (all existing tests):**
- Each test file tests a single class or module in isolation
- Direct construction of objects, method calls, and assertion checks
- Examples: `TestSpatialVector.cpp` (arithmentic, cross product, dot product), `TestRotation.cpp` (constructors, accessors, inverse)
- Mathematical property tests verify invariants

**Integration Tests:**
- `TestInverseDynamics.cpp` — tests RNEA on 1-link, 2-link, 3-link chains, branching trees, edge cases
- `TestForwardDynamics.cpp` — tests ABA on 1-link, 2-link, 3-link chains, branching trees, gravity
- `TestDynamicsConsistency.cpp` — round-trip tests connecting RNEA and ABA
- `TestSpatialOperations.cpp` — cross product and inertia transform using multiple classes
- These tests exercise multiple classes together but are still in the same test process

**E2E Tests:**
- None. No integration with external systems.

## Common Patterns

**Arrange-Act-Assert Comments:**
```cpp
// Arrange
MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
MotionVector v2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));

// Act
SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);

// Assert
EXPECT_NEAR(result.getAngular()[0], 0.0, EPSILON);
```
Used in `TestSpatialOperations.cpp` and `TestSpatialUtils.cpp`. Not used consistently in older test files.

**Doxygen on Tests:**
- Every `TEST()` or `TEST_F()` gets a `@brief` and optional `@details` block
- Test function comments explain what invariant or scenario is being verified
- References to literature: `@see Featherstone 2008, Chapter 2`

**Inline Calculation Comments:**
```cpp
// ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
// ω1×v2 = (1,0,0)×(0,0,1) = (0,-1,0)
// v1×ω2 = (0,1,0)×(0,1,0) = (0,0,0)
// Linear = (0,-1,0) + (0,0,0) = (0,-1,0)
// Result: [(0,0,1); (0,-1,0)]
```

**Error Testing:**
```cpp
// Exception tests
TEST(LowerTriangularTest, InvalidDimensions)
{
    LowerTriangular L(3);
    EXPECT_THROW(L.multiply(LowerTriangular(4)), std::invalid_argument);
}
```
Exception testing used only in `TestLowerTriangular.cpp`.

**Loop-Based Matrix Verification:**
```cpp
// Full matrix element-by-element check
for (int i = 0; i < 3; ++i)
{
    for (int j = 0; j < 3; ++j)
    {
        if (i == j)
            EXPECT_DOUBLE_EQ(rot(i, j), 1.0);
        else
            EXPECT_DOUBLE_EQ(rot(i, j), 0.0);
    }
}
```

**Label-Based Requirements References:**
```cpp
// ===== LTR-01: Packed Storage Tests =====
// ===== LTR-02: Dense Matrix Multiplication Tests =====
// ===== LTR-03: Vector Multiplication Tests =====
// INR-01: Default constructor should create massless body
// VER-01: Default constructor should create zero inertia
```
Labels like `LTR-01`, `INR-01`, `VER-01` appear in section comments to trace test cases back to requirements.

---

*Testing analysis: 2026-06-05*
