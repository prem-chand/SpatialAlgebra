# Testing Patterns

**Analysis Date:** 2026-05-17

## Test Framework

**Runner:**
- Google Test (GTest) — `gtest` library
- Version: As provided by Homebrew (`brew install googletest`)
- Config: Test executables defined in `CMakeLists.txt` (lines 30-132)
- CTest integration via `enable_testing()` (line 27) and `add_test()` calls (lines 76-132)

**Build Commands:**
```bash
cmake -B build                           # Configure
cmake --build build                      # Build tests
cd build && ctest --output-on-failure    # Run all registered tests
# Run individual test executables:
build/TestSpatialVector
build/TestPluckerTransform
build/TestRotation
```

**Assertion Library:**
- GTest built-in: `EXPECT_DOUBLE_EQ`, `EXPECT_NEAR`, `EXPECT_THROW`, `EXPECT_TRUE`, `EXPECT_GT`, `EXPECT_FALSE`, `EXPECT_EQ`
- No external assertion library (e.g., Catch2, doctest not used)

## Test File Organization

**Location:**
- All test files in `tests/` directory:
```
tests/
├── TestSpatialVector.cpp       (609 lines)
├── TestPluckerTransform.cpp    (946 lines)
├── TestRotation.cpp            (413 lines)
├── TestLowerTriangular.cpp     (450 lines)
├── TestSpatialUtils.cpp        (335 lines)
├── TestRigidBodyInertia.cpp    (376 lines)
├── TestArticulatedBodyInertia.cpp (531 lines)
├── TestForwardDynamics.cpp     (253 lines)
├── TestInverseDynamics.cpp     (192 lines)
├── TestDynamicsConsistency.cpp (266 lines)
├── TestSpatialOperations.cpp   (320 lines)
```

**Naming:**
- `Test<ClassName>.cpp` for each class-specific test suite
- `TestDynamicsConsistency.cpp` for cross-component integration tests
- `TestSpatialUtils.cpp` covers both `SpatialUtils.h` free functions and `SpatialOperations.h`

**Structure within each test file:**
1. Includes (header under test, `<gtest/gtest.h>`, Eigen headers, optional standard headers)
2. `using namespace` declarations
3. Tolerance constant definition (typically `EPSILON` or `TOLERANCE`)
4. Section-separator comments (`// =====` style)
5. Test cases (Doxygen-documented TEST blocks)
6. `main()` function with `InitGoogleTest` + `RUN_ALL_TESTS`

## Test Structure

**Suite Organization:**
- Class-scoped test suites use `TEST(SuiteName, TestName)`:
```cpp
TEST(TestSpatialVector, Constructor) { ... }
TEST(TestSpatialVector, Getters) { ... }
TEST(TestSpatialVector, Addition) { ... }
```
- Functional-group suites use descriptive names:
```cpp
TEST(TransformMotionTest, IdentityTransform) { ... }
TEST(TransformMotionTest, PureRotation) { ... }
TEST(TransformMotionTest, Property_Linearity) { ... }
```

**Test Fixture Classes** (used in 4 of 11 test files):
```cpp
class TestCrossProductMotion : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};
```
Then: `TEST_F(TestCrossProductMotion, SimpleRotationVectors) { ... }`

Fixtures are used in these files:
- `tests/TestPluckerTransform.cpp` — `TestInverse`, `TestMultiply`, `TestPrint`
- `tests/TestSpatialUtils.cpp` — `TestSkew`, `TestDot`, `TestCross`, `TestSpatialOperations`
- `tests/TestSpatialOperations.cpp` — `TestCrossProductMotion`, `TestCrossProductForce`, `TestTransformInertia`

**Main Entry Point** (every test file):
```cpp
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

**Test Categories:**
- **Constructor tests:** Verify default and parameterized constructors store correct values
- **Getter tests:** Verify accessors return stored data correctly
- **Arithmetic tests:** Add, subtract, scale operations with known values
- **Property tests:** Mathematical invariants (commutativity, anti-commutativity, associativity, distributivity)
- **Textbook example tests:** Known problems from Featherstone 2008 (e.g., `FeatherstoneExample2_1`, `WrenchExample`)
- **Edge case tests:** Zero input, identity transform, pure rotation/pure translation
- **Linearity/proportionality tests:** Double input → double output
- **Round-trip tests:** Apply transform then inverse returns original
- **Consistency tests:** Forward dynamics + inverse dynamics consistency (cross-component)

## Floating-Point Tolerance

Throughout the tested codebase, two patterns exist:

**Pattern 1 — `constexpr double EPSILON`** (preferred in newer tests):
```cpp
constexpr double EPSILON = 1e-10;
```
Used in: `TestPluckerTransform.cpp`, `TestSpatialOperations.cpp`, `TestForwardDynamics.cpp`, `TestInverseDynamics.cpp`, `TestDynamicsConsistency.cpp`

**Pattern 2 — `const double TOLERANCE`** (in older tests):
```cpp
const double TOLERANCE = 1e-10;
```
Used in: `TestRotation.cpp`, `TestLowerTriangular.cpp`, `TestRigidBodyInertia.cpp`, `TestArticulatedBodyInertia.cpp`

Both patterns use `1e-10`. The `DynamicsConsistency` test uses a looser `1e-8`.

**Assertion Usage:**
- `EXPECT_DOUBLE_EQ(a, b)` — Exact double equality (used when no computation error, e.g., getter return values)
- `EXPECT_NEAR(a, b, EPSILON)` — Floating-point approximate equality (preferred for computed results)
- `EXPECT_THROW(expr, exception_type)` — Exception testing (e.g., `EXPECT_THROW(L(0,1)=5.0, std::invalid_argument)`)
- `EXPECT_TRUE(condition)` / `EXPECT_FALSE(condition)` — Boolean conditions
- `EXPECT_GT(a, b)` / `EXPECT_LT(a, b)` / `EXPECT_EQ(a, b)` — Comparison assertions
- `EXPECT_NE(output.find("Rotation"), std::string::npos)` — String matching (in `TestPluckerTransform.cpp`)
- `.norm()` comparison for vector results: `EXPECT_NEAR((a - b).norm(), 0.0, EPSILON)`

## Mocking

**Framework:** None used.

**What to Mock:**
- No mocking framework (Mockito, GMock, etc.) is used or installed
- No interfaces or abstract classes designed for mocking exist
- `ForwardDynamics.h` uses `ArticulatedBodyInertia` and `PluckerTransform` directly (no dependency injection)

**What NOT to Mock:**
- All classes are concrete types with value semantics
- Tests create real instances with test data
- Mathematical verification relies on known expected values, not mock expectations

## Fixtures and Factories

**Test Data Pattern:**
No dedicated test data files or factory functions. Each test creates its own data inline:

```cpp
TEST(TestSpatialVector, Addition)
{
    SpatialVector v1(Vector3d(1.0, 2.0, 3.0), Vector3d(4.0, 5.0, 6.0));
    SpatialVector v2(Vector3d(2.0, 3.0, 4.0), Vector3d(5.0, 6.0, 7.0));
    SpatialVector sum = v1 + v2;
    EXPECT_DOUBLE_EQ(sum.getAngular()[0], 3.0);
    // ...
}
```

**Helper Functions** (defined in test file, not in a shared header):
```cpp
// TestSpatialOperations.cpp:
LowerTriangular createIdentityInertia() {
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    return LowerTriangular::fromFullMatrix(identity);
}

LowerTriangular createDiagonalInertia(double value) {
    Eigen::Matrix3d diagonal = Eigen::Matrix3d::Identity() * value;
    return LowerTriangular::fromFullMatrix(diagonal);
}
```

**Location:** Helper functions are file-local in test files — no shared test utility directory.

## Coverage

**Requirements:** None enforced.

- No coverage tools (gcov, lcov, CodeCov) configured
- No coverage CI step
- No coverage threshold defined

**View Coverage:**
```bash
# Not configured — would need gcov + manual setup
```

## Test Types

**Unit Tests:**
- Scope: Each test file covers a single class or header
- Approach: White-box testing with known mathematical formulas
- Pattern: Create objects with known values, assert expected results

**Integration Tests:**
- `tests/TestDynamicsConsistency.cpp` — Cross-component test verifying forward dynamics (ABA) and inverse dynamics (RNEA) produce consistent results
- Tests round-trip: `ABA(RNEA(tau)) ≈ tau` and `RNEA(ABA(qddot)) ≈ qddot`
- Tests multi-link serial chains (3-link) and branching trees (Y-configuration)
- Tests PluckerTransform usage within dynamics algorithms

**E2E Tests:** Not used.

**Python Tests:** Not used. `robot_dynamics/rnea.py` is a standalone implementation with no test coverage.

## Common Patterns

**Arrange-Act-Assert Comments:**
Many test files annotate each test with `// Arrange`, `// Act`, `// Assert` sections:
```cpp
TEST(TestCrossProductMotion, SimpleRotationVectors) {
    // Arrange: pure angular velocities about X and Y axes
    MotionVector v1(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    MotionVector v2(Vector3d(0, 1, 0), Vector3d(0, 0, 0));
    
    // Act
    SpatialVector result = SpatialOperations::crossProductMotion(v1, v2);
    
    // Assert: ω1×ω2 = (1,0,0)×(0,1,0) = (0,0,1)
    EXPECT_NEAR(result.getAngular()[2], 1.0, EPSILON);
}
```

**Property-Based Tests:**
Tests verify mathematical invariants rather than specific values:
```
CrossProductAntiCommutativity     a×b = -(b×a)
DotProductCommutativity           a·b = b·a
CrossProductDistributivity        a×(b+c) = a×b + a×c
ScalarMultiplicationProperty      (k*a)×b = k*(a×b)
MassConservation                  mass unchanged by Plücker transform
InverseRoundTrip                  X^(-1)(X(v)) = v
LinearityProperty                 f(a+b) = f(a) + f(b)
```

**Manual Math Verification in Comments:**
Tests explicitly show the expected calculation step-by-step:
```cpp
// Expected:
// ω' = R*ω = [0,0,1]
// v' = R*(v - r×ω) = R*([1,0,0] - [1,0,0]×[0,0,1])
//    = R*([1,0,0] - [0,-1,0]) = R*[1,1,0]
// R*[1,1,0] = [-1,1,0] (90° rotation)
```

**Dynamics Test Pattern:**
Forward/Inverse dynamics tests use a struct-based robot model setup:
```cpp
ForwardDynamics fd;
Link link;
link.parent = -1;
link.X = PluckerTransform(Rotation(Eigen::Matrix3d::Identity()), Vector3d::Zero());
link.I = RigidBodyInertia(1.0, Vector3d::Zero(), lt::Identity(3));
link.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
link.q = 0.0;
link.qdot = 0.0;
link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());
fd.links.push_back(link);
```

## CMake Test Registration

Each test executable is registered in `CMakeLists.txt` with:
1. `add_executable(TestName tests/TestName.cpp)`
2. `target_link_libraries(TestName SpatialAlgebra GTest::GTest GTest::Main)`
3. `add_test(NAME TestName COMMAND TestName)`

Currently registered test executables (10 total):
`TestSpatialVector`, `TestPluckerTransform`, `TestRotation`, `TestLowerTriangular`, `TestSpatialUtils`, `TestRigidBodyInertia`, `TestArticulatedBodyInertia`, `TestForwardDynamics`, `TestSpatialOperations`, `TestInverseDynamics`, `TestDynamicsConsistency`

## Known Issues

- `TestSpatialVector.cpp` uses `.eval()` on `.cross()` results (line 45 in `SpatialVector.cpp`), which is an anti-pattern — Eigen cross product returns an expression template that does not require `.eval()`
- No test for `src/main.cpp` — the example program is not tested
- No test for Python `robot_dynamics/rnea.py`
- `TestSpatialVector.cpp` property tests use `EXPECT_DOUBLE_EQ` for cross product results, which may be fragile for computed values; newer tests correctly use `EXPECT_NEAR`
- Test tolerance `EPSILON = 1e-10` is uniformly applied but some operations (dynamics round-trip) use looser `1e-8`

---

*Testing analysis: 2026-05-17*
