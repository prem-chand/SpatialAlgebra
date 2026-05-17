cmake_minimum_required(VERSION 3.10)
project(SpatialAlgebra)

set(CMAKE_CXX_STANDARD 17)

#  Use system default compiler (removed hardcoded g++)

# Include directories
include_directories(include)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find GTest with FetchContent fallback
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/release-1.12.1.zip
    )
    FetchContent_MakeAvailable(googletest)
    include_directories(${googletest_SOURCE_DIR}/googletest/include)
endif()

# Source files
file(GLOB SOURCES "src/*.cpp")

# Remove empty stub files from compilation
list(REMOVE_ITEM SOURCES
    "${CMAKE_SOURCE_DIR}/src/RigidBodyInertia.cpp"
    "${CMAKE_SOURCE_DIR}/src/ArticulatedBodyInertia.cpp"
)

# Build library
add_library(SpatialAlgebra ${SOURCES})

# Link Eigen3
target_link_libraries(SpatialAlgebra Eigen3::Eigen)

# Enable code coverage (gcov)
option(ENABLE_COVERAGE "Enable coverage flags for CI" OFF)
if(ENABLE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# Build tests
enable_testing()

# Add test executables
add_executable(TestSpatialVector tests/TestSpatialVector.cpp)
add_executable(TestPluckerTransform tests/TestPluckerTransform.cpp)
add_executable(TestRotation tests/TestRotation.cpp)
add_executable(TestLowerTriangular tests/TestLowerTriangular.cpp)
add_executable(TestSpatialUtils tests/TestSpatialUtils.cpp)

# Link test executables with Google Test and the main library
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

target_link_libraries(TestRotation
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

target_link_libraries(TestLowerTriangular
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

target_link_libraries(TestSpatialUtils
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_executable(TestRigidBodyInertia tests/TestRigidBodyInertia.cpp)

target_link_libraries(TestRigidBodyInertia
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

# Add tests to CTest
add_test(NAME TestSpatialVector COMMAND TestSpatialVector)
add_test(NAME TestPluckerTransform COMMAND TestPluckerTransform)
add_test(NAME TestRotation COMMAND TestRotation)
add_test(NAME TestLowerTriangular COMMAND TestLowerTriangular)
add_test(NAME TestSpatialUtils COMMAND TestSpatialUtils)

add_test(NAME TestRigidBodyInertia COMMAND TestRigidBodyInertia)

add_executable(TestArticulatedBodyInertia tests/TestArticulatedBodyInertia.cpp)

target_link_libraries(TestArticulatedBodyInertia
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestArticulatedBodyInertia COMMAND TestArticulatedBodyInertia)

add_executable(TestForwardDynamics tests/TestForwardDynamics.cpp)

target_link_libraries(TestForwardDynamics
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestForwardDynamics COMMAND TestForwardDynamics)

add_executable(TestSpatialOperations tests/TestSpatialOperations.cpp)

target_link_libraries(TestSpatialOperations
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestSpatialOperations COMMAND TestSpatialOperations)

add_executable(TestInverseDynamics tests/TestInverseDynamics.cpp)

target_link_libraries(TestInverseDynamics
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestInverseDynamics COMMAND TestInverseDynamics)

add_executable(TestDynamicsConsistency tests/TestDynamicsConsistency.cpp)

target_link_libraries(TestDynamicsConsistency
    SpatialAlgebra
    GTest::GTest
    GTest::Main
)

add_test(NAME TestDynamicsConsistency COMMAND TestDynamicsConsistency)

# Build examples
add_subdirectory(examples)
