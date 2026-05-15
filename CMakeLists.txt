cmake_minimum_required(VERSION 3.10)
project(SpatialAlgebra)

set(CMAKE_CXX_STANDARD 17)

#  set gcc compiler
set(CMAKE_CXX_COMPILER "g++")

# Include directories
include_directories(include)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find GTest
find_package(GTest REQUIRED)

# Source files
file(GLOB SOURCES "src/*.cpp")

# Build library
add_library(SpatialAlgebra ${SOURCES})

# Link Eigen3
target_link_libraries(SpatialAlgebra Eigen3::Eigen)

# Build tests
enable_testing()

# Add test executables
add_executable(TestSpatialVector tests/TestSpatialVector.cpp)
add_executable(TestPluckerTransform tests/TestPluckerTransform.cpp)
add_executable(TestRotation tests/TestRotation.cpp)
add_executable(TestLowerTriangular tests/TestLowerTriangular.cpp)

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

# Add tests to CTest
add_test(NAME TestSpatialVector COMMAND TestSpatialVector)
add_test(NAME TestPluckerTransform COMMAND TestPluckerTransform)
add_test(NAME TestRotation COMMAND TestRotation)
add_test(NAME TestLowerTriangular COMMAND TestLowerTriangular)
