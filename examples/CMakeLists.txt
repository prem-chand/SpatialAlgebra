# Examples CMakeLists.txt
# Build usage examples for SpatialAlgebra library

set(CMAKE_CXX_STANDARD 17)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Basic vectors example
add_executable(example_vectors basic_vectors.cpp)
target_link_libraries(example_vectors SpatialAlgebra Eigen3::Eigen)

# Transforms example
add_executable(example_transforms transforms.cpp)
target_link_libraries(example_transforms SpatialAlgebra Eigen3::Eigen)

# Inertia example
add_executable(example_inertia inertia.cpp)
target_link_libraries(example_inertia SpatialAlgebra Eigen3::Eigen)

# Dynamics example (ABA)
add_executable(example_dynamics dynamics.cpp)
target_link_libraries(example_dynamics SpatialAlgebra Eigen3::Eigen)
