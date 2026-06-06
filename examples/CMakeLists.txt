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

# Robot example: 2-link Z-Z planar arm
add_executable(example_robot_2link example_robot_2link.cpp)
target_link_libraries(example_robot_2link SpatialAlgebra Eigen3::Eigen)

# Robot example: 3-link Z-Y-Z spatial arm
add_executable(example_robot_3link example_robot_3link.cpp)
target_link_libraries(example_robot_3link SpatialAlgebra Eigen3::Eigen)
