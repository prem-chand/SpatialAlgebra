# Common shared utilities library for benchmarks
# Sources are auto-discovered via GLOB — Plan 02 adds model_factory.cpp and random_state.cpp
file(GLOB COMMON_SOURCES *.cpp)
add_library(bench_common STATIC ${COMMON_SOURCES})
target_link_libraries(bench_common PUBLIC SpatialAlgebra Eigen3::Eigen)
target_include_directories(bench_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
