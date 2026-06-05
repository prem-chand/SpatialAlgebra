cmake_minimum_required(VERSION 3.19)
project(SpatialAlgebraBenchmarks LANGUAGES CXX)

# Suppress Google Benchmark's internal test/install targets
# Without this, FetchContent pulls in GTest as transitive dependency,
# conflicting with the project's own GTest setup (Pitfall 1, RESEARCH.md)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

# Fetch Google Benchmark v1.9.5
include(FetchContent)
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
)
FetchContent_MakeAvailable(googlebenchmark)

# Explicit C++17 — Google Benchmark v1.9.5 requires it to build
set(CMAKE_CXX_STANDARD 17)

# Shared utilities library
add_subdirectory(common)

# Domain benchmark directories
add_subdirectory(aba)
add_subdirectory(rnea)
add_subdirectory(core)

# Collect all benchmark sources from subdirectories
# Each subdirectory CMakeLists.txt sets a PARENT_SCOPE variable
set(BENCH_SOURCES ${ABA_SOURCES} ${RNEA_SOURCES} ${CORE_SOURCES})

# CMake 4.x requires at least one source per target. If no sources exist yet
# (Plans 02-03 add them later), use a placeholder to satisfy CMake.
if(NOT BENCH_SOURCES)
    set(BENCH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/_init.cpp)
endif()

# Unified benchmark executable
# bench_all.cpp is created in Plan 03 (Wave 2) and added via target_sources()
add_executable(bench_all ${BENCH_SOURCES})

# bench_all.cpp is excluded from Plan 01's BENCH_SOURCES because the file doesn't
# exist until Plan 03. Adding it there would cause CMake 4.x "Cannot find source file"
# errors. Instead, append it via target_sources() after the file exists.
target_sources(bench_all PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bench_all.cpp)

# Add include path so subdirectory stubs can resolve #include "common/model_factory.h"
target_include_directories(bench_all PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(bench_all
    PRIVATE
    SpatialAlgebra
    benchmark::benchmark
    Eigen3::Eigen
    bench_common
)

# Benchmark compilation profile per D-12/D-13/D-14:
# -O3 -DNDEBUG applied per-target (not as CMAKE_BUILD_TYPE override)
# -flto on both compile and link per D-14
target_compile_options(bench_all PRIVATE -O3 -DNDEBUG -flto)
target_link_options(bench_all PRIVATE -flto)
