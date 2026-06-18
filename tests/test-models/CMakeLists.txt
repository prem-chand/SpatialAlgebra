# Zero-dependency test model library.
# test_models INTERFACE target links only Eigen3::Eigen — no SpatialAlgebra.
# sa_test_adapter is the sole bridge between test-models and libSpatialAlgebra.
#
# Phase 20 (Plan 01) creates the CMake infrastructure and header-only types.
# Plan 04 adds sa_adapter.cpp and updates compile_smoke_test.cpp.
# Plan 02-03 add the 11 chains/ model factory headers.

# ---------------------------------------------------------------------------
# INTERFACE library: test_models (Eigen3 only — zero SpatialAlgebra dependency)
# Header-only. Consumers include robot_model.h and robot_solver.h.
# ---------------------------------------------------------------------------
add_library(test_models INTERFACE)
target_include_directories(test_models INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_models INTERFACE Eigen3::Eigen)

# ---------------------------------------------------------------------------
# Compiled adapter: sa_test_adapter (bridges test_models ↔ SpatialAlgebra)
# sa_adapter.cpp is created in Plan 04. Target declared now so CMake configure
# succeeds; CMake evaluates the source list at build time, not configure time.
# ---------------------------------------------------------------------------
add_library(sa_test_adapter sa_adapter.cpp)
target_link_libraries(sa_test_adapter PUBLIC test_models SpatialAlgebra)
target_include_directories(sa_test_adapter PRIVATE ${CMAKE_SOURCE_DIR}/include)

# ---------------------------------------------------------------------------
# Pinocchio cross-validation adapter: pinocchio_adapter (bridges test_models ↔ Pinocchio)
# Guarded by SA_BUILD_PINOCCHIO_BENCHMARKS (default OFF, set in root CMakeLists.txt).
# Requires Pinocchio C++ library from https://stack-of-tasks.github.io/pinocchio/.
# The pinocchio::pinocchio IMPORTED target is created in root CMakeLists.txt via
# find_path/find_library (NOT find_package) to bypass pinocchioConfig.cmake's
# find_package(Boost REQUIRED COMPONENTS system) call, which fails on Homebrew
# Boost 1.89.0 because boost_system is header-only (no libboost_system.dylib).
# ---------------------------------------------------------------------------
if(SA_BUILD_PINOCCHIO_BENCHMARKS)
    add_library(pinocchio_test_adapter pinocchio_adapter.cpp)
    target_link_libraries(pinocchio_test_adapter PUBLIC test_models pinocchio::pinocchio)
    target_include_directories(pinocchio_test_adapter PRIVATE ${CMAKE_SOURCE_DIR}/include)

    # Cross-validation comparison benchmark (GTest)
    add_executable(bench_pinocchio_comparison
        ${CMAKE_SOURCE_DIR}/tests/bench_pinocchio_comparison.cpp
    )
    target_link_libraries(bench_pinocchio_comparison PRIVATE
        test_models
        sa_test_adapter
        pinocchio_test_adapter
        SpatialAlgebra
        GTest::GTest
        GTest::Main
    )
    target_include_directories(bench_pinocchio_comparison PRIVATE
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME bench_pinocchio_comparison COMMAND bench_pinocchio_comparison)
endif()
# ---------------------------------------------------------------------------
# Smoke test: verifies zero-dependency guarantee compiles
# Links ONLY test_models (NOT SpatialAlgebra). If any SA type leaks into
# the test_models headers, this target fails to build.
# compile_smoke_test.cpp is updated in Plan 04 to include robot_model.h
# and robot_solver.h instead of SpatialAlgebra.h.
# ---------------------------------------------------------------------------
add_executable(compile_smoke_test
    ${CMAKE_SOURCE_DIR}/tests/compile_smoke_test.cpp
)
target_link_libraries(compile_smoke_test PRIVATE test_models)
add_test(NAME compile_smoke_test COMMAND compile_smoke_test)

# ---------------------------------------------------------------------------
# Install rules: standalone installable library (TML-05)
# ---------------------------------------------------------------------------
install(TARGETS test_models EXPORT SpatialAlgebraTargets)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/test-models
    FILES_MATCHING PATTERN "*.h"
)
