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
