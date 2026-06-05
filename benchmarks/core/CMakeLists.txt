# Core microbenchmark sources — auto-discovered via GLOB
# Plan 03 adds bench_core_stub.cpp; Phase 17 replaces stub with real benchmark
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
set(CORE_SOURCES ${SOURCES} PARENT_SCOPE)
