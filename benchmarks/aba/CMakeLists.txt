# ABA benchmark sources — auto-discovered via GLOB
# Plan 03 adds bench_aba_stub.cpp; Phase 17 replaces stub with real benchmark
# Use absolute paths so PARENT_SCOPE resolves correctly from benchmarks/CMakeLists.txt
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
set(ABA_SOURCES ${SOURCES} PARENT_SCOPE)
