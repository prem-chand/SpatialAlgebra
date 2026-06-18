---
quick_id: 260618-hap
slug: fix-cmake-configure-failure-homebrew-boo
description: Fix CMake configure failure — Homebrew Boost 1.89.0 does not ship boost_systemConfig.cmake
date: 2026-06-18
status: complete
commit: 8bbd433
---

# Summary

## Root Cause Found

`find_package(pinocchio REQUIRED)` → pinocchioConfig.cmake → `find_package(Boost REQUIRED COMPONENTS system)` with `Boost_NO_BOOST_CMAKE ON` (find-module mode) → FindBoost.cmake looks for `libboost_system.dylib` → not found (header-only in Boost 1.89.0) → FATAL_ERROR.

The system had a fake `boost_systemConfig.cmake` as a workaround, but it doesn't help because pinocchioConfig.cmake forces find-module mode that ignores config files.

## Fix Applied

The working tree already had the correct fix in CMakeLists.txt (using `find_path`/`find_library` to construct a manual `pinocchio::pinocchio` IMPORTED target, bypassing pinocchioConfig.cmake entirely). This commit finalizes the fix:

1. **CMakeLists.txt** — `find_path`/`find_library` approach confirmed correct, bypasses Boost chain entirely
2. **tests/test-models/CMakeLists.txt** — fixed wrong comment that said `find_package(pinocchio REQUIRED)` was called (it is not)
3. **CLAUDE.md** — documented the Boost/pinocchio workaround to prevent future regression

## Verification

`cmake -B build -DSA_BUILD_PINOCCHIO_BENCHMARKS=ON` configures successfully without the system-level fake `boost_systemConfig.cmake`.
