---
quick_id: 260618-hap
slug: fix-cmake-configure-failure-homebrew-boo
description: Fix CMake configure failure — Homebrew Boost 1.89.0 does not ship boost_systemConfig.cmake
date: 2026-06-18
status: ready
must_haves:
  truths:
    - cmake -B build -DSA_BUILD_PINOCCHIO_BENCHMARKS=ON succeeds without fake system boost_systemConfig.cmake
    - tests/test-models/CMakeLists.txt comment is accurate about how pinocchio is located
    - CLAUDE.md documents the Boost/pinocchio workaround
  artifacts:
    - CMakeLists.txt: find_path/find_library approach for pinocchio (no find_package(pinocchio))
    - tests/test-models/CMakeLists.txt: corrected comment
    - CLAUDE.md: updated with configure workaround note
---

# Quick Task 260618-hap: Fix CMake Boost Configure Failure

## Root Cause

`find_package(pinocchio REQUIRED)` chains into pinocchioConfig.cmake, which calls
`find_package(Boost REQUIRED COMPONENTS system)` with `Boost_NO_BOOST_CMAKE ON` (find-module
mode). FindBoost.cmake then looks for a `libboost_system.a/dylib`, but Boost 1.89.0's
`boost_system` is header-only — no library file — so CMake FATAL_ERRORs.

The fake `/usr/local/lib/cmake/boost_system-1.89.0/boost_systemConfig.cmake` does NOT fix
this because pinocchioConfig.cmake forces find-module mode (bypasses the config file).

## Fix Already in Working Tree

`CMakeLists.txt` working tree already contains the correct fix: uses `find_path/find_library`
directly to locate pinocchio headers and libpinocchio_default, bypassing pinocchioConfig.cmake
entirely so Boost is never needed during cmake configuration.

## Tasks

### Task 1: Fix wrong comment in tests/test-models/CMakeLists.txt

**File:** `tests/test-models/CMakeLists.txt`
**Action:** Update comment that incorrectly says "find_package(pinocchio REQUIRED) is called
from root CMakeLists.txt" — the root CMakeLists.txt uses find_path/find_library, NOT
find_package. Update comment to reflect the find_path/find_library bypass approach.
**Verify:** Comment accurately describes how pinocchio::pinocchio is constructed.
**Done:** grep shows no "find_package(pinocchio REQUIRED)" in comment.

### Task 2: Update CLAUDE.md with Boost/configure workaround note

**File:** `CLAUDE.md`
**Action:** Add a note under the existing "Eigen version mismatch" note in the Build & Test
section explaining the Boost issue. Document that pinocchio is located via find_path/find_library
(not find_package) to avoid pinocchioConfig.cmake's FindBoost.cmake requirement, which fails
because Homebrew Boost 1.89.0's boost_system is header-only (no libboost_system.dylib).
**Verify:** CLAUDE.md has the note.
**Done:** Note visible under build section.

### Task 3: Commit the working tree CMakeLists.txt fix

**Files:** `CMakeLists.txt`, `tests/test-models/CMakeLists.txt`, `CLAUDE.md`
**Action:** Stage and commit these three files with message:
  `fix(cmake): bypass pinocchioConfig.cmake to avoid Boost 1.89.0 header-only system failure`
**Verify:** `git log --oneline -1` shows the commit.
**Done:** Commit hash recorded in SUMMARY.md.
