# Deferred Items — Phase 13-05

## Pre-existing Test Failures

### 1. TestDynamicsConsistency — 3 failing sub-tests
- **File:** `tests/TestDynamicsConsistency.cpp`
- **Failing tests:** `ThreeLinkSerialChain`, `BranchingYConfiguration`, `TwoLinkRoundTripWithGravity`
- **Issue:** ABA qddot values do not match expected inputs for multi-link chains with gravity
- **Root cause:** Pre-existing dynamics algorithm accuracy issue (not caused by phase 13 changes)
- **When to resolve:** Phase 12 (dynamics consistency) or dedicated algorithm bug-fix phase
- **Scope boundary:** Out of scope — our changes are code quality (namespace, OpenMP, CMake, CI)
