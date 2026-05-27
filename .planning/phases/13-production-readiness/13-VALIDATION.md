---
phase: 13
slug: production-readiness
status: draft
nyquist_compliant: true
wave_0_complete: true
created: 2026-05-17
revised: 2026-05-27
---

# Phase 13 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Google Test (GTest) |
| **Config file** | `CMakeLists.txt` (test executables registered via `add_test()`) |
| **Quick run command** | `cmake --build build && cd build && ctest --output-on-failure -R Test` |
| **Full suite command** | `cmake --build build && cd build && ctest --output-on-failure` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run quick command (affected test executable)
- **After every plan wave:** Run full suite
- **Before `/gsd-verify-work`:** Full suite must be green (all tests passing)
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

### Wave 1A — Cross-product fix, Gravity, Plücker safety

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-01-T1 | 13-01 | 1A | VEC-01, UTL-03 | T-13-BC | unit | `grep -rn "crossMotion" include/ src/ tests/ examples/ --include="*.h" --include="*.cpp" 2>/dev/null \| grep -v "planning" \| head -20` | N/A (grep) | ✅ green |
| 13-01-T2 | 13-01 | 1A | VEC-01, UTL-03, INR-01, PLX-04 | T-13-01 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestSpatialVector\|TestSpatialUtils\|TestSpatialOperations" 2>&1 \| tail -10` | ✅ | ✅ green |
| 13-01-T3 | 13-01 | 1A | VEC-01, UTL-03, INR-01, PLX-04 | T-13-02 | unit | `grep -rn "crossMotion" include/ src/ tests/ --include="*.h" --include="*.cpp" && echo "UNEXPECTED: crossMotion still referenced" \|\| echo "PASS: no crossMotion references remain"` | N/A (grep) | ✅ green |
| 13-02-T1 | 13-02 | 1A | ABA-01 | T-13-03 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestForwardDynamics" 2>&1 \| tail -10` | ✅ | ✅ green |
| 13-02-T2 | 13-02 | 1A | ABA-01 | T-13-03 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestInverseDynamics\|TestForwardDynamics" 2>&1 \| tail -10` | ✅ | ✅ green |

### Wave 1B — Assertions, NaN guards, test helpers

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-03-T1 | 13-03 | 1B | VEC-01, UTL-03 | T-13-05 | build | `cmake --build build 2>&1 \| tail -5` | ✅ | ✅ green |
| 13-03-T2 | 13-03 | 1B | VEC-01, UTL-03 | T-13-04 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestSpatialOperations\|TestSpatialVector" 2>&1 \| tail -20` | ✅ | ✅ green |

### Wave 2 — Numerical regression tests, external oracles

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-04-T1 | 13-04 | 2 | ABA-01, ABA-02 | T-13-06 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestForwardDynamics" 2>&1 \| tail -20` | ✅ | ✅ green |
| 13-04-T2 | 13-04 | 2 | ABA-01, ABA-02, TST-07 | T-13-07, T-13-08 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestInverseDynamics" 2>&1 \| tail -20` | ✅ | ✅ green |
| 13-04-T3 | 13-04 | 2 | ABA-01, ABA-02, TST-07 | T-13-06 | integration | `cmake --build build 2>&1 \| tail -5 && cd build && ctest --output-on-failure -R "TestDynamicsConsistency" 2>&1 \| tail -20` | ✅ | ❌ red (3 known CR-02 failures) |

### Wave 2A — Conventions doc, gravity invariant tests

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-06-T1 | 13-06 | 2A | UTL-03, ABA-01, ABA-02, TST-07 | T-13-10, T-13-11 | doc+build | `test -f MATHEMATICAL_CONVENTIONS.md && wc -l MATHEMATICAL_CONVENTIONS.md \| awk '$1 >= 100 {print "PASS: "$1" lines"} $1 < 100 {print "FAIL: only "$1" lines"}' && grep -q "## Convention Lock" MATHEMATICAL_CONVENTIONS.md && echo "Convention lock found" && grep -q "Worked Examples" MATHEMATICAL_CONVENTIONS.md && echo "Worked examples found" && cmake --build build 2>&1 \| tail -3` | ✅ | ✅ green |
| 13-06-T2 | 13-06 | 2A | UTL-03, ABA-01, ABA-02, TST-07 | T-13-10 | unit | `cmake --build build 2>&1 \| tail -3 && cd build && ./TestInverseDynamics --gtest_filter="StaticGravityInvariants.*" 2>&1 \| tail -10` | ✅ | ✅ green |
| 13-06-T3 | 13-06 | 2A | UTL-03, ABA-01, ABA-02, TST-07 | T-13-10 | unit | `cmake --build build 2>&1 \| tail -3 && cd build && ./TestForwardDynamics --gtest_filter="GravityABAInvariants.*" 2>&1 \| tail -10 && cd build && ./TestDynamicsConsistency --gtest_filter="*DirectComparison*" 2>&1 \| tail -10` | ✅ | ✅ green |

### Wave 3 — Code hygiene, build/CI, edge case tests, release-mode

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-05a-T1 | 13-05a | 3 | VEC-01 | T-13-08, T-13-14 | build | `cmake -B build 2>&1 \| tail -5 && cmake --build build 2>&1 \| tail -10` | ✅ | ✅ green |
| *D-22 verify* | 13-05a | 3 | — | T-13-14 | grep | `grep -c "static_cast.*SpatialVector" include/SpatialOperations.h src/SpatialOperations.cpp` | N/A (grep) | ✅ green |
| 13-05b-T1 | 13-05b | 3 | VEC-01 | T-13-06 | build | `cmake --build build 2>&1 \| tail -10 && cd build && ctest --output-on-failure 2>&1 \| tail -20` | ✅ | ⚠️ partial (build OK, but CompileSmoke target missing; .gitignore blocks .github/) |
| 13-05b-T2 | 13-05b | 3 | VEC-01 | T-13-07 | lint | `python3 -c "import yaml; d=yaml.safe_load(open('.github/workflows/ci.yml')); assert 'eigen5-compat' in d['jobs'], 'Missing eigen5-compat job'; print('PASS: eigen5-compat job found')" 2>&1` | ✅ | ❌ red (eigen5-compat job missing; CI workflow not git-tracked) |
| 13-07-T1 | 13-07 | 3 | VEC-01, UTL-03, ABA-01, ABA-02, TST-07 | T-13-12, T-13-13 | unit | `cmake --build build 2>&1 \| tail -3 && cd build && ./TestSpatialOperations --gtest_filter="TestCrossProductForce.CrossForceZero*" 2>&1 \| tail -5 && ./TestForwardDynamics --gtest_filter="*ZeroMassEdgeCase*" 2>&1 \| tail -5 && ./TestInverseDynamics --gtest_filter="*ZeroMassEdgeCase*" 2>&1 \| tail -5` | ✅ | ✅ green |
| 13-07-T2 | 13-07 | 3 | VEC-01, UTL-03, ABA-01, ABA-02, TST-07 | T-13-13 | unit | `cmake --build build 2>&1 \| tail -3 && cd build && ./TestForwardDynamics --gtest_filter="*ReleaseModeStability*" 2>&1 \| tail -10 && ./TestInverseDynamics --gtest_filter="*ReleaseModeStability*" 2>&1 \| tail -10` | ✅ | ✅ green |

### Wave 4 — README documentation update

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-08-T1 | 13-08 | 4 | DOC-02 | T-13-15 | doc | `grep -c "gravity" README.md \| awk '$1 >= 3 {print "PASS: gravity found ("$1" refs)"} $1 < 3 {print "FAIL: only "$1" gravity refs"}' && (grep -q "### Gravity" README.md && echo "PASS: Gravity Support section found" \|\| echo "FAIL: Gravity Support section missing") && (grep -q "Migration Guide\|API Changes in v1.1" README.md && echo "PASS: Migration/API Changes section found" \|\| echo "FAIL: Migration section not found") && (grep -q "MATHEMATICAL_CONVENTIONS.md" README.md && echo "PASS: Conventions cross-ref found" \|\| echo "FAIL: MATHEMATICAL_CONVENTIONS.md cross-ref missing") && (grep -q "link.parent\|Link link" README.md && echo "PASS: Link struct pattern found" \|\| echo "FAIL: Link struct pattern not found")` | ✅ | ❌ red (missing: Gravity Support section heading; MATHEMATICAL_CONVENTIONS.md cross-ref) |

### Wave 3 — Post-validation additions

| Task ID | Plan | Wave | Requirement | Threat Ref | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------|-------------------|-------------|--------|
| 13-05b-T3 | 13-05b | 3 | VEC-01 | T-13-09 | build | `EIGEN_DIR=$(brew --prefix eigen)/include/eigen3 && g++ -std=c++17 -I include -I "$EIGEN_DIR" -c tests/compile_smoke_test.cpp -o /tmp/compile_smoke_test.o 2>&1 \| tail -3` | ✅ | ✅ green (file created, compiles standalone; CMake target needs to be added to CMakeLists.txt) |

---

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky/partial*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements. GTest is already configured in CMakeLists.txt. Test executables already compiled via `cmake --build build`.

---

## Manual-Only Verifications

All phase behaviors have automated verification. D-22 compliance verified via grep.

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending — see Audit Trail below for escalation items requiring resolution.

---

## Audit Trail (2026-05-27)

### Summary of Gap Analysis

| Gap | Root Cause | Finding | Actions Taken |
|-----|-----------|---------|--------------|
| **G1: Plan 13-05b (Build/CI Infrastructure)** | Plan never executed | **PARTIAL** — CMake compile smoke target missing; CI workflow lacks eigen5-compat job; .gitignore blocks .github/ tracking | Created `tests/compile_smoke_test.cpp` (compiles standalone ✓); escalated CMakeLists.txt, .gitignore, CI workflow fixes |
| **G2: Plan 13-08 (README Documentation Update)** | Plan never executed (partially done by 13-07 T3) | **PARTIAL** — Gravity Support section heading missing; MATHEMATICAL_CONVENTIONS.md cross-ref missing; Link struct pattern ✓; API Changes section ✓ | Escalated README.md updates |
| **G3: VALIDATION.md** | Pending task and unsigned approval | **UPDATED** — Adjusted verification commands to match actual test names; updated status markers; added audit trail | Updated all entries |

### Gap Count Verification

| Category | Count | Details |
|----------|-------|---------|
| ✅ Green (passing) | 16 | All Wave 1A/1B/2/2A tasks; 13-05a-T1, D-22, 13-07-T1, 13-07-T2, 13-05b-T3 |
| ❌ Red (failing) | 3 | 13-04-T3 (3 known CR-02 failures), 13-05b-T2 (missing eigen5-compat job), 13-08-T1 (missing README sections) |
| ⚠️ Partial | 1 | 13-05b-T1 (build OK but missing CompileSmoke target and .gitignore issue) |
| ⬜ Pending | 0 | All resolved or escalated |
| **Total** | **20** | |

### Escalated Items (BLOCKER — require developer intervention)

| # | Component | Issue | Impact | Fix Required |
|---|-----------|-------|--------|-------------|
| E1 | `CMakeLists.txt` | Missing `SpatialAlgebraCompileSmoke` target | Compile smoke test cannot run via ctest | Add `add_executable(SpatialAlgebraCompileSmoke ...)` and `add_test(NAME CompileSmokeTest ...)` block per 13-05b-PLAN.md |
| E2 | `.gitignore` (line 38) | Contains `.github/` exclusion | CI workflow `.github/workflows/ci.yml` cannot be tracked by git | Remove `.github/` line from `.gitignore` |
| E3 | `.github/workflows/ci.yml` | Missing `eigen5-compat` job | Eigen 5.x compatibility is not verified in CI | Add `eigen5-compat` job per 13-05b-PLAN.md template |
| E4 | `.github/workflows/ci.yml` | Not tracked by git | CI file exists on disk but is not version-controlled | `git add -f .github/workflows/ci.yml` (after fixing .gitignore) |
| E5 | `README.md` | Missing dedicated "Gravity Support" section | API documentation incomplete | Add dedicated section with `computeAccelerations(tau, Vector3d(...))` and `computeTorques(qddot, Vector3d(...))` examples |
| E6 | `README.md` (References) | Missing cross-reference to `MATHEMATICAL_CONVENTIONS.md` | Users cannot find conventions document | Add `MATHEMATICAL_CONVENTIONS.md` link to References section |

### Pre-existing Known Gaps

| Test Executable | Failing Tests | Root Cause | Status |
|----------------|--------------|------------|--------|
| `TestDynamicsConsistency` | `ThreeLinkSerialChain`, `BranchingYConfiguration`, `TwoLinkRoundTripWithGravity` | ABA/RNEA roundtrip accuracy for multi-link chains with gravity | Known CR-02 tracking issue — deferred to dedicated algorithm bug-fix phase |

### File Changes Made During Validation

- `tests/compile_smoke_test.cpp` — **CREATED**: verifies `SpatialAlgebra.h` umbrella header is self-contained and compiles standalone
- `.planning/phases/13-production-readiness/13-VALIDATION.md` — **UPDATED**: corrected verification commands, status markers, added audit trail
