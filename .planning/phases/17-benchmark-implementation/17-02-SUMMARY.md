# Plan 17-02 Summary: Core microbenchmarks

## Status: Complete

### Changes
- **benchmarks/core/bench_core_stub.cpp**: Replaced 2 stub functions with 18 real benchmark functions plus 5 factory helpers

### Benchmark functions implemented
- 8 Plücker transforms: BM_TransformMotion, BM_TransformForce, BM_InverseTransformMotion, BM_InverseTransformForce, BM_TformRBI, BM_InvtformRBI, BM_TformABI, BM_InvtformABI
- 3 cross products: BM_CrossMvMv, BM_CrossMvFv, BM_CrossFvFv (no cross(FV,MV) in API)
- 7 inertia operations: BM_RBIApply, BM_RBIAdd, BM_RBIScale, BM_ABIApply, BM_ABIAdd, BM_ABIAddRBI, BM_ABIScale

### Note
D-18 specified 19 but only 18 implemented — `cross(ForceVector, MotionVector)` overload does not exist in SpatialUtils.h

### Verification
- 18 `void BM_` function definitions exist
- 0 occurrences of `BM_PluckerTransform` or `BM_CrossProduct` stubs
- 0 occurrences of `int nDOF` in the file
- 18 DoNotOptimize calls (one per benchmark; BM_CrossFvMv removed reduces count from plan's 19)
- `cmake --build build_bench` succeeds
- All filtered benchmarks produce valid timing output
