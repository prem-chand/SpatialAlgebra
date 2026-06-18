/**
 * @file bench_core_stub.cpp
 * @brief Core microbenchmarks — Plücker transforms, cross products, and inertia operations
 * @details Implements 18 individual benchmark functions covering the full spatial algebra
 *          API surface (per D-18). Each operation gets its own benchmark function with
 *          pre-constructed test objects and DoNotOptimize to prevent dead-code elimination.
 *
 *          Plücker transforms (8): transformMotion, transformForce, inverseTransformMotion,
 *          inverseTransformForce, tformRBI, invtformRBI, tformABI, invtformABI
 *
 *          Cross products (3): cross(MV, MV), cross(MV, FV), cross(FV, FV)
 *
 *          Inertia operations (7): RBI apply/add/scale, ABI apply/add/addRBI/scale
 */
#include <benchmark/benchmark.h>
#include "PluckerTransform.h"
#include "MotionVector.h"
#include "ForceVector.h"
#include "SpatialUtils.h"
#include "RigidBodyInertia.h"
#include "ArticulatedBodyInertia.h"
#include "LowerTriangular.h"

using namespace SpatialAlgebra;

// ── Factory helpers ──────────────────────────────────────────────────────

static PluckerTransform make_test_transform() {
    return PluckerTransform(
        Rotation(Eigen::AngleAxisd(0.5, Vector3d::UnitZ())),
        Vector3d(0.3, 0.1, 0.7)
    );
}

static MotionVector make_test_mv() {
    return MotionVector(Vector3d(0.2, 0.5, 0.8), Vector3d(0.1, 0.3, 0.6));
}

static ForceVector make_test_fv() {
    return ForceVector(Vector3d(0.4, 0.9, 0.2), Vector3d(0.7, 0.1, 0.5));
}

static RigidBodyInertia make_test_rbi() {
    return RigidBodyInertia(1.5, Vector3d(0.1, 0.2, 0.3), LowerTriangular::Identity(3));
}

static ArticulatedBodyInertia make_test_abi() {
    return ArticulatedBodyInertia(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );
}

// ── 8 Plücker transform benchmarks ──────────────────────────────────────

void BM_TransformMotion(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        MotionVector result = X.transformMotion(mv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_TransformForce(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ForceVector fv = make_test_fv();
    for (auto _ : state) {
        ForceVector result = X.transformForce(fv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_InverseTransformMotion(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        MotionVector result = X.inverseTransformMotion(mv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_InverseTransformForce(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ForceVector fv = make_test_fv();
    for (auto _ : state) {
        ForceVector result = X.inverseTransformForce(fv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_TformRBI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        RigidBodyInertia result = X.tformRBI(rbi);
        benchmark::DoNotOptimize(result);
    }
}

void BM_InvtformRBI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        RigidBodyInertia result = X.invtformRBI(rbi);
        benchmark::DoNotOptimize(result);
    }
}

void BM_TformABI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ArticulatedBodyInertia abi = make_test_abi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = X.tformABI(abi);
        benchmark::DoNotOptimize(result);
    }
}

void BM_InvtformABI(benchmark::State& state) {
    PluckerTransform X = make_test_transform();
    ArticulatedBodyInertia abi = make_test_abi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = X.invtformABI(abi);
        benchmark::DoNotOptimize(result);
    }
}

// ── 3 Cross product benchmarks ──────────────────────────────────────────

void BM_CrossMvMv(benchmark::State& state) {
    MotionVector mv1 = make_test_mv();
    MotionVector mv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));
    for (auto _ : state) {
        MotionVector result = cross(mv1, mv2);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossMvFv(benchmark::State& state) {
    MotionVector mv = make_test_mv();
    ForceVector fv = make_test_fv();
    for (auto _ : state) {
        ForceVector result = cross(mv, fv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_CrossFvFv(benchmark::State& state) {
    ForceVector fv1 = make_test_fv();
    ForceVector fv2(Vector3d(0.5, 0.3, 0.1), Vector3d(0.9, 0.7, 0.4));
    for (auto _ : state) {
        ForceVector result = cross(fv1, fv2);
        benchmark::DoNotOptimize(result);
    }
}

// ── 7 Inertia operation benchmarks ──────────────────────────────────────

void BM_RBIApply(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        ForceVector result = rbi.apply(mv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_RBIAdd(benchmark::State& state) {
    RigidBodyInertia rbi1 = make_test_rbi();
    RigidBodyInertia rbi2(2.0, Vector3d(0.05, 0.1, 0.15), LowerTriangular::Identity(3));
    for (auto _ : state) {
        RigidBodyInertia result = rbi1 + rbi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_RBIScale(benchmark::State& state) {
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        RigidBodyInertia result = rbi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIApply(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    MotionVector mv = make_test_mv();
    for (auto _ : state) {
        ForceVector result = abi.apply(mv);
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIAdd(benchmark::State& state) {
    ArticulatedBodyInertia abi1 = make_test_abi();
    ArticulatedBodyInertia abi2(
        LowerTriangular::Identity(3),
        Eigen::Matrix3d::Identity(),
        LowerTriangular::Identity(3)
    );
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi1 + abi2;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIAddRBI(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    RigidBodyInertia rbi = make_test_rbi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi + rbi;
        benchmark::DoNotOptimize(result);
    }
}

void BM_ABIScale(benchmark::State& state) {
    ArticulatedBodyInertia abi = make_test_abi();
    for (auto _ : state) {
        ArticulatedBodyInertia result = abi * 2.0;
        benchmark::DoNotOptimize(result);
    }
}
