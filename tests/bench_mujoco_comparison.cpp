#include <gtest/gtest.h>
#include "mujoco_adapter.h"
#include "sa_adapter.h"
#include "mjcf_parser.h"

#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace Eigen;
using namespace test_models;

namespace {

struct ModelSpec {
    std::string name;
    std::string path;
};

const std::vector<ModelSpec> MODELS = {
    {"Go1",  std::string(MODELS_ROOT) + "/menagerie/unitree_go1/go1.xml"},
    {"Go2",  std::string(MODELS_ROOT) + "/menagerie/unitree_go2/go2.xml"},
    {"H1",   std::string(MODELS_ROOT) + "/menagerie/unitree_h1/h1.xml"},
    {"G1",   std::string(MODELS_ROOT) + "/menagerie/unitree_g1/g1.xml"},
};

// Random state within joint limits
VectorXd randomQ(const RobotModel& model, std::mt19937& rng) {
    int n = model.getDOF();
    VectorXd q(n);
    // Use reasonable joint range: [-π, π] for revolutes, [-0.5, 0.5] for prismatics
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    for (int i = 0; i < n; i++) q(i) = dist(rng);
    return q;
}

VectorXd randomV(int n, std::mt19937& rng) {
    VectorXd v(n);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < n; i++) v(i) = dist(rng);
    return v;
}

} // anonymous namespace

// =========================================================================
// Accuracy benchmark — compare SA vs MuJoCo across 100 random configs
// =========================================================================
TEST(MujocoAccuracy, AllModels) {
    const int TRIALS = 100;
    const Vector3d GRAVITY(0, 0, -9.81);

    std::mt19937 rng(42);

    for (const auto& spec : MODELS) {
        std::cout << "\n=== " << spec.name << " ===\n";

        RobotModel model;
        EXPECT_NO_THROW(model = buildRobotModelFromMJCF(spec.path));
        int n = model.getDOF();
        std::cout << "  DOF: " << n << "\n";

        auto sa = SpatialAlgebraAdapter::build(model);
        auto mj = MujocoAdapter::build(model);

        double tau_max_err = 0.0, tau_sum_err = 0.0;
        double qddot_max_err = 0.0, qddot_sum_err = 0.0;
        double grav_max_err = 0.0, grav_sum_err = 0.0;
        double mass_err = 0.0;

        double sa_time_total = 0.0, mj_time_total = 0.0;

        for (int trial = 0; trial < TRIALS; trial++) {
            VectorXd q = randomQ(model, rng);
            VectorXd qdot = randomV(n, rng);
            VectorXd qddot = randomV(n, rng);

            // ---- Inverse dynamics ----
            auto t0 = std::chrono::high_resolution_clock::now();
            sa->setState(q, qdot);
            VectorXd tau_sa = sa->computeTorques(qddot, GRAVITY);
            auto t1 = std::chrono::high_resolution_clock::now();
            mj->setState(q, qdot);
            VectorXd tau_mj = mj->computeTorques(qddot, GRAVITY);
            auto t2 = std::chrono::high_resolution_clock::now();

            sa_time_total += std::chrono::duration<double, std::micro>(t1 - t0).count();
            mj_time_total += std::chrono::duration<double, std::micro>(t2 - t1).count();

            double tau_err = (tau_sa - tau_mj).norm() / (tau_sa.norm() + 1e-10);
            tau_max_err = std::max(tau_max_err, tau_err);
            tau_sum_err += tau_err;

            // ---- Forward dynamics ----
            sa->setState(q, qdot);
            VectorXd qddot_sa = sa->computeAccelerations(tau_sa, GRAVITY);
            mj->setState(q, qdot);
            VectorXd qddot_mj = mj->computeAccelerations(tau_mj, GRAVITY);

            double qddot_err = (qddot_sa - qddot_mj).norm() / (qddot_sa.norm() + 1e-10);
            qddot_max_err = std::max(qddot_max_err, qddot_err);
            qddot_sum_err += qddot_err;

            // ---- Gravity torques ----
            VectorXd g_sa = sa->computeGravityTorques(GRAVITY);
            VectorXd g_mj = mj->computeGravityTorques(GRAVITY);
            double g_err = (g_sa - g_mj).norm() / (g_sa.norm() + 1e-10);
            grav_max_err = std::max(grav_max_err, g_err);
            grav_sum_err += g_err;

            // ---- Mass matrix (once) ----
            if (trial == 0) {
                MatrixXd H_sa = sa->computeMassMatrix();
                MatrixXd H_mj = mj->computeMassMatrix();
                mass_err = (H_sa - H_mj).norm() / H_sa.norm();
            }
        }

        double tau_mean = tau_sum_err / TRIALS;
        double qddot_mean = qddot_sum_err / TRIALS;
        double grav_mean = grav_sum_err / TRIALS;

        std::cout << "  Inverse dynamics (tau):\n"
                  << "    max_rel_err = " << tau_max_err << "\n"
                  << "    mean_rel_err = " << tau_mean << "\n";
        std::cout << "  Forward dynamics (qddot):\n"
                  << "    max_rel_err = " << qddot_max_err << "\n"
                  << "    mean_rel_err = " << qddot_mean << "\n";
        std::cout << "  Mass matrix:\n"
                  << "    rel_err = " << mass_err << "\n";
        std::cout << "  Gravity torques:\n"
                  << "    max_rel_err = " << grav_max_err << "\n"
                  << "    mean_rel_err = " << grav_mean << "\n";
        std::cout << "  Avg computeTorques time:\n"
                  << "    SA = " << (sa_time_total / TRIALS) << " μs\n"
                  << "    MuJoCo = " << (mj_time_total / TRIALS) << " μs\n";

        const double errorLimit = 1e-6;
        EXPECT_LT(tau_max_err, errorLimit);
        EXPECT_LT(qddot_max_err, errorLimit);
        EXPECT_LT(grav_max_err, errorLimit);
        EXPECT_LT(mass_err, errorLimit);
    }
}

// =========================================================================
// Efficiency benchmark — detailed timing
// =========================================================================
TEST(MujocoEfficiency, Timing) {
    const Vector3d GRAVITY(0, 0, -9.81);
    std::mt19937 rng(99);

    std::cout << "\n=== Efficiency Benchmark ===\n";
    std::cout << "Model       DOF  RNEA_SA(μs)  RNEA_MJ(μs)  ABA_SA(μs)  ABA_MJ(μs)  MassMat(μs)\n";

    for (const auto& spec : MODELS) {
        auto model = buildRobotModelFromMJCF(spec.path);
        int n = model.getDOF();
        auto sa = SpatialAlgebraAdapter::build(model);
        auto mj = MujocoAdapter::build(model);

        VectorXd q = randomQ(model, rng);
        VectorXd qdot = randomV(n, rng);
        VectorXd qddot = randomV(n, rng);

        const int WARMUP = 100;
        const int TRIALS = 1000;

        // Warmup
        for (int i = 0; i < WARMUP; i++) {
            sa->setState(q, qdot); sa->computeTorques(qddot, GRAVITY);
            mj->setState(q, qdot); mj->computeTorques(qddot, GRAVITY);
        }

        // SA RNEA
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            sa->setState(q, qdot);
            sa->computeTorques(qddot, GRAVITY);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double sa_rnea = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // MuJoCo RNEA
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            mj->setState(q, qdot);
            mj->computeTorques(qddot, GRAVITY);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double mj_rnea = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // SA ABA
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            sa->setState(q, qdot);
            sa->computeAccelerations(qddot, GRAVITY);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double sa_aba = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // MuJoCo ABA (mj_forward)
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            mj->setState(q, qdot);
            mj->computeAccelerations(qddot, GRAVITY);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double mj_aba = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // SA Mass matrix (RNEA-column method, O(n²))
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            sa->setState(q, qdot);
            sa->computeMassMatrix();
        }
        t1 = std::chrono::high_resolution_clock::now();
        double sa_mass = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        printf("%-12s %3d  %10.1f  %10.1f  %10.1f  %10.1f  %10.1f\n",
               spec.name.c_str(), n, sa_rnea, mj_rnea, sa_aba, mj_aba, sa_mass);
    }
}

// =========================================================================
// Scaling benchmark — synthetic N-link chains
// =========================================================================
TEST(MujocoEfficiency, Scaling) {
    const Vector3d GRAVITY(0, 0, -9.81);
    std::mt19937 rng(123);

    std::cout << "\n=== Scaling Benchmark ===\n";
    std::cout << "DOF  RNEA_SA(μs)  RNEA_MJ(μs)  MassMat_SA(μs)\n";

    for (int dof : {6, 12, 24, 48}) {
        // Build N-link serial chain
        RobotModel model;
        for (int i = 0; i < dof; i++) {
            JointSpec js;
            js.parent = (i == 0) ? -1 : (i - 1);
            js.parentToJoint = Matrix4d::Identity();
            js.parentToJoint(0, 3) = 0.5;  // 0.5m link length
            js.jointAxis = Vector3d::UnitZ();
            js.type = JointType::REVOLUTE;
            js.mass = 1.0;
            js.com = Vector3d(0.25, 0, 0);  // COM at midpoint
            js.inertia = Matrix3d::Identity() * 0.01;  // small rotational inertia
            model.joints.push_back(js);
        }

        auto sa = SpatialAlgebraAdapter::build(model);
        auto mj = MujocoAdapter::build(model);

        VectorXd q = VectorXd::Zero(dof);
        VectorXd qdot = VectorXd::Zero(dof);
        VectorXd qddot = VectorXd::Ones(dof);

        const int WARMUP = 50;
        const int TRIALS = 500;

        // Warmup
        for (int i = 0; i < WARMUP; i++) {
            sa->setState(q, qdot); sa->computeTorques(qddot, GRAVITY);
            mj->setState(q, qdot); mj->computeTorques(qddot, GRAVITY);
        }

        // SA RNEA
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            sa->setState(q, qdot);
            sa->computeTorques(qddot, GRAVITY);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double sa_rnea = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // MuJoCo RNEA
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            mj->setState(q, qdot);
            mj->computeTorques(qddot, GRAVITY);
        }
        t1 = std::chrono::high_resolution_clock::now();
        double mj_rnea = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        // SA Mass matrix (O(n²) — warmup not needed separately)
        t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < TRIALS; i++) {
            sa->setState(q, qdot);
            sa->computeMassMatrix();
        }
        t1 = std::chrono::high_resolution_clock::now();
        double sa_mass = std::chrono::duration<double, std::micro>(t1 - t0).count() / TRIALS;

        printf("%3d  %11.1f  %11.1f  %14.1f\n",
               dof, sa_rnea, mj_rnea, sa_mass);
    }
}
