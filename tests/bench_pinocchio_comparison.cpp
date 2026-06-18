/**
 * @file bench_pinocchio_comparison.cpp
 * @brief Cross-validation benchmark comparing SpatialAlgebra vs Pinocchio
 * @details Builds identical kinematic chains through both the SpatialAlgebra
 *          adapter and the Pinocchio adapter, then compares computed dynamics
 *          quantities per-joint with relative error thresholds.
 *
 *          Test models covered (PCC-03):
 *          1. Single-link Z-revolute (mass=1, COM=origin)
 *          2. Two-link serial Z-revolute (1m X spacing)
 *          3. Three-link serial Z-revolute (1m X spacing)
 *          4. Single-link gravity test (X-axis, COM offset)
 *          5. Branching Y configuration
 *
 *          Quantities compared:
 *          - computeTorques (RNEA)
 *          - computeAccelerations (ABA)
 *          - computeMassMatrix
 *          - computeGravityTorques
 *          - getJointTransform
 *          - getLinkCOM
 *
 *          Tolerance threshold: per-joint relative error < 1e-6 (PCC-04 SC-4).
 *
 * @see RobotSolver, SpatialAlgebraAdapter, PinocchioAdapter
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <random>

#include "robot_model.h"
#include "robot_solver.h"
#include "sa_adapter.h"
#include "pinocchio_adapter.h"

// Chain model factories
#include "chains/inverse_dynamics.h"
#include "chains/forward_dynamics.h"
#include "chains/consistency.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace test_models {
namespace {

using Vector3d = Eigen::Vector3d;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

// =========================================================================
// Tolerance
// =========================================================================

/// Relative error threshold for per-joint comparisons
constexpr double kRelTol = 1e-6;

/// Absolute error threshold for near-zero values
constexpr double kAbsTol = 1e-10;

// =========================================================================
// Error comparison helpers
// =========================================================================

/**
 * @brief Compute relative error between computed and reference values
 * @param computed Value from Pinocchio adapter
 * @param reference Value from SpatialAlgebra adapter
 * @return |computed - reference| / max(|reference|, 1)
 */
double relError(double computed, double reference) {
    double denom = std::max(std::abs(reference), 1.0);
    return std::abs(computed - reference) / denom;
}

/**
 * @brief Check vector relative error
 * @param computed Computed vector
 * @param reference Reference vector
 * @param[out] maxErr Maximum relative error found
 * @param[out] meanErr Mean relative error
 * @return true if all components within tolerance
 */
bool checkVectorRelError(const VectorXd& computed,
                          const VectorXd& reference,
                          double& maxErr,
                          double& meanErr) {
    if (computed.size() != reference.size()) return false;

    maxErr = 0.0;
    double sumErr = 0.0;
    int n = computed.size();

    for (int i = 0; i < n; ++i) {
        double err = relError(computed[i], reference[i]);
        maxErr = std::max(maxErr, err);
        sumErr += err;
    }
    meanErr = sumErr / static_cast<double>(n);
    return maxErr < kRelTol;
}

/**
 * @brief Check matrix relative error
 */
bool checkMatrixRelError(const MatrixXd& computed,
                          const MatrixXd& reference,
                          double& maxErr,
                          double& meanErr) {
    if (computed.rows() != reference.rows() ||
        computed.cols() != reference.cols()) return false;

    maxErr = 0.0;
    double sumErr = 0.0;
    int count = computed.rows() * computed.cols();

    for (int r = 0; r < computed.rows(); ++r) {
        for (int c = 0; c < computed.cols(); ++c) {
            double err = relError(computed(r, c), reference(r, c));
            maxErr = std::max(maxErr, err);
            sumErr += err;
        }
    }
    meanErr = sumErr / static_cast<double>(count);
    return maxErr < kRelTol;
}

/**
 * @brief Check 4×4 transform relative error
 */
bool checkTransformRelError(const Eigen::Matrix4d& computed,
                             const Eigen::Matrix4d& reference,
                             double& maxErr,
                             double& meanErr) {
    maxErr = 0.0;
    double sumErr = 0.0;
    int count = 0;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            // Skip bottom row which is always [0,0,0,1]
            if (r == 3) continue;
            double err = relError(computed(r, c), reference(r, c));
            maxErr = std::max(maxErr, err);
            sumErr += err;
        }
        count += 4;
    }
    meanErr = sumErr / static_cast<double>(count - 4); // skip bottom row
    return maxErr < kRelTol;
}

// =========================================================================
// Random state generation
// =========================================================================

VectorXd randomVector(int n, double minVal, double maxVal, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(minVal, maxVal);
    VectorXd v(n);
    for (int i = 0; i < n; ++i) v[i] = dist(rng);
    return v;
}

VectorXd randomPositions(int n, unsigned seed = 42) {
    return randomVector(n, -M_PI, M_PI, seed);
}

VectorXd randomVelocities(int n, unsigned seed = 123) {
    return randomVector(n, -5.0, 5.0, seed);
}

VectorXd randomAccelerations(int n, unsigned seed = 456) {
    return randomVector(n, -10.0, 10.0, seed);
}

VectorXd randomTorques(int n, unsigned seed = 789) {
    return randomVector(n, -50.0, 50.0, seed);
}

// =========================================================================
// Test fixture — runs all checks for a given model
// =========================================================================

struct ComparisonResult {
    std::string modelName;
    int dof;
    bool pass;

    // Per-check results
    struct CheckResult {
        std::string name;
        bool pass;
        double maxRelErr;
        double meanRelErr;
    };
    std::vector<CheckResult> checks;

    void print() const {
        std::cout << "  " << (pass ? "PASS" : "FAIL") << " | "
                  << std::setw(35) << std::left << modelName
                  << " | DOF=" << dof << std::endl;
        for (const auto& c : checks) {
            std::cout << "       " << (c.pass ? "✓" : "✗") << "  "
                      << std::setw(28) << std::left << c.name
                      << " max_err=" << std::scientific << std::setprecision(2)
                      << c.maxRelErr
                      << " mean_err=" << c.meanRelErr
                      << std::endl;
        }
    }
};

ComparisonResult compareAdapters(
    const std::string& name,
    std::unique_ptr<RobotSolver> saSolver,
    std::unique_ptr<RobotSolver> pinSolver) {

    ComparisonResult result;
    result.modelName = name;
    result.dof = saSolver->getDOF();
    result.pass = true;

    // Generate random joint state
    int n = result.dof;
    VectorXd q = randomPositions(n, 42);
    VectorXd qdot = randomVelocities(n, 123);
    VectorXd qddot = randomAccelerations(n, 456);
    VectorXd tau = randomTorques(n, 789);
    Vector3d gravity(0.0, 0.0, -9.81);

    // Set state on both adapters
    saSolver->setState(q, qdot);
    pinSolver->setState(q, qdot);

    // --- 1. computeTorques ---
    {
        VectorXd tauSa = saSolver->computeTorques(qddot, gravity);
        VectorXd tauPin = pinSolver->computeTorques(qddot, gravity);
        double maxErr, meanErr;
        bool ok = checkVectorRelError(tauPin, tauSa, maxErr, meanErr);
        result.checks.push_back({"computeTorques", ok, maxErr, meanErr});
        result.pass &= ok;
    }

    // --- 2. computeAccelerations ---
    {
        VectorXd qddotSa = saSolver->computeAccelerations(tau, gravity);
        VectorXd qddotPin = pinSolver->computeAccelerations(tau, gravity);
        double maxErr, meanErr;
        bool ok = checkVectorRelError(qddotPin, qddotSa, maxErr, meanErr);
        result.checks.push_back({"computeAccelerations", ok, maxErr, meanErr});
        result.pass &= ok;
    }

    // --- 3. computeMassMatrix ---
    {
        MatrixXd MSa = saSolver->computeMassMatrix();
        MatrixXd MPin = pinSolver->computeMassMatrix();
        double maxErr, meanErr;
        bool ok = checkMatrixRelError(MPin, MSa, maxErr, meanErr);
        result.checks.push_back({"computeMassMatrix", ok, maxErr, meanErr});
        result.pass &= ok;
    }

    // --- 4. computeGravityTorques ---
    {
        VectorXd gSa = saSolver->computeGravityTorques(gravity);
        VectorXd gPin = pinSolver->computeGravityTorques(gravity);
        double maxErr, meanErr;
        bool ok = checkVectorRelError(gPin, gSa, maxErr, meanErr);
        result.checks.push_back({"computeGravityTorques", ok, maxErr, meanErr});
        result.pass &= ok;
    }

    // --- 5. getJointTransform (each link) ---
    {
        double overallMaxErr = 0.0;
        double overallMeanErr = 0.0;
        bool ok = true;
        for (int i = 0; i < n; ++i) {
            Eigen::Matrix4d TSa = saSolver->getJointTransform(i);
            Eigen::Matrix4d TPin = pinSolver->getJointTransform(i);
            double maxErr, meanErr;
            bool checkOk = checkTransformRelError(TPin, TSa, maxErr, meanErr);
            ok &= checkOk;
            overallMaxErr = std::max(overallMaxErr, maxErr);
            overallMeanErr += meanErr;
        }
        overallMeanErr /= static_cast<double>(n);
        result.checks.push_back({"getJointTransform", ok, overallMaxErr, overallMeanErr});
        result.pass &= ok;
    }

    // --- 6. getLinkCOM (each link) ---
    {
        double overallMaxErr = 0.0;
        double overallMeanErr = 0.0;
        bool ok = true;
        for (int i = 0; i < n; ++i) {
            Vector3d comSa = saSolver->getLinkCOM(i);
            Vector3d comPin = pinSolver->getLinkCOM(i);
            double maxErr = 0.0;
            double meanErr = 0.0;
            // Manual check for Vector3d
            for (int j = 0; j < 3; ++j) {
                double err = relError(comPin[j], comSa[j]);
                maxErr = std::max(maxErr, err);
                meanErr += err;
            }
            meanErr /= 3.0;
            ok &= (maxErr < kRelTol);
            overallMaxErr = std::max(overallMaxErr, maxErr);
            overallMeanErr += meanErr;
        }
        overallMeanErr /= static_cast<double>(n);
        result.checks.push_back({"getLinkCOM", ok, overallMaxErr, overallMeanErr});
        result.pass &= ok;
    }

    return result;
}

// =========================================================================
// GTest integration
// =========================================================================

class PinocchioComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Empty — models constructed per test
    }

    /**
     * @brief Run comparison for a model and record GTest assertions
     */
    ComparisonResult runComparison(
        const std::string& name,
        std::unique_ptr<RobotSolver> saSolver,
        std::unique_ptr<RobotSolver> pinSolver) {

        ComparisonResult result = compareAdapters(
            name, std::move(saSolver), std::move(pinSolver));

        result.print();

        // GTest assertions for each check
        for (const auto& c : result.checks) {
            EXPECT_TRUE(c.pass)
                << name << " :: " << c.name
                << " — max_rel_err=" << c.maxRelErr
                << " (> " << kRelTol << ")";
        }

        return result;
    }
};

// =========================================================================
// Individual test cases
// =========================================================================

TEST_F(PinocchioComparisonTest, SingleLinkChain) {
    auto model = makeSingleLinkChain();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("SingleLinkChain",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, TwoLinkSerialChain) {
    auto model = makeTwoLinkSerialChain();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("TwoLinkSerialChain",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, ThreeLinkSerialChain) {
    auto model = makeThreeLinkSerialChain();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("ThreeLinkSerialChain",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, SingleLinkXAxisWithCOMOffset) {
    auto model = makeSingleLinkXAxisWithCOMOffset();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("SingleLinkXAxisCOMOffset",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, BranchingYConfiguration) {
    auto model = makeBranchingYConfiguration();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("BranchingYConfiguration",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, TwoLinkWithCOMOffset) {
    auto model = makeTwoLinkWithCOMOffset();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("TwoLinkWithCOMOffset",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, ThreeLinkSerialNonZeroCOM) {
    auto model = makeThreeLinkSerialChainNonZeroCOM();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("ThreeLinkSerialNonZeroCOM",
                   std::move(saPin.first), std::move(saPin.second));
}

TEST_F(PinocchioComparisonTest, TwoLinkWithGravity) {
    auto model = makeTwoLinkSerialChain();
    auto saPin = std::make_pair(
        SpatialAlgebraAdapter::build(model),
        PinocchioAdapter::build(model)
    );
    runComparison("TwoLinkWithGravity",
                   std::move(saPin.first), std::move(saPin.second));
}

// =========================================================================
// Entry point
// =========================================================================

} // anonymous namespace
} // namespace test_models

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
