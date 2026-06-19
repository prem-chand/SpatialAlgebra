/**
 * @file TestCompositeRigidBodyInertia.cpp
 * @brief Tests for the Composite Rigid Body Algorithm
 * @details Validates CRBA mass matrix against the existing RNEA-column method
 *          for serial chains with random configurations. Also tests edge cases:
 *          empty chain, single link, and known kinematic configurations.
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "CompositeRigidBodyInertia.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include "LowerTriangular.h"
#include "SpatialUtils.h"

using namespace Eigen;
using namespace SpatialAlgebra;

namespace {

// =========================================================================
// Helpers: build serial chain test fixtures
// =========================================================================

/**
 * @brief Build an N-link serial chain in the X-Y plane
 * @param n Number of links
 * @param linkLength Length of each link (m)
 * @param linkMass Mass of each link (kg)
 * @return Populated ForwardDynamics with kinematics and inertias set
 *
 * Each link is a REVOLUTE joint about Z at the end of the previous link.
 * Link i rotates about local Z. The joint axis is always Z.
 * Transform X_init places link at end of previous link along X.
 */
ForwardDynamics buildSerialChain(int n, double linkLength = 0.5, double linkMass = 1.0) {
    ForwardDynamics fd;

    for (int i = 0; i < n; ++i) {
        Link link;

        // Parent index: -1 for base, i-1 otherwise
        link.parent = (i == 0) ? -1 : (i - 1);

        // Parent-to-link transform: translate along X by linkLength then rotate about Z
        Matrix3d R = Matrix3d::Identity();
        Vector3d t(linkLength, 0.0, 0.0);
        link.X = PluckerTransform(Rotation(R), t);
        link.X0 = link.X;

        // Joint motion subspace: revolute about local Z
        link.S = MotionVector(Vector3d::UnitZ(), Vector3d::Zero());

        // Link inertia: point mass at link center
        Vector3d com(linkLength * 0.5, 0.0, 0.0);
        // Small rotational inertia about COM (thin rod about center, Y-axis = mL²/12)
        double I_yy = linkMass * linkLength * linkLength / 12.0;
        lt inertia(3);
        inertia(0, 0) = I_yy;  // I_xx ~ 0 for thin rod along X
        inertia(1, 1) = I_yy;  // I_yy
        inertia(2, 2) = 0.001; // I_zz small (point mass approximation)

        link.I = RigidBodyInertia(linkMass, com, inertia);

        link.type = JointType::REVOLUTE;
        link.jointAxis = Vector3d::UnitZ();
        link.q = 0.0;
        link.qdot = 0.0;
        link.qddot = 0.0;

        fd.links.push_back(link);
    }

    return fd;
}

/**
 * @brief Compute mass matrix using RNEA-column method (reference)
 * @param fd ForwardDynamics with links populated
 * @return n×n mass matrix
 *
 * This is the established method used in SpatialAlgebraAdapter::computeMassMatrix().
 * For each DOF j, set unit acceleration and compute torque column via RNEA.
 */
MatrixXd computeMassMatrixRNEA(ForwardDynamics& fd) {
    int n = static_cast<int>(fd.links.size());

    // Build an InverseDynamics solver from the FD link parameters
    InverseDynamics id;
    for (int i = 0; i < n; ++i) {
        InverseDynamicsLink link;
        const auto& src = fd.links[static_cast<size_t>(i)];
        link.parent = src.parent;
        link.X = src.X;
        link.I = src.I;
        link.S = src.S;
        link.q = src.q;
        link.qdot = 0.0;  // Zero velocity — mass matrix is position-only
        link.qddot = 0.0;
        link.type = src.type;
        link.jointAxis = src.jointAxis;
        link.X0 = src.X0;
        id.links.push_back(link);
    }

    MatrixXd H = MatrixXd::Zero(n, n);

    for (int j = 0; j < n; ++j) {
        VectorXd e_j = VectorXd::Zero(n);
        e_j(j) = 1.0;

        VectorXd tau_col = id.computeTorques(e_j, Vector3d::Zero());
        H.col(j) = tau_col;
    }

    return H;
}

/**
 * @brief Compute relative Frobenius norm difference between two matrices
 */
double relativeError(const MatrixXd& A, const MatrixXd& B) {
    return (A - B).norm() / (A.norm() + 1e-15);
}

} // anonymous namespace

// =========================================================================
// Tests
// =========================================================================

TEST(CompositeRigidBodyInertiaTest, EmptyChain) {
    ForwardDynamics fd;
    MatrixXd H = computeCompositeRigidBodyInertia(fd.links);
    EXPECT_EQ(H.rows(), 0);
    EXPECT_EQ(H.cols(), 0);
}

TEST(CompositeRigidBodyInertiaTest, SingleLink) {
    ForwardDynamics fd = buildSerialChain(1);
    fd.updateKinematics();

    MatrixXd H_crba = computeCompositeRigidBodyInertia(fd.links);
    MatrixXd H_rnea = computeMassMatrixRNEA(fd);

    EXPECT_EQ(H_crba.rows(), 1);
    EXPECT_EQ(H_crba.cols(), 1);
    EXPECT_GT(H_crba(0, 0), 0.0);  // Mass matrix is positive definite

    double err = relativeError(H_crba, H_rnea);
    EXPECT_LT(err, 1e-12) << "Single link mass matrix mismatch: " << err;
}

TEST(CompositeRigidBodyInertiaTest, TwoLinkChain) {
    ForwardDynamics fd = buildSerialChain(2);
    fd.updateKinematics();

    MatrixXd H_crba = computeCompositeRigidBodyInertia(fd.links);
    MatrixXd H_rnea = computeMassMatrixRNEA(fd);

    double err = relativeError(H_crba, H_rnea);
    EXPECT_LT(err, 1e-12) << "Two-link mass matrix mismatch: " << err;

    // Verify symmetry
    EXPECT_NEAR((H_crba - H_crba.transpose()).norm(), 0.0, 1e-15);
    EXPECT_NEAR(H_crba(0, 1), H_crba(1, 0), 1e-15);
}

TEST(CompositeRigidBodyInertiaTest, SixLinkChain) {
    ForwardDynamics fd = buildSerialChain(6);
    fd.updateKinematics();

    MatrixXd H_crba = computeCompositeRigidBodyInertia(fd.links);
    MatrixXd H_rnea = computeMassMatrixRNEA(fd);

    double err = relativeError(H_crba, H_rnea);
    EXPECT_LT(err, 1e-12) << "Six-link mass matrix mismatch: " << err;

    // Verify symmetry
    EXPECT_NEAR((H_crba - H_crba.transpose()).norm(), 0.0, 1e-14);
}

TEST(CompositeRigidBodyInertiaTest, NonZeroConfig) {
    // Test at a non-zero joint configuration
    ForwardDynamics fd = buildSerialChain(4);

    // Set arbitrary joint angles
    fd.links[0].q = 0.5;
    fd.links[1].q = -0.3;
    fd.links[2].q = 1.2;
    fd.links[3].q = 0.8;

    fd.updateKinematics();

    MatrixXd H_crba = computeCompositeRigidBodyInertia(fd.links);
    MatrixXd H_rnea = computeMassMatrixRNEA(fd);

    double err = relativeError(H_crba, H_rnea);
    EXPECT_LT(err, 1e-12) << "Non-zero config mass matrix mismatch: " << err;
}

TEST(CompositeRigidBodyInertiaTest, RandomConfigurations) {
    const int NUM_TESTS = 20;

    for (int trial = 0; trial < NUM_TESTS; ++trial) {
        ForwardDynamics fd = buildSerialChain(5);

        // Random joint positions in [-π, π]
        for (int i = 0; i < 5; ++i) {
            fd.links[static_cast<size_t>(i)].q =
                (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 * M_PI - M_PI;
        }

        fd.updateKinematics();

        MatrixXd H_crba = computeCompositeRigidBodyInertia(fd.links);
        MatrixXd H_rnea = computeMassMatrixRNEA(fd);

        double err = relativeError(H_crba, H_rnea);
        EXPECT_LT(err, 1e-12) << "Trial " << trial << ": mass matrix mismatch: " << err;
    }
}

TEST(CompositeRigidBodyInertiaTest, PositiveDefinite) {
    // The mass matrix should be positive definite at any non-singular config
    ForwardDynamics fd = buildSerialChain(6);

    // Random configuration
    for (size_t i = 0; i < fd.links.size(); ++i) {
        fd.links[i].q = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 * M_PI - M_PI;
    }

    fd.updateKinematics();

    MatrixXd H = computeCompositeRigidBodyInertia(fd.links);

    // Check all eigenvalues are positive
    SelfAdjointEigenSolver<MatrixXd> solver(H);
    const VectorXd& eigenvalues = solver.eigenvalues();
    for (int i = 0; i < eigenvalues.size(); ++i) {
        EXPECT_GT(eigenvalues(i), 0.0)
            << "Eigenvalue " << i << " = " << eigenvalues(i) << " should be positive";
    }
}

// =========================================================================
// Symmetry and structural properties
// =========================================================================

TEST(CompositeRigidBodyInertiaTest, Symmetry) {
    ForwardDynamics fd = buildSerialChain(8);
    for (size_t i = 0; i < fd.links.size(); ++i) {
        fd.links[i].q = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 * M_PI - M_PI;
    }
    fd.updateKinematics();

    MatrixXd H = computeCompositeRigidBodyInertia(fd.links);

    // Frobenius norm of (H - H^T) should be near zero
    double asym = (H - H.transpose()).norm();
    EXPECT_LT(asym, 1e-14) << "Mass matrix is not symmetric: ||H - H^T|| = " << asym;
}

TEST(CompositeRigidBodyInertiaTest, VariesWithConfiguration) {
    // Create a 2-link chain with a prismatic joint to verify H(q) changes.
    // Prismatic joints change COM distances from joint axes, making H
    // clearly configuration-dependent.
    ForwardDynamics fd;
    for (int i = 0; i < 2; ++i) {
        Link link;
        link.parent = (i == 0) ? -1 : (i - 1);
        Matrix3d R = Matrix3d::Identity();
        Vector3d t(0.5, 0.0, 0.0);
        link.X = PluckerTransform(Rotation(R), t);
        link.X0 = link.X;
        // First joint revolute Z, second prismatic X
        if (i == 0) {
            link.S = MotionVector(Vector3d::UnitZ(), Vector3d::Zero());
            link.type = JointType::REVOLUTE;
            link.jointAxis = Vector3d::UnitZ();
        } else {
            link.S = MotionVector(Vector3d::Zero(), Vector3d::UnitX());
            link.type = JointType::PRISMATIC;
            link.jointAxis = Vector3d::UnitX();
        }
        Vector3d com(0.25, 0.0, 0.0);
        Matrix3d J = Matrix3d::Zero();
        J(0,0) = 0.001; J(1,1) = 0.001; J(2,2) = 0.001;
        link.I = RigidBodyInertia(1.0, com, LowerTriangular::fromFullMatrix(J));
        link.q = 0.0; link.qdot = 0.0; link.qddot = 0.0;
        fd.links.push_back(link);
    }
    fd.links[1].q = 0.0;
    fd.updateKinematics();
    MatrixXd H0 = computeCompositeRigidBodyInertia(fd.links);

    fd.links[1].q = 2.0;  // Extend prismatic joint
    fd.updateKinematics();
    MatrixXd H1 = computeCompositeRigidBodyInertia(fd.links);

    double diff = (H0 - H1).norm();
    EXPECT_GT(diff, 1.0) << "Mass matrix should change with prismatic joint: diff = " << diff;
}
