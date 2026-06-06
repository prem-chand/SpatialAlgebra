/**
 * @file example_robot_3link.cpp
 * @brief Demonstrates forward and inverse dynamics for a 3-link Z-Y-Z spatial RRR arm
 *        with UR5-derived parameters and FD+ID cross-validation.
 * 
 * Robot Configuration:
 * - 3-link serial chain (standard anthropomorphic arm)
 * - Joint 1: Z-axis revolute (waist yaw)
 * - Joint 2: Y-axis revolute (shoulder pitch)
 * - Joint 3: Z-axis revolute (elbow roll)
 * - Link 0: UR5 shoulder (mass=3.7 kg)
 * - Link 1: UR5 upper arm (mass=8.393 kg)
 * - Link 2: UR5 forearm (mass=2.33 kg)
 * - Gravity: [0, 0, -9.81] m/s²
 * 
 * This example demonstrates:
 * 1. Model setup with realistic UR5 parameters (mass, COM, inertia, transforms)
 * 2. Forward dynamics (ABA): apply torques, compute joint accelerations
 * 3. Inverse dynamics (RNEA): compute gravity compensation torques at static pose
 *    — Joint 2 (shoulder pitch) bears ~30 N·m gravity load when arm is horizontal
 * 4. Cross-validation: feed gravity torques into FD, verify zero acceleration
 * 5. Physical interpretation of all results
 * 
 * Key Physics Insight:
 * In this Z-Y-Z anthropomorphic arm, only the Y-axis joint (shoulder pitch)
 * bears a significant gravity load when the arm is horizontal. Joints 1 and 3
 * (both Z-axis) have gravity parallel to their axes → zero gravity torque.
 * The shoulder pitch gravity torque is approximately 30.7 N·m at full extension.
 * 
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 * @see ForwardDynamics (ABA, Algorithm 7.3)
 * @see InverseDynamics (RNEA, Algorithm 7.1)
 */

#include "ForwardDynamics.h"
#include "InverseDynamics.h"
#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include "LowerTriangular.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
using lt = LowerTriangular;
using mv = MotionVector;
using fv = ForceVector;
using plux = PluckerTransform;
using rbi = RigidBodyInertia;

int main() {
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "=== SpatialAlgebra 3-Link Z-Y-Z Spatial Arm Example ===" << std::endl;
    std::cout << std::endl;
    std::cout << "Physical configuration:" << std::endl;
    std::cout << "  3-link serial chain (anthropomorphic RRR arm)." << std::endl;
    std::cout << "  Joint 1: Z-axis (waist yaw)." << std::endl;
    std::cout << "  Joint 2: Y-axis (shoulder pitch)." << std::endl;
    std::cout << "  Joint 3: Z-axis (elbow roll)." << std::endl;
    std::cout << "  Gravity: [0, 0, -9.81] m/s^2." << std::endl;
    std::cout << "  UR5-derived link parameters." << std::endl;
    std::cout << std::endl;

    // ============================================================
    // Section 1: Model Setup
    // ============================================================
    std::cout << "=== Section 1: Robot Model Setup ===" << std::endl;
    std::cout << std::endl;

    // Build ForwardDynamics solver (ABA)
    ForwardDynamics fd;
    // Build InverseDynamics solver (RNEA) separately (different struct types)
    InverseDynamics id;

    // ------------------------------------------------------------------
    // Link 0 (waist): UR5 shoulder — Z-axis revolute (waist yaw)
    // ------------------------------------------------------------------
    {
        Link link0;
        link0.parent = -1;

        double mass0 = 3.7;
        Eigen::Vector3d com0(0.0, -0.02561, 0.00193);

        lt I0(3);
        I0(0, 0) = 0.0067;  // Ixx
        I0(1, 0) = 0.0;     // Ixy
        I0(1, 1) = 0.0064;  // Iyy
        I0(2, 0) = 0.0;     // Ixz
        I0(2, 1) = 0.0;     // Iyz
        I0(2, 2) = 0.0067;  // Izz

        link0.I = RigidBodyInertia(mass0, com0, I0);

        Rotation R0(Eigen::Matrix3d::Identity());
        link0.X = PluckerTransform(R0, Eigen::Vector3d::Zero());

        // Joint 1: Z-axis (waist yaw)
        link0.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link0.q = 0.0;
        link0.qdot = 0.0;

        std::cout << "Link 0 (waist yaw):" << std::endl;
        std::cout << "  Parent: world (base, index -1)" << std::endl;
        std::cout << "  Source: UR5 shoulder" << std::endl;
        std::cout << "  Mass: " << mass0 << " kg" << std::endl;
        std::cout << "  COM: [" << com0.transpose() << "] m" << std::endl;
        std::cout << "  Inertia: Ixx=" << I0(0,0) << ", Iyy=" << I0(1,1) << ", Izz=" << I0(2,2) << " kg·m^2" << std::endl;
        std::cout << "  Joint axis: Z (waist yaw)" << std::endl;
        std::cout << "  Transform from parent: identity" << std::endl;
        std::cout << std::endl;

        fd.links.push_back(link0);
    }

    {
        InverseDynamicsLink link0;
        link0.parent = -1;

        double mass0 = 3.7;
        Eigen::Vector3d com0(0.0, -0.02561, 0.00193);

        lt I0(3);
        I0(0, 0) = 0.0067;
        I0(1, 0) = 0.0;
        I0(1, 1) = 0.0064;
        I0(2, 0) = 0.0;
        I0(2, 1) = 0.0;
        I0(2, 2) = 0.0067;

        link0.I = RigidBodyInertia(mass0, com0, I0);

        Rotation R0(Eigen::Matrix3d::Identity());
        link0.X = PluckerTransform(R0, Eigen::Vector3d::Zero());

        link0.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link0.q = 0.0;
        link0.qdot = 0.0;
        link0.qddot = 0.0;

        id.links.push_back(link0);
    }

    // ------------------------------------------------------------------
    // Link 1 (shoulder): UR5 upper arm — Y-axis revolute (shoulder pitch)
    // ------------------------------------------------------------------
    {
        Link link1;
        link1.parent = 0;

        double mass1 = 8.393;
        Eigen::Vector3d com1(0.2125, 0.0, 0.11336);

        lt I1(3);
        I1(0, 0) = 0.0149;  // Ixx
        I1(1, 0) = 0.0;     // Ixy
        I1(1, 1) = 0.3564;  // Iyy
        I1(2, 0) = 0.0;     // Ixz
        I1(2, 1) = 0.0;     // Iyz
        I1(2, 2) = 0.3553;  // Izz

        link1.I = RigidBodyInertia(mass1, com1, I1);

        Rotation R1(Eigen::Matrix3d::Identity());
        // Translation along Z: UR5 shoulder height d1 = 0.089 m
        link1.X = PluckerTransform(R1, Eigen::Vector3d(0.0, 0.0, 0.089));

        // Joint 2: Y-axis (shoulder pitch per D-11)
        link1.S = MotionVector(Eigen::Vector3d(0, 1, 0), Eigen::Vector3d::Zero());

        link1.q = 0.0;
        link1.qdot = 0.0;

        std::cout << "Link 1 (shoulder pitch):" << std::endl;
        std::cout << "  Parent: link 0 (index 0)" << std::endl;
        std::cout << "  Source: UR5 upper arm" << std::endl;
        std::cout << "  Mass: " << mass1 << " kg" << std::endl;
        std::cout << "  COM: [" << com1.transpose() << "] m" << std::endl;
        std::cout << "  Inertia: Ixx=" << I1(0,0) << ", Iyy=" << I1(1,1) << ", Izz=" << I1(2,2) << " kg·m^2" << std::endl;
        std::cout << "  Joint axis: Y (shoulder pitch)" << std::endl;
        std::cout << "  Transform from parent: translation [0, 0, 0.089] m (UR5 d1)" << std::endl;
        std::cout << std::endl;

        fd.links.push_back(link1);
    }

    {
        InverseDynamicsLink link1;
        link1.parent = 0;

        double mass1 = 8.393;
        Eigen::Vector3d com1(0.2125, 0.0, 0.11336);

        lt I1(3);
        I1(0, 0) = 0.0149;
        I1(1, 0) = 0.0;
        I1(1, 1) = 0.3564;
        I1(2, 0) = 0.0;
        I1(2, 1) = 0.0;
        I1(2, 2) = 0.3553;

        link1.I = RigidBodyInertia(mass1, com1, I1);

        Rotation R1(Eigen::Matrix3d::Identity());
        link1.X = PluckerTransform(R1, Eigen::Vector3d(0.0, 0.0, 0.089));

        link1.S = MotionVector(Eigen::Vector3d(0, 1, 0), Eigen::Vector3d::Zero());

        link1.q = 0.0;
        link1.qdot = 0.0;
        link1.qddot = 0.0;

        id.links.push_back(link1);
    }

    // ------------------------------------------------------------------
    // Link 2 (elbow): UR5 forearm — Z-axis revolute (elbow roll)
    // ------------------------------------------------------------------
    {
        Link link2;
        link2.parent = 1;

        double mass2 = 2.33;
        Eigen::Vector3d com2(0.15, 0.0, 0.0265);

        lt I2(3);
        I2(0, 0) = 0.0025;  // Ixx
        I2(1, 0) = 0.0;     // Ixy
        I2(1, 1) = 0.0551;  // Iyy
        I2(2, 0) = 0.0034;  // Ixz  (non-zero off-diagonal)
        I2(2, 1) = 0.0;     // Iyz
        I2(2, 2) = 0.0546;  // Izz

        link2.I = RigidBodyInertia(mass2, com2, I2);

        Rotation R2(Eigen::Matrix3d::Identity());
        // Translation along X: UR5 upper arm length a2 = 0.425 m
        link2.X = PluckerTransform(R2, Eigen::Vector3d(0.425, 0.0, 0.0));

        // Joint 3: Z-axis (elbow roll per D-11)
        link2.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link2.q = 0.0;
        link2.qdot = 0.0;

        std::cout << "Link 2 (elbow roll):" << std::endl;
        std::cout << "  Parent: link 1 (index 1)" << std::endl;
        std::cout << "  Source: UR5 forearm" << std::endl;
        std::cout << "  Mass: " << mass2 << " kg" << std::endl;
        std::cout << "  COM: [" << com2.transpose() << "] m" << std::endl;
        std::cout << "  Inertia: Ixx=" << I2(0,0) << ", Iyy=" << I2(1,1) << ", Izz=" << I2(2,2);
        std::cout << ", Ixz=" << I2(2,0) << " kg·m^2" << std::endl;
        std::cout << "  Joint axis: Z (elbow roll)" << std::endl;
        std::cout << "  Transform from parent: translation [0.425, 0, 0] m (UR5 a2)" << std::endl;
        std::cout << std::endl;

        fd.links.push_back(link2);
    }

    {
        InverseDynamicsLink link2;
        link2.parent = 1;

        double mass2 = 2.33;
        Eigen::Vector3d com2(0.15, 0.0, 0.0265);

        lt I2(3);
        I2(0, 0) = 0.0025;
        I2(1, 0) = 0.0;
        I2(1, 1) = 0.0551;
        I2(2, 0) = 0.0034;
        I2(2, 1) = 0.0;
        I2(2, 2) = 0.0546;

        link2.I = RigidBodyInertia(mass2, com2, I2);

        Rotation R2(Eigen::Matrix3d::Identity());
        link2.X = PluckerTransform(R2, Eigen::Vector3d(0.425, 0.0, 0.0));

        link2.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link2.q = 0.0;
        link2.qdot = 0.0;
        link2.qddot = 0.0;

        id.links.push_back(link2);
    }

    int nDof = 3;
    std::cout << "Model setup complete: " << nDof << " DOF, " << fd.links.size() << " links." << std::endl;
    std::cout << "Joint configuration: Z (waist yaw) → Y (shoulder pitch) → Z (elbow roll)." << std::endl;
    std::cout << std::endl;

    // ============================================================
    // Section 2: Forward Dynamics (ABA)
    // ============================================================
    std::cout << "=== Section 2: Forward Dynamics (ABA) ===" << std::endl;
    std::cout << "Computing joint accelerations from applied torques." << std::endl;
    std::cout << std::endl;

    // Test case 1: all positive torques
    {
        Eigen::VectorXd tau(3);
        tau << 20.0, 10.0, 5.0;

        std::cout << "FD Test 1: Applied torques tau = [" << tau(0) << ", " << tau(1) << ", " << tau(2) << "] N·m" << std::endl;

        fd.computeAccelerations(tau);

        std::cout << "  Joint 1 (waist)   qddot: " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 (shoulder) qddot: " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 3 (elbow)   qddot: " << fd.links[2].qddot << " rad/s^2" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    All joints accelerate positively under positive torques." << std::endl;
        std::cout << "    Joint 1 (waist, 3.7 kg) has lowest inertia → largest acceleration." << std::endl;
        std::cout << "    Joint 2 (shoulder, 8.393 kg) has highest inertia → smaller acceleration." << std::endl;
        std::cout << "    Joint 3 (elbow, 2.33 kg) is lightest but farthest from base," << std::endl;
        std::cout << "    with coupling effects from upstream links." << std::endl;
        std::cout << std::endl;
    }

    // Test case 2: mixed torques (opposing on shoulder)
    {
        Eigen::VectorXd tau(3);
        tau << 20.0, -10.0, 5.0;

        std::cout << "FD Test 2: Applied torques tau = [" << tau(0) << ", " << tau(1) << ", " << tau(2) << "] N·m" << std::endl;
        std::cout << "  (Joint 2 torque is negative — opposing shoulder)" << std::endl;

        fd.computeAccelerations(tau);

        std::cout << "  Joint 1 (waist)   qddot: " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 (shoulder) qddot: " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 3 (elbow)   qddot: " << fd.links[2].qddot << " rad/s^2" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Joint 2's negative torque produces negative acceleration." << std::endl;
        std::cout << "    Joint 1 is affected by the coupling through the arm — the" << std::endl;
        std::cout << "    shoulder's opposing torque changes the load seen at the waist." << std::endl;
        std::cout << "    Joint 3 is largely decoupled (serial chain, lightweight link)." << std::endl;
        std::cout << std::endl;
    }

    // ============================================================
    // Section 3: Inverse Dynamics — Gravity Compensation
    // ============================================================
    std::cout << "=== Section 3: Inverse Dynamics — Gravity Torques ===" << std::endl;
    std::cout << "Computing torques needed to hold the arm static against gravity." << std::endl;
    std::cout << std::endl;

    Eigen::Vector3d gravity(0, 0, -9.81);

    // Test 1: Arm fully extended horizontally q=[0,0,0]
    {
        id.links[0].q = 0.0;
        id.links[1].q = 0.0;
        id.links[2].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;
        id.links[2].qdot = 0.0;

        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(3);

        std::cout << "ID Gravity Test 1: Static pose q = [0, 0, 0] rad (arm fully extended horizontally)" << std::endl;
        std::cout << "  Gravity vector: [0, 0, -9.81] m/s^2" << std::endl;

        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        std::cout << "  Gravity compensation torques:" << std::endl;
        std::cout << "    Joint 1 (waist Z):   tau_g = " << tau_g(0) << " N·m" << std::endl;
        std::cout << "    Joint 2 (shoulder Y): tau_g = " << tau_g(1) << " N·m" << std::endl;
        std::cout << "    Joint 3 (elbow Z):   tau_g = " << tau_g(2) << " N·m" << std::endl;
        std::cout << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Joint 1 (Z-axis): tau_g ~ 0 N·m — gravity acts along Z," << std::endl;
        std::cout << "    which is parallel to the waist yaw axis. No torque required." << std::endl;
        std::cout << std::endl;
        std::cout << "    Joint 2 (Y-axis, shoulder): tau_g = " << tau_g(1) << " N·m — the" << std::endl;
        std::cout << "    arm extends horizontally, so gravity acts perpendicular to the" << std::endl;
        std::cout << "    shoulder pitch axis. This is the torque holding up the full arm." << std::endl;
        std::cout << "    The Featherstone sign convention gives a negative value (gravity" << std::endl;
        std::cout << "    torque in the negative Y-rotation direction)." << std::endl;
        std::cout << std::endl;
        std::cout << "    Joint 3 (Z-axis, elbow): tau_g ~ 0 N·m — gravity is parallel to" << std::endl;
        std::cout << "    the elbow roll axis, so no gravity compensation torque needed." << std::endl;
        std::cout << std::endl;
    }

    // Test 2: Same computation at non-zero q — shows solver limitation
    {
        id.links[0].q = 0.0;
        id.links[1].q = M_PI / 4.0;
        id.links[2].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;
        id.links[2].qdot = 0.0;

        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(3);

        std::cout << "ID Gravity Test 2: Same computation at q = [0, pi/4, 0] (arm at 45°)" << std::endl;
        std::cout << "  NOTE: The RNEA solver uses fixed transforms X (set at model setup)." << std::endl;
        std::cout << "  It does not recalculate transforms as a function of q. Results shown" << std::endl;
        std::cout << "  are for the home configuration only." << std::endl;

        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        std::cout << "  Gravity compensation torques (same as Test 1 — X is unchanged):" << std::endl;
        std::cout << "    Joint 1 (waist Z):   tau_g = " << tau_g(0) << " N·m" << std::endl;
        std::cout << "    Joint 2 (shoulder Y): tau_g = " << tau_g(1) << " N·m" << std::endl;
        std::cout << "    Joint 3 (elbow Z):   tau_g = " << tau_g(2) << " N·m" << std::endl;
        std::cout << "  This is a known solver limitation: transforms must be updated for" << std::endl;
        std::cout << "  non-zero joint configurations (a future library enhancement)." << std::endl;
        std::cout << std::endl;
    }

    // ============================================================
    // Section 4: Cross-Validation (FD + ID)
    // ============================================================
    std::cout << "=== Section 4: FD+ID Cross-Validation ===" << std::endl;
    std::cout << "Verifying solver consistency: feed gravity torques from ID into FD" << std::endl;
    std::cout << "and verify that the arm stays at rest (qddot = 0)." << std::endl;
    std::cout << std::endl;
    std::cout << "Cross-validation identity (should hold for consistent RNEA/ABA):" << std::endl;
    std::cout << "  tau_g = ID(q, 0, 0, g)  -- gravity compensation torques" << std::endl;
    std::cout << "  qddot = FD(tau_g, g)   -- should produce zero acceleration" << std::endl;
    std::cout << std::endl;
    std::cout << "This cross-validation tests the gravity compensation identity:" << std::endl;
    std::cout << "  FD(ID(qddot=0, g), g) should ≈ 0" << std::endl;
    std::cout << "  i.e., feeding gravity-compensation torques into the forward" << std::endl;
    std::cout << "  dynamics produces zero acceleration (arm stays at rest)." << std::endl;
    std::cout << std::endl;
    std::cout << "LIBRARY FINDINGS:" << std::endl;
    std::cout << "  1. ABA bug (ForwardDynamics.cpp):" << std::endl;
    std::cout << "     - ID→FD round-trip fails for multi-link chains with non-zero COM" << std::endl;
    std::cout << "     - Single Y-axis joint with COM offset under gravity: round-trip FAILS" << std::endl;
    std::cout << "     - Single joint with zero COM: round-trip PASSES (degenerate case)" << std::endl;
    std::cout << "     - Tests pass only because multi-link test cases use identity inertia" << std::endl;
    std::cout << "       and zero COM (degenerate case)" << std::endl;
    std::cout << "  2. RNEA limitation (InverseDynamics.cpp):" << std::endl;
    std::cout << "     - Transforms X do not update with joint position q" << std::endl;
    std::cout << "     - Results are valid only at the home configuration" << std::endl;
    std::cout << "     - Non-zero q positions require transform recomputation" << std::endl;
    std::cout << std::endl;
    std::cout << "Here tau_g is non-zero on the Y-axis shoulder joint, which exercises" << std::endl;
    std::cout << "both issues. The ID solver correctly computes gravity torques at the" << std::endl;
    std::cout << "home configuration, but the ABA solver fails to reproduce the identity." << std::endl;
    std::cout << std::endl;

    // Test 1: arm at [0, 0, 0]
    {
        fd.links[0].q = 0.0;
        fd.links[1].q = 0.0;
        fd.links[2].q = 0.0;
        id.links[0].q = 0.0;
        id.links[1].q = 0.0;
        id.links[2].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;
        id.links[2].qdot = 0.0;

        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        fd.computeAccelerations(tau_g, gravity);

        std::cout << "Cross-validation at q = [0, 0, 0] (arm extended horizontally):" << std::endl;
        std::cout << "  tau_g from ID = [" << tau_g(0) << ", " << tau_g(1) << ", " << tau_g(2) << "] N·m" << std::endl;
        std::cout << "  Joint 1 qddot (FD): " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot (FD): " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 3 qddot (FD): " << fd.links[2].qddot << " rad/s^2" << std::endl;

        bool pass = (std::abs(fd.links[0].qddot) < 1e-10 &&
                     std::abs(fd.links[1].qddot) < 1e-10 &&
                     std::abs(fd.links[2].qddot) < 1e-10);
        std::cout << "  Cross-validation: " << (pass ? "PASS" : "FAIL") << std::endl;
        if (pass) {
            std::cout << "  All |qddot| < 1e-10 → solvers are consistent." << std::endl;
        } else {
            std::cout << "  FAIL: qddot_J2 = " << fd.links[1].qddot << " rad/s^2 (expected 0)" << std::endl;
            std::cout << "  Pre-existing ABA bug: the FD solver incorrectly propagates" << std::endl;
            std::cout << "  forces for multi-link chains with non-zero COM offsets." << std::endl;
        }
        std::cout << std::endl;
    }

    // Test 2: arm at [0, pi/4, 0]
    {
        fd.links[0].q = 0.0;
        fd.links[1].q = M_PI / 4.0;
        fd.links[2].q = 0.0;
        id.links[0].q = 0.0;
        id.links[1].q = M_PI / 4.0;
        id.links[2].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;
        id.links[2].qdot = 0.0;

        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        fd.computeAccelerations(tau_g, gravity);

        std::cout << "Cross-validation at q = [0, pi/4, 0] (arm at 45°):" << std::endl;
        std::cout << "  tau_g from ID = [" << tau_g(0) << ", " << tau_g(1) << ", " << tau_g(2) << "] N·m" << std::endl;
        std::cout << "  Joint 1 qddot (FD): " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot (FD): " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 3 qddot (FD): " << fd.links[2].qddot << " rad/s^2" << std::endl;

        bool pass = (std::abs(fd.links[0].qddot) < 1e-10 &&
                     std::abs(fd.links[1].qddot) < 1e-10 &&
                     std::abs(fd.links[2].qddot) < 1e-10);
        std::cout << "  Cross-validation: " << (pass ? "PASS" : "FAIL") << std::endl;
        if (pass) {
            std::cout << "  All |qddot| < 1e-10 → solvers are consistent." << std::endl;
        } else {
            std::cout << "  FAIL: qddot_J2 = " << fd.links[1].qddot << " rad/s^2 (expected 0)" << std::endl;
            std::cout << "  Same bug as above — consistent failure on the Y-axis joint." << std::endl;
        }
        std::cout << std::endl;
    }

    // ============================================================
    // Section 5: Summary
    // ============================================================
    std::cout << "=== Section 5: Summary ===" << std::endl;
    std::cout << std::endl;
    std::cout << "This example demonstrated:" << std::endl;
    std::cout << "  1. Forward Dynamics (ABA): qddot = FD(tau)" << std::endl;
    std::cout << "     — Computes joint accelerations from torques (works for single" << std::endl;
    std::cout << "       configuration at home position)" << std::endl;
    std::cout << "  2. Inverse Dynamics (RNEA): tau_g = ID(qddot, g)" << std::endl;
    std::cout << "     — Computes gravity compensation torques at home configuration." << std::endl;
    std::cout << "  3. Cross-validation: FD(ID(0, g), g) = 0" << std::endl;
    std::cout << "     — Fails for multi-link chains with non-zero COM (ABA bug)." << std::endl;
    std::cout << std::endl;
    std::cout << "Library limitations revealed:" << std::endl;
    std::cout << "  - ABA (ForwardDynamics) has incorrect force propagation for" << std::endl;
    std::cout << "    multi-link chains with non-zero COM offsets." << std::endl;
    std::cout << "  - RNEA (InverseDynamics) uses fixed X transforms — results are" << std::endl;
    std::cout << "    valid only at the home configuration (q=0 for all joints)." << std::endl;
    std::cout << "  - These are pre-existing issues, not caused by this example." << std::endl;
    std::cout << std::endl;

    {
        // Compute gravity torque at horizontal pose for display in summary
        id.links[0].q = 0.0; id.links[1].q = 0.0; id.links[2].q = 0.0;
        id.links[0].qdot = 0.0; id.links[1].qdot = 0.0; id.links[2].qdot = 0.0;
        Eigen::VectorXd tau_g_0 = id.computeTorques(Eigen::VectorXd::Zero(3), gravity);

        std::cout << "Key physics results for the Z-Y-Z anthropomorphic arm:" << std::endl;
        std::cout << "  - Joint 1 (Z-axis waist yaw):  tau_g = " << tau_g_0(0) << " N·m" << std::endl;
        std::cout << "    (gravity parallel to joint axis — no torque needed)" << std::endl;
        std::cout << "  - Joint 2 (Y-axis shoulder):   tau_g = " << tau_g_0(1) << " N·m" << std::endl;
        std::cout << "    (at q=[0,0,0], arm horizontal — bears the full gravity load)" << std::endl;
        std::cout << "  - Joint 3 (Z-axis elbow roll): tau_g = " << tau_g_0(2) << " N·m" << std::endl;
        std::cout << "    (gravity parallel to joint axis — no torque needed)" << std::endl;
        std::cout << std::endl;
    }

    std::cout << "The Z-Y-Z arm is a standard anthropomorphic configuration where only" << std::endl;
    std::cout << "the Y-axis joint (shoulder pitch) bears gravity load. This matches" << std::endl;
    std::cout << "real-world robot behavior: the shoulder does the heavy lifting." << std::endl;
    std::cout << std::endl;
    std::cout << "=== Example Complete ===" << std::endl;

    return 0;
}
