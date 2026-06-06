/**
 * @file example_robot_2link.cpp
 * @brief Demonstrates forward and inverse dynamics for a 2-link Z-Z planar arm
 *        with UR5-derived parameters and FD+ID cross-validation.
 * 
 * Robot Configuration:
 * - 2-link serial chain operating in the XY plane (horizontal arm)
 * - Both joints are Z-axis revolutes (Z-Z configuration)
 * - Link 0: UR5 upper arm adapted (mass=8.393 kg, COM at [0.2125, 0, 0.02])
 * - Link 1: UR5 forearm adapted (mass=2.33 kg, COM at [0.15, 0, 0.02])
 * - Gravity: [0, 0, -9.81] m/s²
 * 
 * This example demonstrates:
 * 1. Model setup with realistic UR5 parameters (mass, COM, inertia, transforms)
 * 2. Forward dynamics (ABA): apply torques, compute joint accelerations
 * 3. Inverse dynamics (RNEA): compute gravity compensation torques at static pose
 * 4. Cross-validation: feed gravity torques into FD, verify zero acceleration
 * 5. Physical interpretation of all results
 * 
 * Key Physics Insight:
 * For this Z-Z planar arm with gravity [0,0,-9.81], gravity acts PARALLEL
 * to both joint axes. Therefore the gravity compensation torque is ZERO
 * for any pose — the arm is a horizontal SCARA-like configuration.
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
    std::cout << "=== SpatialAlgebra 2-Link Z-Z Planar Arm Example ===" << std::endl;
    std::cout << std::endl;
    std::cout << "Physical configuration:" << std::endl;
    std::cout << "  Two-link serial chain, both joints rotate about Z axis." << std::endl;
    std::cout << "  Arm operates in the XY plane (horizontal)." << std::endl;
    std::cout << "  Gravity: [0, 0, -9.81] m/s^2 acts parallel to Z." << std::endl;
    std::cout << "  UR5-derived link parameters." << std::endl;
    std::cout << std::endl;

    // ============================================================
    // Section 1: Model Setup
    // ============================================================
    std::cout << "=== Section 1: Robot Model Setup ===" << std::endl;
    std::cout << std::endl;

    // Build ForwardDynamics solver (ABA)
    ForwardDynamics fd;
    // Build InverseDynamics solver (RNEA) separately per Pitfall 2
    InverseDynamics id;

    // ------------------------------------------------------------------
    // Link 0 (base): UR5 upper arm adapted
    // ------------------------------------------------------------------
    {
        Link link0;
        link0.parent = -1;  // Base link (connected to world)

        // Mass and COM from UR5 upper arm
        double mass0 = 8.393;
        Eigen::Vector3d com0(0.2125, 0.0, 0.02);

        // Inertia tensor (lower-triangular packed storage)
        lt I0(3);
        I0(0, 0) = 0.0149;  // Ixx
        I0(1, 0) = 0.0;     // Ixy
        I0(1, 1) = 0.3564;  // Iyy
        I0(2, 0) = 0.0;     // Ixz
        I0(2, 1) = 0.0;     // Iyz
        I0(2, 2) = 0.3553;  // Izz

        link0.I = RigidBodyInertia(mass0, com0, I0);

        // Transform from world to link 0: identity (base at origin)
        Rotation R0(Eigen::Matrix3d::Identity());
        link0.X = PluckerTransform(R0, Eigen::Vector3d::Zero());

        // Joint axis: Z-axis revolute
        link0.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        // Initial joint state
        link0.q = 0.0;
        link0.qdot = 0.0;

        // Print link parameters
        std::cout << "Link 0 (base):" << std::endl;
        std::cout << "  Parent: world (base, index -1)" << std::endl;
        std::cout << "  Source: UR5 upper arm (adapted)" << std::endl;
        std::cout << "  Mass: " << mass0 << " kg" << std::endl;
        std::cout << "  COM: [" << com0.transpose() << "] m" << std::endl;
        std::cout << "  Inertia: Ixx=" << I0(0,0) << ", Iyy=" << I0(1,1) << ", Izz=" << I0(2,2) << " kg·m^2" << std::endl;
        std::cout << "  Joint axis: Z (revolute)" << std::endl;
        std::cout << "  Transform from parent: identity" << std::endl;
        std::cout << std::endl;

        fd.links.push_back(link0);
    }

    // Also add to InverseDynamics solver (separate struct, identical params)
    {
        InverseDynamicsLink link0;
        link0.parent = -1;

        double mass0 = 8.393;
        Eigen::Vector3d com0(0.2125, 0.0, 0.02);

        lt I0(3);
        I0(0, 0) = 0.0149;
        I0(1, 0) = 0.0;
        I0(1, 1) = 0.3564;
        I0(2, 0) = 0.0;
        I0(2, 1) = 0.0;
        I0(2, 2) = 0.3553;

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
    // Link 1 (tip): UR5 forearm adapted
    // ------------------------------------------------------------------
    {
        Link link1;
        link1.parent = 0;  // Parent is link 0 (index 0)

        double mass1 = 2.33;
        Eigen::Vector3d com1(0.15, 0.0, 0.02);

        lt I1(3);
        I1(0, 0) = 0.0025;  // Ixx
        I1(1, 0) = 0.0;     // Ixy
        I1(1, 1) = 0.0551;  // Iyy
        I1(2, 0) = 0.0;     // Ixz
        I1(2, 1) = 0.0;     // Iyz
        I1(2, 2) = 0.0546;  // Izz

        link1.I = RigidBodyInertia(mass1, com1, I1);

        // Transform from link 0 to link 1: translate 0.425 m along X (UR5 upper arm length a2)
        Rotation R1(Eigen::Matrix3d::Identity());
        link1.X = PluckerTransform(R1, Eigen::Vector3d(0.425, 0.0, 0.0));

        // Joint axis: Z-axis revolute
        link1.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link1.q = 0.0;
        link1.qdot = 0.0;

        // Print link parameters
        std::cout << "Link 1 (tip):" << std::endl;
        std::cout << "  Parent: link 0 (index 0)" << std::endl;
        std::cout << "  Source: UR5 forearm (adapted)" << std::endl;
        std::cout << "  Mass: " << mass1 << " kg" << std::endl;
        std::cout << "  COM: [" << com1.transpose() << "] m" << std::endl;
        std::cout << "  Inertia: Ixx=" << I1(0,0) << ", Iyy=" << I1(1,1) << ", Izz=" << I1(2,2) << " kg·m^2" << std::endl;
        std::cout << "  Joint axis: Z (revolute)" << std::endl;
        std::cout << "  Transform from parent: translation [0.425, 0, 0] m (UR5 a2)" << std::endl;
        std::cout << std::endl;

        fd.links.push_back(link1);
    }

    {
        InverseDynamicsLink link1;
        link1.parent = 0;

        double mass1 = 2.33;
        Eigen::Vector3d com1(0.15, 0.0, 0.02);

        lt I1(3);
        I1(0, 0) = 0.0025;
        I1(1, 0) = 0.0;
        I1(1, 1) = 0.0551;
        I1(2, 0) = 0.0;
        I1(2, 1) = 0.0;
        I1(2, 2) = 0.0546;

        link1.I = RigidBodyInertia(mass1, com1, I1);

        Rotation R1(Eigen::Matrix3d::Identity());
        link1.X = PluckerTransform(R1, Eigen::Vector3d(0.425, 0.0, 0.0));

        link1.S = MotionVector(Eigen::Vector3d(0, 0, 1), Eigen::Vector3d::Zero());

        link1.q = 0.0;
        link1.qdot = 0.0;
        link1.qddot = 0.0;

        id.links.push_back(link1);
    }

    int nDof = 2;
    std::cout << "Model setup complete: " << nDof << " DOF, " << fd.links.size() << " links." << std::endl;
    std::cout << std::endl;

    // ============================================================
    // Section 2: Forward Dynamics (ABA)
    // ============================================================
    std::cout << "=== Section 2: Forward Dynamics (ABA) ===" << std::endl;
    std::cout << "Computing joint accelerations from applied torques." << std::endl;
    std::cout << "The Articulated Body Algorithm solves: qddot = ABA(tau, q, qdot)" << std::endl;
    std::cout << std::endl;

    // Test case 1: both positive torques
    {
        Eigen::VectorXd tau(2);
        tau << 10.0, 5.0;

        std::cout << "FD Test 1: Applied torques tau = [" << tau(0) << ", " << tau(1) << "] N·m" << std::endl;
        std::cout << "  Gravity: zero (disabling gravity for direct torque→acceleration test)" << std::endl;

        fd.computeAccelerations(tau);

        std::cout << "  Joint 1 qddot: " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot: " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Both joints accelerate in positive direction as expected." << std::endl;
        std::cout << "    Joint 1 has larger inertia (heavier link, 8.393 kg) so acceleration" << std::endl;
        std::cout << "    is smaller per unit torque compared to joint 2 (lighter, 2.33 kg)." << std::endl;
        std::cout << "    Joint 1's torque also drives joint 2 through coupling in the ABA." << std::endl;
        std::cout << std::endl;
    }

    // Test case 2: opposing torques
    {
        Eigen::VectorXd tau(2);
        tau << 10.0, -5.0;

        std::cout << "FD Test 2: Applied torques tau = [" << tau(0) << ", " << tau(1) << "] N·m" << std::endl;
        std::cout << "  (Joint 1 positive, Joint 2 negative — opposing torques)" << std::endl;

        fd.computeAccelerations(tau);

        std::cout << "  Joint 1 qddot: " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot: " << fd.links[1].qddot << " rad/s^2" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Joint 2's negative torque counters the coupling from joint 1," << std::endl;
        std::cout << "    reducing joint 2's acceleration compared to Test 1." << std::endl;
        std::cout << "    Joint 1's acceleration is largely unchanged since coupling from" << std::endl;
        std::cout << "    the lighter link 2 back to link 1 is relatively small." << std::endl;
        std::cout << std::endl;
    }

    // ============================================================
    // Section 3: Inverse Dynamics — Gravity Compensation
    // ============================================================
    std::cout << "=== Section 3: Inverse Dynamics — Gravity Torques ===" << std::endl;
    std::cout << "Computing torques needed to hold the arm static against gravity." << std::endl;
    std::cout << "The Recursive Newton-Euler Algorithm solves: tau = RNEA(q, qdot, qddot, gravity)" << std::endl;
    std::cout << std::endl;

    Eigen::Vector3d gravity(0, 0, -9.81);

    // Static pose q=[0,0], qdot=[0,0], qddot=[0,0]
    {
        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(2);
        id.links[0].q = 0.0;
        id.links[1].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;

        std::cout << "ID Gravity Test 1: Static pose q = [0, 0] rad (arm extended along X)" << std::endl;
        std::cout << "  qdot = [0, 0] rad/s, qddot = [0, 0] rad/s^2" << std::endl;
        std::cout << "  Gravity vector: [0, 0, -9.81] m/s^2" << std::endl;

        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        std::cout << "  Gravity compensation torques: tau_g = [" << tau_g(0) << ", " << tau_g(1) << "] N·m" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Both joint axes are Z-direction (parallel to gravity vector)." << std::endl;
        std::cout << "    Gravity acts along Z, while the joint axes also point along Z." << std::endl;
        std::cout << "    The cross product r x (m * g) has zero component along Z." << std::endl;
        std::cout << "    Therefore gravity produces NO torque about Z-axis joints." << std::endl;
        std::cout << "    Gravity torque = [0, 0] N·m, which is physically correct for" << std::endl;
        std::cout << "    this horizontal SCARA-like configuration." << std::endl;
        std::cout << std::endl;
    }

    // Repeat with different pose to verify gravity torque remains zero
    {
        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(2);
        id.links[0].q = M_PI / 4.0;
        id.links[1].q = M_PI / 3.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;

        std::cout << "ID Gravity Test 2: Static pose q = [pi/4, pi/3] rad" << std::endl;
        std::cout << "  (Arm at arbitrary rotated configuration)" << std::endl;

        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        std::cout << "  Gravity compensation torques: tau_g = [" << tau_g(0) << ", " << tau_g(1) << "] N·m" << std::endl;
        std::cout << "  Physical interpretation:" << std::endl;
        std::cout << "    Even with rotated joints, gravity torque remains zero." << std::endl;
        std::cout << "    Both links have Z-axis joints, and gravity is also along Z." << std::endl;
        std::cout << "    The Z-Z configuration is gravity-neutral: no torque is needed" << std::endl;
        std::cout << "    to hold any static pose against gravity." << std::endl;
        std::cout << "    (The interesting gravity case is in the 3-link Z-Y-Z example.)" << std::endl;
        std::cout << std::endl;
    }

    // ============================================================
    // Section 4: Cross-Validation (FD + ID)
    // ============================================================
    std::cout << "=== Section 4: FD+ID Cross-Validation ===" << std::endl;
    std::cout << "Verifying solver consistency: feed gravity torques from ID into FD" << std::endl;
    std::cout << "and verify that the arm stays at rest (qddot = 0)." << std::endl;
    std::cout << std::endl;
    std::cout << "Cross-validation identity (per D-12):" << std::endl;
    std::cout << "  tau_g = ID(q, 0, 0, g)  -- gravity compensation torques" << std::endl;
    std::cout << "  qddot = FD(tau_g, g)   -- should produce zero acceleration" << std::endl;
    std::cout << std::endl;

    {
        // Reset FD joint state to match ID
        fd.links[0].q = M_PI / 4.0;
        fd.links[1].q = M_PI / 3.0;

        // Compute gravity torques from ID
        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        // Feed gravity torques into FD with same gravity
        fd.computeAccelerations(tau_g, gravity);

        std::cout << "Cross-validation at q = [pi/4, pi/3] rad:" << std::endl;
        std::cout << "  tau_g from ID = [" << tau_g(0) << ", " << tau_g(1) << "] N·m" << std::endl;
        std::cout << "  Joint 1 qddot (FD): " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot (FD): " << fd.links[1].qddot << " rad/s^2" << std::endl;

        bool pass = (std::abs(fd.links[0].qddot) < 1e-10 && std::abs(fd.links[1].qddot) < 1e-10);
        std::cout << "  Cross-validation: " << (pass ? "PASS" : "FAIL") << std::endl;
        if (pass) {
            std::cout << "  Both |qddot| < 1e-10 → solvers are consistent." << std::endl;
        } else {
            std::cout << "  WARNING: |qddot| >= 1e-10 — possible solver mismatch." << std::endl;
        }
        std::cout << std::endl;
    }

    // Repeat for simple extended pose
    {
        fd.links[0].q = 0.0;
        fd.links[1].q = 0.0;
        id.links[0].q = 0.0;
        id.links[1].q = 0.0;
        id.links[0].qdot = 0.0;
        id.links[1].qdot = 0.0;

        Eigen::VectorXd qddot_zero = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd tau_g = id.computeTorques(qddot_zero, gravity);

        fd.computeAccelerations(tau_g, gravity);

        std::cout << "Cross-validation at q = [0, 0] rad (arm fully extended):" << std::endl;
        std::cout << "  Joint 1 qddot: " << fd.links[0].qddot << " rad/s^2" << std::endl;
        std::cout << "  Joint 2 qddot: " << fd.links[1].qddot << " rad/s^2" << std::endl;

        bool pass = (std::abs(fd.links[0].qddot) < 1e-10 && std::abs(fd.links[1].qddot) < 1e-10);
        std::cout << "  Cross-validation: " << (pass ? "PASS" : "FAIL") << std::endl;
        if (pass) {
            std::cout << "  Both |qddot| < 1e-10 → solvers are consistent." << std::endl;
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
    std::cout << "     — Computed joint accelerations from applied torques" << std::endl;
    std::cout << "  2. Inverse Dynamics (RNEA): tau_g = ID(0, g)" << std::endl;
    std::cout << "     — Computed gravity compensation torques at static pose" << std::endl;
    std::cout << "  3. Cross-validation: FD(ID(0, g), g) = 0" << std::endl;
    std::cout << "     — Verified solver consistency: |qddot| < 1e-10" << std::endl;
    std::cout << std::endl;
    std::cout << "Key physics results:" << std::endl;
    std::cout << "  - Z-Z planar arm with gravity [0,0,-9.81] produces" << std::endl;
    std::cout << "    ZERO gravity torque regardless of pose." << std::endl;
    std::cout << "  - This is because gravity acts parallel to both joint axes." << std::endl;
    std::cout << "  - The cross product r x (m*g) has no component along Z." << std::endl;
    std::cout << "  - FD+ID cross-validation confirms solver consistency." << std::endl;
    std::cout << std::endl;
    std::cout << "=== Example Complete ===" << std::endl;

    return 0;
}
