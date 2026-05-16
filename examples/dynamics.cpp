/**
 * @file dynamics.cpp
 * @brief Demonstrates forward dynamics using the Articulated Body Algorithm (ABA).
 * 
 * This example shows:
 * - Setting up a kinematic chain using Link structures
 * - Configuring link properties (inertia, transform, joint axis)
 * - Computing joint accelerations from applied torques
 * - Understanding the ABA pipeline for articulated body simulation
 * 
 * The Articulated Body Algorithm (Featherstone Algorithm 7.3) computes
 * joint accelerations qdd from joint torques tau in O(n) time for an
 * n-degree-of-freedom articulated body.
 * 
 * Reference: Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 * Springer. Algorithm 7.3.
 */

#include "ForwardDynamics.h"
#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include "LowerTriangular.h"
#include <iostream>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
using lt = LowerTriangular;

int main() {
    std::cout << "=== SpatialAlgebra Forward Dynamics (ABA) Example ===" << std::endl;
    std::cout << std::endl;
    
    // Create a 2-link planar arm
    // Each link is a simple rigid body with uniform properties
    
    ForwardDynamics fd;
    std::cout << "Created forward dynamics solver" << std::endl;
    std::cout << std::endl;
    
    // Link 1: Base link (connected to world)
    // Mass: 1 kg, COM at origin, uniform inertia
    Link link1;
    link1.parent = -1;  // Base link (no parent)
    
    // Set inertia: mass=1kg, COM at origin, diagonal inertia
    lt I1_tensor(3);
    I1_tensor(0, 0) = 0.1; I1_tensor(1, 0) = 0.0; I1_tensor(1, 1) = 0.1;
    I1_tensor(2, 0) = 0.0; I1_tensor(2, 1) = 0.0; I1_tensor(2, 2) = 0.1;
    link1.I = RigidBodyInertia(1.0, Vector3d::Zero(), I1_tensor);
    
    // Transform from world to link 1: identity
    Rotation R1;
    R1.setIdentity();
    link1.X = PluckerTransform(R1, Vector3d::Zero());
    
    // Joint axis: rotation about Z axis (revolute joint)
    link1.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    
    // Initial joint state
    link1.q = 0.0;
    link1.qdot = 0.0;
    
    std::cout << "Link 1 (base):" << std::endl;
    std::cout << "  Mass: 1.0 kg" << std::endl;
    std::cout << "  Joint axis: [0, 0, 1] (rotation about Z)" << std::endl;
    std::cout << "  Parent: base (world)" << std::endl;
    std::cout << std::endl;
    
    fd.links.push_back(link1);
    
    // Link 2: Second link (connected to link 1)
    Link link2;
    link2.parent = 0;  // Parent is link 1 (index 0)
    
    // Same properties as link 1
    lt I2_tensor(3);
    I2_tensor(0, 0) = 0.1; I2_tensor(1, 0) = 0.0; I2_tensor(1, 1) = 0.1;
    I2_tensor(2, 0) = 0.0; I2_tensor(2, 1) = 0.0; I2_tensor(2, 2) = 0.1;
    link2.I = RigidBodyInertia(1.0, Vector3d::Zero(), I2_tensor);
    
    // Transform from link 1 to link 2: translate 0.5m along Z
    Rotation R2;
    R2.setIdentity();
    link2.X = PluckerTransform(R2, Vector3d(0, 0, 0.5));
    
    // Joint axis: rotation about Z axis
    link2.S = MotionVector(Vector3d(0, 0, 1), Vector3d::Zero());
    
    // Initial joint state
    link2.q = 0.0;
    link2.qdot = 0.0;
    
    std::cout << "Link 2:" << std::endl;
    std::cout << "  Mass: 1.0 kg" << std::endl;
    std::cout << "  Joint axis: [0, 0, 1] (rotation about Z)" << std::endl;
    std::cout << "  Parent: link 1" << std::endl;
    std::cout << "  Transform: [0, 0, 0.5] from parent" << std::endl;
    std::cout << std::endl;
    
    fd.links.push_back(link2);
    
    // Apply joint torques
    // Torque on joint 1: 1.0 N·m
    // Torque on joint 2: 0.5 N·m
    Eigen::VectorXd tau(2);
    tau << 1.0, 0.5;
    
    std::cout << "Applied joint torques:" << std::endl;
    std::cout << "  Joint 1: " << tau(0) << " N·m" << std::endl;
    std::cout << "  Joint 2: " << tau(1) << " N·m" << std::endl;
    std::cout << std::endl;
    
    // Compute joint accelerations using Articulated Body Algorithm
    // Results are stored in link.qddot for each link
    fd.computeAccelerations(tau);
    
    std::cout << "Computed joint accelerations (ABA):" << std::endl;
    std::cout << "  Joint 1: " << fd.links[0].qddot << " rad/s²" << std::endl;
    std::cout << "  Joint 2: " << fd.links[1].qddot << " rad/s²" << std::endl;
    std::cout << std::endl;
    
    // Verify results make physical sense
    std::cout << "Physical interpretation:" << std::endl;
    std::cout << "  Positive torque on joint 1 produces positive acceleration" << std::endl;
    std::cout << "  Joint 2 accelerates due to both its own torque and coupling" << std::endl;
    std::cout << "  Coupling between joints affects acceleration distribution" << std::endl;
    std::cout << std::endl;
    
    // Try different torque configuration
    std::cout << "=== Second Test: Opposing Torques ===" << std::endl;
    Eigen::VectorXd tau2(2);
    tau2 << 1.0, -1.0;
    
    std::cout << "Applied torques: [1.0, -1.0] N·m" << std::endl;
    fd.computeAccelerations(tau2);
    
    std::cout << "Resulting accelerations:" << std::endl;
    std::cout << "  Joint 1: " << fd.links[0].qddot << " rad/s²" << std::endl;
    std::cout << "  Joint 2: " << fd.links[1].qddot << " rad/s²" << std::endl;
    std::cout << std::endl;
    
    // Try zero torque (should give zero acceleration if at rest)
    std::cout << "=== Third Test: Zero Torques ===" << std::endl;
    Eigen::VectorXd tau_zero = Eigen::VectorXd::Zero(2);
    
    std::cout << "Applied torques: [0, 0] N·m" << std::endl;
    fd.computeAccelerations(tau_zero);
    
    std::cout << "Resulting accelerations:" << std::endl;
    std::cout << "  Joint 1: " << fd.links[0].qddot << " rad/s²" << std::endl;
    std::cout << "  Joint 2: " << fd.links[1].qddot << " rad/s²" << std::endl;
    std::cout << "  (Expected: zero or near-zero acceleration)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== Example Complete ===" << std::endl;
    std::cout << "The Articulated Body Algorithm computes forward dynamics in O(n) time." << std::endl;
    
    return 0;
}
