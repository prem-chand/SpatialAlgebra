/**
 * @file inertia.cpp
 * @brief Demonstrates rigid body inertia operations.
 * 
 * This example shows:
 * - Creating RigidBodyInertia with mass, center of mass, and inertia tensor
 * - Applying inertia to motion vectors to compute forces
 * - Transforming inertias by Plücker transforms
 * - Physical interpretation of inertia operations
 * 
 * In spatial vector algebra, rigid body inertia maps motion vectors to force
 * vectors: F = I * v. The inertia includes both mass and rotational inertia.
 */

#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include "Rotation.h"
#include "LowerTriangular.h"
#include <iostream>

using namespace SpatialAlgebra;

// Type aliases
using mv = MotionVector;
using fv = ForceVector;
using plux = PluckerTransform;
using rbi = RigidBodyInertia;
using lt = LowerTriangular;

int main() {
    std::cout << "=== SpatialAlgebra Inertia Example ===" << std::endl;
    std::cout << std::endl;
    
    // Create inertia tensor in LowerTriangular packed format
    // For a uniform sphere: I = (2/5) * m * r^2
    // For m=2kg, r=0.5m: I = (2/5) * 2 * 0.25 = 0.2 kg·m²
    double mass = 2.0;  // kg
    double radius = 0.5; // m
    double I_sphere = 0.4 * mass * radius * radius;
    
    // Create LowerTriangular inertia matrix
    // Packed storage: [Ixx, Ixy, Iyy, Ixz, Iyz, Izz]
    lt inertia_tensor(3);
    inertia_tensor(0, 0) = I_sphere;  // Ixx
    inertia_tensor(1, 0) = 0.0;       // Ixy
    inertia_tensor(1, 1) = I_sphere;  // Iyy
    inertia_tensor(2, 0) = 0.0;       // Ixz
    inertia_tensor(2, 1) = 0.0;       // Iyz
    inertia_tensor(2, 2) = I_sphere;  // Izz
    
    // Center of mass at origin
    Vector3d com(0, 0, 0);
    
    // Create rigid body inertia
    rbi I_sphere_body(mass, com, inertia_tensor);
    
    std::cout << "Rigid Body Inertia: Uniform Sphere" << std::endl;
    std::cout << "  Mass: " << mass << " kg" << std::endl;
    std::cout << "  Radius: " << radius << " m" << std::endl;
    std::cout << "  Center of Mass: [0, 0, 0]" << std::endl;
    std::cout << "  Rotational Inertia: " << I_sphere << " kg·m² (all axes)" << std::endl;
    std::cout << "Full inertia:" << std::endl;
    I_sphere_body.print();
    std::cout << std::endl;
    
    // Create a motion vector: pure rotation about Z axis at 1 rad/s
    mv rotation_z(Vector3d(0, 0, 1), Vector3d(0, 0, 0));
    std::cout << "Applied motion: pure rotation about Z at 1 rad/s" << std::endl;
    rotation_z.print();
    std::cout << std::endl;
    
    // Apply inertia to compute the resulting force (torque)
    // F = I * v => torque = I * angular_velocity
    fv torque_result = I_sphere_body.apply(rotation_z);
    std::cout << "Resulting force (torque to produce this rotation):" << std::endl;
    torque_result.print();
    std::cout << "  Expected: [0, 0, " << I_sphere << "] N·m (torque about Z)" << std::endl;
    std::cout << std::endl;
    
    // Create a motion with linear acceleration
    mv linear_x(Vector3d(0, 0, 0), Vector3d(1, 0, 0));
    std::cout << "Applied motion: pure linear acceleration along X at 1 m/s²" << std::endl;
    linear_x.print();
    
    fv force_result = I_sphere_body.apply(linear_x);
    std::cout << "Resulting force:" << std::endl;
    force_result.print();
    std::cout << "  Expected: [" << mass << ", 0, 0] N (F = m*a = 2*1)" << std::endl;
    std::cout << std::endl;
    
    // Create rigid body inertia with offset center of mass
    // This creates coupling between linear and angular components
    Vector3d com_offset(0, 0, 0.5); // COM 0.5m along Z axis
    rbi I_offset(mass, com_offset, inertia_tensor);
    
    std::cout << "Rigid Body Inertia with offset COM:" << std::endl;
    std::cout << "  Mass: " << mass << " kg" << std::endl;
    std::cout << "  Center of Mass: [0, 0, 0.5] m" << std::endl;
    I_offset.print();
    std::cout << std::endl;
    
    // Apply pure linear acceleration to offset inertia
    // This produces both force AND torque due to the offset COM
    fv offset_result = I_offset.apply(linear_x);
    std::cout << "Force from linear acceleration with offset COM:" << std::endl;
    offset_result.print();
    std::cout << "  Note: Linear acceleration produces torque due to moment arm" << std::endl;
    std::cout << std::endl;
    
    // Transform inertia by a Plücker transform
    // This computes the inertia as seen from a different coordinate frame
    Rotation R;
    R.setIdentity();
    Vector3d translation(0, 0, 0.1); // Translate 0.1m along Z
    plux X(R, translation);
    
    rbi I_transformed = X.tformRBI(I_sphere_body);
    std::cout << "Inertia transformed by X(R, [0,0,0.1]):" << std::endl;
    I_transformed.print();
    std::cout << "  Transformed to frame 0.1m away along Z axis" << std::endl;
    std::cout << std::endl;
    
    // Get inertia properties
    std::cout << "Inertia properties:" << std::endl;
    std::cout << "  Mass: " << I_sphere_body.getMass() << " kg" << std::endl;
    std::cout << "  COM: [" << I_sphere_body.getCom().transpose() << "] m" << std::endl;
    
    // Get inertia matrix in lower triangular packed format
    lt I_lt = I_sphere_body.getInertiaMatrixLT();
    std::cout << "  Inertia matrix stored in packed LowerTriangular format" << std::endl;
    std::cout << std::endl;
    
    // Combined motion and force
    mv combined_motion(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    std::cout << "Combined motion (rotation about X, translation along Y):" << std::endl;
    combined_motion.print();
    
    fv combined_force = I_sphere_body.apply(combined_motion);
    std::cout << "Resulting force:" << std::endl;
    combined_force.print();
    std::cout << std::endl;
    
    std::cout << "=== Example Complete ===" << std::endl;
    
    return 0;
}
