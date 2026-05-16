/**
 * @file basic_vectors.cpp
 * @brief Demonstrates basic spatial vector operations with MotionVector and ForceVector.
 * 
 * This example shows:
 * - Creating motion vectors (twists) representing spatial velocity
 * - Creating force vectors (wrenches) representing spatial force
 * - Vector operations: addition, subtraction, scaling
 * - Cross products between spatial vectors
 * - Using type aliases for concise notation
 */

#include "MotionVector.h"
#include "ForceVector.h"
#include "SpatialUtils.h"
#include <iostream>

using namespace SpatialAlgebra;

// Type aliases for concise notation
using mv = MotionVector;
using fv = ForceVector;

int main() {
    std::cout << "=== SpatialAlgebra Basic Vectors Example ===" << std::endl;
    std::cout << std::endl;
    
    // Create a motion vector (twist) with angular velocity [1, 0, 0] rad/s
    // and linear velocity [0, 1, 0] m/s
    // Physical interpretation: rotation about X-axis, translation along Y-axis
    mv twist1(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    std::cout << "Motion Vector 1 (twist):" << std::endl;
    std::cout << "  Angular: [1, 0, 0] rad/s (rotation about X)" << std::endl;
    std::cout << "  Linear:  [0, 1, 0] m/s (translation along Y)" << std::endl;
    twist1.print();
    std::cout << std::endl;
    
    // Create another motion vector
    mv twist2(Vector3d(0, 1, 0), Vector3d(0, 0, 1));
    std::cout << "Motion Vector 2 (twist):" << std::endl;
    std::cout << "  Angular: [0, 1, 0] rad/s (rotation about Y)" << std::endl;
    std::cout << "  Linear:  [0, 0, 1] m/s (translation along Z)" << std::endl;
    twist2.print();
    std::cout << std::endl;
    
    // Vector addition: combine two twists
    mv twist_sum = twist1 + twist2;
    std::cout << "Sum of twist1 + twist2:" << std::endl;
    twist_sum.print();
    std::cout << std::endl;
    
    // Vector subtraction
    mv twist_diff = twist1 - twist2;
    std::cout << "Difference twist1 - twist2:" << std::endl;
    twist_diff.print();
    std::cout << std::endl;
    
    // Scalar multiplication
    mv twist_scaled = twist1 * 2.0;
    std::cout << "twist1 scaled by 2.0:" << std::endl;
    twist_scaled.print();
    std::cout << std::endl;
    
    // Create a force vector (wrench) with torque [0, 0, 5] N·m and force [10, 0, 0] N
    // Physical interpretation: torque about Z-axis, force along X-axis
    fv wrench1(Vector3d(0, 0, 5), Vector3d(10, 0, 0));
    std::cout << "Force Vector 1 (wrench):" << std::endl;
    std::cout << "  Torque: [0, 0, 5] N·m (torque about Z)" << std::endl;
    std::cout << "  Force:  [10, 0, 0] N (force along X)" << std::endl;
    wrench1.print();
    std::cout << std::endl;
    
    // Create another force vector
    fv wrench2(Vector3d(1, 0, 0), Vector3d(0, 5, 0));
    std::cout << "Force Vector 2 (wrench):" << std::endl;
    wrench2.print();
    std::cout << std::endl;
    
    // Force vector addition
    fv wrench_sum = wrench1 + wrench2;
    std::cout << "Sum of wrench1 + wrench2:" << std::endl;
    wrench_sum.print();
    std::cout << std::endl;
    
    // Cross product: motion cross force produces a force
    // Physical interpretation: the effect of moving a force through a velocity field
    fv cross_result = twist1.crossForce(wrench1);
    std::cout << "Cross product twist1 × wrench1 (produces force):" << std::endl;
    cross_result.print();
    std::cout << std::endl;
    
    // Cross product: motion cross motion produces a motion
    mv motion_cross = twist1.crossMotion(twist2);
    std::cout << "Cross product twist1 × twist2 (produces motion):" << std::endl;
    motion_cross.print();
    std::cout << std::endl;
    
    // Dot product between motion and force (power)
    // Returns scalar: power = ω·τ + v·f
    double power = dot(twist1, wrench1);
    std::cout << "Dot product twist1 · wrench1 (power):" << std::endl;
    std::cout << "  Power = " << power << " W" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== Example Complete ===" << std::endl;
    
    return 0;
}
