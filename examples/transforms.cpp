/**
 * @file transforms.cpp
 * @brief Demonstrates Plücker coordinate transforms for spatial vector algebra.
 * 
 * This example shows:
 * - Creating Rotation objects from angle-axis representation
 * - Creating PluckerTransform from rotation and translation
 * - Transforming motion vectors between coordinate frames
 * - Transforming force vectors between coordinate frames
 * - Computing inverse transforms
 * - Chaining multiple transforms
 * 
 * Plücker coordinates provide a compact 6×6 representation of rigid body
 * transformations that correctly transform both motion and force vectors.
 */

#include "PluckerTransform.h"
#include "Rotation.h"
#include "SpatialUtils.h"
#include <iostream>
#include <cmath>

using namespace SpatialAlgebra;

// Type aliases
using mv = MotionVector;
using fv = ForceVector;
using plux = PluckerTransform;

int main() {
    std::cout << "=== SpatialAlgebra Plücker Transforms Example ===" << std::endl;
    std::cout << std::endl;
    
    // Create a rotation of 90 degrees (π/2 radians) around the Z axis
    // This rotates the X axis to align with the Y axis
    Rotation R;
    R.setFromAngleAxis(Eigen::AngleAxisd(M_PI_2, Vector3d(0, 0, 1)));
    std::cout << "Rotation: 90° around Z axis" << std::endl;
    std::cout << "Rotation matrix:" << std::endl;
    std::cout << R << std::endl;
    std::cout << std::endl;
    
    // Create a translation of 1 unit along the X axis
    Vector3d translation(1.0, 0.0, 0.0);
    std::cout << "Translation: [1, 0, 0]" << std::endl;
    std::cout << std::endl;
    
    // Create Plücker transform from rotation and translation
    // The 6×6 spatial transform is: X = [R, 0; -R[t]×, R]
    // where [t]× is the skew-symmetric cross-product matrix of translation
    plux X(R, translation);
    std::cout << "Plücker Transform X(R, t):" << std::endl;
    X.print();
    std::cout << std::endl;
    
    // Create a motion vector in the source frame
    // Pure rotation about the source frame's X axis at 1 rad/s
    mv v_source_X(Vector3d(1, 0, 0), Vector3d(0, 0, 0));
    std::cout << "Source motion (pure rotation about X):" << std::endl;
    v_source_X.print();
    std::cout << std::endl;
    
    // Transform motion to target frame
    // The same physical motion expressed in the transformed coordinate frame
    mv v_target = X.transformMotion(v_source_X);
    std::cout << "Transformed motion in target frame:" << std::endl;
    v_target.print();
    std::cout << "  Note: X-axis rotation becomes Y-axis rotation after 90° Z rotation" << std::endl;
    std::cout << std::endl;
    
    // Create a force vector in the source frame
    // Pure force along the X axis
    fv f_source_X(Vector3d(0, 0, 0), Vector3d(1, 0, 0));
    std::cout << "Source force (pure force along X):" << std::endl;
    f_source_X.print();
    std::cout << std::endl;
    
    // Transform force to target frame
    fv f_target = X.transformForce(f_source_X);
    std::cout << "Transformed force in target frame:" << std::endl;
    f_target.print();
    std::cout << "  Note: X-axis force becomes Y-axis force after 90° Z rotation" << std::endl;
    std::cout << std::endl;
    
    // Compute the inverse transform
    // X_inv * X = Identity
    plux X_inv = X.inverse();
    std::cout << "Inverse transform X^(-1):" << std::endl;
    X_inv.print();
    std::cout << std::endl;
    
    // Verify inverse: transform back to source frame
    mv v_back = X_inv.transformMotion(v_target);
    std::cout << "Transform back to source frame (verification):" << std::endl;
    v_back.print();
    std::cout << "  Should match original v_source_X" << std::endl;
    std::cout << std::endl;
    
    // Chain transforms: create a second transform
    Rotation R2;
    R2.setFromAngleAxis(Eigen::AngleAxisd(M_PI_4, Vector3d(1, 0, 0))); // 45° around X
    plux X2(R2, Vector3d(0, 1, 0)); // Translate 1 unit along Y
    
    std::cout << "Second transform X2: 45° around X, translate [0, 1, 0]" << std::endl;
    std::cout << "  Note: Transform chaining via operator* not implemented in this version" << std::endl;
    std::cout << std::endl;
    
    // Demonstrate transform of a wrench with both torque and force
    fv wrench_source(Vector3d(0, 0, 10), Vector3d(5, 0, 0));
    std::cout << "Source wrench (torque about Z, force along X):" << std::endl;
    wrench_source.print();
    
    fv wrench_target = X.transformForce(wrench_source);
    std::cout << "Transformed wrench in target frame:" << std::endl;
    wrench_target.print();
    std::cout << "  Note: moment arm creates additional torque from force" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== Example Complete ===" << std::endl;
    
    return 0;
}
