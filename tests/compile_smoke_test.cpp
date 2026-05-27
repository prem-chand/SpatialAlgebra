/**
 * @brief Compile smoke test: verifies SpatialAlgebra.h is self-contained
 * @details This file tests that SpatialAlgebra.h can be included standalone
 *          without any prior includes. If this file compiles, the umbrella
 *          header is self-contained and includes all its dependencies in
 *          the correct order.
 */
#include "SpatialAlgebra.h"

int main() {
    // Verify the umbrella header provides access to all types
    using namespace SpatialAlgebra;

    // SpatialVector
    MotionVector mv(Vector3d(1, 0, 0), Vector3d(0, 1, 0));
    ForceVector fv(Vector3d(0, 0, 1), Vector3d(1, 0, 0));

    // Rotation
    Rotation R;
    R.setIdentity();

    // PluckerTransform
    PluckerTransform X(R, Vector3d::Zero());

    // LowerTriangular
    lt L = lt::Identity(3);

    // RigidBodyInertia
    RigidBodyInertia rbi(1.0, Vector3d::Zero(), L);

    // ArticulatedBodyInertia
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    ArticulatedBodyInertia abi(L, H, L);

    // SpatialOperations
    MotionVector result = cross(mv, mv);

    // ForwardDynamics
    ForwardDynamics fd;

    // InverseDynamics
    InverseDynamics id;

    // Verify value semantics work
    double d = mv.getAngular()[0];
    (void)d;

    return 0;
}
