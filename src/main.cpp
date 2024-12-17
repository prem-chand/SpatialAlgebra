#include "SpatialVector.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include <iostream>
#include <Eigen/Geometry> // Include for RotationMatrix

using namespace SpatialAlgebra;
using lt = LowerTriangular;

int main()
{
    // Spatial vector
    SpatialVector motion({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});
    motion.print();

    Rotation E{Eigen::Matrix3d::Identity()};

    // Plücker transform
    PluckerTransform transform(E, Vector3d::Zero());
    SpatialVector transformed = transform.transformMotion(motion);
    transformed.print();

    // Rigid-body inertia
    RigidBodyInertia inertia(10.0, Vector3d::Zero(), lt(Vector6d::Zero(),3));
    inertia.print();

    return 0;
}
