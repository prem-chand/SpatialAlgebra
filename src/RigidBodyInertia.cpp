#include "RigidBodyInertia.h"
#include "PluckerTransform.h"
#include <Eigen/Dense>
#include <iostream>

using namespace SpatialAlgebra;
using namespace Eigen;



// Transform inertia
// RigidBodyInertia RigidBodyInertia::transform(const PluckerTransform &transform) const
// {
//     // Transform the center of mass
//     Vector3d transformedCom = transform.transformPoint(com);

//     // Transform the inertia matrix
//     Matrix6d transformMatrix = transform.getTransformationMatrix();
//     Vector6d transformedInertiaMatrixLT = transformMatrix * inertiaMatrixLT;

//     std::array<double, 3> transformedComArray = {transformedCom[0], transformedCom[1], transformedCom[2]};
//     std::array<double, 6> transformedInertiaMatrixLTArray;
//     for (size_t i = 0; i < 6; ++i)
//     {
//         transformedInertiaMatrixLTArray[i] = transformedInertiaMatrixLT[i];
//     }

//     return RigidBodyInertia(mass, transformedComArray, transformedInertiaMatrixLTArray);
// }

// Basic operations




