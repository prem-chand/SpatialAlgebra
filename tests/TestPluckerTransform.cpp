// write tests for plucker transform

#include "PluckerTransform.h"
#include "SpatialVector.h"

#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace SpatialAlgebra;
using namespace Eigen;

TEST(PluckerTransform, TransformMotion)
{

    // Create a rotation matrix
    Rotation E{Eigen::Matrix3d::Identity()};

    // Create a Plucker transform
    PluckerTransform transform(E, Vector3d::Zero());

    // Create a spatial vector
    SpatialVector motion({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});

    // Transform the spatial vector
    SpatialVector transformed = transform.transformMotion(motion);

    // Check the transformed spatial vector
    EXPECT_DOUBLE_EQ(transformed.getAngular()[0], 1.0);
    EXPECT_DOUBLE_EQ(transformed.getAngular()[1], 2.0);
    EXPECT_DOUBLE_EQ(transformed.getAngular()[2], 3.0);
    EXPECT_DOUBLE_EQ(transformed.getLinear()[0], 4.0);
    EXPECT_DOUBLE_EQ(transformed.getLinear()[1], 5.0);
    EXPECT_DOUBLE_EQ(transformed.getLinear()[2], 6.0);
}

int main()
{
    // create an mv object
    mv m1({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0});

    // create an fv object
    fv f1({2.0, 4.0, 6.0}, {4.0, 5.0, 6.0});

    Rotation E1{Eigen::Matrix3d::Identity()};

    plux p1(E1, Vector3d::Zero());
    plux p2(Rotation(2 * E1), Vector3d::Zero());
    auto m2 = p1.apply(m1);
    auto f2 = p1.apply(f1);

    auto p3 = p1.apply(p2);

    std::cout << "m2: ";
    m2.print();
    std::cout << "f2: ";
    f2.print();
    std::cout << "p3: ";
    p3.print();

    int k = 5;
    auto a = p1.get(k);
    auto b = p1.get(k);
}