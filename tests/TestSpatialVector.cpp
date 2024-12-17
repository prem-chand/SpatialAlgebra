#include "SpatialVector.h"
#include <cassert>

using namespace SpatialAlgebra;

int main()
{
    SpatialVector v1({1, 2, 3}, {4, 5, 6});
    SpatialVector v2 = v1 * 2.0;

    assert(v2.getAngular()[0] == 2);
    assert(v2.getLinear()[1] == 10);

    v1.print();

    return 0;
}
