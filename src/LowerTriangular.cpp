#include "LowerTriangular.h"

using namespace SpatialAlgebra;

Eigen::VectorXd LowerTriangular::operator*(const Eigen::VectorXd& v) const
{
    // Lower-triangular only multiplication (symmetric access via operator())
    int n = 3;
    if (v.size() != n)
        throw std::invalid_argument("Dimension mismatch in LowerTriangular-vector multiplication");

    Eigen::VectorXd result(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i; ++j)
            sum += (*this)(i, j) * v(j);
        result(i) = sum;
    }
    return result;
}

LowerTriangular LowerTriangular::operator*(const LowerTriangular& other) const
{
    LowerTriangular result;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = j; k <= i; ++k)
                sum += (*this)(i, k) * other(k, j);
            result(i, j) = sum;
        }
    return result;
}

Eigen::MatrixXd LowerTriangular::operator*(const Eigen::MatrixXd& other) const
{
    if (3 != other.rows())
        throw std::invalid_argument("Matrix dimension mismatch");

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(3, other.cols());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j <= i; ++j)
            result.row(i) += (*this)(i, j) * other.row(j);
    return result;
}

Eigen::MatrixXd SpatialAlgebra::operator*(const Eigen::MatrixXd& lhs, const LowerTriangular& rhs)
{
    if (lhs.cols() != 3)
        throw std::invalid_argument("Matrix dimension mismatch");

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(lhs.rows(), 3);
    for (int i = 0; i < lhs.rows(); ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = j; k < 3; ++k)
                result(i, j) += lhs(i, k) * rhs(k, j);
    return result;
}

LowerTriangular LowerTriangular::inverse() const
{
    double a = data(0,0), b = data(1,0), c = data(1,1);
    double d = data(2,0), e = data(2,1), f = data(2,2);

    // Forward substitution: L·X = I  → solve L·col_j = e_j
    LowerTriangular result;
    double& x00 = result(0,0);  // col 0
    double& x10 = result(1,0);
    double& x20 = result(2,0);
    if (a == 0.0) throw std::runtime_error("Singular lower triangular matrix");
    x00 = 1.0 / a;
    x10 = -b * x00 / c;
    x20 = -(d * x00 + e * x10) / f;

    double& x11 = result(1,1);  // col 1
    double& x21 = result(2,1);
    x11 = 1.0 / c;
    x21 = -e * x11 / f;

    result(2,2) = 1.0 / f;  // col 2

    return result;
}
