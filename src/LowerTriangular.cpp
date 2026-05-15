#include "LowerTriangular.h"
#include <stdexcept>
#include <cmath>

namespace SpatialAlgebra {

Eigen::VectorXd LowerTriangular::operator*(const Eigen::VectorXd &v) const
{
    if (n != v.size())
        throw std::invalid_argument("Dimension mismatch in LowerTriangular-vector multiplication");
    
    Eigen::VectorXd result(n);
    for (int i = 0; i < n; ++i)
    {
        result(i) = 0.0;
        for (int j = 0; j <= i; ++j)
        {
            result(i) += (*this)(i, j) * v(j);
        }
    }
    return result;
}

LowerTriangular LowerTriangular::inverse() const
{
    LowerTriangular result(n);
    
    // Check for singularity and compute diagonal elements
    for (int i = 0; i < n; ++i)
    {
        double diag = (*this)(i, i);
        if (std::abs(diag) < 1e-15)
            throw std::runtime_error("Matrix is singular (zero diagonal element at index " + std::to_string(i) + ")");
        result(i, i) = 1.0 / diag;
    }
    
    // Compute off-diagonal elements using forward substitution
    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            double sum = 0.0;
            for (int k = j; k < i; ++k)
            {
                sum += (*this)(i, k) * result(k, j);
            }
            result(i, j) = -sum / (*this)(i, i);
        }
    }
    
    return result;
}

} // namespace SpatialAlgebra
