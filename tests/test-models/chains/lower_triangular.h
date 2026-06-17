#pragma once

/**
 * @file lower_triangular.h
 * @brief Test model factories for lower triangular matrix tests
 * @details Defines factory functions that produce RobotModel instances with
 *          specific 3×3 inertia matrices for testing LowerTriangular pack/
 *          unpack operations. Each function constructs a single-link chain
 *          whose inertia field encodes a known lower-triangular matrix.
 *
 *          The LowerTriangular class uses packed storage (1D array with
 *          index mapping: idx = i*(i+1)/2 + j) for memory efficiency.
 *          These models define 3×3 dense matrices that the adapter
 *          converts to packed format for comparison.
 *
 * @see LowerTriangular (include/LowerTriangular.h)
 */

#include <Eigen/Dense>

#include "robot_model.h"

namespace test_models {

/**
 * @brief Identity 3×3 lower triangular inertia for LT identity tests
 * @details Constructs a single-link chain with an identity inertia matrix.
 *          In lower-triangular packed form, this is [1, 0, 1, 0, 0, 1].
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=origin,
 *            inertia=I_3×3, name="link0"
 *
 *          Test usage: verify that pack/unpack round-trips produce the
 *          identity, and that the storage layout is correct.
 *
 * @return RobotModel containing a 1-link chain with identity inertia
 */
inline RobotModel makeIdentityLT()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = Eigen::Matrix3d::Identity();
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Scaled identity inertia for LT scalar operation tests
 * @details Constructs a single-link chain with a diagonal inertia matrix
 *          where all diagonal entries equal the parameter value.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=origin,
 *            inertia=diag(value,value,value), name="link0"
 *
 *          Test usage: verify that scalar operations (multiply, divide) on
 *          LT matrices preserve the diagonal structure and that packed
 *          storage correctly handles uniform diagonal values.
 *
 * @param value Diagonal value (default 2.0). All three diagonal entries
 *              are set to this value.
 *
 * @return RobotModel containing a 1-link chain with scaled identity inertia
 */
inline RobotModel makeDiagonalLT(double value = 2.0)
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();
    link.inertia = value * Eigen::Matrix3d::Identity();
    link.name = "link0";
    model.joints.push_back(link);
    return model;
}

/**
 * @brief Non-trivial lower triangular inertia for LT pack/unpack tests
 * @details Constructs a single-link chain with a dense, non-symmetric 3×3
 *          inertia matrix representing lower-triangular packed values:
 *          @f$ [[1, 0, 0], [2, 3, 0], [4, 5, 6]] @f$
 *
 *          The dense matrix is deliberately non-symmetric — the adapter
 *          symmetrizes when converting to packed LT form (only the lower
 *          triangle is stored). The upper triangle entries are structurally
 *          zero in packed representation.
 *
 *          Chain layout:
 *          - Link 0: root link, Z-revolute, mass=1.0, COM=origin,
 *            inertia=[[1,0,0],[2,3,0],[4,5,6]], name="link0"
 *
 *          Test usage: verify pack/unpack correctness for a matrix with
 *          distinct lower-triangular values. The packed representation
 *          should be [1, 2, 3, 4, 5, 6] (column-major lower triangle).
 *
 * @return RobotModel containing a 1-link chain with arbitrary LT inertia
 */
inline RobotModel makeArbitraryLT()
{
    RobotModel model;
    JointSpec link;
    link.parent = -1;
    link.parentToJoint = Eigen::Matrix4d::Identity();
    link.jointAxis = Eigen::Vector3d::UnitZ();
    link.type = JointType::REVOLUTE;
    link.mass = 1.0;
    link.com = Eigen::Vector3d::Zero();

    // Non-trivial lower triangle: row-major [[1,0,0],[2,3,0],[4,5,6]]
    Eigen::Matrix3d M;
    M << 1.0, 0.0, 0.0,
         2.0, 3.0, 0.0,
         4.0, 5.0, 6.0;
    link.inertia = M;
    link.name = "link0";

    model.joints.push_back(link);
    return model;
}

} // namespace test_models
