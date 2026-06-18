#pragma once

/**
 * @file test_helpers.h
 * @brief Shared test utilities for relative error reporting and edge-case helpers
 * @details Provides macros and utility functions for consistent test output
 *          across all SpatialAlgebra test suites. Includes relative error
 *          reporting (TST-02) alongside traditional EXPECT_NEAR assertions.
 *
 * @see Featherstone, R. (2008). Rigid Body Dynamics Algorithms.
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

namespace test_helpers {

/**
 * @brief Absolute tolerance for near-zero comparisons
 */
constexpr double kAbsTol = 1e-10;

/**
 * @brief Relative tolerance for general comparisons
 */
constexpr double kRelTol = 1e-8;

/**
 * @brief Compute relative error: |computed - expected| / max(|expected|, absTol)
 * @param computed Value from the implementation
 * @param expected Reference value
 * @param absTol Minimum denominator to prevent division by near-zero
 * @return Relative error
 */
inline double relativeError(double computed, double expected,
                             double absTol = kAbsTol) {
    double denom = std::max(std::abs(expected), absTol);
    return std::abs(computed - expected) / denom;
}

/**
 * @brief Report a comparison result with both absolute and relative error
 * @param testName Name of the test or check
 * @param computed Computed value
 * @param expected Expected value
 * @param absTol Absolute tolerance
 * @param relTol Relative tolerance
 * @return true if both absolute and relative errors are within tolerance
 */
inline bool checkNear(double computed, double expected,
                       const std::string& testName = "",
                       double absTol = kAbsTol,
                       double relTol = kRelTol) {
    double absErr = std::abs(computed - expected);
    double relErr = relativeError(computed, expected, absTol);
    bool pass = (absErr < absTol) || (relErr < relTol);

    if (!testName.empty()) {
        std::cout << "  " << (pass ? "✓" : "✗") << " " << testName
                  << "  abs_err=" << std::scientific << std::setprecision(2)
                  << absErr << "  rel_err=" << relErr
                  << std::endl;
    }

    return pass;
}

/**
 * @brief Check a vector component-wise with relative error reporting
 * @param computed Computed vector
 * @param expected Expected vector
 * @param testName Prefix for component names
 * @param absTol Absolute tolerance
 * @param relTol Relative tolerance
 * @return true if all components pass
 */
inline bool checkVector(const Eigen::VectorXd& computed,
                         const Eigen::VectorXd& expected,
                         const std::string& testName = "",
                         double absTol = kAbsTol,
                         double relTol = kRelTol) {
    if (computed.size() != expected.size()) return false;

    bool allPass = true;
    double maxRelErr = 0.0;
    double maxAbsErr = 0.0;

    for (int i = 0; i < computed.size(); ++i) {
        double absErr = std::abs(computed[i] - expected[i]);
        double relErr = relativeError(computed[i], expected[i], absTol);
        maxRelErr = std::max(maxRelErr, relErr);
        maxAbsErr = std::max(maxAbsErr, absErr);
        if (absErr >= absTol && relErr >= relTol) allPass = false;
    }

    if (!testName.empty()) {
        std::cout << "  " << (allPass ? "✓" : "✗") << " " << testName
                  << "  max_abs_err=" << std::scientific << std::setprecision(2)
                  << maxAbsErr << "  max_rel_err=" << maxRelErr
                  << std::endl;
    }

    return allPass;
}

/**
 * @brief GTest assertion that checks both absolute and relative tolerance.
 *        Reports both errors on failure.
 */
#define EXPECT_NEAR_REL(computed, expected, absTol, relTol)                    \
    do {                                                                       \
        double _comp = (computed);                                             \
        double _exp = (expected);                                              \
        double _absErr = std::abs(_comp - _exp);                               \
        double _relErr = relativeError(_comp, _exp, kAbsTol);                  \
        EXPECT_NEAR(_comp, _exp, (absTol))                                     \
            << "  rel_err=" << std::scientific << std::setprecision(2)         \
            << _relErr;                                                        \
    } while (0)

} // namespace test_helpers
