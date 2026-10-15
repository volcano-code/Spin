#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/QuadricTriangulator.h"

#include <Eigen/Core>

#include <string>

namespace spin {

/**
 * First-order coefficient constraints for a fixed Euclidean quadric stratum.
 *
 * The rows are independent linear constraints on a perturbation of the ten
 * normalized polynomial coefficients.  The spherical normalization row is
 * included explicitly.  For rank-deficient strata this is intentionally a
 * lower-dimensional tangent space; it is not a fabricated 10-D Jacobian.
 */
struct QuadricStratumTangent {
    bool defined = false;
    QuadricSurfaceSubtype subtype = QuadricSurfaceSubtype::Plane;
    Eigen::MatrixXd constraintRows;
    /** Orthonormal columns spanning the represented fixed-stratum tangent.
     * This is a reduced coefficient chart; it must not be padded to ten
     * directions when the stratum is rank deficient. */
    Eigen::MatrixXd basis;
    int tangentDimension = 0;
    std::string diagnostic;

    bool contains(const QuadricCoefficients10::Vector& direction,
                  double tolerance = 1e-10) const;
};

QuadricStratumTangent evaluateQuadricStratumTangent(
    const QuadricCoefficients10& quadric);

} // namespace spin
