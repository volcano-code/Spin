#pragma once

#include "spin/QuadricCoefficients.h"

namespace spin {

/** Literal-monomial unit sphere s=Dp/||Dp||; never a native geometry quadric.
 * D=diag(1,1,1,1,1,1,1,2,2,2), dual to independent raw moments.
 */
class DualMonomialCoordinates10 {
public:
    using Vector = QuadricCoefficients10::Vector;
    using Matrix = Eigen::Matrix<double,10,10>;
    static DualMonomialCoordinates10 fromGeometry(const QuadricCoefficients10& p);
    static DualMonomialCoordinates10 fromUnnormalized(const Vector& s);
    const Vector& values() const noexcept { return values_; }
    QuadricCoefficients10 geometry() const;
    /** dp/ds=(I-pp^T)D^-1/||D^-1 s||, including normalization. */
    Matrix geometryJacobian() const;
    Vector projectToTangent(const Vector& direction) const;
    DualMonomialCoordinates10 geodesicUpdate(const Vector& direction, double arcLength) const;
private:
    explicit DualMonomialCoordinates10(const Vector& s) : values_(s) {}
    Vector values_;
};

} // namespace spin
