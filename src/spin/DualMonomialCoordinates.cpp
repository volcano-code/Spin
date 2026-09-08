#include "spin/DualMonomialCoordinates.h"
#include <cmath>
#include <stdexcept>

namespace spin {

DualMonomialCoordinates10 DualMonomialCoordinates10::fromUnnormalized(const Vector& s) {
    const double norm=s.norm();
    if(!s.allFinite() || !std::isfinite(norm) || norm<=0)
        throw std::invalid_argument("dual monomial coordinates must be finite and nonzero");
    return DualMonomialCoordinates10(s/norm);
}

DualMonomialCoordinates10 DualMonomialCoordinates10::fromGeometry(const QuadricCoefficients10& p) {
    Vector s=p.values(); s.tail<3>()*=2;
    return fromUnnormalized(s);
}

QuadricCoefficients10 DualMonomialCoordinates10::geometry() const {
    Vector p=values_; p.tail<3>()*=.5;
    return QuadricCoefficients10::fromUnnormalized(p);
}

DualMonomialCoordinates10::Matrix DualMonomialCoordinates10::geometryJacobian() const {
    Vector ts=values_; ts.tail<3>()*=.5;
    const Vector p=ts.normalized();
    Matrix t=Matrix::Identity(); t.diagonal().tail<3>().setConstant(.5);
    return (Matrix::Identity()-p*p.transpose())*t/ts.norm();
}

DualMonomialCoordinates10::Vector DualMonomialCoordinates10::projectToTangent(const Vector& d) const {
    if(!d.allFinite()) throw std::invalid_argument("dual tangent direction must be finite");
    return d-values_*values_.dot(d);
}

DualMonomialCoordinates10 DualMonomialCoordinates10::geodesicUpdate(const Vector& d, double arcLength) const {
    if(!std::isfinite(arcLength) || arcLength<0)
        throw std::invalid_argument("dual geodesic arc must be finite and nonnegative");
    const Vector tangent=projectToTangent(d);
    const double norm=tangent.norm();
    if(!std::isfinite(norm) || norm<=0)
        throw std::invalid_argument("dual geodesic tangent must be nonzero and finite");
    return fromUnnormalized(std::cos(arcLength)*values_+std::sin(arcLength)*tangent/norm);
}

} // namespace spin
