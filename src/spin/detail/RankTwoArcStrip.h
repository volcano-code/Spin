#pragma once
#include <CGAL/Gmpq.h>
#include <array>
#include <string>
#include <vector>
#include <variant>

namespace spin::detail {
// Coefficient x^i y^j at 5*i+j, total degree at most four.
using RankTwoStripPolynomial=std::array<CGAL::Gmpq,25>;
using RankTwoStripBox=std::array<CGAL::Gmpq,4>; // xlo,xhi,ylo,yhi
// Decode immutable native factor text into the requested exact projection.
// This is only a rational representation, not a geometric certificate.
RankTwoStripPolynomial parseRankTwoStripPolynomial(const std::string&,int projectionAxis);
struct RankTwoRegularStripCertificate {
    std::string inputIdentity;
    int yDerivativeSign=0;
    CGAL::Gmpq derivativeLowerBound,slopeUpperBound,curvatureUpperBound;
    CGAL::Gmpq chordDistanceUpperBound;
    std::vector<int> obstacleSigns;
};
// Internal search outcome: the string is inability to certify, never an
// infeasibility certificate. Input/resource/environment errors still throw.
using RankTwoRegularStripAttempt=std::variant<RankTwoRegularStripCertificate,std::string>;
RankTwoRegularStripAttempt tryCertifyRankTwoRegularStrip(
    const RankTwoStripPolynomial&,const RankTwoStripBox&,
    const std::vector<RankTwoStripPolynomial>&,const CGAL::Gmpq&);
// Sufficient regular-graph certificate, computed entirely in exact rationals.
// Proves a UNIQUE y=f(x) inside the box and a bidirectional 2D chord-distance
// bound for any two actual graph points in it. Also excludes the supplied
// obstacles from the WHOLE box. The caller must bind the graph/points to its
// immutable arrangement and provide all relevant obstacles; neither omitted
// owners, endpoint sectors, 3D lifts nor complete-carrier topology are proved.
// Failure to bound is ARRANGEMENT_UNRESOLVED, never infeasibility or permission
// to discard the region. Subdivision/other charts remain the caller's job.
RankTwoRegularStripCertificate certifyRankTwoRegularStrip(
    const RankTwoStripPolynomial&,const RankTwoStripBox&,
    const std::vector<RankTwoStripPolynomial>& obstacles,
    const CGAL::Gmpq& maximumChordDistance);
} // namespace spin::detail
