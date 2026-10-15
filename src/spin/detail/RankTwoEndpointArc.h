#pragma once
#include "RankTwoArrangementRecords.h"
#include "RankTwoCarrierKernels.h"
#include "RankTwoArcStrip.h"
#include "RankTwoMovingGate.h"

namespace spin {
struct RankTwoEndpointArcCertificate {
    std::string arrangementIdentity,arcIdentity,endpointIdentity,otherEndpointIdentity;
    int projectionAxis=0;
    CGAL::Gmpq fraction;
    std::array<detail::CertifiedAlgebraicReal,2> gateCoordinates; // ORIGINAL chart
    detail::RankTwoRegularStripCertificate strip; // selected (u,v) chart
    std::size_t checkedVertexCount=0;
    std::string identity;
};
// Local endpoint-to-gate arc certificate. The gate's selected projection is
// (1-fraction)*endpoint + fraction*otherEndpoint of the SAME immutable arc.
// projectionAxis=0 uses (x,y), =1 uses (y,x); swapping is an isometry.
// The exact regular strip bounds the whole curve/chord distance, and all
// intervening native arrangement vertices are excluded. Native arc membership
// is independently checked AFTER changing chart; native ordinals are never
// reused as ordinals in the swapped polynomial.
//
// This is NOT a complete incident fan, unrelated-chord exclusion, 3D error
// certificate or complete carrier. Shared endpoints still require common-gate
// fan/adjacency proofs before simultaneous polygonalization.
RankTwoEndpointArcCertificate certifyRankTwoEndpointArc(
    const RankTwoArrangement&,const std::string& arcIdentity,int endpointIndex,
    int projectionAxis,const CGAL::Gmpq& fraction,
    const detail::RankTwoStripBox& chartBox,const CGAL::Gmpq& maximumChordDistance,
    int maximumIsolationBits=512);

struct RankTwoSharedGateArcCertificate {
    std::string arrangementIdentity,arcIdentity,endpointIdentity,gateIdentity;
    int projectionAxis=0;
    std::array<detail::CertifiedAlgebraicReal,2> gateCoordinates;
    detail::RankTwoRegularStripCertificate strip;
    std::size_t checkedVertexCount=0;
    std::string belowRegion,aboveRegion,identity;
};
// Same local arc proof with a shared moving gate; also derives the native
// regions below/above the graph in the selected chart. It accepts an unbounded
// far end only when this finite gate is certified inside the ORIGINAL arc.
RankTwoSharedGateArcCertificate certifyRankTwoEndpointArc(
    const RankTwoArrangement&,const std::string& arcIdentity,const RankTwoMovingGate&,
    const detail::RankTwoStripBox&,const CGAL::Gmpq& maximumChordDistance,
    int maximumIsolationBits=512,
    const std::vector<detail::RankTwoStripPolynomial>& obstacles={});
} // namespace spin
