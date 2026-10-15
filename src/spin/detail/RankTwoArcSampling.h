#pragma once
#include "RankTwoArrangementRecords.h"
#include "RankTwoCarrierKernels.h"

namespace spin {
struct RankTwoArcSample {
    std::string arrangementIdentity,arcIdentity;
    std::array<std::string,2> endpointIdentities;
    CGAL::Gmpq fraction;
    int interiorBranch=-1;
    bool vertical=false;
    std::array<detail::CertifiedAlgebraicReal,2> coordinates;
    std::string identity;
};
// A point on a finite arc of an already certified immutable arrangement.
// This does not certify a replacement chord, endpoint sector or triangulation.
RankTwoArcSample sampleRankTwoArc(const RankTwoArrangement&,const std::string& arcIdentity,
    const CGAL::Gmpq& fraction,int maximumIsolationBits=512);
namespace detail {
// Internal exact fiber bridge. The registry owns/validates the parent DAG;
// this function returns coordinates only, never a published source feature.
// Projection must lie STRICTLY inside the original native monotone arc.
std::array<CertifiedAlgebraicReal,2> evaluateRankTwoArcProjection(
    const RankTwoArrangement&,const std::string& arcIdentity,
    const CertifiedAlgebraicReal& projection,int maximumIsolationBits=512);
}
} // namespace spin
