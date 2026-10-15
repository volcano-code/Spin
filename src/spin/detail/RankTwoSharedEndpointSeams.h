#pragma once
#include "RankTwoArcRegistrySeams.h"
#include "RankTwoEndpointNeighborhood.h"

namespace spin {
struct RankTwoSharedEndpointBinding {
    std::size_t face;
    std::string intervalIdentity,arcIdentity,endpointIdentity;
    std::string neighborhoodIdentity,proofIdentity,gateIdentity,anchorIdentity;
    int projectionAxis;
    CGAL::Gmpq fraction;
    std::size_t pointIndex;
};
struct RankTwoSharedEndpointPoint {
    std::string identity,intervalIdentity,endpointIdentity;
    std::array<detail::CertifiedAlgebraicReal,2> coordinates;
    std::vector<std::size_t> consumers;
    bool masterAlias=false; // both definitions structurally equal one original moving master sample
    std::string definitionPairIdentity;
};
// Both factory-certified definitions are retained. Their exact CURRENT point
// equality is checked, not tolerance-welded. Source motion must later prove
// the two definitions have equal velocities; this class grants no derivative.
class RankTwoSharedEndpointSeams {
public:
    const RankTwoSharedFiniteSeams& master()const;
    const std::string& masterIdentity()const;
    const std::vector<RankTwoSharedEndpointPoint>& points()const;
    const std::vector<RankTwoSharedEndpointBinding>& bindings()const;
    const std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>>& neighborhoods()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoSharedEndpointSeams(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoSharedEndpointSeams publishRankTwoSharedEndpointSeams(
        const PreparedRankTwoSource&,const RankTwoSharedFiniteSeams&,
        const std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>>&,
        const RankTwoSourceLimits&);
};
RankTwoSharedEndpointSeams publishRankTwoSharedEndpointSeams(
    const PreparedRankTwoSource&,const RankTwoSharedFiniteSeams&,
    const std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>>&,
    const RankTwoSourceLimits& limits={});
RankTwoArcRegistry buildRankTwoArcRegistryWithSharedEndpointSeams(
    const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
    const RankTwoSharedEndpointSeams&,std::size_t face,const RankTwoSourceLimits& limits={});
} // namespace spin
