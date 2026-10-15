#pragma once
#include "RankTwoMixedArcRefinement.h"

namespace spin {
struct RankTwoMixedArcInterval {
    std::string arcIdentity;
    std::array<std::string,2> pointIdentities;
};
// A canonical, non-overlapping set of ACTUAL native intervals refined in ONE
// common registry. Every selected ordinary proof is rebuilt after any arc
// changes. Endpoint chords keep their complete neighborhood authorities.
// No claim about other edges, global neighborhood compatibility/isotopy,
// clipping, triangulation or material is granted by this local certificate.
class RankTwoMixedArcIntervals {
public:
    const std::vector<RankTwoMixedArcInterval>& intervals()const;
    const RankTwoArcRegistry& registry()const;
    const std::vector<RankTwoRegularEdge>& proofs()const;
    const std::vector<std::string>& endpointEdges()const;
    const std::vector<RankTwoMixedArcRefinementRound>& rounds()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoMixedArcIntervals(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoMixedArcIntervals refineRankTwoMixedArcIntervalsFromRegistry(
        const RankTwoArrangement&,const RankTwoArcRegistry&,const RankTwoArcSamplingPlan&,
        const std::vector<RankTwoEndpointNeighborhood>&,
        const std::vector<RankTwoMixedArcInterval>&,const CGAL::Gmpq&,
        const RankTwoArcRefinementSettings&,const RankTwoSourceLimits&);
    friend RankTwoMixedArcIntervals refineRankTwoMixedArcIntervals(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
        const std::vector<RankTwoEndpointNeighborhood>&,
        const std::vector<RankTwoMixedArcInterval>&,const CGAL::Gmpq&,
        const RankTwoArcRefinementSettings&,const RankTwoSourceLimits&);
};
// Intervals are sorted canonically. Point IDs must exist in increasing native
// order in the initial common registry; duplicates/edge overlaps are errors.
// Intervals sharing only a vertex are permitted. Unbounded arcs are retained
// elsewhere, not silently clipped or promoted to finite interval coverage.
RankTwoMixedArcIntervals refineRankTwoMixedArcIntervals(
    const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
    const std::vector<RankTwoEndpointNeighborhood>&,
    const std::vector<RankTwoMixedArcInterval>&,
    const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings={},const RankTwoSourceLimits& limits={});
// Same refinement kernel, starting from an actual sealed registry. Every
// round retains its primitive plan, endpoint proofs and common seam authority.
// Selected middle intervals do NOT grant coverage of their omitted endpoints.
RankTwoMixedArcIntervals refineRankTwoMixedArcIntervalsFromRegistry(
    const RankTwoArrangement&,const RankTwoArcRegistry&,const RankTwoArcSamplingPlan&,
    const std::vector<RankTwoEndpointNeighborhood>&,
    const std::vector<RankTwoMixedArcInterval>&,const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings={},const RankTwoSourceLimits& limits={});
} // namespace spin
