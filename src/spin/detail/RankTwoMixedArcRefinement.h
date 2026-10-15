#pragma once
#include "RankTwoArcRegistrySubdivisions.h"
#include "RankTwoArcRefinement.h"

namespace spin {
struct RankTwoMixedArcRefinementRound {
    std::string registryIdentity;
    std::size_t selectedEdges=0,endpointEdges=0,regularEdges=0,splitEdges=0,boxTrials=0;
};
// An explicitly selected native interval covered WITHOUT A GAP by original
// complete endpoint chords and final-registry regular strips. Subdivision
// uses actual gate-relative curve points, never rationalized gate coordinates.
// This is a bridge toward compatible global coverage, NOT a global isotopy,
// clipping, triangulation or material certificate.
class RankTwoMixedArcRefinement {
public:
    const std::string& arcIdentity()const;
    const std::array<std::string,2>& endpointIdentities()const;
    const RankTwoArcRegistry& registry()const;
    const std::vector<RankTwoRegularEdge>& proofs()const;
    const std::vector<std::string>& endpointEdges()const;
    const std::vector<RankTwoMixedArcRefinementRound>& rounds()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoMixedArcRefinement(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoMixedArcRefinement refineRankTwoMixedArcInterval(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
        const std::vector<RankTwoEndpointNeighborhood>&,const std::string&,
        const std::array<std::string,2>&,const CGAL::Gmpq&,
        const RankTwoArcRefinementSettings&,const RankTwoSourceLimits&);
};
// Both point identities must occur in increasing ORIGINAL native order in
// the initial common registry. All initial-plan/native owners remain present;
// only unproved ordinary edges within this explicit interval are subdivided.
// Changed partitions invalidate EVERY regular proof, including accepted ones.
RankTwoMixedArcRefinement refineRankTwoMixedArcInterval(
    const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
    const std::vector<RankTwoEndpointNeighborhood>&,const std::string& arcIdentity,
    const std::array<std::string,2>& intervalPointIdentities,
    const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings={},const RankTwoSourceLimits& limits={});
} // namespace spin
