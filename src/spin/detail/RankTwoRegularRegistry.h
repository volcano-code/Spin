#pragma once
#include "RankTwoArcRefinement.h"

namespace spin {
// All actual finite edges must have final-registry regular strip proofs.
// This path preserves original arcs/owners/regions/holes, but rejects real
// unbounded arcs and cannot substitute for simultaneous multi-factor fans.
// It grants local curve/chord error coverage, NOT global isotopy or triangles.
class RankTwoRegularRegistryRefinement {
public:
    const RankTwoArcRegistry& registry()const;
    const std::vector<RankTwoRegularEdge>& proofs()const;
    const std::vector<RankTwoArcRefinementRound>& rounds()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoRegularRegistryRefinement(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoRegularRegistryRefinement refineRankTwoRegularRegistry(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,const CGAL::Gmpq&,
        const RankTwoArcRefinementSettings&,const RankTwoSourceLimits&);
};
RankTwoRegularRegistryRefinement refineRankTwoRegularRegistry(
    const RankTwoArrangement&,const RankTwoArcSamplingPlan&,const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings=RankTwoArcRefinementSettings{},
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
