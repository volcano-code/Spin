#pragma once
#include "RankTwoRegistryEdge.h"

namespace spin {
struct RankTwoArcRefinementSettings {
    CGAL::Gmpq initialBoxGridStep=CGAL::Gmpq(1,16);
    int maximumBoxHalvings=16;
};
struct RankTwoArcRefinementRound {
    std::string registryIdentity;
    std::size_t selectedEdges=0,certifiedEdges=0,splitEdges=0,boxTrials=0;
};
class RankTwoRegularArcRefinement {
public:
    const std::string& arcIdentity()const;
    const RankTwoArcRegistry& registry()const;
    const std::vector<RankTwoRegularEdge>& proofs()const;
    const std::vector<RankTwoArcRefinementRound>& rounds()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoRegularArcRefinement(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoRegularArcRefinement refineRankTwoArcInterval(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::string&,const CGAL::Gmpq&,const CGAL::Gmpq&,
        const CGAL::Gmpq&,const RankTwoArcRefinementSettings&,const RankTwoSourceLimits&);
};
// Refine an explicitly bounded INTERIOR interval of one actual finite arc.
// The supplied complete shared plan must already contain both fractions.
// Failed ordinary chords are bisected in the same native moving coordinate.
// Every partition change invalidates and rebuilds ALL selected proofs.
// Boxes are chosen with exact value floors, never root-interval midpoints.
// This does NOT cover endpoint fans, other arcs, mutual chord topology or a
// full carrier. Its interval scope remains explicit in the saved attachment.
RankTwoRegularArcRefinement refineRankTwoArcInterval(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::string& arcIdentity,
    const CGAL::Gmpq& beginFraction,const CGAL::Gmpq& endFraction,
    const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings=RankTwoArcRefinementSettings{},
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
