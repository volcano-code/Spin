#pragma once
#include "RankTwoFinitePartition.h"
#include "RankTwoMixedArcIntervals.h"

namespace spin {
// Complete LOCAL approximation coverage of the owned finite native domain.
// Includes Empty/Cut/Full cells and every outer/inner ring. The factory, not
// its caller, selects the arcs. All witnesses consume ONE final registry.
// This is necessary but NOT sufficient for global neighborhood compatibility,
// isotopy, profile precision, strict thickness or a three-dimensional carrier.
class RankTwoFiniteArcCoverage {
public:
    const std::string& arrangementIdentity()const;
    const RankTwoFinitePartition& partition()const;
    const RankTwoMixedArcIntervals& refinement()const;
    std::size_t edgeCount()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoFiniteArcCoverage(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoFiniteArcCoverage certifyRankTwoFiniteArcCoverage(
        const RankTwoArrangement&,const RankTwoFinitePartition&,const RankTwoMixedArcIntervals&);
};
// A selected-interval result cannot be promoted by a boolean or a count:
// compare actual whole native arcs and final shared edges to the partition.
RankTwoFiniteArcCoverage certifyRankTwoFiniteArcCoverage(
    const RankTwoArrangement&,const RankTwoFinitePartition&,const RankTwoMixedArcIntervals&);
RankTwoFiniteArcCoverage refineRankTwoFiniteArcCoverage(
    const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
    const std::vector<RankTwoEndpointNeighborhood>&,const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings={},const RankTwoSourceLimits& limits={});
} // namespace spin
