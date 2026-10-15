#pragma once
#include "RankTwoArcRegistry.h"
#include "RankTwoEndpointFan.h"

namespace spin {
// Adds every gate point from each private, fully certified incident fan to
// the shared registry, and binds its EXACT event-to-gate chord to the proof.
// Ordinary samples may not split that chord: independently sampled incident
// curves would lose the simultaneous common-gate ordering guarantee.
// This covers endpoint edges only, not the remaining regular arc segments.
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
