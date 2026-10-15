#pragma once
#include "RankTwoArcRegistryFans.h"
#include "RankTwoEndpointNeighborhood.h"

namespace spin {
// Bind each factory-owned COMPLETE local star to the actual shared endpoint
// edges. Both graph sectors AND straight boundary rays are required. A sample
// inside a proof-bound chord invalidates it; unbounded incident arcs require
// certified clipping and cannot be silently dropped. This is not a global
// chord embedding, isotopy or triangulation grant. Empty input preserves v1.
RankTwoArcRegistry buildRankTwoArcRegistryWithNeighborhoods(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointNeighborhood>&,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
