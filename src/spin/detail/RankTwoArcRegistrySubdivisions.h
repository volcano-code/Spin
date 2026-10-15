#pragma once
#include "RankTwoArcRegistryNeighborhoods.h"

namespace spin {
// A topologically sorted moving-point DAG, normalized by stable dependency ID.
// Parents must be actual published points on the SAME immutable native arc.
// Only the native monotone coordinate is interpolated; the other coordinate is
// evaluated on its exact original fiber branch, NEVER interpolated on a chord.
// Requests are coordinate-free; unknown parents, duplicate dependencies and
// subdivisions inside a proof-bound endpoint chord are rejected. Every face
// consumes the updated common chain. New samples require NEW strip/topology
// proofs; this API grants none. Empty input preserves historical v1/v2 bytes.
RankTwoArcRegistry buildRankTwoArcRegistryWithSubdivisions(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointNeighborhood>&,
    const std::vector<RankTwoArcSubdivision>&,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
