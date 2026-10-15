#pragma once
#include "RankTwoArcRegistrySubdivisions.h"

namespace spin {
// Rebuild the SAME publication with an extended moving-point DAG. The native
// arrangement, initial sampling plan, endpoint neighborhoods and immutable
// seam import remain authoritative. A local request cannot insert new seam
// points: that requires a new common master consumed by both owners.
// Existing subdivisions cannot disappear. This grants no strip/isotopy proof.
RankTwoArcRegistry rebuildRankTwoArcRegistryWithSubdivisions(
    const RankTwoArrangement&,const RankTwoArcRegistry& authority,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointNeighborhood>&,
    const std::vector<RankTwoArcSubdivision>&,
    const RankTwoSourceLimits& limits={});
} // namespace spin
