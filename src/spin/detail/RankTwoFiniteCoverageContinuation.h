#pragma once
#include "RankTwoFiniteArcCoverage.h"

namespace spin {
// Whole local finite coverage starts from THIS actual factory-owned registry.
// Private seam/endpoint authorities and prior moving dependencies survive all
// rounds. No global compatibility/isotopy/precision/thickness/3D grant follows.
RankTwoFiniteArcCoverage refineRankTwoFiniteArcCoverageFromRegistry(
    const RankTwoArrangement&,const RankTwoArcRegistry&,const RankTwoArcSamplingPlan&,
    const std::vector<RankTwoEndpointNeighborhood>&,const CGAL::Gmpq& maximumChordDistance,
    const RankTwoArcRefinementSettings& settings={},const RankTwoSourceLimits& limits={});
} // namespace spin
