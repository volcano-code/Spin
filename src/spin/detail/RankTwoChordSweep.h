#pragma once
#include "RankTwoChordIncidence.h"
#define RANK_TWO_FINITE_PARTITION_INCIDENCE_V1 1

namespace spin {
struct RankTwoChordSweepSettings {
    std::uint64_t maximumPredicateGroups=1000000;
};
// Exact-coordinate sweep over ALL actual finite points/chords. Sorted closed
// interval windows prove disjointness of excluded pairs; remaining pairs use
// the original exact incidence predicates. Counts describe total coverage,
// while predicateGroupCount records actual work under the explicit budget.
// It does not certify native-curve isotopy, infinite geometry or a carrier.
RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const RankTwoChordSweepSettings& settings=RankTwoChordSweepSettings{});
// ALL and ONLY points/edges of an immutable native finite partition. The
// caller cannot choose a subset to hide a crossing. Full-registry APIs remain
// unchanged; excluded finite arcs are not granted incidence or isotopy.
RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement&,
    const RankTwoFinitePartition&,const RankTwoChordSweepSettings& settings=RankTwoChordSweepSettings{});
} // namespace spin
