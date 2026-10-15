#pragma once
#include "RankTwoArcRegistry.h"

namespace spin {
struct RankTwoExcludedRegion {
    std::size_t registryRegionIndex;
    bool outsideProjection;
    std::vector<std::size_t> negativeNormalHalfplanes;
};
// Exact restriction of the ORIGINAL native DCEL to D>0 and the actual normal
// cone. Empty cells remain present for complement construction. No invented
// truncation coordinates or independent owner sampling. This is not a chord
// isotopy, triangulation, thickness or full carrier certificate.
class RankTwoFinitePartition {
public:
    const std::string& arrangementIdentity()const;
    const std::string& registryIdentity()const;
    const RankTwoArcRegistry& registry()const;
    const std::vector<std::size_t>& regionIndices()const;
    const std::vector<std::size_t>& edgeIndices()const;
    const std::vector<std::size_t>& pointIndices()const;
    const std::vector<RankTwoExcludedRegion>& excludedRegions()const;
    const std::vector<std::string>& excludedArcs()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoFinitePartition(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoFinitePartition buildRankTwoFinitePartition(
        const RankTwoArrangement&,const RankTwoArcRegistry&,const RankTwoSourceLimits&);
};
RankTwoFinitePartition buildRankTwoFinitePartition(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const RankTwoSourceLimits& limits={});
// Rebuild from the original authorities. Untrusted flags, counts and labels
// never grant inclusion, omitted holes, or a positive completion state.
RankTwoFinitePartition reloadRankTwoFinitePartition(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const std::string& saved,const RankTwoSourceLimits& limits={});
} // namespace spin
