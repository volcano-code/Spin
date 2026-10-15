#pragma once
#include "RankTwoMovingInterior.h"

namespace spin {
struct RankTwoInteriorInventorySettings {
    std::size_t maximumPoints=1000000;
    std::uint64_t maximumCandidateChecks=1000000;
    RankTwoInteriorPointSettings pointProof;
};
// ALL active native vertical gaps, not a caller-chosen subset. Every point
// retains the adjacent-event/fiber definition and its full native proof.
// This inventory alone grants neither carrier insertion nor global isotopy.
class RankTwoInteriorInventory {
public:
    const std::string& arrangementIdentity()const;
    std::size_t slabCount()const;
    const std::vector<RankTwoMovingInteriorPoint>& points()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoInteriorInventory(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoInteriorInventory enumerateRankTwoMovingInteriors(
        const RankTwoArrangement&,const RankTwoInteriorInventorySettings&);
};
RankTwoInteriorInventory enumerateRankTwoMovingInteriors(const RankTwoArrangement&,
    const RankTwoInteriorInventorySettings& settings={});
RankTwoInteriorInventory reloadRankTwoMovingInteriors(const RankTwoArrangement&,
    const std::string& saved,const RankTwoInteriorInventorySettings& settings={});
} // namespace spin
