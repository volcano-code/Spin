#pragma once
#include "RankTwoInteriorInventory.h"
#include "RankTwoWholeArrangement.h"

namespace spin {
// Native owner order is fixed by the admitted normal grid. Empty owners are
// retained, not filtered out. This is an inventory, NOT a welded shared mesh.
class RankTwoWholeInteriorInventory {
public:
    const std::string& sourceIdentity()const;
    const std::string& arrangementIdentity()const;
    const std::string& normalIdentity()const;
    const std::vector<RankTwoInteriorInventory>& owners()const;
    std::size_t pointCount()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoWholeInteriorInventory(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoWholeInteriorInventory enumerateRankTwoWholeMovingInteriors(
        const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoInteriorInventorySettings&);
};
// maximumPoints is a TOTAL cap. Candidate/point-proof budgets remain per native
// owner; exceeding any cap rejects the whole composition without partial return.
RankTwoWholeInteriorInventory enumerateRankTwoWholeMovingInteriors(
    const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoInteriorInventorySettings& settings={});
RankTwoWholeInteriorInventory reloadRankTwoWholeMovingInteriors(
    const PreparedRankTwoSource&,const RankTwoArrangement&,const std::string& saved,
    const RankTwoInteriorInventorySettings& settings={});
} // namespace spin
