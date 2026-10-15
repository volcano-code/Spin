#pragma once
#include "RankTwoAlgebraic.h"
#include "RankTwoConstruction.h"

namespace spin {
struct RankTwoSharedInterval {
    std::string identity;
    CORE::Expr lowerParameter,upperParameter; // exact open interval in canonical t
    std::array<std::string,2> sourceArcs;
    std::array<bool,2> sourceReversed{}; // t increases opposite to lexicographic arc
    std::array<int,3> originalSigns{};
    std::vector<RankTwoPolynomialOwner> originalOwners;
};
struct RankTwoSharedEdgeRecord;
const std::vector<RankTwoSharedInterval>& sharedIntervals(const RankTwoSharedEdgeRecord&);
} // namespace spin
