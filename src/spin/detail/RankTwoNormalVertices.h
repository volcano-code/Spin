#pragma once
#include "RankTwoNormalOverlay.h"
#include <optional>

namespace spin {
enum class RankTwoNormalVertexLocation { FinitePositive,Ideal,Negative };
struct RankTwoNormalVertexRecord {
    std::size_t normalVertex=0;
    std::string identity;
    RankTwoNormalVertexLocation location=RankTwoNormalVertexLocation::Ideal;
    // A fixed rational ray has a rational preimage in the rank-two graph.
    // No fake finite point is attached to ideal or negative rays.
    std::optional<std::array<CGAL::Gmpq,2>> point;
    std::array<int,3> originalSigns{}; // F,h,D; meaningful only when point exists
};
const std::vector<RankTwoNormalVertexRecord>& normalVertexRecords(const RankTwoArrangement&);
} // namespace spin
