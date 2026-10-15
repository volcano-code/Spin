#pragma once
#include "RankTwoConstruction.h"

namespace spin {
// Exact stored normal rays, not independently re-normalized unit vectors.
// This is ONE cone overlay, not certification of a whole normal sphere mesh.
struct RankTwoNormalFace {
    std::array<std::size_t,3> vertexIds{};
    std::array<std::array<CGAL::Gmpq,3>,3> rays{};
};

RankTwoArrangement analyzeRankTwoNormalFace(
    const PreparedRankTwoSource&,const RankTwoArrangement& global,
    const RankTwoNormalFace&,const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
