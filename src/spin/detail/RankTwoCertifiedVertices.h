#pragma once
#include "RankTwoConstruction.h"
#include "RankTwoAlgebraic.h"
#include <array>

namespace spin {
// Borrow only values minted by the current local DCEL factory. No caller
// certificate, predicate result or foreign arrangement is accepted. The
// reference lives as long as its immutable arrangement owner.
const std::array<detail::CertifiedAlgebraicReal,2>& certifiedRankTwoVertex(
    const RankTwoArrangement&,const std::string& vertexIdentity);
} // namespace spin
