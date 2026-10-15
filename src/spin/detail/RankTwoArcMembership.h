#pragma once
#include "RankTwoArcSampling.h"

namespace spin::detail {
// Exact finite native-arc membership. Unlike a zero-set query, this checks
// the ORIGINAL interior fiber ordinal and native monotone endpoint bounds.
// Correlated Cartesian coordinate roots are not re-eliminated as independent
// variables. No approximation, tolerance or geometric grant is introduced.
bool isOnRankTwoNativeArc(const RankTwoArrangement&,const std::string& arcIdentity,
    const std::array<CertifiedAlgebraicReal,2>& coordinates);
} // namespace spin::detail
