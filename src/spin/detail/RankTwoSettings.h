#pragma once
#include "spin/QuadricTriangulator.h"
#include "RankTwoWholeArrangement.h"

namespace spin {
// Strict settings-only admission. Returns a complete algebraic partition,
// NOT a finite carrier, closed source mesh or an accepted material object.
RankTwoArrangement analyzeRankTwoArrangement(
    const PreparedRankTwoSource&, const QuadricTriangulationSettings&);
}
