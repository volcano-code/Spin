#pragma once

#include "spin/QuadricMomentEvaluator.h"
#include <stdexcept>

namespace spin::detail {

enum class QuadricExportRole { VariableMaterial, Cavity };

// This is a rejection gate, not a geometry/fabrication certification API.
// Shared by the two production export paths so each can be exercised directly.
inline void requireNonemptyQuadricExportRegion(
    const QuadricMomentEvaluation& evaluation, QuadricExportRole role)
{
    if (evaluation.materialRegion != QuadricMaterialRegion::CertifiedEmpty)
        return;
    if (role == QuadricExportRole::VariableMaterial)
        throw std::runtime_error("certified empty variable region is moment-only; fixed-shell geometry assembly is not an export certificate");
    throw std::runtime_error("certified empty cavity region requires an explicit empty-artifact export protocol");
}

} // namespace spin::detail
