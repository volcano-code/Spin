#pragma once
#include "RankTwoNormalOverlay.h"
#include "RankTwoAlgebraic.h"

namespace spin {
struct RankTwoNormalBoundaryEvent {
    std::array<std::size_t,2> edge{}; // canonical endpoint order
    std::size_t ordinal=0; // among all F/h/D events on the open positive cone edge
    std::ptrdiff_t radiusOrdinal=-1; // among F roots ONLY; matches old edge certificates
    std::vector<std::string> zeroPolynomials;
    std::string identity;
    std::array<detail::AlgebraicRootCertificate,2> point;
    CORE::Expr parameter; // exact t in (1-t)*n_low + t*n_high, never a double
};
// Read-only events recovered from the actual overlay's finite vertices.
// Endpoint normal vertices have their own typed registry (RankTwoNormalVertices),
// never an ordinal in this open-edge list.
const std::vector<RankTwoNormalBoundaryEvent>& normalBoundaryEvents(const RankTwoArrangement&);
} // namespace spin
