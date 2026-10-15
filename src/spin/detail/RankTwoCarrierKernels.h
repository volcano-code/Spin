#pragma once
#include "RankTwoAlgebraic.h"
#include <CGAL/Gmpq.h>
#include <array>

namespace spin::detail {
// The identical representation is an exact equality proof, not proximity.
// Avoid constructing a prohibitively high-degree SubRep(x,x) at a seam.
CORE::Expr rankTwoExactLiftDifference(const CORE::Expr& upper,const CORE::Expr& lower);
// Exact local conditions for a piecewise-affine vertical prism. Edge i is
// (i,(i+1)%3). Boundary classification MUST come from the certified material
// complex; this function does not certify that classification or topology.
void certifyRankTwoTriangleThickness(const std::array<CORE::Expr,3>& lower,
    const std::array<CORE::Expr,3>& upper,
    const std::array<bool,3>& trueMaterialBoundary);

// Exact (1-t)*a+t*b, 0<t<1. This is a numeric algebraic-coordinate bridge,
// not the carrier's provenance DAG: callers must retain endpoint feature IDs
// and t even when a==b. Never average two isolator endpoints as a source point.
CertifiedAlgebraicReal sampleRankTwoAlgebraicCoordinate(
    const CertifiedAlgebraicReal& a,const CertifiedAlgebraicReal& b,
    const CGAL::Gmpq& t,int maximumIsolationBits=512);
// Bounded signed affine ruler, -1<t<1, t!=0. Negative weights are for a
// separately identified reflected endpoint gate, NOT interior arc samples.
// Both signs use outward intervals for the actual selected algebraic roots.
CertifiedAlgebraicReal affineRankTwoAlgebraicCoordinate(
    const CertifiedAlgebraicReal& a,const CertifiedAlgebraicReal& b,
    const CGAL::Gmpq& t,int maximumIsolationBits=512);
AlgebraicReal rankTwoAlgebraicHandle(const CertifiedAlgebraicReal&);
} // namespace spin::detail
