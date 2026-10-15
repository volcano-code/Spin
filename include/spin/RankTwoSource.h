#pragma once

#include "spin/QuadricCoefficients.h"
#include <CGAL/Gmpq.h>
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

namespace spin {

// Same ordering as QuadricCoefficients10: c,bx,by,bz,xx,yy,zz,yz,xz,xy.
// Off-diagonal entries are matrix entries (polynomial cross coefficient / 2).
using RankTwoExactCoefficients = std::array<CGAL::Gmpq,10>;

struct RankTwoSourceLimits {
    int maximumRefinementDepth=24;
    int maximumIsolationBits=512;
    std::size_t maximumArrangementCells=200000;
    std::size_t maximumVertices=1000000;
    double maximumDeviationOverInputScale=0x1p-20;
    double maximumCoordinateErrorOverInputScale=0x1p-44;
};

class RankTwoSourceError : public std::invalid_argument {
public:
    RankTwoSourceError(std::string code,std::string witness);
    const std::string& code() const noexcept { return code_; }
    const std::string& witness() const noexcept { return witness_; }
private:
    std::string code_,witness_;
};

// Factory-owned, exposed as const. All entries refer to the authority polynomial.
struct RankTwoExactGraph {
    std::array<CGAL::Gmpq,3> kernel{};
    std::array<std::array<CGAL::Gmpq,2>,3> plane{};
    CGAL::Gmpq beta{};
    // h(s,t): 1,s,t,s^2,st,t^2 (st stores the full polynomial coefficient).
    std::array<CGAL::Gmpq,6> height{};
    // G(s,t): per component coefficients of 1,s,t.
    std::array<std::array<CGAL::Gmpq,3>,3> gradient{};
    int kernelCoordinate=-1;
    std::string branchIdentity;
    // Coefficient (i,j) is stored at 5*i+j, i,j in [0,4].
    std::array<CGAL::Gmpq,25> spherePolynomial(const CGAL::Gmpq& radius) const;
};

class PreparedRankTwoSource {
public:
    static PreparedRankTwoSource fromExact(
        const RankTwoExactCoefficients&,const std::string& frameIdentity);
    static PreparedRankTwoSource fromRepresented(
        const QuadricCoefficients10&,const std::string& frameIdentity);
    const RankTwoExactCoefficients& exactCoefficients() const;
    const QuadricCoefficients10& normalizedView() const;
    const std::array<CGAL::Gmpq,10>& normalizedViewErrorBounds() const;
    const RankTwoExactGraph& graph() const;
    const std::string& inputHash() const;
    const std::string& geometryHash() const;
    const std::string& numericViewHash() const;
private:
    struct Impl;
    explicit PreparedRankTwoSource(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};

} // namespace spin
