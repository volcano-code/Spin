#pragma once
#include <cstdint>
#include "RankTwoArcSampling.h"

namespace spin {
struct RankTwoInteriorPointRequest {
    std::string regionIdentity;
    std::array<std::string,2> eventIdentities; // adjacent DISTINCT native x events
    std::array<std::string,2> boundaryArcIdentities; // lower/upper native fibers
    CGAL::Gmpq abscissaFraction=CGAL::Gmpq(1,2),ordinateFraction=CGAL::Gmpq(1,2);
};
struct RankTwoInteriorPointSettings {
    int maximumIsolationBits=512;
    std::uint64_t maximumNativeChecks=1000000,maximumFiberRoots=1000000;
};
// A new moving feature, NOT an ArcSample or a fixed rational region witness.
// Private factory validates the full native vertical gap and active material
// signs. Dependencies are retained even when the evaluated point is rational.
// Carrier insertion, triangle/edge thickness and source motion are separate.
class RankTwoMovingInteriorPoint {
public:
    const std::string& arrangementIdentity()const;
    const RankTwoInteriorPointRequest& definition()const;
    const std::array<detail::CertifiedAlgebraicReal,2>& coordinates()const;
    const std::array<detail::CertifiedAlgebraicReal,2>& fiberOrdinates()const;
    RankTwoRegionKind regionKind()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoMovingInteriorPoint(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoMovingInteriorPoint buildRankTwoMovingInteriorPoint(
        const RankTwoArrangement&,const RankTwoInteriorPointRequest&,const RankTwoInteriorPointSettings&);
};
RankTwoMovingInteriorPoint buildRankTwoMovingInteriorPoint(
    const RankTwoArrangement&,const RankTwoInteriorPointRequest&,const RankTwoInteriorPointSettings& settings={});
RankTwoMovingInteriorPoint reloadRankTwoMovingInteriorPoint(
    const RankTwoArrangement&,const RankTwoInteriorPointRequest&,const std::string&,
    const RankTwoInteriorPointSettings& settings={});
} // namespace spin
