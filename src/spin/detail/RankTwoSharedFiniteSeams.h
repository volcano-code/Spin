#pragma once
#include "RankTwoWholeArrangement.h"
#include <cstdint>

namespace spin {
struct RankTwoSharedFiniteSeamSettings {
    std::vector<CGAL::Gmpq> fractions{CGAL::Gmpq(1,2)};
    std::size_t maximumPoints=1000000;
    std::uint64_t maximumChecks=2000000;
    int maximumIsolationBits=512;
};
struct RankTwoSharedFinitePoint {
    std::string identity;
    bool endpoint;
    std::array<detail::CertifiedAlgebraicReal,2> coordinates;
    std::array<std::string,2> parentIdentities;
    CGAL::Gmpq fraction;
};
struct RankTwoSharedFiniteInterval {
    std::string identity;
    std::array<std::size_t,2> edge,incidentFaces;
    std::array<std::string,2> sourceArcs;
    std::array<bool,2> sourceReversed;
    std::array<int,3> originalSigns;
    std::vector<RankTwoPolynomialOwner> originalOwners;
    bool finite;
    std::vector<std::size_t> points; // canonical increasing Gauss t, SAME indices for both consumers
};
// Owns one exact coordinate pair per shared feature, not two rounded face
// copies. Every original edge/interval remains in the authority, including
// Empty and outside-projection parts. This table is not triangulation,
// global isotopy, profile precision, thickness or public carrier permission.
class RankTwoSharedFiniteSeams {
public:
    const RankTwoArrangement& wholeArrangement()const;
    const std::vector<RankTwoSharedFinitePoint>& points()const;
    const std::vector<RankTwoSharedFiniteInterval>& intervals()const;
    std::vector<std::size_t> consumerPoints(const std::string& interval,int side)const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoSharedFiniteSeams(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoSharedFiniteSeams buildRankTwoSharedFiniteSeams(
        const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoSharedFiniteSeamSettings&);
};
RankTwoSharedFiniteSeams buildRankTwoSharedFiniteSeams(const PreparedRankTwoSource&,
    const RankTwoArrangement&,const RankTwoSharedFiniteSeamSettings& settings={});
} // namespace spin
