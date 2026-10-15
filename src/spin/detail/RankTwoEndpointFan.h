#pragma once
#include "RankTwoEndpointArc.h"
#include <memory>
#include <vector>

namespace spin {
enum class RankTwoFanSupportKind { Graph,AtStartBoundary,Excluded };
struct RankTwoFanSupportRecord {
    std::string polynomial;
    std::vector<RankTwoPolynomialOwner> owners;
    RankTwoFanSupportKind kind;
    // Only Excluded uses this sign. Graph/boundary do not claim a sign.
    int exclusionSign=0;
    std::string arcIdentity;
};
class RankTwoEndpointFan {
public:
    const RankTwoMovingGate& gate() const;
    const std::vector<RankTwoSharedGateArcCertificate>& orderedArcs() const;
    const std::vector<RankTwoFanSupportRecord>& supportRecords() const;
    const std::string& canonicalRecord() const;
    const std::string& scientificHash() const;
private:
    struct Impl;
    explicit RankTwoEndpointFan(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoEndpointFan certifyRankTwoEndpointFan(const RankTwoArrangement&,
        const RankTwoMovingGate&,const detail::RankTwoStripBox&,const CGAL::Gmpq&,int);
};
// Certifies ONE half-sector from an actual event to ONE shared moving gate.
// All immutable zero-set factors/normal owners must either have their unique
// graph in this sector, be strictly excluded from the whole rectangle, or
// form the start projection boundary. Gate order plus exact event exclusion
// proves a simultaneous order-preserving deformation to straight chords.
// Region adjacency is checked against native DCEL IDs, not parity.
// This is not a complete finite carrier, hole triangulation or 3D bound.
RankTwoEndpointFan certifyRankTwoEndpointFan(const RankTwoArrangement&,
    const RankTwoMovingGate&,const detail::RankTwoStripBox&,
    const CGAL::Gmpq& maximumChordDistance,int maximumIsolationBits=512);
} // namespace spin
