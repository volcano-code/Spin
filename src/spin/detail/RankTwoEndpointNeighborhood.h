#pragma once
#include "RankTwoEndpointFan.h"

namespace spin {
struct RankTwoEndpointNeighborhoodSettings {
    CGAL::Gmpq maximumChordDistance=CGAL::Gmpq(1,1048576);
    CGAL::Gmpq heightMargin=CGAL::Gmpq(1,100);
    int maximumGateHalvings=16;
    std::size_t maximumCandidateTests=10000;
    int maximumIsolationBits=512;
};
struct RankTwoEndpointBoundaryRay {
    RankTwoMovingGate gate;
    RankTwoSharedGateArcCertificate proof;
    detail::RankTwoStripBox box;
};
// Factory-only proof of a COMPLETE local native endpoint star. Two opposite
// ordered graph sectors are joined with unchanged straight projection-boundary
// rays. Not a global embedding/isotopy, triangulation or material certificate.
class RankTwoEndpointNeighborhood {
public:
    const std::string& arrangementIdentity() const;
    const std::string& endpointIdentity() const;
    int projectionAxis() const;
    const std::vector<RankTwoEndpointFan>& fans() const;
    const std::vector<RankTwoEndpointBoundaryRay>& boundaryRays() const;
    const std::vector<std::string>& incidentArcs() const;
    const std::string& canonicalRecord() const;
    const std::string& scientificHash() const;
private:
    struct Impl;
    explicit RankTwoEndpointNeighborhood(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoEndpointNeighborhood certifyRankTwoEndpointNeighborhood(
        const RankTwoArrangement&,const std::string&,const RankTwoEndpointNeighborhoodSettings&);
};
RankTwoEndpointNeighborhood certifyRankTwoEndpointNeighborhood(
    const RankTwoArrangement&,const std::string& endpointIdentity,
    const RankTwoEndpointNeighborhoodSettings& settings={});
} // namespace spin
