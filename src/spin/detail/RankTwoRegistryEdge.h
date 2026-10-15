#pragma once
#include "RankTwoArcRegistry.h"
#include "RankTwoArcStrip.h"
#include "RankTwoStripMemo.h"

namespace spin {
// A private, actual-registry chord proof. The full immutable arrangement
// supplies EVERY other zero-set factor, including normal-boundary owners.
// This certifies one regular graph segment, not all chords or a 3D mesh.
class RankTwoRegularEdge {
public:
    const std::string& arrangementIdentity()const;
    const std::string& registryIdentity()const;
    const std::string& edgeIdentity()const;
    const std::array<std::string,2>& pointIdentities()const;
    const std::array<std::array<detail::CertifiedAlgebraicReal,2>,2>& endpointCoordinates()const;
    int projectionAxis()const;
    const detail::RankTwoStripBox& box()const;
    const CGAL::Gmpq& maximumChordDistance()const;
    const detail::RankTwoRegularStripCertificate& strip()const;
    const RankTwoZeroFactor& support()const;
    const std::vector<RankTwoZeroFactor>& obstacles()const;
    std::size_t checkedVertexCount()const;
    const std::string& belowRegion()const;
    const std::string& aboveRegion()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoRegularEdge(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const std::string&,int,
        const detail::RankTwoStripBox&,const CGAL::Gmpq&);
    friend RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const std::string&,int,
        const detail::RankTwoStripBox&,const CGAL::Gmpq&,detail::RankTwoRegularStripMemo*);
    friend RankTwoRegularEdge certifyRankTwoRegistryEdgeAtIndex(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const std::string&,std::size_t,int,
        const detail::RankTwoStripBox&,const CGAL::Gmpq&,detail::RankTwoRegularStripMemo*);
    friend std::variant<RankTwoRegularEdge,std::string> tryCertifyRankTwoRegistryEdgeAtIndex(
        const RankTwoArrangement&,const RankTwoArcRegistry&,const std::string&,std::size_t,int,
        const detail::RankTwoStripBox&,const CGAL::Gmpq&,detail::RankTwoRegularStripMemo*);
};
// axis=0 uses (x,y), axis=1 uses the isometric chart (y,x). Shared endpoints
// come from the immutable registry, never regenerated or supplied by caller.
// Whole-strip Py/bracketing/C2/chord bounds and all-factor exclusion are
// checked, then all native events in the chord span. An endpoint fan edge
// must retain its simultaneous fan proof instead of acquiring this one.
RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const std::string& edgeIdentity,int projectionAxis,
    const detail::RankTwoStripBox&,const CGAL::Gmpq& maximumChordDistance);
// Reuse ONLY the pure strip arithmetic. All current registry membership,
// endpoint, native event and authority bindings are checked anew.
RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const std::string& edgeIdentity,int projectionAxis,
    const detail::RankTwoStripBox&,const CGAL::Gmpq& maximumChordDistance,
    detail::RankTwoRegularStripMemo*);
// The refiner already owns an index in the CURRENT immutable publication.
// Verify its bounds and identity, then run the same complete native proof.
// An index is not proof authority and cannot be reused across publications.
RankTwoRegularEdge certifyRankTwoRegistryEdgeAtIndex(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const std::string& edgeIdentity,std::size_t edgeIndex,
    int projectionAxis,const detail::RankTwoStripBox&,
    const CGAL::Gmpq& maximumChordDistance,detail::RankTwoRegularStripMemo*);
// Same fresh native proof. Only mathematical inability returns its original
// witness; invalid authority, resources and exact-environment errors throw.
using RankTwoRegularEdgeAttempt=std::variant<RankTwoRegularEdge,std::string>;
RankTwoRegularEdgeAttempt tryCertifyRankTwoRegistryEdgeAtIndex(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const std::string&,std::size_t,int,
    const detail::RankTwoStripBox&,const CGAL::Gmpq&,detail::RankTwoRegularStripMemo*);
} // namespace spin
