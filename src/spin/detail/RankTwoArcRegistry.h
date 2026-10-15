#pragma once
#include "RankTwoArcSampling.h"
#include <map>
#include <vector>

namespace spin {
class RankTwoEndpointFan;
class RankTwoEndpointNeighborhood;
class RankTwoRegistrySeamImport;
// A plan must name EVERY finite immutable arc. Fractions are exact moving
// dependencies, not supplied coordinates. Order is normalized; duplicates,
// endpoints, unknown arcs and unbounded arcs are rejected.
using RankTwoArcSamplingPlan=std::map<std::string,std::vector<CGAL::Gmpq>>;
struct RankTwoArcSubdivision {
    std::string arcIdentity;
    std::array<std::string,2> pointIdentities;
    CGAL::Gmpq fraction;
};
enum class RankTwoRegistryPointKind { Event,ArcSample,SharedGate,BoundaryRay,ArcSubdivision,SharedSeam };
struct RankTwoEndpointPointBinding {
    std::string neighborhoodIdentity,proofIdentity,gateIdentity,fanIdentity;
};
struct RankTwoRegistryPoint {
    std::string identity,eventIdentity,arcIdentity;
    RankTwoRegistryPointKind kind;
    std::array<detail::CertifiedAlgebraicReal,2> coordinates;
    std::array<std::string,2> endpointIdentities;
    CGAL::Gmpq fraction;
    std::string gateIdentity,endpointFanIdentity;
    int projectionAxis=-1;
    std::string endpointNeighborhoodIdentity,endpointArcIdentity;
    // Additional sealed proofs for the SAME original moving dependency, e.g.
    // the F and D endpoint chords meeting at one original master midpoint.
    // The primary legacy fields stay intact; these are not extra vertices.
    std::vector<RankTwoEndpointPointBinding> additionalEndpointBindings;
};
struct RankTwoRegistryEdge {
    std::string identity,arcIdentity;
    std::array<std::size_t,2> points; // original native increasing direction
    std::vector<RankTwoPolynomialOwner> owners;
    std::string endpointFanIdentity;
    std::string endpointNeighborhoodIdentity;
};
struct RankTwoRegistryArc {
    std::string identity,supportPolynomial,leftRegion,rightRegion;
    bool finite=false,vertical=false;
    int branch=-1;
    std::array<RankTwoArcEndRecord,2> ends;
    std::vector<RankTwoPolynomialOwner> owners;
    std::vector<std::size_t> edges;
};
struct RankTwoRegistryDirectedEdge {std::size_t edgeIndex;bool reversed;};
struct RankTwoRegistryLoop {
    bool finite=false;
    std::vector<RankTwoDirectedArcRecord> sourceArcs;
    std::vector<RankTwoRegistryDirectedEdge> edges;
};
struct RankTwoRegistryRegion {
    std::string identity;
    RankTwoRegionKind kind;
    bool unbounded;
    std::array<int,3> originalSigns;
    std::vector<int> normalSigns;
    std::vector<RankTwoRegistryLoop> outerLoops,innerLoops;
    std::vector<std::size_t> isolatedPoints;
};
// Factory-certified exclusion of an OPEN unbounded incidence from the finite
// D domain. The actual endpoint, native arc/regions and complete star remain.
// A cone-filtered Outside label is NOT sufficient: both ORIGINAL D signs and
// the certified ray gate must be strictly negative, and D must not own the arc.
struct RankTwoOutsideIncidence {
    std::string arcIdentity,endpointIdentity,neighborhoodIdentity;
    std::string endpointArcIdentity,gateIdentity,leftRegion,rightRegion;
};
class RankTwoArcRegistry {
public:
    const std::string& arrangementIdentity()const;
    const std::vector<RankTwoRegistryPoint>& points()const;
    const std::vector<RankTwoRegistryEdge>& edges()const;
    const std::vector<RankTwoRegistryArc>& arcs()const;
    const std::vector<RankTwoRegistryRegion>& regions()const;
    const std::vector<RankTwoEndpointFan>& endpointFans()const;
    const std::vector<RankTwoEndpointNeighborhood>& endpointNeighborhoods()const;
    const std::vector<RankTwoArcSubdivision>& subdivisions()const;
    const std::vector<RankTwoOutsideIncidence>& outsideIncidences()const;
    // Returns shared global point indices, in the ORIGINAL directed DCEL
    // loop order. No closing duplicate. Infinite loops cannot be truncated.
    std::vector<std::size_t> loopVertices(const std::string& regionIdentity,
        bool inner,std::size_t loopIndex)const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoArcRegistry(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoArcRegistry rebuildRankTwoArcRegistryWithSubdivisions(
        const RankTwoArrangement&,const RankTwoArcRegistry&,const RankTwoArcSamplingPlan&,
        const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
        const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,
        const std::vector<RankTwoArcSubdivision>&,const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,
        const std::vector<RankTwoArcSubdivision>&,const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&,
        const RankTwoArcRegistry*);
    friend RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,
        const std::vector<RankTwoArcSubdivision>&,const RankTwoSourceLimits&);
};
// Factory-owned shared samples/edge references and original labels/holes.
// Infinite arcs are retained, not sampled or silently discarded. Their loops
// remain explicitly unresolved until exact finite clipping is implemented.
// This is NOT a chord topology/error proof and NOT triangulation permission:
// endpoint fans and regular strips must cover every edge before consuming
// these edges as non-crossing constraints. Full carrier acceptance is separate.
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
