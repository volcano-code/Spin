#pragma once
#include "RankTwoArcRegistrySubdivisions.h"
#include "RankTwoSharedFiniteSeams.h"
#include <set>

namespace spin {
class RankTwoSharedEndpointSeams;
// Only the immutable-master adapter may construct an import. The existing
// registry builder consumes its actual coordinates; it does not independently
// resample a seam, weld coordinates or grant a curve/isotopy certificate.
class RankTwoRegistrySeamImport {
    std::shared_ptr<const RankTwoSharedFiniteSeams> authority_;
    std::string masterIdentity_,faceAuthority_;
    std::size_t face_=0;
    std::map<std::string,std::vector<RankTwoRegistryPoint>> chains_;
    std::shared_ptr<const RankTwoSharedEndpointSeams> endpointAuthority_;
    std::map<std::string,std::string> endpointBindings_; // actual owner proof -> shared point
    std::map<std::string,std::size_t> endpointCounts_;
    std::set<std::string> endpointAliasProofs_; // same actual affine master dependency, not a coordinate weld
    RankTwoRegistrySeamImport()=default;
    friend RankTwoArcRegistry buildRankTwoArcRegistryWithSharedSeams(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
        const RankTwoSharedFiniteSeams&,std::size_t,const RankTwoSourceLimits&);
    friend RankTwoRegistrySeamImport makeRankTwoRegistrySeamImport(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
        const RankTwoSharedFiniteSeams&,std::size_t);
    friend RankTwoArcRegistry buildRankTwoArcRegistryWithSharedEndpointSeams(
        const RankTwoArrangement&,const RankTwoArcSamplingPlan&,
        const RankTwoSharedEndpointSeams&,std::size_t,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
        const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&);
    friend RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
        const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
        const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
        const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&,const RankTwoArcRegistry*);
};
RankTwoRegistrySeamImport makeRankTwoRegistrySeamImport(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const RankTwoSharedFiniteSeams&,std::size_t);
RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
    const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&);
// Continuation only: no caller-supplied coordinates or global cache. Reuse
// a factory-owned point only after checking the complete frozen authority.
RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
    const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&,const RankTwoArcRegistry*);
// The plan still names every finite native arc. For imported arcs its fractions
// must equal the master's fractions in native increasing order. Partial-arc
// common refinements require certified composition and are explicitly UNKNOWN.
RankTwoArcRegistry buildRankTwoArcRegistryWithSharedSeams(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const RankTwoSharedFiniteSeams&,std::size_t face,
    const RankTwoSourceLimits& limits={});
} // namespace spin
