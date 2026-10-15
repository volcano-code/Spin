#include "detail/RankTwoArcRegistrySeams.h"
#include <algorithm>

namespace spin {
RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const std::vector<RankTwoEndpointNeighborhood>&,const std::vector<RankTwoArcSubdivision>&,
    const RankTwoRegistrySeamImport*,const RankTwoSourceLimits&);
RankTwoRegistrySeamImport makeRankTwoRegistrySeamImport(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const RankTwoSharedFiniteSeams& master,
    std::size_t face){
    const auto& whole=wholeArrangementRecords(master.wholeArrangement());
    if(face>=whole.normalFaces.size()||whole.normalFaces[face].scientificHash()!=a.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","master normal-face authority differs from native registry input");
    RankTwoRegistrySeamImport imported;imported.authority_=std::make_shared<const RankTwoSharedFiniteSeams>(master);
    imported.masterIdentity_=master.scientificHash();imported.faceAuthority_=a.scientificHash();imported.face_=face;
    std::map<std::string,const RankTwoArcRecord*> arcs;
    for(const auto& arc:arrangementRecords(a).arcs)arcs.emplace(arc.identity,&arc);
    for(const auto& interval:master.intervals())if(interval.finite)for(int side=0;side<2;++side)
        if(interval.incidentFaces[side]==face){
            const auto found=arcs.find(interval.sourceArcs[side]);
            if(found==arcs.end())throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","master interval lacks actual native source arc");
            const auto& arc=*found->second;const auto indices=master.consumerPoints(interval.identity,side);
            if(indices.size()<2||!arc.ends[0].vertexIdentity||!arc.ends[1].vertexIdentity
                ||master.points().at(indices.front()).identity!=*arc.ends[0].vertexIdentity
                ||master.points().at(indices.back()).identity!=*arc.ends[1].vertexIdentity)
                throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","master interval is a partial native arc; certified composition required");
            std::vector<RankTwoRegistryPoint> chain;std::vector<CGAL::Gmpq> fractions;
            for(const auto k:indices){const auto& p=master.points().at(k);
                if(p.endpoint)chain.push_back({p.identity,p.identity,{},RankTwoRegistryPointKind::Event,p.coordinates,{},CGAL::Gmpq(0)});
                else {chain.push_back({p.identity,{},arc.identity,RankTwoRegistryPointKind::SharedSeam,p.coordinates,p.parentIdentities,p.fraction});
                    fractions.push_back(interval.sourceReversed[side]?1-p.fraction:p.fraction);}}
            const auto request=plan.find(arc.identity);auto requested=request==plan.end()?std::vector<CGAL::Gmpq>{}:request->second;
            std::sort(requested.begin(),requested.end());
            if(request==plan.end()||requested!=fractions)
                throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","native seam sampling plan differs from published master fractions");
            if(!imported.chains_.emplace(arc.identity,std::move(chain)).second)
                throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","multiple master subintervals need certified native-arc composition");
        }
    return imported;
}
RankTwoArcRegistry buildRankTwoArcRegistryWithSharedSeams(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const RankTwoSharedFiniteSeams& master,
    std::size_t face,const RankTwoSourceLimits& limits){
    const auto imported=makeRankTwoRegistrySeamImport(a,plan,master,face);
    return buildRankTwoArcRegistryImported(a,plan,{},{},{},&imported,limits);
}
} // namespace spin
