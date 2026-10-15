#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcRegistryNeighborhoods.h")
#include "spin/detail/RankTwoArcRegistryNeighborhoods.h"
#include "spin/detail/RankTwoNormalOverlay.h"
namespace {
using Q=CGAL::Gmpq;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
void tests(const char* path){
    using namespace spin;
    const auto p=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"endpoint-neighborhood-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(p,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(p,global,cone);
    std::string event;RankTwoArcSamplingPlan plan;
    for(const auto& v:arrangementRecords(a).vertices)
        if(C::reload(v.coordinates[0]).expression()==1 && C::reload(v.coordinates[1]).expression()==0)event=v.identity;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    require(!event.empty(),"actual endpoint absent");
    const auto n=certifyRankTwoEndpointNeighborhood(a,event);
    const auto legacy=buildRankTwoArcRegistry(a,plan);
    const auto r=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n});
    require(r.endpointNeighborhoods().size()==1 && r.endpointFans().size()==2,
        "complete proof was not retained by shared registry");
    require(r.points().size()==legacy.points().size()+8 && r.edges().size()==legacy.edges().size()+8,
        "some of eight endpoint gates were lost or double-counted");
    require(r.canonicalRecord().find("rank-two-shared-arc-registry-v2\n")==0
        && legacy.canonicalRecord().find("rank-two-shared-arc-registry-v1\n")==0,
        "new evidence changed legacy format or lacks explicit version");
    std::map<std::string,const RankTwoRegistryPoint*> gates;
    std::map<std::string,const RankTwoSharedGateArcCertificate*> proofs;
    for(const auto& f:n.fans())for(const auto& proof:f.orderedArcs())proofs.emplace(proof.arcIdentity,&proof);
    for(const auto& ray:n.boundaryRays())proofs.emplace(ray.proof.arcIdentity,&ray.proof);
    std::size_t rays=0,bound=0;
    for(const auto& point:r.points())if(!point.endpointNeighborhoodIdentity.empty()){
        require(point.endpointNeighborhoodIdentity==n.scientificHash()
            && gates.emplace(point.arcIdentity,&point).second,"gate authority or unique incidence lost");
        const auto& proof=*proofs.at(point.arcIdentity);
        require(point.endpointArcIdentity==proof.identity && point.gateIdentity==proof.gateIdentity
            && point.coordinates[0].expression()==proof.gateCoordinates[0].expression()
            && point.coordinates[1].expression()==proof.gateCoordinates[1].expression(),
            "registry/checker/source coordinates or proof identities disagree");
        if(point.kind==RankTwoRegistryPointKind::BoundaryRay){++rays;
            require(point.endpointFanIdentity.empty() && proof.strip.chordDistanceUpperBound==0,
                "straight ray was relabeled as graph-sector evidence");}
    }
    require(gates.size()==8 && rays==2,"full native star not published");
    for(const auto& edge:r.edges())if(!edge.endpointNeighborhoodIdentity.empty()){
        ++bound;require(edge.endpointNeighborhoodIdentity==n.scientificHash(),"edge lost complete proof identity");
        const auto& x=r.points().at(edge.points[0]);const auto& y=r.points().at(edge.points[1]);
        require((x.eventIdentity==event && y.identity==gates.at(edge.arcIdentity)->identity)
            || (y.eventIdentity==event && x.identity==gates.at(edge.arcIdentity)->identity),
            "certificate attached to a different shared edge");
    }
    require(bound==8,"complete certificate did not bind all native endpoint edges");
    // Every finite directed ring consumes the SAME shared chain, in native order.
    for(const auto& region:r.regions())for(bool inner:{false,true}){
        const auto& loops=inner?region.innerLoops:region.outerLoops;
        for(std::size_t i=0;i<loops.size();++i){
            if(!loops[i].finite){rejects([&]{r.loopVertices(region.identity,inner,i);},
                "ARRANGEMENT_UNRESOLVED","infinity silently clipped");continue;}
            std::vector<RankTwoRegistryDirectedEdge> expected;
            for(const auto& ref:loops[i].sourceArcs){
                const RankTwoRegistryArc* parent=nullptr;for(const auto& arc:r.arcs())if(arc.identity==ref.identity)parent=&arc;
                require(parent && parent->finite,"finite ring lost native arc");
                if(ref.reversed)for(auto e=parent->edges.rbegin();e!=parent->edges.rend();++e)expected.push_back({*e,true});
                else for(auto e:parent->edges)expected.push_back({e,false});
            }
            require(expected.size()==loops[i].edges.size(),"consumer dropped shared gate");
            for(std::size_t j=0;j<expected.size();++j)
                require(expected[j].edgeIndex==loops[i].edges[j].edgeIndex
                    && expected[j].reversed==loops[i].edges[j].reversed,"consumer resampled or reversed native arc");
            require(r.loopVertices(region.identity,inner,i).size()==expected.size(),"finite ring not closed");
        }
    }
    rejects([&]{buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n,n});},"INPUT_IDENTITY_MISMATCH",
        "duplicate complete endpoint accepted");
    rejects([&]{buildRankTwoArcRegistryWithNeighborhoods(global,plan,{n});},"INPUT_IDENTITY_MISMATCH",
        "foreign neighborhood accepted");
    auto split=plan;const auto id=n.incidentArcs().front();
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.identity==id)
        split.at(id).push_back(arc.ends[0].vertexIdentity==event?Q(1,1099511627776L):Q(1099511627775L,1099511627776L));
    rejects([&]{buildRankTwoArcRegistryWithNeighborhoods(a,split,{n});},"ARRANGEMENT_UNRESOLVED",
        "extra sample split proof-bound endpoint chord");
    RankTwoSourceLimits tiny;tiny.maximumVertices=r.points().size()-1;
    rejects([&]{buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n},tiny);},"RESOURCE_LIMIT",
        "complete proof gates ignored vertex budget");
    const auto again=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n});
    require(again.scientificHash()==r.scientificHash() && again.canonicalRecord()==r.canonicalRecord(),
        "shared registry depends on prior root refinement");
    const auto empty=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{});
    require(empty.canonicalRecord()==legacy.canonicalRecord() && empty.scientificHash()==legacy.scientificHash(),
        "empty neighborhood path changed legacy scientific output");
    std::ofstream out(path);require(bool(out),"cannot save binding evidence");
    out<<"rank-two-neighborhood-registry-test-v1\nREGISTRY "<<r.scientificHash()<<'\n'<<r.canonicalRecord();
    out.flush();require(bool(out),"incomplete binding evidence");
    std::cout<<"PASS complete eight-direction neighborhood bound to actual shared edges and original directed consumers; no clipping or global topology grant\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcRegistryNeighborhoods.h")
    if(argc!=2)throw std::invalid_argument("binding evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: complete endpoint neighborhood bound to shared arc registry");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
