#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoRegularRegistry.h")
#include "spin/detail/RankTwoRegularRegistry.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include <cstring>
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& a){
    spin::RankTwoArcSamplingPlan p;
    for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite
            && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
            p.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    return p;
}
void check(const spin::RankTwoArrangement& a,const spin::RankTwoRegularRegistryRefinement& result,
    std::ostream& out,bool annulus){
    const auto& registry=result.registry();
    require(result.proofs().size()==registry.edges().size() && !registry.edges().empty(),
        "partial or duplicate edge coverage granted whole registry refinement");
    std::set<std::string> covered;
    for(const auto& proof:result.proofs()){
        require(proof.arrangementIdentity()==a.scientificHash()
            && proof.registryIdentity()==registry.scientificHash()
            && covered.emplace(proof.edgeIdentity()).second,
            "accepted a stale proof across a changed partition");
        require(proof.strip().chordDistanceUpperBound<=Q(1,32),"whole registry error relaxed");
        const spin::RankTwoRegistryEdge* edge=nullptr;
        for(const auto& e:registry.edges())if(e.identity==proof.edgeIdentity())edge=&e;
        require(edge && proof.pointIdentities()[0]==registry.points()[edge->points[0]].identity
            && proof.pointIdentities()[1]==registry.points()[edge->points[1]].identity,
            "proof does not consume actual exported endpoints");
    }
    require(result.rounds().size()>1 && result.rounds().front().splitEdges>0
        && result.rounds().back().splitEdges==0,
        "frozen coarse registry did not exercise actual subdivision");
    for(const auto& round:result.rounds())require(round.selectedEdges==round.certifiedEdges+round.splitEdges,
        "unaccounted edge silently dropped during refinement");
    std::size_t holes=0,unresolved=0;
    for(const auto& r:registry.regions()){
        holes+=r.innerLoops.size();for(const auto& l:r.outerLoops)unresolved+=!l.finite;
    }
    require(holes>0 && unresolved==1,"native holes or fictitious infinity loop disappeared");
    if(annulus)require(holes>1,"annulus hole topology lost");
    out<<"REGULAR_REGISTRY "<<result.scientificHash()<<'\n'<<result.canonicalRecord();
}
void tests(const char* path){
    using namespace spin;
    std::size_t stripHits=0;
    detail::ScopedRankTwoDiagnostics actualMemo([](const char* phase,void* p)noexcept{
        if(std::strcmp(phase,"regular-strip:memo-hit")==0)++*static_cast<std::size_t*>(p);
    },&stripHits);
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"whole-regular-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));
    const auto p=plan(a);
    const auto result=refineRankTwoRegularRegistry(a,p,Q(1,32));
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"whole-regular-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    const auto other=refineRankTwoRegularRegistry(b,plan(b),Q(1,32));
    require(stripHits>0,"MISSING_FEATURE: actual refiner recomputes identical strip queries");
    auto missing=p;missing.erase(missing.begin());
    rejects([&]{refineRankTwoRegularRegistry(a,missing,Q(1,32));},"INPUT_IDENTITY_MISMATCH","missing native arc accepted");
    RankTwoSourceLimits zero;zero.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoRegularRegistry(a,p,Q(1,32),{},zero);},"RESOURCE_LIMIT","coarse unproved edges granted completion");
    RankTwoSourceLimits small;small.maximumVertices=1;
    rejects([&]{refineRankTwoRegularRegistry(a,p,Q(1,32),{},small);},"RESOURCE_LIMIT","vertex cap ignored");
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto infinite=analyzeRankTwoNormalFace(disk,a,cone);
    rejects([&]{refineRankTwoRegularRegistry(infinite,plan(infinite),Q(1,32));},
        "ARRANGEMENT_UNRESOLVED","infinite native geometry silently truncated");
    std::ofstream out(path);require(bool(out),"cannot save whole registry evidence");
    out<<"rank-two-whole-regular-registry-test-v1\n";
    check(a,result,out,false);check(b,other,out,true);
    out.flush();require(bool(out),"incomplete whole registry evidence");
    std::cout<<"PASS actual all-edge finite regular refinement, original holes, stale partition rebuild and resource rejection\n";
    std::cout<<"STRIP_MEMO_CONSUMED "<<stripHits<<" FRESH_REGISTRY_PROOFS\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoRegularRegistry.h")
    if(argc!=2)throw std::invalid_argument("whole registry evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: automatic refinement proves every finite regular registry edge");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
