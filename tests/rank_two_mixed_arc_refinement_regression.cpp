#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoMixedArcRefinement.h")
#include "spin/detail/RankTwoMixedArcRefinement.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include <set>
namespace {
using Q=CGAL::Gmpq;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"endpoint-neighborhood-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(source,global,cone);
    std::string event;RankTwoArcSamplingPlan plan;
    for(const auto& v:arrangementRecords(a).vertices)
        if(C::reload(v.coordinates[0]).expression()==1 && C::reload(v.coordinates[1]).expression()==0)event=v.identity;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    // Registered small-scale method test, NOT the paper's geometric accuracy.
    // Both local proof policies are frozen before any refinement result exists.
    RankTwoEndpointNeighborhoodSettings ns;ns.maximumChordDistance=Q(1,32);ns.heightMargin=Q(1,8);
    const auto n=certifyRankTwoEndpointNeighborhood(a,event,ns);
    const auto initial=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n});
    std::string arcId,begin;
    for(const auto& edge:initial.edges()){
        bool h=false;for(const auto& owner:edge.owners)h|=owner.polynomial=="h";
        if(!h || !edge.endpointNeighborhoodIdentity.empty())continue;
        const auto& l=initial.points()[edge.points[0]],&r=initial.points()[edge.points[1]];
        if(l.coordinates[1].expression()>0 && r.coordinates[1].expression()>0
            && r.kind==RankTwoRegistryPointKind::SharedGate && l.kind==RankTwoRegistryPointKind::ArcSample){
            require(arcId.empty(),"ambiguous frozen h span");arcId=edge.arcIdentity;begin=l.identity;
        }
    }
    require(!arcId.empty(),"frozen h span absent");
    const auto result=refineRankTwoMixedArcInterval(a,plan,{n},arcId,{begin,event},Q(1,32));
    require(result.arcIdentity()==arcId && result.endpointIdentities()==std::array<std::string,2>{begin,event},
        "requested actual interval identity changed");
    const auto& registry=result.registry();
    require(result.rounds().size()>1 && result.rounds().front().splitEdges>0
        && result.rounds().back().splitEdges==0 && !registry.subdivisions().empty(),
        "coarse regular interval was not actually refined beside the endpoint neighborhood");
    std::set<std::string> covered;
    for(const auto& id:result.endpointEdges())require(covered.emplace(id).second,"duplicate endpoint edge");
    require(result.endpointEdges().size()==1 && !result.proofs().empty(),"fan-only or strip-only coverage");
    for(const auto& p:result.proofs()){
        require(p.registryIdentity()==registry.scientificHash() && p.arrangementIdentity()==a.scientificHash()
            && covered.emplace(p.edgeIdentity()).second,"stale or duplicated final partition proof");
        require(p.strip().chordDistanceUpperBound<=Q(1,32),"regular error relaxed");
    }
    const auto& last=result.rounds().back();
    require(covered.size()==last.selectedEdges && last.selectedEdges==last.endpointEdges+last.regularEdges,
        "selected chain has an unproved gap");
    for(const auto& round:result.rounds())require(round.selectedEdges==round.endpointEdges+round.regularEdges+round.splitEdges,
        "a failed edge disappeared from round accounting");
    std::size_t bound=0;for(const auto& e:registry.edges())bound+=!e.endpointNeighborhoodIdentity.empty();
    require(bound==8,"refinement invalidated another incident endpoint proof");
    for(const auto& p:registry.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision){
        const auto x=p.coordinates[0].expression(),y=p.coordinates[1].expression();
        require((x*x+y*y+y-1).sign()==0 && y>0,"refinement published a chord or wrong native branch");
    }
    RankTwoSourceLimits zero;zero.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoMixedArcInterval(a,plan,{n},arcId,{begin,event},Q(1,32),{},zero);},
        "RESOURCE_LIMIT","unproved coarse mixed interval granted success");
    rejects([&]{refineRankTwoMixedArcInterval(a,plan,{n},arcId,{event,begin},Q(1,32));},
        "INPUT_IDENTITY_MISMATCH","reversed actual interval accepted");
    rejects([&]{refineRankTwoMixedArcInterval(a,plan,{n},arcId,{begin,"undeclared"},Q(1,32));},
        "INPUT_IDENTITY_MISMATCH","undeclared moving endpoint accepted");
    RankTwoSourceLimits tiny;tiny.maximumVertices=1;
    rejects([&]{refineRankTwoMixedArcInterval(a,plan,{n},arcId,{begin,event},Q(1,32),{},tiny);},
        "RESOURCE_LIMIT","mixed vertex budget ignored");
    std::ofstream out(path);require(bool(out),"cannot write mixed interval evidence");
    out<<"rank-two-mixed-arc-refinement-test-v1\nMIXED_INTERVAL "<<result.scientificHash()<<'\n'<<result.canonicalRecord();
    out.flush();require(bool(out),"incomplete mixed interval evidence");
    std::cout<<"PASS actual endpoint-plus-regular chain, gate-relative subdivision, final-partition proof rebuild and bounded failure; no global topology grant\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoMixedArcRefinement.h")
    if(argc!=2)throw std::invalid_argument("mixed interval evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual common endpoint-neighborhood and regular-strip interval coverage");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
