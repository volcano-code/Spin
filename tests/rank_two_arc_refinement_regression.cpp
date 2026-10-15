#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcRefinement.h")
#include "spin/detail/RankTwoArcRefinement.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include <set>
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=spin::detail::CertifiedAlgebraicReal;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool rejected=false;try{f();}catch(const spin::RankTwoSourceError& e){rejected=e.code()==code;}
    require(rejected,why);
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"arc-refinement-annulus-v1");
    const auto arrangement=analyzeRankTwoGlobalArrangement(source,Q(3,2));
    RankTwoArcSamplingPlan plan;std::string chosen;
    for(const auto& arc:arrangementRecords(arrangement).arcs){
        if(arc.ends[0].location!=RankTwoArcEndLocation::Finite
            ||arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        plan[arc.identity]={Q(1,2)};
        const E y=C::reload(arc.interiorWitness[1]).expression();
        bool F=false;for(const auto& owner:arc.owners)F|=owner.polynomial=="F";
        if(!arc.vertical && F && y>0 && y<1){require(chosen.empty(),"ambiguous frozen arc");chosen=arc.identity;}
    }
    require(!chosen.empty(),"frozen upper-inner quartic unavailable");
    const Q begin(31,64),end(33,64),error(1,1048576);
    plan.at(chosen)={begin,end};
    // This wide nonlinear interval cannot meet the locked error with one
    // chord. A no-op subdivision or reused old-registry proof must fail.
    const RankTwoArcRefinementSettings settings{};
    const auto result=refineRankTwoArcInterval(arrangement,plan,chosen,begin,end,error,settings);
    require(result.arcIdentity()==chosen && result.proofs().size()>2,"actual nonlinear interval was not subdivided");
    require(result.rounds().size()>1 && result.rounds().front().splitEdges>0
        && result.rounds().back().splitEdges==0,"failed chords were not actually refined");
    std::set<std::string> proofIds;Q next=begin;
    for(const auto& proof:result.proofs()){
        require(proof.registryIdentity()==result.registry().scientificHash()
            && proof.arrangementIdentity()==arrangement.scientificHash(),"old-registry proof survived a partition change");
        require(proof.strip().chordDistanceUpperBound<=error
            && proof.obstacles().size()+1==arrangement.factors().size(),"unproved chord or missing native factor");
        const RankTwoRegistryEdge* edge=nullptr;
        for(const auto& candidate:result.registry().edges())if(candidate.identity==proof.edgeIdentity())edge=&candidate;
        require(edge && proofIds.insert(edge->identity).second,"foreign or duplicated accepted edge");
        const auto& a=result.registry().points()[edge->points[0]];
        const auto& b=result.registry().points()[edge->points[1]];
        require(a.kind==RankTwoRegistryPointKind::ArcSample && b.kind==RankTwoRegistryPointKind::ArcSample
            && a.fraction==next && b.fraction>a.fraction,"interval coverage has gap/reorder or unowned points");
        next=b.fraction;
    }
    require(next==end,"subdivision omitted the end of its exact moving interval");
    const auto again=refineRankTwoArcInterval(arrangement,plan,chosen,begin,end,error,settings);
    require(result.scientificHash()==again.scientificHash()
        && result.canonicalRecord()==again.canonicalRecord(),"automatic exact boxes depend on previous root refinement");
    RankTwoSourceLimits exhausted;exhausted.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoArcInterval(arrangement,plan,chosen,begin,end,error,settings,exhausted);},
        "RESOURCE_LIMIT","unproved coarse chord escaped its subdivision budget");
    auto bad=plan;bad.erase(chosen);
    rejects([&]{refineRankTwoArcInterval(arrangement,bad,chosen,begin,end,error,settings);},
        "INPUT_IDENTITY_MISMATCH","missing actual arc plan silently repaired");
    rejects([&]{refineRankTwoArcInterval(arrangement,plan,chosen,Q(0),end,error,settings);},
        "INPUT_IDENTITY_MISMATCH","critical endpoint was passed as an ordinary interior interval");
    rejects([&]{refineRankTwoArcInterval(arrangement,plan,"foreign",begin,end,error,settings);},
        "INPUT_IDENTITY_MISMATCH","foreign arc accepted for refinement");
    // Characterize the exact zero-curvature boundary: no subdivision depth
    // must be spent merely because the native line is vertical in (x,y).
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto overlay=analyzeRankTwoNormalFace(source,arrangement,cone);
    RankTwoArcSamplingPlan verticalPlan;std::string vertical;
    for(const auto& arc:arrangementRecords(overlay).arcs){
        if(arc.ends[0].location!=RankTwoArcEndLocation::Finite
            ||arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        verticalPlan[arc.identity]={Q(1,2)};
        bool owner=false;for(const auto& o:arc.owners)owner|=o.polynomial=="NORMAL_EDGE:0:2";
        if(arc.vertical && owner && C::reload(arc.interiorWitness[1]).expression()==E(1)){
            require(vertical.empty(),"ambiguous exact vertical interval");vertical=arc.identity;
        }
    }
    require(!vertical.empty(),"frozen exact vertical interval missing");
    verticalPlan.at(vertical)={begin,end};
    const auto line=refineRankTwoArcInterval(overlay,verticalPlan,vertical,begin,end,error,settings,exhausted);
    require(line.rounds().size()==1 && line.proofs().size()==1
        && line.rounds()[0].splitEdges==0 && line.proofs()[0].projectionAxis()==1
        && line.proofs()[0].strip().chordDistanceUpperBound==0,
        "exact vertical normal edge consumed subdivision or acquired spurious error");
    std::ofstream out(path);require(bool(out),"cannot open refinement evidence");
    out<<"rank-two-regular-interval-test-v1\nREFINEMENT "<<result.scientificHash()<<'\n'
        <<result.canonicalRecord()<<"REFINEMENT "<<line.scientificHash()<<'\n'<<line.canonicalRecord();
    out.flush();require(bool(out),"incomplete refinement evidence");
    std::cout<<"PASS automatic exact boxes, actual dyadic subdivision, final shared-registry proofs and resource rejection\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcRefinement.h")
    if(argc!=2)throw std::invalid_argument("refinement evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: automatic shared regular-arc subdivision with exact strip proofs");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
