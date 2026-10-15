#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoFiniteArcCoverage.h")
#include "spin/detail/RankTwoFiniteArcCoverage.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <chrono>
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& a){
    spin::RankTwoArcSamplingPlan p;
    for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite
            && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
            p.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    return p;
}
void check(const spin::RankTwoArrangement& a,const spin::RankTwoFiniteArcCoverage& c,
    std::ostream& out,const char* name,bool holes){
    const auto& r=c.partition().registry();const auto& mixed=c.refinement();
    require(c.arrangementIdentity()==a.scientificHash()
        && r.scientificHash()==mixed.registry().scientificHash(),"coverage does not bind final common partition");
    std::set<std::string> wanted,covered;
    for(const auto i:c.partition().edgeIndices())wanted.insert(r.edges().at(i).identity);
    for(const auto& p:mixed.proofs()){
        require(p.registryIdentity()==r.scientificHash() && p.arrangementIdentity()==a.scientificHash()
            && p.strip().chordDistanceUpperBound<=Q(1,32),"stale or relaxed actual strip proof");
        require(covered.insert(p.edgeIdentity()).second,"duplicate strip counted as new coverage");
    }
    for(const auto& e:mixed.endpointEdges())require(covered.insert(e).second,"endpoint/strip coverage overlaps");
    require(covered==wanted && c.edgeCount()==wanted.size(),"some finite domain edges lack a proof");
    std::set<std::string> arcIds;
    for(const auto i:c.partition().edgeIndices())arcIds.insert(r.edges().at(i).arcIdentity);
    require(arcIds.size()==mixed.intervals().size(),"native arc ends or whole interval omitted");
    for(const auto& i:mixed.intervals()){
        const auto arc=std::find_if(r.arcs().begin(),r.arcs().end(),[&](const auto& x){return x.identity==i.arcIdentity;});
        require(arc!=r.arcs().end() && arc->finite
            && i.pointIdentities[0]==*arc->ends[0].vertexIdentity
            && i.pointIdentities[1]==*arc->ends[1].vertexIdentity,"coverage stops at an interior sample");
    }
    require(mixed.rounds().size()>1 && mixed.rounds().front().splitEdges>0
        && mixed.rounds().back().splitEdges==0,"actual finite domain subdivision was not exercised");
    std::size_t nativeHoles=0,empty=0;
    for(const auto i:c.partition().regionIndices()){
        const auto& region=r.regions().at(i);nativeHoles+=region.innerLoops.size();
        empty+=region.kind==spin::RankTwoRegionKind::Empty;
    }
    require(!holes || nativeHoles>0,"annulus holes lost in full coverage");
    require(c.canonicalRecord().find("global_isotopy NOT_GRANTED\n")!=std::string::npos
        && c.canonicalRecord().find("strict_thickness NOT_GRANTED\n")!=std::string::npos,
        "local coverage incorrectly promoted to full carrier");
    out<<"CASE "<<name<<" COVERAGE "<<c.scientificHash()<<'\n'<<c.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" all_finite_edges="<<wanted.size()<<" arcs="<<arcIds.size()
        <<" native_holes="<<nativeHoles<<" Empty="<<empty<<'\n';
}
void tests(const char* path){
    using namespace spin;std::ostringstream saved;saved<<"rank-two-finite-arc-coverage-test-v1\n";
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));const auto p=plan(a);
    std::cerr<<"CASE_BEGIN full-disk\n";
    const auto first=refineRankTwoFiniteArcCoverage(a,p,{},Q(1,32));check(a,first,saved,"full-disk",false);
    std::cerr<<"CASE_END full-disk\n";
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));
    std::cerr<<"CASE_BEGIN annulus\n";
    const auto second=refineRankTwoFiniteArcCoverage(b,plan(b),{},Q(1,32));check(b,second,saved,"annulus",true);
    std::cerr<<"CASE_END annulus\n";
    // A genuine successful local proof is NOT full-domain coverage. The
    // production omission this catches is checking only supplied intervals.
    const auto& arc=first.refinement().intervals().front();
    const auto partial=refineRankTwoMixedArcIntervals(a,p,{}, {arc},Q(1,32));
    const auto partialPartition=buildRankTwoFinitePartition(a,partial.registry());
    rejects([&]{certifyRankTwoFiniteArcCoverage(a,partialPartition,partial);},
        "ARRANGEMENT_UNRESOLVED","successful single native arc passed as complete finite domain");
    rejects([&]{certifyRankTwoFiniteArcCoverage(b,first.partition(),first.refinement());},
        "INPUT_IDENTITY_MISMATCH","foreign authority accepted");
    rejects([&]{certifyRankTwoFiniteArcCoverage(a,first.partition(),partial);},
        "INPUT_IDENTITY_MISMATCH","different final registry reused old proof");
    auto missing=p;missing.erase(missing.begin());
    rejects([&]{refineRankTwoFiniteArcCoverage(a,missing,{},Q(1,32));},
        "INPUT_IDENTITY_MISMATCH","missing sampling arc silently filled or omitted");
    RankTwoSourceLimits zero;zero.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoFiniteArcCoverage(a,p,{},Q(1,32),{},zero);},
        "RESOURCE_LIMIT","partial round granted positive coverage");
    // Outside infinite arcs remain under exact native exclusion, not a
    // permission to accept unproved finite endpoint segments.
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto clipped=analyzeRankTwoNormalFace(disk,a,cone);
    const auto cr=buildRankTwoArcRegistry(clipped,plan(clipped));
    const auto cp=buildRankTwoFinitePartition(clipped,cr);
    require(!cp.excludedArcs().empty(),"infinite-exclusion counterexample missing");
    rejects([&]{refineRankTwoFiniteArcCoverage(clipped,plan(clipped),{},Q(1,32),{},zero);},
        "RESOURCE_LIMIT","unproved clipped endpoint promoted to complete coverage");
    std::ofstream out(path);require(bool(out),"cannot save finite coverage evidence");out<<saved.str();out.flush();
    require(bool(out),"incomplete finite coverage evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoFiniteArcCoverage.h")
    if(argc!=2)throw std::invalid_argument("finite coverage evidence path required");
    const auto start=std::chrono::steady_clock::now();
    const auto observe=+[](const char* phase,void* context)noexcept{
        try{
            const auto& origin=*static_cast<const std::chrono::steady_clock::time_point*>(context);
            std::cerr<<"DIAGNOSTIC_NS "<<std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-origin).count()<<' '<<phase<<'\n';
        }catch(...){}
    };
    spin::detail::ScopedRankTwoDiagnostics observer(observe,const_cast<std::chrono::steady_clock::time_point*>(&start));
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: complete native finite-domain endpoint/strip coverage, not selected intervals");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
