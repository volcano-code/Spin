#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoChordSweep.h")
#include "spin/detail/RankTwoChordSweep.h"
#include "spin/detail/RankTwoRegularRegistry.h"
#include "spin/detail/RankTwoRegistryRotation.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include <chrono>
#include <string_view>
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& a,bool quarters,bool skew=false){
    spin::RankTwoArcSamplingPlan p;
    for(const auto& arc:spin::arrangementRecords(a).arcs){
        if(arc.ends[0].location!=spin::RankTwoArcEndLocation::Finite
            ||arc.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        bool domain=false;for(const auto& owner:arc.owners)domain|=owner.polynomial=="D";
        p.emplace(arc.identity,quarters?std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)}:
            std::vector<Q>{skew&&domain?Q(1,128):Q(1,2)});
    }
    return p;
}
void coverage(const spin::RankTwoChordIncidence& c,const spin::RankTwoArcRegistry& r){
    const auto n=r.points().size(),m=r.edges().size();
    require(c.pointPairCount()==n*(n-1)/2 && c.pointEdgeCount()==m*(n-2)
        && c.edgePairCount()==m*(m-1)/2,"sweep silently omitted incidence coverage");
    require(c.registryIdentity()==r.scientificHash(),"sweep did not bind actual final registry");
    require(c.predicateGroupCount()>0 && c.predicateGroupCount()<=1000000,
        "sweep exceeded or falsified fixed predicate budget");
}
void tests(const char* path){
    using namespace spin;
    const auto beginning=std::chrono::steady_clock::now();
    const auto observe=[](const char* phase,void* context) noexcept {
        const std::string_view name(phase);
        if(name.find("chord-sweep:")!=0 && name.find("test:")!=0
            && name.find("regular-registry:round:")!=0)return;
        const auto start=*static_cast<std::chrono::steady_clock::time_point*>(context);
        std::cerr<<"OBSERVATION "<<phase<<" ELAPSED_SECONDS "
            <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
    };
    auto started=beginning;detail::ScopedRankTwoDiagnostics diagnostics(observe,&started);
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"whole-regular-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));
    const auto small=buildRankTwoArcRegistry(a,plan(a,false));
    const auto legacy=certifyRankTwoChordIncidence(a,small);
    const auto quick=certifyRankTwoChordIncidenceSweep(a,small);
    coverage(quick,small);
    require(quick.pointPairCount()==legacy.pointPairCount() && quick.pointEdgeCount()==legacy.pointEdgeCount()
        && quick.edgePairCount()==legacy.edgePairCount(),"exact sweep disagrees with exhaustive small-graph coverage");
    require(certifyRankTwoChordIncidenceSweep(a,small).canonicalRecord()==quick.canonicalRecord(),
        "same-process sweep depends on predicate evaluation history");
    RankTwoChordSweepSettings exhausted;exhausted.maximumPredicateGroups=0;
    rejects([&]{certifyRankTwoChordIncidenceSweep(a,small,exhausted);},"RESOURCE_LIMIT",
        "budget exhaustion returned incomplete positive incidence");
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"whole-regular-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    for(bool skew:{false,true}){
        const auto crossed=buildRankTwoArcRegistry(b,plan(b,false,skew));
        rejects([&]{certifyRankTwoChordIncidenceSweep(b,crossed);},"ARRANGEMENT_UNRESOLVED",
            "sweep candidate pruning hid actual coarse crossing");
    }
    const auto r=refineRankTwoRegularRegistry(a,plan(a,true),Q(1,32));
    detail::rankTwoDiagnostic("test:refined-circle-ready");
    const auto s=refineRankTwoRegularRegistry(b,plan(b,true),Q(1,32));
    detail::rankTwoDiagnostic("test:refined-annulus-ready");
    require(r.registry().edges().size()==24 && s.registry().edges().size()==2180,
        "frozen regular-refinement graph changed");
    const auto first=certifyRankTwoChordIncidenceSweep(a,r.registry());
    detail::rankTwoDiagnostic("test:circle-incidence-ready");
    const auto second=certifyRankTwoChordIncidenceSweep(b,s.registry());
    detail::rankTwoDiagnostic("test:annulus-incidence-ready");
    coverage(first,r.registry());coverage(second,s.registry());
    const auto firstRotation=certifyRankTwoRegistryRotation(a,r.registry(),first);
    detail::rankTwoDiagnostic("test:circle-rotation-ready");
    const auto secondRotation=certifyRankTwoRegistryRotation(b,s.registry(),second);
    detail::rankTwoDiagnostic("test:annulus-rotation-ready");
    std::ofstream regular(std::string(path)+".regular");require(bool(regular),"cannot save same-registry error proofs");
    regular<<"rank-two-whole-regular-registry-test-v1\nREGULAR_REGISTRY "<<r.scientificHash()<<'\n'
        <<r.canonicalRecord()<<"REGULAR_REGISTRY "<<s.scientificHash()<<'\n'<<s.canonicalRecord();
    regular.flush();require(bool(regular),"incomplete same-registry error proofs");
    require(firstRotation.finiteDirectionCount()==48 && secondRotation.finiteDirectionCount()==4360,
        "native rotation did not consume every actual refined chord");
    rejects([&]{certifyRankTwoChordIncidenceSweep(a,s.registry());},"INPUT_IDENTITY_MISMATCH",
        "sweep accepted foreign arrangement");
    rejects([&]{certifyRankTwoRegistryRotation(a,r.registry(),quick);},"INPUT_IDENTITY_MISMATCH",
        "coarse incidence was promoted to a different refined registry");
    std::ofstream out(path);require(bool(out),"cannot save sweep evidence");
    out<<"rank-two-chord-sweep-test-v1\n";
    for(const auto* p:{&quick,&first,&second})out<<"SWEEP "<<p->scientificHash()<<'\n'<<p->canonicalRecord();
    out<<"ROTATION "<<firstRotation.scientificHash()<<'\n'<<firstRotation.canonicalRecord()
        <<"ROTATION "<<secondRotation.scientificHash()<<'\n'<<secondRotation.canonicalRecord();
    out.flush();require(bool(out),"incomplete saved incidence/rotation");
    std::cout<<"PASS complete exact sweep incidence and native rotation on actual 24/2180-edge refined registries; coarse crossings and stale/budget rejection\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoChordSweep.h")
    if(argc!=2)throw std::invalid_argument("sweep evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact complete incidence sweep on refined actual shared chords");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
