#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "spin/detail/RankTwoRegistryContinuation.h"
#include "spin/detail/RankTwoDiagnostics.h"

namespace {
using Q=CGAL::Gmpq;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){bool caught=false;try{f();}
    catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"reuse bypassed immutable authority or resource validation");}
struct Trace {std::string retained;std::size_t reevaluated=0,newEvaluated=0;};
void observe(const char* p,void* raw)noexcept{
    auto& t=*static_cast<Trace*>(raw);const std::string phase(p),suffix=":coordinates-ready";
    if(phase.rfind("registry:subdivision:",0)!=0||phase.size()<suffix.size()
        ||phase.compare(phase.size()-suffix.size(),suffix.size(),suffix)!=0)return;
    if(phase=="registry:subdivision:"+t.retained+suffix)++t.reevaluated;
    else ++t.newEvaluated;
}
void test(const char* path){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"registry-reuse-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(1));RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto base=buildRankTwoArcRegistry(a,plan);const RankTwoRegistryEdge* edge=nullptr;
    for(const auto& e:base.edges())if(base.points().at(e.points[0]).kind==RankTwoRegistryPointKind::ArcSample
        &&base.points().at(e.points[1]).kind==RankTwoRegistryPointKind::ArcSample){edge=&e;break;}
    require(edge,"actual curved primitive interval absent");
    const RankTwoArcSubdivision first{edge->arcIdentity,{base.points().at(edge->points[0]).identity,
        base.points().at(edge->points[1]).identity},Q(1,3)};
    const auto initial=rebuildRankTwoArcRegistryWithSubdivisions(a,base,plan,{}, {first});
    const auto old=std::find_if(initial.points().begin(),initial.points().end(),
        [](const auto& p){return p.kind==RankTwoRegistryPointKind::ArcSubdivision;});
    require(old!=initial.points().end(),"first actual moving dependency absent");
    const RankTwoArcSubdivision second{first.arcIdentity,{old->identity,first.pointIdentities[1]},Q(1,2)};
    Trace trace{old->identity};
    const auto rebuilt=[&]{detail::ScopedRankTwoDiagnostics sink(observe,&trace);
        return rebuildRankTwoArcRegistryWithSubdivisions(a,initial,plan,{}, {first,second});}();
    std::cerr<<"RETAINED_REEVALUATIONS "<<trace.reevaluated<<" NEW_EVALUATIONS "<<trace.newEvaluated<<'\n';
    // A continuation that recalculates an already sealed exact coordinate
    // makes this fail before accepting any speed improvement. New roots must
    // still be actually computed; no global cache or input mutation is allowed.
    require(trace.reevaluated==0&&trace.newEvaluated==1,
        "MISSING_REUSE: continuation reevaluates an already factory-certified moving dependency");
    const auto fresh=buildRankTwoArcRegistryWithSubdivisions(a,plan,{}, {second,first});
    require(rebuilt.canonicalRecord()==fresh.canonicalRecord(),"reuse changed scientific registry vs fresh exact construction");
    const auto again=rebuildRankTwoArcRegistryWithSubdivisions(a,rebuilt,plan,{}, {first,second});
    require(again.canonicalRecord()==rebuilt.canonicalRecord(),"same-state reuse changed scientific output");
    const auto retained=std::find_if(rebuilt.points().begin(),rebuilt.points().end(),[&](const auto& p){return p.identity==old->identity;});
    require(retained!=rebuilt.points().end()&&retained->endpointIdentities==old->endpointIdentities
        &&retained->fraction==old->fraction,"reuse rewrote moving dependency");
    for(int d=0;d<2;++d)require(retained->coordinates[d].certificate().rootIdentity==old->coordinates[d].certificate().rootIdentity,
        "reuse lost original exact coordinate identity");
    auto changed=plan;changed.begin()->second.push_back(Q(1,5));
    rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,changed,{}, {first,second});},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,plan,{}, {second});},"INPUT_IDENTITY_MISMATCH");
    const auto foreign=analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"foreign"),Q(1));
    rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(foreign,initial,plan,{}, {first,second});},"INPUT_IDENTITY_MISMATCH");
    RankTwoSourceLimits tiny;tiny.maximumVertices=1;
    rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,plan,{}, {first,second},tiny);},"RESOURCE_LIMIT");
    RankTwoSourceLimits narrow;narrow.maximumIsolationBits=1;
    rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,plan,{}, {first,second},narrow);},"RESOURCE_LIMIT");
    std::ofstream out(path);require(bool(out),"cannot save actual reuse evidence");
    out<<"rank-two-registry-reuse-test-v1\nINITIAL_BEGIN\n"<<initial.canonicalRecord()<<"INITIAL_END\n"
        <<"FINAL_BEGIN\n"<<rebuilt.canonicalRecord()<<"FINAL_END\nretained_evaluations 0\nnew_evaluations 1\n"
        <<"whole_coverage_isotopy_thickness NOT_GRANTED\n";
    out.flush();require(bool(out),"incomplete reuse output");
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");
    test(argv[1]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
