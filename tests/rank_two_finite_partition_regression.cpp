#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#if __has_include("spin/detail/RankTwoFinitePartition.h")
#include "spin/detail/RankTwoFinitePartition.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
Q eval(const std::array<Q,25>& p,const std::array<Q,2>& w){
    Q answer(0),x(1);for(int i=0;i<5;++i){Q y(1);for(int j=0;j<5;++j){answer+=p[5*i+j]*x*y;y*=w[1];}x*=w[0];}return answer;
}
void check(const spin::RankTwoArrangement& a,const char* name,std::ostream& out,
    bool needHole,bool needInfinite,bool needEmpty){
    using namespace spin;
    RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    const auto r=buildRankTwoArcRegistry(a,plan);
    const auto p=buildRankTwoFinitePartition(a,r);
    require(p.arrangementIdentity()==a.scientificHash() && p.registryIdentity()==r.scientificHash(),"partition identities lost");
    const auto& native=arrangementRecords(a);
    std::set<std::size_t> wanted,actual(p.regionIndices().begin(),p.regionIndices().end());
    bool hole=false,infinity=false,empty=false;std::set<std::size_t> points,edges;
    for(std::size_t i=0;i<r.regions().size();++i){
        const auto& region=r.regions()[i];
        const auto n=std::find_if(native.regions.begin(),native.regions.end(),[&](const auto& x){return x.identity==region.identity;});
        require(n!=native.regions.end(),"native region missing");
        bool selected=eval(native.originalFhD[2],n->witness)>0;
        for(const auto& normal:native.normalHalfplanes)selected=selected && eval(normal,n->witness)>0;
        if(!selected)continue;wanted.insert(i);
        require(!region.unbounded && region.kind!=RankTwoRegionKind::Outside,"unbounded/Outside cell inside finite domain");
        empty=empty || region.kind==RankTwoRegionKind::Empty;hole=hole || !region.innerLoops.empty();
        for(bool inner:{false,true}){
            const auto& loops=inner?region.innerLoops:region.outerLoops;
            for(std::size_t l=0;l<loops.size();++l){
                require(loops[l].finite,"infinite boundary retained in positive ellipse");
                const auto chain=r.loopVertices(region.identity,inner,l);
                require(chain.size()==loops[l].edges.size(),"hole boundary omitted");
                for(const auto& e:loops[l].edges)edges.insert(e.edgeIndex);
                points.insert(chain.begin(),chain.end());
            }
        }
        points.insert(region.isolatedPoints.begin(),region.isolatedPoints.end());
    }
    require(wanted==actual && !actual.empty(),"finite clipping omitted a domain cell");
    require(points==std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end()),"finite point inventory wrong");
    require(edges==std::set<std::size_t>(p.edgeIndices().begin(),p.edgeIndices().end()),"finite edge inventory wrong");
    require(p.excludedRegions().size()+actual.size()==r.regions().size(),"excluded regions have no proof");
    for(const auto& arc:r.arcs())if(!arc.finite){
        infinity=true;require(std::find(p.excludedArcs().begin(),p.excludedArcs().end(),arc.identity)!=p.excludedArcs().end(),
            "unbounded native curve truncated without exclusion witness");
    }
    require(!needHole || hole,"finite hole fixture not covered");
    require(!needInfinite || infinity,"unbounded branch fixture not covered");
    require(!needEmpty || empty,"Empty cells were dropped as if positive material only");
    require(reloadRankTwoFinitePartition(a,r,p.canonicalRecord()).scientificHash()==p.scientificHash(),"finite reload changed identity");
    auto bad=p.canonicalRecord();const auto at=bad.find("INCLUDE_REGION ");require(at!=std::string::npos,"no finite region record");
    bad.erase(at,bad.find('\n',at)-at+1);
    rejects([&]{reloadRankTwoFinitePartition(a,r,bad);},"INPUT_IDENTITY_MISMATCH","missing actual finite region accepted");
    bad=p.canonicalRecord();const auto h=bad.find("LOOP INNER ");
    if(h!=std::string::npos){bad.erase(h,bad.find('\n',h)-h+1);
        rejects([&]{reloadRankTwoFinitePartition(a,r,bad);},"INPUT_IDENTITY_MISMATCH","omitted hole accepted");}
    RankTwoSourceLimits noCells;noCells.maximumArrangementCells=0;
    rejects([&]{buildRankTwoFinitePartition(a,r,noCells);},"RESOURCE_LIMIT","empty clipping budget granted success");
    out<<"CASE "<<name<<" PARTITION "<<p.scientificHash()<<'\n'<<p.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" finite_cells="<<actual.size()<<" edges="<<edges.size()<<" holes="<<hole<<" Empty="<<empty<<'\n';
}
void tests(const char* path){
    using namespace spin;std::ostringstream saved;saved<<"rank-two-finite-partition-test-v1\n";
    const auto full=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-full-disk-v1");
    check(analyzeRankTwoGlobalArrangement(full,Q(1)),"full-disk",saved,false,false,false);
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-annulus-v1");
    const auto global=analyzeRankTwoGlobalArrangement(ann,Q(3,2));
    check(global,"annulus",saved,true,false,false);
    const auto lobes=PreparedRankTwoSource::fromExact(rank_two_cases::twoLobes(),"finite-two-lobes-v1");
    check(analyzeRankTwoGlobalArrangement(lobes,Q(3,2)),"two-lobes",saved,false,true,true);
    const auto equator=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"finite-equator-overlay-v1");
    const auto eq=analyzeRankTwoGlobalArrangement(equator,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto overlay=analyzeRankTwoNormalFace(equator,eq,cone);
    check(overlay,"equator-overlay",saved,false,true,false);
    RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(global).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    const auto registry=buildRankTwoArcRegistry(global,plan);
    rejects([&]{buildRankTwoFinitePartition(overlay,registry);},"INPUT_IDENTITY_MISMATCH","foreign arrangement accepted");
    std::ofstream out(path);require(bool(out),"cannot write finite partition evidence");out<<saved.str();out.flush();
    require(bool(out),"incomplete finite partition evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoFinitePartition.h")
    if(argc!=2)throw std::invalid_argument("finite partition evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact finite domain with all native holes and excluded infinity witnesses");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
