#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoCoveredFiniteMesh.h")
#include "spin/detail/RankTwoCoveredFiniteMesh.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
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
void check(const spin::RankTwoArrangement& a,const spin::RankTwoFiniteArcCoverage& coverage,
    std::ostream& out,const char* name,bool annulus){
    using namespace spin;
    // Catches triangulating a stale/coarse partition while attaching the new
    // curve proofs, filling holes by parity, and losing shared native edges.
    std::cerr<<"MESH_BEGIN "<<name<<'\n';
    const auto mesh=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    const auto& p=mesh.coverage().partition();const auto& r=p.registry();const auto& t=mesh.triangulation();
    require(mesh.coverage().scientificHash()==coverage.scientificHash()
        && t.partitionIdentity()==p.scientificHash() && t.registryIdentity()==r.scientificHash(),
        "covered mesh consumed a stale/unrefined partition");
    std::set<std::size_t> used;std::map<std::size_t,std::map<std::pair<std::size_t,std::size_t>,int>> boundary;
    for(const auto& triangle:t.triangles()){
        used.insert(triangle.points.begin(),triangle.points.end());
        for(std::size_t i=0;i<3;++i){const auto x=triangle.points[i],y=triangle.points[(i+1)%3];
            boundary[triangle.registryRegionIndex][std::minmax(x,y)]+=x<y?1:-1;}
    }
    require(used==std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end()),
        "refined boundary points omitted or anonymous mesh point accepted");
    std::size_t holes=0,empty=0,cut=0,full=0,consumed=0;
    for(const auto& summary:t.regionSummaries()){
        const auto& region=r.regions().at(summary.registryRegionIndex);
        holes+=summary.holes;empty+=region.kind==RankTwoRegionKind::Empty;
        cut+=region.kind==RankTwoRegionKind::Cut;full+=region.kind==RankTwoRegionKind::Full;
        require(summary.holes==region.innerLoops.size(),"regional mesh lost a native hole");
        std::map<std::pair<std::size_t,std::size_t>,int> wanted;
        const auto addLoop=[&](const RankTwoRegistryLoop& loop){
            for(const auto& ref:loop.edges){const auto& edge=r.edges().at(ref.edgeIndex);
                auto x=edge.points[0],y=edge.points[1];if(ref.reversed)std::swap(x,y);
                wanted[std::minmax(x,y)]+=x<y?1:-1;++consumed;
            }
        };
        for(const auto& loop:region.outerLoops)addLoop(loop);
        for(const auto& loop:region.innerLoops)addLoop(loop);
        auto actual=boundary.at(summary.registryRegionIndex);
        for(auto i=actual.begin();i!=actual.end();)if(i->second==0)i=actual.erase(i);else ++i;
        for(auto i=wanted.begin();i!=wanted.end();)if(i->second==0)i=wanted.erase(i);else ++i;
        require(actual==wanted,"triangles do not consume the same directed covered edges");
    }
    if(annulus)require(holes==3 && empty==1 && cut==2 && full==1,
        "annulus labels/holes replaced by parity filling");
    else {
        require(holes==0 && empty==0 && cut==0 && full==1 && used.size()==24,
            "full disk did not use actual error-controlled boundary");
        // The native D=0 equator is a true boundary, but boundary-only
        // triangulation still has zero lifted thickness throughout every
        // active triangle. Actual exact coordinates, not artificial zeros.
        for(const auto& triangle:t.triangles()){
            std::array<CORE::Expr,3> lo,hi;
            for(std::size_t i=0;i<3;++i){const auto& c=r.points().at(triangle.points[i]).coordinates;
                const auto D=CORE::Expr(1)-c[0].expression()*c[0].expression()-c[1].expression()*c[1].expression();
                require(D.sign()==0,"full disk vertex is not native equator");
                hi[i]=CORE::sqrt(D);lo[i]=-hi[i];}
            rejects([&]{detail::certifyRankTwoTriangleThickness(lo,hi,{true,true,true});},
                "SOURCE_TOPOLOGY_INVALID","zero-thickness actual disk mesh granted carrier");
        }
    }
    require(!t.triangles().empty() && consumed>=coverage.edgeCount(),"actual covered mesh not built");
    rejects([&]{triangulateRankTwoCoveredFiniteDomain(a,coverage,{}, {}, {0,1000000});},
        "RESOURCE_LIMIT","zero triangle budget granted covered mesh");
    out<<"CASE "<<name<<" COVERED_MESH "<<mesh.scientificHash()<<'\n'
        <<mesh.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" points="<<used.size()<<" triangles="<<t.triangles().size()
        <<" covered_edges="<<coverage.edgeCount()<<" holes="<<holes<<" labels="<<empty<<'/'<<cut<<'/'<<full<<'\n';
    std::cerr<<"MESH_END "<<name<<'\n';
}
void tests(const char* path){
    using namespace spin;std::ostringstream out;out<<"rank-two-covered-finite-mesh-test-v1\n";
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));
    const auto first=refineRankTwoFiniteArcCoverage(a,plan(a),{},Q(1,32));check(a,first,out,"full-disk",false);
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));
    const auto second=refineRankTwoFiniteArcCoverage(b,plan(b),{},Q(1,32));check(b,second,out,"annulus",true);
    rejects([&]{triangulateRankTwoCoveredFiniteDomain(b,first);},
        "INPUT_IDENTITY_MISMATCH","foreign coverage granted mesh");
    std::ofstream saved(path);require(bool(saved),"cannot save actual covered mesh");saved<<out.str();saved.flush();
    require(bool(saved),"incomplete covered mesh write");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoCoveredFiniteMesh.h")
    if(argc!=2)throw std::invalid_argument("covered mesh evidence path required");
    const auto start=std::chrono::steady_clock::now();
    const auto observe=+[](const char* phase,void* context)noexcept{
        try{const auto& t=*static_cast<const std::chrono::steady_clock::time_point*>(context);
            std::cerr<<"DIAGNOSTIC_NS "<<std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-t).count()<<' '<<phase<<'\n';}catch(...){}
    };
    spin::detail::ScopedRankTwoDiagnostics observation(observe,const_cast<std::chrono::steady_clock::time_point*>(&start));
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual common triangulation of the complete error-controlled finite domain");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
