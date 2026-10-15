#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoFiniteCellTopology.h")
#include "spin/detail/RankTwoFiniteCellTopology.h"
#include "fixtures/rank_two_source_cases.h"
#include <sstream>
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
            &&arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
            p.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    return p;
}
void tests(const char* path){
    using namespace spin;using namespace spin::detail;
    const std::vector<std::array<std::size_t,3>> disk{{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
    const std::vector<RankTwoComplexBoundary> boundary{{false,{0,1,2,3}}};
    const auto local=inspectRankTwoTriangleComplex(disk,boundary,1000000);
    require(local.components.size()==1&&local.links.size()==5&&local.components[0].holes==0,
        "actual disk links/component boundary not inspected");
    // Same oriented boundary sums and total Euler checks are not enough:
    // two disks touching only at one vertex are a non-manifold pinch.
    rejects([&]{inspectRankTwoTriangleComplex({{0,1,2},{0,3,4}},
        {{false,{0,1,2}},{false,{0,3,4}}},1000000);},"SOURCE_TOPOLOGY_INVALID","pinched vertex granted manifold");
    rejects([&]{inspectRankTwoTriangleComplex(disk,{{false,{0,3,2,1}}},1000000);},
        "SOURCE_TOPOLOGY_INVALID","reversed boundary granted correspondence");
    rejects([&]{inspectRankTwoTriangleComplex(disk,{},1000000);},
        "SOURCE_TOPOLOGY_INVALID","omitted native boundary accepted");
    rejects([&]{inspectRankTwoTriangleComplex({{0,1,2},{0,1,3}},
        {{false,{0,1,2}},{false,{0,1,3}}},1000000);},"SOURCE_TOPOLOGY_INVALID","same-direction shared edge accepted");
    rejects([&]{inspectRankTwoTriangleComplex(disk,boundary,0);},"RESOURCE_LIMIT","unbounded topology inspection");
    std::ostringstream out;out<<"rank-two-finite-cell-topology-test-v1\n";
    const auto save=[&](const RankTwoArrangement& a,const char* name,bool annulus){
        const auto c=refineRankTwoFiniteArcCoverage(a,plan(a),{},Q(1,32));
        const auto m=triangulateRankTwoCoveredFiniteDomain(a,c);
        const auto t=certifyRankTwoFiniteCellTopology(a,m);
        require(t.coveredMesh().scientificHash()==m.scientificHash(),"topology used different actual mesh");
        require(t.canonicalRecord().find("global_isotopy NOT_GRANTED\n")!=std::string::npos,
            "finite manifold inspection upgraded to isotopy");
        rejects([&]{certifyRankTwoFiniteCellTopology(a,m,{0});},"RESOURCE_LIMIT","zero topology budget granted");
        require(t.regionCount()==(annulus?4u:1u)&&t.nativeHoleCount()==(annulus?3u:0u),
            "actual native cells/holes omitted");
        out<<"CASE "<<name<<" CELL_TOPOLOGY "<<t.scientificHash()<<'\n'<<t.canonicalRecord()<<"CASE_END\n";
        std::cout<<"PASS "<<name<<" cells="<<t.regionCount()<<" holes="<<t.nativeHoleCount()<<'\n';
        return m;
    };
    const auto d=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(d,Q(1));const auto mesh=save(a,"full-disk",false);
    const auto s=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(s,Q(3,2));save(b,"annulus",true);
    rejects([&]{certifyRankTwoFiniteCellTopology(b,mesh);},"INPUT_IDENTITY_MISMATCH","foreign cell authority accepted");
    std::ofstream file(path);require(bool(file),"cannot save finite cell topology evidence");
    file<<out.str();file.flush();require(bool(file),"incomplete topology evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoFiniteCellTopology.h")
    if(argc!=2)throw std::invalid_argument("cell topology evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual finite vertex links and oriented native cell correspondence");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
