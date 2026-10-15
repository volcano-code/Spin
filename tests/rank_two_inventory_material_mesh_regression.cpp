#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoInventoryMaterialMesh.h")
#include "spin/detail/RankTwoInventoryMaterialMesh.h"
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong failure class");return;}throw std::runtime_error("foreign inventory accepted");}
void run(const char* path){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
        plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto coverage=refineRankTwoFiniteArcCoverage(a,plan,{},Q(1,32));
    const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    RankTwoInteriorInventorySettings frozen;frozen.maximumPoints=99;
    const auto inventory=enumerateRankTwoMovingInteriors(a,frozen);
    const auto foreignSource=PreparedRankTwoSource::fromExact(source.exactCoefficients(),"foreign-inventory-material");
    const auto foreignA=analyzeRankTwoGlobalArrangement(foreignSource,Q(1));
    const auto foreignInventory=enumerateRankTwoMovingInteriors(foreignA,frozen);
    rejects([&]{buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,foreignInventory);},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{buildRankTwoAutomaticMaterialMeshFromInventory(foreignSource,a,base,inventory);},"INPUT_IDENTITY_MISMATCH");
    const auto wrongRadius=analyzeRankTwoGlobalArrangement(source,Q(2));
    const auto wrongInventory=enumerateRankTwoMovingInteriors(wrongRadius,frozen);
    rejects([&]{buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,wrongInventory);},"INPUT_IDENTITY_MISMATCH");
    RankTwoStrictMaterialMeshSettings noRepair;noRepair.maximumNewEdgePoints=0;
    rejects([&]{buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,inventory,{},noRepair);},"RESOURCE_LIMIT");
    const auto result=buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,inventory);
    require(result.interiors().scientificHash()==inventory.scientificHash()
        &&result.interiors().canonicalRecord()==inventory.canonicalRecord(),"inventory silently re-enumerated or policy replaced");
    require(result.baseMeshIdentity()==base.scientificHash(),"actual covered mesh replaced");
    const auto& initial=result.interiorMesh();const auto& final=result.strictMaterial().mesh();
    require(initial.points().size()==base.coverage().partition().registry().points().size()+inventory.points().size()
        &&initial.triangles().size()==base.triangulation().triangles().size()+2*inventory.points().size(),"inventory not actually inserted");
    std::set<std::string> ids;for(const auto& p:final.points())ids.insert(p.identity);
    for(const auto& p:inventory.points())require(ids.count(p.scientificHash()),"inventory point lost in final material");
    require(result.strictMaterial().lifts().allActiveStrictPassed(),"incomplete thickness validation");
    std::ofstream out(path);require(bool(out),"cannot save supplied-inventory material");
    out<<"rank-two-inventory-material-test-v1\nCASE full-disk AUTOMATIC "<<result.scientificHash()<<'\n'
       <<result.canonicalRecord()<<"CASE_END\n";out.flush();require(bool(out),"truncated supplied-inventory material");
    std::cout<<"PASS frozen inventory material points="<<inventory.points().size()<<" repairs="<<result.strictMaterial().insertedEdgePoints().size()<<std::endl;
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoInventoryMaterialMesh.h")
    if(argc!=2)throw std::invalid_argument("inventory material evidence path required");run(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual immutable owner inventory to strict material mesh");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
