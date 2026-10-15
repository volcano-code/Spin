#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoAutomaticMaterialMesh.h")
#include "spin/detail/RankTwoAutomaticMaterialMesh.h"
#include "fixtures/rank_two_source_cases.h"
#include <set>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong failure class");return;}throw std::runtime_error("incomplete automatic material accepted");}
void run(const char* path){using namespace spin;
    std::ofstream out(path);require(bool(out),"cannot save automatic material evidence");
    out<<"rank-two-automatic-material-mesh-test-v1\n";
    for(bool annulus:{false,true}){
        const auto source=PreparedRankTwoSource::fromExact(annulus?rank_two_cases::annulus():RankTwoExactCoefficients{2,0,0,1,1,1,0,0,0,0},
            annulus?"finite-coverage-annulus-v1":"finite-coverage-disk-v1");
        const auto a=analyzeRankTwoGlobalArrangement(source,annulus?Q(3,2):Q(1));
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
        const auto coverage=refineRankTwoFiniteArcCoverage(a,plan,{},Q(1,32));
        const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);
        RankTwoAutomaticMaterialMeshSettings limited;limited.interiors.maximumPoints=annulus?12:0;
        rejects([&]{buildRankTwoAutomaticMaterialMesh(source,a,base,limited);},"RESOURCE_LIMIT");
        const auto foreign=PreparedRankTwoSource::fromExact(source.exactCoefficients(),"foreign-automatic-source");
        rejects([&]{buildRankTwoAutomaticMaterialMesh(foreign,a,base);},"INPUT_IDENTITY_MISMATCH");
        if(!annulus){RankTwoAutomaticMaterialMeshSettings noRepair;noRepair.thickness.maximumNewEdgePoints=0;
            rejects([&]{buildRankTwoAutomaticMaterialMesh(source,a,base,noRepair);},"RESOURCE_LIMIT");}
        const auto result=buildRankTwoAutomaticMaterialMesh(source,a,base);
        const auto& inventory=result.interiors();const auto& initial=result.interiorMesh();
        const auto& strict=result.strictMaterial();const auto& final=strict.mesh();
        const std::size_t n=annulus?13:1;
        require(result.baseMeshIdentity()==base.scientificHash()&&inventory.points().size()==n,"foreign base or incomplete native inventory");
        require(initial.points().size()==base.coverage().partition().registry().points().size()+n
            &&initial.triangles().size()==base.triangulation().triangles().size()+2*n,"not ALL native interiors actually inserted");
        std::set<std::string> ids;for(const auto& p:final.points())ids.insert(p.identity);
        for(const auto& p:inventory.points())require(ids.count(p.scientificHash()),"complete inventory point absent from final material");
        require(strict.parentMeshIdentity()==initial.scientificHash()&&strict.lifts().meshIdentity()==final.scientificHash()
            &&strict.lifts().allActiveStrictPassed(),"automatic composition detached from actual material");
        const auto repairs=strict.insertedEdgePoints().size();
        require(final.points().size()==initial.points().size()+repairs
            &&final.triangles().size()==initial.triangles().size()+2*repairs,"incomplete BOTH-sided actual edge repair");
        std::size_t holes=0;std::set<std::size_t> active;
        for(const auto& r:final.regionSummaries())holes+=r.holes;
        for(const auto& t:strict.lifts().triangles())if(t.active){require(t.strictPassed,"invalid active triangle hidden");active.insert(t.registryRegionIndex);}
        require(holes==(annulus?3u:0u)&&active.size()==(annulus?3u:1u),"native rings or active owners omitted");
        out<<"CASE "<<(annulus?"annulus":"full-disk")<<" AUTOMATIC "<<result.scientificHash()<<'\n'
            <<result.canonicalRecord()<<"CASE_END\n";
        std::cout<<"PASS automatic material "<<(annulus?"annulus":"disk")<<" native_interiors="<<n
            <<" repairs="<<repairs<<" final_points="<<final.points().size()<<" final_triangles="<<final.triangles().size()<<std::endl;
    }
    out.flush();require(bool(out),"incomplete automatic material evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoAutomaticMaterialMesh.h")
    if(argc!=2)throw std::invalid_argument("automatic material evidence path required");run(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: ALL automatic native interior points to actual strict material");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
