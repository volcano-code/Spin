#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoInteriorInventory.h")
#include "spin/detail/RankTwoInteriorInventory.h"
#include "spin/detail/RankTwoStrictMaterialMesh.h"
#include "fixtures/rank_two_source_cases.h"
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong rejection class");return;}throw std::runtime_error("incomplete inventory accepted");}
void tests(const char* path){using namespace spin;
    std::ofstream out(path);require(bool(out),"cannot create inventory evidence");
    out<<"rank-two-interior-inventory-test-v1\n";
    for(bool annulus:{false,true}){
        const auto source=PreparedRankTwoSource::fromExact(annulus?rank_two_cases::annulus():RankTwoExactCoefficients{2,0,0,1,1,1,0,0,0,0},
            annulus?"finite-coverage-annulus-v1":"finite-coverage-disk-v1");
        const Q radius=annulus?Q(3,2):Q(1);
        const auto a=analyzeRankTwoGlobalArrangement(source,radius);
        const auto inventory=enumerateRankTwoMovingInteriors(a);
        require(inventory.arrangementIdentity()==a.scientificHash(),"foreign native authority");
        // Hand-derived native x-events and ALL vertical gaps, not just one
        // representative per face. Annulus: 0,1,3,5,3,1,0 across7slabs.
        require(inventory.slabCount()==(annulus?7u:1u)&&inventory.points().size()==(annulus?13u:1u),"active vertical gaps omitted or duplicated");
        std::map<std::pair<int,int>,std::size_t> counts;std::set<std::string> regions,pointIds;
        for(const auto& p:inventory.points()){
            const auto& r=p.definition();
            require(r.abscissaFraction==Q(1,2)&&r.ordinateFraction==Q(1,2),"inventory used an answer-selected fraction");
            require(p.arrangementIdentity()==a.scientificHash()&&pointIds.insert(p.scientificHash()).second,"foreign or duplicate internal point");
            const auto x=p.coordinates()[0].expression(),y=p.coordinates()[1].expression();
            const E h=annulus?x*x+y*y-2:-x*x-y*y-2;
            const E R=annulus?E(3)/2:E(1),D=R*R-x*x-y*y;
            require(D.sign()>0,"automatic point outside native disk");
            const E upper=CORE::sqrt(D),lower=p.regionKind()==RankTwoRegionKind::Full?-upper:h;
            require((upper-lower).sign()>0,"automatic active point lacks positive thickness");
            counts[{static_cast<int>(p.regionKind()),h.sign()}]++;regions.insert(r.regionIdentity);
        }
        require(counts[{3,-1}]==1&&(!annulus||(counts[{2,-1}]==4&&counts[{2,1}]==8)),"annulus Full/Cut sign regions or opposite-y branches omitted");
        require(regions.size()==(annulus?3u:1u),"active native face removed");
        RankTwoInteriorInventorySettings limited;limited.maximumPoints=inventory.points().size()-1;
        rejects([&]{enumerateRankTwoMovingInteriors(a,limited);},"RESOURCE_LIMIT");
        RankTwoInteriorInventorySettings noChecks;noChecks.maximumCandidateChecks=0;
        rejects([&]{enumerateRankTwoMovingInteriors(a,noChecks);},"RESOURCE_LIMIT");
        auto altered=inventory.canonicalRecord();altered+="caller_passed true\n";
        rejects([&]{reloadRankTwoMovingInteriors(a,altered);},"INPUT_IDENTITY_MISMATCH");
        out<<"CASE "<<(annulus?"annulus":"full-disk")<<" INVENTORY "<<inventory.scientificHash()<<'\n'<<inventory.canonicalRecord()<<"CASE_END\n";
        if(!annulus){
            RankTwoArcSamplingPlan plan;
            for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
            const auto covered=refineRankTwoFiniteArcCoverage(a,plan,{},Q(1,32));
            const auto base=triangulateRankTwoCoveredFiniteDomain(a,covered);
            const auto inserted=insertRankTwoMovingInteriorNodes(a,base,inventory.points());
            const auto strict=completeRankTwoMaterialThickness(source,a,inserted);
            require(strict.lifts().allActiveStrictPassed()&&!strict.mesh().triangles().empty(),"automatic native inventory did not enter actual material construction");
            out<<"DISK_STRICT_BEGIN\n"<<strict.canonicalRecord()<<"DISK_STRICT_END\n";
        }
        std::cout<<"PASS native inventory "<<(annulus?"annulus":"disk")<<" points="<<inventory.points().size()<<" slabs="<<inventory.slabCount()<<" active_regions="<<regions.size()<<'\n';
    }
    out.flush();require(bool(out),"incomplete inventory evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoInteriorInventory.h")
    if(argc!=2)throw std::invalid_argument("inventory evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: factory-owned ALL active native vertical-gap inventory");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
