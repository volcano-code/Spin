#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoRegistryRotation.h")
#include "spin/detail/RankTwoRegistryRotation.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& arrangement){
    spin::RankTwoArcSamplingPlan result;
    for(const auto& a:spin::arrangementRecords(arrangement).arcs){
        if(a.ends[0].location!=spin::RankTwoArcEndLocation::Finite
            || a.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        bool domain=false;for(const auto& owner:a.owners)domain|=owner.polynomial=="D";
        result[a.identity]={Q(1,2)};
        if(domain){result[a.identity].clear();for(int k=1;k<32;++k)result[a.identity].push_back(Q(k,32));}
    }
    return result;
}
void tests(const char* path){
    using namespace spin;
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"native-rotation-disk-v1");
    const auto diskArrangement=analyzeRankTwoGlobalArrangement(disk,Q(1));
    const auto diskRegistry=buildRankTwoArcRegistry(diskArrangement,plan(diskArrangement));
    const auto diskIncidence=certifyRankTwoChordIncidence(diskArrangement,diskRegistry);
    const auto diskProof=certifyRankTwoRegistryRotation(diskArrangement,diskRegistry,diskIncidence);
    require(diskProof.pointCount()==diskRegistry.points().size()
        && diskProof.finiteDirectionCount()==2*diskRegistry.edges().size()
        && diskProof.maximumFiniteDegree()==2 && diskProof.omittedInfiniteDirectionCount()==0,
        "circle finite-star coverage or native rotation was lost");
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"native-rotation-annulus-v1");
    const auto annulusArrangement=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    const auto annulusRegistry=buildRankTwoArcRegistry(annulusArrangement,plan(annulusArrangement));
    const auto annulusIncidence=certifyRankTwoChordIncidence(annulusArrangement,annulusRegistry);
    const auto annulusProof=certifyRankTwoRegistryRotation(annulusArrangement,annulusRegistry,annulusIncidence);
    require(annulusProof.finiteDirectionCount()==2*annulusRegistry.edges().size()
        && annulusProof.omittedInfiniteDirectionCount()==0,"annulus native stars omitted a component");
    // Frozen nontrivial star: D circle plus perpendicular native normal
    // boundaries. The origin has four actual finite directions; circle/line
    // events also retain native infinite directions rather than clipping them.
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto crossing=analyzeRankTwoNormalFace(disk,diskArrangement,cone);
    const auto crossingRegistry=buildRankTwoArcRegistry(crossing,plan(crossing));
    const auto crossingIncidence=certifyRankTwoChordIncidence(crossing,crossingRegistry);
    const auto crossingProof=certifyRankTwoRegistryRotation(crossing,crossingRegistry,crossingIncidence);
    require(crossingProof.maximumFiniteDegree()==4
        && crossingProof.omittedInfiniteDirectionCount()>0
        && crossingProof.pointCount()==crossingRegistry.points().size(),
        "perpendicular four-valent native star was reduced to a degree-two regression");
    const auto again=certifyRankTwoRegistryRotation(crossing,crossingRegistry,crossingIncidence);
    require(again.canonicalRecord()==crossingProof.canonicalRecord()
        && again.scientificHash()==crossingProof.scientificHash(),"rotation depends on prior predicates");
    rejects([&]{certifyRankTwoRegistryRotation(crossing,crossingRegistry,diskIncidence);},
        "INPUT_IDENTITY_MISMATCH","rotation reused foreign incidence authority");
    auto otherPlan=plan(diskArrangement);otherPlan.begin()->second.push_back(Q(1,3));
    const auto otherRegistry=buildRankTwoArcRegistry(diskArrangement,otherPlan);
    rejects([&]{certifyRankTwoRegistryRotation(diskArrangement,otherRegistry,diskIncidence);},
        "INPUT_IDENTITY_MISMATCH","changed sampling reused old incidence result");
    RankTwoRegistryRotationSettings empty;empty.maximumRotationSteps=0;
    rejects([&]{certifyRankTwoRegistryRotation(crossing,crossingRegistry,crossingIncidence,empty);},
        "RESOURCE_LIMIT","exhausted native rotation budget returned partial success");
    std::ofstream out(path);require(bool(out),"cannot create rotation evidence");
    out<<"rank-two-registry-rotation-test-v1\n";
    for(const auto* proof:{&diskProof,&annulusProof,&crossingProof})
        out<<"ROTATION "<<proof->scientificHash()<<'\n'<<proof->canonicalRecord();
    out.flush();require(bool(out),"incomplete rotation evidence");
    std::cout<<"PASS actual circle/annulus/four-valent native finite-star rotation, retained infinity and stale-authority rejection\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoRegistryRotation.h")
    if(argc!=2)throw std::invalid_argument("rotation evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual finite-star rotation agrees with authoritative native DCEL");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
