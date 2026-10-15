#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoChordIncidence.h")
#include "spin/detail/RankTwoChordIncidence.h"
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& arrangement,bool skew=false,bool refineDomain=false){
    spin::RankTwoArcSamplingPlan result;
    for(const auto& a:spin::arrangementRecords(arrangement).arcs){
        if(a.ends[0].location!=spin::RankTwoArcEndLocation::Finite
            ||a.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        bool domain=false;for(const auto& owner:a.owners)domain|=owner.polynomial=="D";
        // Frozen before implementation: skewing the enclosing circle's
        // interior sample makes a long chord cut across the inner polygons.
        // Actual native curves remain disjoint; a DCEL pass cannot detect it.
        result[a.identity]={skew&&domain?Q(1,128):Q(1,2)};
        // A separate method input, NOT an upgrade of the coarse registry:
        // fixed 1/32 grid on each enclosing-circle arc, fixed before its
        // first run; all original events/arcs and holes remain unchanged.
        if(refineDomain&&domain){
            result[a.identity].clear();
            for(int k=1;k<32;++k)result[a.identity].push_back(Q(k,32));
        }
    }
    return result;
}
void tests(const char* path){
    using namespace spin;
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"chord-incidence-disk-v1");
    const auto diskArrangement=analyzeRankTwoGlobalArrangement(disk,Q(1));
    const auto diskRegistry=buildRankTwoArcRegistry(diskArrangement,plan(diskArrangement));
    const auto diskProof=certifyRankTwoChordIncidence(diskArrangement,diskRegistry);
    require(diskProof.registryIdentity()==diskRegistry.scientificHash()
        && diskProof.edgePairCount()==6 && diskProof.pointEdgeCount()==8
        && diskProof.pointPairCount()==6,"actual four-edge disk incidence omitted a pair");
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"chord-incidence-annulus-v1");
    const auto arrangement=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    const auto coarse=buildRankTwoArcRegistry(arrangement,plan(arrangement));
    // Measured true F/D coarse intersection: all-midpoint DCEL was never
    // an incidence-positive fixture. Retain it as a stronger negative.
    rejects([&]{certifyRankTwoChordIncidence(arrangement,coarse);},"ARRANGEMENT_UNRESOLVED",
        "all-midpoint native annulus was silently upgraded to noncrossing chords");
    const auto registry=buildRankTwoArcRegistry(arrangement,plan(arrangement,false,true));
    const auto proof=[&]{try{return certifyRankTwoChordIncidence(arrangement,registry);}
        catch(const RankTwoSourceError& e){
            std::cerr<<"FROZEN_ANNULUS_INPUT "<<registry.scientificHash()<<'\n'
                <<registry.canonicalRecord();throw;
        }}();
    require(proof.arrangementIdentity()==arrangement.scientificHash()
        && proof.registryIdentity()==registry.scientificHash(),"chord geometry detached from native authority");
    const auto n=registry.points().size(),m=registry.edges().size();
    require(proof.edgePairCount()==m*(m-1)/2 && proof.pointPairCount()==n*(n-1)/2
        && proof.pointEdgeCount()==m*(n-2),"quadratic incidence coverage was silently sampled");
    const auto again=certifyRankTwoChordIncidence(arrangement,registry);
    require(again.canonicalRecord()==proof.canonicalRecord() && again.scientificHash()==proof.scientificHash(),
        "incidence depends on previous predicate evaluation");
    const auto crossed=buildRankTwoArcRegistry(arrangement,plan(arrangement,true));
    rejects([&]{certifyRankTwoChordIncidence(arrangement,crossed);},"ARRANGEMENT_UNRESOLVED",
        "native arrangement pass hid crossing actual coarse chords");
    rejects([&]{certifyRankTwoChordIncidence(diskArrangement,registry);},"INPUT_IDENTITY_MISMATCH",
        "chord certificate accepted a foreign arrangement");
    RankTwoChordIncidenceSettings exhausted;exhausted.maximumIncidenceTests=0;
    rejects([&]{certifyRankTwoChordIncidence(arrangement,registry,exhausted);},"RESOURCE_LIMIT",
        "pair budget exhaustion returned a partial positive certificate");
    std::ofstream out(path);require(bool(out),"cannot open incidence evidence");
    out<<"rank-two-chord-incidence-test-v1\nINCIDENCE "<<diskProof.scientificHash()<<'\n'
        <<diskProof.canonicalRecord()<<"INCIDENCE "<<proof.scientificHash()<<'\n'<<proof.canonicalRecord();
    out.flush();require(bool(out),"incomplete incidence evidence");
    std::cout<<"PASS actual disk/annulus chord incidences, shared endpoints, holes, crossing coarse-chord rejection and complete pair budgets\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoChordIncidence.h")
    if(argc!=2)throw std::invalid_argument("incidence evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact all-pair incidence of actual shared registry chords");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
