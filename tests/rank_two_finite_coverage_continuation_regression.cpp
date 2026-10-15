#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include "spin/detail/RankTwoFiniteArcCoverage.h"
#include "spin/detail/RankTwoRegistryContinuation.h"
#ifdef RANK_TWO_VERIFY_REFINEMENT_REUSE
#include "spin/detail/RankTwoDiagnostics.h"
#endif
#if __has_include("spin/detail/RankTwoFiniteCoverageContinuation.h")
#include "spin/detail/RankTwoFiniteCoverageContinuation.h"
#endif

namespace {
using Q=CGAL::Gmpq;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){bool caught=false;try{f();}
    catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"finite continuation lost identity/resource validation");}
void tests(const char* path){using namespace spin;
    const auto s=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(s,Q(1));
    RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto base=buildRankTwoArcRegistry(a,plan);
    const RankTwoRegistryEdge* seedEdge=nullptr;
    for(const auto& edge:base.edges())if(base.points()[edge.points[0]].kind==RankTwoRegistryPointKind::ArcSample
        &&base.points()[edge.points[1]].kind==RankTwoRegistryPointKind::ArcSample){seedEdge=&edge;break;}
    require(seedEdge,"actual finite curved seed interval absent");
    const RankTwoArcSubdivision seed{seedEdge->arcIdentity,{base.points()[seedEdge->points[0]].identity,
        base.points()[seedEdge->points[1]].identity},Q(1,3)};
    const auto initial=rebuildRankTwoArcRegistryWithSubdivisions(a,base,plan,{}, {seed});
    const RankTwoRegistryPoint* point=nullptr;
    for(const auto& p:initial.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision)point=&p;
    require(point&&initial.subdivisions().size()==1,"frozen one-third moving dependency absent");
    std::ofstream out(path);require(bool(out),"cannot save continuation evidence");
    out<<"rank-two-finite-coverage-continuation-test-v1\nSEED "<<point->identity
       <<"\nINITIAL_BEGIN\n"<<initial.canonicalRecord()<<"INITIAL_END\n";out.flush();
#if __has_include("spin/detail/RankTwoFiniteCoverageContinuation.h")
#ifdef RANK_TWO_VERIFY_REFINEMENT_REUSE
    std::map<std::string,std::size_t> evaluations;
    detail::ScopedRankTwoDiagnostics sink([](const char* p,void* raw)noexcept{
        const std::string phase(p),prefix="registry:subdivision:",suffix=":coordinates-ready";
        if(phase.rfind(prefix,0)==0&&phase.size()>prefix.size()+suffix.size()
            &&phase.compare(phase.size()-suffix.size(),suffix.size(),suffix)==0)
            ++(*static_cast<std::map<std::string,std::size_t>*>(raw))[phase.substr(prefix.size(),phase.size()-prefix.size()-suffix.size())];
    },&evaluations);
#endif
    const auto result=refineRankTwoFiniteArcCoverageFromRegistry(a,initial,plan,{},Q(1,32));
#else
    // Exercise the real old handoff. It reconstructs a sampling-only registry
    // and loses a pre-existing moving dependency before whole-domain coverage.
    const auto result=refineRankTwoFiniteArcCoverage(a,plan,{},Q(1,32));
#endif
    const auto& r=result.partition().registry();
    const auto retained=std::find_if(r.points().begin(),r.points().end(),[&](const auto& p){return p.identity==point->identity;});
    require(retained!=r.points().end(),"MISSING_FEATURE: whole finite coverage discarded the actual input registry dependency");
    require(retained->endpointIdentities==seed.pointIdentities&&retained->fraction==Q(1,3),"retained dependency was rewritten");
    for(int k=0;k<2;++k)require(retained->coordinates[k].certificate().rootIdentity==point->coordinates[k].certificate().rootIdentity,
        "finite continuation silently resampled a different root");
    require(result.refinement().rounds().front().registryIdentity==initial.scientificHash(),"coverage did not start from actual registry");
    require(result.refinement().rounds().back().splitEdges==0,"uncompleted coverage accepted");
    std::set<std::string> covered,actual;
    for(const auto& proof:result.refinement().proofs()){
        require(proof.registryIdentity()==r.scientificHash()&&proof.strip().chordDistanceUpperBound<=Q(1,32),"obsolete final witness");
        require(covered.insert(proof.edgeIdentity()).second,"duplicate strip witness");
    }
    for(const auto& edge:result.refinement().endpointEdges())require(covered.insert(edge).second,"overlapping endpoint witness");
    for(const auto index:result.partition().edgeIndices())actual.insert(r.edges()[index].identity);
    require(covered==actual&&!actual.empty(),"not all actual finite edges proved");
    require(result.refinement().registry().scientificHash()==r.scientificHash(),"proofs and partition use different final maps");
#ifdef RANK_TWO_VERIFY_REFINEMENT_REUSE
    std::size_t repeated=0;for(const auto& entry:evaluations)if(entry.second!=1)++repeated;
    std::cerr<<"REFINEMENT_NEW_COORDINATES "<<evaluations.size()<<" REPEATED "<<repeated<<'\n';
    require(!evaluations.count(point->identity)&&repeated==0&&evaluations.size()+initial.subdivisions().size()==r.subdivisions().size(),
        "MISSING_REFINEMENT_REUSE: whole refiner rebuilds already certified subdivision roots in later rounds");
    out<<"REFINEMENT_REUSE new_coordinates "<<evaluations.size()<<" repeated 0\n";
#endif
    out<<"COVERAGE_BEGIN\n"<<result.canonicalRecord()<<"COVERAGE_END\n";
#if __has_include("spin/detail/RankTwoFiniteCoverageContinuation.h")
    auto changed=plan;changed.begin()->second.push_back(Q(1,5));
    rejects([&]{refineRankTwoFiniteArcCoverageFromRegistry(a,initial,changed,{},Q(1,32));},"INPUT_IDENTITY_MISMATCH");
    const auto other=analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"other-frame"),Q(1));
    rejects([&]{refineRankTwoFiniteArcCoverageFromRegistry(other,initial,plan,{},Q(1,32));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{refineRankTwoFiniteArcCoverageFromRegistry(a,initial,plan,{},Q(0));},"INPUT_IDENTITY_MISMATCH");
    RankTwoSourceLimits noDepth;noDepth.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoFiniteArcCoverageFromRegistry(a,initial,plan,{},Q(1,1048576),{},noDepth);},"RESOURCE_LIMIT");
#endif
    out.flush();require(bool(out),"incomplete continuation evidence");
    std::cout<<"PASS actual one-third moving dependency retained; whole finite coverage "
             <<actual.size()<<" edges, "<<r.subdivisions().size()<<" moving subdivisions\n";
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");
    tests(argv[1]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
