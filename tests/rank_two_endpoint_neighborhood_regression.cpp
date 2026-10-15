#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoEndpointNeighborhood.h")
#include "spin/detail/RankTwoEndpointNeighborhood.h"
#include "spin/detail/RankTwoNormalOverlay.h"
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
void tests(const char* path){
    using namespace spin;
    // Existing equator/tangency fixture, not selected from successful outputs.
    const auto p=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"endpoint-neighborhood-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(p,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};
    cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(p,global,cone);
    std::string event;
    for(const auto& v:arrangementRecords(a).vertices)
        if(C::reload(v.coordinates[0]).expression()==1 && C::reload(v.coordinates[1]).expression()==0)event=v.identity;
    require(!event.empty(),"equator endpoint absent");
    const auto n=certifyRankTwoEndpointNeighborhood(a,event);
    require(n.arrangementIdentity()==a.scientificHash() && n.endpointIdentity()==event,
        "neighborhood lost actual arrangement/event authority");
    require(n.projectionAxis()==1 && n.fans().size()==2 && n.boundaryRays().size()==2,
        "equator needs both three-graph half-sectors and both fixed boundary rays");
    require(n.incidentArcs().size()==8,"complete eight-direction native endpoint star was not retained");
    std::set<std::string> covered;
    for(const auto& f:n.fans()){
        require(f.orderedArcs().size()==3 && f.supportRecords().size()==6,
            "incident factor or normal owner was omitted from a half-sector");
        for(const auto& c:f.orderedArcs())require(covered.insert(c.arcIdentity).second,
            "same native arc counted twice across half-sectors");
    }
    for(const auto& r:n.boundaryRays()){
        require(r.proof.endpointIdentity==event && r.proof.arrangementIdentity==a.scientificHash()
            && r.proof.strip.chordDistanceUpperBound==0,
            "unchanged straight boundary ray is not actually certified");
        require(covered.insert(r.proof.arcIdentity).second,"boundary ray double-counted");
    }
    require(covered==std::set<std::string>(n.incidentArcs().begin(),n.incidentArcs().end()),
        "neighborhood passed while some native incidence remained uncovered");
    // Catch an evaluated/rejected geometric candidate silently disappearing
    // from the recorded budget ledger (including pre-arc height rejections).
    std::istringstream recorded(n.canonicalRecord());std::string line;
    std::size_t declared=0,listed=0;
    while(std::getline(recorded,line)){
        if(line.rfind("CANDIDATES ",0)==0)declared=std::stoull(line.substr(11));
        if(line.rfind("TRIAL ",0)==0)++listed;
    }
    require(declared==listed && declared>0,"evaluated candidate missing from endpoint budget ledger");
    const auto again=certifyRankTwoEndpointNeighborhood(a,event);
    require(again.canonicalRecord()==n.canonicalRecord() && again.scientificHash()==n.scientificHash(),
        "same-process proof depends on prior root refinement");
    rejects([&]{certifyRankTwoEndpointNeighborhood(a,"missing");},"INPUT_IDENTITY_MISMATCH",
        "unknown endpoint got a coverage certificate");
    RankTwoEndpointNeighborhoodSettings zero;zero.maximumCandidateTests=0;
    rejects([&]{certifyRankTwoEndpointNeighborhood(a,event,zero);},"RESOURCE_LIMIT",
        "zero candidate budget returned partial positive coverage");
    RankTwoEndpointNeighborhoodSettings coarse;coarse.maximumGateHalvings=1;
    rejects([&]{certifyRankTwoEndpointNeighborhood(a,event,coarse);},"ARRANGEMENT_UNRESOLVED",
        "unfinished coarse half-sectors received complete coverage");
    RankTwoEndpointNeighborhoodSettings invalid;invalid.maximumChordDistance=-1;
    rejects([&]{certifyRankTwoEndpointNeighborhood(a,event,invalid);},"INPUT_IDENTITY_MISMATCH",
        "negative error contract was accepted");
    std::ofstream out(path);require(bool(out),"cannot save native neighborhood evidence");
    out<<"rank-two-endpoint-neighborhood-test-v1\nNEIGHBORHOOD "<<n.scientificHash()<<'\n'<<n.canonicalRecord();
    out.flush();require(bool(out),"incomplete neighborhood evidence");
    std::cout<<"PASS complete equator endpoint neighborhood: all eight native incidences, both ordered half-sectors, two unchanged straight rays, strict error and budget rejection\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoEndpointNeighborhood.h")
    if(argc!=2)throw std::invalid_argument("neighborhood evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: complete actual endpoint neighborhood, including start-boundary rays");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
