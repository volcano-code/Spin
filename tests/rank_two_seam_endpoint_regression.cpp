#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <set>
#include "spin/detail/RankTwoArcRegistrySeams.h"
#include "spin/detail/RankTwoEndpointNeighborhood.h"
#include "spin/detail/RankTwoRegistryContinuation.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "fixtures/rank_two_source_cases.h"
#if __has_include("spin/detail/RankTwoSharedEndpointSeams.h")
#include "spin/detail/RankTwoSharedEndpointSeams.h"
#endif
#if __has_include("spin/detail/RankTwoEndpointRegistry.h")
#include "spin/detail/RankTwoEndpointRegistry.h"
#endif
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){bool caught=false;try{f();}
    catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"foreign/incomplete endpoint authority or budget was not rejected");}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto a=analyzeRankTwoArrangement(source,Q(3,2),octa());
    const auto master=buildRankTwoSharedFiniteSeams(source,a);
    const auto& whole=wholeArrangementRecords(a);
    std::string endpoint;
    for(const auto& p:master.points())if(p.endpoint&&p.coordinates[1].expression()==0
        &&p.coordinates[0].expression()>CORE::Expr(0)&&p.coordinates[0].expression()<CORE::Expr(1))endpoint=p.identity;
    require(!endpoint.empty(),"frozen inner quartic/normal-edge endpoint absent");
    std::ofstream saved(path);require(bool(saved),"cannot save endpoint evidence");
    saved<<"rank-two-seam-endpoint-test-v1\nENDPOINT "<<endpoint<<"\nMASTER_BEGIN\n"
         <<master.canonicalRecord()<<"MASTER_END\n";saved.flush();
    std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>> neighborhoods;
    // Positive x preimage belongs to the -x/north shared normal edge, whose
    // owners are faces 1 and 2. Input is frozen, not chosen from successful runs.
    for(const std::size_t face:{1,2}){
        std::cerr<<"ENDPOINT_BEGIN FACE "<<face<<'\n';
        const auto n=certifyRankTwoEndpointNeighborhood(whole.normalFaces.at(face),endpoint);
        require(n.fans().size()==2&&n.boundaryRays().size()==2&&n.incidentArcs().size()==4,
            "inner circle/normal-line complete four-direction star lost");
        saved<<"FACE "<<face<<" NEIGHBORHOOD "<<n.scientificHash()<<'\n'<<n.canonicalRecord()
             <<"NEIGHBORHOOD_END\n";saved.flush();
        neighborhoods.emplace_back(face,n);std::cerr<<"ENDPOINT_END FACE "<<face<<'\n';
    }
#if __has_include("spin/detail/RankTwoSharedEndpointSeams.h")
    const auto publication=publishRankTwoSharedEndpointSeams(source,master,neighborhoods);
    require(publication.bindings().size()==4,"both actual boundary rays on both owners must be published");
    require(publication.points().size()==2,"identical shared gate dependencies must own only two points");
    saved<<"PUBLICATION_BEGIN\n"<<publication.canonicalRecord()<<"PUBLICATION_END\n";
    require(publication.masterIdentity()==master.scientificHash(),"publication lost immutable original master");
    // Exact gate membership and shared identity must be checked before import;
    // no local graph proof or 3D acceptance is inferred from publication alone.
    for(const auto& point:publication.points()){
        require(point.coordinates[1].expression()==0,"shared gate left actual normal seam");
        require(point.consumers.size()==2,"point lacks either adjacent owner proof");
    }
    auto one=neighborhoods;one.pop_back();
    rejects([&]{publishRankTwoSharedEndpointSeams(source,master,one);},"ARRANGEMENT_UNRESOLVED");
    auto duplicate=neighborhoods;duplicate.push_back(duplicate.front());
    rejects([&]{publishRankTwoSharedEndpointSeams(source,master,duplicate);},"INPUT_IDENTITY_MISMATCH");
    auto wrong=neighborhoods;wrong[0].first=0;
    rejects([&]{publishRankTwoSharedEndpointSeams(source,master,wrong);},"INPUT_IDENTITY_MISMATCH");
    const auto foreign=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"foreign-frame-v1");
    rejects([&]{publishRankTwoSharedEndpointSeams(foreign,master,neighborhoods);},"INPUT_IDENTITY_MISMATCH");
    RankTwoSourceLimits tiny;tiny.maximumVertices=1;
    rejects([&]{publishRankTwoSharedEndpointSeams(source,master,neighborhoods,tiny);},"RESOURCE_LIMIT");
    auto disagree=neighborhoods;RankTwoEndpointNeighborhoodSettings narrower;narrower.heightMargin=Q(1,200);
    disagree[1].second=certifyRankTwoEndpointNeighborhood(whole.normalFaces[2],endpoint,narrower);
    rejects([&]{publishRankTwoSharedEndpointSeams(source,master,disagree);},"ARRANGEMENT_UNRESOLVED");
#if __has_include("spin/detail/RankTwoEndpointRegistry.h")
    for(const auto& item:neighborhoods){const auto face=item.first;const auto& n=item.second;
        const auto& local=whole.normalFaces[face];RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(local).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        const auto r=buildRankTwoArcRegistryWithSharedEndpointSeams(local,plan,publication,face);
        const auto partition=buildRankTwoFinitePartition(local,r);
        require(r.endpointNeighborhoods().size()==1&&r.endpointFans().size()==2,"actual registry lost complete neighborhood");
        require(r.arcs().size()==arrangementRecords(local).arcs.size(),"native outside or infinite arc dropped");
        std::map<std::string,std::size_t> ids;for(std::size_t k=0;k<r.points().size();++k)ids.emplace(r.points()[k].identity,k);
        std::set<std::string> bound;std::size_t shared=0;
        for(const auto& b:publication.bindings())if(b.face==face){const auto& p=publication.points()[b.pointIndex];
            require(ids.count(p.identity),"registry ignored common gate and independently resampled");
            const auto& actual=r.points()[ids.at(p.identity)];
            require(actual.endpointNeighborhoodIdentity==n.scientificHash()&&actual.endpointArcIdentity==b.proofIdentity,
                "common coordinate lost actual owner proof binding");
            for(int d=0;d<2;++d)require(actual.coordinates[d].expression()==p.coordinates[d].expression()
                &&actual.coordinates[d].certificate().rootIdentity==p.coordinates[d].certificate().rootIdentity,
                "consumer copied another common expression");++shared;
        }
        for(std::size_t k=0;k<r.edges().size();++k){const auto& e=r.edges()[k];if(e.endpointNeighborhoodIdentity.empty())continue;
            require(e.endpointNeighborhoodIdentity==n.scientificHash()&&bound.insert(e.arcIdentity).second,"duplicate/foreign endpoint chord");
            require(r.points()[e.points[0]].eventIdentity==endpoint||r.points()[e.points[1]].eventIdentity==endpoint,
                "endpoint authority assigned to wrong edge");
        }
        require(bound==std::set<std::string>(n.incidentArcs().begin(),n.incidentArcs().end())&&shared==2,
            "actual four-direction star or dual shared gates omitted");
        std::size_t masterConsumers=0;
        for(const auto& interval:master.intervals())if(interval.finite)for(int side=0;side<2;++side)
            if(interval.incidentFaces[side]==face){++masterConsumers;
                auto chain=master.consumerPoints(interval.identity,side);
                std::vector<std::string> expected;for(auto k:chain)expected.push_back(master.points()[k].identity);
                const auto arc=std::find_if(r.arcs().begin(),r.arcs().end(),[&](const auto& v){return v.identity==interval.sourceArcs[side];});
                require(arc!=r.arcs().end(),"augmented import lost original master arc");
                std::vector<std::string> actual;
                for(const auto e:arc->edges){const auto& edge=r.edges()[e];
                    if(actual.empty())actual.push_back(r.points()[edge.points[0]].identity);
                    actual.push_back(r.points()[edge.points[1]].identity);
                    require(std::find(partition.edgeIndices().begin(),partition.edgeIndices().end(),e)!=partition.edgeIndices().end(),
                        "actual finite partition omitted imported seam");
                }
                std::vector<std::string> base;for(const auto& id:actual)
                    if(std::find(expected.begin(),expected.end(),id)!=expected.end())base.push_back(id);
                require(base==expected,"augmentation omitted/reversed old master dependency");
                std::size_t extras=0;for(const auto& b:publication.bindings())if(b.face==face&&b.intervalIdentity==interval.identity)++extras;
                require(actual.size()==expected.size()+extras,"local or duplicate seam sample inserted outside publication");
            }
        require(masterConsumers==8,"some finite/outside-labeled master interval lost");
        const auto again=rebuildRankTwoArcRegistryWithSubdivisions(local,r,plan,{n},{},{});
        require(again.canonicalRecord()==r.canonicalRecord(),"continuation lost augmented sealed master");
        saved<<"IMPORTED_FACE "<<face<<" REGISTRY "<<r.scientificHash()<<"\nREGISTRY_BEGIN\n"<<r.canonicalRecord()
             <<"REGISTRY_END\nPARTITION_BEGIN\n"<<partition.canonicalRecord()<<"PARTITION_END\n";
        auto split=plan;for(const auto& ray:n.boundaryRays())split.at(ray.proof.arcIdentity).push_back(Q(1,4096));
        rejects([&]{buildRankTwoArcRegistryWithSharedEndpointSeams(local,split,publication,face);},"INPUT_IDENTITY_MISMATCH");
        rejects([&]{buildRankTwoArcRegistryWithSharedEndpointSeams(local,plan,publication,(face+1)%8);},"INPUT_IDENTITY_MISMATCH");
        rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(local,r,plan,{},{});},"INPUT_IDENTITY_MISMATCH");
    }
    std::cout<<"PASS two common gate points/four owner bindings consumed by both actual registries and complete four-direction stars\n";
#else
    throw std::runtime_error("MISSING_FEATURE: shared endpoint publication consumed by actual registry");
#endif
#else
    throw std::runtime_error("MISSING_FEATURE: proof-bound common endpoint seam publication");
#endif
    saved.flush();require(bool(saved),"incomplete endpoint evidence");
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("endpoint evidence path required");
    tests(argv[1]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
