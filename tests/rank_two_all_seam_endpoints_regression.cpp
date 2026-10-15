#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include "spin/detail/RankTwoEndpointRegistry.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "spin/detail/RankTwoRegistryContinuation.h"
#include "fixtures/rank_two_source_cases.h"

namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
void tests(const char* path){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto a=analyzeRankTwoArrangement(source,Q(3,2),octa());
    const auto master=buildRankTwoSharedFiniteSeams(source,a);
    const auto& whole=wholeArrangementRecords(a);
    std::ofstream out(path);require(bool(out),"cannot save all endpoint evidence");
    out<<"rank-two-all-seam-endpoints-test-v1\nMASTER_BEGIN\n"<<master.canonicalRecord()<<"MASTER_END\n";out.flush();
    std::map<std::string,std::set<std::size_t>> wanted;
    // Both F circles, all four axes. Membership is the exact original sphere/
    // quadric intersection equation, not a selection of successful endpoints.
    for(const auto& point:master.points())if(point.endpoint){
        const auto x=point.coordinates[0].expression(),y=point.coordinates[1].expression();
        const auto u=x*x+y*y;
        if((x==0||y==0)&&4*u*u-12*u+7==0){
            auto& owners=wanted[point.identity];
            for(const auto& interval:master.intervals())if(interval.finite)
                for(const auto index:interval.points)if(master.points()[index].identity==point.identity)
                    owners.insert(interval.incidentFaces.begin(),interval.incidentFaces.end());
        }
    }
    require(wanted.size()==8,"original two-circle/eight shared F endpoints not preserved");
    std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>> neighborhoods;
    for(const auto& endpoint:wanted){
        require(endpoint.second.size()==2,"actual common endpoint lacks two owners");
        for(const auto face:endpoint.second){
            std::cerr<<"ALL_ENDPOINT_BEGIN "<<endpoint.first<<" FACE "<<face<<'\n';
            const auto n=certifyRankTwoEndpointNeighborhood(whole.normalFaces.at(face),endpoint.first);
            require(n.fans().size()==2&&n.boundaryRays().size()==2&&n.incidentArcs().size()==4,
                "actual circle/line four-direction star incomplete");
            neighborhoods.emplace_back(face,n);
            out<<"ENDPOINT "<<endpoint.first<<" FACE "<<face<<" NEIGHBORHOOD_BEGIN\n"
               <<n.canonicalRecord()<<"NEIGHBORHOOD_END\n";out.flush();
            std::cerr<<"ALL_ENDPOINT_END "<<endpoint.first<<" FACE "<<face<<'\n';
        }
    }
    const auto publication=publishRankTwoSharedEndpointSeams(source,master,neighborhoods);
    require(publication.points().size()==16&&publication.bindings().size()==32,
        "all common endpoint gates were not published together");
    require(std::count_if(publication.points().begin(),publication.points().end(),[](const auto& p){return p.masterAlias;})==4,
        "four original moving midpoint aliases were not retained structurally");
    for(const auto& p:publication.points())if(p.masterAlias){
        const auto original=std::find_if(master.points().begin(),master.points().end(),[&](const auto& q){return q.identity==p.identity;});
        require(original!=master.points().end()&&!original->endpoint,"structural alias invented or replaced original sample");
        for(const auto index:p.consumers){const auto& b=publication.bindings().at(index);
            const bool forward=b.endpointIdentity==original->parentIdentities[0]&&b.anchorIdentity==original->parentIdentities[1]
                &&b.fraction==original->fraction;
            const bool reversed=b.endpointIdentity==original->parentIdentities[1]&&b.anchorIdentity==original->parentIdentities[0]
                &&b.fraction==Q(1)-original->fraction;
            require(forward||reversed,"mere current coordinate equality was treated as moving definition equality");
        }
    }
    out<<"PUBLICATION_BEGIN\n"<<publication.canonicalRecord()<<"PUBLICATION_END\n";out.flush();
    std::size_t neighborhoodsUsed=0,bindingsUsed=0;
    for(std::size_t face=0;face<4;++face){
        const auto& local=whole.normalFaces.at(face);RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(local).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        std::cerr<<"ALL_IMPORT_BEGIN FACE "<<face<<'\n';
        const auto registry=buildRankTwoArcRegistryWithSharedEndpointSeams(local,plan,publication,face);
        const auto partition=buildRankTwoFinitePartition(local,registry);
        require(registry.endpointNeighborhoods().size()==4,"owner lost one of its four F endpoint stars");
        neighborhoodsUsed+=registry.endpointNeighborhoods().size();
        std::vector<RankTwoEndpointNeighborhood> localNeighborhoods;
        for(const auto& n:neighborhoods)if(n.first==face)localNeighborhoods.push_back(n.second);
        const auto again=rebuildRankTwoArcRegistryWithSubdivisions(local,registry,plan,localNeighborhoods,{},{});
        require(again.canonicalRecord()==registry.canonicalRecord(),"multi-endpoint continuation lost common authority");
        std::map<std::string,std::size_t> points;
        for(std::size_t k=0;k<registry.points().size();++k)points.emplace(registry.points()[k].identity,k);
        for(const auto& b:publication.bindings())if(b.face==face){
            const auto& expected=publication.points().at(b.pointIndex);
            require(points.count(expected.identity),"common gate was dropped or resampled by owner");
            const auto& actual=registry.points().at(points.at(expected.identity));
            require(actual.kind==RankTwoRegistryPointKind::SharedSeam&&actual.endpointArcIdentity==b.proofIdentity
                &&actual.endpointNeighborhoodIdentity==b.neighborhoodIdentity,"common gate lost actual owner proof");
            for(int d=0;d<2;++d)require(actual.coordinates[d].certificate().rootIdentity==expected.coordinates[d].certificate().rootIdentity,
                "owner uses different common root");++bindingsUsed;
        }
        out<<"IMPORTED_FACE "<<face<<" REGISTRY_BEGIN\n"<<registry.canonicalRecord()
           <<"REGISTRY_END\nPARTITION_BEGIN\n"<<partition.canonicalRecord()<<"PARTITION_END\n";out.flush();
        std::cerr<<"ALL_IMPORT_END FACE "<<face<<'\n';
    }
    require(neighborhoodsUsed==16&&bindingsUsed==32,"some required owner evidence omitted");
    // This sample coincides with an inner endpoint ray CURRENTLY, but the ray
    // uses another native anchor with the same projection, not the master
    // origin. The two definitions must not be aliased by root equality alone.
    RankTwoSharedFiniteSeamSettings coincident;coincident.fractions={Q(1,128),Q(127,128)};
    const auto different=buildRankTwoSharedFiniteSeams(source,a,coincident);
    bool collision=false;
    for(const auto& p:publication.points())if(!p.masterAlias)
        for(const auto& q:different.points())if(!q.endpoint
            &&p.coordinates[0].expression()==q.coordinates[0].expression()
            &&p.coordinates[1].expression()==q.coordinates[1].expression())collision=true;
    require(collision,"frozen negative input does not actually produce a current point collision");
    bool rejected=false;
    try{publishRankTwoSharedEndpointSeams(source,different,neighborhoods);}
    catch(const RankTwoSourceError& e){rejected=e.code()=="ARRANGEMENT_UNRESOLVED"
        &&std::string(e.what()).find("different moving definitions")!=std::string::npos;}
    require(rejected,"equal coordinates with different actual anchors were welded");
    out<<"all_eight_F_endpoint_stars PASSED\nall_four_actual_owner_imports PASSED\n"
       <<"full_finite_coverage NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\nstrict_thickness NOT_GRANTED\n";
    out.flush();require(bool(out),"incomplete all endpoint evidence");
    std::cout<<"PASS all8 shared F endpoints/16 actual neighborhoods/32 owner bindings in four real finite partitions\n";
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");
    tests(argv[1]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
