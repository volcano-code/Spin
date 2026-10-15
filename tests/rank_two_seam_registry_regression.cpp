#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcRegistrySeams.h")
#include "spin/detail/RankTwoArcRegistrySeams.h"
#include "spin/detail/RankTwoFiniteTriangulation.h"
#include "spin/detail/RankTwoChordSweep.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* s){if(!x)throw std::runtime_error(s);}
template<class F>void rejects(F f,const char* code){bool caught=false;try{f();}
    catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"foreign authority or exhausted budget was not rejected");}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& a){
    spin::RankTwoArcSamplingPlan p;for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite&&arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
            p.emplace(arc.identity,std::vector<Q>{Q(1,2)});return p;}
void check(const spin::PreparedRankTwoSource& s,Q radius,const char* name,bool mesh,std::ostream& out){
    using namespace spin;const auto a=analyzeRankTwoArrangement(s,radius,octa());
    const auto m=buildRankTwoSharedFiniteSeams(s,a);const auto& whole=wholeArrangementRecords(a);
    std::size_t consumers=0,triangles=0;out<<"CASE "<<name<<"\nMASTER_BEGIN\n"<<m.canonicalRecord()<<"MASTER_END\n";
    for(std::size_t face=0;face<whole.normalFaces.size();++face){
        const auto& local=whole.normalFaces[face];const auto p=plan(local);
        const auto r=buildRankTwoArcRegistryWithSharedSeams(local,p,m,face);
        const auto domain=buildRankTwoFinitePartition(local,r);
        require(r.arcs().size()==arrangementRecords(local).arcs.size(),"import dropped original outside or infinite arc");
        std::map<std::string,std::size_t> indices;for(std::size_t k=0;k<r.points().size();++k)indices.emplace(r.points()[k].identity,k);
        for(const auto& interval:m.intervals())if(interval.finite)for(int side=0;side<2;++side)if(interval.incidentFaces[side]==face){
            ++consumers;const auto master=m.consumerPoints(interval.identity,side);std::vector<std::size_t> chain;
            for(auto k:master){const auto& mp=m.points().at(k);require(indices.count(mp.identity),"actual face registry omitted master feature");
                auto index=indices.at(mp.identity);chain.push_back(index);const auto& point=r.points()[index];
                for(int d=0;d<2;++d)require(point.coordinates[d].certificate().rootIdentity==mp.coordinates[d].certificate().rootIdentity
                    &&point.coordinates[d].expression()==mp.coordinates[d].expression(),"consumer reconstructed different master coordinate");
                if(!mp.endpoint)require(point.kind==RankTwoRegistryPointKind::SharedSeam,"seam was independently sampled");}
            const RankTwoRegistryArc* arc=nullptr;for(const auto& row:r.arcs())if(row.identity==interval.sourceArcs[side])arc=&row;
            require(arc&&arc->edges.size()+1==chain.size(),"master subdivision not actually used by native arc edges");
            for(std::size_t k=1;k<chain.size();++k){const auto edge=arc->edges[k-1];
                require(r.edges()[edge].points==std::array<std::size_t,2>{chain[k-1],chain[k]},"actual seam edge has independently generated endpoint");
                require(std::find(domain.edgeIndices().begin(),domain.edgeIndices().end(),edge)!=domain.edgeIndices().end(),"finite partition omitted actual shared edge");}
        }
        out<<"FACE "<<face<<"\nREGISTRY_BEGIN\n"<<r.canonicalRecord()<<"REGISTRY_END\nPARTITION_BEGIN\n"<<domain.canonicalRecord()<<"PARTITION_END\n";
        if(mesh){const auto incidence=certifyRankTwoChordIncidenceSweep(local,r);
            const auto rotation=certifyRankTwoRegistryRotation(local,r,incidence);
            const auto t=triangulateRankTwoFinitePartition(local,domain,incidence,rotation);
            std::map<std::array<std::size_t,2>,std::size_t> uses;
            for(const auto& triangle:t.triangles())for(int k=0;k<3;++k){std::array<std::size_t,2> edge{triangle.points[k],triangle.points[(k+1)%3]};
                if(edge[0]>edge[1])std::swap(edge[0],edge[1]);++uses[edge];}
            for(const auto& interval:m.intervals())if(interval.finite)for(int side=0;side<2;++side)if(interval.incidentFaces[side]==face){
                auto chain=m.consumerPoints(interval.identity,side);for(std::size_t k=1;k<chain.size();++k){std::array<std::size_t,2> edge{
                    indices.at(m.points()[chain[k-1]].identity),indices.at(m.points()[chain[k]].identity)};
                    if(edge[0]>edge[1])std::swap(edge[0],edge[1]);require(uses[edge]==1,"triangulation did not consume actual master boundary exactly once");}}
            triangles+=t.triangles().size();out<<"MESH_BEGIN\n"<<t.canonicalRecord()<<"MESH_END\n";
        }
        auto missing=p;if(!missing.empty()){missing.erase(missing.begin());rejects([&]{buildRankTwoArcRegistryWithSharedSeams(local,missing,m,face);},"INPUT_IDENTITY_MISMATCH");}
        for(const auto& interval:m.intervals())if(interval.finite)for(int side=0;side<2;++side)if(interval.incidentFaces[side]==face){
            auto incompatible=p;incompatible.at(interval.sourceArcs[side])={Q(1,3)};
            rejects([&]{buildRankTwoArcRegistryWithSharedSeams(local,incompatible,m,face);},"INPUT_IDENTITY_MISMATCH");}
        RankTwoSourceLimits tiny;tiny.maximumVertices=1;rejects([&]{buildRankTwoArcRegistryWithSharedSeams(local,p,m,face,tiny);},"RESOURCE_LIMIT");
        rejects([&]{buildRankTwoArcRegistryWithSharedSeams(local,p,m,(face+1)%8);},"INPUT_IDENTITY_MISMATCH");
    }
    require(consumers==(mesh?8:32),"shared finite dual-consumer coverage incomplete");
    // Each quadrant has origin, two axis ends, two axis midpoints and one
    // circle point: six boundary vertices imply four triangles, four owners.
    if(mesh)require(triangles==16,"actual quadrant triangulation coverage incorrect");
    out<<"SUMMARY CONSUMERS "<<consumers<<" TRIANGLES "<<triangles<<"\nCASE_END\n";
    std::cout<<"PASS "<<name<<" consumers="<<consumers<<" actual_triangles="<<triangles<<'\n';
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcRegistrySeams.h")
    if(argc!=2)throw std::invalid_argument("actual registry evidence path required");std::ostringstream saved;
    saved<<"rank-two-seam-registry-test-v1\n";
    check(spin::PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1"),Q(3,2),"annulus",false,saved);
    check(spin::PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"shared-seams-full-disk-v1"),Q(1),"full-disk",true,saved);
    std::ofstream out(argv[1]);require(bool(out),"cannot save actual registry result");out<<saved.str();out.flush();require(bool(out),"incomplete result");return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: immutable master seams consumed by actual per-face arc registries and finite mesh");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
