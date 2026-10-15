#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoRegistryContinuation.h")
#include "spin/detail/RankTwoRegistryContinuation.h"
#include "spin/detail/RankTwoArcRegistrySeams.h"
#include "spin/detail/RankTwoMixedArcIntervals.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "spin/detail/RankTwoArcMembership.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){bool caught=false;try{f();}
    catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"continuation accepted a foreign authority or local seam modification");}
void checkMaster(const spin::RankTwoArcRegistry& r,const spin::RankTwoSharedFiniteSeams& m,std::size_t face){
    using namespace spin;std::map<std::string,std::size_t> ids;
    for(std::size_t k=0;k<r.points().size();++k)ids.emplace(r.points()[k].identity,k);
    for(const auto& i:m.intervals())if(i.finite)for(int side=0;side<2;++side)if(i.incidentFaces[side]==face){
        const auto chain=m.consumerPoints(i.identity,side);const RankTwoRegistryArc* arc=nullptr;
        for(const auto& row:r.arcs())if(row.identity==i.sourceArcs[side])arc=&row;
        require(arc&&arc->edges.size()+1==chain.size(),"refinement changed authoritative seam subdivision");
        for(std::size_t k=0;k<chain.size();++k){const auto& p=m.points().at(chain[k]);
            require(ids.count(p.identity),"refinement discarded immutable master point");
            const auto& actual=r.points().at(ids.at(p.identity));
            for(int d=0;d<2;++d)require(actual.coordinates[d].expression()==p.coordinates[d].expression()
                &&actual.coordinates[d].certificate().rootIdentity==p.coordinates[d].certificate().rootIdentity,
                "refinement independently reconstructed master coordinate");
            if(!p.endpoint)require(actual.kind==RankTwoRegistryPointKind::SharedSeam,"master sample identity lost");
            if(k)require(r.edges().at(arc->edges[k-1]).points==std::array<std::size_t,2>{
                ids.at(m.points().at(chain[k-1]).identity),ids.at(p.identity)},"actual native seam edge differs");
        }
    }
}
void test(const char* path){using namespace spin;
    const auto s=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto normals=PreparedRankTwoNormalMesh::fromExact(
        {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
        {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});
    const auto whole=analyzeRankTwoArrangement(s,Q(3,2),normals);
    const auto master=buildRankTwoSharedFiniteSeams(s,whole);
    const auto& faces=wholeArrangementRecords(whole).normalFaces;
    std::ostringstream out;out<<"rank-two-seam-refinement-test-v1\nMASTER_BEGIN\n"<<master.canonicalRecord()<<"MASTER_END\n";
    std::size_t changed=0;
    // Same frozen annulus and octahedron. Every north owner keeps its actual
    // master radial seams while ALL curved middle intervals are refined.
    for(std::size_t face=0;face<4;++face){const auto& a=faces.at(face);std::set<std::string> seams;
        for(const auto& i:master.intervals())if(i.finite)for(int side=0;side<2;++side)
            if(i.incidentFaces[side]==face)seams.insert(i.sourceArcs[side]);
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,seams.count(arc.identity)?std::vector<Q>{Q(1,2)}:std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
        const auto initial=buildRankTwoArcRegistryWithSharedSeams(a,plan,master,face);checkMaster(initial,master,face);
        const auto domain=buildRankTwoFinitePartition(a,initial);std::set<std::string> domainArcs;
        for(auto edge:domain.edgeIndices())domainArcs.insert(initial.edges().at(edge).arcIdentity);
        std::vector<RankTwoMixedArcInterval> intervals;
        for(const auto& arc:initial.arcs())if(arc.finite&&domainArcs.count(arc.identity)&&!seams.count(arc.identity)){
            std::array<std::string,2> ends;
            for(const auto& p:initial.points())if(p.arcIdentity==arc.identity&&p.kind==RankTwoRegistryPointKind::ArcSample){
                if(p.fraction==Q(1,4))ends[0]=p.identity;if(p.fraction==Q(3,4))ends[1]=p.identity;}
            require(!ends[0].empty()&&!ends[1].empty(),"frozen curved interval omitted");intervals.push_back({arc.identity,ends});}
        // Four analytic circles, but the global quartic's inner critical
        // abscissa ALSO splits the outer F quarter into two native branches.
        require(intervals.size()==5,"annulus has five curved native intervals per north owner");
        const auto refined=refineRankTwoMixedArcIntervalsFromRegistry(a,initial,plan,{},intervals,Q(1,32));
        const auto& r=refined.registry();checkMaster(r,master,face);
        require(refined.rounds().size()>1&&!r.subdivisions().empty(),"only initialized, did not refine actual curved arcs");
        require(r.arcs().size()==initial.arcs().size(),"continuation discarded outside/native arcs");
        std::set<std::string> covered;
        for(const auto& p:refined.proofs())require(p.registryIdentity()==r.scientificHash()
            &&p.strip().chordDistanceUpperBound<=Q(1,32)&&covered.insert(p.edgeIdentity()).second,
            "curve proof reused obsolete registry or relaxed error");
        require(refined.rounds().back().splitEdges==0&&covered.size()==refined.rounds().back().selectedEdges,
            "selected curved interval has an unproved edge");
        for(const auto& p:r.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision){
            require(!seams.count(p.arcIdentity)&&detail::isOnRankTwoNativeArc(a,p.arcIdentity,p.coordinates),
                "refinement split a master seam or left native curve");++changed;}
        const auto rebuilt=rebuildRankTwoArcRegistryWithSubdivisions(a,r,plan,{},r.subdivisions());
        require(rebuilt.canonicalRecord()==r.canonicalRecord(),"sealed continuation changed actual scientific output");
        rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,r,plan,{},{});},"INPUT_IDENTITY_MISMATCH");
        rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(faces.at((face+1)%4),initial,plan,{},{});},"INPUT_IDENTITY_MISMATCH");
        auto wrong=plan;wrong.at(*seams.begin())={Q(1,3)};
        rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,wrong,{},{});},"INPUT_IDENTITY_MISMATCH");
        const auto seam=std::find_if(initial.arcs().begin(),initial.arcs().end(),[&](const auto& arc){return seams.count(arc.identity);});
        const auto& edge=initial.edges().at(seam->edges.front());
        rejects([&]{rebuildRankTwoArcRegistryWithSubdivisions(a,initial,plan,{},{{seam->identity,
            {initial.points().at(edge.points[0]).identity,initial.points().at(edge.points[1]).identity},Q(1,2)}});},"ARRANGEMENT_UNRESOLVED");
        out<<"FACE "<<face<<"\nINITIAL_BEGIN\n"<<initial.canonicalRecord()<<"INITIAL_END\nREFINEMENT_BEGIN\n"
            <<refined.canonicalRecord()<<"REFINEMENT_END\n";
        std::cout<<"PASS owner="<<face<<" rounds="<<refined.rounds().size()<<" curved_edges="<<covered.size()<<'\n';
    }
    require(changed>0,"no actual moving subdivisions published");out<<"SCOPE CURVED_MIDDLE_INTERVALS_ONLY\nGLOBAL_ISOTOPY NOT_GRANTED\n";
    std::ofstream saved(path);require(bool(saved),"cannot save actual refinement");saved<<out.str();saved.flush();require(bool(saved),"incomplete science output");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoRegistryContinuation.h")
    if(argc!=2)throw std::invalid_argument("evidence output required");test(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: immutable seam authority preserved through actual curve refinement rounds");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
