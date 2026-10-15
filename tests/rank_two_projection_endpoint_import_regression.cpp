#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include "spin/detail/RankTwoEndpointRegistry.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "spin/detail/RankTwoRegistryContinuation.h"
#include "fixtures/rank_two_source_cases.h"

namespace {
using Q=CGAL::Gmpq;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
void test(const char* path){using namespace spin;
    // Catches rejection of a complete D star because its OUTSIDE native ray
    // is unbounded, and any fix that instead drops that ray or clips it.
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto whole=analyzeRankTwoArrangement(source,Q(3,2),octa());
    const auto master=buildRankTwoSharedFiniteSeams(source,whole);
    std::ofstream out(path);require(bool(out),"cannot save projection endpoint import");
    out<<"rank-two-projection-endpoint-import-v1\nMASTER_BEGIN\n"<<master.canonicalRecord()<<"MASTER_END\n";
    std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>> neighborhoods;
    for(std::size_t face=0;face<4;++face){const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        const auto base=buildRankTwoArcRegistryWithSharedSeams(a,plan,master,face);
        const auto domain=buildRankTwoFinitePartition(a,base);
        std::set<std::string> required;
        for(const auto e:domain.edgeIndices()){
            const auto arc=std::find_if(arrangementRecords(a).arcs.begin(),arrangementRecords(a).arcs.end(),
                [&](const auto& p){return p.identity==base.edges().at(e).arcIdentity;});
            require(arc!=arrangementRecords(a).arcs.end(),"finite edge lost actual native arc");
            for(const auto& end:arc->ends)if(end.vertexIdentity)required.insert(*end.vertexIdentity);
        }
        for(const auto& event:arrangementRecords(a).vertices){
            if(!required.count(event.identity))continue;
            const auto x=detail::CertifiedAlgebraicReal::reload(event.coordinates[0]).expression();
            const auto y=detail::CertifiedAlgebraicReal::reload(event.coordinates[1]).expression();
            if((x!=0&&y!=0)||4*(x*x+y*y)!=9)continue;
            const auto n=certifyRankTwoEndpointNeighborhood(a,event.identity);
            require(n.incidentArcs().size()==4&&n.fans().size()==2&&n.boundaryRays().size()==2,
                "original projection extreme lost complete four-direction native star");
            neighborhoods.emplace_back(face,n);
        }
    }
    require(neighborhoods.size()==8,"D extreme census changed");
    const auto publication=publishRankTwoSharedEndpointSeams(source,master,neighborhoods);
    require(publication.points().size()==4&&publication.bindings().size()==8,
        "D inward seam gate not shared by both actual owners");
    out<<"PUBLICATION_BEGIN\n"<<publication.canonicalRecord()<<"PUBLICATION_END\n";out.flush();
    std::size_t stars=0,excluded=0;
    for(std::size_t face=0;face<4;++face){const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        std::cerr<<"PROJECTION_IMPORT_BEGIN FACE "<<face<<'\n';
        const auto registry=buildRankTwoArcRegistryWithSharedEndpointSeams(a,plan,publication,face);
        const auto partition=buildRankTwoFinitePartition(a,registry);
        require(registry.arcs().size()==arrangementRecords(a).arcs.size()
            &&registry.regions().size()==arrangementRecords(a).regions.size(),"outside native DCEL records were deleted");
        std::vector<RankTwoEndpointNeighborhood> local;
        for(const auto& p:neighborhoods)if(p.first==face)local.push_back(p.second);
        require(registry.endpointNeighborhoods().size()==2,"D owner lost complete endpoint authority");
        for(const auto& n:registry.endpointNeighborhoods()){++stars;
            std::set<std::string> finiteBound;
            for(const auto& edge:registry.edges())if(edge.endpointNeighborhoodIdentity==n.scientificHash())finiteBound.insert(edge.arcIdentity);
            require(finiteBound.size()==3,"D star has wrong finite incidence count");
            for(const auto& ray:n.boundaryRays())if(ray.gate.fraction()<0){++excluded;
                const auto arc=std::find_if(registry.arcs().begin(),registry.arcs().end(),[&](const auto& p){return p.identity==ray.proof.arcIdentity;});
                require(arc!=registry.arcs().end()&&!arc->finite&&arc->edges.empty(),"outward original ray was removed or artificially clipped");
                for(const auto& point:registry.points())require(point.arcIdentity!=arc->identity,
                    "outside ray was inserted as an uncertified finite carrier point");
                for(const auto& region:registry.regions())if(region.identity==arc->leftRegion||region.identity==arc->rightRegion)
                    require(region.originalSigns[2]==-1,"outside exclusion used cone classification rather than original D sign");
            }
        }
        const auto again=rebuildRankTwoArcRegistryWithSubdivisions(a,registry,plan,local,{},{});
        require(again.canonicalRecord()==registry.canonicalRecord(),"continuation changed outside-incidence authority");
        bool rejected=false;try{rebuildRankTwoArcRegistryWithSubdivisions(a,registry,plan,{}, {},{});}
        catch(const RankTwoSourceError& e){rejected=e.code()=="INPUT_IDENTITY_MISMATCH";}
        require(rejected,"continuation dropped complete D star certificate");
        out<<"IMPORTED_FACE "<<face<<" REGISTRY_BEGIN\n"<<registry.canonicalRecord()
           <<"REGISTRY_END\nPARTITION_BEGIN\n"<<partition.canonicalRecord()<<"PARTITION_END\n";out.flush();
    }
    require(stars==8&&excluded==8,"some actual outside incidence omitted");
    out<<"D_STARS "<<stars<<" OUTSIDE_NATIVE_RAYS "<<excluded
       <<"\nwhole_finite_coverage NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nstrict_thickness NOT_GRANTED\n";
    std::cout<<"PASS eight complete D stars imported into four actual finite partitions; eight outside rays retained\n";
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");test(argv[1]);return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
