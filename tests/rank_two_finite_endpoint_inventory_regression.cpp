#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include "spin/detail/RankTwoEndpointRegistry.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "fixtures/rank_two_source_cases.h"

namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
void tests(const char* path){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto whole=analyzeRankTwoArrangement(source,Q(3,2),octa());
    const auto master=buildRankTwoSharedFiniteSeams(source,whole);
    std::ofstream out(path);require(bool(out),"cannot save finite endpoint inventory");
    out<<"rank-two-finite-endpoint-inventory-v1\nMASTER_BEGIN\n"<<master.canonicalRecord()<<"MASTER_END\n";
    std::size_t required=0,passed=0,unknown=0,extremeStars=0,reflectedGates=0;
    for(std::size_t face=0;face<4;++face){
        const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        const auto registry=buildRankTwoArcRegistryWithSharedSeams(a,plan,master,face);
        const auto partition=buildRankTwoFinitePartition(a,registry);
        std::set<std::string> endpoints;
        for(const auto index:partition.edgeIndices()){
            const auto& edge=registry.edges().at(index);
            const auto arc=std::find_if(arrangementRecords(a).arcs.begin(),arrangementRecords(a).arcs.end(),
                [&](const auto& v){return v.identity==edge.arcIdentity;});
            require(arc!=arrangementRecords(a).arcs.end(),"finite edge lost native arc");
            for(const auto& end:arc->ends)if(end.vertexIdentity)endpoints.insert(*end.vertexIdentity);
        }
        require(!endpoints.empty(),"finite boundary endpoint inventory empty");
        out<<"FACE "<<face<<" REQUIRED "<<endpoints.size()<<"\nPARTITION_BEGIN\n"
           <<partition.canonicalRecord()<<"PARTITION_END\n";out.flush();
        for(const auto& id:endpoints){++required;
            out<<"ENDPOINT "<<id<<" FACE "<<face<<'\n';
            std::cerr<<"FINITE_ENDPOINT_BEGIN "<<id<<" FACE "<<face<<'\n';
            try{
                const auto n=certifyRankTwoEndpointNeighborhood(a,id);
                const auto event=std::find_if(arrangementRecords(a).vertices.begin(),arrangementRecords(a).vertices.end(),
                    [&](const auto& v){return v.identity==id;});
                require(event!=arrangementRecords(a).vertices.end(),"certified star lost actual native event");
                const auto x=detail::CertifiedAlgebraicReal::reload(event->coordinates[0]).expression();
                const auto y=detail::CertifiedAlgebraicReal::reload(event->coordinates[1]).expression();
                const bool extreme=(x==0||y==0)&&4*(x*x+y*y)==9;
                std::size_t reflected=0;std::set<std::string> incident;
                for(const auto& fan:n.fans())for(const auto& p:fan.orderedArcs())
                    require(incident.insert(p.arcIdentity).second,"duplicated native fan direction");
                for(const auto& ray:n.boundaryRays()){
                    require(incident.insert(ray.proof.arcIdentity).second,"duplicated native ray direction");
                    if(ray.gate.fraction()<0){++reflected;++reflectedGates;
                        const auto anchor=std::find_if(arrangementRecords(a).vertices.begin(),arrangementRecords(a).vertices.end(),
                            [&](const auto& v){return v.identity==ray.gate.anchorIdentity();});
                        require(anchor!=arrangementRecords(a).vertices.end(),"reflected ruler has no actual anchor");
                        const int axis=ray.gate.projectionAxis();
                        const auto start=detail::CertifiedAlgebraicReal::reload(event->coordinates[axis]);
                        const auto parent=detail::CertifiedAlgebraicReal::reload(anchor->coordinates[axis]);
                        const auto weight=ray.gate.fraction();std::ostringstream q;q<<weight;
                        const CORE::Expr f(CORE::BigRat(q.str().c_str()));
                        require(ray.gate.coordinate().expression()==(1-f)*start.expression()+f*parent.expression(),
                            "reflected gate uses wrong sign or a rounded parent");
                        require((ray.gate.coordinate().expression()-start.expression()).sign()
                            ==-(parent.expression()-start.expression()).sign(),"reflected gate is not outside actual extreme");
                    }
                }
                require(incident==std::set<std::string>(n.incidentArcs().begin(),n.incidentArcs().end()),
                    "partial local coverage counted as a complete native star");
                if(extreme){++extremeStars;require(reflected==1&&incident.size()==4,
                    "D extreme omitted its actual unbounded outward incidence");}
                else require(reflected==0,"reflection changed a nonextreme neighborhood");
                ++passed;out<<"NEIGHBORHOOD_BEGIN\n"<<n.canonicalRecord()<<"NEIGHBORHOOD_END\n";
                std::cerr<<"FINITE_ENDPOINT_PASSED "<<id<<" FACE "<<face<<'\n';
            }catch(const RankTwoSourceError& e){
                // Keep every unsupported star, not only the first failed event.
                ++unknown;out<<"UNKNOWN "<<e.code()<<' '<<e.what()<<'\n';
                std::cerr<<"FINITE_ENDPOINT_UNKNOWN "<<id<<" FACE "<<face<<' '<<e.code()<<' '<<e.what()<<'\n';
            }
            out.flush();require(bool(out),"truncated endpoint inventory");
        }
    }
    out<<"TOTAL_REQUIRED "<<required<<" PASSED "<<passed<<" UNKNOWN "<<unknown
       <<" EXTREME_STARS "<<extremeStars<<" REFLECTED_GATES "<<reflectedGates<<'\n';out.flush();
    require(required>16,"inventory did not include non-F endpoints");
    require(unknown==0,"MISSING_FEATURE: actual finite boundary contains unsupported native endpoint stars");
    require(required==40&&extremeStars==8&&reflectedGates==8,"frozen whole finite endpoint census changed");
    const auto& a=wholeArrangementRecords(whole).normalFaces.front();
    const auto& vertices=arrangementRecords(a).vertices;
    const auto endpoint=vertices.front().identity,anchor=vertices.back().identity;
    const auto rejects=[&](auto f,const char* code){bool caught=false;try{f();}
        catch(const RankTwoSourceError& e){caught=e.code()==code;}require(caught,"invalid reflected gate accepted");};
    rejects([&]{RankTwoMovingGate::fromReflectedVertices(a,endpoint,anchor,0,Q(0));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{RankTwoMovingGate::fromReflectedVertices(a,endpoint,anchor,0,Q(1));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{RankTwoMovingGate::fromReflectedVertices(a,endpoint,anchor,0,Q(-1,2));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{RankTwoMovingGate::fromReflectedVertices(a,"foreign",anchor,0,Q(1,2));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{RankTwoMovingGate::fromVertices(a,endpoint,anchor,0,Q(-1,2));},"INPUT_IDENTITY_MISMATCH");
    out<<"reflected_gate_validation PASSED\nfinite_registry_import NOT_GRANTED\nwhole_carrier NOT_GRANTED\n";out.flush();
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");
    tests(argv[1]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
