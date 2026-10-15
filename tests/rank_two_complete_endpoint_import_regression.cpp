#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
#define RANK_TWO_VERIFY_COMPLETE_FINITE_MESH 1
#include "spin/detail/RankTwoInventoryMaterialMesh.h"
#include "spin/detail/RankTwoWholeInteriorInventory.h"
#endif
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_MESH
#define RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE 1
#include "spin/detail/RankTwoCoveredFiniteMesh.h"
#include <cstring>
#endif
#include "spin/detail/RankTwoEndpointRegistry.h"
#include "spin/detail/RankTwoFinitePartition.h"
#include "spin/detail/RankTwoRegistryContinuation.h"
#ifdef RANK_TWO_VERIFY_COVERAGE_TIMING
#include "rank_two_phase_profile.h"
#include "spin/detail/RankTwoDiagnostics.h"
#define PROFILE_STAGE(name) profile.stage(name)
#else
#define PROFILE_STAGE(name) do{}while(false)
#endif
#ifdef RANK_TWO_VERIFY_ENDPOINT_CONSUMPTION
#include "spin/detail/RankTwoMixedArcIntervals.h"
#endif
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE
#include "spin/detail/RankTwoFiniteCoverageContinuation.h"
#include "spin/detail/RankTwoDiagnostics.h"
#endif
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
spin::PreparedRankTwoNormalMesh octa(){return spin::PreparedRankTwoNormalMesh::fromExact(
    {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
    {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}});}
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_MESH
spin::RankTwoCoveredFiniteMesh verifyCoveredMesh(const spin::RankTwoArrangement& a,const spin::RankTwoFiniteArcCoverage& coverage,
    std::size_t face,std::ostream& out){
    using namespace spin;
    std::cerr<<"COMPLETE_FINITE_MESH_BEGIN FACE "<<face<<'\n';
    const auto mesh=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    const auto& p=coverage.partition();const auto& r=p.registry();const auto& t=mesh.triangulation();
    require(mesh.coverage().scientificHash()==coverage.scientificHash()
        &&t.partitionIdentity()==p.scientificHash()&&t.registryIdentity()==r.scientificHash(),
        "complete mesh consumed stale coverage or an independently rebuilt owner registry");
    const auto n=p.pointIndices().size(),m=p.edgeIndices().size();
    require(mesh.incidence().partitionIdentity()==p.scientificHash()
        &&mesh.rotation().partitionIdentity()==p.scientificHash()
        &&mesh.incidence().pointPairCount()==n*(n-1)/2
        &&mesh.incidence().pointEdgeCount()==m*(n-2)
        &&mesh.incidence().edgePairCount()==m*(m-1)/2
        &&mesh.rotation().pointCount()==n&&mesh.rotation().finiteDirectionCount()==2*m,
        "complete finite mesh silently omitted native partition predicates or consumed excluded coarse arcs");
    using Edge=std::pair<std::size_t,std::size_t>;
    std::set<std::size_t> used,regions;
    std::map<std::size_t,std::map<Edge,int>> boundaries;
    for(const auto& triangle:t.triangles()){
        used.insert(triangle.points.begin(),triangle.points.end());regions.insert(triangle.registryRegionIndex);
        for(int j=0;j<3;++j){const auto x=triangle.points[j],y=triangle.points[(j+1)%3];
            boundaries[triangle.registryRegionIndex][std::minmax(x,y)]+=x<y?1:-1;}
    }
    require(used==std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end())
        &&regions==std::set<std::size_t>(p.regionIndices().begin(),p.regionIndices().end()),
        "complete mesh omitted a native point/cell or introduced an unowned one");
    std::set<std::size_t> summaries;
    for(const auto& summary:t.regionSummaries()){
        require(summaries.insert(summary.registryRegionIndex).second,"duplicate native region summary");
        const auto& cell=r.regions().at(summary.registryRegionIndex);
        require(summary.holes==cell.innerLoops.size()&&summary.connectedComponents==cell.outerLoops.size(),
            "complete mesh lost a native hole or component");
        std::map<Edge,int> wanted;
        const auto add=[&](const RankTwoRegistryLoop& loop){for(const auto& ref:loop.edges){
            const auto& e=r.edges().at(ref.edgeIndex);auto x=e.points[0],y=e.points[1];
            if(ref.reversed)std::swap(x,y);wanted[std::minmax(x,y)]+=x<y?1:-1;}};
        for(const auto& loop:cell.outerLoops)add(loop);for(const auto& loop:cell.innerLoops)add(loop);
        auto actual=boundaries.at(summary.registryRegionIndex);
        for(auto i=actual.begin();i!=actual.end();)if(!i->second)i=actual.erase(i);else ++i;
        for(auto i=wanted.begin();i!=wanted.end();)if(!i->second)i=wanted.erase(i);else ++i;
        require(actual==wanted,"complete triangles do not consume the final directed native boundary");
    }
    require(summaries==regions&&!t.triangles().empty(),"incomplete covered-mesh region census");
    // Save actual triangles and their final authority, without duplicating the
    // multi-GB coverage/registry record that was written immediately above.
    out<<"COMPLETE_FINITE_MESH_BEGIN FACE "<<face<<'\n'<<"MESH_AUTHORITY "
       <<a.scientificHash()<<' '<<coverage.scientificHash()<<' '<<mesh.scientificHash()<<' '
       <<p.scientificHash()<<' '<<r.scientificHash()<<' '<<mesh.incidence().scientificHash()<<' '
       <<mesh.rotation().scientificHash()<<' '<<t.scientificHash()<<'\n';
    out<<"MESH_PARTITION_SCOPE "<<mesh.incidence().partitionIdentity()<<' '
       <<mesh.rotation().partitionIdentity()<<' '<<n<<' '<<m<<' '
       <<mesh.incidence().pointPairCount()<<' '<<mesh.incidence().pointEdgeCount()<<' '
       <<mesh.incidence().edgePairCount()<<' '<<mesh.rotation().finiteDirectionCount()<<' '
       <<mesh.rotation().omittedInfiniteDirectionCount()<<' '
       <<mesh.rotation().omittedOutsidePartitionDirectionCount()<<'\n';
    for(const auto& triangle:t.triangles())out<<"TRIANGLE "<<triangle.registryRegionIndex<<' '
        <<triangle.points[0]<<' '<<triangle.points[1]<<' '<<triangle.points[2]<<'\n';
    for(const auto& s:t.regionSummaries())out<<"REGION_SUMMARY "<<s.registryRegionIndex<<' '
        <<s.triangleCount<<' '<<s.vertexCount<<' '<<s.edgeCount<<' '<<s.connectedComponents<<' '
        <<s.holes<<' '<<s.eulerCharacteristic<<'\n';
    out<<"COMPLETE_FINITE_MESH_END FACE "<<face<<'\n';out.flush();
    require(bool(out),"incomplete covered-mesh evidence write");
    std::cerr<<"COMPLETE_FINITE_MESH_PASS FACE "<<face<<" POINTS "<<used.size()
        <<" TRIANGLES "<<t.triangles().size()<<" CELLS "<<regions.size()<<'\n';
    return mesh;
}
#endif
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
struct MaterialSeamConsumer {
    std::size_t face;
    std::vector<std::string> points;
    std::vector<std::array<std::string,2>> roots;
    std::vector<CORE::Expr> lower,upper;
    std::vector<int> direction;
    bool active;
};
using MaterialSeamTable=std::map<std::string,std::vector<MaterialSeamConsumer>>;
void verifyOwnerMaterial(const spin::PreparedRankTwoSource& source,const spin::RankTwoArrangement& a,
    const spin::RankTwoCoveredFiniteMesh& base,const spin::RankTwoInteriorInventory& inventory,
    const spin::RankTwoSharedFiniteSeams& master,std::size_t face,MaterialSeamTable& consumers,std::ostream& out){
    using namespace spin;
    require(inventory.points().size()==6,"original native owner inventory was replaced by global/coarse census");
    std::cerr<<"COMPLETE_OWNER_MATERIAL_BEGIN FACE "<<face<<" INVENTORY "<<inventory.scientificHash()<<std::endl;
    // Observe only this new composition stage. The existing whole-coverage
    // observer intentionally filters these events; retain the ACTUAL failing
    // affine input/edge instead of rerunning a guessed synthetic replacement.
    detail::ScopedRankTwoDiagnostics materialDiagnostic([](const char* phase,void*)noexcept{
        if(std::strncmp(phase,"strict-material-mesh:",sizeof("strict-material-mesh:")-1)==0
            ||std::strncmp(phase,"material-lift:",sizeof("material-lift:")-1)==0
            ||std::strncmp(phase,"weighted-coordinate:",sizeof("weighted-coordinate:")-1)==0)
            std::cerr<<"OWNER_MATERIAL_DIAGNOSTIC "<<phase<<'\n';
    },nullptr);
    const auto result=buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,inventory);
    const auto& initial=result.interiorMesh();const auto& strict=result.strictMaterial();const auto& mesh=strict.mesh();
    const auto& registry=base.coverage().partition().registry();const auto& lifts=strict.lifts();
    require(result.interiors().scientificHash()==inventory.scientificHash()&&result.baseMeshIdentity()==base.scientificHash(),
        "original immutable owner inventory or covered mesh replaced");
    require(initial.points().size()==registry.points().size()+6&&initial.triangles().size()==base.triangulation().triangles().size()+12,
        "all original six owner points not inserted in original exact mesh");
    const auto repairs=strict.insertedEdgePoints().size();
    require(mesh.points().size()==initial.points().size()+repairs&&mesh.triangles().size()==initial.triangles().size()+2*repairs,
        "incomplete all-edge material repair");
    require(lifts.allActiveStrictPassed(),"original owner active material failed strict thickness");
    std::set<std::string> ids;for(const auto& p:mesh.points())ids.insert(p.identity);
    for(const auto& p:inventory.points())require(ids.count(p.scientificHash()),"original inventory point lost after thickness repair");
    using Pair=std::pair<std::size_t,std::size_t>;
    std::map<Pair,std::vector<std::pair<std::size_t,int>>> adjacent;
    for(const auto& t:mesh.triangles())for(int j=0;j<3;++j){const auto u=t.points[j],v=t.points[(j+1)%3];
        adjacent[std::minmax(u,v)].push_back({t.registryRegionIndex,u<v?1:-1});}
    for(const auto index:base.coverage().partition().edgeIndices()){
        const auto& e=registry.edges().at(index);
        require(adjacent.count(std::minmax(e.points[0],e.points[1])),"original native constraint disappeared from actual material");
    }
    for(const auto& seam:master.intervals())if(seam.finite)for(int side=0;side<2;++side)if(seam.incidentFaces[side]==face){
        const auto arc=std::find_if(registry.arcs().begin(),registry.arcs().end(),[&](const auto& v){return v.identity==seam.sourceArcs[side];});
        require(arc!=registry.arcs().end()&&!arc->edges.empty(),"material omitted actual shared native arc");
        MaterialSeamConsumer c;c.face=face;c.active=false;
        std::vector<std::size_t> indices;std::optional<std::size_t> activeRegion;
        for(const auto edge:arc->edges){const auto& e=registry.edges().at(edge);
            if(indices.empty())indices.push_back(e.points[0]);
            require(indices.back()==e.points[0],"shared native edges not one literal ordered chain");indices.push_back(e.points[1]);
            const auto found=adjacent.find(std::minmax(e.points[0],e.points[1]));
            require(found!=adjacent.end(),"shared constraint was independently resampled during material construction");
            int directed=0;std::size_t activeCount=0;
            for(const auto& q:found->second){const auto kind=registry.regions().at(q.first).kind;
                if(kind==RankTwoRegionKind::Full||kind==RankTwoRegionKind::Cut){++activeCount;
                    require(!activeRegion||*activeRegion==q.first,"shared interval crossed original active native cell");activeRegion=q.first;
                    directed=q.second*(e.points[0]<e.points[1]?1:-1);}}
            require(activeCount<=1,"shared normal constraint is internal to more than one local active owner");
            c.direction.push_back(directed*(seam.sourceReversed[side]?-1:1));
        }
        if(seam.sourceReversed[side]){std::reverse(indices.begin(),indices.end());std::reverse(c.direction.begin(),c.direction.end());}
        c.active=bool(activeRegion);
        for(const auto i:indices){const auto& p=mesh.points().at(i);const auto& old=registry.points().at(i);
            require(p.identity==old.identity,"material changed literal common point identity");c.points.push_back(p.identity);
            c.roots.push_back({p.coordinates[0].certificate().rootIdentity,p.coordinates[1].certificate().rootIdentity});
            if(c.active){c.lower.push_back(lifts.lower(*activeRegion,i));c.upper.push_back(lifts.upper(i));}}
        if(c.active)for(std::size_t i=1;i<indices.size();++i)
            require((c.upper[i-1]-c.lower[i-1]+c.upper[i]-c.lower[i]).sign()>0,
                "local normal-boundary exemption hid zero thickness at a globally INTERNAL shared edge");
        out<<"MATERIAL_SHARED_CONSUMER "<<seam.identity<<' '<<face<<' '<<side<<' '<<c.active<<' '<<c.points.size()<<'\n';
        for(std::size_t i=0;i<c.points.size();++i)out<<"MATERIAL_SHARED_POINT "<<i<<' '<<c.points[i]<<' '<<c.roots[i][0]<<' '<<c.roots[i][1]<<'\n';
        for(std::size_t i=0;i<c.direction.size();++i)out<<"MATERIAL_SHARED_DIRECTION "<<i<<' '<<c.direction[i]<<'\n';
        consumers[seam.identity].push_back(std::move(c));
    }
    out<<"COMPLETE_OWNER_MATERIAL_BEGIN FACE "<<face<<" AUTOMATIC "<<result.scientificHash()<<'\n'
       <<result.canonicalRecord()<<"COMPLETE_OWNER_MATERIAL_END FACE "<<face<<'\n';out.flush();require(bool(out),"incomplete original owner material evidence");
    std::cerr<<"COMPLETE_OWNER_MATERIAL_PASS FACE "<<face<<" POINTS "<<mesh.points().size()
        <<" TRIANGLES "<<mesh.triangles().size()<<" REPAIRS "<<repairs<<std::endl;
}
void verifyMaterialSeams(const spin::RankTwoSharedFiniteSeams& master,const MaterialSeamTable& consumers,std::ostream& out){
    std::size_t finite=0,active=0;
    for(const auto& seam:master.intervals())if(seam.finite){++finite;
        const auto found=consumers.find(seam.identity);require(found!=consumers.end()&&found->second.size()==2,"shared material interval lacks BOTH actual consumers");
        const auto& a=found->second[0];const auto& b=found->second[1];
        require(a.face!=b.face&&std::set<std::size_t>{a.face,b.face}==std::set<std::size_t>{seam.incidentFaces[0],seam.incidentFaces[1]},"wrong incident material consumers");
        require(a.points==b.points&&a.roots==b.roots&&a.direction.size()==b.direction.size()&&a.active==b.active,"shared material nodes/roots/labels disagree");
        if(a.active){++active;for(std::size_t i=0;i<a.points.size();++i)
            require(a.lower[i]==b.lower[i]&&a.upper[i]==b.upper[i],"same shared source point has different material lifts");}
        for(std::size_t i=0;i<a.direction.size();++i)require(a.direction[i]==-b.direction[i],"shared actual material boundary directions are not opposite");
    }
    require(finite==16&&consumers.size()==finite&&active==12,"shared material finite/active/empty census changed");
    out<<"COMPLETE_SHARED_MATERIAL_COMPATIBILITY 16 ACTIVE 12 EMPTY 4\n";
}
#endif
void test(const char* path){using namespace spin;
#ifdef RANK_TWO_VERIFY_COVERAGE_TIMING
    rank_two_test::PhaseProfile profile(std::cerr);
    detail::ScopedRankTwoDiagnostics sourceDiagnostic([](const char* phase,void* p)noexcept{
        static_cast<rank_two_test::PhaseProfile*>(p)->event(phase);
    },&profile);
#endif
    PROFILE_STAGE("harness:source-prepare");
    // A subset of separately passing F/D stars is not a whole import. This
    // test demands EVERY original endpoint of actual finite boundary edges.
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    PROFILE_STAGE("harness:arrangement");
    const auto whole=analyzeRankTwoArrangement(source,Q(3,2),octa());
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
    const auto inventory=enumerateRankTwoWholeMovingInteriors(source,whole);
    require(inventory.owners().size()==8&&inventory.pointCount()==24,"incomplete original normal-owner inventory");
    MaterialSeamTable materialConsumers;
#endif
    PROFILE_STAGE("harness:master-seams");
    const auto master=buildRankTwoSharedFiniteSeams(source,whole);
    PROFILE_STAGE("harness:endpoint-neighborhoods");
    std::ofstream out(path);require(bool(out),"cannot save complete endpoint import");
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
    out<<"rank-two-complete-owner-material-mesh-v1\n";
    out<<"WHOLE_INVENTORY_BEGIN\n"<<inventory.canonicalRecord()<<"WHOLE_INVENTORY_END\n";
#elif defined(RANK_TWO_VERIFY_COMPLETE_FINITE_MESH)
    out<<"rank-two-complete-finite-mesh-v1\n";
#elif defined(RANK_TWO_VERIFY_ENDPOINT_CONSUMPTION)
    out<<"rank-two-complete-endpoint-binding-v1\n";
#elif defined(RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE)
    out<<"rank-two-complete-finite-coverage-v1\n";
#else
    out<<"rank-two-complete-endpoint-import-v1\n";
#endif
    out<<"MASTER_BEGIN\n"<<master.canonicalRecord()<<"MASTER_END\n";
    std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>> neighborhoods;
    std::vector<RankTwoArcSamplingPlan> plans;
    std::vector<std::set<std::string>> required;
    for(std::size_t face=0;face<4;++face){const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        RankTwoArcSamplingPlan plan;
        for(const auto& arc:arrangementRecords(a).arcs)
            if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
                plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
        const auto base=buildRankTwoArcRegistryWithSharedSeams(a,plan,master,face);const auto domain=buildRankTwoFinitePartition(a,base);
        std::set<std::string> wanted;
        for(const auto e:domain.edgeIndices()){
            const auto arc=std::find_if(arrangementRecords(a).arcs.begin(),arrangementRecords(a).arcs.end(),
                [&](const auto& p){return p.identity==base.edges().at(e).arcIdentity;});
            require(arc!=arrangementRecords(a).arcs.end(),"finite edge lost actual native arc");
            for(const auto& end:arc->ends)if(end.vertexIdentity)wanted.insert(*end.vertexIdentity);
        }
        require(wanted.size()==10,"original complete finite endpoint census changed");
        for(const auto& event:wanted){const auto n=certifyRankTwoEndpointNeighborhood(a,event);
            neighborhoods.emplace_back(face,n);out<<"ENDPOINT "<<event<<" FACE "<<face<<" NEIGHBORHOOD_BEGIN\n"
                <<n.canonicalRecord()<<"NEIGHBORHOOD_END\n";out.flush();
            const auto vertex=std::find_if(arrangementRecords(a).vertices.begin(),arrangementRecords(a).vertices.end(),
                [&](const auto& v){return v.identity==event;});
            if(detail::CertifiedAlgebraicReal::reload(vertex->coordinates[0]).expression()==0
                &&detail::CertifiedAlgebraicReal::reload(vertex->coordinates[1]).expression()==0){
                require(n.fans().size()==2,"normal vertex lost its two actual graph sectors");
                for(const auto& fan:n.fans()){
                    require(fan.orderedArcs().size()==1,"normal-vertex graph sector census changed");
                    const auto arc=std::find_if(arrangementRecords(a).arcs.begin(),arrangementRecords(a).arcs.end(),
                        [&](const auto& p){return p.identity==fan.orderedArcs().front().arcIdentity;});
                    require(arc!=arrangementRecords(a).arcs.end(),"normal-vertex fan lost actual native arc");
                    const auto other=arc->ends[0].vertexIdentity==event?arc->ends[1].vertexIdentity:arc->ends[0].vertexIdentity;
                    require(bool(other),"normal-vertex shared arc lacks other original finite endpoint");
                    std::cerr<<"NORMAL_VERTEX_ANCHOR FACE "<<face<<" ACTUAL "<<fan.gate().anchorIdentity()<<" EXPECTED "<<*other<<'\n';
                    require(fan.gate().anchorIdentity()==*other,
                        "normal-vertex fan ruler uses unrelated off-seam native event instead of actual adjacent shared endpoint");
                }
            }
        }
        required.push_back(wanted);plans.push_back(plan);
    }
    require(neighborhoods.size()==40,"some required native endpoint omitted before publication");
    std::cerr<<"COMPLETE_PUBLICATION_BEGIN forty actual required native stars\n";
    PROFILE_STAGE("harness:endpoint-publication");
    const auto publication=publishRankTwoSharedEndpointSeams(source,master,neighborhoods);
    out<<"PUBLICATION_BEGIN\n"<<publication.canonicalRecord()<<"PUBLICATION_END\n";out.flush();
    std::size_t used=0,outside=0,localPrimitiveAliases=0;
    for(std::size_t face=0;face<4;++face){const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        PROFILE_STAGE("harness:import-and-validate");
        std::cerr<<"COMPLETE_IMPORT_BEGIN FACE "<<face<<'\n';
        const auto registry=buildRankTwoArcRegistryWithSharedEndpointSeams(a,plans[face],publication,face);
        const auto domain=buildRankTwoFinitePartition(a,registry);
        // A local cone-excluded incidence still belongs to a complete star.
        // Its fan may name exactly an EXISTING native primitive dependency.
        // Preserve that dependency instead of inserting a second zero chord.
        const auto original=buildRankTwoArcRegistryWithSharedSeams(a,plans[face],master,face);
        for(const auto& item:neighborhoods)if(item.first==face)for(const auto& fan:item.second.fans()){
            const auto& g=fan.gate();
            for(const auto& proof:fan.orderedArcs()){
                if(std::any_of(publication.bindings().begin(),publication.bindings().end(),
                    [&](const auto& b){return b.face==face&&b.proofIdentity==proof.identity;}))continue;
                const auto old=std::find_if(original.points().begin(),original.points().end(),[&](const auto& p){
                    return p.kind==RankTwoRegistryPointKind::ArcSample&&p.arcIdentity==proof.arcIdentity
                        &&((p.endpointIdentities[0]==g.endpointIdentity()&&p.endpointIdentities[1]==g.anchorIdentity()
                            &&p.fraction==g.fraction())
                           ||(p.endpointIdentities[1]==g.endpointIdentity()&&p.endpointIdentities[0]==g.anchorIdentity()
                            &&p.fraction==Q(1)-g.fraction()));});
                if(old==original.points().end())continue;
                const auto actual=std::find_if(registry.points().begin(),registry.points().end(),
                    [&](const auto& p){return p.identity==old->identity;});
                require(actual!=registry.points().end()&&actual->kind==RankTwoRegistryPointKind::ArcSample
                    &&actual->endpointIdentities==old->endpointIdentities&&actual->fraction==old->fraction,
                    "local fan replaced the original literal primitive dependency");
                require(std::any_of(actual->additionalEndpointBindings.begin(),actual->additionalEndpointBindings.end(),
                    [&](const auto& b){return b.neighborhoodIdentity==item.second.scientificHash()
                        &&b.proofIdentity==proof.identity&&b.gateIdentity==g.identity()&&b.fanIdentity==fan.scientificHash();}),
                    "local primitive alias omitted its complete star proof");
                for(int d=0;d<2;++d)require(actual->coordinates[d].certificate().rootIdentity==old->coordinates[d].certificate().rootIdentity,
                    "local primitive alias independently resampled original native root");
                ++localPrimitiveAliases;
            }
        }
        require(registry.endpointNeighborhoods().size()==10,"common finite registry lost required complete star");
        require(registry.arcs().size()==arrangementRecords(a).arcs.size(),"whole import omitted native arc");
        std::vector<RankTwoEndpointNeighborhood> local;std::set<std::string> actual;
        for(const auto& n:registry.endpointNeighborhoods()){
            require(actual.insert(n.endpointIdentity()).second,"duplicated complete star");local.push_back(n);++used;
            std::set<std::string> covered;
            for(const auto& e:registry.edges())if(e.endpointNeighborhoodIdentity==n.scientificHash())covered.insert(e.arcIdentity);
            for(const auto& x:registry.outsideIncidences())if(x.neighborhoodIdentity==n.scientificHash()){
                require(covered.insert(x.arcIdentity).second,"outside incidence also has finite edge");++outside;}
            require(covered==std::set<std::string>(n.incidentArcs().begin(),n.incidentArcs().end()),"whole star authority is only partial");
        }
        require(actual==required[face],"whole import changed original endpoint set");
        const auto again=rebuildRankTwoArcRegistryWithSubdivisions(a,registry,plans[face],local,{},{});
        require(again.canonicalRecord()==registry.canonicalRecord(),"whole import continuation changed shared authority");
        for(const auto& b:publication.bindings())if(b.face==face){const auto& p=publication.points().at(b.pointIndex);
            const auto actualPoint=std::find_if(registry.points().begin(),registry.points().end(),[&](const auto& q){return q.identity==p.identity;});
            require(actualPoint!=registry.points().end(),"owner omitted or independently resampled shared endpoint gate");
            for(int d=0;d<2;++d)require(actualPoint->coordinates[d].certificate().rootIdentity==p.coordinates[d].certificate().rootIdentity,
                "owner rounded shared common root independently");}
        out<<"IMPORTED_FACE "<<face<<" REGISTRY_BEGIN\n"<<registry.canonicalRecord()<<"REGISTRY_END\nPARTITION_BEGIN\n"
           <<domain.canonicalRecord()<<"PARTITION_END\n";out.flush();
#ifdef RANK_TWO_VERIFY_ENDPOINT_CONSUMPTION
        // Consume ALL actual finite-domain endpoint chords through the real
        // mixed refiner. Primary metadata alone loses the second master proof
        // and every local primitive alias; an import-only test cannot catch it.
        std::vector<RankTwoMixedArcInterval> intervals;std::set<std::string> wantedEdges;
        for(const auto k:domain.edgeIndices()){
            const auto& e=registry.edges().at(k);if(e.endpointNeighborhoodIdentity.empty())continue;
            intervals.push_back({e.arcIdentity,{registry.points().at(e.points[0]).identity,
                registry.points().at(e.points[1]).identity}});wantedEdges.insert(e.identity);
        }
        require(!wantedEdges.empty(),"fixture has no actual endpoint chords to consume");
        std::cerr<<"COMPLETE_ENDPOINT_CONSUMPTION_BEGIN FACE "<<face<<" EDGES "<<wantedEdges.size()<<'\n';
        PROFILE_STAGE("harness:endpoint-consumption");
#ifdef RANK_TWO_VERIFY_COVERAGE_TIMING
        detail::ScopedRankTwoDiagnostics diagnostic([](const char* phase,void* p)noexcept{
            static_cast<rank_two_test::PhaseProfile*>(p)->event(phase);
        },&profile);
#endif
        const auto consumed=refineRankTwoMixedArcIntervalsFromRegistry(a,registry,plans[face],local,intervals,Q(1,1048576));
        PROFILE_STAGE("harness:validate-consumed-output");
        require(consumed.registry().canonicalRecord()==registry.canonicalRecord(),"endpoint consumption changed sealed registry");
        require(std::set<std::string>(consumed.endpointEdges().begin(),consumed.endpointEdges().end())==wantedEdges
            &&consumed.proofs().empty()&&consumed.rounds().size()==1&&consumed.rounds().front().splitEdges==0,
            "mixed refiner did not consume each actual sealed endpoint chord exactly once");
        out<<"ENDPOINT_CONSUMPTION_BEGIN\n"<<consumed.canonicalRecord()<<"ENDPOINT_CONSUMPTION_END\n";out.flush();
#endif
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE
        // A refiner that proves only endpoint chords or selected successful
        // regular intervals must fail this test. The real whole factory
        // chooses ALL actual finite-domain arcs from the same sealed authority.
        std::cerr<<"COMPLETE_FINITE_COVERAGE_BEGIN FACE "<<face<<" EDGES "<<domain.edgeIndices().size()<<'\n';
        PROFILE_STAGE("harness:whole-coverage");
#ifdef RANK_TWO_VERIFY_COVERAGE_TIMING
        detail::ScopedRankTwoDiagnostics diagnostic([](const char* phase,void* p)noexcept{
            static_cast<rank_two_test::PhaseProfile*>(p)->event(phase);
        },&profile);
#else
        detail::ScopedRankTwoDiagnostics diagnostic([](const char* phase,void*)noexcept{
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_MESH
            // Observation only: do not serialize tens of millions of root
            // events. Preserve product predicates and their original budgets.
            if(std::strncmp(phase,"covered-finite-mesh:",sizeof("covered-finite-mesh:")-1)!=0
                &&std::strncmp(phase,"chord-sweep:",sizeof("chord-sweep:")-1)!=0
                &&std::strncmp(phase,"finite-cdt:insert-",sizeof("finite-cdt:insert-")-1)!=0
                &&std::strncmp(phase,"finite-cdt:constraint-",sizeof("finite-cdt:constraint-")-1)!=0)return;
#endif
            std::cerr<<"WHOLE_COVERAGE_DIAGNOSTIC "<<phase<<'\n';
        },nullptr);
#endif
        const auto covered=refineRankTwoFiniteArcCoverageFromRegistry(a,registry,plans[face],local,Q(1,1048576));
        PROFILE_STAGE("harness:validate-coverage-output");
        const auto& final=covered.partition().registry();const auto& refinement=covered.refinement();
        require(refinement.rounds().front().registryIdentity==registry.scientificHash(),
            "whole coverage did not consume original complete common registry");
        require(final.endpointNeighborhoods().size()==10&&final.outsideIncidences().size()==registry.outsideIncidences().size(),
            "whole refinement discarded original stars or certified outside incidences");
        std::set<std::string> actualEdges,provedEdges;
        for(const auto i:covered.partition().edgeIndices())actualEdges.insert(final.edges().at(i).identity);
        for(const auto& proof:refinement.proofs()){
            require(proof.registryIdentity()==final.scientificHash()&&proof.strip().chordDistanceUpperBound<=Q(1,1048576),
                "ordinary whole witness is obsolete or exceeds locked precision");
            require(provedEdges.insert(proof.edgeIdentity()).second,"duplicate ordinary whole witness");
        }
        for(const auto& id:refinement.endpointEdges())require(provedEdges.insert(id).second,"overlapping whole endpoint witness");
        require(!refinement.proofs().empty()&&!refinement.endpointEdges().empty()&&actualEdges==provedEdges,
            "whole coverage omitted ordinary or endpoint finite-domain chords");
        require(refinement.rounds().back().splitEdges==0,"incomplete whole refinement was accepted");
        for(const auto& old:registry.points()){
            const auto p=std::find_if(final.points().begin(),final.points().end(),[&](const auto& q){return q.identity==old.identity;});
            require(p!=final.points().end(),"whole refinement lost an original common dependency");
            for(int d=0;d<2;++d)require(p->coordinates[d].certificate().rootIdentity==old.coordinates[d].certificate().rootIdentity,
                "whole refinement independently resampled an original common root");
        }
        out<<"WHOLE_FINITE_COVERAGE_BEGIN\n"<<covered.canonicalRecord()<<"WHOLE_FINITE_COVERAGE_END\n";out.flush();
        std::cerr<<"COMPLETE_FINITE_COVERAGE_PASS FACE "<<face<<" EDGES "<<actualEdges.size()
            <<" ORDINARY "<<refinement.proofs().size()<<" ENDPOINT "<<refinement.endpointEdges().size()<<'\n';
#ifdef RANK_TWO_VERIFY_COMPLETE_FINITE_MESH
        const auto actualMesh=verifyCoveredMesh(a,covered,face,out);
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
        verifyOwnerMaterial(source,a,actualMesh,inventory.owners().at(face),master,face,materialConsumers,out);
#endif
#endif
#endif
    }
    require(used==40&&outside==8&&localPrimitiveAliases==4,"whole endpoint or literal local alias census changed during actual import");
#ifdef RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH
    verifyMaterialSeams(master,materialConsumers,out);
#endif
    out<<"COMPLETE_ENDPOINT_IMPORT 40 OUTSIDE_NATIVE_RAYS 8 LOCAL_PRIMITIVE_ALIASES 4\nregular_strip_coverage NOT_GRANTED\n"
       <<"global_isotopy NOT_GRANTED\nstrict_thickness NOT_GRANTED\n";
#ifdef RANK_TWO_VERIFY_COVERAGE_TIMING
    // A separate diagnostic method replay, AFTER saved endpoint science.
    // Freeze first finite nonvertical native arc in each original owner and
    // its midpoint before executing. This is not ordinary full coverage.
    PROFILE_STAGE("harness:source-projection-replay");
    for(std::size_t face=0;face<4;++face){
        const auto& a=wholeArrangementRecords(whole).normalFaces.at(face);
        const auto& arcs=arrangementRecords(a).arcs;
        const auto arc=std::find_if(arcs.begin(),arcs.end(),[](const auto& v){
            return !v.vertical&&v.ends[0].location==RankTwoArcEndLocation::Finite
                &&v.ends[1].location==RankTwoArcEndLocation::Finite;
        });
        require(arc!=arcs.end(),"actual source has no finite nonvertical projection replay");
        const auto s=sampleRankTwoArc(a,arc->identity,Q(1,2));
        const auto projected=detail::evaluateRankTwoArcProjection(a,arc->identity,s.coordinates[0]);
        require(detail::rankTwoAlgebraicHandle(projected[0]).compare(detail::rankTwoAlgebraicHandle(s.coordinates[0]))==CGAL::EQUAL
            &&detail::rankTwoAlgebraicHandle(projected[1]).compare(detail::rankTwoAlgebraicHandle(s.coordinates[1]))==CGAL::EQUAL,
            "actual source projection and sampled native fiber disagree");
    }
    // Diagnostic-only actual common-registry conversion. Select by immutable
    // input degree, not runtime or success. Keep every original seam/star;
    // replay eight left children of ONE ordinary domain edge, not full error
    // coverage. The product/geometry budgets and arithmetic stay unchanged.
    PROFILE_STAGE("harness:source-subdivision-replay");
    const auto& replayArrangement=wholeArrangementRecords(whole).normalFaces.at(0);
    auto replay=buildRankTwoArcRegistryWithSharedEndpointSeams(replayArrangement,plans[0],publication,0);
    const auto replayDomain=buildRankTwoFinitePartition(replayArrangement,replay);
    std::set<std::string> publishedSeamArcs;
    for(const auto& seam:master.intervals())if(seam.finite)
        for(int side=0;side<2;++side)if(seam.incidentFaces[side]==0)
            publishedSeamArcs.insert(seam.sourceArcs[side]);
    const RankTwoRegistryEdge* chosen=nullptr;std::size_t chosenProduct=0;
    for(const auto index:replayDomain.edgeIndices()){
        const auto& edge=replay.edges().at(index);
        if(!edge.endpointNeighborhoodIdentity.empty()||!edge.endpointFanIdentity.empty())continue;
        if(publishedSeamArcs.count(edge.arcIdentity))continue;
        const auto& arcs=arrangementRecords(replayArrangement).arcs;
        const auto arc=std::find_if(arcs.begin(),arcs.end(),[&](const auto& v){return v.identity==edge.arcIdentity;});
        require(arc!=arcs.end(),"diagnostic subdivision has no native arc");
        const int axis=arc->vertical?1:0;
        const auto degree=[&](std::size_t point){return replay.points().at(point).coordinates[axis].certificate().coefficients.size()-1;};
        const auto product=degree(edge.points[0])*degree(edge.points[1]);
        if(!chosen||product>chosenProduct||(product==chosenProduct&&edge.identity<chosen->identity)){
            chosen=&edge;chosenProduct=product;
        }
    }
    require(chosen,"source diagnostic has no ordinary finite domain edge");
    const auto arcId=chosen->arcIdentity,first=replay.points().at(chosen->points[0]).identity;
    auto right=replay.points().at(chosen->points[1]).identity;
    std::cerr<<"SOURCE_SUBDIVISION_REPLAY_BEGIN FACE 0 EDGE "<<chosen->identity
        <<" PARENT_DEGREE_PRODUCT "<<chosenProduct<<" LEFT_CHILDREN 8 NOT_FULL_COVERAGE\n";
    std::vector<RankTwoEndpointNeighborhood> localReplay;
    for(const auto& n:neighborhoods)if(n.first==0)localReplay.push_back(n.second);
    auto subdivisions=replay.subdivisions();
    for(int depth=0;depth<8;++depth){
        const std::array<std::string,2> parents{first,right};
        subdivisions.push_back({arcId,parents,Q(1,2)});
        auto next=rebuildRankTwoArcRegistryWithSubdivisions(replayArrangement,replay,plans[0],localReplay,subdivisions);
        const auto point=std::find_if(next.points().begin(),next.points().end(),[&](const auto& v){
            return v.kind==RankTwoRegistryPointKind::ArcSubdivision&&v.arcIdentity==arcId
                &&v.endpointIdentities==parents&&v.fraction==Q(1,2);
        });
        require(point!=next.points().end(),"actual diagnostic conversion omitted requested moving point");
        right=point->identity;replay=std::move(next);
        std::cerr<<"SOURCE_SUBDIVISION_REPLAY_STEP "<<depth<<" REGISTRY "<<replay.scientificHash()<<'\n';
    }
    profile.finish();
#endif
}
}
int main(int argc,char** argv){try{if(argc!=2)throw std::invalid_argument("evidence path required");test(argv[1]);return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
