#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoStrictMaterialMesh.h")
#include "spin/detail/RankTwoStrictMaterialMesh.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
using Pair=std::pair<std::size_t,std::size_t>;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong failure class");return;}throw std::runtime_error("invalid strict mesh accepted");}
spin::RankTwoInteriorPointRequest request(const spin::RankTwoArrangement& a,spin::RankTwoRegionKind kind,int sign){
    using namespace spin;const auto& n=arrangementRecords(a);RankTwoInteriorPointRequest r;
    const RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:n.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0&&(!left||x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0&&(!right||x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left&&right,"native adjacent events missing");r.eventIdentities={left->identity,right->identity};
    for(const auto& c:n.regions)if(c.kind==kind&&c.originalSigns[1]==sign){require(r.regionIdentity.empty(),"ambiguous region");r.regionIdentity=c.identity;}
    const auto x=[&](const std::string& id){for(const auto& v:n.vertices)if(v.identity==id)return C::reload(v.coordinates[0]).expression();throw std::runtime_error("native event absent");};
    for(const auto& arc:n.arcs){if(arc.vertical||arc.ends[0].location!=RankTwoArcEndLocation::Finite||arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        if(!(x(*arc.ends[0].vertexIdentity)<0&&x(*arc.ends[1].vertexIdentity)>0))continue;
        if(kind==RankTwoRegionKind::Cut&&C::reload(arc.interiorWitness[1]).expression().sign()<=0)continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous lower fiber");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous upper fiber");r.boundaryArcIdentities[1]=arc.identity;}}
    return r;
}
std::set<Pair> zeroEdges(const spin::RankTwoInteriorMesh& m,const spin::RankTwoMaterialLifts& l){std::set<Pair> out;
    for(const auto& w:l.triangles())for(int j=0;j<3;++j)if(w.zeroInternalEdges[j]){const auto& t=m.triangles().at(w.triangleIndex);out.insert(std::minmax(t.points[j],t.points[(j+1)%3]));}return out;}
std::size_t failed(const spin::RankTwoMaterialLifts& l){std::size_t n=0;for(const auto& w:l.triangles())n+=w.active&&!w.strictPassed;return n;}
void run(const char* path,bool annulus){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(annulus?rank_two_cases::annulus():RankTwoExactCoefficients{2,0,0,1,1,1,0,0,0,0},
        annulus?"finite-coverage-annulus-v1":"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,annulus?Q(3,2):Q(1));RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite&&arc.ends[1].location==RankTwoArcEndLocation::Finite)
        plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto coverage=refineRankTwoFiniteArcCoverage(a,plan,{},Q(1,32));const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    std::vector<RankTwoMovingInteriorPoint> inside;
    inside.push_back(buildRankTwoMovingInteriorPoint(a,request(a,RankTwoRegionKind::Full,-1)));
    if(annulus)for(int sign:{-1,1})inside.push_back(buildRankTwoMovingInteriorPoint(a,request(a,RankTwoRegionKind::Cut,sign)));
    const auto parent=insertRankTwoMovingInteriorNodes(a,base,inside);
    const auto before=buildRankTwoMaterialLifts(source,a,parent);const auto inventory=zeroEdges(parent,before);
    std::cout<<"INITIAL_CENSUS points="<<parent.points().size()<<" triangles="<<parent.triangles().size()
        <<" zero_edges="<<inventory.size()<<" failed="<<failed(before)<<" parent="<<parent.scientificHash()<<std::endl;
    // The historical annulus has 232 edges on its historical 2180-node
    // registry. Current unchanged sampling/error inputs construct 716 native
    // registry nodes; diagnostic10 measured ALL 272 edges and 292 failed faces.
    // Do not select the old 232-edge subset, or change the old regression.
    require(inventory.size()==(annulus?272u:21u)&&failed(before)==(annulus?292u:24u),
        "current frozen parent did not exercise the complete thickness inventory");
    std::size_t initialActive=0,initialEmpty=0,initialHoles=0;
    std::set<std::size_t> initialActiveOwners;
    for(const auto& w:before.triangles()){
        initialActive+=w.active;initialEmpty+=!w.active;
        if(w.active)initialActiveOwners.insert(w.registryRegionIndex);
    }
    for(const auto& r:parent.regionSummaries())initialHoles+=r.holes;
    require(initialActiveOwners.size()==(annulus?3u:1u)&&initialHoles==(annulus?3u:0u)
        &&(annulus?initialEmpty>0:initialEmpty==0),"native material owners, holes or Empty domain omitted");
    RankTwoStrictMaterialMeshSettings limited;limited.maximumNewEdgePoints=inventory.size()-1;
    rejects([&]{completeRankTwoMaterialThickness(source,a,parent,limited);},"RESOURCE_LIMIT");
    const auto foreign=PreparedRankTwoSource::fromExact(source.exactCoefficients(),"foreign-strict-frame");
    rejects([&]{completeRankTwoMaterialThickness(foreign,a,parent);},"INPUT_IDENTITY_MISMATCH");
    RankTwoStrictMaterialMeshSettings noPredicates;noPredicates.lifts.maximumPredicates=0;
    rejects([&]{completeRankTwoMaterialThickness(source,a,parent,noPredicates);},"RESOURCE_LIMIT");
    const auto result=completeRankTwoMaterialThickness(source,a,parent);
    const auto& mesh=result.mesh();const auto& lifts=result.lifts();
    require(result.parentMeshIdentity()==parent.scientificHash()&&lifts.meshIdentity()==mesh.scientificHash(),"strict grant detached from actual mesh");
    require(result.insertedEdgePoints().size()==inventory.size()&&mesh.points().size()==parent.points().size()+inventory.size()
        &&mesh.triangles().size()==parent.triangles().size()+2*inventory.size(),"factory did not consume every actual zero edge on BOTH sides");
    std::set<Pair> inserted;for(const auto& node:result.insertedEdgePoints()){
        require(node.parentMeshIdentity()==parent.scientificHash()&&node.definition().fraction==Q(1,2),"repair changed parent or frozen weight");
        inserted.insert(std::minmax(node.definition().pointIndices[0],node.definition().pointIndices[1]));}
    require(inserted==inventory&&zeroEdges(mesh,lifts).empty()&&lifts.allActiveStrictPassed()&&failed(lifts)==0,"partial repair promoted to strict material mesh");
    for(const auto& w:lifts.triangles())if(w.active){const auto& t=mesh.triangles().at(w.triangleIndex);std::array<E,3> lo,hi;
        for(int j=0;j<3;++j){lo[j]=lifts.lower(t.registryRegionIndex,t.points[j]);hi[j]=lifts.upper(t.points[j]);}
        detail::certifyRankTwoTriangleThickness(lo,hi,w.materialBoundary);}
    require(parent.regionSummaries().size()==mesh.regionSummaries().size(),"native region removed");
    std::size_t active=0,empty=0;for(const auto& w:lifts.triangles()){active+=w.active;empty+=!w.active;}
    require(active==initialActive+2*inventory.size()&&empty==initialEmpty,
        "actual material/Empty native census changed beyond both-sided edge subdivision");
    for(std::size_t i=0;i<mesh.regionSummaries().size();++i){const auto& p=parent.regionSummaries()[i];const auto& q=mesh.regionSummaries()[i];
        require(p.registryRegionIndex==q.registryRegionIndex&&p.connectedComponents==q.connectedComponents&&p.holes==q.holes
            &&p.eulerCharacteristic==q.eulerCharacteristic,"strict repair changed holes or owner topology");}
    // Complete input is a no-op, not another repair selected by stale parent indices.
    RankTwoStrictMaterialMeshSettings noNodes;noNodes.maximumNewEdgePoints=0;
    const auto again=completeRankTwoMaterialThickness(source,a,mesh,noNodes);
    require(again.insertedEdgePoints().empty()&&again.mesh().scientificHash()==mesh.scientificHash(),"already-strict mesh resampled or remeshed");
    require(parent.scientificHash()==result.parentMeshIdentity()&&!before.allActiveStrictPassed(),"original negative object overwritten");
    std::ofstream out(path);require(bool(out),"cannot save strict evidence");
    out<<"rank-two-strict-material-mesh-test-v1\nCASE "<<(annulus?"annulus":"full-disk")<<" INITIAL_ZERO_EDGES "<<inventory.size()
        <<" INITIAL_FAILED "<<failed(before)<<" FINAL_FAILED "<<failed(lifts)<<" ACTIVE "<<active<<" EMPTY "<<empty<<'\n'
        <<result.canonicalRecord()<<"CASE_END\n";out.flush();require(bool(out),"incomplete evidence");
    std::cout<<"PASS auto strict material mesh edges="<<inventory.size()<<" active="<<active<<" empty="<<empty<<"; native rings/negative budget/authority/no-op retained\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoStrictMaterialMesh.h")
    if(argc!=2&&!(argc==3&&std::string(argv[2])=="annulus"))throw std::invalid_argument("evidence path and optional annulus required");run(argv[1],argc==3);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: factory-owned complete zero-edge repair and strict material mesh grant");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
