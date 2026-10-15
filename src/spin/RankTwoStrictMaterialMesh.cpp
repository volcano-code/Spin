#include "detail/RankTwoStrictMaterialMesh.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Pair=std::pair<std::size_t,std::size_t>;
[[noreturn]]void resource(const char* why){throw RankTwoSourceError("RESOURCE_LIMIT",why);}
}
struct RankTwoStrictMaterialMesh::Impl {
    RankTwoInteriorMesh parent,finalMesh;
    RankTwoMaterialLifts initialLifts,finalLifts;
    std::vector<RankTwoMovingMeshEdgePoint> nodes;
    std::string record,hash;
    Impl(const RankTwoInteriorMesh& p,const RankTwoMaterialLifts& first,
        const RankTwoInteriorMesh& m,const RankTwoMaterialLifts& last,
        std::vector<RankTwoMovingMeshEdgePoint> n)
        :parent(p),finalMesh(m),initialLifts(first),finalLifts(last),nodes(std::move(n)){}
};
RankTwoStrictMaterialMesh::RankTwoStrictMaterialMesh(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoStrictMaterialMesh::parentMeshIdentity()const{return impl_->parent.scientificHash();}
const RankTwoInteriorMesh& RankTwoStrictMaterialMesh::mesh()const{return impl_->finalMesh;}
const RankTwoMaterialLifts& RankTwoStrictMaterialMesh::lifts()const{return impl_->finalLifts;}
const std::vector<RankTwoMovingMeshEdgePoint>& RankTwoStrictMaterialMesh::insertedEdgePoints()const{return impl_->nodes;}
const std::string& RankTwoStrictMaterialMesh::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoStrictMaterialMesh::scientificHash()const{return impl_->hash;}

RankTwoStrictMaterialMesh completeRankTwoMaterialThickness(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const RankTwoInteriorMesh& parent,
    const RankTwoStrictMaterialMeshSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(!settings.maximumInventoryChecks)resource("positive complete thickness inventory budget required");
    // These operations validate the original source/arrangement/mesh authority
    // and derive boundary exemptions from native adjacency, not caller flags.
    const auto first=buildRankTwoMaterialLifts(source,a,parent,settings.lifts);
    std::set<Pair> zeroEdges;
    std::uint64_t checks=0;
    const auto consume=[&]{if(checks==settings.maximumInventoryChecks)
        resource("complete thickness inventory budget exhausted");++checks;};
    std::size_t initialFailed=0;
    for(const auto& w:first.triangles()){
        consume();if(!w.active)continue;
        initialFailed+=!w.strictPassed;
        const auto& t=parent.triangles().at(w.triangleIndex);
        for(int j=0;j<3;++j){consume();if(w.zeroInternalEdges[j])
            zeroEdges.insert(std::minmax(t.points[j],t.points[(j+1)%3]));}
    }
    // Freeze ALL repairs before construction. An exhausted budget cannot
    // return a partially repaired object with the other failed cells omitted.
    if(zeroEdges.size()>settings.maximumNewEdgePoints)resource("all zero internal edges exceed registered new-point budget");
    std::vector<RankTwoMovingMeshEdgePoint> nodes;nodes.reserve(zeroEdges.size());
    for(const auto& edge:zeroEdges){
        RankTwoMovingMeshEdgeRequest request;request.pointIndices={edge.first,edge.second};
        request.fraction=CGAL::Gmpq(1,2);
        detail::rankTwoDiagnostic(("strict-material-mesh:edge:"+std::to_string(edge.first)+":"+std::to_string(edge.second)).c_str());
        // Each point is proved on the same frozen ACTUAL parent edge and in
        // its native analytic cell. Not a root-interval or rounded midpoint.
        nodes.push_back(buildRankTwoMovingMeshEdgePoint(a,parent,request,settings.edgePoints));
    }
    const auto mesh=nodes.empty()?parent:insertRankTwoMovingEdgeNodes(a,parent,nodes,settings.insertion);
    const auto last=nodes.empty()?first:buildRankTwoMaterialLifts(source,a,mesh,settings.lifts);
    if(!last.allActiveStrictPassed()){
        for(const auto& w:last.triangles())if(w.active&&!w.strictPassed)
            throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",
                "strict thickness remains unproved at actual triangle "+std::to_string(w.triangleIndex)
                +"; certified native interior or shared-constraint repair still required");
        throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID","no active material triangles for strict thickness grant");
    }
    auto out=std::make_shared<RankTwoStrictMaterialMesh::Impl>(parent,first,mesh,last,std::move(nodes));
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-strict-material-mesh-v1\nAUTHORITY "<<source.inputHash()<<' '<<a.scientificHash()<<' '
        <<parent.scientificHash()<<' '<<mesh.scientificHash()<<' '<<first.scientificHash()<<' '<<last.scientificHash()<<'\n'
        <<"POLICY ALL_ACTUAL_ZERO_INTERNAL_EDGES FIXED_WEIGHT 1/2 ONE_FROZEN_PARENT\n"
        <<"BUDGET "<<settings.maximumNewEdgePoints<<' '<<settings.maximumInventoryChecks<<" USED "<<checks<<'\n';
    for(const auto& n:out->nodes){const auto& r=n.definition();
        record<<"REPAIR_EDGE "<<r.pointIndices[0]<<' '<<r.pointIndices[1]<<' '<<r.fraction<<' '<<n.scientificHash()<<'\n';}
    record<<"SUMMARY INITIAL_FAILED "<<initialFailed<<" INSERTED_EDGES "<<out->nodes.size()<<" FINAL_FAILED 0\n"
        <<"INITIAL_LIFTS_BEGIN\n"<<first.canonicalRecord()<<"INITIAL_LIFTS_END\n"
        <<"FINAL_LIFTS_BEGIN\n"<<last.canonicalRecord()<<"FINAL_LIFTS_END\n"
        <<"local_actual_strict_thickness PASSED\nglobal_isotopy NOT_GRANTED\n"
        <<"closed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-strict-material-mesh-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoStrictMaterialMesh(std::move(out));
}
} // namespace spin
