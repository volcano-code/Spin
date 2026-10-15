#include "detail/RankTwoInventoryMaterialMesh.h"
#include "verification/Sha256.h"
#include <locale>
#include <sstream>

namespace spin {
struct RankTwoAutomaticMaterialMesh::Impl {
    RankTwoInteriorInventory inventory;
    RankTwoInteriorMesh interior;
    RankTwoStrictMaterialMesh strict;
    std::string record,hash;
    Impl(RankTwoInteriorInventory i,RankTwoInteriorMesh m,RankTwoStrictMaterialMesh s)
        :inventory(std::move(i)),interior(std::move(m)),strict(std::move(s)){}
};
RankTwoAutomaticMaterialMesh::RankTwoAutomaticMaterialMesh(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoAutomaticMaterialMesh::baseMeshIdentity()const{return impl_->interior.baseMesh().scientificHash();}
const RankTwoInteriorInventory& RankTwoAutomaticMaterialMesh::interiors()const{return impl_->inventory;}
const RankTwoInteriorMesh& RankTwoAutomaticMaterialMesh::interiorMesh()const{return impl_->interior;}
const RankTwoStrictMaterialMesh& RankTwoAutomaticMaterialMesh::strictMaterial()const{return impl_->strict;}
const std::string& RankTwoAutomaticMaterialMesh::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoAutomaticMaterialMesh::scientificHash()const{return impl_->hash;}
RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMesh(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const RankTwoCoveredFiniteMesh& base,
    const RankTwoAutomaticMaterialMeshSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(source.inputHash()!=a.sourceIdentity()||base.coverage().arrangementIdentity()!=a.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","automatic material requires SAME source, arrangement and actual covered mesh");
    const auto inventory=enumerateRankTwoMovingInteriors(a,settings.interiors);
    return buildRankTwoAutomaticMaterialMeshFromInventory(source,a,base,inventory,settings.insertion,settings.thickness);
}
RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMeshFromInventory(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const RankTwoCoveredFiniteMesh& base,const RankTwoInteriorInventory& inventory,
    const RankTwoInteriorMeshSettings& insertion,const RankTwoStrictMaterialMeshSettings& thickness){
    detail::requireExactAlgebraicEnvironment();
    if(source.inputHash()!=a.sourceIdentity() || base.coverage().arrangementIdentity()!=a.scientificHash()
        || inventory.arrangementIdentity()!=a.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","supplied native inventory requires SAME source, arrangement and actual covered mesh");
    // Excluded native active owners are rejected by exact insertion. Never
    // drop their points to make a partial covered domain look complete.
    const auto interior=insertRankTwoMovingInteriorNodes(a,base,inventory.points(),insertion);
    const auto strict=completeRankTwoMaterialThickness(source,a,interior,thickness);
    auto out=std::make_shared<RankTwoAutomaticMaterialMesh::Impl>(inventory,interior,strict);
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-automatic-material-mesh-v1\nAUTHORITY "<<source.inputHash()<<' '<<a.scientificHash()<<' '
        <<base.scientificHash()<<' '<<inventory.scientificHash()<<' '<<interior.scientificHash()<<' '<<strict.scientificHash()<<'\n'
        <<"POLICY ALL_NATIVE_GAPS THEN_EXACT_INSERTION THEN_ALL_ZERO_EDGE_REPAIR\n"
        <<"INVENTORY_BEGIN\n"<<inventory.canonicalRecord()<<"INVENTORY_END\n"
        <<"STRICT_BEGIN\n"<<strict.canonicalRecord()<<"STRICT_END\n"
        <<"automatic_local_material_composition PASSED\nglobal_isotopy NOT_GRANTED\n"
        <<"closed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-automatic-material-mesh-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoAutomaticMaterialMesh(std::move(out));
}
} // namespace spin
