#include "detail/RankTwoWholeInteriorInventory.h"
#include "verification/Sha256.h"
#include <locale>
#include <sstream>

namespace spin {
namespace {
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
}
struct RankTwoWholeInteriorInventory::Impl {
    std::string source,authority,normal,record,hash;
    std::vector<RankTwoInteriorInventory> owners;
    std::size_t points=0;
};
RankTwoWholeInteriorInventory::RankTwoWholeInteriorInventory(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoWholeInteriorInventory::sourceIdentity()const{return impl_->source;}
const std::string& RankTwoWholeInteriorInventory::arrangementIdentity()const{return impl_->authority;}
const std::string& RankTwoWholeInteriorInventory::normalIdentity()const{return impl_->normal;}
const std::vector<RankTwoInteriorInventory>& RankTwoWholeInteriorInventory::owners()const{return impl_->owners;}
std::size_t RankTwoWholeInteriorInventory::pointCount()const{return impl_->points;}
const std::string& RankTwoWholeInteriorInventory::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoWholeInteriorInventory::scientificHash()const{return impl_->hash;}

RankTwoWholeInteriorInventory enumerateRankTwoWholeMovingInteriors(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const RankTwoInteriorInventorySettings& settings){
    detail::requireExactAlgebraicEnvironment();
    const auto& whole=wholeArrangementRecords(a); // Reject a local DCEL.
    if(whole.global.sourceIdentity()!=source.inputHash())invalid("whole inventory source/frame authority differs");
    if(whole.normalFaces.size()!=whole.normals.faces().size())invalid("whole inventory normal owner set incomplete");
    if(!settings.maximumCandidateChecks)throw RankTwoSourceError("RESOURCE_LIMIT","positive per-owner inventory budget required");
    auto out=std::make_shared<RankTwoWholeInteriorInventory::Impl>();
    out->source=source.inputHash();out->authority=a.scientificHash();out->normal=whole.normals.inputHash();
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-whole-native-interior-inventory-v1\nAUTHORITY "<<out->source<<' '<<out->authority<<' '<<out->normal
        <<' '<<whole.global.scientificHash()<<' '<<whole.global.truncationRadius()
        <<"\nPOLICY ALL_NATIVE_NORMAL_OWNERS INCLUDING_EMPTY IN_NATIVE_FACE_ORDER\nNORMAL_BEGIN\n"
        <<whole.normals.canonicalRecord()<<"NORMAL_END\n";
    std::size_t activeOwners=0,slabs=0;
    for(std::size_t face=0;face<whole.normalFaces.size();++face){
        const auto& local=whole.normalFaces[face];
        if(local.sourceIdentity()!=out->source||local.truncationRadius()!=whole.global.truncationRadius())
            invalid("normal owner source/radius authority differs");
        if(arrangementRecords(local).normalHalfplanes.size()!=3)invalid("normal owner lacks its complete cone constraints");
        auto perOwner=settings;perOwner.maximumPoints=settings.maximumPoints-out->points;
        auto inventory=enumerateRankTwoMovingInteriors(local,perOwner);
        if(inventory.points().size()>perOwner.maximumPoints)invalid("native inventory exceeded remaining total cap");
        out->points+=inventory.points().size();slabs+=inventory.slabCount();
        activeOwners+=!inventory.points().empty();
        const auto& f=whole.normals.faces()[face];
        record<<"OWNER_BEGIN "<<face<<' '<<local.scientificHash()<<' '<<inventory.scientificHash()<<' '
            <<inventory.points().size()<<' '<<inventory.slabCount()<<' '<<f[0]<<' '<<f[1]<<' '<<f[2]<<'\n'
            <<inventory.canonicalRecord()<<"OWNER_END\n";
        out->owners.push_back(std::move(inventory));
    }
    record<<"COUNTS OWNERS "<<out->owners.size()<<" ACTIVE "<<activeOwners<<" POINTS "<<out->points<<" SLABS "<<slabs
        <<"\nBUDGET TOTAL_POINTS "<<settings.maximumPoints<<" PER_OWNER_CANDIDATE_CHECKS "<<settings.maximumCandidateChecks
        <<"\nwhole_native_owner_inventory PASSED\nshared_mesh_compatibility NOT_GRANTED\ncarrier_insertion NOT_GRANTED\n"
        <<"global_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-whole-native-interior-inventory-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoWholeInteriorInventory(std::move(out));
}
RankTwoWholeInteriorInventory reloadRankTwoWholeMovingInteriors(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const std::string& saved,const RankTwoInteriorInventorySettings& settings){
    const auto fresh=enumerateRankTwoWholeMovingInteriors(source,a,settings);
    if(fresh.canonicalRecord()!=saved)invalid("saved whole inventory differs from ALL native owner reconstruction");return fresh;
}
} // namespace spin
