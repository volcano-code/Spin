#include "detail/RankTwoMovingMeshEdge.h"
#include "detail/RankTwoCarrierKernels.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=detail::CertifiedAlgebraicReal;
template<class T>std::string text(const T& x){std::ostringstream s;s.imbue(std::locale::classic());s<<x;return s.str();}
E exact(const Q& q){return E(CORE::BigRat(CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator()))));}
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
int sign(const std::array<Q,25>& coefficients,const std::array<C,2>& xy){
    std::array<E,5> x,y;x[0]=y[0]=1;
    for(int k=1;k<5;++k){x[k]=x[k-1]*xy[0].expression();y[k]=y[k-1]*xy[1].expression();}
    E value(0);for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(coefficients[i*5+j]!=0)
        value+=exact(coefficients[i*5+j])*x[i]*y[j];
    return value.sign();
}
void coord(std::ostream& s,const C& value){const auto& c=value.certificate();
    s<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& k:c.coefficients)s<<k<<',';s<<'\n';}
}
struct RankTwoMovingMeshEdgePoint::Impl {
    std::string authority,parent,record,hash;
    RankTwoMovingMeshEdgeRequest request;
    std::array<C,2> coordinates;
    std::size_t owner;
    Impl(std::string a,std::string p,RankTwoMovingMeshEdgeRequest r,std::array<C,2> xy,std::size_t o)
        :authority(std::move(a)),parent(std::move(p)),request(std::move(r)),coordinates(std::move(xy)),owner(o){}
};
RankTwoMovingMeshEdgePoint::RankTwoMovingMeshEdgePoint(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoMovingMeshEdgePoint::arrangementIdentity()const{return impl_->authority;}
const std::string& RankTwoMovingMeshEdgePoint::parentMeshIdentity()const{return impl_->parent;}
const RankTwoMovingMeshEdgeRequest& RankTwoMovingMeshEdgePoint::definition()const{return impl_->request;}
const std::array<C,2>& RankTwoMovingMeshEdgePoint::coordinates()const{return impl_->coordinates;}
std::size_t RankTwoMovingMeshEdgePoint::registryRegionIndex()const{return impl_->owner;}
const std::string& RankTwoMovingMeshEdgePoint::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoMovingMeshEdgePoint::scientificHash()const{return impl_->hash;}
RankTwoMovingMeshEdgePoint buildRankTwoMovingMeshEdgePoint(const RankTwoArrangement& a,
    const RankTwoInteriorMesh& parent,const RankTwoMovingMeshEdgeRequest& request,
    const RankTwoMovingMeshEdgeSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(!settings.maximumNativeChecks)throw RankTwoSourceError("RESOURCE_LIMIT","positive moving-edge proof budget required");
    if(request.fraction<=0 || request.fraction>=1 || request.pointIndices[0]==request.pointIndices[1]
        || settings.maximumIsolationBits<1 || settings.maximumIsolationBits>512)
        invalid("distinct edge endpoints, strict rational fraction and valid isolation budget required");
    const auto& partition=parent.baseMesh().coverage().partition();const auto& registry=partition.registry();
    if(a.scientificHash()!=partition.arrangementIdentity())invalid("moving edge authority differs from parent mesh");
    auto r=request;
    if(r.pointIndices[0]>=parent.points().size() || r.pointIndices[1]>=parent.points().size())invalid("moving edge endpoint absent");
    if(r.pointIndices[0]>r.pointIndices[1]){std::swap(r.pointIndices[0],r.pointIndices[1]);r.fraction=1-r.fraction;}
    std::uint64_t checks=0;const auto consume=[&]{if(checks==settings.maximumNativeChecks)
        throw RankTwoSourceError("RESOURCE_LIMIT","moving-edge proof budget exhausted");++checks;};
    const auto edge=std::make_pair(r.pointIndices[0],r.pointIndices[1]);
    for(const auto i:partition.edgeIndices()){consume();const auto& e=registry.edges().at(i);
        if(std::make_pair(std::min(e.points[0],e.points[1]),std::max(e.points[0],e.points[1]))==edge)topology("moving edge is an original boundary/interface constraint");}
    std::vector<const RankTwoFiniteTriangle*> incident;
    for(const auto& t:parent.triangles()){consume();for(int j=0;j<3;++j)
        if(std::make_pair(std::min(t.points[j],t.points[(j+1)%3]),std::max(t.points[j],t.points[(j+1)%3]))==edge)incident.push_back(&t);}
    if(incident.size()!=2 || incident[0]->registryRegionIndex!=incident[1]->registryRegionIndex)
        topology("moving edge needs exactly two actual incident faces with the same native owner");
    const auto owner=incident[0]->registryRegionIndex;
    if(std::find(partition.regionIndices().begin(),partition.regionIndices().end(),owner)==partition.regionIndices().end())
        invalid("moving edge owner excluded from covered domain");
    const auto& cell=registry.regions().at(owner);
    if(cell.kind!=RankTwoRegionKind::Full && cell.kind!=RankTwoRegionKind::Cut)topology("moving edge owner is not active material");
    const auto& first=parent.points()[r.pointIndices[0]];const auto& second=parent.points()[r.pointIndices[1]];
    const std::array<C,2> xy{
        detail::sampleRankTwoAlgebraicCoordinate(first.coordinates[0],second.coordinates[0],r.fraction,settings.maximumIsolationBits),
        detail::sampleRankTwoAlgebraicCoordinate(first.coordinates[1],second.coordinates[1],r.fraction,settings.maximumIsolationBits)};
    const auto& native=arrangementRecords(a);
    std::vector<std::array<Q,25>> polys(native.originalFhD.begin(),native.originalFhD.end());
    polys.insert(polys.end(),native.normalHalfplanes.begin(),native.normalHalfplanes.end());
    std::ostringstream proof;proof.imbue(std::locale::classic());
    for(std::size_t i=0;i<polys.size();++i){consume();const auto actual=sign(polys[i],xy);
        const auto expected=i<3?cell.originalSigns[i]:cell.normalSigns.at(i-3);
        if(actual!=expected)topology("moving edge point disagrees with original native material signs");
        proof<<"POLYNOMIAL "<<i;for(const auto& q:polys[i])proof<<' '<<q;
        proof<<"\nPOINT_SIGN "<<i<<' '<<actual<<'\n';}
    if(sign(native.originalFhD[2],xy)<=0 || (cell.kind==RankTwoRegionKind::Cut && sign(native.originalFhD[0],xy)>=0)
        || (cell.kind==RankTwoRegionKind::Full && (sign(native.originalFhD[0],xy)<=0 || sign(native.originalFhD[1],xy)>=0)))
        topology("moving edge point has no strict active thickness");
    auto out=std::make_shared<RankTwoMovingMeshEdgePoint::Impl>(a.scientificHash(),parent.scientificHash(),r,xy,owner);
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-moving-mesh-edge-point-v1\nAUTHORITY "<<a.scientificHash()<<' '<<parent.scientificHash()<<' '<<owner<<'\n'
        <<"EDGE_BLEND "<<r.pointIndices[0]<<' '<<first.identity<<' '<<r.pointIndices[1]<<' '<<second.identity<<' '<<r.fraction<<'\n'
        <<"BUDGET "<<settings.maximumIsolationBits<<' '<<settings.maximumNativeChecks<<" USED "<<checks<<"\nPOINT_X\n";
    coord(record,xy[0]);record<<"POINT_Y\n";coord(record,xy[1]);record<<proof.str()
        <<"actual_internal_edge PASSED\nstrict_active_point PASSED\ncarrier_insertion NOT_GRANTED\n"
        <<"source_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-moving-mesh-edge-point-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoMovingMeshEdgePoint(std::move(out));
}
} // namespace spin
