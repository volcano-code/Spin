#include "detail/RankTwoMovingInterior.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <locale>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using C=detail::CertifiedAlgebraicReal;using E=CORE::Expr;
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;using P=Traits::Polynomial_2;
template<class T>std::string text(const T& x){std::ostringstream s;s.imbue(std::locale::classic());s<<x;return s.str();}
E exact(const Q& q){return E(CORE::BigRat(CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator()))));}
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
[[noreturn]]void unresolved(const char* why){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",why);}
void coord(std::ostream& s,const C& value){const auto& c=value.certificate();
    s<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& k:c.coefficients)s<<k<<',';s<<'\n';}
P polynomial(const std::array<Q,25>& coefficients){
    CGAL::Gmpz denominator(1);for(const auto& q:coefficients)denominator*=q.denominator();
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);P p(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(coefficients[i*5+j]!=0){
        const Q integer=coefficients[i*5+j]*Q(denominator);
        if(integer.denominator()!=1)invalid("nonintegral positive polynomial scaling");
        p+=P(CORE::BigInt(text(integer.numerator())))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }return p;
}
int sign(const std::array<Q,25>& coefficients,const std::array<C,2>& xy){
    std::array<E,5> x,y;x[0]=y[0]=1;
    for(int k=1;k<5;++k){x[k]=x[k-1]*xy[0].expression();y[k]=y[k-1]*xy[1].expression();}
    E p(0);for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(coefficients[i*5+j]!=0)p+=exact(coefficients[i*5+j])*x[i]*y[j];
    return p.sign();
}
}
struct RankTwoMovingInteriorPoint::Impl {
    std::string authority,record,hash;
    RankTwoInteriorPointRequest request;
    std::array<C,2> coordinates,ordinates;
    RankTwoRegionKind kind;
    Impl(std::string a,RankTwoInteriorPointRequest r,std::array<C,2> xy,std::array<C,2> y,RankTwoRegionKind k)
        :authority(std::move(a)),request(std::move(r)),coordinates(std::move(xy)),ordinates(std::move(y)),kind(k){}
};
RankTwoMovingInteriorPoint::RankTwoMovingInteriorPoint(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoMovingInteriorPoint::arrangementIdentity()const{return impl_->authority;}
const RankTwoInteriorPointRequest& RankTwoMovingInteriorPoint::definition()const{return impl_->request;}
const std::array<C,2>& RankTwoMovingInteriorPoint::coordinates()const{return impl_->coordinates;}
const std::array<C,2>& RankTwoMovingInteriorPoint::fiberOrdinates()const{return impl_->ordinates;}
RankTwoRegionKind RankTwoMovingInteriorPoint::regionKind()const{return impl_->kind;}
const std::string& RankTwoMovingInteriorPoint::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoMovingInteriorPoint::scientificHash()const{return impl_->hash;}
RankTwoMovingInteriorPoint buildRankTwoMovingInteriorPoint(const RankTwoArrangement& a,
    const RankTwoInteriorPointRequest& r,const RankTwoInteriorPointSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(r.abscissaFraction<=0 || r.abscissaFraction>=1 || r.ordinateFraction<=0 || r.ordinateFraction>=1
        || settings.maximumIsolationBits<1 || settings.maximumIsolationBits>512)
        invalid("interior rational fractions and valid isolation budget required");
    if(!settings.maximumNativeChecks || !settings.maximumFiberRoots)
        throw RankTwoSourceError("RESOURCE_LIMIT","positive moving-interior proof budgets required");
    std::uint64_t checks=0,roots=0;
    const auto consume=[&]{if(checks==settings.maximumNativeChecks)
        throw RankTwoSourceError("RESOURCE_LIMIT","moving-interior native check budget exhausted");++checks;};
    const auto& native=arrangementRecords(a);const int bits=settings.maximumIsolationBits;
    const auto event=[&](const std::string& id)->const RankTwoVertexRecord&{
        for(const auto& v:native.vertices){consume();if(v.identity==id)return v;}invalid("interior event absent from native arrangement");};
    const auto arc=[&](const std::string& id)->const RankTwoArcRecord&{
        for(const auto& v:native.arcs){consume();if(v.identity==id)return v;}invalid("interior boundary arc absent from native arrangement");};
    const RankTwoRegionRecord* cell=nullptr;
    for(const auto& v:native.regions){consume();if(v.identity==r.regionIdentity)cell=&v;}
    if(!cell)invalid("interior region absent from native arrangement");
    if(cell->kind!=RankTwoRegionKind::Cut && cell->kind!=RankTwoRegionKind::Full)
        topology("interior node must belong to an active material cell");
    const auto left=C::reload(event(r.eventIdentities[0]).coordinates[0],bits);
    const auto right=C::reload(event(r.eventIdentities[1]).coordinates[0],bits);
    if(!(left.expression()<right.expression()))invalid("native event abscissae not strictly ordered");
    for(const auto& v:native.vertices){consume();const auto x=C::reload(v.coordinates[0],bits);
        if(left.expression()<x.expression() && x.expression()<right.expression())
            unresolved("interior slab crosses an actual native critical/intersection event");}
    const auto x=detail::sampleRankTwoAlgebraicCoordinate(left,right,r.abscissaFraction,bits);
    const auto& lowerArc=arc(r.boundaryArcIdentities[0]);const auto& upperArc=arc(r.boundaryArcIdentities[1]);
    if(lowerArc.vertical || upperArc.vertical || lowerArc.leftRegion!=r.regionIdentity || upperArc.rightRegion!=r.regionIdentity)
        topology("native lower/upper sides do not enclose the declared material cell");
    const auto lower=detail::evaluateRankTwoArcProjection(a,lowerArc.identity,x,bits);
    const auto upper=detail::evaluateRankTwoArcProjection(a,upperArc.identity,x,bits);
    if(!(lower[1].expression()<upper[1].expression()))topology("native fiber ordinates not strictly ordered");
    const auto y=detail::sampleRankTwoAlgebraicCoordinate(lower[1],upper[1],r.ordinateFraction,bits);
    const std::array<C,2> xy{x,y};
    std::vector<std::array<Q,25>> polys(native.originalFhD.begin(),native.originalFhD.end());
    polys.insert(polys.end(),native.normalHalfplanes.begin(),native.normalHalfplanes.end());
    std::ostringstream proof;proof.imbue(std::locale::classic());
    Traits traits;
    for(std::size_t index=0;index<polys.size();++index){consume();
        proof<<"POLYNOMIAL "<<index;for(const auto& q:polys[index])proof<<' '<<q;proof<<'\n';
        const auto p=polynomial(polys[index]);int count=0;
        if(!CGAL::is_zero(p) && CGAL::degree(p)>0){
            const auto curve=traits.construct_curve_2_object()(p);
            const auto fiber=curve.status_line_at_exact_x(detail::rankTwoAlgebraicHandle(x));
            if(fiber.covers_line())throw RankTwoSourceError("DERIVATIVE_BRANCH_EVENT","native vertical fiber is a whole zero component");
            count=fiber.number_of_events();
            for(int k=0;k<count;++k){
                if(roots==settings.maximumFiberRoots)throw RankTwoSourceError("RESOURCE_LIMIT","moving-interior fiber-root budget exhausted");
                ++roots;const auto value=C::fromRoot(fiber.algebraic_real_2(k).y(),bits);
                if(lower[1].expression()<value.expression() && value.expression()<upper[1].expression())
                    unresolved("another native factor lies inside the requested vertical gap");
                proof<<"FIBER_ROOT "<<index<<' '<<k<<'\n';coord(proof,value);
            }
        }
        proof<<"FIBER_COUNT "<<index<<' '<<count<<'\n';
        const auto actual=sign(polys[index],xy);
        const auto expected=index<3?cell->originalSigns[index]:cell->normalSigns.at(index-3);
        if(actual!=expected)topology("moving point signs disagree with original native region");
        proof<<"POINT_SIGN "<<index<<' '<<actual<<'\n';
    }
    if(sign(native.originalFhD[2],xy)<=0 || (cell->kind==RankTwoRegionKind::Cut && sign(native.originalFhD[0],xy)>=0)
        || (cell->kind==RankTwoRegionKind::Full && (sign(native.originalFhD[0],xy)<=0 || sign(native.originalFhD[1],xy)>=0)))
        topology("moving interior node lacks strict active-material sign conditions");
    auto out=std::make_shared<RankTwoMovingInteriorPoint::Impl>(a.scientificHash(),r,xy,std::array<C,2>{lower[1],upper[1]},cell->kind);
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-moving-interior-point-v1\nAUTHORITY "<<a.scientificHash()<<' '<<r.regionIdentity<<' '<<static_cast<int>(cell->kind)
        <<"\nX_BLEND "<<r.eventIdentities[0]<<' '<<r.eventIdentities[1]<<' '<<r.abscissaFraction
        <<"\nLOWER_FIBER "<<lowerArc.identity<<' '<<lowerArc.branch
        <<"\nUPPER_FIBER "<<upperArc.identity<<' '<<upperArc.branch
        <<"\nY_BLEND "<<lower[1].certificate().rootIdentity<<' '<<upper[1].certificate().rootIdentity<<' '<<r.ordinateFraction
        <<"\nBUDGET "<<bits<<' '<<settings.maximumNativeChecks<<' '<<settings.maximumFiberRoots<<" USED "<<checks<<' '<<roots
        <<"\nPOINT_X\n";coord(record,x);record<<"POINT_Y\n";coord(record,y);
    record<<"LOWER_Y\n";coord(record,lower[1]);record<<"UPPER_Y\n";coord(record,upper[1]);
    record<<proof.str()<<"NATIVE_BEGIN\n"<<a.canonicalRecord()<<"NATIVE_END\n"
        <<"native_vertical_gap PASSED\nstrict_active_point PASSED\ncarrier_insertion NOT_GRANTED\n"
        <<"triangle_and_internal_edge_thickness NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-moving-interior-point-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoMovingInteriorPoint(std::move(out));
}
RankTwoMovingInteriorPoint reloadRankTwoMovingInteriorPoint(const RankTwoArrangement& a,
    const RankTwoInteriorPointRequest& r,const std::string& saved,const RankTwoInteriorPointSettings& settings){
    const auto point=buildRankTwoMovingInteriorPoint(a,r,settings);
    if(point.canonicalRecord()!=saved)invalid("saved moving interior definition differs from actual native reconstruction");
    return point;
}
} // namespace spin
