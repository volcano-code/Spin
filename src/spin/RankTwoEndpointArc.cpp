#include "detail/RankTwoEndpointArc.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <locale>
#include <optional>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=detail::CertifiedAlgebraicReal;
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using P=Traits::Polynomial_2;
using S=detail::RankTwoStripPolynomial;
[[noreturn]] void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]] void unresolved(const char* why){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",why);}
E expression(const Q& q){std::ostringstream s;s<<q;return E(CORE::BigRat(s.str().c_str()));}
S coefficients(const std::string& record){
    S out{};std::istringstream input(record);input.imbue(std::locale::classic());
    std::string term;
    while(std::getline(input,term,';')){
        const auto comma=term.find(','),colon=term.find(':');
        if(comma==std::string::npos || colon==std::string::npos || comma>=colon)
            invalid("invalid immutable supporting polynomial");
        const int i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        if(i<0 || j<0 || i+j>4)invalid("endpoint supporting degree outside contract");
        out[5*i+j]+=Q(term.substr(colon+1));
    }
    return out;
}
P polynomial(const S& p){
    // Immutable supporting records have integer coefficients.
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);
    P out(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0){
        const auto& q=p[5*i+j];if(q.denominator()!=1)invalid("noninteger support record");
        std::ostringstream s;s<<q.numerator();
        out+=P(CORE::BigInt(s.str()))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    return out;
}
bool onSupport(const P& p,const std::array<C,2>& xy){
    auto& k=Traits::CKvA_2::instance().kernel();
    const auto point=k.construct_algebraic_real_2_object()(
        detail::rankTwoAlgebraicHandle(xy[0]),detail::rankTwoAlgebraicHandle(xy[1]));
    return k.is_zero_at_2_object()(CGAL::make_square_free(p),point);
}
}

struct RankTwoMovingGate::Impl {
    std::string arrangement,endpoint,anchor,identity;
    int axis;
    Q fraction;
    C coordinate;
    Impl(std::string a,std::string e,std::string n,int x,Q f,C c,std::string h):
        arrangement(std::move(a)),endpoint(std::move(e)),anchor(std::move(n)),
        identity(std::move(h)),axis(x),fraction(std::move(f)),coordinate(std::move(c)){}
};
RankTwoMovingGate::RankTwoMovingGate(std::shared_ptr<const Impl> impl):impl_(std::move(impl)){}
const std::string& RankTwoMovingGate::arrangementIdentity()const{return impl_->arrangement;}
const std::string& RankTwoMovingGate::endpointIdentity()const{return impl_->endpoint;}
const std::string& RankTwoMovingGate::anchorIdentity()const{return impl_->anchor;}
int RankTwoMovingGate::projectionAxis()const{return impl_->axis;}
const Q& RankTwoMovingGate::fraction()const{return impl_->fraction;}
const C& RankTwoMovingGate::coordinate()const{return impl_->coordinate;}
const std::string& RankTwoMovingGate::identity()const{return impl_->identity;}

RankTwoMovingGate RankTwoMovingGate::fromVertices(const RankTwoArrangement& arrangement,
    const std::string& endpoint,const std::string& anchor,int axis,const Q& fraction,int bits){
    detail::requireExactAlgebraicEnvironment();
    if((axis!=0 && axis!=1) || fraction<=0 || fraction>=1 || bits<1 || bits>512)
        invalid("valid gate axis, interior fraction and isolation budget required");
    return fromVertexWeight(arrangement,endpoint,anchor,axis,fraction,bits);
}
RankTwoMovingGate RankTwoMovingGate::fromReflectedVertices(const RankTwoArrangement& arrangement,
    const std::string& endpoint,const std::string& anchor,int axis,const Q& fraction,int bits){
    detail::requireExactAlgebraicEnvironment();
    if((axis!=0 && axis!=1) || fraction<=0 || fraction>=1 || bits<1 || bits>512)
        invalid("valid reflected gate axis, interior ruler fraction and isolation budget required");
    return fromVertexWeight(arrangement,endpoint,anchor,axis,-fraction,bits);
}
RankTwoMovingGate RankTwoMovingGate::fromVertexWeight(const RankTwoArrangement& arrangement,
    const std::string& endpoint,const std::string& anchor,int axis,const Q& fraction,int bits){
    const auto& records=arrangementRecords(arrangement);
    const auto vertex=[&](const std::string& id)->const RankTwoVertexRecord&{
        for(const auto& v:records.vertices)if(v.identity==id)return v;
        invalid("shared gate event absent from arrangement");
    };
    const C a=C::reload(vertex(endpoint).coordinates[axis],bits);
    const C b=C::reload(vertex(anchor).coordinates[axis],bits);
    if(a.expression()==b.expression())unresolved("shared gate has coincident event projections");
    C coordinate=detail::affineRankTwoAlgebraicCoordinate(a,b,fraction,bits);
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<(fraction>0?"rank-two-moving-gate-v1\n":"rank-two-reflected-moving-gate-v1\n")
        <<arrangement.scientificHash()<<'\n'
        <<endpoint<<'\n'<<anchor<<'\n'<<axis<<' '<<fraction<<'\n'
        <<coordinate.certificate().rootIdentity<<'\n';
    verification::Sha256 hash;const auto bytes=record.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return RankTwoMovingGate(std::make_shared<const Impl>(arrangement.scientificHash(),
        endpoint,anchor,axis,fraction,std::move(coordinate),hash.finish()));
}

RankTwoSharedGateArcCertificate certifyRankTwoEndpointArc(const RankTwoArrangement& arrangement,
    const std::string& id,const RankTwoMovingGate& gate,
    const detail::RankTwoStripBox& box,const Q& distance,int bits,const std::vector<S>& obstacles){
    detail::requireExactAlgebraicEnvironment();
    if(gate.arrangementIdentity()!=arrangement.scientificHash() || bits<1 || bits>512)
        invalid("shared gate authority or isolation budget mismatch");
    const int axis=gate.projectionAxis();
    const auto& records=arrangementRecords(arrangement);
    const RankTwoArcRecord* arc=nullptr;
    for(const auto& a:records.arcs)if(a.identity==id){arc=&a;break;}
    if(!arc)invalid("endpoint arc absent from arrangement");
    int endIndex=-1;
    for(int i=0;i<2;++i)if(arc->ends[i].vertexIdentity==gate.endpointIdentity())endIndex=i;
    if(endIndex<0)unresolved("shared gate endpoint is not incident to arc");
    const auto vertex=[&](int i)->const RankTwoVertexRecord&{
        if(!arc->ends[i].vertexIdentity)invalid("arc has no finite endpoint at requested index");
        for(const auto& v:records.vertices)if(v.identity==*arc->ends[i].vertexIdentity)return v;
        invalid("endpoint vertex absent from arrangement");
    };
    const auto& start=vertex(endIndex);
    const std::array<C,2> startXY{C::reload(start.coordinates[0],bits),C::reload(start.coordinates[1],bits)};
    const S original=coefficients(arc->supportPolynomial);
    S chart{};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)chart[axis?5*j+i:5*i+j]=original[5*i+j];
    const auto strip=detail::certifyRankTwoRegularStrip(chart,box,obstacles,distance);
    const auto inBox=[&](const std::array<C,2>& xy){
        return xy[axis].expression()>=expression(box[0]) && xy[axis].expression()<=expression(box[1])
            && xy[1-axis].expression()>=expression(box[2]) && xy[1-axis].expression()<=expression(box[3]);
    };
    const P nativePolynomial=polynomial(original);
    if(!inBox(startXY) || !onSupport(nativePolynomial,startXY))
        unresolved("actual endpoint not on certified strip graph");
    const C moving=gate.coordinate();
    if(moving.expression()<expression(box[0]) || moving.expression()>expression(box[1]))
        unresolved("moving gate lies outside certified strip");
    Traits traits;
    const auto chartCurve=traits.construct_curve_2_object()(polynomial(chart));
    const auto fiber=chartCurve.status_line_at_exact_x(detail::rankTwoAlgebraicHandle(moving));
    if(fiber.covers_line())unresolved("gate fiber is an entire support component");
    std::optional<C> height;
    for(int branch=0;branch<fiber.number_of_events();++branch){
        const auto value=fiber.algebraic_real_2(branch).y();
        // Select by the whole-strip root proof, NOT the original arc ordinal.
        const auto certified=C::fromRoot(value,bits);
        if(certified.expression()>=expression(box[2]) && certified.expression()<=expression(box[3])){
            if(height)unresolved("multiple graph roots inside certified gate");
            height=certified;
        }
    }
    if(!height)unresolved("certified gate has no supported root");
    std::array<C,2> xy=axis?std::array<C,2>{*height,moving}:std::array<C,2>{moving,*height};
    if(!inBox(xy) || !onSupport(nativePolynomial,xy))
        unresolved("gate failed original support membership");
    const auto withinNativeRange=[&](int coordinate){
        for(int i=0;i<2;++i){
            if(arc->ends[i].location==RankTwoArcEndLocation::Finite){
                const E bound=C::reload(vertex(i).coordinates[coordinate],bits).expression();
                if(i==0?xy[coordinate].expression()<=bound:xy[coordinate].expression()>=bound)return false;
            }else{
                const auto expected=arc->vertical?(i==0?RankTwoArcEndLocation::BottomInfinity:
                    RankTwoArcEndLocation::TopInfinity):(i==0?RankTwoArcEndLocation::LeftInfinity:
                    RankTwoArcEndLocation::RightInfinity);
                if(arc->ends[i].location!=expected){
                    if(!arc->vertical && arc->ends[i].finiteAbscissa
                        && (arc->ends[i].location==RankTwoArcEndLocation::BottomInfinity
                            || arc->ends[i].location==RankTwoArcEndLocation::TopInfinity)){
                        const E bound=C::reload(*arc->ends[i].finiteAbscissa,bits).expression();
                        if(i==0?xy[0].expression()<=bound:xy[0].expression()>=bound)return false;
                    }else invalid("inconsistent native arc boundary");
                }
            }
        }
        return true;
    };
    if(arc->vertical){
        if(xy[0].expression()!=startXY[0].expression() || !withinNativeRange(1))
            unresolved("gate not in original vertical arc");
    }else{
        if(!withinNativeRange(0))
            unresolved("gate not in original arc abscissa interval");
        const auto nativeCurve=traits.construct_curve_2_object()(nativePolynomial);
        const auto nativeFiber=nativeCurve.status_line_at_exact_x(detail::rankTwoAlgebraicHandle(xy[0]));
        if(nativeFiber.covers_line() || arc->branch<0 || arc->branch>=nativeFiber.number_of_events()
            || nativeFiber.algebraic_real_2(arc->branch).y().compare(
                detail::rankTwoAlgebraicHandle(xy[1]))!=CGAL::EQUAL)
            unresolved("gate is on a different original fiber branch");
    }
    // The regular graph minus its start is connected. Since it contains the
    // native gate and contains no other arrangement vertex, it belongs to that
    // same open native arc, including after an isometric chart swap. All
    // vertices count (not just intersections); include gate, exclude start.
    std::size_t checked=0;
    const E a=startXY[axis].expression(),g=moving.expression();
    for(const auto& v:records.vertices){
        ++checked;
        if(v.identity==start.identity)continue;
        const C u=C::reload(v.coordinates[axis],bits);
        const E t=u.expression();
        if(!(g>a?(t>a && t<=g):(t>=g && t<a)))continue;
        const C w=C::reload(v.coordinates[1-axis],bits);
        const std::array<C,2> point=axis?std::array<C,2>{w,u}:std::array<C,2>{u,w};
        if(inBox(point) && onSupport(nativePolynomial,point))
            unresolved("another native event lies between endpoint and gate");
    }
    const int nativeDirection=(xy[arc->vertical?1:0].expression()-startXY[arc->vertical?1:0].expression()).sign();
    const int gateDirection=(g-a).sign();
    if(!nativeDirection || !gateDirection)unresolved("gate direction is not separated from zero");
    const bool aboveIsNativeLeft=(axis==0?1:-1)*nativeDirection*gateDirection>0;
    const std::string above=aboveIsNativeLeft?arc->leftRegion:arc->rightRegion;
    const std::string below=aboveIsNativeLeft?arc->rightRegion:arc->leftRegion;
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"rank-two-shared-gate-arc-v1\n"<<arrangement.scientificHash()<<'\n'
        <<arc->identity<<'\n'<<start.identity<<'\n'<<gate.identity()<<'\n'
        <<strip.inputIdentity<<'\n'<<below<<'\n'<<above<<'\n';
    for(const auto& v:xy)identity<<v.certificate().rootIdentity<<'\n';
    verification::Sha256 hash;const auto bytes=identity.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return {arrangement.scientificHash(),arc->identity,start.identity,gate.identity(),
        axis,std::move(xy),strip,checked,below,above,hash.finish()};
}

RankTwoEndpointArcCertificate certifyRankTwoEndpointArc(const RankTwoArrangement& arrangement,
    const std::string& id,int endIndex,int axis,const Q& fraction,
    const detail::RankTwoStripBox& box,const Q& distance,int bits){
    if((endIndex!=0 && endIndex!=1) || (axis!=0 && axis!=1)
        || fraction<=0 || fraction>=1 || bits<1 || bits>512)
        invalid("valid endpoint, axis, interior fraction and isolation budget required");
    const auto& records=arrangementRecords(arrangement);
    const RankTwoArcRecord* arc=nullptr;
    for(const auto& a:records.arcs)if(a.identity==id){arc=&a;break;}
    if(!arc)invalid("endpoint arc absent from arrangement");
    for(const auto& end:arc->ends)
        if(end.location!=RankTwoArcEndLocation::Finite || !end.vertexIdentity)
            unresolved("endpoint chart needs finite certified endpoints");
    const auto& start=*arc->ends[endIndex].vertexIdentity;
    const auto& other=*arc->ends[1-endIndex].vertexIdentity;
    const auto gate=RankTwoMovingGate::fromVertices(arrangement,start,other,axis,fraction,bits);
    auto proof=certifyRankTwoEndpointArc(arrangement,id,gate,box,distance,bits);
    // Preserve the prior local-arc identity; the new shared-gate identity is
    // separately used by fans. Old evidence never gains the new fan scope.
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"rank-two-endpoint-arc-chart-v1\n"<<arrangement.scientificHash()<<'\n'
        <<arc->identity<<'\n'<<start<<'\n'<<other<<'\n'
        <<axis<<' '<<fraction<<'\n'<<proof.strip.inputIdentity<<'\n';
    for(const auto& v:proof.gateCoordinates)identity<<v.certificate().rootIdentity<<'\n';
    verification::Sha256 hash;const auto bytes=identity.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return {arrangement.scientificHash(),arc->identity,start,other,axis,fraction,
        std::move(proof.gateCoordinates),proof.strip,proof.checkedVertexCount,hash.finish()};
}
} // namespace spin
