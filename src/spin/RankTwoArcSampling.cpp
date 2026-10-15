#include "detail/RankTwoArcSampling.h"
#include "detail/RankTwoArcMembership.h"
#include "detail/RankTwoDiagnostics.h"
#include "detail/RankTwoNativeFiber.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <locale>
#include <sstream>

namespace spin {
namespace {
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using P=Traits::Polynomial_2;
using C=detail::CertifiedAlgebraicReal;
[[noreturn]] void invalid(const char* message){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",message);}
P polynomial(const std::string& record){
    // Only a private, constructed/replayed arrangement's record reaches here.
    // Keep its WHOLE supporting polynomial: factor replacement may change
    // the native interior fiber index even for the same geometric component.
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);
    P out(0);std::istringstream input(record);input.imbue(std::locale::classic());
    std::string term;
    while(std::getline(input,term,';')){
        const auto comma=term.find(','),colon=term.find(':');
        if(comma==std::string::npos || colon==std::string::npos || comma>=colon)
            invalid("invalid supporting polynomial term");
        const int i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        if(i<0 || j<0 || i>4 || j>4)invalid("supporting polynomial degree outside arrangement contract");
        out+=P(CORE::BigInt(term.substr(colon+1)))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    if(CGAL::is_zero(out))invalid("zero supporting polynomial");
    return out;
}
}

RankTwoArcSample sampleRankTwoArc(const RankTwoArrangement& arrangement,const std::string& id,
    const CGAL::Gmpq& fraction,int bits){
    detail::requireExactAlgebraicEnvironment();
    if(fraction<=0 || fraction>=1 || bits<1 || bits>512)invalid("interior sample fraction and valid isolation budget required");
    const auto& records=arrangementRecords(arrangement);
    const RankTwoArcRecord* arc=nullptr;
    for(const auto& candidate:records.arcs)if(candidate.identity==id){arc=&candidate;break;}
    if(!arc)invalid("sample arc absent from arrangement");
    for(const auto& end:arc->ends)
        if(end.location!=RankTwoArcEndLocation::Finite || !end.vertexIdentity)
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","unbounded arc requires certified finite clipping before sampling");
    const auto endpoint=[&](int i)->const RankTwoVertexRecord&{
        for(const auto& v:records.vertices)if(v.identity==*arc->ends[i].vertexIdentity)return v;
        invalid("sample endpoint absent from arrangement");
    };
    const auto& left=endpoint(0);const auto& right=endpoint(1);
    detail::rankTwoDiagnostic("arc-sample:endpoint-reload");
    const int axis=arc->vertical?1:0;
    const auto a=C::reload(left.coordinates[axis],bits),b=C::reload(right.coordinates[axis],bits);
    if(!(a.expression()<b.expression()))invalid("sample arc endpoints are not ordered on their monotone axis");
    detail::rankTwoDiagnostic("arc-sample:weighted-coordinate");
    const auto moving=detail::sampleRankTwoAlgebraicCoordinate(a,b,fraction,bits);
    detail::rankTwoDiagnostic("arc-sample:weighted-coordinate-complete");
    const auto other=[&](){
        if(arc->vertical){
            auto x=C::reload(left.coordinates[0],bits);
            if(x.expression()!=C::reload(right.coordinates[0],bits).expression())invalid("vertical arc endpoints have different abscissae");
            return x;
        }
        Traits traits;
        const auto curve=traits.construct_curve_2_object()(polynomial(arc->supportPolynomial));
        const auto alpha=detail::rankTwoAlgebraicHandle(moving);
        // An interval representative section can lie at a DIFFERENT x. Use
        // the exact algebraic abscissa and the saved INTERIOR branch ordinal;
        // endpoint branch ordinals are not interchangeable with this one.
        detail::rankTwoDiagnostic("arc-sample:exact-fiber");
        const auto line=curve.status_line_at_exact_x(alpha);
        if(line.covers_line() || arc->branch<0 || arc->branch>=line.number_of_events())
            throw RankTwoSourceError("DERIVATIVE_BRANCH_EVENT","arc interior sample has no valid isolated fiber branch");
        const auto point=line.algebraic_real_2(arc->branch);
        detail::rankTwoDiagnostic("arc-sample:y-handle");
        std::optional<detail::RankTwoNativeRootSet> proof;
        const auto y=detail::rankTwoNativeFiberOrdinate(point,bits,&proof);
        detail::rankTwoDiagnostic("arc-sample:y-certificate");
        const auto result=proof?C::fromRoot(y,*proof,bits):C::fromRoot(y,bits);
        detail::rankTwoDiagnostic("arc-sample:y-certificate-complete");
        return result;
    }();
    std::array<C,2> coordinates=arc->vertical?std::array<C,2>{other,moving}:std::array<C,2>{moving,other};
    std::ostringstream canonical;canonical.imbue(std::locale::classic());
    canonical<<"rank-two-moving-arc-sample-v1\n"<<arrangement.scientificHash()<<'\n'
        <<arc->identity<<'\n'<<left.identity<<'\n'<<right.identity<<'\n'
        <<fraction.numerator()<<'/'<<fraction.denominator()<<'\n'<<arc->vertical<<' '<<arc->branch<<'\n';
    verification::Sha256 hash;const auto bytes=canonical.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return {arrangement.scientificHash(),arc->identity,{left.identity,right.identity},fraction,
        arc->branch,arc->vertical,std::move(coordinates),hash.finish()};
}
namespace detail {
bool isOnRankTwoNativeArc(const RankTwoArrangement& arrangement,const std::string& id,
    const std::array<C,2>& xy){
    requireExactAlgebraicEnvironment();
    const auto& records=arrangementRecords(arrangement);const RankTwoArcRecord* arc=nullptr;
    for(const auto& a:records.arcs)if(a.identity==id)arc=&a;
    if(!arc)invalid("membership arc absent from native arrangement");
    for(const auto& end:arc->ends)if(end.location!=RankTwoArcEndLocation::Finite || !end.vertexIdentity)
        throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","native membership requires certified finite bounds");
    const auto endpoint=[&](int k)->const RankTwoVertexRecord&{
        for(const auto& v:records.vertices)if(v.identity==*arc->ends[k].vertexIdentity)return v;
        invalid("native membership endpoint absent");
    };
    const int axis=arc->vertical?1:0;
    const auto u=rankTwoAlgebraicHandle(xy[axis]);
    const auto left=rankTwoAlgebraicHandle(C::reload(endpoint(0).coordinates[axis]));
    const auto right=rankTwoAlgebraicHandle(C::reload(endpoint(1).coordinates[axis]));
    if(left.compare(right)!=CGAL::SMALLER)invalid("native membership endpoints not monotone");
    const auto lower=u.compare(left),upper=u.compare(right);
    if(lower==CGAL::SMALLER || upper==CGAL::LARGER)return false;
    if(lower==CGAL::EQUAL || upper==CGAL::EQUAL){
        const auto& v=endpoint(lower==CGAL::EQUAL?0:1);
        // A zero at the same abscissa on another branch is NOT this endpoint.
        for(int k=0;k<2;++k)if(rankTwoAlgebraicHandle(xy[k]).compare(
            rankTwoAlgebraicHandle(C::reload(v.coordinates[k])))!=CGAL::EQUAL)return false;
        return true;
    }
    if(arc->vertical)return rankTwoAlgebraicHandle(xy[0]).compare(
        rankTwoAlgebraicHandle(C::reload(endpoint(0).coordinates[0])))==CGAL::EQUAL;
    Traits traits;const auto curve=traits.construct_curve_2_object()(polynomial(arc->supportPolynomial));
    const auto line=curve.status_line_at_exact_x(u);
    if(line.covers_line() || arc->branch<0 || arc->branch>=line.number_of_events())
        throw RankTwoSourceError("DERIVATIVE_BRANCH_EVENT","native membership has no valid interior fiber branch");
    return rankTwoNativeFiberMatches(line.algebraic_real_2(arc->branch),xy[1]);
}
std::array<C,2> evaluateRankTwoArcProjection(const RankTwoArrangement& arrangement,
    const std::string& id,const C& projection,int bits){
    requireExactAlgebraicEnvironment();
    if(bits<1 || bits>512)invalid("valid projection isolation budget required");
    const auto& records=arrangementRecords(arrangement);const RankTwoArcRecord* arc=nullptr;
    for(const auto& a:records.arcs)if(a.identity==id)arc=&a;
    if(!arc)invalid("projection arc absent from arrangement");
    for(const auto& end:arc->ends)if(end.location!=RankTwoArcEndLocation::Finite || !end.vertexIdentity)
        throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","subdivision requires certified finite clipping");
    const auto vertex=[&](int end)->const RankTwoVertexRecord&{
        for(const auto& v:records.vertices)if(v.identity==*arc->ends[end].vertexIdentity)return v;
        invalid("projection arc endpoint absent from arrangement");
    };
    const int axis=arc->vertical?1:0;
    const auto lower=C::reload(vertex(0).coordinates[axis],bits),upper=C::reload(vertex(1).coordinates[axis],bits);
    rankTwoDiagnostic("arc-projection:order-check");
    if(!(lower.expression()<projection.expression() && projection.expression()<upper.expression()))
        invalid("projection outside original native open arc interval");
    if(arc->vertical){
        const auto x=C::reload(vertex(0).coordinates[0],bits);
        if(x.expression()!=C::reload(vertex(1).coordinates[0],bits).expression())
            invalid("vertical projection arc has inconsistent abscissae");
        return {x,projection};
    }
    rankTwoDiagnostic("arc-projection:construct-curve");
    Traits traits;const auto curve=traits.construct_curve_2_object()(polynomial(arc->supportPolynomial));
    rankTwoDiagnostic("arc-projection:exact-fiber");
    const auto line=curve.status_line_at_exact_x(rankTwoAlgebraicHandle(projection));
    if(line.covers_line() || arc->branch<0 || arc->branch>=line.number_of_events())
        throw RankTwoSourceError("DERIVATIVE_BRANCH_EVENT","subdivision has no valid native isolated fiber branch");
    rankTwoDiagnostic("arc-projection:y-handle");
    std::optional<RankTwoNativeRootSet> proof;
    const auto native=rankTwoNativeFiberOrdinate(line.algebraic_real_2(arc->branch),bits,&proof);
    const auto y=proof?C::fromRoot(native,*proof,bits):C::fromRoot(native,bits);
    rankTwoDiagnostic("arc-projection:complete");
    requireExactAlgebraicEnvironment();return {projection,y};
}
}
} // namespace spin
