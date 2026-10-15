#include "detail/RankTwoEndpointFan.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <algorithm>
#include <locale>
#include <map>
#include <optional>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=detail::CertifiedAlgebraicReal;
using S=detail::RankTwoStripPolynomial;
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using P=Traits::Polynomial_2;
[[noreturn]] void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]] void unresolved(const char* why){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",why);}
E expression(const Q& q){std::ostringstream s;s<<q;return E(CORE::BigRat(s.str().c_str()));}
S coefficients(const std::string& record){
    S out{};std::istringstream input(record);input.imbue(std::locale::classic());std::string term;
    while(std::getline(input,term,';')){
        const auto comma=term.find(','),colon=term.find(':');
        if(comma==std::string::npos || colon==std::string::npos || comma>=colon)
            invalid("malformed immutable fan support");
        const int i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        if(i<0 || j<0 || i+j>4)invalid("fan support degree outside contract");
        out[5*i+j]+=Q(term.substr(colon+1));
    }
    return out;
}
P polynomial(const S& p){
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);P out(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0){
        if(p[5*i+j].denominator()!=1)invalid("noninteger immutable fan support");
        std::ostringstream s;s<<p[5*i+j].numerator();
        out+=P(CORE::BigInt(s.str()))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    if(CGAL::is_zero(out))invalid("zero immutable fan support");
    return out;
}
S swapped(const S& p,int axis){
    S out{};for(int i=0;i<5;++i)for(int j=0;j<5;++j)out[axis?5*j+i:5*i+j]=p[5*i+j];return out;
}
bool boundaryPlane(const S& p){
    if(p[5]==0)return false;
    for(std::size_t i=0;i<p.size();++i)if(i!=0 && i!=5 && p[i]!=0)return false;
    return true; // p(u,v)=a*u+b, a!=0; actual start membership is checked.
}
void coordinateRecord(std::ostream& out,const C& value){
    const auto& c=value.certificate();
    out<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& coefficient:c.coefficients)out<<coefficient<<',';
    out<<'\n';
}
}

struct RankTwoEndpointFan::Impl {
    std::optional<RankTwoMovingGate> gate;
    std::vector<RankTwoSharedGateArcCertificate> arcs;
    std::vector<RankTwoFanSupportRecord> supports;
    std::string record,hash;
};
RankTwoEndpointFan::RankTwoEndpointFan(std::shared_ptr<const Impl> impl):impl_(std::move(impl)){}
const RankTwoMovingGate& RankTwoEndpointFan::gate()const{return *impl_->gate;}
const std::vector<RankTwoSharedGateArcCertificate>& RankTwoEndpointFan::orderedArcs()const{return impl_->arcs;}
const std::vector<RankTwoFanSupportRecord>& RankTwoEndpointFan::supportRecords()const{return impl_->supports;}
const std::string& RankTwoEndpointFan::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoEndpointFan::scientificHash()const{return impl_->hash;}

RankTwoEndpointFan certifyRankTwoEndpointFan(const RankTwoArrangement& arrangement,
    const RankTwoMovingGate& gate,const detail::RankTwoStripBox& box,const Q& distance,int bits){
    detail::requireExactAlgebraicEnvironment();
    if(gate.arrangementIdentity()!=arrangement.scientificHash() || bits<1 || bits>512)
        invalid("fan gate authority or isolation budget mismatch");
    if(box[0]>=box[1] || box[2]>=box[3] || distance<=0)
        invalid("ordered fan rectangle and positive accuracy required");
    const auto& records=arrangementRecords(arrangement);
    const RankTwoVertexRecord* start=nullptr;
    for(const auto& v:records.vertices)if(v.identity==gate.endpointIdentity()){start=&v;break;}
    if(!start)invalid("fan endpoint absent from immutable arrangement");
    const int axis=gate.projectionAxis();
    const std::array<C,2> startXY{C::reload(start->coordinates[0],bits),C::reload(start->coordinates[1],bits)};
    if(startXY[axis].expression()<expression(box[0]) || startXY[axis].expression()>expression(box[1])
        || startXY[1-axis].expression()<expression(box[2]) || startXY[1-axis].expression()>expression(box[3])
        || gate.coordinate().expression()<expression(box[0]) || gate.coordinate().expression()>expression(box[1]))
        unresolved("fan endpoint or common gate outside rectangle");
    auto& k=Traits::CKvA_2::instance().kernel();
    const auto point=k.construct_algebraic_real_2_object()(
        detail::rankTwoAlgebraicHandle(startXY[0]),detail::rankTwoAlgebraicHandle(startXY[1]));
    struct Factor {RankTwoZeroFactor original;S chart;RankTwoFanSupportKind kind;};
    std::vector<Factor> factors;
    std::vector<S> obstacles;
    std::vector<std::size_t> excludedIndices;
    std::map<std::string,bool> supportIdentities;
    // This list is owned by the private, fully constructed arrangement.
    // The caller does not select a subset of arcs or normal owners.
    for(const auto& factor:arrangement.factors()){
        if(!supportIdentities.emplace(factor.polynomial,true).second || factor.owners.empty())
            invalid("duplicate or unowned immutable fan factor");
        const S original=coefficients(factor.polynomial),chart=swapped(original,axis);
        const bool containsStart=k.is_zero_at_2_object()(CGAL::make_square_free(polynomial(original)),point);
        const auto kind=!containsStart?RankTwoFanSupportKind::Excluded:
            boundaryPlane(chart)?RankTwoFanSupportKind::AtStartBoundary:RankTwoFanSupportKind::Graph;
        if(kind==RankTwoFanSupportKind::Excluded){excludedIndices.push_back(factors.size());obstacles.push_back(chart);}
        factors.push_back({factor,chart,kind});
    }
    if(factors.empty())invalid("fan authority has no zero-set factors");
    auto out=std::make_shared<RankTwoEndpointFan::Impl>();
    out->gate=gate;
    std::vector<int> exclusionSigns;
    for(const auto& factor:factors){
        RankTwoFanSupportRecord saved{factor.original.polynomial,factor.original.owners,factor.kind,0,{}};
        if(factor.kind==RankTwoFanSupportKind::Graph){
            // Prove complete root/obstacle coverage on the rectangle FIRST.
            // Incident factors meet at the start, so they require the ordered
            // fan proof below instead of an impossible positive separation.
            const auto strip=detail::certifyRankTwoRegularStrip(factor.chart,box,obstacles,distance);
            if(exclusionSigns.empty())exclusionSigns=strip.obstacleSigns;
            else if(exclusionSigns!=strip.obstacleSigns)invalid("inconsistent whole-box obstacle signs");
            std::optional<RankTwoSharedGateArcCertificate> selected;
            for(const auto& arc:records.arcs){
                if(arc.supportPolynomial!=factor.original.polynomial)continue;
                bool incident=false;
                for(const auto& end:arc.ends)incident|=end.vertexIdentity==start->identity;
                if(!incident)continue;
                try{
                    auto proof=certifyRankTwoEndpointArc(arrangement,arc.identity,gate,box,distance,bits,obstacles);
                    if(selected)throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",
                        "common graph gate belongs to multiple native incident arcs");
                    if(proof.strip.inputIdentity!=strip.inputIdentity)
                        invalid("arc and fan use different strip/obstacle proofs");
                    selected=std::move(proof);
                }catch(const RankTwoSourceError& e){
                    if(e.code()!="ARRANGEMENT_UNRESOLVED")throw;
                }
            }
            if(!selected)unresolved("incident graph cannot be bound to its original native arc");
            saved.arcIdentity=selected->arcIdentity;
            out->arcs.push_back(std::move(*selected));
        }
        out->supports.push_back(std::move(saved));
    }
    if(out->arcs.empty())unresolved("no regular incident graph in this fan chart");
    if(exclusionSigns.size()!=excludedIndices.size())invalid("incomplete nonincident factor exclusion");
    for(std::size_t i=0;i<excludedIndices.size();++i)
        out->supports[excludedIndices[i]].exclusionSign=exclusionSigns[i];
    const int heightAxis=1-axis;
    std::sort(out->arcs.begin(),out->arcs.end(),[&](const auto& a,const auto& b){
        return a.gateCoordinates[heightAxis].expression()<b.gateCoordinates[heightAxis].expression();
    });
    for(std::size_t i=1;i<out->arcs.size();++i){
        const auto& a=out->arcs[i-1];const auto& b=out->arcs[i];
        if(!(a.gateCoordinates[heightAxis].expression()<b.gateCoordinates[heightAxis].expression()))
            throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID","incident graphs meet again at common gate");
        if(a.aboveRegion!=b.belowRegion)
            throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID","native region adjacency differs from ordered fan");
    }
    // Each pair has no intervening arrangement vertex by the actual-arc
    // certificates. Hence its graphs stay strictly ordered in the open
    // half-sector. Common-gate chords have the same strict order, and the
    // simultaneous interpolation homotopy stays ordered for every u!=start.
    std::ostringstream scientific,record;
    scientific.imbue(std::locale::classic());record.imbue(std::locale::classic());
    scientific<<"rank-two-common-endpoint-fan-v1\n"<<arrangement.scientificHash()<<'\n'
        <<gate.identity()<<'\n';
    for(const auto& b:box)scientific<<b<<' ';scientific<<'\n'<<distance<<'\n';
    record<<"FAN_AUTHORITY "<<arrangement.scientificHash()<<' '<<gate.identity()<<' '
        <<gate.endpointIdentity()<<' '<<gate.anchorIdentity()<<' '<<axis<<' '<<gate.fraction()<<'\n';
    record<<"FAN_BOX";for(const auto& b:box)record<<' '<<b;record<<" ERROR "<<distance<<'\n';
    coordinateRecord(record,gate.coordinate());
    for(const auto& support:out->supports){
        const char* kind=support.kind==RankTwoFanSupportKind::Graph?"GRAPH":
            support.kind==RankTwoFanSupportKind::AtStartBoundary?"START_BOUNDARY":"EXCLUDED";
        scientific<<kind<<' '<<support.polynomial<<' '<<support.arcIdentity<<' '
            <<support.exclusionSign<<'\n';
        record<<"SUPPORT "<<kind<<' '<<support.polynomial<<' '<<support.owners.size();
        for(const auto& owner:support.owners){
            scientific<<"OWNER "<<owner.polynomial<<' '<<owner.multiplicity<<'\n';
            record<<' '<<owner.polynomial<<':'<<owner.multiplicity;
        }
        if(support.kind==RankTwoFanSupportKind::Excluded)record<<" SIGN "<<support.exclusionSign;
        if(support.kind==RankTwoFanSupportKind::Graph)record<<" ARC "<<support.arcIdentity;
        record<<'\n';
    }
    for(const auto& arc:out->arcs){
        scientific<<arc.identity<<'\n'<<arc.strip.inputIdentity<<'\n';
        for(const auto& v:arc.gateCoordinates)scientific<<v.certificate().rootIdentity<<'\n';
        record<<"ORDERED_ARC "<<arc.arcIdentity<<' '<<arc.identity<<' '<<arc.gateIdentity<<' '
            <<arc.belowRegion<<' '<<arc.aboveRegion<<' '<<arc.strip.inputIdentity<<' '
            <<arc.strip.chordDistanceUpperBound<<' '<<arc.checkedVertexCount<<'\n';
        for(const auto& v:arc.gateCoordinates)coordinateRecord(record,v);
    }
    verification::Sha256 hash;const auto bytes=scientific.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    out->hash=hash.finish();out->record=record.str();
    detail::requireExactAlgebraicEnvironment();
    return RankTwoEndpointFan(std::move(out));
}
} // namespace spin
