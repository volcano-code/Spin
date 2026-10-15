#include "detail/RankTwoEndpointNeighborhood.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <algorithm>
#include <locale>
#include <map>
#include <optional>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using Z=CGAL::Gmpz;
using C=detail::CertifiedAlgebraicReal;
using S=detail::RankTwoStripPolynomial;
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using P=Traits::Polynomial_2;
[[noreturn]] void invalid(const char* s){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",s);}
[[noreturn]] void unresolved(const char* s){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",s);}
template<class T>std::string text(const T& x){std::ostringstream s;s.imbue(std::locale::classic());s<<x;return s.str();}
CORE::BigRat native(const Q& x){return CORE::BigRat(text(x).c_str());}
CORE::Expr expression(const Q& x){return CORE::Expr(native(x));}
S coefficients(const std::string& record){
    S out{};std::istringstream input(record);std::string term;
    while(std::getline(input,term,';')){
        const auto comma=term.find(','),colon=term.find(':');
        if(comma==std::string::npos || colon==std::string::npos || comma>=colon)invalid("bad native support");
        const int i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        if(i<0 || j<0 || i+j>4)invalid("unsupported native support degree");
        out[5*i+j]+=Q(term.substr(colon+1));
    }return out;
}
P polynomial(const S& p){
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);P out(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0){
        if(p[5*i+j].denominator()!=1)invalid("noninteger native polynomial");
        out+=P(CORE::BigInt(text(p[5*i+j].numerator())))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }return out;
}
bool boundary(const S& p,int axis){
    const std::size_t linear=axis?1:5;
    if(p[linear]==0)return false;
    for(std::size_t i=0;i<p.size();++i)if(i!=0 && i!=linear && p[i]!=0)return false;
    return true;
}
S derivative(const S& p,int axis){
    S out{};for(int i=0;i<5;++i)for(int j=0;j<5;++j){
        if(!axis && i)out[5*(i-1)+j]+=i*p[5*i+j];
        if(axis && j)out[5*i+j-1]+=j*p[5*i+j];
    }return out;
}
// Native bounds are only search hints. Exact comparison selects a unique
// FIXED rational grid cell (or an exact grid point), independent of refinement
// history. Neither bound nor midpoint is substituted for a curve coordinate.
std::array<Q,2> enclosure(const C& value,const Q& grid,int bits){
    auto root=detail::rankTwoAlgebraicHandle(value);
    if(root.is_rational()){const Q q(text(root.rational()));return {q,q};}
    for(int attempt=0;attempt<4*bits;++attempt){
        const Q low=Q(text(root.low()))/grid,high=Q(text(root.high()))/grid;
        Z a,b;mpz_fdiv_q(a.mpz(),low.numerator().mpz(),low.denominator().mpz());
        mpz_fdiv_q(b.mpz(),high.numerator().mpz(),high.denominator().mpz());
        if(a==b)return {Q(a)*grid,Q(a+1)*grid};
        if(b==a+1 && root.compare(native(Q(b)*grid))==CGAL::EQUAL)
            return {Q(b)*grid,Q(b)*grid};
        root.refine();
    }throw RankTwoSourceError("RESOURCE_LIMIT","canonical neighborhood enclosure exhausted");
}
void coord(std::ostream& out,const C& v){const auto& c=v.certificate();
    out<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& a:c.coefficients)out<<a<<',';out<<'\n';
}
void gateRecord(std::ostream& out,const RankTwoMovingGate& g){
    out<<"GATE "<<g.identity()<<' '<<g.endpointIdentity()<<' '<<g.anchorIdentity()<<' '
       <<g.projectionAxis()<<' '<<g.fraction()<<'\n';coord(out,g.coordinate());
}
} // namespace
struct RankTwoEndpointNeighborhood::Impl {
    std::string arrangement,endpoint,record,hash;
    int axis=0;
    std::vector<RankTwoEndpointFan> fans;
    std::vector<RankTwoEndpointBoundaryRay> rays;
    std::vector<std::string> incident;
};
RankTwoEndpointNeighborhood::RankTwoEndpointNeighborhood(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoEndpointNeighborhood::arrangementIdentity()const{return impl_->arrangement;}
const std::string& RankTwoEndpointNeighborhood::endpointIdentity()const{return impl_->endpoint;}
int RankTwoEndpointNeighborhood::projectionAxis()const{return impl_->axis;}
const std::vector<RankTwoEndpointFan>& RankTwoEndpointNeighborhood::fans()const{return impl_->fans;}
const std::vector<RankTwoEndpointBoundaryRay>& RankTwoEndpointNeighborhood::boundaryRays()const{return impl_->rays;}
const std::vector<std::string>& RankTwoEndpointNeighborhood::incidentArcs()const{return impl_->incident;}
const std::string& RankTwoEndpointNeighborhood::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoEndpointNeighborhood::scientificHash()const{return impl_->hash;}

RankTwoEndpointNeighborhood certifyRankTwoEndpointNeighborhood(const RankTwoArrangement& a,
    const std::string& id,const RankTwoEndpointNeighborhoodSettings& s){
    detail::requireExactAlgebraicEnvironment();
    if(s.maximumChordDistance<=0 || s.heightMargin<=0 || s.maximumGateHalvings<1
        || s.maximumGateHalvings>64 || s.maximumIsolationBits<1 || s.maximumIsolationBits>512)
        invalid("invalid complete endpoint neighborhood settings");
    const auto& records=arrangementRecords(a);
    const RankTwoVertexRecord* start=nullptr;
    for(const auto& v:records.vertices)if(v.identity==id)start=&v;
    if(!start)invalid("endpoint missing from actual arrangement");
    const std::array<C,2> xy{C::reload(start->coordinates[0],s.maximumIsolationBits),
                           C::reload(start->coordinates[1],s.maximumIsolationBits)};
    std::vector<const RankTwoArcRecord*> incident;
    for(const auto& arc:records.arcs){
        int count=0;for(const auto& end:arc.ends)if(end.vertexIdentity==id)++count;
        if(count>1)unresolved("loop has two endpoint incidences; local graph contract unavailable");
        if(count)incident.push_back(&arc);
    }
    std::sort(incident.begin(),incident.end(),[](auto x,auto y){return x->identity<y->identity;});
    if(incident.empty())unresolved("isolated event has no graph-star neighborhood");
    std::map<std::string,std::array<C,2>> vertices;
    for(const auto& v:records.vertices)vertices.emplace(v.identity,std::array<C,2>{
        C::reload(v.coordinates[0],s.maximumIsolationBits),C::reload(v.coordinates[1],s.maximumIsolationBits)});
    auto& kernel=Traits::CKvA_2::instance().kernel();
    const auto point=kernel.construct_algebraic_real_2_object()(
        detail::rankTwoAlgebraicHandle(xy[0]),detail::rankTwoAlgebraicHandle(xy[1]));
    std::size_t tests=0;std::ostringstream trials,failures;trials.imbue(std::locale::classic());
    failures.imbue(std::locale::classic());
    const auto consume=[&](){
        if(tests==s.maximumCandidateTests)throw RankTwoSourceError("RESOURCE_LIMIT","endpoint candidate budget exhausted");
        ++tests;
    };
    const auto anchor=[&](int axis,int direction)->std::string{
        std::string best;
        for(const auto& v:vertices){
            const auto& value=v.second[axis].expression();
            if(direction<0?!(value<xy[axis].expression()):!(value>xy[axis].expression()))continue;
            if(best.empty() || (direction<0?value>vertices.at(best)[axis].expression():
                   value<vertices.at(best)[axis].expression()))best=v.first;
        }return best;
    };
    const auto fanAnchor=[&](int axis,int direction)->std::string{
        // A sole straight shared graph incidence already has an authoritative
        // opposite native endpoint. Use that dependency, not an unrelated
        // event whose projection happens to coincide at the current source.
        // This preserves the existing no-coordinate-welding publication rule.
        std::string selected;std::size_t matches=0;
        for(const auto* arc:incident){
            if(!std::any_of(arc->owners.begin(),arc->owners.end(),[](const auto& owner){
                return owner.polynomial.rfind("NORMAL_EDGE:",0)==0;}))continue;
            const auto p=coefficients(arc->supportPolynomial);bool line=true;
            for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(i+j>1&&p[5*i+j]!=0)line=false;
            if(!line||boundary(p,axis))continue;
            const auto other=arc->ends[0].vertexIdentity==id?arc->ends[1].vertexIdentity:arc->ends[0].vertexIdentity;
            if(!other||!vertices.count(*other))continue;
            const auto delta=vertices.at(*other)[axis].expression()-xy[axis].expression();
            if(delta.sign()!=direction)continue;
            selected=*other;++matches;
        }
        return matches==1?selected:anchor(axis,direction);
    };
    for(int axis=0;axis<2;++axis){
        bool regular=true;
        for(const auto& f:a.factors()){
            const S p=coefficients(f.polynomial);
            if(!kernel.is_zero_at_2_object()(CGAL::make_square_free(polynomial(p)),point) || boundary(p,axis))continue;
            const P d=polynomial(derivative(p,1-axis));
            if(CGAL::is_zero(d) || kernel.is_zero_at_2_object()(d,point)){regular=false;break;}
        }
        if(!regular){trials<<"CHART "<<axis<<" NONREGULAR\n";continue;}
        auto out=std::make_shared<RankTwoEndpointNeighborhood::Impl>();
        out->arrangement=a.scientificHash();out->endpoint=id;out->axis=axis;
        std::array<detail::RankTwoStripBox,2> fanBoxes;
        bool complete=true;
        for(int side=0;side<2;++side){
            const int direction=side?1:-1;const auto other=fanAnchor(axis,direction);
            if(other.empty()){complete=false;break;}
            bool found=false;Q fraction(1);
            for(int half=1;half<=s.maximumGateHalvings;++half){
                fraction/=2;consume();
                const auto gate=RankTwoMovingGate::fromVertices(a,id,other,axis,fraction,s.maximumIsolationBits);
                const Q grid=s.heightMargin*fraction/8;
                const auto u=enclosure(xy[axis],grid,s.maximumIsolationBits);
                const auto g=enclosure(gate.coordinate(),grid,s.maximumIsolationBits);
                const auto v=enclosure(xy[1-axis],s.heightMargin/8,s.maximumIsolationBits);
                const detail::RankTwoStripBox box{std::min(u[0],g[0]),std::max(u[1],g[1]),
                    v[0]-s.heightMargin,v[1]+s.heightMargin};
                trials<<"TRIAL FAN "<<axis<<' '<<direction<<' '<<half<<' '<<gate.identity();
                for(const auto& b:box)trials<<' '<<b;
                try{
                    out->fans.push_back(certifyRankTwoEndpointFan(a,gate,box,s.maximumChordDistance,s.maximumIsolationBits));
                    fanBoxes[side]=box;trials<<" ACCEPT\n";found=true;break;
                }catch(const RankTwoSourceError& e){
                    trials<<' '<<e.code()<<'\n';if(e.code()!="ARRANGEMENT_UNRESOLVED")throw;
                    failures<<"FAN "<<axis<<' '<<direction<<' '<<half<<' '<<e.witness()<<'\n';
                }
            }
            if(!found){complete=false;break;}
        }
        if(!complete)continue;
        std::set<std::string> covered;
        for(const auto& f:out->fans)for(const auto& arc:f.orderedArcs())
            if(!covered.insert(arc.arcIdentity).second)unresolved("overlapping endpoint half-sectors");
        for(const auto* arc:incident){
            out->incident.push_back(arc->identity);
            if(covered.count(arc->identity))continue;
            if(!boundary(coefficients(arc->supportPolynomial),axis)){complete=false;break;}
            bool found=false;
            for(int direction:{-1,1}){
                auto other=anchor(1-axis,direction);
                const bool reflected=other.empty();
                if(reflected)other=anchor(1-axis,-direction);
                if(other.empty())continue;Q fraction(1);
                for(int half=1;half<=s.maximumGateHalvings;++half){
                    fraction/=2;consume();
                    const auto gate=reflected
                        ?RankTwoMovingGate::fromReflectedVertices(a,id,other,1-axis,fraction,s.maximumIsolationBits)
                        :RankTwoMovingGate::fromVertices(a,id,other,1-axis,fraction,s.maximumIsolationBits);
                    const auto& q=gate.coordinate().expression();
                    if(q<expression(std::max(fanBoxes[0][2],fanBoxes[1][2])) ||
                       q>expression(std::min(fanBoxes[0][3],fanBoxes[1][3]))){
                        trials<<"TRIAL RAY "<<arc->identity<<' '<<direction<<' '<<half<<' '
                            <<gate.identity()<<" SKIP_OUTSIDE_HEIGHT\n";
                        continue;
                    }
                    const Q grid=s.heightMargin*fraction/8;
                    const auto u=enclosure(xy[1-axis],grid,s.maximumIsolationBits);
                    const auto g=enclosure(gate.coordinate(),grid,s.maximumIsolationBits);
                    const auto v=enclosure(xy[axis],s.heightMargin/8,s.maximumIsolationBits);
                    const detail::RankTwoStripBox box{std::min(u[0],g[0]),std::max(u[1],g[1]),
                        v[0]-s.heightMargin,v[1]+s.heightMargin};
                    trials<<"TRIAL RAY "<<arc->identity<<' '<<direction<<' '<<half<<' '<<gate.identity();
                    for(const auto& b:box)trials<<' '<<b;
                    try{
                        auto proof=certifyRankTwoEndpointArc(a,arc->identity,gate,box,s.maximumChordDistance,s.maximumIsolationBits);
                        if(proof.strip.chordDistanceUpperBound!=0)invalid("start-boundary line is not unchanged");
                        out->rays.push_back({gate,std::move(proof),box});covered.insert(arc->identity);
                        trials<<" ACCEPT\n";found=true;break;
                    }catch(const RankTwoSourceError& e){
                        trials<<' '<<e.code()<<'\n';if(e.code()!="ARRANGEMENT_UNRESOLVED")throw;
                        failures<<"RAY "<<arc->identity<<' '<<direction<<' '<<half<<' '<<e.witness()<<'\n';
                    }
                }if(found)break;
            }
            if(!found){complete=false;break;}
        }
        if(!complete || covered!=std::set<std::string>(out->incident.begin(),out->incident.end()))continue;
        std::ostringstream record;record.imbue(std::locale::classic());
        const bool reflected=std::any_of(out->rays.begin(),out->rays.end(),[](const auto& r){return r.gate.fraction()<0;});
        record<<(reflected?"rank-two-complete-endpoint-neighborhood-v2\nAUTHORITY ":
            "rank-two-complete-endpoint-neighborhood-v1\nAUTHORITY ")<<out->arrangement<<' '<<id<<' '<<axis
            <<"\nSETTINGS "<<s.maximumChordDistance<<' '<<s.heightMargin<<' '<<s.maximumGateHalvings
            <<' '<<s.maximumCandidateTests<<' '<<s.maximumIsolationBits<<"\nCANDIDATES "<<tests<<'\n';
        coord(record,xy[0]);coord(record,xy[1]);
        record<<"NATIVE_BEGIN\n"<<a.canonicalRecord()<<"NATIVE_END\n";
        for(const auto* arc:incident)record<<"INCIDENT "<<arc->identity<<' '<<arc->supportPolynomial<<'\n';
        record<<trials.str();
        for(const auto& f:out->fans)record<<"FAN_BEGIN "<<f.scientificHash()<<'\n'<<f.canonicalRecord()<<"FAN_END\n";
        for(const auto& r:out->rays){
            record<<"RAY_BEGIN "<<r.proof.arcIdentity<<' '<<r.proof.identity<<'\n';gateRecord(record,r.gate);
            record<<"RAY_BOX";for(const auto& b:r.box)record<<' '<<b;
            record<<"\nRAY_PROOF "<<r.proof.strip.inputIdentity<<' '<<r.proof.strip.chordDistanceUpperBound
                <<' '<<r.proof.belowRegion<<' '<<r.proof.aboveRegion<<' '<<r.proof.checkedVertexCount<<'\n';
            coord(record,r.proof.gateCoordinates[0]);coord(record,r.proof.gateCoordinates[1]);record<<"RAY_END\n";
        }
        out->record=record.str();const auto bytes="rank-two-endpoint-neighborhood-certificate-v1\n"+out->record;
        verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
        out->hash=hash.finish();detail::requireExactAlgebraicEnvironment();
        return RankTwoEndpointNeighborhood(std::move(out));
    }
    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
        "cannot cover complete native endpoint star within frozen neighborhood contract\n"+failures.str()+trials.str());
}
} // namespace spin
