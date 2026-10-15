#include "detail/RankTwoRegistryEdge.h"
#include "detail/RankTwoDiagnostics.h"
#include "detail/RankTwoArcMembership.h"
#include "detail/RankTwoCertifiedVertices.h"
#include "verification/Sha256.h"
#include <locale>
#include <map>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=detail::CertifiedAlgebraicReal;
using S=detail::RankTwoStripPolynomial;
[[noreturn]] void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
E expression(const Q& q){std::ostringstream out;out.imbue(std::locale::classic());out<<q;return E(CORE::BigRat(out.str().c_str()));}
S coefficients(const std::string& saved){
    return detail::parseRankTwoStripPolynomial(saved,0);
}
S chart(const S& p,int axis){
    S out{};for(int i=0;i<5;++i)for(int j=0;j<5;++j)out[axis?5*j+i:5*i+j]=p[5*i+j];return out;
}
void coordinateRecord(std::ostream& out,const C& value){
    const auto& c=value.certificate();out<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& coefficient:c.coefficients)out<<coefficient<<',';out<<'\n';
}
}
struct RankTwoRegularEdge::Impl {
    std::string arrangement,registry,edge,below,above,record,hash;
    std::array<std::string,2> points;
    std::array<std::array<C,2>,2> coordinates;
    int axis;
    detail::RankTwoStripBox box;
    Q error;
    detail::RankTwoRegularStripCertificate strip;
    RankTwoZeroFactor support;
    std::vector<RankTwoZeroFactor> obstacles;
    std::size_t checked=0;
    Impl(const RankTwoArcRegistry& r,const RankTwoRegistryEdge& e,int a,
        const detail::RankTwoStripBox& b,const Q& limit):
        arrangement(r.arrangementIdentity()),registry(r.scientificHash()),edge(e.identity),
        points{r.points()[e.points[0]].identity,r.points()[e.points[1]].identity},
        coordinates{r.points()[e.points[0]].coordinates,r.points()[e.points[1]].coordinates},
        axis(a),box(b),error(limit){}
};
RankTwoRegularEdge::RankTwoRegularEdge(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoRegularEdge::arrangementIdentity()const{return impl_->arrangement;}
const std::string& RankTwoRegularEdge::registryIdentity()const{return impl_->registry;}
const std::string& RankTwoRegularEdge::edgeIdentity()const{return impl_->edge;}
const std::array<std::string,2>& RankTwoRegularEdge::pointIdentities()const{return impl_->points;}
const std::array<std::array<C,2>,2>& RankTwoRegularEdge::endpointCoordinates()const{return impl_->coordinates;}
int RankTwoRegularEdge::projectionAxis()const{return impl_->axis;}
const detail::RankTwoStripBox& RankTwoRegularEdge::box()const{return impl_->box;}
const Q& RankTwoRegularEdge::maximumChordDistance()const{return impl_->error;}
const detail::RankTwoRegularStripCertificate& RankTwoRegularEdge::strip()const{return impl_->strip;}
const RankTwoZeroFactor& RankTwoRegularEdge::support()const{return impl_->support;}
const std::vector<RankTwoZeroFactor>& RankTwoRegularEdge::obstacles()const{return impl_->obstacles;}
std::size_t RankTwoRegularEdge::checkedVertexCount()const{return impl_->checked;}
const std::string& RankTwoRegularEdge::belowRegion()const{return impl_->below;}
const std::string& RankTwoRegularEdge::aboveRegion()const{return impl_->above;}
const std::string& RankTwoRegularEdge::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoRegularEdge::scientificHash()const{return impl_->hash;}

RankTwoRegularEdgeAttempt tryCertifyRankTwoRegistryEdgeAtIndex(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const std::string& edgeId,std::size_t edgeIndex,int axis,
    const detail::RankTwoStripBox& box,const Q& error,detail::RankTwoRegularStripMemo* memo){
    detail::requireExactAlgebraicEnvironment();
    if((axis!=0 && axis!=1) || registry.arrangementIdentity()!=arrangement.scientificHash())
        invalid("regular edge requires matching arrangement and valid chart");
    if(edgeIndex>=registry.edges().size() || registry.edges()[edgeIndex].identity!=edgeId)
        invalid("actual shared edge absent from immutable registry");
    const auto* edge=&registry.edges()[edgeIndex];
    if(!edge->endpointFanIdentity.empty())invalid("simultaneous fan edge cannot be replaced by independent strip");
    const auto& records=arrangementRecords(arrangement);
    const RankTwoArcRecord* arc=nullptr;
    for(const auto& a:records.arcs)if(a.identity==edge->arcIdentity){arc=&a;break;}
    if(!arc)invalid("registry edge does not belong to the actual native arc");
    const auto& a=registry.points()[edge->points[0]].coordinates,
        &b=registry.points()[edge->points[1]].coordinates;
    const E u0=a[axis].expression(),u1=b[axis].expression();
    if(u0==u1)return std::string("zero chord projection in requested regular chart");
    const auto inside=[&](const std::array<C,2>& xy){
        return xy[axis].expression()>=expression(box[0]) && xy[axis].expression()<=expression(box[1])
            && xy[1-axis].expression()>=expression(box[2]) && xy[1-axis].expression()<=expression(box[3]);
    };
    std::map<std::string,bool> supports;std::vector<S> obstacles;
    std::vector<std::string> obstacleNames;
    const RankTwoZeroFactor* supportFactor=nullptr;
    const auto parsed=[&](const std::string& saved){
        return memo?memo->parsedChart(saved,axis):chart(coefficients(saved),axis);
    };
    for(const auto& factor:arrangement.factors()){
        if(factor.owners.empty() || !supports.emplace(factor.polynomial,true).second)
            invalid("duplicate or unowned immutable regular-edge factor");
        if(factor.polynomial==arc->supportPolynomial)supportFactor=&factor;
        else if(memo)obstacleNames.push_back(factor.polynomial);
        else obstacles.push_back(parsed(factor.polynomial));
    }
    if(!supportFactor)invalid("actual arc supporting factor omitted from arrangement authority");
    auto stripAttempt=memo?memo->tryCertifyChart(supportFactor->polynomial,axis,box,obstacleNames,error):
        detail::tryCertifyRankTwoRegularStrip(parsed(supportFactor->polynomial),box,obstacles,error);
    if(auto* witness=std::get_if<std::string>(&stripAttempt))return std::move(*witness);
    auto strip=std::move(std::get<detail::RankTwoRegularStripCertificate>(stripAttempt));
    detail::rankTwoDiagnostic("regular-edge:strip-ready");
    if(!inside(a) || !inside(b) || !detail::isOnRankTwoNativeArc(arrangement,arc->identity,a)
        || !detail::isOnRankTwoNativeArc(arrangement,arc->identity,b))
        return std::string("actual shared endpoints are not on the certified regular graph");
    detail::rankTwoDiagnostic("regular-edge:endpoints-ready");
    // The unique implicit graph plus absence of all other zero-set factors
    // excludes curve intersections on the ENTIRE strip/chord homotopy. Also
    // check every native event in its actual span, not just sampled points.
    const E lower=u0<u1?u0:u1,upper=u0<u1?u1:u0;
    std::size_t checked=0;
    for(const auto& v:records.vertices){
        ++checked;
        if(v.identity==registry.points()[edge->points[0]].eventIdentity
            || v.identity==registry.points()[edge->points[1]].eventIdentity)continue;
        const auto& xy=certifiedRankTwoVertex(arrangement,v.identity);
        const C& u=xy[axis];
        if(u.expression()<lower || u.expression()>upper)continue;
        if(inside(xy))return std::string("a native event intervenes in the actual regular chord span");
    }
    // Most candidate boxes fail. Only materialize immutable proof metadata
    // AFTER the same fresh native/strip/event checks have actually passed.
    // These factors and coordinates still belong to the current const owner;
    // no geometric result, endpoint or proof authority is cached here.
    auto out=std::make_shared<RankTwoRegularEdge::Impl>(registry,*edge,axis,box,error);
    out->strip=std::move(strip);out->checked=checked;out->support=*supportFactor;
    for(const auto& factor:arrangement.factors())
        if(factor.polynomial!=arc->supportPolynomial)out->obstacles.push_back(factor);
    const bool aboveLeft=(axis==0?1:-1)*(u1-u0).sign()>0;
    out->above=aboveLeft?arc->leftRegion:arc->rightRegion;
    out->below=aboveLeft?arc->rightRegion:arc->leftRegion;
    std::ostringstream scientific,record;scientific.imbue(std::locale::classic());record.imbue(std::locale::classic());
    scientific<<"rank-two-registry-regular-edge-v1\n"<<out->arrangement<<'\n'<<out->registry<<'\n'
        <<edgeId<<'\n'<<arc->identity<<'\n'<<axis<<'\n'<<out->strip.inputIdentity<<'\n'
        <<out->below<<'\n'<<out->above<<'\n';
    record<<"EDGE_AUTHORITY "<<out->arrangement<<' '<<out->registry<<' '<<edgeId<<' '<<arc->identity
        <<' '<<axis<<' '<<out->below<<' '<<out->above<<'\n';
    for(int i=0;i<2;++i){
        scientific<<out->points[i]<<'\n';record<<"EDGE_POINT "<<out->points[i]<<'\n';
        for(const auto& c:out->coordinates[i]){scientific<<c.certificate().rootIdentity<<'\n';coordinateRecord(record,c);}
    }
    record<<"STRIP_BOX";for(const auto& c:box)record<<' '<<c;record<<" ERROR "<<error<<'\n';
    const auto saveFactor=[&](const RankTwoZeroFactor& factor,const char* kind){
        record<<kind<<' '<<factor.polynomial<<' '<<factor.owners.size();
        for(const auto& owner:factor.owners)record<<' '<<owner.polynomial<<':'<<owner.multiplicity;
        record<<'\n'<<"CHART_COEFFICIENTS";for(const auto& c:parsed(factor.polynomial))record<<' '<<c;
        record<<'\n';
    };
    saveFactor(out->support,"GRAPH");for(const auto& other:out->obstacles)saveFactor(other,"OBSTACLE");
    record<<"STRIP_PROOF "<<out->strip.inputIdentity<<' '<<out->strip.yDerivativeSign
        <<' '<<out->strip.derivativeLowerBound<<' '<<out->strip.slopeUpperBound<<' '
        <<out->strip.curvatureUpperBound<<' '<<out->strip.chordDistanceUpperBound;
    for(const auto sign:out->strip.obstacleSigns)record<<' '<<sign;
    record<<'\n'<<"CHECKED_VERTICES "<<out->checked<<'\n'
        <<"full_chord_coverage NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    verification::Sha256 hash;const auto bytes=scientific.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    out->hash=hash.finish();out->record=record.str();
    detail::requireExactAlgebraicEnvironment();return RankTwoRegularEdge(std::move(out));
}
RankTwoRegularEdge certifyRankTwoRegistryEdgeAtIndex(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const std::string& edgeId,std::size_t edgeIndex,int axis,
    const detail::RankTwoStripBox& box,const Q& error,detail::RankTwoRegularStripMemo* memo){
    auto result=tryCertifyRankTwoRegistryEdgeAtIndex(arrangement,registry,edgeId,edgeIndex,axis,box,error,memo);
    if(auto* proof=std::get_if<RankTwoRegularEdge>(&result))return std::move(*proof);
    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",std::get<std::string>(result));
}
RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const std::string& edgeId,int axis,
    const detail::RankTwoStripBox& box,const Q& error,detail::RankTwoRegularStripMemo* memo){
    detail::requireExactAlgebraicEnvironment();
    if((axis!=0 && axis!=1) || registry.arrangementIdentity()!=arrangement.scientificHash())
        invalid("regular edge requires matching arrangement and valid chart");
    for(std::size_t i=0;i<registry.edges().size();++i)if(registry.edges()[i].identity==edgeId)
        return certifyRankTwoRegistryEdgeAtIndex(arrangement,registry,edgeId,i,axis,box,error,memo);
    invalid("actual shared edge absent from immutable registry");
}
RankTwoRegularEdge certifyRankTwoRegistryEdge(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const std::string& edgeId,int axis,
    const detail::RankTwoStripBox& box,const Q& error){
    return certifyRankTwoRegistryEdge(arrangement,registry,edgeId,axis,box,error,nullptr);
}
} // namespace spin
