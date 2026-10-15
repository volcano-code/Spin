#include "detail/RankTwoSharedFiniteSeams.h"
#include "detail/RankTwoCarrierKernels.h"
#include "detail/RankTwoArrangementRecords.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=detail::CertifiedAlgebraicReal;
template<class T>std::string text(const T& x){std::ostringstream o;o.imbue(std::locale::classic());o<<x;return o.str();}
std::string fraction(const Q& x){return text(x.numerator())+"/"+text(x.denominator());}
E expression(const Q& x){return E(CORE::BigRat(CORE::BigInt(text(x.numerator())),CORE::BigInt(text(x.denominator()))));}
std::string hash(const std::string& x){verification::Sha256 s;s.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());return s.finish();}
[[noreturn]]void invalid(const std::string& x){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",x);}
[[noreturn]]void unresolved(const std::string& x){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",x);}
struct Budget {std::uint64_t used=0,limit;void tick(){if(used==limit)throw RankTwoSourceError("RESOURCE_LIMIT","shared seam check budget exhausted");++used;}};
E evaluate(const std::array<Q,25>& p,const E& x,const E& y){
    std::array<E,5> xx{E(1),x,x*x,x*x*x,x*x*x*x},yy{E(1),y,y*y,y*y*y,y*y*y*y};E result(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0)result+=expression(p[5*i+j])*xx[i]*yy[j];
    return result;
}
C rational(const Q& q,int bits){detail::AlgebraicKernel kernel;
    return C::fromRoot(kernel.construct_algebraic_real_1_object()(CORE::BigRat(
        CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator())))),bits);}
struct Endpoint {std::string id;std::array<C,2> xy;};
Endpoint endpoint(const RankTwoWholeArrangementRecords& w,const RankTwoSharedEdgeRecord& edge,const E& t,int bits){
    for(const auto& event:edge.events)if(event.parameter==t)
        return {event.identity,{C::reload(event.point[0],bits),C::reload(event.point[1],bits)}};
    if(t==E(0)||t==E(1)){
        const auto index=edge.edge[t==E(0)?0:1];
        const auto& vertex=w.normalVertices.at(index);
        if(vertex.normalVertex!=index||vertex.location!=RankTwoNormalVertexLocation::FinitePositive||!vertex.point)
            unresolved("finite seam end has no finite normal preimage");
        return {vertex.identity,{rational((*vertex.point)[0],bits),rational((*vertex.point)[1],bits)}};
    }
    unresolved("finite shared interval end is absent from actual normal events");
}
const RankTwoArcRecord& arc(const RankTwoArrangementRecords& r,const std::string& id){
    for(const auto& a:r.arcs)if(a.identity==id)return a;
    unresolved("shared interval names absent native arc");
}
const RankTwoVertexRecord& vertex(const RankTwoArrangementRecords& r,const std::string& id){
    for(const auto& v:r.vertices)if(v.identity==id)return v;
    unresolved("native seam bound names absent vertex");
}
void membership(const RankTwoArrangementRecords& r,const RankTwoArcRecord& a,
    const std::array<C,2>& xy,const std::string& owner,int bits,Budget& budget){
    budget.tick();
    if(std::none_of(a.owners.begin(),a.owners.end(),[&](const auto& o){return o.polynomial==owner;}))
        unresolved("shared consumer lost actual native straight normal owner");
    // The actual NORMAL_EDGE owner certifies the straight-line component of
    // the supporting polynomial, even when that polynomial has other factors.
    // Bound its ORIGINAL monotone branch, not just an arbitrary zero set.
    const int axis=a.vertical?1:0;const E& u=xy[axis].expression();
    for(int k=0;k<2;++k)if(a.ends[k].location==RankTwoArcEndLocation::Finite){
        if(!a.ends[k].vertexIdentity)unresolved("finite native seam end lacks identity");
        const auto& v=vertex(r,*a.ends[k].vertexIdentity);const auto bound=C::reload(v.coordinates[axis],bits);
        const int sign=(u-bound.expression()).sign();
        if((k==0&&sign<0)||(k==1&&sign>0))unresolved("master point outside actual native seam bounds");
        if(sign==0)for(int d=0;d<2;++d)
            if(xy[d].expression()!=C::reload(v.coordinates[d],bits).expression())
                unresolved("master point uses another branch at native endpoint");
    }
    if(a.vertical){const auto& v=vertex(r,*a.ends[0].vertexIdentity);
        if(xy[0].expression()!=C::reload(v.coordinates[0],bits).expression())unresolved("master point left native vertical line");}
}
E parameter(const PreparedRankTwoSource& s,const RankTwoWholeArrangementRecords& w,
    const std::array<std::size_t,2>& edge,const std::array<C,2>& xy){
    const auto& n=w.normals.rays()[edge[0]];const auto& m=w.normals.rays()[edge[1]];
    std::array<E,3> g;
    for(int d=0;d<3;++d){const auto& p=s.graph().gradient[d];g[d]=expression(p[0])+expression(p[1])*xy[0].expression()+expression(p[2])*xy[1].expression();}
    for(int i=0;i<3;++i)for(int j=i+1;j<3;++j){const Q det=n[i]*m[j]-n[j]*m[i];if(det==0)continue;
        const E alpha=(g[i]*expression(m[j])-g[j]*expression(m[i]))/expression(det);
        const E beta=(expression(n[i])*g[j]-expression(n[j])*g[i])/expression(det);
        for(int d=0;d<3;++d)if(alpha*expression(n[d])+beta*expression(m[d])!=g[d])
            unresolved("master point not on actual shared Gauss line");
        if(alpha.sign()<0||beta.sign()<0||(alpha+beta).sign()<=0)unresolved("master point outside actual positive cone");
        return beta/(alpha+beta);
    }
    unresolved("dependent normal edge in admitted mesh");
}
void root(std::ostream& o,const char* axis,const C& x){const auto& c=x.certificate();
    o<<axis<<" ROOT "<<c.rootIdentity<<" INDEX "<<c.rootIndex<<" POLY ";
    for(std::size_t i=0;i<c.coefficients.size();++i)o<<(i?",":"")<<c.coefficients[i];
    o<<" LOWER "<<c.lower<<" UPPER "<<c.upper<<'\n';}
}
struct RankTwoSharedFiniteSeams::Impl {
    RankTwoArrangement whole;std::vector<RankTwoSharedFinitePoint> points;
    std::vector<RankTwoSharedFiniteInterval> intervals;std::string record,science;
    explicit Impl(const RankTwoArrangement& a):whole(a){}
};
RankTwoSharedFiniteSeams::RankTwoSharedFiniteSeams(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoArrangement& RankTwoSharedFiniteSeams::wholeArrangement()const{return impl_->whole;}
const std::vector<RankTwoSharedFinitePoint>& RankTwoSharedFiniteSeams::points()const{return impl_->points;}
const std::vector<RankTwoSharedFiniteInterval>& RankTwoSharedFiniteSeams::intervals()const{return impl_->intervals;}
const std::string& RankTwoSharedFiniteSeams::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoSharedFiniteSeams::scientificHash()const{return impl_->science;}
std::vector<std::size_t> RankTwoSharedFiniteSeams::consumerPoints(const std::string& id,int side)const{
    if(side!=0&&side!=1)invalid("seam consumer must be one of two incident faces");
    for(const auto& i:impl_->intervals)if(i.identity==id){auto result=i.points;
        if(i.sourceReversed[side])std::reverse(result.begin(),result.end());return result;}
    invalid("unknown shared interval consumer");
}
RankTwoSharedFiniteSeams buildRankTwoSharedFiniteSeams(const PreparedRankTwoSource& s,
    const RankTwoArrangement& a,const RankTwoSharedFiniteSeamSettings& settings){
    detail::requireExactAlgebraicEnvironment();const auto& w=wholeArrangementRecords(a);
    if(w.global.sourceIdentity()!=s.inputHash())invalid("shared seams refer to another source");
    if(settings.maximumIsolationBits<1||settings.maximumIsolationBits>512)invalid("invalid shared seam isolation budget");
    if(!settings.maximumPoints||!settings.maximumChecks)throw RankTwoSourceError("RESOURCE_LIMIT","positive seam budgets required");
    Q last(0);for(const auto& f:settings.fractions){if(f<=last||f>=1)invalid("strictly increasing interior fractions required");last=f;}
    if(settings.fractions.empty())invalid("at least one finite interval sample required");
    Budget budget{0,settings.maximumChecks};auto out=std::make_shared<RankTwoSharedFiniteSeams::Impl>(a);
    std::map<std::string,RankTwoSharedFinitePoint> points;
    std::map<std::string,std::vector<std::string>> intervalPoints;
    const auto add=[&](RankTwoSharedFinitePoint p){budget.tick();const auto found=points.find(p.identity);
        if(found!=points.end()){
            for(int d=0;d<2;++d)if(found->second.coordinates[d].expression()!=p.coordinates[d].expression())
                unresolved("same shared feature has unequal exact coordinates");
            if(found->second.endpoint!=p.endpoint||found->second.parentIdentities!=p.parentIdentities||found->second.fraction!=p.fraction)
                unresolved("same shared feature has unequal dependency");
        }else {if(points.size()==settings.maximumPoints)throw RankTwoSourceError("RESOURCE_LIMIT","shared seam point budget exceeded");
            points.emplace(p.identity,std::move(p));}
    };
    for(const auto& edge:w.sharedEdges){budget.tick();std::string authority=s.geometryHash()+";R="+fraction(w.global.truncationRadius());
        for(const auto index:edge.edge){authority+=";V="+std::to_string(index);for(const auto& q:w.normals.rays()[index])authority+=";"+fraction(q);}
        const std::string owner="NORMAL_EDGE:"+std::to_string(edge.edge[0])+":"+std::to_string(edge.edge[1]);
        for(const auto& part:edge.intervals){budget.tick();
            if(part.originalSigns[2]==0)unresolved("entire straight interval on projection ellipse unsupported");
            RankTwoSharedFiniteInterval row{part.identity,edge.edge,edge.incidentFaces,part.sourceArcs,part.sourceReversed,
                part.originalSigns,part.originalOwners,part.originalSigns[2]>0,{}};
            if(row.finite){
                const auto left=endpoint(w,edge,part.lowerParameter,settings.maximumIsolationBits);
                const auto right=endpoint(w,edge,part.upperParameter,settings.maximumIsolationBits);
                add({left.id,true,left.xy,{"",""},Q(0)});add({right.id,true,right.xy,{"",""},Q(0)});
                std::vector<std::string> ids{left.id};
                for(const auto& f:settings.fractions){const std::string id=hash("rank-two-shared-affine-point-v1;"+authority+";"+
                    part.identity+";"+left.id+";"+right.id+";"+fraction(f));
                    std::array<C,2> xy{detail::sampleRankTwoAlgebraicCoordinate(left.xy[0],right.xy[0],f,settings.maximumIsolationBits),
                        detail::sampleRankTwoAlgebraicCoordinate(left.xy[1],right.xy[1],f,settings.maximumIsolationBits)};
                    add({id,false,std::move(xy),{left.id,right.id},f});ids.push_back(id);}
                ids.push_back(right.id);E previous=part.lowerParameter;
                for(std::size_t p=0;p<ids.size();++p){const auto& xy=points.at(ids[p]).coordinates;
                    const E t=parameter(s,w,edge.edge,xy);budget.tick();
                    if((p==0&&t!=part.lowerParameter)||(p+1==ids.size()&&t!=part.upperParameter)
                        ||(p>0&&t<=previous))unresolved("master point Gauss order disagrees with actual shared interval");
                    previous=t;
                    for(int side=0;side<2;++side){const auto& r=arrangementRecords(w.normalFaces.at(edge.incidentFaces[side]));
                        membership(r,arc(r,part.sourceArcs[side]),xy,owner,settings.maximumIsolationBits,budget);
                        for(int k=0;k<3;++k){budget.tick();const int sign=evaluate(r.originalFhD[k],xy[0].expression(),xy[1].expression()).sign();
                            if(p>0&&p+1<ids.size()&&sign!=part.originalSigns[k])unresolved("shared affine sample changed original polynomial sign");
                            if(k==2&&sign<0)unresolved("master finite point outside projection domain");}
                    }
                }
                if(!intervalPoints.emplace(part.identity,std::move(ids)).second)unresolved("duplicate shared interval identity");
            }
            out->intervals.push_back(std::move(row));
        }
    }
    std::map<std::string,std::size_t> indices;
    for(auto& p:points){indices.emplace(p.first,out->points.size());out->points.push_back(std::move(p.second));}
    for(auto& i:out->intervals)if(i.finite)for(const auto& id:intervalPoints.at(i.identity))i.points.push_back(indices.at(id));
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-shared-finite-seams-v1\nAUTHORITY "<<s.inputHash()<<' '<<a.scientificHash()<<' '<<w.normals.inputHash()
        <<"\nRADIUS "<<fraction(w.global.truncationRadius())<<"\nFRACTIONS";
    for(const auto& f:settings.fractions)record<<' '<<fraction(f);record<<"\nPOINTS "<<out->points.size()<<'\n';
    for(const auto& p:out->points){record<<"POINT "<<p.identity<<" ENDPOINT "<<p.endpoint<<" PARENTS "
        <<(p.parentIdentities[0].empty()?"-":p.parentIdentities[0])<<' '<<(p.parentIdentities[1].empty()?"-":p.parentIdentities[1])
        <<" FRACTION "<<fraction(p.fraction)<<'\n';root(record,"X",p.coordinates[0]);root(record,"Y",p.coordinates[1]);}
    for(const auto& edge:w.sharedEdges)record<<"EDGE "<<edge.edge[0]<<' '<<edge.edge[1]<<" INTERVALS "<<edge.intervals.size()<<'\n';
    record<<"INTERVALS "<<out->intervals.size()<<'\n';
    for(const auto& i:out->intervals){record<<"INTERVAL "<<i.identity<<" EDGE "<<i.edge[0]<<' '<<i.edge[1]<<" FINITE "<<i.finite
        <<" SIGNS "<<i.originalSigns[0]<<' '<<i.originalSigns[1]<<' '<<i.originalSigns[2]<<" OWNERS";
        for(const auto& o:i.originalOwners)record<<' '<<o.polynomial<<':'<<o.multiplicity;
        record<<"\nPOINT_INDICES";for(const auto p:i.points)record<<' '<<p;record<<'\n';
        for(int side=0;side<2;++side){const auto& r=arrangementRecords(w.normalFaces[i.incidentFaces[side]]);const auto& v=arc(r,i.sourceArcs[side]);
            record<<"CONSUMER "<<side<<" FACE "<<i.incidentFaces[side]<<" ARC "<<i.sourceArcs[side]<<" REVERSED "<<i.sourceReversed[side]
                <<" LEFT "<<v.leftRegion<<" RIGHT "<<v.rightRegion<<" POINTS";
            auto ids=i.points;if(i.sourceReversed[side])std::reverse(ids.begin(),ids.end());for(const auto p:ids)record<<' '<<p;record<<'\n';}
    }
    record<<"shared_point_publication PASSED\nper_face_mesh_import NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\n"
        <<"strict_thickness NOT_GRANTED\nfull_carrier NOT_GRANTED\nWHOLE_BEGIN\n"<<a.canonicalRecord()<<"WHOLE_END\n";
    out->record=record.str();out->science=hash(out->record);return RankTwoSharedFiniteSeams(std::move(out));
}
} // namespace spin
