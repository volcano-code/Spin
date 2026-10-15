#include "detail/RankTwoMaterialLifts.h"
#include "detail/RankTwoCarrierKernels.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using Pair=std::pair<std::size_t,std::size_t>;
template<class T>std::string text(const T& x){std::ostringstream s;s.imbue(std::locale::classic());s<<x;return s.str();}
E exact(const Q& q){return E(CORE::BigRat(CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator()))));}
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
bool active(RankTwoRegionKind k){return k==RankTwoRegionKind::Full || k==RankTwoRegionKind::Cut;}
Pair pair(std::size_t a,std::size_t b){return std::minmax(a,b);}
E evaluate(const std::array<Q,25>& p,const std::array<detail::CertifiedAlgebraicReal,2>& xy){
    std::array<E,5> x,y;x[0]=y[0]=1;
    for(int i=1;i<5;++i){x[i]=x[i-1]*xy[0].expression();y[i]=y[i-1]*xy[1].expression();}
    E v(0);for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0)v+=exact(p[5*i+j])*x[i]*y[j];return v;
}
struct PointLift {E height,sphere,negativeSphere;unsigned zeroMask=0;int hSign=0,dSign=0;};
}
struct RankTwoMaterialLifts::Impl {
    PreparedRankTwoSource source;
    RankTwoArrangement arrangement;
    RankTwoInteriorMesh mesh;
    std::map<std::size_t,PointLift> points;
    std::vector<RankTwoLiftedTriangle> triangles;
    std::string record,hash;
    bool allPassed=false;
    Impl(const PreparedRankTwoSource& s,const RankTwoArrangement& a,const RankTwoInteriorMesh& m):source(s),arrangement(a),mesh(m){}
};
RankTwoMaterialLifts::RankTwoMaterialLifts(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoMaterialLifts::sourceIdentity()const{return impl_->source.inputHash();}
const std::string& RankTwoMaterialLifts::meshIdentity()const{return impl_->mesh.scientificHash();}
const std::vector<RankTwoLiftedTriangle>& RankTwoMaterialLifts::triangles()const{return impl_->triangles;}
const E& RankTwoMaterialLifts::lower(std::size_t r,std::size_t p)const{
    const auto k=impl_->mesh.baseMesh().coverage().partition().registry().regions().at(r).kind;
    if(k==RankTwoRegionKind::Full)return impl_->points.at(p).negativeSphere;
    if(k==RankTwoRegionKind::Cut)return impl_->points.at(p).height;
    topology("inactive region has no material lower lift");
}
const E& RankTwoMaterialLifts::upper(std::size_t p)const{return impl_->points.at(p).sphere;}
bool RankTwoMaterialLifts::allActiveStrictPassed()const{return impl_->allPassed;}
const std::string& RankTwoMaterialLifts::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoMaterialLifts::scientificHash()const{return impl_->hash;}
RankTwoMaterialLifts buildRankTwoMaterialLifts(const PreparedRankTwoSource& source,
    const RankTwoArrangement& a,const RankTwoInteriorMesh& mesh,const RankTwoMaterialLiftSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    detail::rankTwoDiagnostic("material-lift:begin");
    if(!settings.maximumPredicates)throw RankTwoSourceError("RESOURCE_LIMIT","positive material-lift predicate budget required");
    const auto& partition=mesh.baseMesh().coverage().partition();const auto& registry=partition.registry();
    if(source.inputHash()!=a.sourceIdentity() || a.scientificHash()!=partition.arrangementIdentity()
        || registry.arrangementIdentity()!=a.scientificHash())invalid("material lifts require original source, arrangement and actual mesh");
    const auto& graph=source.graph();const auto& native=arrangementRecords(a);const Q radius=a.truncationRadius();
    Q kernelNorm=0;for(const auto& q:graph.kernel)kernelNorm+=q*q;
    if(kernelNorm<=0 || graph.beta<=0)invalid("material lift chart is not an oriented rank-two source");
    std::array<Q,25> height{},projection{};const int powers[6]={0,5,1,10,6,2};
    for(int i=0;i<6;++i)height[powers[i]]=graph.height[i];projection[0]=radius*radius;
    for(int i=0;i<3;++i){
        if(graph.kernel[i]*graph.plane[i][0]!=0 || graph.kernel[i]*graph.plane[i][1]!=0){
            // Individual products need not vanish for a non-axis-aligned chart.
            Q first=0,second=0;for(int j=0;j<3;++j){first+=graph.kernel[j]*graph.plane[j][0];second+=graph.kernel[j]*graph.plane[j][1];}
            if(first!=0 || second!=0)invalid("source chart is not orthogonal to its kernel");
        }
        projection[10]-=graph.plane[i][0]*graph.plane[i][0];
        projection[6]-=2*graph.plane[i][0]*graph.plane[i][1];projection[2]-=graph.plane[i][1]*graph.plane[i][1];
    }
    if(native.originalFhD[0]!=graph.spherePolynomial(radius) || native.originalFhD[1]!=height || native.originalFhD[2]!=projection)
        invalid("native source polynomials differ from the lift chart");
    std::uint64_t predicates=0;const auto consume=[&]{if(predicates==settings.maximumPredicates)
        throw RankTwoSourceError("RESOURCE_LIMIT","material-lift predicate budget exhausted");++predicates;};
    auto out=std::make_shared<RankTwoMaterialLifts::Impl>(source,a,mesh);
    std::map<std::size_t,unsigned> nativeZeros;std::map<Pair,std::size_t> constraints;
    std::map<Pair,std::vector<std::size_t>> uses;
    std::set<std::size_t> usedPoints;
    // Zero equations come only from immutable source-feature incidence, not
    // caller flags or approximate evaluation of independent conjugate roots.
    for(std::size_t i=0;i<registry.edges().size();++i){consume();const auto& e=registry.edges()[i];unsigned mask=0;
        for(const auto& o:e.owners)if(o.multiplicity>0){if(o.polynomial=="F")mask|=1;if(o.polynomial=="h")mask|=2;if(o.polynomial=="D")mask|=4;}
        for(const auto p:e.points)nativeZeros[p]|=mask;
    }
    for(const auto i:partition.edgeIndices()){consume();const auto& e=registry.edges().at(i);
        if(!constraints.emplace(pair(e.points[0],e.points[1]),i).second)topology("duplicate actual native constraint");}
    for(std::size_t i=0;i<mesh.triangles().size();++i){consume();const auto& t=mesh.triangles()[i];
        if(std::find(partition.regionIndices().begin(),partition.regionIndices().end(),t.registryRegionIndex)==partition.regionIndices().end())
            invalid("material triangle refers to an excluded native owner");
        for(int j=0;j<3;++j){usedPoints.insert(t.points[j]);uses[pair(t.points[j],t.points[(j+1)%3])].push_back(i);}}
    for(const auto p:usedPoints){consume();const auto& point=mesh.points().at(p);PointLift lift;
        detail::rankTwoDiagnostic(("material-lift:point:"+std::to_string(p)+":begin").c_str());
        if(p<registry.points().size()){
            if(point.identity!=registry.points()[p].identity)invalid("actual mesh point no longer matches native feature");
            lift.zeroMask=nativeZeros[p];
        }
        lift.height=(lift.zeroMask&2)?E(0):evaluate(height,point.coordinates);lift.hSign=lift.height.sign();
        detail::rankTwoDiagnostic(("material-lift:point:"+std::to_string(p)+":height-sign").c_str());
        const E d=(lift.zeroMask&4)?E(0):evaluate(projection,point.coordinates);lift.dSign=d.sign();
        detail::rankTwoDiagnostic(("material-lift:point:"+std::to_string(p)+":projection-sign").c_str());
        if(lift.dSign<0)topology("actual material mesh point lies outside sphere projection");
        if(lift.zeroMask&4)lift.sphere=E(0);
        else if(lift.zeroMask&1)lift.sphere=lift.hSign<0?-lift.height:lift.height;
        else lift.sphere=sqrt(d/exact(kernelNorm));
        // On F=0,h<0, the Full lower lift IS the same source h expression
        // consumed by the adjacent Cut owner. Reuse it, not -(-h) DAGs.
        lift.negativeSphere=(lift.zeroMask&1) && lift.hSign<0?lift.height:-lift.sphere;
        out->points.emplace(p,std::move(lift));
        detail::rankTwoDiagnostic(("material-lift:point:"+std::to_string(p)+":end").c_str());
    }
    std::map<Pair,bool> boundary;
    for(const auto& [edge,incident]:uses){consume();if(incident.empty() || incident.size()>2)topology("actual mesh edge is non-manifold");
        const auto found=constraints.find(edge);std::size_t activeCount=0;
        for(const auto i:incident)activeCount+=active(registry.regions().at(mesh.triangles()[i].registryRegionIndex).kind);
        if(found!=constraints.end()){
            const auto& e=registry.edges()[found->second];const auto arc=std::find_if(registry.arcs().begin(),registry.arcs().end(),[&](const auto& x){return x.identity==e.arcIdentity;});
            if(arc==registry.arcs().end())invalid("constraint source arc absent");
            std::set<std::string> expected,actual;
            for(const auto r:partition.regionIndices()){const auto& c=registry.regions()[r];if(c.identity==arc->leftRegion || c.identity==arc->rightRegion)expected.insert(c.identity);}
            for(const auto i:incident)actual.insert(registry.regions().at(mesh.triangles()[i].registryRegionIndex).identity);
            if(actual!=expected || actual.size()!=incident.size())topology("actual constraint sides differ from original native adjacency");
        }else if(incident.size()!=2 || mesh.triangles()[incident[0]].registryRegionIndex!=mesh.triangles()[incident[1]].registryRegionIndex)
            topology("nonconstraint edge does not have two same-owner faces");
        boundary.emplace(edge,activeCount==1);
    }
    std::size_t activeCount=0,failed=0;
    for(std::size_t i=0;i<mesh.triangles().size();++i){consume();const auto& t=mesh.triangles()[i];RankTwoLiftedTriangle w;
        detail::rankTwoDiagnostic(("material-lift:triangle:"+std::to_string(i)+":begin").c_str());
        w.triangleIndex=i;w.registryRegionIndex=t.registryRegionIndex;const auto kind=registry.regions()[t.registryRegionIndex].kind;
        w.active=active(kind);
        for(int j=0;j<3;++j)w.materialBoundary[j]=boundary.at(pair(t.points[j],t.points[(j+1)%3]));
        if(w.active){++activeCount;std::array<E,3> thickness;
            for(int j=0;j<3;++j){const auto& point=out->points.at(t.points[j]);const E lo=kind==RankTwoRegionKind::Full?point.negativeSphere:point.height;
                thickness[j]=detail::rankTwoExactLiftDifference(point.sphere,lo);if(thickness[j].sign()<0)topology("actual active corner has negative thickness");}
            w.sumThicknessSign=(thickness[0]+thickness[1]+thickness[2]).sign();
            w.strictPassed=w.sumThicknessSign>0;
            for(int j=0;j<3;++j)if(!w.materialBoundary[j]){
                w.zeroInternalEdges[j]=(thickness[j]+thickness[(j+1)%3]).sign()==0;
                if(w.zeroInternalEdges[j])w.strictPassed=false;}
            failed+=!w.strictPassed;
        }out->triangles.push_back(w);
        detail::rankTwoDiagnostic(("material-lift:triangle:"+std::to_string(i)+":end").c_str());
    }
    out->allPassed=activeCount>0 && failed==0;
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-material-lifts-v1\nAUTHORITY "<<source.inputHash()<<' '<<a.scientificHash()<<' '<<mesh.scientificHash()<<"\nCOEFFICIENTS";
    for(const auto& q:source.exactCoefficients())record<<' '<<q;
    record<<"\nRADIUS "<<radius<<" KERNEL_NORM "<<kernelNorm<<"\nKERNEL";for(const auto& q:graph.kernel)record<<' '<<q;
    for(int i=0;i<3;++i)record<<"\nPLANE "<<i<<' '<<graph.plane[i][0]<<' '<<graph.plane[i][1];
    record<<"\nHEIGHT";for(const auto& q:graph.height)record<<' '<<q;
    for(std::size_t i=0;i<3;++i){record<<"\nPOLYNOMIAL "<<i;for(const auto& q:native.originalFhD[i])record<<' '<<q;}
    record<<"\nBUDGET "<<settings.maximumPredicates<<" USED "<<predicates<<'\n';
    for(const auto& [p,l]:out->points)record<<"POINT_LIFT "<<p<<' '<<l.zeroMask<<' '<<l.hSign<<' '<<l.dSign<<' '
        <<((l.zeroMask&4)?"D_ZERO":(l.zeroMask&1)?"F_ABS_HEIGHT":"POSITIVE_SQRT_D_OVER_K")<<'\n';
    for(const auto& [edge,i]:constraints)record<<"NATIVE_CONSTRAINT "<<edge.first<<' '<<edge.second<<' '<<i<<' '<<registry.edges()[i].identity<<'\n';
    for(const auto& w:out->triangles){record<<"TRIANGLE_LIFT "<<w.triangleIndex<<' '<<w.registryRegionIndex<<' '<<w.active;
        for(const auto b:w.materialBoundary)record<<' '<<b;record<<' '<<w.strictPassed<<' '<<w.sumThicknessSign;
        for(const auto b:w.zeroInternalEdges)record<<' '<<b;record<<'\n';}
    record<<"SUMMARY "<<usedPoints.size()<<' '<<activeCount<<' '<<failed<<"\nall_active_strict "<<out->allPassed
        <<"\nclosed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\nMESH_BEGIN\n"
        <<mesh.canonicalRecord()<<"MESH_END\n";
    out->record=record.str();const auto bytes="rank-two-material-lifts-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoMaterialLifts(std::move(out));
}
} // namespace spin
