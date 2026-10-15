#include "detail/RankTwoConstruction.h"
#include "detail/RankTwoNormalOverlay.h"
#include "detail/RankTwoNormalEvents.h"
#include "detail/RankTwoNormalVertices.h"
#include "detail/RankTwoWholeArrangement.h"
#include "detail/RankTwoArrangementRecords.h"
#include "detail/RankTwoAlgebraic.h"
#include "detail/RankTwoCertifiedVertices.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <algorithm>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using Z=CGAL::Gmpz;
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using P=Traits::Polynomial_2;
using PT=CGAL::Polynomial_traits_d<P>;
using AR=Traits::Algebraic_real_1;
using Bound=Traits::Bound;
using Coefficients=std::array<Q,25>;
using Dcel=CGAL::Arr_extended_dcel<Traits,std::size_t,std::size_t,std::size_t>;
using Arrangement=CGAL::Arrangement_2<Traits,Dcel>;
template<class T> std::string text(const T& x) {
    std::ostringstream out;out.imbue(std::locale::classic());out<<x;return out.str();
}
std::string fraction(const Q& x){return text(x.numerator())+"/"+text(x.denominator());}
std::string sha(const std::string& s) {
    verification::Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());return hash.finish();
}
[[noreturn]] void unresolved(const char* reason) {
    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",reason);
}
void checkBudget(const Arrangement& a,const RankTwoSourceLimits& limits) {
    if(a.number_of_faces()+a.number_of_vertices()+a.number_of_edges()>limits.maximumArrangementCells)
        throw RankTwoSourceError("RESOURCE_LIMIT","global arrangement cell budget exceeded");
}
P integerPolynomial(const Coefficients& input) {
    Z den(1);
    for(const auto& q:input)mpz_lcm(den.mpz(),den.mpz(),q.denominator().mpz());
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);
    P result(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(input[5*i+j]!=0){
        const Z v=input[5*i+j].numerator()*(den/input[5*i+j].denominator());
        result+=P(CORE::BigInt(text(v)))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    return result;
}
P primitive(const P& p) {
    if(CGAL::is_zero(p))return p;
    Z gcd=0,leading=0;
    for(int j=0;j<=p.degree();++j)for(int i=0;i<=p[j].degree();++i){
        const Z c(text(p[j][i]));
        mpz_gcd(gcd.mpz(),gcd.mpz(),c.mpz());
        if(c!=0)leading=c;
    }
    if(leading<0)gcd=-gcd;
    return p/P(CORE::BigInt(text(gcd)));
}
std::string polynomialIdentity(const P& p) {
    const auto c=primitive(p);std::ostringstream out;out.imbue(std::locale::classic());
    for(int j=0;j<=c.degree();++j)for(int i=0;i<=c[j].degree();++i)
        if(c[j][i]!=0)out<<i<<','<<j<<':'<<c[j][i]<<';';
    return out.str();
}
Q evaluate(const Coefficients& p,const Q& x,const Q& y) {
    Q result=0;
    for(int i=4;i>=0;--i){
        Q row=0;for(int j=4;j>=0;--j)row=row*y+p[5*i+j];
        result=result*x+row;
    }
    return result;
}
CORE::Expr expression(const Q& q){return CORE::Expr(CORE::BigRat(
    CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator()))));}
CORE::Expr evaluate(const Coefficients& p,const CORE::Expr& x,const CORE::Expr& y) {
    CORE::Expr result(0);
    for(int i=4;i>=0;--i){
        CORE::Expr row(0);for(int j=4;j>=0;--j)row=row*y+expression(p[5*i+j]);
        result=result*x+row;
    }
    return result;
}
int sign(const Q& v){return v>0?1:(v<0?-1:0);}
Bound bound(const Q& v) {
    return Bound(CORE::BigInt(text(v.numerator())),CORE::BigInt(text(v.denominator())));
}
void uniqueRoots(std::vector<AR>& roots) {
    std::sort(roots.begin(),roots.end(),[](const AR& a,const AR& b){return a.compare(b)==CGAL::SMALLER;});
    roots.erase(std::unique(roots.begin(),roots.end(),[](const AR& a,const AR& b){
        return a.compare(b)==CGAL::EQUAL;}),roots.end());
}
Z floorScaledRoot(AR root,const Z& denominator,int bits) {
    // A floor is an exact value property. In contrast, bound_between and
    // isolator endpoints depend on how much an earlier overlay refined a
    // shared CGAL root, and must not choose our scientific sign witness.
    for(int step=0;step<=4*bits;++step){
        const Q l=Q(text(root.low()))*Q(denominator),h=Q(text(root.high()))*Q(denominator);
        Z lo,hi;
        mpz_fdiv_q(lo.mpz(),l.numerator().mpz(),l.denominator().mpz());
        mpz_fdiv_q(hi.mpz(),h.numerator().mpz(),h.denominator().mpz());
        if(lo==hi)return lo;
        if(hi==lo+1)
            return root.compare(bound(Q(hi,denominator)))==CGAL::SMALLER?lo:hi;
        root.refine();
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","cannot certify slab floor within isolation budget");
}
std::vector<Bound> gaps(std::vector<AR> roots,int bits) {
    uniqueRoots(roots);
    if(roots.empty())return {Bound(0)};
    std::vector<Bound> result{bound(Q(floorScaledRoot(roots.front(),Z(1),bits)-1))};
    for(std::size_t i=1;i<roots.size();++i){
        Z denominator=1;bool found=false;
        for(int precision=0;precision<=bits;++precision,denominator*=2){
            const auto candidate=bound(Q(floorScaledRoot(roots[i-1],denominator,bits)+1,denominator));
            if(roots[i].compare(candidate)==CGAL::LARGER){
                result.push_back(candidate);found=true;break;
            }
        }
        if(!found)throw RankTwoSourceError("RESOURCE_LIMIT","dyadic slab witness budget exceeded");
    }
    result.push_back(bound(Q(floorScaledRoot(roots.back(),Z(1),bits)+1)));
    return result;
}
using Cell=RankTwoRegionKind;
struct Face {
    bool certified=false,unbounded=false;
    Q x,y;
    std::array<int,3> signs{};
    std::array<int,3> coneSigns{};
    Cell cell=Cell::Outside;
    std::string id;
};
std::string rootRecord(const AR& root,int bits,std::set<std::string>& proofs) {
    const auto real=detail::CertifiedAlgebraicReal::fromRoot(root,bits);
    const auto& c=real.certificate();std::ostringstream out;
    proofs.insert("ROOT_PROOF "+c.rootIdentity+" "+c.lower+" "+c.upper);
    out<<c.rootIdentity<<' '<<c.rootIndex<<' ';
    for(const auto& v:c.coefficients)out<<v<<',';
    return out.str();
}
std::string rootRecord(const detail::AlgebraicRootCertificate& c,std::set<std::string>& proofs) {
    proofs.insert("ROOT_PROOF "+c.rootIdentity+" "+c.lower+" "+c.upper);
    std::string record=c.rootIdentity+" "+std::to_string(c.rootIndex)+" ";
    for(const auto& v:c.coefficients)record+=v+",";
    return record;
}
std::string pointRecord(const Traits::Point_2& p,int bits,std::set<std::string>& proofs) {
    // Point_2::y() in this pinned release returns a reference to a temporary.
    // xy().y() is the correct value-returning coordinate accessor.
    const AR y=p.xy().y();
    return rootRecord(p.x(),bits,proofs)+" | "+rootRecord(y,bits,proofs);
}
std::size_t componentCount(const Arrangement& a,const std::vector<Face>& faces,Cell kind) {
    std::vector<std::size_t> parent(faces.size());std::iota(parent.begin(),parent.end(),0);
    auto root=[&](std::size_t n){while(parent[n]!=n)n=parent[n];return n;};
    for(auto edge=a.edges_begin();edge!=a.edges_end();++edge){
        const auto left=edge->face()->data(),right=edge->twin()->face()->data();
        if(faces[left].cell==kind && faces[right].cell==kind)parent[root(left)]=root(right);
    }
    std::set<std::size_t> components;
    for(std::size_t i=0;i<faces.size();++i)if(faces[i].cell==kind)components.insert(root(i));
    return components.size();
}
std::size_t boundaryLoops(const Arrangement& a,const std::vector<Face>& faces,Cell kind) {
    using H=Arrangement::Halfedge_const_handle;
    std::set<const void*> visited;std::size_t count=0;
    for(auto edge=a.halfedges_begin();edge!=a.halfedges_end();++edge){
        if(faces[edge->face()->data()].cell!=kind || faces[edge->twin()->face()->data()].cell==kind)
            continue;
        H start=edge;
        if(visited.count(&*start))continue;
        ++count;H current=start;
        do {
            if(!visited.insert(&*current).second)unresolved("selected-region boundary is not a cycle");
            H next=current->next();std::size_t steps=0;
            while(faces[next->twin()->face()->data()].cell==kind){
                if(++steps>a.number_of_halfedges())unresolved("region boundary traversal did not close");
                next=next->twin()->next();
            }
            current=next;
        }while(current!=start);
    }
    return count;
}
} // namespace

struct RankTwoArrangement::Impl {
    // Own traits longer than its DCEL. No native handles enter serialized IDs.
    Traits traits;
    // Arrangement copies retain the original traits pointer in pinned CGAL.
    // The shared global owner therefore outlives every overlay and its DCEL.
    std::shared_ptr<const Impl> parent;
    Arrangement arrangement;
    std::array<Coefficients,3> original{};
    std::vector<Traits::Curve_2> curves;
    std::vector<AR> events;
    std::vector<Coefficients> halfplanes;
    std::vector<std::string> normalRecords;
    std::optional<RankTwoNormalFace> normalFace;
    std::vector<RankTwoNormalBoundaryEvent> normalEvents;
    std::vector<RankTwoNormalVertexRecord> normalVertices;
    RankTwoArrangementRecords structured;
    // Factory values, not predicate outcomes. Each overlay builds its own map
    // from its actual DCEL; serialized reload still replays construction.
    std::map<std::string,std::array<detail::CertifiedAlgebraicReal,2>> certifiedVertices;
    std::shared_ptr<const RankTwoWholeArrangementRecords> whole;
    std::string inputHash;
    Q radius;
    std::vector<Face> faces;
    std::vector<RankTwoZeroFactor> factors;
    std::size_t cuts=0,cutLoops=0,fulls=0,isolated=0;
    std::string record,hash;
    Impl():arrangement(&traits){}
    explicit Impl(std::shared_ptr<const Impl> base):parent(std::move(base)),
        arrangement(parent->arrangement),original(parent->original),
        curves(parent->curves),events(parent->events),inputHash(parent->inputHash),
        radius(parent->radius),factors(parent->factors){}
    static void finalize(const std::shared_ptr<Impl>&,const PreparedRankTwoSource&,
                         const Q&,const RankTwoSourceLimits&);
    static void collectEvents(const std::shared_ptr<Impl>&,const PreparedRankTwoSource&,
        std::vector<std::string>&,std::map<const void*,std::string>&,
        const std::map<const void*,std::array<detail::CertifiedAlgebraicReal,2>>&);
    void requireLocal() const {
        if(whole)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","local DCEL accessor used on whole-grid arrangement");
    }
};
RankTwoArrangement::RankTwoArrangement(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
std::size_t RankTwoArrangement::faceCount()const{impl_->requireLocal();return impl_->faces.size();}
std::size_t RankTwoArrangement::certifiedFaceCount()const{
    impl_->requireLocal();
    return std::count_if(impl_->faces.begin(),impl_->faces.end(),[](const Face& f){return f.certified;});
}
std::size_t RankTwoArrangement::cutComponents()const{impl_->requireLocal();return impl_->cuts;}
std::size_t RankTwoArrangement::cutBoundaryLoops()const{impl_->requireLocal();return impl_->cutLoops;}
std::size_t RankTwoArrangement::fullComponents()const{impl_->requireLocal();return impl_->fulls;}
std::size_t RankTwoArrangement::isolatedPoints()const{impl_->requireLocal();return impl_->isolated;}
const std::vector<RankTwoZeroFactor>& RankTwoArrangement::factors()const{impl_->requireLocal();return impl_->factors;}
const std::string& RankTwoArrangement::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoArrangement::scientificHash()const{return impl_->hash;}
const std::string& RankTwoArrangement::sourceIdentity()const{impl_->requireLocal();return impl_->inputHash;}
const CGAL::Gmpq& RankTwoArrangement::truncationRadius()const{impl_->requireLocal();return impl_->radius;}
const std::vector<RankTwoNormalBoundaryEvent>& normalBoundaryEvents(const RankTwoArrangement& a){
    a.impl_->requireLocal();
    return a.impl_->normalEvents;
}
const std::vector<RankTwoNormalVertexRecord>& normalVertexRecords(const RankTwoArrangement& a){
    a.impl_->requireLocal();
    return a.impl_->normalVertices;
}
const RankTwoArrangementRecords& arrangementRecords(const RankTwoArrangement& a){a.impl_->requireLocal();return a.impl_->structured;}
const std::array<detail::CertifiedAlgebraicReal,2>& certifiedRankTwoVertex(
    const RankTwoArrangement& a,const std::string& id){
    detail::requireExactAlgebraicEnvironment();
    a.impl_->requireLocal();
    const auto found=a.impl_->certifiedVertices.find(id);
    if(found==a.impl_->certifiedVertices.end())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","native event absent from current local arrangement");
    return found->second;
}
const RankTwoWholeArrangementRecords& wholeArrangementRecords(const RankTwoArrangement& a){
    if(!a.impl_->whole)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","whole-grid accessor used on local DCEL");
    return *a.impl_->whole;
}
RankTwoArrangement RankTwoArrangement::reload(const PreparedRankTwoSource& source,const Q& radius,
    const std::string& saved,const RankTwoSourceLimits& limits) {
    const std::string header="rank-two-global-arrangement-v2\nINPUT "+source.inputHash()
        +"\nRADIUS "+fraction(radius)+"\n";
    if(saved.compare(0,header.size(),header)!=0)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","partition authority/radius/schema differs");
    // This early global adapter replays the exact construction rather than
    // trusting counts, hashes or a caller's acceptance field. Later carrier
    // reload consumes the structured, checked construction attachment.
    auto result=analyzeRankTwoGlobalArrangement(source,radius,limits);
    if(result.canonicalRecord()!=saved)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","saved global DCEL/sign evidence differs");
    return result;
}

RankTwoArrangement analyzeRankTwoGlobalArrangement(const PreparedRankTwoSource& source,
    const Q& radius,const RankTwoSourceLimits& limits) {
    detail::requireExactAlgebraicEnvironment();
    if(radius<=0 || limits.maximumArrangementCells==0 || limits.maximumIsolationBits<1
        || limits.maximumIsolationBits>512)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","invalid global partition settings");
    auto out=std::make_shared<RankTwoArrangement::Impl>();
    out->inputHash=source.inputHash();out->radius=radius;
    auto& arr=out->arrangement;
    const auto& graph=source.graph();
    std::array<Coefficients,3> original{};
    original[0]=graph.spherePolynomial(radius);
    const int powers[6]={0,5,1,10,6,2};
    for(int i=0;i<6;++i)original[1][powers[i]]=graph.height[i];
    original[2][0]=radius*radius;
    for(int i=0;i<3;++i){
        original[2][10]-=graph.plane[i][0]*graph.plane[i][0];
        original[2][6]-=2*graph.plane[i][0]*graph.plane[i][1];
        original[2][2]-=graph.plane[i][1]*graph.plane[i][1];
    }
    const char* names[3]={"F","h","D"};
    std::vector<P> zeroSets;
    std::array<P,3> polys;
    for(int i=0;i<3;++i){
        polys[i]=integerPolynomial(original[i]);
        std::vector<std::pair<P,int>> factors;
        PT::Square_free_factorize_up_to_constant_factor{}(polys[i],std::back_inserter(factors));
        for(const auto& f:factors)if(CGAL::total_degree(f.first)>0)zeroSets.push_back(primitive(f.first));
    }
    // Split shared factors before insertion. Every final factor retains its
    // multiplicity in each ORIGINAL polynomial, not the square-free sign.
    bool changed=true;
    while(changed){
        changed=false;
        std::map<std::string,P> unique;
        for(const auto& p:zeroSets)unique.emplace(polynomialIdentity(p),p);
        zeroSets.clear();for(const auto& entry:unique)zeroSets.push_back(entry.second);
        for(std::size_t i=0;i<zeroSets.size()&&!changed;++i)for(std::size_t j=i+1;j<zeroSets.size();++j){
            const P g=primitive(PT::Gcd_up_to_constant_factor{}(zeroSets[i],zeroSets[j]));
            if(CGAL::total_degree(g)==0)continue;
            const P left=primitive(PT::Integral_division_up_to_constant_factor{}(zeroSets[i],g));
            const P right=primitive(PT::Integral_division_up_to_constant_factor{}(zeroSets[j],g));
            zeroSets.erase(zeroSets.begin()+j);zeroSets.erase(zeroSets.begin()+i);
            zeroSets.push_back(g);
            if(CGAL::total_degree(left)>0)zeroSets.push_back(left);
            if(CGAL::total_degree(right)>0)zeroSets.push_back(right);
            changed=true;break;
        }
    }
    std::vector<Traits::Curve_2> curves;
    std::vector<AR> events;
    for(const auto& p:zeroSets){
        RankTwoZeroFactor factor;factor.polynomial=polynomialIdentity(p);
        for(int i=0;i<3;++i){
            P remaining=polys[i];int multiplicity=0;
            while(CGAL::total_degree(remaining)>0){
                const P g=primitive(PT::Gcd_up_to_constant_factor{}(remaining,p));
                if(g!=p)break;
                remaining=PT::Integral_division_up_to_constant_factor{}(remaining,p);++multiplicity;
            }
            if(multiplicity)factor.owners.push_back({names[i],multiplicity});
        }
        out->factors.push_back(std::move(factor));
        curves.push_back(out->traits.construct_curve_2_object()(p));
        auto& curve=curves.back();
        CGAL::insert(arr,curve);checkBudget(arr,limits);
        for(int i=0;i<curve.number_of_status_lines_with_event();++i)
            events.push_back(curve.status_line_at_event(i).x());
    }
    out->original=original;out->curves=curves;out->events=events;
    RankTwoArrangement::Impl::finalize(out,source,radius,limits);
    return RankTwoArrangement(std::move(out));
}

void RankTwoArrangement::Impl::finalize(const std::shared_ptr<Impl>& out,
    const PreparedRankTwoSource& source,const Q& radius,const RankTwoSourceLimits& limits) {
    detail::rankTwoDiagnostic("finalize-faces:begin");
    auto& arr=out->arrangement;
    const auto& original=out->original;
    const auto& curves=out->curves;
    auto events=out->events;
    const char* names[3]={"F","h","D"};
    checkBudget(arr,limits);
    out->structured.originalFhD=out->original;
    out->structured.normalHalfplanes=out->halfplanes;
    std::size_t faceIndex=0;
    for(auto face=arr.faces_begin();face!=arr.faces_end();++face){
        face->set_data(faceIndex++);out->faces.push_back(Face{});
        out->faces.back().unbounded=face->is_unbounded();
    }
    for(auto vertex=arr.vertices_begin();vertex!=arr.vertices_end();++vertex){
        if(!vertex->is_at_open_boundary())events.push_back(vertex->point().x());
        if(vertex->is_isolated())++out->isolated;
    }
    CGAL::Arr_naive_point_location<Arrangement> locate(arr);
    // The event set contains all critical, intersection, isolated and vertical
    // asymptote abscissae. One exact rational slice per complementary interval
    // intersects every two-dimensional cell. This is not a sampling grid.
    for(const auto& x:gaps(events,limits.maximumIsolationBits)){
        std::vector<AR> ordinates;
        for(const auto& curve:curves){
            const auto line=curve.status_line_at_exact_x(x);
            if(line.covers_line())unresolved("slab abscissa lies on a vertical zero component");
            for(int j=0;j<line.number_of_events();++j)
                ordinates.push_back(line.algebraic_real_2(j).y());
        }
        for(const auto& y:gaps(ordinates,limits.maximumIsolationBits)){
            const auto located=locate.locate(Traits::Point_2(x,y));
            const auto* face=CGAL::Arr_point_location_result<Arrangement>{}
                .template assign<Arrangement::Face_const_handle>(located);
            if(!face)unresolved("slab witness lies on a zero stratum");
            auto& f=out->faces[(*face)->data()];
            if(f.certified)continue;
            f.x=Q(text(x));f.y=Q(text(y));
            for(int i=0;i<3;++i)f.signs[i]=sign(evaluate(original[i],f.x,f.y));
            if(!f.signs[0]||!f.signs[1]||!f.signs[2])unresolved("open witness has a zero original sign");
            f.cell=f.signs[2]<0?Cell::Outside:(f.signs[0]<0?Cell::Cut:
                (f.signs[1]<0?Cell::Full:Cell::Empty));
            for(std::size_t i=0;i<out->halfplanes.size();++i){
                f.coneSigns[i]=sign(evaluate(out->halfplanes[i],f.x,f.y));
                if(!f.coneSigns[i])unresolved("open witness has zero normal halfplane sign");
                if(f.coneSigns[i]<0)f.cell=Cell::Outside;
            }
            f.certified=true;
            f.id=sha("rank-two-face-v1;"+fraction(f.x)+";"+fraction(f.y));
        }
    }
    for(const auto& face:out->faces)if(!face.certified)unresolved("DCEL face lacks certified slab witness");
    out->cuts=componentCount(arr,out->faces,Cell::Cut);
    out->fulls=componentCount(arr,out->faces,Cell::Full);
    out->cutLoops=boundaryLoops(arr,out->faces,Cell::Cut);
    detail::rankTwoDiagnostic("finalize-faces:end");

    std::vector<std::string> records=out->normalRecords;
    std::set<std::string> proofs;
    for(int i=0;i<3;++i){
        std::string r=std::string("ORIGINAL ")+names[i]+" ";
        for(const auto& c:original[i])r+=fraction(c)+";";
        records.push_back(r);
    }
    for(const auto& f:out->factors){
        std::string r="FACTOR "+f.polynomial+" OWNERS ";
        for(const auto& owner:f.owners)r+=owner.polynomial+":"+std::to_string(owner.multiplicity)+";";
        records.push_back(r);
    }
    // Read each DCEL point once: every incident arc references the same frozen
    // root certificate, even if later comparisons simplify a CGAL root handle.
    std::map<const void*,std::string> vertexIds;
    std::map<const void*,std::string> pointRecords;
    std::map<const void*,RankTwoVertexRecord> pointData;
    // Event collection and the immutable local owner consume these same
    // certified coordinates, never a cross-face cache or a second root chosen
    // after a native handle has been simplified.
    std::map<const void*,std::array<detail::CertifiedAlgebraicReal,2>> certifiedPoints;
    detail::rankTwoDiagnostic("finalize-vertices:begin");
    for(auto vertex=arr.vertices_begin();vertex!=arr.vertices_end();++vertex){
        if(vertex->is_at_open_boundary())continue;
        RankTwoVertexRecord data;
        const auto& coordinate=certifiedPoints.emplace(&*vertex,
            std::array<detail::CertifiedAlgebraicReal,2>{
                detail::CertifiedAlgebraicReal::fromRoot(vertex->point().x(),limits.maximumIsolationBits),
                detail::CertifiedAlgebraicReal::fromRoot(vertex->point().xy().y(),limits.maximumIsolationBits)}).first->second;
        data.coordinates={coordinate[0].certificate(),coordinate[1].certificate()};
        data.isolated=vertex->is_isolated();
        const auto record=rootRecord(data.coordinates[0],proofs)+" | "+rootRecord(data.coordinates[1],proofs);
        const auto id=sha("rank-two-global-vertex-v1;"+record);
        vertexIds[&*vertex]=id;
        pointRecords[&*vertex]=record;
        pointData.emplace(&*vertex,std::move(data));
    }
    detail::rankTwoDiagnostic("finalize-vertices:end");
    detail::rankTwoDiagnostic("finalize-events:begin");
    collectEvents(out,source,records,vertexIds,certifiedPoints);
    detail::rankTwoDiagnostic("finalize-events:end");
    for(const auto& entry:pointRecords){
        records.push_back("VERTEX "+vertexIds.at(entry.first)+" "+entry.second);
        auto data=pointData.at(entry.first);data.identity=vertexIds.at(entry.first);
        out->certifiedVertices.emplace(data.identity,certifiedPoints.at(entry.first));
        out->structured.vertices.push_back(std::move(data));
    }
    std::map<const void*,std::string> edgeIds;
    detail::rankTwoDiagnostic("finalize-arcs:begin");
    for(auto edge=arr.edges_begin();edge!=arr.edges_end();++edge){
        const auto& arc=edge->curve();
        const bool forward=edge->direction()==CGAL::ARR_LEFT_TO_RIGHT;
        RankTwoArcRecord data;
        data.supportPolynomial=polynomialIdentity(arc.curve().polynomial_2());
        data.vertical=arc.is_vertical();if(!data.vertical)data.branch=arc.arcno();
        std::string record="ARC "+polynomialIdentity(arc.curve().polynomial_2())+" ";
        record+=arc.is_vertical()?"VERTICAL":"BRANCH "+std::to_string(arc.arcno());
        for(const auto end:{CGAL::ARR_MIN_END,CGAL::ARR_MAX_END}){
            auto& endpoint=data.ends[end==CGAL::ARR_MIN_END?0:1];
            record+=" END "+std::to_string(static_cast<int>(arc.location(end)))+" ";
            if(arc.is_finite(end)){
                const bool minimum=end==CGAL::ARR_MIN_END;
                const auto vertex=(minimum==forward)?edge->source():edge->target();
                record+=vertexIds.at(&*vertex);
                endpoint.vertexIdentity=vertexIds.at(&*vertex);
            }
            else {
                switch(arc.location(end)){
                    case CGAL::ARR_LEFT_BOUNDARY:endpoint.location=RankTwoArcEndLocation::LeftInfinity;break;
                    case CGAL::ARR_RIGHT_BOUNDARY:endpoint.location=RankTwoArcEndLocation::RightInfinity;break;
                    case CGAL::ARR_TOP_BOUNDARY:endpoint.location=RankTwoArcEndLocation::TopInfinity;break;
                    case CGAL::ARR_BOTTOM_BOUNDARY:endpoint.location=RankTwoArcEndLocation::BottomInfinity;break;
                    default:unresolved("unclassified nonfinite arc end");
                }
                if(arc.location(end)==CGAL::ARR_TOP_BOUNDARY || arc.location(end)==CGAL::ARR_BOTTOM_BOUNDARY){
                    endpoint.finiteAbscissa=detail::CertifiedAlgebraicReal::fromRoot(
                        arc.curve_end_x(end),limits.maximumIsolationBits).certificate();
                    record+=rootRecord(*endpoint.finiteAbscissa,proofs);
                }
            }
        }
        const auto id=sha(record);
        data.identity=id;
        edgeIds[&*edge]=(forward?"+":"-")+id;
        edgeIds[&*edge->twin()]=(forward?"-":"+")+id;
        data.leftRegion=out->faces[(forward?edge:edge->twin())->face()->data()].id;
        data.rightRegion=out->faces[(forward?edge->twin():edge)->face()->data()].id;
        record+=" ID "+id+" LEFT "+data.leftRegion+" RIGHT "+data.rightRegion;
        // One certified INTERIOR point per actual arc identifies its owners,
        // even when its support polynomial contains several components.
        std::vector<AR> bounds;
        bool hasLower=false,hasUpper=false;
        for(const auto end:{CGAL::ARR_MIN_END,CGAL::ARR_MAX_END}){
            const bool finiteBound=arc.is_vertical()?arc.is_finite(end):
                (arc.location(end)!=CGAL::ARR_LEFT_BOUNDARY && arc.location(end)!=CGAL::ARR_RIGHT_BOUNDARY);
            if(!finiteBound)continue;
            if(end==CGAL::ARR_MIN_END)hasLower=true;else hasUpper=true;
            if(arc.is_vertical()){
                const auto v=((end==CGAL::ARR_MIN_END)==forward)?edge->source():edge->target();
                bounds.push_back(v->point().xy().y());
            }else bounds.push_back(arc.curve_end_x(end));
        }
        const auto choices=gaps(bounds,limits.maximumIsolationBits);
        const Bound coordinate=hasLower?(hasUpper?choices.at(1):choices.back()):choices.front();
        AR x,y;
        if(arc.is_vertical()){x=arc.x();y=AR(coordinate);}
        else {
            x=AR(coordinate);
            if(!arc.is_in_x_range_interior(x))unresolved("arc witness abscissa is not interior");
            const auto section=arc.curve().status_line_at_exact_x(coordinate);
            if(section.covers_line())unresolved("arc witness has a vertical support component");
            y=section.algebraic_real_2(arc.arcno(x)).y();
        }
        const auto X=detail::CertifiedAlgebraicReal::fromRoot(x,limits.maximumIsolationBits);
        const auto Y=detail::CertifiedAlgebraicReal::fromRoot(y,limits.maximumIsolationBits);
        data.interiorWitness={X.certificate(),Y.certificate()};
        auto atPoint=[&](const P& p){
            CORE::Expr value(0);
            for(int j=p.degree();j>=0;--j){
                CORE::Expr row(0);for(int i=p[j].degree();i>=0;--i)row=row*X.expression()+CORE::Expr(p[j][i]);
                value=value*Y.expression()+row;
            }
            return value.sign();
        };
        const P support=primitive(arc.curve().polynomial_2());
        if(atPoint(support)!=0)unresolved("arc witness not on actual support");
        auto addOwner=[&](const Coefficients& c,const std::string& owner){
            P remaining=integerPolynomial(c);int multiplicity=0;
            while(CGAL::total_degree(remaining)>0 && atPoint(remaining)==0){
                const P factor=primitive(PT::Gcd_up_to_constant_factor{}(remaining,support));
                if(CGAL::total_degree(factor)==0)unresolved("arc owner has no common support factor");
                remaining=PT::Integral_division_up_to_constant_factor{}(remaining,factor);++multiplicity;
            }
            if(multiplicity)data.owners.push_back({owner,multiplicity});
        };
        for(int i=0;i<3;++i)addOwner(original[i],names[i]);
        if(out->normalFace)for(int i=0;i<3;++i){
            const auto a=out->normalFace->vertexIds[(i+1)%3],b=out->normalFace->vertexIds[(i+2)%3];
            addOwner(out->halfplanes[i],"NORMAL_EDGE:"+std::to_string(std::min(a,b))+":"+std::to_string(std::max(a,b)));
        }
        if(data.owners.empty())unresolved("arc lost every zero-set owner");
        record+=" OWNERS ";for(const auto& owner:data.owners)record+=owner.polynomial+":"+std::to_string(owner.multiplicity)+";";
        record+=" WITNESS "+rootRecord(data.interiorWitness[0],proofs)+" | "+rootRecord(data.interiorWitness[1],proofs);
        records.push_back(record);
        out->structured.arcs.push_back(std::move(data));
    }
    detail::rankTwoDiagnostic("finalize-arcs:end");
    detail::rankTwoDiagnostic("finalize-records:begin");
    for(auto face=arr.faces_begin();face!=arr.faces_end();++face){
        const auto& f=out->faces[face->data()];
        RankTwoRegionRecord data;
        data.identity=f.id;data.kind=f.cell;data.witness={f.x,f.y};
        data.originalSigns=f.signs;data.unbounded=f.unbounded;
        if(!out->halfplanes.empty())data.normalSigns.assign(f.coneSigns.begin(),f.coneSigns.end());
        std::string record="FACE "+f.id+" SAMPLE "+fraction(f.x)+" "+fraction(f.y)
            +" TYPE "+std::to_string(static_cast<int>(f.cell))+" SIGNS ";
        for(int s:f.signs)record+=std::to_string(s)+",";
        if(!out->halfplanes.empty()){
            record+=" CONE_SIGNS ";for(int s:f.coneSigns)record+=std::to_string(s)+",";
        }
        record+=f.unbounded?" UNBOUNDED":" BOUNDED";
        std::vector<std::string> loops;
        auto appendLoop=[&](auto begin,const char* kind){
            std::vector<std::string> edges;auto current=begin;
            do {
                edges.push_back(current->is_fictitious()?"AT_INFINITY":edgeIds.at(&*current));++current;
            }while(current!=begin);
            // Cycles have no distinguished start; lexicographically minimal
            // rotation canonicalizes them without changing their orientation.
            std::vector<std::string> best=edges;
            for(std::size_t i=1;i<edges.size();++i){
                std::rotate(edges.begin(),edges.begin()+1,edges.end());
                if(edges<best)best=edges;
            }
            std::string loop=kind;for(const auto& e:best)loop+=" "+e;loops.push_back(loop);
            std::vector<RankTwoDirectedArcRecord> refs;
            for(const auto& e:best)refs.push_back(e=="AT_INFINITY"?
                RankTwoDirectedArcRecord{"",false,true}:RankTwoDirectedArcRecord{e.substr(1),e[0]=='-',false});
            (std::string(kind)=="OUTER"?data.outerLoops:data.innerLoops).push_back(std::move(refs));
        };
        for(auto c=face->outer_ccbs_begin();c!=face->outer_ccbs_end();++c)appendLoop(*c,"OUTER");
        for(auto c=face->inner_ccbs_begin();c!=face->inner_ccbs_end();++c)appendLoop(*c,"INNER");
        std::sort(loops.begin(),loops.end());for(const auto& loop:loops)record+=" LOOP "+loop;
        for(auto v=face->isolated_vertices_begin();v!=face->isolated_vertices_end();++v){
            records.push_back("ISOLATED "+f.id+" "+vertexIds.at(&*v));
            data.isolatedVertices.push_back(vertexIds.at(&*v));
        }
        records.push_back(record);
        out->structured.regions.push_back(std::move(data));
    }
    std::sort(records.begin(),records.end());
    auto byId=[](const auto& a,const auto& b){return a.identity<b.identity;};
    std::sort(out->structured.vertices.begin(),out->structured.vertices.end(),byId);
    std::sort(out->structured.arcs.begin(),out->structured.arcs.end(),byId);
    std::sort(out->structured.regions.begin(),out->structured.regions.end(),byId);
    out->record=(out->parent?"rank-two-normal-overlay-v2\n":"rank-two-global-arrangement-v2\n")
        +std::string("INPUT ")+source.inputHash()+"\nRADIUS "+fraction(radius)+"\n";
    if(out->parent)out->record+="GLOBAL "+out->parent->hash+"\n";
    for(const auto& r:records)out->record+=r+"\n";
    out->hash=sha(out->record);
    // Refining an isolator changes proof evidence, not scientific root IDs.
    for(const auto& proof:proofs)out->record+=proof+"\n";
    detail::requireExactAlgebraicEnvironment();
    detail::rankTwoDiagnostic("finalize-records:end");
}

RankTwoArrangement reloadRankTwoNormalFace(const PreparedRankTwoSource& source,
    const RankTwoArrangement& global,const RankTwoNormalFace& face,const std::string& saved,
    const RankTwoSourceLimits& limits) {
    const auto result=analyzeRankTwoNormalFace(source,global,face,limits);
    if(saved!=result.canonicalRecord())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","saved normal cone/topology/event evidence differs");
    return result;
}

RankTwoArrangement analyzeRankTwoNormalFace(const PreparedRankTwoSource& source,
    const RankTwoArrangement& global,const RankTwoNormalFace& input,
    const RankTwoSourceLimits& limits) {
    detail::rankTwoDiagnostic("overlay:begin");
    detail::requireExactAlgebraicEnvironment();
    if(global.impl_->whole || global.impl_->parent || global.impl_->inputHash!=source.inputHash()
        || limits.maximumArrangementCells==0 || limits.maximumIsolationBits<1
        || limits.maximumIsolationBits>512)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","invalid global authority/overlay settings");
    auto face=input;
    // Same face under any local vertex permutation has the same construction
    // order and evidence; each edge still retains its canonical endpoint IDs.
    for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)if(face.vertexIds[j]<face.vertexIds[i]){
        std::swap(face.vertexIds[j],face.vertexIds[i]);std::swap(face.rays[j],face.rays[i]);
    }
    if(face.vertexIds[0]==face.vertexIds[1] || face.vertexIds[1]==face.vertexIds[2])
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","duplicate normal vertex identity");
    auto cross=[](const auto& a,const auto& b){return std::array<Q,3>{
        a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};};
    auto dot=[](const auto& a,const auto& b){Q sum=0;for(int i=0;i<3;++i)sum+=a[i]*b[i];return sum;};
    const Q det=dot(face.rays[0],cross(face.rays[1],face.rays[2]));
    if(det==0)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","singular normal triangle");
    detail::rankTwoDiagnostic("overlay-copy:begin");
    auto out=std::make_shared<RankTwoArrangement::Impl>(global.impl_);
    detail::rankTwoDiagnostic("overlay-copy:end");
    out->normalFace=face;
    std::set<std::string> known;
    for(const auto& curve:out->curves)known.insert(polynomialIdentity(curve.polynomial_2()));
    std::string faceRecord="NORMAL_FACE";
    for(int i=0;i<3;++i){
        faceRecord+=" "+std::to_string(face.vertexIds[i])+":";
        for(const auto& q:face.rays[i])faceRecord+=fraction(q)+",";
    }
    out->normalRecords.push_back(faceRecord);
    detail::rankTwoDiagnostic("overlay-lines:begin");
    for(int i=0;i<3;++i){
        const int j=(i+1)%3,k=(i+2)%3;
        const auto cof=cross(face.rays[j],face.rays[k]);
        Coefficients line{};
        for(int d=0;d<3;++d)for(int power=0;power<3;++power)
            line[power==0?0:(power==1?5:1)]+=cof[d]*source.graph().gradient[d][power]/det;
        // Invertible N and regular G ensure no coordinate is identically zero.
        const P p=primitive(integerPolynomial(line));
        if(CGAL::is_zero(p))unresolved("normal cone coordinate identically zero");
        out->halfplanes.push_back(line);
        const auto lo=std::min(face.vertexIds[j],face.vertexIds[k]);
        const auto hi=std::max(face.vertexIds[j],face.vertexIds[k]);
        const std::string owner="NORMAL_EDGE:"+std::to_string(lo)+":"+std::to_string(hi);
        std::string r=owner+" HALFPLANE ";for(const auto& q:line)r+=fraction(q)+";";
        out->normalRecords.push_back(r);
        if(CGAL::total_degree(p)==0)continue; // constant signed cone constraint
        const std::string identity=polynomialIdentity(p);
        auto existing=std::find_if(out->factors.begin(),out->factors.end(),
            [&](const RankTwoZeroFactor& f){return f.polynomial==identity;});
        if(existing!=out->factors.end())existing->owners.push_back({owner,1});
        else out->factors.push_back({identity,{{owner,1}}});
        if(!known.insert(identity).second)continue;
        out->curves.push_back(out->traits.construct_curve_2_object()(p));
        auto& curve=out->curves.back();
        CGAL::insert(out->arrangement,curve);checkBudget(out->arrangement,limits);
        for(int e=0;e<curve.number_of_status_lines_with_event();++e)
            out->events.push_back(curve.status_line_at_event(e).x());
    }
    detail::rankTwoDiagnostic("overlay-lines:end");
    RankTwoArrangement::Impl::finalize(out,source,out->radius,limits);
    detail::requireExactAlgebraicEnvironment();
    detail::rankTwoDiagnostic("overlay:end");
    return RankTwoArrangement(std::move(out));
}

void RankTwoArrangement::Impl::collectEvents(const std::shared_ptr<Impl>& out,
    const PreparedRankTwoSource& source,std::vector<std::string>& records,
    std::map<const void*,std::string>& vertexIds,
    const std::map<const void*,std::array<detail::CertifiedAlgebraicReal,2>>& certifiedPoints) {
    if(!out->normalFace)return;
    const auto& face=*out->normalFace;
    struct Candidate { RankTwoNormalBoundaryEvent event;const void* vertex; };
    struct Point {
        const void* vertex;
        Traits::Point_2 actual;
        std::array<detail::CertifiedAlgebraicReal,2> coordinate;
    };
    std::vector<Point> points;
    detail::rankTwoDiagnostic("events-certificates:begin");
    for(auto v=out->arrangement.vertices_begin();v!=out->arrangement.vertices_end();++v){
        if(v->is_at_open_boundary())continue;
        points.push_back({&*v,v->point(),certifiedPoints.at(&*v)});
    }
    detail::rankTwoDiagnostic("events-certificates:end");
    detail::rankTwoDiagnostic("events-normal-vertices:begin");
    for(int index=0;index<3;++index){
        RankTwoNormalVertexRecord v;v.normalVertex=face.vertexIds[index];
        const auto& ray=face.rays[index];const auto& g=source.graph();
        Q projection=0;
        for(int d=0;d<3;++d)projection+=g.kernel[d]*ray[d];
        std::string key=source.geometryHash()+";R="+fraction(out->radius)+";"+std::to_string(v.normalVertex)+":";
        for(const auto& q:ray)key+=fraction(q)+",";
        v.identity=sha("rank-two-normal-vertex-v1;"+key);
        if(projection>0){
            v.location=RankTwoNormalVertexLocation::FinitePositive;
            // u.G=beta fixes the positive scale of the ray. G(v) is an
            // injective affine map of the graph plane: solve a nonsingular
            // 2x2 minor, then verify ALL THREE equations exactly.
            const Q lambda=g.beta/projection;
            std::array<Q,3> rhs{};
            for(int d=0;d<3;++d)rhs[d]=lambda*ray[d]-g.gradient[d][0];
            std::array<Q,2> p{};bool found=false;
            for(int i=0;i<3 && !found;++i)for(int j=i+1;j<3;++j){
                const Q det=g.gradient[i][1]*g.gradient[j][2]-g.gradient[j][1]*g.gradient[i][2];
                if(det==0)continue;
                p[0]=(rhs[i]*g.gradient[j][2]-rhs[j]*g.gradient[i][2])/det;
                p[1]=(g.gradient[i][1]*rhs[j]-g.gradient[j][1]*rhs[i])/det;
                found=true;break;
            }
            if(!found)unresolved("normal vertex has singular graph gradient map");
            for(int d=0;d<3;++d)if(g.gradient[d][1]*p[0]+g.gradient[d][2]*p[1]!=rhs[d])
                unresolved("normal vertex preimage failed gradient identity");
            v.point=p;
            for(int z=0;z<3;++z)v.originalSigns[z]=sign(evaluate(out->original[z],p[0],p[1]));
            std::size_t matches=0;
            for(const auto& actual:points)
                if(actual.coordinate[0].expression()==expression(p[0])
                    && actual.coordinate[1].expression()==expression(p[1])){
                    vertexIds.at(actual.vertex)=v.identity;++matches;
                }
            if(matches!=1)unresolved("finite normal vertex missing or duplicated in actual DCEL");
        }else v.location=projection==0?RankTwoNormalVertexLocation::Ideal:RankTwoNormalVertexLocation::Negative;
        std::string row="NORMAL_VERTEX "+v.identity+" INDEX "+std::to_string(v.normalVertex)+" LOCATION ";
        row+=v.location==RankTwoNormalVertexLocation::FinitePositive?"FINITE_POSITIVE":
             (v.location==RankTwoNormalVertexLocation::Ideal?"IDEAL":"NEGATIVE");
        if(v.point){
            row+=" POINT "+fraction((*v.point)[0])+" "+fraction((*v.point)[1])+" SIGNS ";
            for(int s:v.originalSigns)row+=std::to_string(s)+",";
        }
        records.push_back(row);out->normalVertices.push_back(std::move(v));
    }
    detail::rankTwoDiagnostic("events-normal-vertices:end");
    detail::rankTwoDiagnostic("events-normal-edges:begin");
    for(int low=0;low<3;++low)for(int high=low+1;high<3;++high){
        const auto& n=face.rays[low];const auto& m=face.rays[high];
        int i=-1,j=-1;Q det;
        for(int a=0;a<3 && i<0;++a)for(int b=a+1;b<3;++b){
            det=n[a]*m[b]-n[b]*m[a];if(det!=0){i=a;j=b;break;}
        }
        if(i<0)unresolved("dependent normal edge rays");
        // With independent n,m, g lies in span(n,m) iff (n x m).g=0.
        // Collapse this affine polynomial over Q BEFORE evaluating algebraic
        // coordinates. The two Cramer equations below are identities; asking
        // CORE to re-prove their expanded cancellations at every point builds
        // unnecessarily expensive root-separation expressions.
        std::array<Q,3> membership{};
        for(int d=0;d<3;++d){
            const Q cross=n[(d+1)%3]*m[(d+2)%3]-n[(d+2)%3]*m[(d+1)%3];
            for(int power=0;power<3;++power)
                membership[power]+=cross*source.graph().gradient[d][power];
        }
        Coefficients membershipCoefficients{};
        membershipCoefficients[0]=membership[0];membershipCoefficients[5]=membership[1];
        membershipCoefficients[1]=membership[2];
        const P membershipPolynomial=integerPolynomial(membershipCoefficients);
        const bool constantMembership=membership[1]==0 && membership[2]==0;
        if(constantMembership){
            if(membership[0]==0)unresolved("normal edge contains the entire graph gradient plane");
            continue;
        }
        // A common edge has one free coordinate. Use its identical rational
        // affine map on both incident faces, instead of a quotient of two
        // independently reconstructed (but algebraically dependent) roots.
        // The exact parameter/point equality checks below remain mandatory.
        const int freeCoordinate=membership[2]!=0?0:1;
        const int dependent=1-freeCoordinate;
        std::array<Q,2> offset{},direction{};direction[freeCoordinate]=1;
        offset[dependent]=-membership[0]/membership[dependent+1];
        direction[dependent]=-membership[freeCoordinate+1]/membership[dependent+1];
        std::array<std::array<Q,2>,3> along{};
        for(int d=0;d<3;++d){
            const auto& g=source.graph().gradient[d];
            along[d][0]=g[0]+g[1]*offset[0]+g[2]*offset[1];
            along[d][1]=g[1]*direction[0]+g[2]*direction[1];
        }
        std::array<Q,2> alpha{},beta{};
        for(int k=0;k<2;++k){
            alpha[k]=(along[i][k]*m[j]-along[j][k]*m[i])/det;
            beta[k]=(n[i]*along[j][k]-n[j]*along[i][k])/det;
            for(int d=0;d<3;++d)
                if(alpha[k]*n[d]+beta[k]*m[d]!=along[d][k])
                    unresolved("normal edge parameter failed exact span identity");
        }
        auto& algebraicKernel=Traits::CKvA_2::instance().kernel();
        std::vector<Candidate> candidates;
        for(const auto& point:points){
            detail::rankTwoDiagnostic("event-membership:begin");
            // Preserve the native algebraic relation between x and y. This
            // exact zero predicate has no finite-precision UNKNOWN-as-zero
            // overload. Nonconstant affine polynomials are square-free;
            // constants are handled before entering curve analysis.
            const bool onPlane=algebraicKernel.is_zero_at_2_object()(membershipPolynomial,point.actual.xy());
            detail::rankTwoDiagnostic("event-membership:end");
            if(!onPlane)continue;
            detail::rankTwoDiagnostic("event-positive:begin");
            const auto& coordinate=point.coordinate[freeCoordinate].expression();
            const CORE::Expr a=expression(alpha[0])+expression(alpha[1])*coordinate;
            const CORE::Expr b=expression(beta[0])+expression(beta[1])*coordinate;
            const bool positive=a.sign()>0 && b.sign()>0;
            detail::rankTwoDiagnostic("event-positive:end");
            if(!positive)continue;
            RankTwoNormalBoundaryEvent event;
            event.edge={face.vertexIds[low],face.vertexIds[high]};
            const char* names[3]={"F","h","D"};
            detail::rankTwoDiagnostic("event-original-zero:begin");
            for(int z=0;z<3;++z)if(evaluate(out->original[z],point.coordinate[0].expression(),
                    point.coordinate[1].expression()).sign()==0)event.zeroPolynomials.push_back(names[z]);
            detail::rankTwoDiagnostic("event-original-zero:end");
            if(event.zeroPolynomials.empty())continue;
            event.parameter=b/(a+b);
            for(int d=0;d<2;++d)event.point[d]=point.coordinate[d].certificate();
            candidates.push_back({std::move(event),point.vertex});
        }
        std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b){
            return a.event.parameter<b.event.parameter;});
        std::string edgeKey=source.geometryHash()+";R="+fraction(out->radius)+";";
        for(int index:{low,high}){
            edgeKey+=std::to_string(face.vertexIds[index])+":";
            for(const auto& q:face.rays[index])edgeKey+=fraction(q)+",";
        }
        std::ptrdiff_t radiusOrdinal=0;
        for(std::size_t index=0;index<candidates.size();++index){
            auto& c=candidates[index];auto& e=c.event;e.ordinal=index;
            if(index && e.parameter==candidates[index-1].event.parameter)
                unresolved("duplicate exact normal-edge event");
            if(std::find(e.zeroPolynomials.begin(),e.zeroPolynomials.end(),"F")!=e.zeroPolynomials.end())
                e.radiusOrdinal=radiusOrdinal++;
            std::string kind;
            for(const auto& z:e.zeroPolynomials)kind+=z+",";
            e.identity=sha("rank-two-normal-boundary-event-v1;"+edgeKey+";"+kind+";"+std::to_string(index));
            vertexIds.at(c.vertex)=e.identity;
            records.push_back("NORMAL_EVENT "+e.identity+" EDGE "+std::to_string(e.edge[0])+" "
                +std::to_string(e.edge[1])+" ORDINAL "+std::to_string(e.ordinal)
                +" RADIUS_ORDINAL "+std::to_string(e.radiusOrdinal)+" ZERO "+kind);
            out->normalEvents.push_back(std::move(e));
        }
    }
    detail::rankTwoDiagnostic("events-normal-edges:end");
}

namespace {
struct EdgeArcInterval {
    CORE::Expr lower,upper;
    std::string arc;
    bool reversed=false;
    std::array<int,3> signs{};
    std::vector<RankTwoPolynomialOwner> owners;
};
std::vector<EdgeArcInterval> edgeArcIntervals(const PreparedRankTwoSource& source,
    const PreparedRankTwoNormalMesh& mesh,const std::array<std::size_t,2>& edge,
    const RankTwoArrangement& overlay,const RankTwoSourceLimits& limits){
    const auto& n=mesh.rays()[edge[0]];const auto& m=mesh.rays()[edge[1]];
    const auto& graph=source.graph();const auto& records=arrangementRecords(overlay);
    const std::array<Q,3> normal={n[1]*m[2]-n[2]*m[1],n[2]*m[0]-n[0]*m[2],n[0]*m[1]-n[1]*m[0]};
    std::array<Q,3> line{};
    for(int d=0;d<3;++d)for(int p=0;p<3;++p)line[p]+=normal[d]*graph.gradient[d][p];
    const std::string owner="NORMAL_EDGE:"+std::to_string(edge[0])+":"+std::to_string(edge[1]);
    if(line[1]==0 && line[2]==0){
        if(line[0]==0)unresolved("normal edge contains the entire graph gradient plane");
        return {}; // the supporting plane never meets the affine Gauss image
    }
    const bool vertical=line[2]==0;
    // Parametrize the actual supporting line by x, or by y for a vertical
    // line. Clip in THIS coordinate before applying its fractional t map;
    // a+b can vanish outside the positive cone and must not reorder it.
    std::array<Q,2> offset{},direction{};
    if(vertical){offset[0]=-line[0]/line[1];direction[1]=1;}
    else {offset[1]=-line[0]/line[2];direction={1,-line[1]/line[2]};}
    std::array<std::array<Q,2>,3> G{};
    for(int d=0;d<3;++d){
        G[d][0]=graph.gradient[d][0]+graph.gradient[d][1]*offset[0]+graph.gradient[d][2]*offset[1];
        G[d][1]=graph.gradient[d][1]*direction[0]+graph.gradient[d][2]*direction[1];
    }
    std::array<Q,2> a{},b{};bool found=false;
    for(int i=0;i<3 && !found;++i)for(int j=i+1;j<3;++j){
        const Q det=n[i]*m[j]-n[j]*m[i];if(det==0)continue;
        for(int k=0;k<2;++k){a[k]=(G[i][k]*m[j]-G[j][k]*m[i])/det;b[k]=(n[i]*G[j][k]-n[j]*G[i][k])/det;}
        found=true;break;
    }
    if(!found)unresolved("dependent shared normal edge");
    for(int d=0;d<3;++d)for(int k=0;k<2;++k)
        if(a[k]*n[d]+b[k]*m[d]!=G[d][k])unresolved("shared normal line failed exact span identity");
    const Q derivative=b[1]*a[0]-b[0]*a[1];
    if(derivative==0)unresolved("shared normal parameter is constant");
    std::map<std::string,const RankTwoVertexRecord*> vertices;
    for(const auto& v:records.vertices)if(!vertices.emplace(v.identity,&v).second)
        unresolved("duplicate actual vertex identity in normal overlay");
    std::vector<EdgeArcInterval> result;
    for(const auto& arc:records.arcs){
        if(std::none_of(arc.owners.begin(),arc.owners.end(),[&](const auto& o){return o.polynomial==owner;}))continue;
        if(arc.vertical!=vertical)unresolved("normal owner disagrees with actual straight-arc direction");
        std::array<std::optional<CORE::Expr>,2> ends;
        for(int i=0;i<2;++i){
            const auto& e=arc.ends[i];
            if(e.vertexIdentity){
                const auto& v=*vertices.at(*e.vertexIdentity);
                ends[i]=detail::CertifiedAlgebraicReal::reload(v.coordinates[vertical?1:0],limits.maximumIsolationBits).expression();
            }else {
                const auto expected=vertical?(i?RankTwoArcEndLocation::TopInfinity:RankTwoArcEndLocation::BottomInfinity):
                    (i?RankTwoArcEndLocation::RightInfinity:RankTwoArcEndLocation::LeftInfinity);
                if(e.location!=expected)unresolved("unexpected ideal endpoint on normal line");
            }
        }
        bool empty=false;
        for(const auto& coefficient:{a,b}){
            if(coefficient[1]==0){if(coefficient[0]<=0)empty=true;continue;}
            const CORE::Expr root=expression(-coefficient[0]/coefficient[1]);
            const int side=coefficient[1]>0?0:1;
            if(!ends[side] || (side==0?root>*ends[0]:root<*ends[1]))ends[side]=root;
        }
        if(empty || (ends[0] && ends[1] && *ends[0]>=*ends[1]))continue;
        auto parameter=[&](const std::optional<CORE::Expr>& endpoint){
            if(endpoint){
                const CORE::Expr aa=expression(a[0])+expression(a[1])**endpoint;
                const CORE::Expr bb=expression(b[0])+expression(b[1])**endpoint;
                if(aa.sign()<0 || bb.sign()<0 || (aa+bb).sign()<=0)
                    unresolved("clipped normal interval endpoint outside positive closure");
                return bb/(aa+bb);
            }
            if(a[1]+b[1]==0)unresolved("unbounded positive interval has no finite normal parameter limit");
            return expression(b[1]/(a[1]+b[1]));
        };
        EdgeArcInterval part;part.lower=parameter(ends[0]);part.upper=parameter(ends[1]);
        part.reversed=derivative<0;if(part.reversed)std::swap(part.lower,part.upper);
        if(part.lower<CORE::Expr(0) || part.upper>CORE::Expr(1) || part.lower>=part.upper)
            unresolved("invalid shared positive normal parameter interval");
        part.arc=arc.identity;
        const auto X=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[0],limits.maximumIsolationBits);
        const auto Y=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[1],limits.maximumIsolationBits);
        const char* names[3]={"F","h","D"};
        for(int k=0;k<3;++k){
            part.signs[k]=evaluate(records.originalFhD[k],X.expression(),Y.expression()).sign();
            const auto match=std::find_if(arc.owners.begin(),arc.owners.end(),[&](const auto& o){return o.polynomial==names[k];});
            if((part.signs[k]==0)!=(match!=arc.owners.end()))unresolved("shared arc original sign/owner disagreement");
            if(match!=arc.owners.end())part.owners.push_back(*match);
        }
        result.push_back(std::move(part));
    }
    std::sort(result.begin(),result.end(),[](const auto& x,const auto& y){return x.lower<y.lower;});
    Q un=0,um=0;for(int d=0;d<3;++d){un+=graph.kernel[d]*n[d];um+=graph.kernel[d]*m[d];}
    if(un<=0 && um<=0){if(!result.empty())unresolved("ideal/negative normal edge has finite positive intervals");return {};}
    const CORE::Expr first=un<0?expression(un/(un-um)):CORE::Expr(0);
    const CORE::Expr last=um<0?expression(un/(un-um)):CORE::Expr(1);
    CORE::Expr position=first;
    for(const auto& part:result){
        if(part.lower!=position)unresolved("gap or overlap in actual shared normal edge cover");
        position=part.upper;
    }
    if(result.empty() || position!=last)unresolved("actual shared arcs do not cover positive Gauss edge");
    return result;
}

std::vector<RankTwoSharedInterval> compareSharedIntervals(const PreparedRankTwoSource& source,const Q& radius,
    const PreparedRankTwoNormalMesh& mesh,const RankTwoSharedEdgeRecord& edge,
    const std::vector<RankTwoArrangement>& faces,const RankTwoSourceLimits& limits){
    const auto left=edgeArcIntervals(source,mesh,edge.edge,faces[edge.incidentFaces[0]],limits);
    const auto right=edgeArcIntervals(source,mesh,edge.edge,faces[edge.incidentFaces[1]],limits);
    if(left.empty()!=right.empty())unresolved("one incident face lost the positive shared edge");
    if(left.empty())return {};
    std::size_t i=0,j=0;CORE::Expr current=left.front().lower;
    std::vector<RankTwoSharedInterval> result;
    while(i<left.size() && j<right.size()){
        const auto& a=left[i];const auto& b=right[j];
        if(a.lower>current || b.lower>current)unresolved("shared arc refinement has a gap");
        const CORE::Expr next=a.upper<b.upper?a.upper:b.upper;
        if(next<=current || a.signs!=b.signs || a.owners.size()!=b.owners.size())
            unresolved("incident shared arc signs/owners differ");
        for(std::size_t o=0;o<a.owners.size();++o)
            if(a.owners[o].polynomial!=b.owners[o].polynomial || a.owners[o].multiplicity!=b.owners[o].multiplicity)
                unresolved("incident shared arc original multiplicity differs");
        const std::string key=source.geometryHash()+";R="+fraction(radius)+";"+std::to_string(edge.edge[0])+":"
            +std::to_string(edge.edge[1])+";"+std::to_string(result.size());
        result.push_back({sha("rank-two-shared-interval-v1;"+key),current,next,{a.arc,b.arc},
            {a.reversed,b.reversed},a.signs,a.owners});
        current=next;if(a.upper==next)++i;if(b.upper==next)++j;
    }
    if(i!=left.size() || j!=right.size())unresolved("incident shared interval covers end differently");
    return result;
}
} // namespace

const std::vector<RankTwoSharedInterval>& sharedIntervals(const RankTwoSharedEdgeRecord& edge){return edge.intervals;}

RankTwoArrangement analyzeRankTwoArrangement(const PreparedRankTwoSource& source,const Q& radius,
    const PreparedRankTwoNormalMesh& mesh,const RankTwoSourceLimits& limits) {
    detail::requireExactAlgebraicEnvironment();
    if(mesh.rays().size()>limits.maximumVertices)
        throw RankTwoSourceError("RESOURCE_LIMIT","whole-grid normal vertex budget exceeded");
    auto global=analyzeRankTwoGlobalArrangement(source,radius,limits);
    auto data=std::make_shared<RankTwoWholeArrangementRecords>(RankTwoWholeArrangementRecords{global,mesh,{},{},{},0});
    auto consume=[&](const RankTwoArrangement& a){
        const auto& r=arrangementRecords(a);
        for(const auto count:{r.vertices.size(),r.arcs.size(),r.regions.size()}){
            if(count>limits.maximumArrangementCells-data->totalArrangementCells)
                throw RankTwoSourceError("RESOURCE_LIMIT","whole-grid aggregate arrangement cell budget exceeded");
            data->totalArrangementCells+=count;
        }
    };
    consume(global);
    std::map<std::array<std::size_t,2>,std::vector<std::size_t>> incident;
    std::map<std::size_t,RankTwoNormalVertexRecord> vertexRegistry;
    std::vector<std::size_t> vertexUse(mesh.rays().size(),0),expectedUse(mesh.rays().size(),0);
    for(std::size_t f=0;f<mesh.faces().size();++f){
        const auto& tri=mesh.faces()[f];
        for(int j=0;j<3;++j){
            ++expectedUse[tri[j]];
            const auto pair=std::minmax(tri[j],tri[(j+1)%3]);
            incident[{pair.first,pair.second}].push_back(f);
        }
        auto localLimits=limits;
        localLimits.maximumArrangementCells=limits.maximumArrangementCells-data->totalArrangementCells;
        if(!localLimits.maximumArrangementCells)
            throw RankTwoSourceError("RESOURCE_LIMIT","whole-grid aggregate arrangement cell budget exhausted");
        data->normalFaces.push_back(analyzeRankTwoNormalFace(source,global,mesh.face(f),localLimits));
        const auto& overlay=data->normalFaces.back();consume(overlay);
        const auto& vertices=normalVertexRecords(overlay);
        if(vertices.size()!=3)unresolved("normal face lost a vertex record");
        for(const auto& v:vertices){
            if(std::find(tri.begin(),tri.end(),v.normalVertex)==tri.end())
                unresolved("normal vertex record belongs to another face");
            ++vertexUse.at(v.normalVertex);
            const auto inserted=vertexRegistry.emplace(v.normalVertex,v);
            if(!inserted.second){
                const auto& old=inserted.first->second;
                if(old.identity!=v.identity || old.location!=v.location || old.point!=v.point
                    || old.originalSigns!=v.originalSigns)
                    unresolved("incident faces disagree on exact normal vertex");
            }
        }
    }
    if(vertexUse!=expectedUse || vertexRegistry.size()!=mesh.rays().size())
        unresolved("whole normal vertex registry lacks an incident face");
    for(const auto& entry:vertexRegistry)data->normalVertices.push_back(entry.second);
    if(incident.size()!=mesh.edges().size())unresolved("whole shared edge registry size mismatch");
    detail::rankTwoDiagnostic("atlas-shared-edges:begin");
    for(const auto& e:mesh.edges()){
        detail::rankTwoDiagnostic("atlas-shared-edge:begin");
        const auto& owners=incident.at(e);
        if(owners.size()!=2 || owners[0]==owners[1])unresolved("normal edge is not two-sided");
        std::array<std::vector<RankTwoNormalBoundaryEvent>,2> events;
        for(int side=0;side<2;++side)
            for(const auto& event:normalBoundaryEvents(data->normalFaces.at(owners[side])))
                if(event.edge==e)events[side].push_back(event);
        if(events[0].size()!=events[1].size())unresolved("incident faces disagree on shared edge event count");
        for(std::size_t k=0;k<events[0].size();++k){
            const auto& a=events[0][k];const auto& b=events[1][k];
            detail::rankTwoDiagnostic("shared-event-identity:begin");
            if(a.identity!=b.identity || a.ordinal!=b.ordinal || a.radiusOrdinal!=b.radiusOrdinal
                || a.zeroPolynomials!=b.zeroPolynomials || a.parameter!=b.parameter)
                unresolved("incident faces disagree on shared edge event identity");
            detail::rankTwoDiagnostic("shared-event-identity:end");
            detail::rankTwoDiagnostic("shared-event-point:begin");
            for(int d=0;d<2;++d){
                const auto x=detail::CertifiedAlgebraicReal::reload(a.point[d],limits.maximumIsolationBits);
                const auto y=detail::CertifiedAlgebraicReal::reload(b.point[d],limits.maximumIsolationBits);
                if(x.expression()!=y.expression())unresolved("incident faces disagree on exact shared point");
            }
            detail::rankTwoDiagnostic("shared-event-point:end");
        }
        RankTwoSharedEdgeRecord shared{e,{owners[0],owners[1]},std::move(events[0]),{}};
        detail::rankTwoDiagnostic("shared-intervals:begin");
        shared.intervals=compareSharedIntervals(source,radius,mesh,shared,data->normalFaces,limits);
        detail::rankTwoDiagnostic("shared-intervals:end");
        data->sharedEdges.push_back(std::move(shared));
        detail::rankTwoDiagnostic("atlas-shared-edge:end");
    }
    detail::rankTwoDiagnostic("atlas-shared-edges:end");
    auto result=std::make_shared<RankTwoArrangement::Impl>();result->whole=data;
    std::string record="rank-two-whole-arrangement-v2\nINPUT "+source.inputHash()+"\nRADIUS "+fraction(radius)
        +"\nNORMAL_MESH "+mesh.inputHash()+"\nGLOBAL "+global.scientificHash()+"\nCELLS "
        +std::to_string(data->totalArrangementCells)+"\n";
    for(std::size_t i=0;i<data->normalFaces.size();++i)
        record+="NORMAL_FACE "+std::to_string(i)+" "+data->normalFaces[i].scientificHash()+"\n";
    for(const auto& e:data->sharedEdges){
        record+="SHARED_EDGE "+std::to_string(e.edge[0])+" "+std::to_string(e.edge[1])+" FACES "
            +std::to_string(e.incidentFaces[0])+" "+std::to_string(e.incidentFaces[1])+" EVENTS";
        for(const auto& event:e.events)record+=" "+event.identity;
        record+='\n';
        for(const auto& interval:e.intervals){
            record+="SHARED_INTERVAL "+interval.identity+" ARCS "+interval.sourceArcs[0]+" "+interval.sourceArcs[1]
                +" REVERSED "+std::to_string(interval.sourceReversed[0])+" "+std::to_string(interval.sourceReversed[1])+" SIGNS";
            for(int sign:interval.originalSigns)record+=" "+std::to_string(sign);
            record+=" OWNERS";for(const auto& o:interval.originalOwners)record+=" "+o.polynomial+":"+std::to_string(o.multiplicity);
            record+='\n';
        }
    }
    for(const auto& v:data->normalVertices)
        record+="SHARED_VERTEX "+std::to_string(v.normalVertex)+" "+v.identity+"\n";
    result->hash=sha(record);result->record=record+"NORMAL_ADMISSION_BEGIN\n"+mesh.canonicalRecord()
        +"NORMAL_ADMISSION_END\nGLOBAL_BEGIN\n"+global.canonicalRecord()+"GLOBAL_END\n";
    for(std::size_t i=0;i<data->normalFaces.size();++i)
        result->record+="OVERLAY_BEGIN "+std::to_string(i)+"\n"+data->normalFaces[i].canonicalRecord()+"OVERLAY_END\n";
    detail::requireExactAlgebraicEnvironment();return RankTwoArrangement(std::move(result));
}

RankTwoArrangement reloadRankTwoArrangement(const PreparedRankTwoSource& source,const Q& radius,
    const PreparedRankTwoNormalMesh& mesh,const std::string& saved,const RankTwoSourceLimits& limits){
    const auto result=analyzeRankTwoArrangement(source,radius,mesh,limits);
    if(result.canonicalRecord()!=saved)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","whole-grid construction replay mismatch");
    return result;
}
} // namespace spin
