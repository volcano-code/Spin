#include "detail/RankTwoArcStrip.h"
#include "detail/RankTwoStripMemo.h"
#include "detail/RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <sstream>
#include <map>
#include <tuple>

namespace spin::detail {
namespace {
using Q=CGAL::Gmpq;
using P=RankTwoStripPolynomial;
struct Interval {Q lo,hi;};
Q power(Q x,int n){Q out(1);while(n-->0)out*=x;return out;}
int choose(int n,int k){
    if(k<0 || k>n)return 0;
    int out=1;for(int i=1;i<=k;++i)out=out*(n-i+1)/i;return out;
}
// Exact convex-hull enclosure for the COMPLETE horizontal edge. Mapping
// u=a+(b-a)t and degree-elevating to Bernstein degree four avoids expanded
// monomial cancellation without replacing the edge by corner samples.
Interval horizontalRange(const P& p,const Interval& x,const Q& y){
    std::array<Q,5> monomial{},shifted{},bernstein{};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)monomial[i]+=p[5*i+j]*power(y,j);
    for(int k=0;k<5;++k)for(int i=k;i<5;++i)
        shifted[k]+=monomial[i]*Q(choose(i,k))*power(x.lo,i-k)*power(x.hi-x.lo,k);
    for(int j=0;j<5;++j)for(int k=0;k<=j;++k)
        bernstein[j]+=shifted[k]*Q(choose(j,k),choose(4,k));
    return {*std::min_element(bernstein.begin(),bernstein.end()),
        *std::max_element(bernstein.begin(),bernstein.end())};
}
Interval power(const Interval& x,int n){
    if(n==0)return {Q(1),Q(1)};
    const Q a=power(x.lo,n),b=power(x.hi,n);
    if(n%2)return {a,b};
    return {x.lo<=0 && x.hi>=0?Q(0):std::min(a,b),std::max(a,b)};
}
Interval multiply(const Interval& a,const Interval& b){
    const std::array<Q,4> p{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    return {*std::min_element(p.begin(),p.end()),*std::max_element(p.begin(),p.end())};
}
Interval evaluate(const P& p,const Interval& x,const Interval& y){
    Interval out{Q(0),Q(0)};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0){
        const auto term=multiply({p[5*i+j],p[5*i+j]},multiply(power(x,i),power(y,j)));
        out.lo+=term.lo;out.hi+=term.hi;
    }
    return out;
}
// Exact tensor Bernstein convex hull over the ENTIRE rectangle. Used only
// when the monomial enclosure cannot certify obstacle exclusion. No corners,
// samples, floating approximation or tolerance can establish the sign.
Interval tensorBernsteinRange(const P& p,const Interval& x,const Interval& y){
    std::array<Q,25> shifted{},controls{};
    for(int k=0;k<5;++k)for(int l=0;l<5;++l)
        for(int i=k;i<5;++i)for(int j=l;j<5;++j)if(p[5*i+j]!=0)
            shifted[5*k+l]+=p[5*i+j]*Q(choose(i,k))*Q(choose(j,l))
                *power(x.lo,i-k)*power(y.lo,j-l)*power(x.hi-x.lo,k)*power(y.hi-y.lo,l);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)
        for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)
            controls[5*i+j]+=shifted[5*k+l]*Q(choose(i,k),choose(4,k))*Q(choose(j,l),choose(4,l));
    return {*std::min_element(controls.begin(),controls.end()),*std::max_element(controls.begin(),controls.end())};
}
P derivative(const P& p,int axis){
    P out{};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j){
        const int n=axis==0?i:j;
        if(n)out[5*(i-(axis==0))+(j-(axis==1))]+=Q(n)*p[5*i+j];
    }
    return out;
}
int strictSign(const Interval& v){return v.lo>0?1:v.hi<0?-1:0;}
Q absoluteUpper(const Interval& v){return std::max(CGAL::abs(v.lo),CGAL::abs(v.hi));}
RankTwoRegularStripCertificate requireCertificate(RankTwoRegularStripAttempt result){
    if(auto* proof=std::get_if<RankTwoRegularStripCertificate>(&result))return std::move(*proof);
    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",std::get<std::string>(result));
}
void validate(const P& p){
    bool nonzero=false;
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(p[5*i+j]!=0){
        nonzero=true;
        if(i+j>4)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","strip polynomial degree exceeds four");
    }
    if(!nonzero)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","zero strip polynomial");
}
}

P parseRankTwoStripPolynomial(const std::string& saved,int axis){
    if(axis!=0 && axis!=1)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","regular edge requires valid chart");
    P out{};std::istringstream input(saved);input.imbue(std::locale::classic());std::string term;
    while(std::getline(input,term,';')){
        const auto comma=term.find(','),colon=term.find(':');
        if(comma==std::string::npos || colon==std::string::npos || comma>=colon)
            throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","malformed immutable regular-edge support");
        const int i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        if(i<0 || j<0 || i+j>4)
            throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","regular-edge support degree outside contract");
        out[axis?5*j+i:5*i+j]+=Q(term.substr(colon+1));
    }
    return out;
}

RankTwoRegularStripAttempt tryCertifyRankTwoRegularStrip(const P& p,const RankTwoStripBox& box,
    const std::vector<P>& obstacles,const Q& maximumDistance){
    validate(p);for(const auto& other:obstacles)validate(other);
    if(box[0]>=box[1] || box[2]>=box[3] || maximumDistance<=0)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","ordered nonempty strip and positive chord error required");
    const Interval x{box[0],box[1]},y{box[2],box[3]};
    const P px=derivative(p,0),py=derivative(p,1);
    const auto dy=evaluate(py,x,y);const int sign=strictSign(dy);
    if(!sign)return std::string("vertical derivative not separated from zero on whole strip");
    // Prove existence over each COMPLETE horizontal edge, not just corners.
    // Fixed strict Py then proves uniqueness and the implicit C2 graph.
    if(strictSign(horizontalRange(p,x,y.lo))!=-sign
        ||strictSign(horizontalRange(p,x,y.hi))!=sign)
        return std::string("horizontal edges do not bracket one graph for every abscissa");
    RankTwoRegularStripCertificate out;out.yDerivativeSign=sign;
    out.derivativeLowerBound=sign>0?dy.lo:-dy.hi;
    const Q m=out.derivativeLowerBound;
    out.slopeUpperBound=absoluteUpper(evaluate(px,x,y))/m;
    const Q k=out.slopeUpperBound;
    out.curvatureUpperBound=(absoluteUpper(evaluate(derivative(px,0),x,y))
        +Q(2)*absoluteUpper(evaluate(derivative(px,1),x,y))*k
        +absoluteUpper(evaluate(derivative(py,1),x,y))*k*k)/m;
    const Q width=x.hi-x.lo;
    // Lagrange interpolation remainder; matching abscissae gives both
    // directed Euclidean distances, not only curve-to-chord sampling error.
    out.chordDistanceUpperBound=out.curvatureUpperBound*width*width/Q(8);
    if(out.chordDistanceUpperBound>maximumDistance)
        return std::string("regular graph chord bound exceeds requested accuracy");
    for(const auto& other:obstacles){
        int otherSign=strictSign(evaluate(other,x,y));
        if(!otherSign)otherSign=strictSign(tensorBernsteinRange(other,x,y));
        if(!otherSign)return std::string("another curve cannot be excluded from whole strip");
        out.obstacleSigns.push_back(otherSign);
    }
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-regular-strip-v1\n";
    for(const auto& c:p)record<<c<<' ';record<<'\n';
    for(const auto& c:box)record<<c<<' ';record<<'\n'<<maximumDistance<<'\n';
    record<<obstacles.size()<<'\n';
    for(const auto& other:obstacles){for(const auto& c:other)record<<c<<' ';record<<'\n';}
    verification::Sha256 hash;const auto bytes=record.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    out.inputIdentity=hash.finish();return out;
}
RankTwoRegularStripCertificate certifyRankTwoRegularStrip(const P& p,const RankTwoStripBox& box,
    const std::vector<P>& obstacles,const Q& maximumDistance){
    return requireCertificate(tryCertifyRankTwoRegularStrip(p,box,obstacles,maximumDistance));
}

struct RankTwoRegularStripMemo::Impl {
    using Context=std::tuple<P,std::vector<P>>;
    struct ContextLess {
        using is_transparent=void;
        template<class A,class B>bool operator()(const A& a,const B& b)const {
            if(std::get<0>(a)<std::get<0>(b))return true;
            if(std::get<0>(b)<std::get<0>(a))return false;
            return std::get<1>(a)<std::get<1>(b);
        }
    };
    using Result=RankTwoRegularStripAttempt;
    using Query=std::tuple<RankTwoStripBox,Q>;
    using Queries=std::map<Query,Result>;
    using Contexts=std::map<Context,Queries,ContextLess>;
    using NamedContext=std::tuple<std::string,int,std::vector<std::string>>;
    struct NamedContextLess {
        using is_transparent=void;
        template<class A,class B>bool operator()(const A& a,const B& b)const {
            if(std::get<0>(a)<std::get<0>(b))return true;
            if(std::get<0>(b)<std::get<0>(a))return false;
            if(std::get<1>(a)!=std::get<1>(b))return std::get<1>(a)<std::get<1>(b);
            return std::get<2>(a)<std::get<2>(b);
        }
    };
    // Intern exact contexts once; lookups compare references, never rebuild
    // or copy 25*factor-count rationals for each of the million box trials.
    Contexts contexts;
    std::map<std::pair<std::string,int>,P> parsed;
    // Text identities resolve to owned coefficient contexts, not cached
    // geometric/authority checks. Map insertion never invalidates iterators.
    std::map<NamedContext,Contexts::iterator,NamedContextLess> named;
    Contexts::iterator intern(const P& p,const std::vector<P>& obstacles){
        auto context=contexts.find(std::tie(p,obstacles));
        if(context==contexts.end())
            context=contexts.emplace(std::make_tuple(p,obstacles),Queries{}).first;
        return context;
    }
    Result attempt(Contexts::iterator context,const RankTwoStripBox& box,const Q& error){
        auto key=std::make_tuple(box,error);
        const auto found=context->second.find(key);
        if(found!=context->second.end()){
            rankTwoDiagnostic("regular-strip:memo-hit");
            return found->second;
        }
        rankTwoDiagnostic("regular-strip:memo-miss");
        // Validation exceptions escape before insertion. Only the actual
        // exact mathematical attempt may become an owned reusable result.
        auto result=tryCertifyRankTwoRegularStrip(std::get<0>(context->first),box,
            std::get<1>(context->first),error);
        context->second.emplace(std::move(key),result);return result;
    }
};
RankTwoRegularStripMemo::RankTwoRegularStripMemo():impl_(std::make_unique<Impl>()){}
RankTwoRegularStripMemo::~RankTwoRegularStripMemo()=default;
const P& RankTwoRegularStripMemo::parsedChart(const std::string& saved,int axis){
    const auto key=std::make_pair(saved,axis);
    const auto found=impl_->parsed.find(key);
    if(found!=impl_->parsed.end())return found->second;
    // Parse before insertion. Malformed inputs never become cached witnesses.
    const auto p=parseRankTwoStripPolynomial(saved,axis);
    return impl_->parsed.emplace(key,p).first->second;
}
RankTwoRegularStripCertificate RankTwoRegularStripMemo::certify(const P& p,
    const RankTwoStripBox& box,const std::vector<P>& obstacles,const Q& error){
    return requireCertificate(impl_->attempt(impl_->intern(p,obstacles),box,error));
}
RankTwoRegularStripCertificate RankTwoRegularStripMemo::certifyChart(const std::string& saved,int axis,
    const RankTwoStripBox& box,const std::vector<std::string>& obstacles,const Q& error){
    return requireCertificate(tryCertifyChart(saved,axis,box,obstacles,error));
}
RankTwoRegularStripAttempt RankTwoRegularStripMemo::tryCertifyChart(const std::string& saved,int axis,
    const RankTwoStripBox& box,const std::vector<std::string>& obstacles,const Q& error){
    if(axis!=0 && axis!=1)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","regular edge requires valid chart");
    auto named=impl_->named.find(std::tie(saved,axis,obstacles));
    if(named==impl_->named.end()){
        const auto& p=parsedChart(saved,axis);std::vector<P> decoded;
        decoded.reserve(obstacles.size());
        for(const auto& other:obstacles)decoded.push_back(parsedChart(other,axis));
        const auto context=impl_->intern(p,decoded);
        named=impl_->named.emplace(std::make_tuple(saved,axis,obstacles),context).first;
    }
    return impl_->attempt(named->second,box,error);
}
} // namespace spin::detail
