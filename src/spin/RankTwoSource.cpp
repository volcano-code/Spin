#include "spin/RankTwoSource.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <locale>
#include <sstream>
#include <utility>

namespace {
using Q=CGAL::Gmpq;
using Z=CGAL::Gmpz;
using Matrix=std::array<std::array<Q,3>,3>;
Q absolute(const Q& q) { return q<0?-q:q; }
int sign(const Q& q) { return q>0?1:(q<0?-1:0); }
std::string sha(const std::string& content) {
    verification::Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(content.data()),content.size());
    return hash.finish();
}
Q cofactor(const Matrix& a,int row,int column) {
    std::array<int,2> rows{},columns{};
    for(int i=0,r=0,c=0;i<3;++i){
        if(i!=row)rows[r++]=i;
        if(i!=column)columns[c++]=i;
    }
    const Q minor=a[rows[0]][columns[0]]*a[rows[1]][columns[1]]
                 -a[rows[0]][columns[1]]*a[rows[1]][columns[0]];
    return (row+column)%2?-minor:minor;
}
spin::RankTwoExactCoefficients primitive(const spin::RankTwoExactCoefficients& input) {
    Z denominator(1),divisor(0);
    for(const auto& q:input)
        mpz_lcm(denominator.mpz(),denominator.mpz(),q.denominator().mpz());
    std::array<Z,10> integers{};
    for(int i=0;i<10;++i){
        integers[i]=input[i].numerator()*(denominator/input[i].denominator());
        mpz_gcd(divisor.mpz(),divisor.mpz(),integers[i].mpz());
    }
    // This function is only called after nonzero rank-two admission.
    spin::RankTwoExactCoefficients result{};
    for(int i=0;i<10;++i)result[i]=Q(integers[i]/divisor);
    return result;
}
std::string identity(const spin::RankTwoExactCoefficients& input,
                     const std::string& frame,const char* version) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out<<version<<';'<<frame.size()<<':'<<frame<<';';
    for(const auto& q:input)out<<q.numerator()<<'/'<<q.denominator()<<';';
    return sha(out.str());
}
spin::RankTwoExactGraph makeGraph(const spin::RankTwoExactCoefficients& q) {
    const Matrix a{{{{q[4],q[9],q[8]}},{{q[9],q[5],q[7]}},{{q[8],q[7],q[6]}}}};
    Q determinant=0;
    for(int j=0;j<3;++j)determinant+=a[0][j]*cofactor(a,0,j);
    if(determinant!=0)
        throw spin::RankTwoSourceError("EXACT_STRATUM_MISMATCH","quadratic determinant is nonzero");
    std::array<Q,3> minors{};
    for(int k=0;k<3;++k)minors[k]=cofactor(a,k,k);
    int chosen=-1;Q largest=0;
    // Principal block order xy,xz,yz breaks ties deterministically.
    for(int k:{2,1,0})if(absolute(minors[k])>largest){
        chosen=k;largest=absolute(minors[k]);
    }
    if(chosen<0)
        throw spin::RankTwoSourceError("EXACT_STRATUM_MISMATCH","quadratic rank is less than two");
    spin::RankTwoExactGraph g;
    g.kernelCoordinate=chosen;
    for(int i=0;i<3;++i){
        g.kernel[i]=cofactor(a,chosen,i)/minors[chosen];
        g.beta+=g.kernel[i]*q[i+1];
    }
    if(g.beta==0)
        throw spin::RankTwoSourceError("EXACT_STRATUM_MISMATCH","kernel linear coefficient is zero");
    const bool reversed=g.beta<0;
    if(reversed){g.beta=-g.beta;for(auto& v:g.kernel)v=-v;}
    for(int axis=0,column=0;axis<3;++axis)if(axis!=chosen){
        g.plane[axis][column]=1;
        g.plane[chosen][column]=-g.kernel[axis]/g.kernel[chosen];
        ++column;
    }
    g.height[0]=-q[0]/g.beta;
    for(int i=0;i<3;++i){
        g.gradient[i][0]=q[i+1];
        for(int s=0;s<2;++s){
            g.height[1+s]-=q[i+1]*g.plane[i][s]/g.beta;
            for(int j=0;j<3;++j)g.gradient[i][1+s]+=2*a[i][j]*g.plane[j][s];
        }
        for(int j=0;j<3;++j){
            g.height[3]-=g.plane[i][0]*a[i][j]*g.plane[j][0]/g.beta;
            g.height[4]-=2*g.plane[i][0]*a[i][j]*g.plane[j][1]/g.beta;
            g.height[5]-=g.plane[i][1]*a[i][j]*g.plane[j][1]/g.beta;
        }
    }
    std::ostringstream branch;
    branch.imbue(std::locale::classic());
    branch<<"rank-two-orthogonal-chart-v1;"<<chosen<<';'<<reversed<<';';
    for(int i=0;i<3;++i){
        branch<<sign(minors[i])<<';';
        for(int j=i+1;j<3;++j)branch<<sign(absolute(minors[i])-absolute(minors[j]))<<';';
    }
    g.branchIdentity=sha(branch.str());
    return g;
}

struct NumericView {
    spin::QuadricCoefficients10 coefficients;
    std::array<Q,10> errors;
};
NumericView makeView(const spin::RankTwoExactCoefficients& input) {
    Q maximum=0;
    for(const auto& v:input)maximum=std::max(maximum,absolute(v));
    spin::RankTwoExactCoefficients scaled{};
    spin::QuadricCoefficients10::Vector view;
    Q normSquared=0;
    for(int i=0;i<10;++i){
        scaled[i]=input[i]/maximum;
        normSquared+=scaled[i]*scaled[i];
        view[i]=CGAL::to_double(scaled[i]);
    }
    auto normalized=spin::QuadricCoefficients10::fromUnnormalized(view);
    // sqrt(normSquared) is in [1,4]. Exact dyadic refinement encloses the
    // intended normalization, including error in the stored double norm.
    Q lower=1,upper=4;
    for(int bit=0;bit<128;++bit){
        const Q mid=(lower+upper)/2;
        if(mid*mid<=normSquared)lower=mid;else upper=mid;
    }
    std::array<Q,10> errors{};
    for(int i=0;i<10;++i){
        Q lo=scaled[i]/upper,hi=scaled[i]/lower;
        if(lo>hi)std::swap(lo,hi);
        const Q represented(normalized.values()[i]);
        errors[i]=std::max(absolute(represented-lo),absolute(represented-hi));
    }
    return {std::move(normalized),std::move(errors)};
}
std::string viewIdentity(const NumericView& view,const std::string& inputHash) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out<<"rank-two-positive-unit-view-maxscaled-dyadic128-v1;"<<inputHash<<';';
    for(int i=0;i<10;++i){
        std::uint64_t bits=0;
        const double value=view.coefficients.values()[i];
        static_assert(sizeof(bits)==sizeof(value),"binary64 numeric view required");
        std::memcpy(&bits,&value,sizeof(bits));
        out<<bits<<':'<<view.errors[i].numerator()<<'/'<<view.errors[i].denominator()<<';';
    }
    return sha(out.str());
}
} // namespace

namespace spin {

RankTwoSourceError::RankTwoSourceError(std::string code,std::string witness)
    : std::invalid_argument(code+": "+witness),code_(std::move(code)),witness_(std::move(witness)) {}

struct PreparedRankTwoSource::Impl {
    RankTwoExactCoefficients exact;
    RankTwoExactGraph graph;
    NumericView numeric;
    std::string inputHash,geometryHash,numericViewHash;
    Impl(RankTwoExactCoefficients q,const std::string& frame)
        : exact(std::move(q)),graph(makeGraph(exact)),numeric(makeView(exact)),
          inputHash(identity(exact,frame,"rank-two-source-input-v1")),
          geometryHash(identity(primitive(exact),frame,"rank-two-source-geometry-v1")),
          numericViewHash(viewIdentity(numeric,inputHash)) {}
};

PreparedRankTwoSource::PreparedRankTwoSource(std::shared_ptr<const Impl> impl)
    : impl_(std::move(impl)) {}

PreparedRankTwoSource PreparedRankTwoSource::fromExact(
    const RankTwoExactCoefficients& input,const std::string& frame) {
    if(frame.empty() || frame.find('\0')!=std::string::npos)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","nonempty NUL-free frame identity required");
    auto canonical=input;
    for(auto& q:canonical){
        if(q.denominator()==0)
            throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","zero rational denominator");
        q=Q(q.numerator(),q.denominator());
    }
    return PreparedRankTwoSource(std::make_shared<const Impl>(std::move(canonical),frame));
}

PreparedRankTwoSource PreparedRankTwoSource::fromRepresented(
    const QuadricCoefficients10& input,const std::string& frame) {
    RankTwoExactCoefficients exact{};
    for(int i=0;i<10;++i)exact[i]=Q(input.values()[i]);
    return fromExact(exact,frame);
}

const RankTwoExactCoefficients& PreparedRankTwoSource::exactCoefficients() const {return impl_->exact;}
const QuadricCoefficients10& PreparedRankTwoSource::normalizedView() const {return impl_->numeric.coefficients;}
const std::array<Q,10>& PreparedRankTwoSource::normalizedViewErrorBounds() const {return impl_->numeric.errors;}
const RankTwoExactGraph& PreparedRankTwoSource::graph() const {return impl_->graph;}
const std::string& PreparedRankTwoSource::inputHash() const {return impl_->inputHash;}
const std::string& PreparedRankTwoSource::geometryHash() const {return impl_->geometryHash;}
const std::string& PreparedRankTwoSource::numericViewHash() const {return impl_->numericViewHash;}

std::array<Q,25> RankTwoExactGraph::spherePolynomial(const Q& radius) const {
    if(radius<=0)throw std::invalid_argument("sphere radius must be positive");
    std::array<Q,25> result{};
    result[0]=-radius*radius;
    Q uSquared=0;
    for(int i=0;i<3;++i){
        uSquared+=kernel[i]*kernel[i];
        result[10]+=plane[i][0]*plane[i][0];
        result[6]+=2*plane[i][0]*plane[i][1];
        result[2]+=plane[i][1]*plane[i][1];
    }
    const int x[6]={0,1,0,2,1,0},y[6]={0,0,1,0,1,2};
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)
        result[5*(x[i]+x[j])+y[i]+y[j]]+=uSquared*height[i]*height[j];
    return result;
}
} // namespace spin
