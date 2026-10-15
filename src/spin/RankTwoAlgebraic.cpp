#include "detail/RankTwoAlgebraic.h"
#include "detail/RankTwoDiagnostics.h"
#include "detail/RankTwoSquareFreeProof.h"
#include "spin/RankTwoSource.h"
#include "verification/Sha256.h"
#include <CGAL/CORE/poly/Sturm.h>
#include <algorithm>
#include <locale>
#include <map>
#include <optional>
#include <sstream>
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
#include <flint.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_poly_factor.h>
#include <fmpz_mpoly.h>
#include <nmod_poly.h>
#endif

namespace spin::detail {
namespace {
using Q=CGAL::Gmpq;
using Z=CGAL::Gmpz;
template<class T> std::string text(const T& value) {
    std::ostringstream out;out.imbue(std::locale::classic());out<<value;return out.str();
}
std::string fraction(const Q& value) {
    return text(value.numerator())+"/"+text(value.denominator());
}
[[noreturn]] void invalid(const char* why) {
    throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);
}
bool canonicalInteger(const std::string& s) {
    std::size_t i=(!s.empty()&&s[0]=='-')?1:0;
    if(i==s.size() || s.size()>100000)return false;
    if(s[i]=='0')return s.size()==1;
    if(s[i]<'1'||s[i]>'9')return false;
    return std::all_of(s.begin()+i+1,s.end(),[](char c){return c>='0'&&c<='9';});
}
Q rational(const std::string& s) {
    const auto slash=s.find('/');
    if(slash==std::string::npos || !canonicalInteger(s.substr(0,slash))
        || !canonicalInteger(s.substr(slash+1)))invalid("malformed rational isolator");
    const Z numerator(s.substr(0,slash)),denominator(s.substr(slash+1));
    if(denominator<=0)invalid("nonpositive rational denominator");
    const Q value(numerator,denominator);
    if(fraction(value)!=s)invalid("noncanonical rational isolator");
    return value;
}
CORE::BigRat coreRational(const Q& q) {
    return CORE::BigRat(CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator())));
}
CORE::BigFloat dyadic(const Q& q,int maximumBits) {
    const auto denominator=q.denominator();
    if(mpz_popcount(denominator.mpz())!=1)invalid("isolator endpoint is not dyadic");
    const auto power=mpz_scan1(denominator.mpz(),0);
    if(power>static_cast<unsigned long>(maximumBits))
        throw RankTwoSourceError("RESOURCE_LIMIT","isolator dyadic bits exceed bound");
    const CORE::BigFloat result(CORE::BigInt(text(q.numerator())),-static_cast<long>(power));
    if(!result.isExact() || result.BigRatValue()!=coreRational(q))
        invalid("dyadic conversion changed value");
    return result;
}
std::vector<std::string> canonicalCoefficients(const AlgebraicPolynomial& p) {
    std::vector<Z> coefficients;Z divisor=0;
    for(int i=0;i<=p.degree();++i){
        coefficients.emplace_back(text(p[i]));
        mpz_gcd(divisor.mpz(),divisor.mpz(),coefficients.back().mpz());
    }
    if(p.degree()<1 || divisor==0)invalid("nonconstant root polynomial required");
    if(coefficients.back()<0)divisor=-divisor;
    std::vector<std::string> out;
    for(const auto& c:coefficients)out.push_back(text(c/divisor));
    return out;
}
std::string rootHash(const AlgebraicRootCertificate& r) {
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<r.schema<<';'<<r.rootIndex<<';';
    for(const auto& c:r.coefficients)out<<c.size()<<':'<<c<<';';
    verification::Sha256 hash;const auto bytes=out.str();
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return hash.finish();
}
std::vector<CORE::BigInt> readCoefficients(const AlgebraicRootCertificate& r) {
    if(r.schema!="rank-two-algebraic-root-v1" || r.coefficients.size()<2
       || r.coefficients.size()>257)invalid("root schema or polynomial degree invalid");
    std::vector<CORE::BigInt> out;
    for(const auto& c:r.coefficients){
        if(!canonicalInteger(c))invalid("malformed polynomial coefficient");
        out.emplace_back(c);
    }
    const AlgebraicPolynomial p(out.begin(),out.end());
    if(canonicalCoefficients(p)!=r.coefficients)invalid("noncanonical root polynomial");
    if(rootHash(r)!=r.rootIdentity)invalid("root identity mismatch");
    return out;
}
// Only the private native factory may provide this proof. Its sorted complete
// root set proves the frozen polynomial ordinal and strict unique isolator;
// public reload always passes nullptr and independently revalidates both.
struct VerifiedRoot {
    CORE::Expr expression;
    AlgebraicReal handle;
};
// Private proof data, never persisted or caller-authored. It may be carried
// synchronously by the opaque internally minted complete root-set attachment.
struct VerifiedFactors {
    AlgebraicPolynomial original;
    std::vector<AlgebraicPolynomial> polynomials;
};
AlgebraicPolynomial completeNativeResultantWithComponents(
    const AlgebraicBivariatePolynomial&,const AlgebraicBivariatePolynomial&,
    std::vector<AlgebraicPolynomial>*);
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
// Factorization only supplies candidates. CGAL/CORE independently prove the
// product, selected native root and strict isolator; the saved certificate is
// NEVER replaced by the factor's polynomial or its different root ordinal.
struct NativePolynomial {
    fmpz_poly_t value;
    NativePolynomial(){fmpz_poly_init(value);}
    ~NativePolynomial(){fmpz_poly_clear(value);}
    NativePolynomial(const NativePolynomial&)=delete;
};
struct NativeFactors {
    fmpz_poly_factor_t value;
    NativeFactors(){fmpz_poly_factor_init(value);}
    ~NativeFactors(){fmpz_poly_factor_clear(value);}
    NativeFactors(const NativeFactors&)=delete;
};
struct NativeModPolynomial {
    nmod_poly_t value;
    NativeModPolynomial(){nmod_poly_init(value,rankTwoSquareFreePrime);}
    ~NativeModPolynomial(){nmod_poly_clear(value);}
    NativeModPolynomial(const NativeModPolynomial&)=delete;
};
CORE::BigInt nativeInteger(const fmpz* value) {
    char* raw=fmpz_get_str(nullptr,10,value);
    if(!raw)invalid("native factor coefficient allocation failed");
    const std::unique_ptr<char,decltype(&flint_free)> owned(raw,&flint_free);
    return CORE::BigInt(raw);
}
struct NativeResultantContext {
    fmpz_mpoly_ctx_t value;
    NativeResultantContext(){fmpz_mpoly_ctx_init(value,2,ORD_LEX);}
    ~NativeResultantContext(){fmpz_mpoly_ctx_clear(value);}
    NativeResultantContext(const NativeResultantContext&)=delete;
};
struct NativeMultivariatePolynomial {
    fmpz_mpoly_t value;const fmpz_mpoly_ctx_struct* context;
    explicit NativeMultivariatePolynomial(const NativeResultantContext& c):context(c.value){
        fmpz_mpoly_init(value,context);}
    ~NativeMultivariatePolynomial(){fmpz_mpoly_clear(value,context);}
    NativeMultivariatePolynomial(const NativeMultivariatePolynomial&)=delete;
};
struct NativeInteger {
    fmpz_t value;
    NativeInteger(){fmpz_init(value);}
    ~NativeInteger(){fmpz_clear(value);}
    NativeInteger(const NativeInteger&)=delete;
};
void transferResultantOperand(NativeMultivariatePolynomial& out,const AlgebraicBivariatePolynomial& p){
    NativeInteger coefficient;
    for(int i=0;i<=p.degree();++i)for(int j=0;j<=p[i].degree();++j){
        const auto decimal=text(p[i][j]);
        if(fmpz_set_str(coefficient.value,decimal.c_str(),10))invalid("resultant operand integer transfer failed");
        const ulong exponents[2]={static_cast<ulong>(i),static_cast<ulong>(j)};
        fmpz_mpoly_set_coeff_fmpz_ui(out.value,coefficient.value,exponents,out.context);
    }
    if(!fmpz_mpoly_is_canonical(out.value,out.context))invalid("resultant operand is not canonical");
    // Check every original coefficient, including explicit zeros. X0 is the
    // OUTER nested power, not CGAL shift(...,0)'s inner coordinate.
    for(int i=0;i<=p.degree();++i)for(int j=0;j<=p[i].degree();++j){
        const ulong exponents[2]={static_cast<ulong>(i),static_cast<ulong>(j)};
        fmpz_mpoly_get_coeff_fmpz_ui(coefficient.value,out.value,exponents,out.context);
        if(nativeInteger(coefficient.value)!=p[i][j])invalid("complete resultant operand transfer mismatch");
    }
}
std::optional<AlgebraicPolynomial> directNativeResultant(
    const AlgebraicBivariatePolynomial& a,const AlgebraicBivariatePolynomial& b,long long bound){
    NativeResultantContext context;
    NativeMultivariatePolynomial aa(context),bb(context),result(context);
    transferResultantOperand(aa,a);transferResultantOperand(bb,b);
    if(fmpz_mpoly_resultant(result.value,aa.value,bb.value,0,context.value)!=1)return std::nullopt;
    if(!fmpz_mpoly_is_canonical(result.value,context.value))invalid("raw resultant is not canonical");
    for(slong term=0;term<fmpz_mpoly_length(result.value,context.value);++term){
        if(!fmpz_mpoly_term_exp_fits_ui(result.value,term,context.value))invalid("raw resultant exponent does not fit");
        ulong exponents[2];fmpz_mpoly_get_term_exp_ui(exponents,result.value,term,context.value);
        if(exponents[0]!=0||exponents[1]>static_cast<ulong>(bound))
            invalid("raw resultant did not eliminate the original outer variable");
    }
    NativePolynomial univariate;
    if(fmpz_mpoly_get_fmpz_poly(univariate.value,result.value,1,context.value)!=1)
        invalid("raw resultant is not univariate");
    std::vector<CORE::BigInt> coefficients;
    for(slong i=0;i<=std::max<slong>(0,fmpz_poly_degree(univariate.value));++i){
        NativeInteger coefficient;fmpz_poly_get_coeff_fmpz(coefficient.value,univariate.value,i);
        coefficients.push_back(nativeInteger(coefficient.value));
    }
    return AlgebraicPolynomial(coefficients.begin(),coefficients.end());
}
std::optional<AlgebraicPolynomial> completeProductResultant(
    const AlgebraicBivariatePolynomial& a,const AlgebraicBivariatePolynomial& b,long long bound,
    std::vector<AlgebraicPolynomial>* components=nullptr){
    const auto pureOuter=[](const AlgebraicBivariatePolynomial& p){
        for(int i=0;i<=p.degree();++i)if(p[i].degree()!=0)return false;return true;};
    const bool right=b.degree()>0&&pureOuter(b),left=!right&&a.degree()>0&&pureOuter(a);
    if(!right&&!left)return std::nullopt;
    const auto& curve=right?a:b;const auto& univariate=right?b:a;
    if(curve.degree()<1||curve.degree()>4)return std::nullopt;
    std::vector<CORE::BigInt> originalCoefficients;NativePolynomial input;NativeInteger coefficient;
    for(int i=0;i<=univariate.degree();++i){
        originalCoefficients.push_back(univariate[i][0]);const auto decimal=text(univariate[i][0]);
        if(fmpz_set_str(coefficient.value,decimal.c_str(),10))invalid("complete product input parse failed");
        fmpz_poly_set_coeff_fmpz(input.value,i,coefficient.value);
    }
    for(int i=0;i<=univariate.degree();++i){
        fmpz_poly_get_coeff_fmpz(coefficient.value,input.value,i);
        if(nativeInteger(coefficient.value)!=univariate[i][0])invalid("complete product input roundtrip failed");
    }
    const AlgebraicPolynomial original(originalCoefficients.begin(),originalCoefficients.end());
    NativeFactors factors;fmpz_poly_factor(factors.value,input.value);
    const auto content=nativeInteger(&factors.value->c);
    AlgebraicPolynomial reconstructed(content),raw(CGAL::ipower(content,curve.degree()));
    int innerDegree=0;for(int i=0;i<=curve.degree();++i)innerDegree=std::max(innerDegree,curve[i].degree());
    if(factors.value->num<1)invalid("complete product lost all factors");
    for(slong k=0;k<factors.value->num;++k){
        const auto degree=fmpz_poly_degree(factors.value->p+k),multiplicity=factors.value->exp[k];
        if(degree<1||degree>original.degree()||multiplicity<1||multiplicity>original.degree())
            invalid("complete product factor degree/multiplicity invalid");
        std::vector<CORE::BigInt> coefficients;
        for(slong i=0;i<=degree;++i){
            fmpz_poly_get_coeff_fmpz(coefficient.value,factors.value->p+k,i);
            coefficients.push_back(nativeInteger(coefficient.value));
        }
        const AlgebraicPolynomial factor(coefficients.begin(),coefficients.end());
        if(factor.degree()!=degree)invalid("complete product factor lost leading coefficient");
        reconstructed*=CGAL::ipower(factor,static_cast<int>(multiplicity));
        if(reconstructed.degree()>original.degree())invalid("complete product exceeds original degree");
        const auto operand=CGAL::Polynomial_traits_d<AlgebraicBivariatePolynomial>::Swap()(
            AlgebraicBivariatePolynomial(factor),0,1);
        const auto part=directNativeResultant(curve,operand,degree*innerDegree);
        if(!part){rankTwoDiagnostic("native-fiber:raw-resultant-full-product-direct-fallback");return std::nullopt;}
        raw*=CGAL::ipower(*part,static_cast<int>(multiplicity));
        if(components)components->push_back(*part);
    }
    // This is a COMPLETE raw resultant identity, not a selected factor/root
    // optimization. Signed content, ALL factors (including complex-only ones)
    // and each multiplicity must reconstruct the original operand exactly.
    if(reconstructed!=original)invalid("complete product differs from original operand");
    if(left&&(a.degree()%2)&&(b.degree()%2))raw=-raw;
    if(raw.degree()>bound)invalid("complete raw product exceeds original degree bound");
    rankTwoDiagnostic("native-fiber:raw-resultant-full-product-verified");return raw;
}
bool certifySquareFreeCandidate(const AlgebraicPolynomial& candidate){
    NativeModPolynomial q,derivative,g,s,t;
    std::vector<std::uint64_t> coefficients;
    for(int i=0;i<=candidate.degree();++i){
        fmpz_t c;fmpz_init(c);const auto decimal=text(candidate[i]);
        const int bad=fmpz_set_str(c,decimal.c_str(),10);
        if(bad){fmpz_clear(c);invalid("square-free modular coefficient parse failed");}
        const auto residue=fmpz_fdiv_ui(c,rankTwoSquareFreePrime);fmpz_clear(c);
        coefficients.push_back(residue);nmod_poly_set_coeff_ui(q.value,i,residue);
    }
    if(coefficients.back()==0||nmod_poly_degree(q.value)!=candidate.degree())return false;
    nmod_poly_derivative(derivative.value,q.value);
    nmod_poly_xgcd(g.value,s.value,t.value,q.value,derivative.value);
    if(nmod_poly_degree(g.value)!=0||nmod_poly_get_coeff_ui(g.value,0)!=1)return false;
    const auto copy=[](const nmod_poly_struct* polynomial){
        std::vector<std::uint64_t> out;
        for(slong i=0;i<nmod_poly_length(polynomial);++i)
            out.push_back(nmod_poly_get_coeff_ui(polynomial,i));
        return out;
    };
    return rankTwoVerifySquareFreeModulo(coefficients,copy(s.value),copy(t.value));
}
AlgebraicPolynomial authenticatedSquareFree(const AlgebraicPolynomial& original){
    if(std::string(FLINT_VERSION)!="3.3.1"||flint_get_num_threads()!=1)
        invalid("square-free backend version or thread contract mismatch");
    NativePolynomial input;
    for(int i=0;i<=original.degree();++i){
        fmpz_t coefficient;fmpz_init(coefficient);const auto decimal=text(original[i]);
        const int bad=fmpz_set_str(coefficient,decimal.c_str(),10);
        if(bad){fmpz_clear(coefficient);invalid("square-free coefficient parse failed");}
        fmpz_poly_set_coeff_fmpz(input.value,i,coefficient);fmpz_clear(coefficient);
    }
    NativeFactors factors;rankTwoDiagnostic("native-fiber:square-free-begin");
    fmpz_poly_factor_squarefree(factors.value,input.value);
    AlgebraicPolynomial product(nativeInteger(&factors.value->c)),candidate(1);
    for(slong i=0;i<factors.value->num;++i){
        const auto degree=fmpz_poly_degree(factors.value->p+i),exponent=factors.value->exp[i];
        if(degree<1||exponent<1||degree>original.degree()||exponent>original.degree())
            invalid("square-free decomposition degree/multiplicity invalid");
        std::vector<CORE::BigInt> coefficients;std::vector<std::string> encoded;
        for(slong j=0;j<=degree;++j){
            coefficients.push_back(nativeInteger(factors.value->p[i].coeffs+j));
            encoded.push_back(text(coefficients.back()));
        }
        const AlgebraicPolynomial factor(coefficients.begin(),coefficients.end());
        if(canonicalCoefficients(factor)!=encoded)invalid("square-free group not primitive positive");
        product*=CGAL::ipower(factor,static_cast<int>(exponent));candidate*=factor;
        if(product.degree()>original.degree())invalid("square-free product exceeds original degree");
    }
    if(product!=original||candidate.degree()<1)invalid("square-free complete product mismatch");
    rankTwoDiagnostic("native-fiber:square-free-product");
    // The product of primitive positive groups is primitive positive (Gauss).
    // Exact product plus square-free candidate proves ALL distinct factors
    // of original occur exactly once, including complex-only factors.
    rankTwoDiagnostic("native-fiber:square-free-modular");
    if(!certifySquareFreeCandidate(candidate)){
        // A bad prime or unproved candidate is not evidence of geometry
        // failure. Retain the SAME original exact arithmetic operation.
        rankTwoDiagnostic("native-fiber:square-free-original-fallback");
        return CGAL::make_square_free(original);
    }
    rankTwoDiagnostic("native-fiber:square-free-certified");
    return candidate;
}
VerifiedFactors authenticateFactors(const std::vector<CORE::BigInt>& original,bool emit=true) {
    // Do not modify global thread or CORE arithmetic settings.
    if(std::string(FLINT_VERSION)!="3.3.1" || flint_get_num_threads()!=1)
        invalid("native factor backend version or thread contract mismatch");
    NativePolynomial input;
    for(std::size_t i=0;i<original.size();++i){
        const auto c=text(original[i]);
        // CORE BigInt already owns the exact GMP value; no decimal rounding.
        fmpz_t coefficient;fmpz_init(coefficient);
        const int status=fmpz_set_str(coefficient,c.c_str(),10);
        if(status){fmpz_clear(coefficient);invalid("native factor coefficient parse failed");}
        fmpz_poly_set_coeff_fmpz(input.value,static_cast<slong>(i),coefficient);
        fmpz_clear(coefficient);
    }
    NativeFactors factors;
    if(emit)rankTwoDiagnostic("algebraic-root:native-factor");
    fmpz_poly_factor(factors.value,input.value);
    const AlgebraicPolynomial native(original.begin(),original.end());
    AlgebraicPolynomial product(nativeInteger(&factors.value->c));
    VerifiedFactors proof{native,{}};
    for(slong i=0;i<factors.value->num;++i){
        const auto degree=fmpz_poly_degree(factors.value->p+i);
        const auto multiplicity=factors.value->exp[i];
        if(degree<1 || degree>native.degree() || multiplicity<1
            || multiplicity>native.degree())invalid("invalid native factor degree/multiplicity");
        std::vector<CORE::BigInt> coefficients;
        for(slong j=0;j<=degree;++j)
            coefficients.push_back(nativeInteger(factors.value->p[i].coeffs+j));
        const AlgebraicPolynomial factor(coefficients.begin(),coefficients.end());
        for(slong j=0;j<multiplicity;++j)product=product*factor;
        if(product.degree()>native.degree())invalid("native factor product degree exceeds original");
        proof.polynomials.push_back(factor);
    }
    if(product!=native || proof.polynomials.empty())invalid("native factor product mismatch");
    if(emit)rankTwoDiagnostic("algebraic-root:native-factor-product-verified");
    return proof;
}
VerifiedFactors authenticateCompleteResultantComponents(const AlgebraicPolynomial& original,
    const std::vector<AlgebraicPolynomial>& components){
    // Components come only from the synchronous complete resultant producer,
    // AFTER exact signed input/product verification. Every raw component is
    // still independently factorized and reconstructs its own FULL input.
    rankTwoDiagnostic("algebraic-root:native-factor");
    std::map<std::vector<std::string>,AlgebraicPolynomial> all;
    for(const auto& part:components){
        if(part==AlgebraicPolynomial(0))invalid("zero component in nonzero resultant proof");
        if(part.degree()==0)continue; // nonzero constants have no roots
        std::vector<CORE::BigInt> coefficients;
        for(int i=0;i<=part.degree();++i)coefficients.push_back(part[i]);
        const auto local=authenticateFactors(coefficients,false);
        if(local.original!=part)invalid("resultant component proof input mismatch");
        for(const auto& factor:local.polynomials){
            const auto key=canonicalCoefficients(factor);
            std::vector<CORE::BigInt> normalized;for(const auto& c:key)normalized.emplace_back(c);
            all.emplace(key,AlgebraicPolynomial(normalized.begin(),normalized.end()));
        }
    }
    AlgebraicPolynomial product(1);VerifiedFactors proof{original,{}};
    for(const auto& entry:all){product*=entry.second;proof.polynomials.push_back(entry.second);
        if(product.degree()>original.degree())invalid("complete component product exceeds original SF degree");}
    // Deduplication is of exactly equal primitive polynomials, never roots,
    // real-only factors, interval guesses or a reduced original certificate.
    if(product!=original||proof.polynomials.empty())invalid("complete resultant components differ from original SF polynomial");
    rankTwoDiagnostic("algebraic-root:native-factor-product-verified");
    rankTwoDiagnostic("algebraic-root:resultant-component-product-verified");
    return proof;
}
std::optional<VerifiedFactors> authenticateCompleteComponentSquareFree(
    const std::vector<AlgebraicPolynomial>& components){
    // ONLY the synchronous complete raw resultant producer calls this path.
    // It already reconstructed its signed original operand with ALL factors
    // and multiplicities. Authenticate every complete component once; their
    // distinct primitive-positive factor product is a candidate for the SAME
    // full square-free eliminant, including complex-only factors.
    rankTwoDiagnostic("algebraic-root:native-factor");
    std::map<std::vector<std::string>,AlgebraicPolynomial> all;
    for(const auto& part:components){
        if(part==AlgebraicPolynomial(0))return std::nullopt;
        if(part.degree()==0)continue;
        std::vector<CORE::BigInt> coefficients;
        for(int i=0;i<=part.degree();++i)coefficients.push_back(part[i]);
        const auto local=authenticateFactors(coefficients,false);
        if(local.original!=part)invalid("resultant component proof input mismatch");
        for(const auto& factor:local.polynomials){
            const auto key=canonicalCoefficients(factor);
            std::vector<CORE::BigInt> normalized;
            for(const auto& c:key)normalized.emplace_back(c);
            all.emplace(key,AlgebraicPolynomial(normalized.begin(),normalized.end()));
        }
    }
    VerifiedFactors proof{AlgebraicPolynomial(1),{}};
    for(const auto& entry:all){
        proof.original*=entry.second;proof.polynomials.push_back(entry.second);
    }
    if(proof.polynomials.empty())return std::nullopt;
    // Exact complete component products cover every original factor. The
    // independent degree-preserving modular Bezout proof excludes repeated
    // factors even if the candidate factorization were not irreducible.
    // A bad prime/degree loss retains the original whole-polynomial path.
    if(!certifySquareFreeCandidate(proof.original))return std::nullopt;
    rankTwoDiagnostic("algebraic-root:native-factor-product-verified");
    rankTwoDiagnostic("algebraic-root:resultant-component-product-verified");
    rankTwoDiagnostic("native-fiber:component-square-free-certified");
    return proof;
}
std::vector<AlgebraicReal> completeDistinctFactorRoots(const VerifiedFactors& proof) {
    AlgebraicKernel kernel;std::vector<AlgebraicReal> roots;
    for(const auto& factor:proof.polynomials){
        const auto phase="algebraic-root:isolate-factor-degree-"+std::to_string(factor.degree());
        rankTwoDiagnostic(phase.c_str());
        // false requests CGAL's own multiplicity handling; do not assume
        // unverified irreducibility or square-freeness from the backend.
        kernel.solve_1_object()(factor,false,std::back_inserter(roots));
    }
    rankTwoDiagnostic("algebraic-root:sort-distinct-factors");
    std::sort(roots.begin(),roots.end(),[](const auto& a,const auto& b){
        return a.compare(b)==CGAL::SMALLER;
    });
    roots.erase(std::unique(roots.begin(),roots.end(),[](const auto& a,const auto& b){
        return a.compare(b)==CGAL::EQUAL;
    }),roots.end());
    if(roots.size()>static_cast<std::size_t>(proof.original.degree()))
        invalid("factor distinct root count exceeds original degree");
    // Exact product proof plus complete root isolation of EVERY factor proves
    // the complete original distinct-real set. Content/multiplicities do not
    // alter this set. No interval midpoints or approximate sort/deduplication.
    const auto phase="algebraic-root:factor-root-set-"+std::to_string(roots.size());
    rankTwoDiagnostic(phase.c_str());
    return roots;
}
CORE::Polynomial<CORE::BigInt> authenticatedExpressionPolynomial(
    const std::vector<CORE::BigInt>& original,const AlgebraicReal& root,
    const CORE::BigFloat& lower,const CORE::BigFloat& upper,
    const VerifiedFactors* nativeFactorProof) {
    std::optional<VerifiedFactors> owned;
    if(!nativeFactorProof){owned.emplace(authenticateFactors(original));nativeFactorProof=&*owned;}
    const AlgebraicPolynomial native(original.begin(),original.end());
    if(nativeFactorProof->original!=native)invalid("factor proof belongs to different original polynomial");
    AlgebraicKernel kernel;std::vector<CORE::BigInt> selected;int matches=0;
    for(const auto& factor:nativeFactorProof->polynomials){
        if(kernel.sign_at_1_object()(factor,root)==CGAL::ZERO){
            ++matches;selected.clear();
            for(int i=0;i<=factor.degree();++i)selected.push_back(factor[i]);
        }
    }
    if(matches!=1)invalid("native factor selected root mismatch");
    const CORE::Polynomial<CORE::BigInt> result(selected);
    if(CORE::sign(result.evalExactSign(lower))==0
        || CORE::sign(result.evalExactSign(upper))==0)
        invalid("selected factor endpoint is a root");
    CORE::Sturm<CORE::BigInt> sturm(result);
    if(sturm.numberOfRoots(lower,upper)!=1)
        invalid("selected factor lacks unique strict root");
    const auto phase="algebraic-root:expression-degree-"+std::to_string(selected.size()-1);
    rankTwoDiagnostic(phase.c_str());
    return result;
}
#endif
VerifiedRoot reconstruct(const AlgebraicRootCertificate& r,int maximumBits,
    const AlgebraicReal* nativeProof=nullptr,const CORE::Expr* affineProof=nullptr,
    const VerifiedFactors* factorProof=nullptr) {
    requireExactAlgebraicEnvironment();
    if(maximumBits<1 || maximumBits>512)invalid("isolation budget outside supported range");
    if(affineProof && !nativeProof)invalid("affine proof requires private native root proof");
    if(factorProof && !nativeProof)invalid("factor proof requires private native root proof");
    const auto coefficients=readCoefficients(r);
    const CORE::Polynomial<CORE::BigInt> p(coefficients);
    const Q lower=rational(r.lower),upper=rational(r.upper);
    if(lower>=upper)invalid("strict ordered isolator required");
    const auto L=dyadic(lower,maximumBits),U=dyadic(upper,maximumBits);
    // CORE's rootOf interval constructor uses v.front() before checking size.
    // Verify its complete preconditions *before* constructing any root DAG.
    if(CORE::sign(p.evalExactSign(L))==0 || CORE::sign(p.evalExactSign(U))==0)
        invalid("isolator endpoint is a root");
    if(nativeProof){
        if(nativeProof->compare(coreRational(lower))!=CGAL::LARGER
            || nativeProof->compare(coreRational(upper))!=CGAL::SMALLER)
            invalid("native proof left strict certificate isolator");
        rankTwoDiagnostic("algebraic-root:native-isolator-proof");
    }else{
        rankTwoDiagnostic("algebraic-root:reload-sturm");
        CORE::Sturm<CORE::BigInt> sturm(p);
        if(sturm.numberOfRoots(L,U)!=1)invalid("isolator must contain exactly one distinct root");
        const int below=sturm.numberOfRootsBelow(L);
        if(below<0 || static_cast<std::size_t>(below)!=r.rootIndex)invalid("root index mismatch");
    }
    if(rankTwoDiagnosticTarget.sink){
        std::ostringstream capture;capture.imbue(std::locale::classic());
        capture<<"algebraic-root:validated-certificate "<<r.rootIdentity<<' '<<r.rootIndex
            <<' '<<r.lower<<' '<<r.upper<<' ';
        for(const auto& c:r.coefficients)capture<<c<<',';
        rankTwoDiagnostic(capture.str().c_str());
    }
    // The native interval constructor requires a square-free polynomial.
    // Preserve the original certificate and test its substitution on the
    // square-free root handle; repeated roots remain accepted as before.
    const AlgebraicPolynomial native(coefficients.begin(),coefficients.end());
    AlgebraicKernel kernel;
    const auto nativeRoot=nativeProof?*nativeProof:
        kernel.construct_algebraic_real_1_object()(kernel.make_square_free_1_object()(native),
            coreRational(lower),coreRational(upper));
    rankTwoDiagnostic("algebraic-root:reload-native-substitution");
    if(kernel.sign_at_1_object()(native,nativeRoot)!=CGAL::ZERO)
        invalid("native root failed original polynomial substitution");
    rankTwoDiagnostic(affineProof?"algebraic-root:proved-affine-expression":
        "algebraic-root:reload-expression");
    const CORE::Expr result=affineProof?*affineProof:CORE::rootOf(
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
        authenticatedExpressionPolynomial(coefficients,nativeRoot,L,U,factorProof),
#else
        p,
#endif
        CORE::BFInterval(L,U));
    // Pinned CORE ConstPolyRep stores this defining polynomial and isolates
    // its unique root in [L,U]. Together with Sturm and native sign-at above,
    // that is an exact substitution proof, not a numerical zero shortcut.
    // Re-expanding p(rootOf(p)) as a composite DAG can make the exact zero
    // separation bound grow catastrophically (actual degree-40 regression).
    rankTwoDiagnostic("algebraic-root:reload-interval-check");
    if(!(CORE::Expr(coreRational(lower))<result)
        || !(result<CORE::Expr(coreRational(upper))))
        invalid("reconstructed expression failed exact polynomial/interval check");
    requireExactAlgebraicEnvironment();
    rankTwoDiagnostic("algebraic-root:reload-complete");
    return {result,nativeRoot};
}
} // namespace

struct RankTwoNativeRootSet::Impl {
    std::vector<std::string> coefficients;
    std::optional<VerifiedFactors> factors;
    std::vector<AlgebraicReal> roots;
    std::shared_ptr<const RankTwoNativeFiberWitness> nativeFiberWitness;
    explicit Impl(const AlgebraicPolynomial& p,std::optional<VerifiedFactors> prepared=std::nullopt)
        :coefficients(canonicalCoefficients(p)) {
        std::vector<CORE::BigInt> original;
        for(const auto& c:coefficients)original.emplace_back(c);
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
        if(prepared){
            const AlgebraicPolynomial canonical(original.begin(),original.end());
            if(prepared->original!=canonical)invalid("private component proof belongs to different original polynomial");
            factors=std::move(prepared);
        }else factors.emplace(authenticateFactors(original));
        roots=completeDistinctFactorRoots(*factors);
#else
        const AlgebraicPolynomial canonical(original.begin(),original.end());
        AlgebraicKernel{}.solve_1_object()(canonical,false,std::back_inserter(roots));
#endif
    }
};
RankTwoNativeRootSet::RankTwoNativeRootSet(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
RankTwoNativeRootSet RankTwoNativeRootSet::fromPolynomial(const AlgebraicPolynomial& p) {
    requireExactAlgebraicEnvironment();
    if(p.degree()<1||p.degree()>256)invalid("complete native root polynomial degree outside supported range");
    return RankTwoNativeRootSet(std::make_shared<const Impl>(p));
}
RankTwoNativeRootSet RankTwoNativeRootSet::fromResultant(const AlgebraicBivariatePolynomial& a,
    const AlgebraicBivariatePolynomial& b,AlgebraicPolynomial& original){
    requireExactAlgebraicEnvironment();std::vector<AlgebraicPolynomial> components;
    const auto raw=completeNativeResultantWithComponents(a,b,&components);
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    if(raw.degree()>0&&!components.empty()){
        auto proof=authenticateCompleteComponentSquareFree(components);
        if(proof){
            original=proof->original;
            if(original.degree()<1||original.degree()>256)
                invalid("complete native root polynomial degree outside supported range");
            rankTwoDiagnostic("native-fiber:complete-roots");
            return RankTwoNativeRootSet(std::make_shared<const Impl>(original,std::move(proof)));
        }
        rankTwoDiagnostic("native-fiber:component-square-free-original-fallback");
    }
#endif
    original=rankTwoSquareFreeNativePolynomial(raw);
    if(original.degree()<1||original.degree()>256)invalid("complete native root polynomial degree outside supported range");
    rankTwoDiagnostic("native-fiber:complete-roots");
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    if(!components.empty()){
        const auto encoded=canonicalCoefficients(original);std::vector<CORE::BigInt> coefficients;
        for(const auto& c:encoded)coefficients.emplace_back(c);
        const AlgebraicPolynomial canonical(coefficients.begin(),coefficients.end());
        auto proof=authenticateCompleteResultantComponents(canonical,components);
        return RankTwoNativeRootSet(std::make_shared<const Impl>(original,std::move(proof)));
    }
#endif
    // Unsupported or failed structural construction retains SAME complete
    // original factor/root path. No backend or geometry approximation switch.
    return fromPolynomial(original);
}
RankTwoNativeRootSet RankTwoNativeRootSet::withNativeFiberWitness(
    std::shared_ptr<const RankTwoNativeFiberWitness> witness) const {
    requireExactAlgebraicEnvironment();
    if(!impl_||!witness)invalid("missing private native elimination proof");
    auto copy=std::make_shared<Impl>(*impl_);copy->nativeFiberWitness=std::move(witness);
    return RankTwoNativeRootSet(std::move(copy));
}
const std::vector<AlgebraicReal>& RankTwoNativeRootSet::roots() const {
    requireExactAlgebraicEnvironment();
    if(!impl_)invalid("missing scoped complete root proof");
    return impl_->roots;
}

void requireExactAlgebraicEnvironment() {
    if(CORE::get_static_EscapePrec()!=CORE::extLong::getPosInfty()
       || CORE::get_static_EscapePrecFlag()<0)
        throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","CORE escape precision is finite or tainted");
}
struct CertifiedAlgebraicReal::Impl {
    AlgebraicRootCertificate certificate;
    VerifiedRoot value;
    std::shared_ptr<const RankTwoNativeFiberWitness> nativeFiberWitness;
    Impl(AlgebraicRootCertificate r,int bits,const AlgebraicReal* nativeProof=nullptr,
        const CORE::Expr* affineProof=nullptr,const VerifiedFactors* factorProof=nullptr,
        std::shared_ptr<const RankTwoNativeFiberWitness> witness=nullptr)
        :certificate(std::move(r)),value(reconstruct(certificate,bits,nativeProof,affineProof,factorProof)),
         nativeFiberWitness(std::move(witness)) {}
};
CertifiedAlgebraicReal::CertifiedAlgebraicReal(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
CertifiedAlgebraicReal CertifiedAlgebraicReal::reload(const AlgebraicRootCertificate& r,int bits) {
    return CertifiedAlgebraicReal(std::make_shared<const Impl>(r,bits));
}
CertifiedAlgebraicReal CertifiedAlgebraicReal::fromRoot(const AlgebraicReal& root,int bits) {
    return fromRootWithExpression(root,bits,nullptr);
}
CertifiedAlgebraicReal CertifiedAlgebraicReal::fromRoot(const AlgebraicReal& root,
    const RankTwoNativeRootSet& proof,int bits) {
    return fromRootWithExpression(root,bits,nullptr,&proof);
}
CertifiedAlgebraicReal CertifiedAlgebraicReal::fromProvenAffineRoot(const AlgebraicReal& root,
    const CertifiedAlgebraicReal& a,const CertifiedAlgebraicReal& b,const Q& t,int bits,
    const RankTwoNativeRootSet* sharedProof) {
    // The sampler has proved that this exact combination is the sole root in
    // its outward interval of the all-conjugate resultant. Factor selection
    // retains that root. Construct the expression from certified parents, not
    // from an approximate coordinate or a caller-supplied expression.
    const CORE::Expr affine=CORE::Expr(coreRational(1-t))*a.expression()
        +CORE::Expr(coreRational(t))*b.expression();
    return fromRootWithExpression(root,bits,&affine,sharedProof);
}
CertifiedAlgebraicReal CertifiedAlgebraicReal::fromRootWithExpression(
    const AlgebraicReal& root,int bits,const CORE::Expr* affineProof,
    const RankTwoNativeRootSet* sharedProof) {
    requireExactAlgebraicEnvironment();
    if(bits<1 || bits>512)invalid("isolation budget outside supported range");
    rankTwoDiagnostic("algebraic-root:certify");
    // Freeze before compare/refine: CGAL can simplify shared root handles.
    const AlgebraicPolynomial frozen=root.polynomial();
    if(rankTwoDiagnosticTarget.sink){
        const auto phase="algebraic-root:frozen-degree-"+std::to_string(frozen.degree());
        rankTwoDiagnostic(phase.c_str());
    }
    AlgebraicRootCertificate r;r.coefficients=canonicalCoefficients(frozen);
    std::vector<AlgebraicReal> roots;
    rankTwoDiagnostic("algebraic-root:solve-frozen");
    const VerifiedFactors* factorProof=nullptr;
    std::shared_ptr<const RankTwoNativeFiberWitness> nativeFiberWitness;
    std::optional<RankTwoNativeRootSet> ownedProof;
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    // Exact canonical polynomial equality, NOT a digest or selected-root value.
    // A conjugate membership check may simplify a native root definition; its
    // old proof is then inapplicable and must be independently reconstructed.
    if(sharedProof&&sharedProof->impl_&&sharedProof->impl_->coefficients==r.coefficients){
        rankTwoDiagnostic("algebraic-root:complete-proof-reused");
        nativeFiberWitness=sharedProof->impl_->nativeFiberWitness;
    }else{
        if(sharedProof)rankTwoDiagnostic("algebraic-root:complete-proof-mismatch");
        ownedProof.emplace(RankTwoNativeRootSet::fromPolynomial(frozen));sharedProof=&*ownedProof;
    }
    factorProof=&*sharedProof->impl_->factors;
    roots=sharedProof->impl_->roots;
#else
    AlgebraicKernel kernel;kernel.solve_1_object()(frozen,false,std::back_inserter(roots));
#endif
    rankTwoDiagnostic("algebraic-root:find-ordinal");
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    // A normal copy shares CGAL's rep/related-root list. Cross-factor compare
    // can learn a gcd and simplify every related caller root, changing later
    // feature certificates. Construct a detached exact comparison operand;
    // its rational value or original strict isolator identifies the SAME root.
    AlgebraicKernel ordinalKernel;
    const AlgebraicReal ordinalRoot=root.is_rational()
        ?ordinalKernel.construct_algebraic_real_1_object()(root.rational())
        :ordinalKernel.construct_algebraic_real_1_object()(frozen,root.low(),root.high());
#else
    const AlgebraicReal& ordinalRoot=root;
#endif
    bool found=false;
    for(std::size_t i=0;i<roots.size();++i)if(roots[i].compare(ordinalRoot)==CGAL::EQUAL){
        r.rootIndex=i;found=true;break;
    }
    if(!found)invalid("root not in its frozen defining polynomial");
    r.rootIdentity=rootHash(r);
    rankTwoDiagnostic("algebraic-root:enclose");
    AlgebraicReal isolated=root;
    // Enclose the exact handle with outward dyadics, then validate all roots.
    // No interval midpoint is used as an approximate root.
    for(int precision=std::min(4,bits);;precision=std::min(bits,2*precision)){
        Z denominator(1);mpz_mul_2exp(denominator.mpz(),denominator.mpz(),precision);
        const Q step(Z(1),denominator);
        int refinements=0;
        while(Q(text(isolated.high()))-Q(text(isolated.low()))>step){
            if(++refinements>4*bits)
                throw RankTwoSourceError("RESOURCE_LIMIT","root refinement budget exhausted");
            isolated.refine();
        }
        const Q low=Q(text(isolated.low()))*Q(denominator);
        const Q high=Q(text(isolated.high()))*Q(denominator);
        Z lo,hi;
        mpz_fdiv_q(lo.mpz(),low.numerator().mpz(),low.denominator().mpz());
        mpz_cdiv_q(hi.mpz(),high.numerator().mpz(),high.denominator().mpz());
        const Q L(lo-1,denominator),U(hi+1,denominator);
        bool unique=true;
        for(std::size_t i=0;i<roots.size();++i)if(i!=r.rootIndex
            && roots[i].compare(coreRational(L))!=CGAL::SMALLER
            && roots[i].compare(coreRational(U))!=CGAL::LARGER){unique=false;break;}
        if(unique){
            r.lower=fraction(L);r.upper=fraction(U);
            // Ordinal was proved against ALL roots of frozen above. Every
            // other root is outside the closed dyadic interval, and the
            // selected handle is strictly inside. Reuse this native proof,
            // not the public untrusted reload path. A sampler-proved affine
            // expression avoids reconstructing the same root as a new DAG;
            // ordinary roots and public reload still use CORE rootOf.
            return CertifiedAlgebraicReal(std::make_shared<const Impl>(r,bits,&root,affineProof,
                factorProof,std::move(nativeFiberWitness)));
        }
        if(precision==bits)break;
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","cannot certify distinct root within dyadic budget");
}
const CORE::Expr& CertifiedAlgebraicReal::expression() const {
    requireExactAlgebraicEnvironment();return impl_->value.expression;
}
const AlgebraicRootCertificate& CertifiedAlgebraicReal::certificate() const {return impl_->certificate;}
std::shared_ptr<const RankTwoNativeFiberWitness> CertifiedAlgebraicReal::nativeFiberWitness() const {
    requireExactAlgebraicEnvironment();return impl_->nativeFiberWitness;
}
AlgebraicReal rankTwoAlgebraicHandle(const CertifiedAlgebraicReal& value){
    requireExactAlgebraicEnvironment();
    // Both native construction and PUBLIC reload already validated this
    // same selected handle against the frozen polynomial and strict isolator.
    // Return a value copy,not an untrusted root or a process-global cache.
    return value.impl_->value.handle;
}
AlgebraicPolynomial rankTwoAuthenticatedParentComponent(const CertifiedAlgebraicReal& value){
    requireExactAlgebraicEnvironment();
    const auto& certificate=value.certificate();
    const auto coefficients=readCoefficients(certificate);
    const AlgebraicPolynomial original(coefficients.begin(),coefficients.end());
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    // The product contains EVERY original factor and multiplicity, including
    // complex-only ones. Only a native immutable value enters; no factor or
    // root-set proof is supplied by the caller.
    const auto proof=authenticateFactors(coefficients);
    if(proof.original!=original)invalid("parent component proof belongs to another polynomial");
    AlgebraicKernel kernel;
    // CGAL copies share related-root lists: membership can simplify all
    // related roots. Reconstruct from the SAME certified strict isolator to
    // keep this operation isolated from the parent's source authority.
    const auto detached=kernel.construct_algebraic_real_1_object()(
        kernel.make_square_free_1_object()(original),
        coreRational(rational(certificate.lower)),coreRational(rational(certificate.upper)));
    std::optional<AlgebraicPolynomial> selected;
    for(const auto& factor:proof.polynomials){
        if(kernel.sign_at_1_object()(factor,detached)==CGAL::ZERO){
            if(selected)invalid("parent root belongs to multiple authenticated components");
            selected=factor;
        }
    }
    if(!selected)invalid("parent root missing from complete authenticated factor product");
    if(rankTwoDiagnosticTarget.sink){
        std::ostringstream record;record.imbue(std::locale::classic());
        record<<"weighted-coordinate:authenticated-parent "<<certificate.rootIdentity<<' '
            <<certificate.rootIndex<<' '<<original.degree()<<' '<<selected->degree()<<' ';
        for(int i=0;i<=selected->degree();++i)record<<(*selected)[i]<<',';
        rankTwoDiagnostic(record.str().c_str());
    }
    requireExactAlgebraicEnvironment();return *selected;
#else
    return original;
#endif
}
std::vector<AlgebraicReal> rankTwoCompleteNativeRealRoots(const AlgebraicPolynomial& p){
    requireExactAlgebraicEnvironment();
    if(p.degree()<1 || p.degree()>256)invalid("complete native root polynomial degree outside supported range");
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    return RankTwoNativeRootSet::fromPolynomial(p).roots();
#else
    std::vector<AlgebraicReal> roots;
    AlgebraicKernel{}.solve_1_object()(p,false,std::back_inserter(roots));
    return roots;
#endif
}
AlgebraicPolynomial rankTwoSquareFreeNativePolynomial(const AlgebraicPolynomial& p){
    requireExactAlgebraicEnvironment();
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    // A quartic resultant of an existing <=256-degree coordinate has degree
    // <=1024. Other legacy inputs retain their original exact operation.
    if(p.degree()>=1&&p.degree()<=1024)return authenticatedSquareFree(p);
#endif
    return CGAL::make_square_free(p);
}
namespace {
AlgebraicPolynomial completeNativeResultantWithComponents(
    const AlgebraicBivariatePolynomial& a,const AlgebraicBivariatePolynomial& b,
    std::vector<AlgebraicPolynomial>* components){
    requireExactAlgebraicEnvironment();
    if(components)components->clear();
    const auto original=[&]{rankTwoDiagnostic("native-fiber:raw-resultant-original-fallback");
        return CGAL::resultant(a,b);};
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    const auto innerDegree=[](const AlgebraicBivariatePolynomial& p){
        int maximum=0;for(int i=0;i<=p.degree();++i)maximum=std::max(maximum,p[i].degree());return maximum;};
    const auto da=innerDegree(a),db=innerDegree(b);
    const long long bound=static_cast<long long>(a.degree())*db+static_cast<long long>(b.degree())*da;
    // The measured ordinary scope is quartic versus <=256-degree abscissa.
    // Do not infer geometry failure for any unsupported/zero/constant input.
    if(a==AlgebraicBivariatePolynomial(0)||b==AlgebraicBivariatePolynomial(0)
        ||(a.degree()==0&&b.degree()==0)||a.degree()>256||b.degree()>256
        ||da>256||db>256||std::min(a.degree(),b.degree())>4||bound>1024)return original();
    if(std::string(FLINT_VERSION)!="3.3.1"||flint_get_num_threads()!=1)
        invalid("resultant backend version or thread contract mismatch");
    rankTwoDiagnostic("native-fiber:raw-resultant-begin");
    auto out=completeProductResultant(a,b,bound,components);
    if(!out){if(components)components->clear();out=directNativeResultant(a,b,bound);}
    if(!out)return original();
    rankTwoDiagnostic("native-fiber:raw-resultant-complete");return *out;
#else
    return original();
#endif
}
} // namespace
AlgebraicPolynomial rankTwoCompleteNativeResultant(
    const AlgebraicBivariatePolynomial& a,const AlgebraicBivariatePolynomial& b){
    return completeNativeResultantWithComponents(a,b,nullptr);
}
} // namespace spin::detail
