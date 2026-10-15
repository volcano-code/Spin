#include "detail/RankTwoCarrierKernels.h"
#include "detail/RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include <CGAL/Polynomial_traits_d.h>
#include <algorithm>
#include <locale>
#include <sstream>
#include <vector>

namespace spin::detail {
namespace {
using Q=CGAL::Gmpq;
using Z=CGAL::Gmpz;
using I=CORE::BigInt;
using P=AlgebraicPolynomial;
using P2=CGAL::Polynomial<P>;
template<class T> std::string text(const T& x){
    std::ostringstream out;out.imbue(std::locale::classic());out<<x;return out.str();
}
CORE::BigRat bound(const Q& x){return CORE::BigRat(I(text(x.numerator())),I(text(x.denominator())));}
P polynomial(const CertifiedAlgebraicReal& value){
    // The private adapter has already validated polynomial, root ordinal,
    // strict isolator and hash. No untrusted certificate enters this path.
    std::vector<I> coefficients;
    for(const auto& x:value.certificate().coefficients)coefficients.emplace_back(x);
    return P(coefficients.begin(),coefficients.end());
}
void diagnosticRoot(std::ostream& out,const CertifiedAlgebraicReal& value){
    const auto& c=value.certificate();out<<c.rootIdentity<<' '<<c.rootIndex<<' '
        <<c.lower<<' '<<c.upper<<' ';
    for(const auto& coefficient:c.coefficients)out<<coefficient<<',';
}
[[noreturn]] void resource(const char* why){throw RankTwoSourceError("RESOURCE_LIMIT",why);}
}

CORE::Expr rankTwoExactLiftDifference(const CORE::Expr& upper,const CORE::Expr& lower){
    requireExactAlgebraicEnvironment();
    return upper.Rep()==lower.Rep()?CORE::Expr(0):upper-lower;
}

void certifyRankTwoTriangleThickness(const std::array<CORE::Expr,3>& lower,
    const std::array<CORE::Expr,3>& upper,const std::array<bool,3>& boundary){
    requireExactAlgebraicEnvironment();
    std::array<CORE::Expr,3> d;
    for(std::size_t i=0;i<3;++i){
        d[i]=rankTwoExactLiftDifference(upper[i],lower[i]);
        if(d[i].sign()<0)throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",
            "negative carrier thickness at corner "+std::to_string(i));
    }
    if((d[0]+d[1]+d[2]).sign()<=0)
        throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID","zero thickness throughout active carrier triangle");
    for(std::size_t i=0;i<3;++i)if(!boundary[i] && (d[i]+d[(i+1)%3]).sign()<=0)
        throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",
            "zero thickness along internal carrier edge "+std::to_string(i));
    requireExactAlgebraicEnvironment();
}

CertifiedAlgebraicReal sampleRankTwoAlgebraicCoordinate(
    const CertifiedAlgebraicReal& a,const CertifiedAlgebraicReal& b,const Q& t,int bits){
    requireExactAlgebraicEnvironment();
    if(t<=0 || t>=1 || bits<1 || bits>512)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","interior rational weight and valid isolation budget required");
    return affineRankTwoAlgebraicCoordinate(a,b,t,bits);
}

CertifiedAlgebraicReal affineRankTwoAlgebraicCoordinate(
    const CertifiedAlgebraicReal& a,const CertifiedAlgebraicReal& b,const Q& t,int bits){
    requireExactAlgebraicEnvironment();
    if(t<=-1 || t>=1 || t==0 || bits<1 || bits>512)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","bounded nonzero affine weight and valid isolation budget required");
    rankTwoDiagnostic("weighted-coordinate:capture-input");
    if(rankTwoDiagnosticTarget.sink){
        std::ostringstream capture;capture.imbue(std::locale::classic());
        capture<<"weighted-coordinate:input "<<t<<' '<<bits<<' ';
        diagnosticRoot(capture,a);capture<<' ';diagnosticRoot(capture,b);
        rankTwoDiagnostic(capture.str().c_str());
    }
    const auto finished=[](CertifiedAlgebraicReal result){
        rankTwoDiagnostic("weighted-coordinate:capture-output");
        if(rankTwoDiagnosticTarget.sink){
            std::ostringstream capture;capture.imbue(std::locale::classic());
            capture<<"weighted-coordinate:output ";diagnosticRoot(capture,result);
            rankTwoDiagnostic(capture.str().c_str());
        }
        rankTwoDiagnostic("weighted-coordinate:complete");return result;
    };
    rankTwoDiagnostic("weighted-coordinate:dependency-check");
    auto left=rankTwoAlgebraicHandle(a),right=rankTwoAlgebraicHandle(b);
    if(left.compare(right)==CGAL::EQUAL)return finished(CertifiedAlgebraicReal::fromRoot(left,bits));
    if(left.compare(-right)==CGAL::EQUAL){
        // Preserve a PROVED dependency instead of treating conjugate endpoint
        // values as independent eliminant variables. This is the SAME moving
        // coordinate, not a lower-resolution or approximate sampling rule.
        rankTwoDiagnostic("weighted-coordinate:certified-antipodal");
        const Q scale=1-Q(2)*t;
        AlgebraicKernel kernel;
        if(scale==0)return finished(CertifiedAlgebraicReal::fromRoot(
            kernel.construct_algebraic_real_1_object()(CORE::BigRat(0)),bits));
        const P p=polynomial(a);
        const I N(text(scale.numerator())),D(text(scale.denominator()));
        std::vector<I> transformed;
        for(int i=0;i<=p.degree();++i)
            transformed.push_back(p[i]*CGAL::ipower(D,i)*CGAL::ipower(N,p.degree()-i));
        P q(transformed.begin(),transformed.end());
        q=CGAL::Polynomial_traits_d<P>::Make_square_free()(q);
        // Use the certificate's STRICT bounds; a native rational handle may
        // have low()==high(). Signed scaling reverses the interval as needed.
        Q lo=scale*Q(a.certificate().lower),hi=scale*Q(a.certificate().upper);
        if(lo>hi)std::swap(lo,hi);
        return finished(CertifiedAlgebraicReal::fromRoot(
            kernel.construct_algebraic_real_1_object()(q,bound(lo),bound(hi)),bits));
    }
    rankTwoDiagnostic("weighted-coordinate:resultant-build");
    P p=polynomial(a),q=polynomial(b);
    if(p.degree()>256/q.degree()){
        // Preserve the legacy path and bytes for its supported operands.
        // For over-bound parents, authenticate ALL factors of each original
        // certificate, then their SAME selected values' exact membership.
        // This defines only new elimination operands, NOT replacement parent
        // certificates or a fake complete proof of the original polynomials.
        p=rankTwoAuthenticatedParentComponent(a);
        q=rankTwoAuthenticatedParentComponent(b);
    }
    if(p.degree()>256/q.degree())resource("weighted-coordinate resultant degree exceeds bound");
    const I B(text(t.numerator())),C(text(t.denominator())),A=C-B;
    // Nested polynomial: z is inner variable 0, u is eliminated variable 1.
    // P2(p) would mean p(z), NOT p(u). Construct p(u) explicitly.
    const P2 z=CGAL::shift(P2(1),1,0),u=CGAL::shift(P2(1),1,1);
    P2 pu(0);
    for(int i=p.degree();i>=0;--i)pu=pu*u+P2(p[i]);
    const P2 v=P2(C)*z-P2(A)*u;
    P2 homogeneous(0);
    for(int i=0;i<=q.degree();++i)
        homogeneous+=P2(q[i]*CGAL::ipower(B,q.degree()-i))*CGAL::ipower(v,i);
    P resultant=CGAL::Polynomial_traits_d<P2>::Resultant()(pu,homogeneous);
    if(CGAL::is_zero(resultant) || resultant.degree()<1)
        throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","weighted-coordinate eliminant is invalid");
    rankTwoDiagnostic("weighted-coordinate:resultant-square-free");
    resultant=CGAL::Polynomial_traits_d<P>::Make_square_free()(resultant);
    AlgebraicKernel kernel;std::vector<AlgebraicReal> roots;
    rankTwoDiagnostic("weighted-coordinate:solve-all-roots");
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    // Keep the same complete ORIGINAL proof alive through selection. The
    // private factory exact-matches its complete canonical coefficients;
    // native conjugate checks may rewrite selected's polynomial, in which
    // case the existing independent fresh-proof fallback remains mandatory.
    const auto completeProof=RankTwoNativeRootSet::fromPolynomial(resultant);
    roots=completeProof.roots();
#else
    kernel.solve_1_object()(resultant,true,std::back_inserter(roots));
#endif
    rankTwoDiagnostic("weighted-coordinate:select-isolator");
    const Q oneMinus=1-t;
    for(int precision=std::min(4,bits);;precision=std::min(bits,2*precision)){
        Z den(1);mpz_mul_2exp(den.mpz(),den.mpz(),precision);
        const Q step(Z(1),den);
        for(auto* value:{&left,&right}){
            int count=0;
            while(Q(text(value->high()))-Q(text(value->low()))>step){
                if(++count>4*bits)resource("weighted-coordinate input refinement budget exhausted");
                value->refine();
            }
        }
        // Outward bounds on the true affine combination, not a guessed root.
        const Q lo=(oneMinus*Q(text(left.low()))+t*Q(text(t>0?right.low():right.high())))*Q(den);
        const Q hi=(oneMinus*Q(text(left.high()))+t*Q(text(t>0?right.high():right.low())))*Q(den);
        Z floor,ceil;
        mpz_fdiv_q(floor.mpz(),lo.numerator().mpz(),lo.denominator().mpz());
        mpz_cdiv_q(ceil.mpz(),hi.numerator().mpz(),hi.denominator().mpz());
        const auto L=bound(Q(floor-1,den)),U=bound(Q(ceil+1,den));
        std::size_t count=0,index=0;
        for(std::size_t i=0;i<roots.size();++i)
            if(roots[i].compare(L)==CGAL::LARGER && roots[i].compare(U)==CGAL::SMALLER){++count;index=i;}
        if(count==1){
            auto selected=roots[index];
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
            if(!selected.is_rational()){
                // Factor roots prove the COMPLETE set, not the defining
                // polynomial/ordinal of the public affine certificate.
                // Authenticate a strict unique ORIGINAL-polynomial isolator
                // before reconstructing that same root and applying the
                // existing conjugate checks. No selected-factor substitution.
                bool reconstructed=false;
                for(int refinement=0;refinement<=4*bits;++refinement){
                    const auto lower=selected.low(),upper=selected.high();
                    bool unique=lower<upper&&resultant.sign_at(lower)!=CGAL::ZERO
                        &&resultant.sign_at(upper)!=CGAL::ZERO;
                    for(std::size_t i=0;i<roots.size()&&unique;++i)if(i!=index
                        &&roots[i].compare(lower)!=CGAL::SMALLER&&roots[i].compare(upper)!=CGAL::LARGER)unique=false;
                    if(unique){
                        rankTwoDiagnostic("weighted-coordinate:original-root");
                        selected=kernel.construct_algebraic_real_1_object()(resultant,lower,upper);
                        reconstructed=true;break;
                    }
                    if(refinement==4*bits)break;
                    selected.refine();
                }
                if(!reconstructed)resource("original affine root isolation exceeds fixed budget");
            }
#endif
            rankTwoDiagnostic("weighted-coordinate:conjugate-check");
            if(p==q){
                // The all-conjugate eliminant also admits b=a and b=-a,
                // regardless of the actual selected endpoint roots. Query
                // exact membership BEFORE freezing: CGAL splits by gcd and
                // keeps the factor containing this SAME isolated real root,
                // including when membership is false. Never assume that an
                // affine combination cannot coincide with an ancestor root.
                (void)selected.is_root_of(p);
                const Q scale=1-Q(2)*t;
                if(scale==0){
                    (void)selected.is_root_of(P(I(0),I(1)));
                }else{
                    const I N(text(scale.numerator())),D(text(scale.denominator()));
                    std::vector<I> transformed;
                    for(int i=0;i<=p.degree();++i)
                        transformed.push_back(p[i]*CGAL::ipower(D,i)*CGAL::ipower(N,p.degree()-i));
                    const P scaled(transformed.begin(),transformed.end());
                    (void)selected.is_root_of(scaled);
                }
                rankTwoDiagnostic("weighted-coordinate:selected-conjugate-factors");
            }
            return finished(CertifiedAlgebraicReal::fromProvenAffineRoot(selected,a,b,t,bits,
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
                &completeProof
#else
                nullptr
#endif
                ));
        }
        if(count==0)throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","weighted-coordinate interval lost all eliminant roots");
        if(precision==bits)break;
    }
    resource("weighted coordinate not isolated within fixed bit budget");
}
} // namespace spin::detail
