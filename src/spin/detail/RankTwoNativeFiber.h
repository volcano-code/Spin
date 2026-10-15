#pragma once
#include "RankTwoAlgebraic.h"
#include "RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <algorithm>
#include <locale>
#include <optional>
#include <sstream>

namespace spin::detail {
using RankTwoFiberTraits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using RankTwoNativeFiberPoint=RankTwoFiberTraits::Algebraic_kernel_d_2::Algebraic_real_2;

// Private in-memory provenance, never public certificate fields or a global
// cache. ALL original operands are frozen by value before any native handle
// comparison can rewrite its defining polynomial.
class RankTwoNativeFiberWitness {
    friend struct RankTwoNativeFiberAccess;
    std::vector<std::vector<std::string>> curve_;
    std::vector<std::string> abscissa_;
    AlgebraicPolynomial complete_;
    RankTwoNativeFiberWitness(std::vector<std::vector<std::string>> curve,
        std::vector<std::string> abscissa,const AlgebraicPolynomial& complete)
        :curve_(std::move(curve)),abscissa_(std::move(abscissa)),complete_(complete){}
};
inline AlgebraicReal rankTwoNativeFiberOrdinate(const RankTwoNativeFiberPoint&,int,
    std::optional<RankTwoNativeRootSet>*);
struct RankTwoNativeFiberAccess {
private:
    friend AlgebraicReal rankTwoNativeFiberOrdinate(const RankTwoNativeFiberPoint&,int,
        std::optional<RankTwoNativeRootSet>*);
    template<class Poly> static std::vector<std::string> coefficients(const Poly& p){
        std::vector<std::string> out;
        for(int i=0;i<=p.degree();++i){std::ostringstream s;s.imbue(std::locale::classic());s<<p[i];out.push_back(s.str());}
        return out;
    }
    static std::vector<std::vector<std::string>> curveCoefficients(const RankTwoFiberTraits::Polynomial_2& p){
        std::vector<std::vector<std::string>> out;
        for(int i=0;i<=p.degree();++i)out.push_back(coefficients(p[i]));return out;
    }
    // ONLY the actual getter may bind the already computed ORIGINAL resultant.
    static RankTwoNativeRootSet bind(const RankTwoNativeRootSet& proof,
        const RankTwoFiberTraits::Polynomial_2& f,const AlgebraicPolynomial& u,
        const AlgebraicPolynomial& complete){
        auto witness=std::shared_ptr<const RankTwoNativeFiberWitness>(new RankTwoNativeFiberWitness(
            curveCoefficients(f),coefficients(u),complete));
        return proof.withNativeFiberWitness(std::move(witness));
    }
public:
    static std::shared_ptr<const RankTwoNativeFiberWitness> witness(const CertifiedAlgebraicReal& expected){
        return expected.nativeFiberWitness();
    }
    static bool matches(const RankTwoNativeFiberWitness& witness,
        const RankTwoFiberTraits::Polynomial_2& f,const AlgebraicPolynomial& u){
        return witness.curve_==curveCoefficients(f)&&witness.abscissa_==coefficients(u);
    }
    static const AlgebraicPolynomial& complete(const RankTwoNativeFiberWitness& witness){return witness.complete_;}
};

// Arithmetic-equivalent ordinary-fiber adapter. The actual native point's
// guaranteed ordinate bounds select ONE root from the authenticated COMPLETE
// original resultant root set. Bounds are never substituted as coordinates.
// Critical and low-degree cases retain the original native getter verbatim.
inline AlgebraicReal rankTwoNativeFiberOrdinate(const RankTwoNativeFiberPoint& point,int maximumBits=512,
    std::optional<RankTwoNativeRootSet>* scopedProof=nullptr){
    requireExactAlgebraicEnvironment();
    if(scopedProof)scopedProof->reset();
    if(maximumBits<1 || maximumBits>512)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","fiber isolation budget outside supported range");
#ifndef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    return point.y();
#else
    const auto curve=point.curve();
    const auto line=curve.status_line_at_exact_x(point.x());
    // Only the measured high-degree ordinary path is optimized. This is not
    // unsupported geometry: all other cases use the SAME original provider.
    if(line.is_event() || point.x().degree()<=8)return point.y();
    using P=RankTwoFiberTraits::Polynomial_2;
    using PT=CGAL::Polynomial_traits_d<P>;
    rankTwoDiagnostic("native-fiber:resultant");
    const auto abscissa=point.x().polynomial();
    const P f=curve.primitive_polynomial_2(),r(abscissa);
    AlgebraicPolynomial original;
    const auto complete=RankTwoNativeRootSet::fromResultant(PT::Swap()(f,0,1),PT::Swap()(r,0,1),original);
    const auto proof=RankTwoNativeFiberAccess::bind(complete,
        f,abscissa,original);
    auto roots=proof.roots();
    const auto count="native-fiber:factor-root-set-"+std::to_string(roots.size());rankTwoDiagnostic(count.c_str());
    if(roots.empty())throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","ordinary native fiber has no real roots");
    rankTwoDiagnostic("native-fiber:select-branch");
    const auto approximate=curve.kernel()->approximate_absolute_y_2_object();
    int bits=std::min(16,maximumBits);auto bounds=approximate(point,bits);
    std::vector<std::size_t> candidates;
    for(std::size_t i=0;i<roots.size();++i)
        if(!(roots[i].high()<bounds.first || bounds.second<roots[i].low()))candidates.push_back(i);
    while(candidates.size()>1){
        if(bits==maximumBits)throw RankTwoSourceError("RESOURCE_LIMIT","ordinary fiber branch isolation budget exhausted");
        bits=std::min(maximumBits,bits*2);bounds=approximate(point,bits);
        std::vector<std::size_t> remaining;
        for(auto i:candidates)if(!(roots[i].high()<bounds.first || bounds.second<roots[i].low()))remaining.push_back(i);
        candidates=std::move(remaining);
    }
    if(candidates.empty())throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","native branch bounds miss complete resultant root set");
    const std::size_t selected=candidates.front();auto chosen=roots[selected];
    // Preserve rational-root definition exactly as the original getter does.
    if(chosen.is_rational())return point.y();
    AlgebraicKernel kernel;
    for(int step=0;step<=4*maximumBits;++step){
        const auto lo=chosen.low(),hi=chosen.high();bool unique=lo<hi;
        if(original.sign_at(lo)==CGAL::ZERO || original.sign_at(hi)==CGAL::ZERO)unique=false;
        for(std::size_t i=0;i<roots.size() && unique;++i)if(i!=selected
            && roots[i].compare(lo)!=CGAL::SMALLER && roots[i].compare(hi)!=CGAL::LARGER)unique=false;
        if(unique){
            // Returning a factor root would change the frozen certificate and
            // original ordinal. Reconstruct the ORIGINAL defining polynomial
            // in this independently isolated selected-root interval instead.
            rankTwoDiagnostic("native-fiber:original-root");
            auto result=kernel.construct_algebraic_real_1_object()(original,lo,hi);
            if(scopedProof)*scopedProof=proof;
            return result;
        }
        if(step==4*maximumBits)break;
        chosen.refine();
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","cannot isolate original native fiber root within budget");
#endif
}
// Exact ordinary-fiber membership without enumerating a root set already
// certified by expected. Only a validated arithmetic object is accepted, not
// caller-authored certificate fields. Its ORIGINAL polynomial and unique
// strict isolator remain immutable; guaranteed bounds are never coordinates.
inline bool rankTwoNativeFiberMatches(const RankTwoNativeFiberPoint& point,
    const CertifiedAlgebraicReal& expected){
    requireExactAlgebraicEnvironment();
    const auto original=[&]{rankTwoDiagnostic("native-fiber:membership-original-fallback");
        return rankTwoNativeFiberOrdinate(point).compare(rankTwoAlgebraicHandle(expected))==CGAL::EQUAL;};
#ifndef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
    return original();
#else
    const auto curve=point.curve();
    const auto line=curve.status_line_at_exact_x(point.x());
    if(line.is_event() || point.x().degree()<=8)return original();
    const auto& certificate=expected.certificate();
    const AlgebraicKernel::Bound lower(CORE::BigRat(certificate.lower.c_str()));
    const AlgebraicKernel::Bound upper(CORE::BigRat(certificate.upper.c_str()));
    const auto approximate=curve.kernel()->approximate_absolute_y_2_object();
    int bits=16;auto bounds=approximate(point,bits);
    const auto capture=[&](const char* label,const AlgebraicPolynomial* polynomial=nullptr){
        if(!rankTwoDiagnosticTarget.sink)return;
        std::ostringstream out;out.imbue(std::locale::classic());
        out<<"native-fiber:membership-bound-proof "<<label<<' '<<certificate.rootIdentity
            <<' '<<certificate.rootIndex<<' '<<certificate.lower<<' '<<certificate.upper
            <<' '<<bits<<' '<<bounds.first<<' '<<bounds.second;
        if(polynomial){out<<' '<<polynomial->degree()+1;
            for(int i=0;i<=polynomial->degree();++i)out<<' '<<(*polynomial)[i];}
        rankTwoDiagnostic(out.str().c_str());
    };
    while(!(lower<bounds.first && bounds.second<upper)){
        // The native branch value is in its GUARANTEED bounds; the validated
        // expected root is strictly in its own isolator. Disjointness proves
        // inequality independently of any polynomial-representation choice.
        if(bounds.second<=lower || upper<=bounds.first){
            capture("disjoint");rankTwoDiagnostic("native-fiber:membership-disjoint");return false;
        }
        if(bits==512){rankTwoDiagnostic("native-fiber:membership-isolator-unresolved");return original();}
        bits=std::min(512,bits*2);bounds=approximate(point,bits);
    }
    using P=RankTwoFiberTraits::Polynomial_2;
    using PT=CGAL::Polynomial_traits_d<P>;
    const auto abscissa=point.x().polynomial();
    const P f=curve.primitive_polynomial_2(),r(abscissa);
    const auto witness=RankTwoNativeFiberAccess::witness(expected);
    const auto complete=[&]{
        if(witness&&RankTwoNativeFiberAccess::matches(*witness,f,abscissa)){
            rankTwoDiagnostic("native-fiber:elimination-proof-reused");
            return RankTwoNativeFiberAccess::complete(*witness);
        }
        if(witness)rankTwoDiagnostic("native-fiber:elimination-proof-mismatch");
        rankTwoDiagnostic("native-fiber:membership-resultant");
        return rankTwoSquareFreeNativePolynomial(rankTwoCompleteNativeResultant(PT::Swap()(f,0,1),PT::Swap()(r,0,1)));
    }();
    bool same=complete.degree()+1==static_cast<int>(certificate.coefficients.size());
    for(int i=0;i<=complete.degree()&&same;++i){
        std::ostringstream out;out.imbue(std::locale::classic());out<<complete[i];
        same=out.str()==certificate.coefficients[i];
    }
    if(!same){rankTwoDiagnostic("native-fiber:membership-polynomial-mismatch");return original();}
    // Actual native value is a root of this ORIGINAL eliminant. Both roots
    // now lie in the validated isolator containing EXACTLY ONE distinct real
    // root of that same polynomial. No interval overlap/zero-set shortcut.
    capture("contained",&complete);rankTwoDiagnostic("native-fiber:membership-certified");return true;
#endif
}
} // namespace spin::detail
