#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoAlgebraic.h")
#include "spin/detail/RankTwoAlgebraic.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include "rank_two_degree40_certificate.h"
#include <CGAL/Cartesian.h>
#include <algorithm>
#include <vector>
#include <fstream>

namespace {
void require(bool ok,const char* why) {if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f) {
    try {f();} catch(const std::invalid_argument&) {return;}
    throw std::runtime_error("invalid root certificate accepted");
}
using namespace spin::detail;
void authenticatedFactorSelection() {
    // A reducible defining polynomial must not inflate the arithmetic DAG.
    // Its public certificate still identifies sqrt(2) as ordinal 2 of p,
    // not ordinal 1 of the selected irreducible factor x^2-2.
    int factorDegree=0;
    ScopedRankTwoDiagnostics observe([](const char* phase,void* context) noexcept {
        const std::string value(phase),prefix="algebraic-root:expression-degree-";
        if(value.rfind(prefix,0)==0)
            *static_cast<int*>(context)=std::atoi(value.c_str()+prefix.size());
    },&factorDegree);
    AlgebraicRootCertificate certificate;
    certificate.coefficients={"6","0","-5","0","1"};
    certificate.rootIndex=2;certificate.lower="11/8";certificate.upper="23/16";
    std::string bytes=certificate.schema+";2;";
    for(const auto& c:certificate.coefficients)
        bytes+=std::to_string(c.size())+":"+c+";";
    verification::Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    certificate.rootIdentity=hash.finish();
    const auto root=CertifiedAlgebraicReal::reload(certificate);
    require(factorDegree==2,"reducible root still uses full-degree arithmetic polynomial");
    require((root.expression()*root.expression()-2).sign()==0 && root.expression().sign()>0,
        "factor selection changed root/conjugate");
    require(root.certificate().coefficients==certificate.coefficients
        && root.certificate().rootIndex==2
        && root.certificate().rootIdentity==certificate.rootIdentity,
        "factor selection rewrote original root identity");
    // The same original polynomial's adjacent conjugate must select x^2-3.
    certificate.rootIndex=3;certificate.lower="27/16";certificate.upper="7/4";
    bytes=certificate.schema+";3;";
    for(const auto& c:certificate.coefficients)
        bytes+=std::to_string(c.size())+":"+c+";";
    verification::Sha256 nextHash;
    nextHash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    certificate.rootIdentity=nextHash.finish();
    const auto adjacent=CertifiedAlgebraicReal::reload(certificate);
    require((adjacent.expression()*adjacent.expression()-3).sign()==0
        && adjacent.expression()>root.expression(),"adjacent factor/conjugate confused");
    std::cout<<"PASS authenticated minimal factor, original ordinal/identity and adjacent conjugate\n";
}
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
void scopedCompleteRootProof() {
    // Hand-derived p=(x^2-2)(x^2-3); ONE minted complete proof may serve
    // adjacent original roots. A proof of x^2-2 is not a proof of p, even
    // though sqrt(2) has the same value. Public reload must still use Sturm.
    const AlgebraicPolynomial p(6,0,-5,0,1);
    const auto proof=RankTwoNativeRootSet::fromPolynomial(p);
    struct Events {int factor=0,reuse=0,mismatch=0,sturm=0;};
    Events events;
    ScopedRankTwoDiagnostics observe([](const char* phase,void* context) noexcept {
        auto& e=*static_cast<Events*>(context);const std::string name(phase);
        if(name=="algebraic-root:native-factor")++e.factor;
        if(name=="algebraic-root:complete-proof-reused")++e.reuse;
        if(name=="algebraic-root:complete-proof-mismatch")++e.mismatch;
        if(name=="algebraic-root:reload-sturm")++e.sturm;
    },&events);
    AlgebraicKernel kernel;
    const auto a=kernel.construct_algebraic_real_1_object()(p,CORE::BigRat(11,8),CORE::BigRat(23,16));
    const auto b=kernel.construct_algebraic_real_1_object()(p,CORE::BigRat(27,16),CORE::BigRat(7,4));
    const auto first=CertifiedAlgebraicReal::fromRoot(a,proof);
    const auto next=CertifiedAlgebraicReal::fromRoot(b,proof);
    require(events.factor==0&&events.reuse==2&&events.mismatch==0,
        "scoped complete proof redundantly recomputed or confused original roots");
    const std::vector<std::string> original={"6","0","-5","0","1"};
    require(first.certificate().coefficients==original&&first.certificate().rootIndex==2
        &&next.certificate().coefficients==original&&next.certificate().rootIndex==3,
        "scoped proof substituted a factor or adjacent ordinal");
    require((first.expression()*first.expression()-2).sign()==0&&first.expression()>0
        &&(next.expression()*next.expression()-3).sign()==0&&next.expression()>first.expression(),
        "scoped proof changed actual conjugates");
    require(a.polynomial()==p&&b.polynomial()==p,
        "scoped ordinal comparison mutated caller original polynomials");
    const auto smaller=RankTwoNativeRootSet::fromPolynomial(AlgebraicPolynomial(-2,0,1));
    const int before=events.factor;
    const auto sameValue=CertifiedAlgebraicReal::fromRoot(a,smaller);
    require(events.mismatch==1&&events.factor==before+1&&events.reuse==2
        &&sameValue.certificate().coefficients==original&&sameValue.certificate().rootIndex==2,
        "same-value factor proof was reused for a different original polynomial");
    const auto reload=CertifiedAlgebraicReal::reload(first.certificate());
    require(events.sturm==1&&(reload.expression()-first.expression()).sign()==0,
        "scoped proof bypassed public reload independent validation");
    auto wrong=first.certificate();wrong.rootIndex=3;
    std::string bytes=wrong.schema+";3;";
    for(const auto& c:wrong.coefficients)bytes+=std::to_string(c.size())+":"+c+";";
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    wrong.rootIdentity=hash.finish();rejects([&]{CertifiedAlgebraicReal::reload(wrong);});
    require(events.sturm==2,"rehash wrong ordinal bypassed independent public proof");
    std::cout<<"PASS scoped complete proof: adjacent original ordinals, mismatch fallback, external Sturm\n";
}
void ordinalComparisonIsolation() {
    // CGAL copies share a mutable rep/related-root list. Comparing a factor
    // root with the caller's native root can learn a gcd and rewrite shared
    // defining polynomials before enclosure, changing later feature IDs.
    AlgebraicKernel kernel;
    const AlgebraicPolynomial p(CORE::BigInt(6),CORE::BigInt(0),CORE::BigInt(-5),
        CORE::BigInt(0),CORE::BigInt(1));
    const auto native=kernel.construct_algebraic_real_1_object()(p,
        CORE::BigRat(11,8),CORE::BigRat(23,16));
    struct Observation {const AlgebraicReal* input;int degree=-1;};
    Observation state{&native,-1};
    const auto certified=[&]{
        ScopedRankTwoDiagnostics observe([](const char* phase,void* context) noexcept {
            auto& out=*static_cast<Observation*>(context);
            if(std::string(phase)=="algebraic-root:enclose")out.degree=out.input->polynomial().degree();
        },&state);
        return CertifiedAlgebraicReal::fromRoot(native);
    }();
    require(state.degree==4,"factor ordinal comparison mutated caller native definition before enclosure");
    require(certified.certificate().rootIndex==2 && certified.certificate().coefficients==
        std::vector<std::string>({"6","0","-5","0","1"}),"detached ordinal proof changed original root");
    require((certified.expression()*certified.expression()-2).sign()==0,
        "detached ordinal proof changed actual root");
    std::cout<<"PASS factor ordinal isolation: caller definition unchanged at enclosure\n";
}
void completeFactorRootSet() {
    // Hand-derived p=(x^2-2)(x^2-3)(x^2+1)(x^2-1). Six real roots:
    // -sqrt(3),-sqrt(2),-1,1,sqrt(2),sqrt(3). The complex-only factor
    // contributes none. Dropping a factor or using its LOCAL ordinal fails.
    const std::vector<std::string> original={"-6","0","5","0","5","0","-5","0","1"};
    std::vector<CORE::BigInt> coefficients;
    for(const auto& c:original)coefficients.emplace_back(c);
    const AlgebraicPolynomial p(coefficients.begin(),coefficients.end());
    const CORE::BigRat intervals[][2]={{CORE::BigRat(-7,4),CORE::BigRat(-27,16)},
        {CORE::BigRat(-23,16),CORE::BigRat(-11,8)},
        {CORE::BigRat(-17,16),CORE::BigRat(-15,16)},
        {CORE::BigRat(15,16),CORE::BigRat(17,16)},
        {CORE::BigRat(11,8),CORE::BigRat(23,16)},
        {CORE::BigRat(27,16),CORE::BigRat(7,4)}};
    const CORE::Expr expected[]={-CORE::sqrt(CORE::Expr(3)),-CORE::sqrt(CORE::Expr(2)),
        CORE::Expr(-1),CORE::Expr(1),CORE::sqrt(CORE::Expr(2)),CORE::sqrt(CORE::Expr(3))};
    struct Observed {int maximumDegree=0;int roots=-1;};
    AlgebraicKernel kernel;
    for(std::size_t i=0;i<6;++i){
        Observed events;
        const auto value=[&]{
            ScopedRankTwoDiagnostics observation([](const char* phase,void* context) noexcept {
                auto& out=*static_cast<Observed*>(context);const std::string event(phase);
                const std::string degree="algebraic-root:isolate-factor-degree-";
                const std::string roots="algebraic-root:factor-root-set-";
                if(event.rfind(degree,0)==0)
                    out.maximumDegree=std::max(out.maximumDegree,std::atoi(event.c_str()+degree.size()));
                if(event.rfind(roots,0)==0)out.roots=std::atoi(event.c_str()+roots.size());
            },&events);
            const auto native=kernel.construct_algebraic_real_1_object()(p,intervals[i][0],intervals[i][1]);
            return CertifiedAlgebraicReal::fromRoot(native);
        }();
        require(events.roots==6 && events.maximumDegree==2,
            "native factory still isolates original full-degree root set instead of authenticated factors");
        require(value.certificate().coefficients==original && value.certificate().rootIndex==i,
            "factor root union changed original polynomial or distinct-real ordinal");
        require((value.expression()-expected[i]).sign()==0,"factor root union changed actual selected root");
        const auto reload=CertifiedAlgebraicReal::reload(value.certificate());
        require((reload.expression()-expected[i]).sign()==0,"external original-polynomial proof disagrees");
    }
    std::cout<<"PASS authenticated complete root union: six original ordinals, complex factor, external proofs\n";
}
void actualDegree144() {
    // This is the original LAST degree-144 certificate captured by the full
    // failed coverage job, not its precomputed answer/selected factor.
    std::ifstream input(RANK_TWO_DEGREE144_FIXTURE);
    int count=0;input>>count;
    require(count==145,"actual degree-144 fixture count mismatch");
    std::vector<CORE::BigInt> coefficients;
    std::vector<std::string> original;
    for(int i=0;i<count;++i){
        std::string value;input>>value;
        require(!value.empty(),"truncated actual degree-144 fixture");
        original.push_back(value);coefficients.emplace_back(value);
    }
    std::string extra;require(!(input>>extra),"extra degree-144 fixture data");
    int selectedDegree=0;
    ScopedRankTwoDiagnostics observe([](const char* phase,void* context) noexcept {
        const std::string value(phase),prefix="algebraic-root:expression-degree-";
        if(value.rfind(prefix,0)==0)
            *static_cast<int*>(context)=std::atoi(value.c_str()+prefix.size());
    },&selectedDegree);
    AlgebraicKernel kernel;
    const AlgebraicPolynomial p(coefficients.begin(),coefficients.end());
    const auto native=kernel.construct_algebraic_real_1_object()(
        kernel.make_square_free_1_object()(p),CORE::BigRat(-66769,65536),CORE::BigRat(-33383,32768));
    const auto root=CertifiedAlgebraicReal::fromRoot(native);
    require(root.certificate().coefficients==original && root.certificate().rootIndex==27
        && root.certificate().rootIdentity==
            "e48e09980ed607940710c91c3febf094bc2875be878dcb66ec34b30466de2aa1",
        "factor selection changed original captured root certificate");
    require(selectedDegree==16,"actual high-degree root did not use authenticated smaller factor");
    require(CORE::Expr(CORE::BigRat(-66769,65536))<root.expression()
        && root.expression()<CORE::Expr(CORE::BigRat(-33383,32768)),
        "actual high-degree expression changed conjugate");
    std::cout<<"PASS actual degree-144 native root, degree-16 expression and original ordinal27 identity\n";
}
#endif
void nativeProofReuse() {
    // Catch redundant independent Sturm reconstruction on a root already
    // certified by the native factory. External reload must STILL execute it.
    int sturmPasses=0;
    ScopedRankTwoDiagnostics observe([](const char* phase,void* context) noexcept {
        if(std::string(phase)=="algebraic-root:reload-sturm")
            ++*static_cast<int*>(context);
    },&sturmPasses);
    AlgebraicKernel kernel;std::vector<AlgebraicReal> roots;
    kernel.solve_1_object()(AlgebraicPolynomial(-2,0,1),true,std::back_inserter(roots));
    const auto positive=CertifiedAlgebraicReal::fromRoot(roots.at(1));
    require(sturmPasses==0,"native proven root redundantly rebuilds certificate Sturm sequence");
    require((positive.expression()*positive.expression()-2).sign()==0
        && positive.expression().sign()>0,"native proof reuse changed selected root");
    const auto copy=CertifiedAlgebraicReal::reload(positive.certificate());
    require(sturmPasses==1,"untrusted reload bypassed its independent Sturm proof");
    require((copy.expression()-positive.expression()).sign()==0,"native/reloaded root mismatch");
    auto wrong=positive.certificate();wrong.rootIndex=0;
    // Rehash the false ordinal, so rejection must be mathematical, not a stale
    // digest. Root identity encoding is a public saved-certificate contract.
    std::string bytes=wrong.schema+";0;";
    for(const auto& c:wrong.coefficients)bytes+=std::to_string(c.size())+":"+c+";";
    verification::Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    wrong.rootIdentity=hash.finish();
    rejects([&]{CertifiedAlgebraicReal::reload(wrong);});
    require(sturmPasses==2,"rehash wrong ordinal was not independently checked");
    std::cout<<"PASS native proof reuse, full external reload and rehashed wrong ordinal rejection\n";
}
void degree40() {
    const auto certificate=rankTwoDegree40Certificate();
    const auto root=CertifiedAlgebraicReal::reload(certificate);
    // Hand-frozen strict isolator/ordinal from the actual source arc, not the
    // reconstructor's computed expected value. Must retain exact arithmetic.
    require(root.certificate().rootIdentity==
        "3cdc02ce6cdf06ad292fbf46b43e4e4ae58e6795a47b0f006ee8bc3997e9be57",
        "degree-40 source root identity changed");
    require(CORE::Expr(9)/16<root.expression() && root.expression()<CORE::Expr(3)/4,
        "degree-40 source root left its independently frozen isolator");
    auto corrupt=certificate;corrupt.lower="1/1";corrupt.upper="17/16";
    rejects([&]{CertifiedAlgebraicReal::reload(corrupt);});
    requireExactAlgebraicEnvironment();
    std::cout<<"PASS actual degree-40 arc root, strict isolator and corrupt interval rejection\n";
}
void tests() {
    AlgebraicKernel kernel;
    const AlgebraicPolynomial polynomial(-2,0,1);
    std::vector<AlgebraicReal> roots;
    kernel.solve_1_object()(polynomial,true,std::back_inserter(roots));
    require(roots.size()==2,"sqrt(2) root count");
    const auto positive=CertifiedAlgebraicReal::fromRoot(roots[1]);
    const auto negative=CertifiedAlgebraicReal::fromRoot(roots[0]);
    const auto& x=positive.expression();
    require((x*x-2).sign()==0 && x.sign()>0,"root expression not exact");
    require((x+negative.expression()).sign()==0,"opposite root selected");
    require((x-CORE::sqrt(CORE::Expr(2))).sign()==0,"equal algebraic expressions differ");
    const auto copy=CertifiedAlgebraicReal::reload(positive.certificate());
    require((copy.expression()-x).sign()==0,"reload changed the root");
    require(copy.certificate().rootIdentity==positive.certificate().rootIdentity,
            "reload changed root identity");
    auto refined=positive.certificate();
    refined.lower="11/8";refined.upper="23/16";
    require(CertifiedAlgebraicReal::reload(refined).certificate().rootIdentity
            ==positive.certificate().rootIdentity,"isolation refinement changed feature identity");
    auto invalid=positive.certificate();
    invalid.lower="3/1";invalid.upper="4/1";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);}); // no roots: must reject before CORE.
    invalid.lower="-2/1";invalid.upper="2/1";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);}); // two roots
    invalid=positive.certificate();invalid.rootIndex=0;
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    invalid=positive.certificate();invalid.rootIdentity[0]^=1;
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    invalid=positive.certificate();invalid.lower="4/3";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);}); // not dyadic
    invalid=positive.certificate();invalid.lower="not-a-number";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    invalid=positive.certificate();invalid.coefficients={"0","0","0"};
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    invalid=positive.certificate();invalid.lower="2/1";invalid.upper="1/1";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    invalid=positive.certificate();invalid.lower="1/0";
    rejects([&]{CertifiedAlgebraicReal::reload(invalid);});
    const auto previousEscape=CORE::get_static_EscapePrec();
    CORE::get_static_EscapePrec()=CORE::extLong(64);
    try {
        rejects([&]{CertifiedAlgebraicReal::reload(positive.certificate());});
        rejects([&]{positive.expression();});
    }catch(...){CORE::get_static_EscapePrec()=previousEscape;throw;}
    CORE::get_static_EscapePrec()=previousEscape;

    const auto rational=kernel.construct_algebraic_real_1_object()(CORE::BigRat(1,3));
    const auto third=CertifiedAlgebraicReal::fromRoot(rational);
    require((3*third.expression()-1).sign()==0,"non-dyadic rational root rounded");
    // The interval root has even multiplicity in the original polynomial.
    // Native reconstruction must square-free only its internal handle, not
    // change the saved polynomial or reject a valid prior certificate.
    AlgebraicRootCertificate repeated;
    repeated.coefficients={"4","0","-4","0","1"};
    repeated.rootIndex=1;repeated.lower="1/1";repeated.upper="2/1";
    repeated.rootIdentity="efa23c9f29b8ac40d14f6cf5e99c6383ae853af01de4aec052ac3429236edecb";
    const auto twice=CertifiedAlgebraicReal::reload(repeated);
    require((twice.expression()*twice.expression()-2).sign()==0
        && twice.certificate().coefficients==repeated.coefficients,
        "even-multiplicity original root lost during exact reload");
    // Distinct rational roots below binary64 separation: exact arithmetic and
    // isolation must not fuse them. The test is independent of decimal output.
    const CORE::BigInt denominator("1267650600228229401496703205376"); // 2^100
    const AlgebraicPolynomial first(CORE::BigInt(-1),CORE::BigInt(1));
    const AlgebraicPolynomial second(-(denominator+1),denominator);
    roots.clear();kernel.solve_1_object()(first*second,true,std::back_inserter(roots));
    require(roots.size()==2,"near roots merged");
    const auto a=CertifiedAlgebraicReal::fromRoot(roots[0]);
    const auto b=CertifiedAlgebraicReal::fromRoot(roots[1]);
    require((b.expression()-a.expression()).sign()>0,"near-root arithmetic lost order");
    using K=CGAL::Cartesian<CORE::Expr>;
    const auto midpoint=(a.expression()+b.expression())/2;
    const K::Point_2 p(a.expression(),0),q(b.expression(),0),r(midpoint,CORE::Expr(1)/denominator);
    require(CGAL::orientation(p,q,r)==CGAL::LEFT_TURN,"exact Cartesian predicate lost sign");
    require(CGAL::orientation(p,q,K::Point_2(midpoint,0))==CGAL::COLLINEAR,
            "exact Cartesian predicate failed zero");
    std::cout<<"PASS certified roots, arithmetic, invalid isolators, reload and exact predicates\n";
}
}
#endif
int main(int argc,char** argv) {
    try {
#if __has_include("spin/detail/RankTwoAlgebraic.h")
        if(argc==2 && std::string(argv[1])=="--degree40")degree40();
        else {
#ifdef MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR
            authenticatedFactorSelection();
            scopedCompleteRootProof();
            ordinalComparisonIsolation();
            completeFactorRootSet();
            actualDegree144();
#endif
            nativeProofReuse();tests();degree40();
        }
        return 0;
#else
        throw std::runtime_error("MISSING_FEATURE: certified algebraic arithmetic adapter");
#endif
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
