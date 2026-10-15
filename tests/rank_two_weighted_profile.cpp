// Fresh-process diagnostic of an independently captured ACTUAL affine call.
// This does not change the production sampler or prove whole coverage.
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include <CGAL/Polynomial_traits_d.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace spin::detail;
using Q=CGAL::Gmpq;
using P=AlgebraicPolynomial;
using P2=CGAL::Polynomial<P>;
template<class T> std::string text(const T& x){std::ostringstream out;out<<x;return out.str();}
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
long long ns(clockid_t clock){timespec t{};require(clock_gettime(clock,&t)==0,"clock unavailable");return t.tv_sec*1000000000LL+t.tv_nsec;}
struct Span{
    const char* label;long long wall=ns(CLOCK_MONOTONIC),cpu=ns(CLOCK_THREAD_CPUTIME_ID);
    ~Span(){std::cout<<"WEIGHTED_PHASE "<<label<<' '<<ns(CLOCK_MONOTONIC)-wall<<' '<<ns(CLOCK_THREAD_CPUTIME_ID)-cpu<<'\n';}
};
AlgebraicRootCertificate read(std::istream& in){
    AlgebraicRootCertificate c;int n=0;in>>n;require(n>=2&&n<=257,"captured degree invalid");
    for(int i=0;i<n;++i){std::string value;in>>value;require(!value.empty(),"truncated coefficient");c.coefficients.push_back(value);}
    in>>c.rootIndex>>c.lower>>c.upper>>c.rootIdentity;require(bool(in),"truncated captured root");return c;
}
P polynomial(const AlgebraicRootCertificate& c){std::vector<CORE::BigInt> v;for(const auto& x:c.coefficients)v.emplace_back(x);return P(v.begin(),v.end());}
void emit(const char* label,const P& p){std::cout<<"WEIGHTED_POLYNOMIAL "<<label<<' '<<p.degree()+1;for(int i=0;i<=p.degree();++i)std::cout<<' '<<p[i];std::cout<<'\n';}
void phase(const char* value,void*) noexcept {std::cout<<"WEIGHTED_ROOT_PHASE "<<value<<'\n';}
struct ActualObserver {
    bool enumerating=false;int authenticatedUnions=0;
    int factorizations=0,rootSets=0,proofReuses=0,proofMismatches=0;
    static void observe(const char* value,void* context) noexcept {
        auto& self=*static_cast<ActualObserver*>(context);const std::string name(value);
        if(name=="weighted-coordinate:solve-all-roots")self.enumerating=true;
        if(name=="weighted-coordinate:select-isolator")self.enumerating=false;
        if(self.enumerating&&name=="algebraic-root:factor-root-set-72")++self.authenticatedUnions;
        if(name=="algebraic-root:native-factor")++self.factorizations;
        if(name.rfind("algebraic-root:factor-root-set-",0)==0)++self.rootSets;
        if(name=="algebraic-root:complete-proof-reused")++self.proofReuses;
        if(name=="algebraic-root:complete-proof-mismatch")++self.proofMismatches;
        phase(value,nullptr);
    }
};
#ifdef RANK_TWO_EXPECT_AUTHENTICATED_WEIGHTED_ROOTS
void actualConjugateMismatch(){
    // Hand-derived roots of x^4-10x^2+1 are ±sqrt(3)±sqrt(2).
    // Averaging the TWO positive roots gives sqrt(3). The complete
    // all-pair mean polynomial is x(x^2-2)(x^2-3)(x^4-10x^2+1).
    // Existing conjugate checks exclude the quartic and x factors, so the
    // caller freezes (x^2-2)(x^2-3), ordinal3, NOT the original degree9.
    AlgebraicKernel kernel;const P p(1,0,-10,0,1);
    const auto a=CertifiedAlgebraicReal::fromRoot(kernel.construct_algebraic_real_1_object()(
        p,CORE::BigRat(1,4),CORE::BigRat(1,2)));
    const auto b=CertifiedAlgebraicReal::fromRoot(kernel.construct_algebraic_real_1_object()(
        p,CORE::BigRat(3),CORE::BigRat(13,4)));
    ActualObserver observer;
    const auto result=[&]{ScopedRankTwoDiagnostics observe(&ActualObserver::observe,&observer);
        return sampleRankTwoAlgebraicCoordinate(a,b,Q(1,2),512);}();
    const auto& c=result.certificate();
    require(c.coefficients==std::vector<std::string>({"6","0","-5","0","1"})&&c.rootIndex==3,
        "weighted conjugate checks changed original reduced polynomial or ordinal");
    const auto expected=kernel.construct_algebraic_real_1_object()(P(-3,0,1),CORE::BigRat(1),CORE::BigRat(2));
    require(rankTwoAlgebraicHandle(result).compare(expected)==CGAL::EQUAL,
        "weighted conjugate mismatch changed exact positive sqrt(3)");
    require(observer.factorizations==2&&observer.rootSets==2
        &&observer.proofReuses==0&&observer.proofMismatches==1,
        "weighted conjugate rewrite did not invalidate complete-resultant proof");
    std::cout<<"WEIGHTED_CONJUGATE_MISMATCH "<<c.rootIdentity<<" 4 3 "
        <<observer.factorizations<<' '<<observer.rootSets<<' '<<observer.proofReuses<<' '
        <<observer.proofMismatches<<'\n';
}
#endif
int main(int argc,char** argv){try{
    require(argc==2,"actual affine input required");std::ifstream in(argv[1]);std::string schema,fraction;int bits=0;
    in>>schema>>fraction>>bits;require(schema=="RANK_TWO_ACTUAL_AFFINE_V1"&&bits==512,"wrong actual affine scope");
    const auto ac=read(in),bc=read(in),expected=read(in);std::string extra;require(!(in>>extra),"extra actual input");
    const Q t(fraction);require(t>0&&t<1,"captured interior fraction invalid");
    const auto a=CertifiedAlgebraicReal::reload(ac,bits),b=CertifiedAlgebraicReal::reload(bc,bits);
    const P p=polynomial(ac),q=polynomial(bc);
    require(p.degree()*q.degree()==128,"not the frozen actual degree-product128 input");
    const auto raw=[&]{Span timing{"original-raw-resultant"};
        const CORE::BigInt B(text(t.numerator())),C(text(t.denominator())),A=C-B;
        const P2 z=CGAL::shift(P2(1),1,0),u=CGAL::shift(P2(1),1,1);P2 pu(0),homogeneous(0);
        for(int i=p.degree();i>=0;--i)pu=pu*u+P2(p[i]);
        const P2 v=P2(C)*z-P2(A)*u;
        for(int i=0;i<=q.degree();++i)homogeneous+=P2(q[i]*CGAL::ipower(B,q.degree()-i))*CGAL::ipower(v,i);
        return CGAL::Polynomial_traits_d<P2>::Resultant()(pu,homogeneous);}();
    const auto sf=[&]{Span timing{"original-square-free"};return CGAL::make_square_free(raw);}();
    const auto sfCandidate=[&]{Span timing{"authenticated-square-free"};ScopedRankTwoDiagnostics observer(&phase,nullptr);return rankTwoSquareFreeNativePolynomial(raw);}();
    require(sf==sfCandidate,"complete square-free coefficients differ");emit("raw",raw);emit("square-free",sf);
    AlgebraicKernel kernel;
    const auto original=[&]{Span timing{"original-all-real-roots"};std::vector<AlgebraicReal> v;kernel.solve_1_object()(sf,true,std::back_inserter(v));return v;}();
    const auto factored=[&]{Span timing{"authenticated-all-real-roots"};ScopedRankTwoDiagnostics observer(&phase,nullptr);return rankTwoCompleteNativeRealRoots(sf);}();
    require(original.size()==factored.size()&&!original.empty(),"complete real-root count differs");
    for(std::size_t i=0;i<original.size();++i)require(original[i].compare(factored[i])==CGAL::EQUAL,"complete root order/value differs");
    std::cout<<"WEIGHTED_ROOT_SET "<<sf.degree()<<' '<<original.size()<<'\n';
    // Actual production call remains authoritative. Candidate factor roots
    // are NOT submitted as substituted polynomial/ordinal certificates.
    ActualObserver actualObserver;
    const auto actual=[&]{Span timing{"actual-affine-caller"};ScopedRankTwoDiagnostics observer(&ActualObserver::observe,&actualObserver);return sampleRankTwoAlgebraicCoordinate(a,b,t,bits);}();
    const auto& c=actual.certificate();
    require(c.rootIdentity==expected.rootIdentity&&c.rootIndex==expected.rootIndex&&c.coefficients==expected.coefficients,
        "actual affine caller differs from original captured root identity");
    require(rankTwoAlgebraicHandle(actual).compare(rankTwoAlgebraicHandle(CertifiedAlgebraicReal::reload(expected,bits)))==CGAL::EQUAL,
        "actual affine root differs from captured parent image");
#ifdef RANK_TWO_EXPECT_AUTHENTICATED_WEIGHTED_ROOTS
    std::cout<<"WEIGHTED_ACTUAL_PROOF_COUNTS "<<actualObserver.factorizations<<' '
        <<actualObserver.rootSets<<' '<<actualObserver.proofReuses<<' '
        <<actualObserver.proofMismatches<<'\n';
    require(actualObserver.authenticatedUnions==1,
        "actual weighted caller still lacks authenticated complete root enumeration");
    require(actualObserver.factorizations==1&&actualObserver.rootSets==1
        &&actualObserver.proofReuses==1&&actualObserver.proofMismatches==0,
        "actual weighted factory repeats certified complete root proof");
    actualConjugateMismatch();
#endif
    std::cout<<"WEIGHTED_REPLAY_SCIENCE "<<ac.rootIdentity<<' '<<bc.rootIdentity<<' '<<fraction<<' '<<c.rootIdentity<<' '<<c.rootIndex<<'\n';
    std::cout<<"WEIGHTED_RESULT_ISOLATOR "<<c.lower<<' '<<c.upper<<'\n';
    std::cout<<"WEIGHTED_REPLAY_SCOPE ACTUAL_CAPTURE_AND_TEST_ONLY_ROOT_COMPARISON_NOT_WHOLE_COVERAGE\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
