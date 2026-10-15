#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoStripMemo.h")
#include "spin/detail/RankTwoStripMemo.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include <cstring>
#include <variant>
namespace {
using Q=CGAL::Gmpq;using P=spin::detail::RankTwoStripPolynomial;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
struct Counts {int hit=0,miss=0;};
template<class M>
auto attemptChart(M& m,const std::string& p,int axis,const spin::detail::RankTwoStripBox& b,
    const std::vector<std::string>& o,const Q& e,int)->decltype(m.tryCertifyChart(p,axis,b,o,e)){
    return m.tryCertifyChart(p,axis,b,o,e);
}
template<class M>
std::variant<spin::detail::RankTwoRegularStripCertificate,std::string> attemptChart(M&,
    const std::string&,int,const spin::detail::RankTwoStripBox&,
    const std::vector<std::string>&,const Q&,long){
    throw std::runtime_error("MISSING_FEATURE: typed exact strip mathematical rejection");
}
template<class F>std::string failure(F f,const char* code){
    try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"failure code changed");return e.witness();}
    throw std::runtime_error("invalid or unproved strip accepted");
}
// Break caught: omitting any source/chart/ordered-obstacle/query component
// aliases a different pure rational problem, or creates a fabricated result.
// ADL/SFINAE keeps the missing implementation an actual runtime RED.
template<class M>
auto savedChart(M& memo,const std::string& p,int axis,const spin::detail::RankTwoStripBox& box,
    const std::vector<std::string>& obstacles,const Q& error,int)
    ->decltype(memo.certifyChart(p,axis,box,obstacles,error)){
    return memo.certifyChart(p,axis,box,obstacles,error);
}
template<class M>
spin::detail::RankTwoRegularStripCertificate savedChart(M&,const std::string&,int,
    const spin::detail::RankTwoStripBox&,const std::vector<std::string>&,const Q&,long){
    throw std::runtime_error("MISSING_FEATURE: exact saved-chart strip context with full ordered input binding");
}
void savedChartTests(){using namespace spin::detail;
    Counts counts;ScopedRankTwoDiagnostics sink([](const char* s,void* p)noexcept{
        auto& c=*static_cast<Counts*>(p);
        c.hit+=std::strcmp(s,"regular-strip:memo-hit")==0;
        c.miss+=std::strcmp(s,"regular-strip:memo-miss")==0;
    },&counts);
    RankTwoRegularStripMemo memo;
    const RankTwoStripBox box{Q(-1),Q(1),Q(-1),Q(2)};
    const std::string p="0,1:1;2,0:-1;",safe="0,0:-3;0,1:1;",safe2="0,0:-4;0,1:1;";
    const auto call=[&](const std::string& poly,int axis,const RankTwoStripBox& b,
        const std::vector<std::string>& obs,const Q& e){return savedChart(memo,poly,axis,b,obs,e,0);};
    const auto first=call(p,0,box,{safe,safe2},Q(1));
    const auto repeated=call(p,0,box,{safe,safe2},Q(1));
    require(first.derivativeLowerBound==1 && first.slopeUpperBound==2
        && first.curvatureUpperBound==2 && first.chordDistanceUpperBound==1
        && first.obstacleSigns==std::vector<int>({-1,-1})
        && first.inputIdentity==repeated.inputIdentity && counts.hit==1 && counts.miss==1,
        "saved source context changed analytic parabola or repeated query");
    // Actual input equivalent under additive term decoding: may share ONLY
    // the identical fully interned rational problem, never a textual subset.
    const auto alias=call("0,1:2;0,1:-1;2,0:-1;",0,box,{safe,safe2},Q(1));
    require(alias.inputIdentity==first.inputIdentity && counts.hit==2 && counts.miss==1,
        "equivalent decoded context did not share the same full rational query");
    const auto reordered=call(p,0,box,{safe2,safe},Q(1));
    require(reordered.inputIdentity!=first.inputIdentity,"ordered obstacles lost from source identity");
    failure([&]{call(p,1,box,{safe,safe2},Q(1));},"ARRANGEMENT_UNRESOLVED");
    const auto rotated=call("1,0:1;0,2:-1;",1,box,{"0,0:-3;1,0:1;","0,0:-4;1,0:1;"},Q(1));
    require(rotated.inputIdentity==first.inputIdentity,"exact chart transposition changed the rational problem");
    const auto reversed=call("0,1:-7;2,0:7;",0,box,{safe,safe2},Q(1));
    require(reversed.yDerivativeSign==-1 && reversed.inputIdentity!=first.inputIdentity,
        "signed support scaling incorrectly reused source proof");
    failure([&]{call(p,0,box,{"0,0:-1/2;0,1:1;2,0:-1;"},Q(1));},"ARRANGEMENT_UNRESOLVED");
    const auto why=failure([&]{call(p,0,box,{safe,safe2},Q(999,1000));},"ARRANGEMENT_UNRESOLVED");
    require(failure([&]{call(p,0,box,{safe,safe2},Q(999,1000));},"ARRANGEMENT_UNRESOLVED")==why,
        "saved query failure witness changed");
    failure([&]{call(p,0,{Q(-1),Q(1),Q(2),Q(3)},{safe,safe2},Q(1));},"ARRANGEMENT_UNRESOLVED");
    for(int repeat=0;repeat<2;++repeat){
        failure([&]{call("0,0:0;",0,box,{},Q(1));},"INPUT_IDENTITY_MISMATCH");
        failure([&]{call(p,2,box,{},Q(1));},"INPUT_IDENTITY_MISMATCH");
        failure([&]{call("5,0:1;",0,box,{},Q(1));},"INPUT_IDENTITY_MISMATCH");
        failure([&]{call("malformed",0,box,{},Q(1));},"INPUT_IDENTITY_MISMATCH");
        failure([&]{call(p,0,box,{},Q(0));},"INPUT_IDENTITY_MISMATCH");
    }
    // Mutation/lifetime: the memo owns its key and coefficients, not caller
    // mutable strings. Returning to the old source must retrieve old math.
    std::vector<std::string> mutableObstacles{safe,safe2};
    call(p,0,box,mutableObstacles,Q(1));mutableObstacles[0]="0,0:-1/2;0,1:1;2,0:-1;";
    failure([&]{call(p,0,box,mutableObstacles,Q(1));},"ARRANGEMENT_UNRESOLVED");
    require(call(p,0,box,{safe,safe2},Q(1)).inputIdentity==first.inputIdentity,
        "caller mutation corrupted retained pure source context");
    RankTwoRegularStripMemo fresh;
    require(savedChart(fresh,p,0,box,{safe,safe2},Q(1),0).inputIdentity==first.inputIdentity,
        "new memo inherited wrong source context");
    std::cout<<"PASS actual saved-chart full input binding, axes, ordered obstacles, aliases and invalid inputs\n";
}
void tests(const char* path){using namespace spin::detail;
    Counts counts;ScopedRankTwoDiagnostics sink([](const char* s,void* p)noexcept{
        auto& c=*static_cast<Counts*>(p);
        c.hit+=std::strcmp(s,"regular-strip:memo-hit")==0;
        c.miss+=std::strcmp(s,"regular-strip:memo-miss")==0;
    },&counts);
    RankTwoRegularStripMemo memo;
    P p{};p[1]=1;p[10]=-1;const RankTwoStripBox box{Q(-1),Q(1),Q(-1),Q(2)};
    const auto first=memo.certify(p,box,{},Q(1));
    const auto second=memo.certify(p,box,{},Q(1));
    // Hand-derived y=x²: m=1, slope=2, curvature=2, bidirectional error=1.
    require(first.derivativeLowerBound==1 && second.slopeUpperBound==2
        && second.curvatureUpperBound==2 && second.chordDistanceUpperBound==1
        && first.inputIdentity==second.inputIdentity,"memo changed exact mathematical proof");
    require(counts.hit==1 && counts.miss==1,"identical strip was recomputed instead of reused");
    const auto tight=[&]{memo.certify(p,box,{},Q(999,1000));};
    const auto why=failure(tight,"ARRANGEMENT_UNRESOLVED");
    require(counts.miss==2 && counts.hit==1,"changed accuracy reused old success");
    require(failure(tight,"ARRANGEMENT_UNRESOLVED")==why && counts.hit==2,
        "reused failure lost its original diagnosis");
    P crossing=p;crossing[0]=Q(-1,2);
    failure([&]{memo.certify(p,box,{crossing},Q(1));},"ARRANGEMENT_UNRESOLVED");
    require(counts.miss==3,"changed obstacle reused incomplete obstacle set");
    P safe{};safe[1]=1;safe[0]=-3;
    const auto other=memo.certify(p,box,{safe},Q(1));
    require(other.obstacleSigns==std::vector<int>{-1} && counts.miss==4,"safe obstacle context misbound");
    P changed=p;changed[1]=0;changed[2]=1;
    failure([&]{memo.certify(changed,box,{},Q(1));},"ARRANGEMENT_UNRESOLVED");
    require(counts.miss==5,"changed support reused old derivative");
    failure([&]{memo.certify(p,{Q(-1),Q(1),Q(2),Q(3)},{},Q(1));},"ARRANGEMENT_UNRESOLVED");
    require(counts.miss==6,"changed box reused old bracketing");
    P secondSafe{};secondSafe[1]=1;secondSafe[0]=-4;
    const auto ordered=memo.certify(p,box,{safe,secondSafe},Q(1));
    const auto reordered=memo.certify(p,box,{secondSafe,safe},Q(1));
    require(ordered.obstacleSigns==std::vector<int>({-1,-1}) && ordered.inputIdentity!=reordered.inputIdentity
        && counts.miss==8,"ordered source obstacles were silently canonicalized");
    P reversed=p;for(auto& c:reversed)c*=-7;
    const auto reverse=memo.certify(reversed,box,{},Q(1));
    require(reverse.yDerivativeSign==-1 && reverse.chordDistanceUpperBound==1 && counts.miss==9,
        "support rescaling reused wrong certificate sign/identity");
    for(int k=0;k<2;++k)failure([&]{memo.certify(P{},box,{},Q(1));},"INPUT_IDENTITY_MISMATCH");
    require(counts.miss==11,"invalid input was persisted as a geometric rejection");
    failure([&]{memo.certify(p,box,{},Q(0));},"INPUT_IDENTITY_MISMATCH");
    require(counts.miss==12,"invalid zero accuracy reused an old proof");
    RankTwoRegularStripMemo fresh;
    const auto independent=fresh.certify(p,box,{},Q(1));
    require(counts.miss==13 && independent.inputIdentity==first.inputIdentity,"new memo inherited caller state");
    std::ofstream out(path);require(bool(out),"cannot save memo evidence");
    out<<"rank-two-strip-memo-regression-v1\nPARABOLA "<<first.inputIdentity<<" 1 2 2 1\n"
       <<"COUNTS "<<counts.hit<<' '<<counts.miss<<"\nFAILURE "<<why
       <<"\nfull_registry_proof NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    require(bool(out),"memo evidence write failed");
    savedChartTests();
    // Break caught: a rejected strip is silently promoted, changes witness,
    // uses a separate cache, or input errors become mathematical rejections.
    RankTwoRegularStripMemo typed;
    const std::string raw="0,1:1;2,0:-1;";
    const auto typedProof=attemptChart(typed,raw,0,box,{},Q(1),0);
    require(typedProof.index()==0 && std::get<0>(typedProof).inputIdentity==first.inputIdentity,
        "typed strip path changed exact accepted certificate");
    const int misses=counts.miss,hits=counts.hit;
    for(int i=0;i<2;++i){
        const auto rejected=attemptChart(typed,raw,0,box,{},Q(999,1000),0);
        require(rejected.index()==1 && std::get<1>(rejected)=="regular graph chord bound exceeds requested accuracy",
            "typed rejection lost mathematical cause or granted a certificate");
    }
    require(counts.miss==misses+1 && counts.hit==hits+1,"typed path did not share owned query store");
    require(failure([&]{typed.certifyChart(raw,0,box,{},Q(999,1000));},"ARRANGEMENT_UNRESOLVED")
        =="regular graph chord bound exceeds requested accuracy","legacy exception contract changed");
    for(int i=0;i<2;++i)
        failure([&]{attemptChart(typed,raw,0,box,{},Q(0),0);},"INPUT_IDENTITY_MISMATCH");
    failure([&]{attemptChart(typed,raw,2,box,{},Q(1),0);},"INPUT_IDENTITY_MISMATCH");
    failure([&]{attemptChart(typed,"0,0:0;",0,box,{},Q(1),0);},"INPUT_IDENTITY_MISMATCH");
    std::cout<<"PASS typed mathematical rejection and unchanged legacy/input error semantics\n";
    std::cout<<"PASS exact strip reuse, changed support/box/error/obstacles, signed scaling and invalid input\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoStripMemo.h")
    if(argc!=2)throw std::invalid_argument("memo evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact strip-only memo with full input identity and fresh caller scope");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
