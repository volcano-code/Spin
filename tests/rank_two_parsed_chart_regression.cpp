#include "spin/detail/RankTwoStripMemo.h"
#include "spin/RankTwoSource.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
namespace {
using Q=CGAL::Gmpq;using P=spin::detail::RankTwoStripPolynomial;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F>void rejected(F f){
    try{f();}catch(const spin::RankTwoSourceError& e){
        require(e.code()=="INPUT_IDENTITY_MISMATCH","wrong invalid-input category");return;
    }throw std::runtime_error("malformed chart became cached authority");
}
template<class M>auto check(M& memo,const char* path,int)
    ->decltype(memo.parsedChart(std::string{},0),void()){
    const std::string raw="1,0:1/2;0,1:-3/7;2,1:5/11;0,0:2;0,0:-1;";
    P original{};original[0]=1;original[5]=Q(1,2);original[1]=Q(-3,7);original[11]=Q(5,11);
    P transpose{};transpose[0]=1;transpose[1]=Q(1,2);transpose[5]=Q(-3,7);transpose[7]=Q(5,11);
    const P& first=memo.parsedChart(raw,0);
    require(first==original && memo.parsedChart(raw,1)==transpose,"wrong exact chart or cross term");
    for(int n=0;n<500;++n){
        require(&memo.parsedChart(raw,0)==&first,"identical immutable context reparsed");
        require(memo.parsedChart("0,0:"+std::to_string(n+2)+";",n%2)[0]==n+2,
            "changed coefficient/axis reused stale result");
    }
    require(first==original,"interning other factors invalidated stored coefficients");
    P opposite{};opposite[0]=-1;opposite[1]=-1;
    require(memo.parsedChart("0,0:-1;0,1:-1;",0)==opposite,"material sign lost");
    for(int n=0;n<2;++n){
        rejected([&]{memo.parsedChart("0,0:1;",2);});
        rejected([&]{memo.parsedChart("5,0:1;",0);});
        rejected([&]{memo.parsedChart("0,1;",0);});
    }
    require(memo.parsedChart("0,1:1;",0)[1]==1,"failed parse poisoned later context");
    const auto& p=memo.parsedChart("0,1:1;2,0:-1;",0);
    const auto proof=memo.certify(p,{Q(-1),Q(1),Q(-1),Q(2)},{},Q(1));
    require(proof.derivativeLowerBound==1 && proof.slopeUpperBound==2
        && proof.curvatureUpperBound==2 && proof.chordDistanceUpperBound==1,
        "parsed exact coefficients changed downstream strip proof");
    M fresh;require(fresh.parsedChart(raw,0)==original,"fresh instance acquired wrong state");
    std::ofstream out(path);require(bool(out),"cannot save exact chart evidence");
    out<<"rank-two-parsed-chart-regression-v1\n";
    for(int axis=0;axis<2;++axis){out<<"FACTOR "<<axis<<' '<<std::quoted(raw);
        for(const auto& q:memo.parsedChart(raw,axis))out<<' '<<q;out<<'\n';}
    out<<"PROOF "<<proof.inputIdentity<<" 1 2 2 1\nfull_carrier_completion NOT_GRANTED\n";
    require(bool(out),"evidence write failed");
    std::cout<<"PASS exact parsing, chart identity/lifetime, changed coefficients, invalid input and downstream proof\n";
}
template<class M>void check(M&,const char*,long){
    throw std::runtime_error("MISSING_FEATURE: scoped immutable factor/chart parsing reuse");
}
}
int main(int argc,char** argv){try{
    require(argc==2,"evidence path required");spin::detail::RankTwoRegularStripMemo memo;
    check(memo,argv[1],0);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
