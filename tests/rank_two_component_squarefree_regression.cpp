// Regression: reconstruct the SAME complete square-free eliminant from ALL
// synchronous components without repeating whole-polynomial decomposition.
// Bad modular primes must still use the original exact fallback.
#include "spin/detail/RankTwoAlgebraic.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include <CGAL/Polynomial_traits_d.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace spin::detail;
using P=AlgebraicBivariatePolynomial;
using PT=CGAL::Polynomial_traits_d<P>;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
struct Counts{int component=0,whole=0,fallback=0;};
void observe(const char* phase,void* context) noexcept {
    auto& c=*static_cast<Counts*>(context);
    if(std::strcmp(phase,"native-fiber:component-square-free-certified")==0)++c.component;
    if(std::strcmp(phase,"native-fiber:square-free-begin")==0)++c.whole;
    if(std::strcmp(phase,"native-fiber:component-square-free-original-fallback")==0)++c.fallback;
}
std::map<std::string,std::string> records;
void evaluate(const std::string& label,const P& a,const P& b,bool eligible){
    const auto raw=rankTwoCompleteNativeResultant(a,b);
    const auto expected=CGAL::make_square_free(raw); // original, independent operation
    AlgebraicPolynomial actual;Counts counts;
    const auto roots=[&]{ScopedRankTwoDiagnostics scope(observe,&counts);
        return RankTwoNativeRootSet::fromResultant(a,b,actual);}();
    require(actual==expected,"complete original polynomial changed");
    const auto control=RankTwoNativeRootSet::fromPolynomial(expected);
    require(roots.roots().size()==control.roots().size(),"complete root count changed");
    for(std::size_t i=0;i<roots.roots().size();++i)
        require(roots.roots()[i].compare(control.roots()[i])==CGAL::EQUAL,"complete root order changed");
    if(eligible){
        require(counts.component==1,"MISSING_FEATURE: complete component square-free proof not consumed");
        require(counts.whole==0,"repeated whole polynomial square-free decomposition");
        require(counts.fallback==0,"ordinary complete component unexpectedly fell back");
    }else{
        require(counts.component==0,"bad modular prime incorrectly accepted");
        require(counts.fallback==1&&counts.whole==1,"original exact fallback missing");
    }
    std::ostringstream out;out<<actual.degree()+1;
    for(int i=0;i<=actual.degree();++i)out<<'\t'<<actual[i];
    out<<'\t'<<roots.roots().size();
    // Certify and PUBLICLY reload the actual first/last roots. Interior roots
    // were compared against the complete independent control set above.
    for(std::size_t i=0;i<roots.roots().size();++i){
        if(i!=0&&i+1!=roots.roots().size())continue;
        const auto value=CertifiedAlgebraicReal::fromRoot(roots.roots()[i],roots);
        const auto loaded=CertifiedAlgebraicReal::reload(value.certificate());
        require(rankTwoAlgebraicHandle(loaded).compare(roots.roots()[i])==CGAL::EQUAL,
            "public reload no longer matches original root");
        const auto& c=value.certificate();out<<'\t'<<c.rootIdentity<<'\t'<<c.rootIndex;
    }
    records.emplace(label,out.str());
    std::cout<<"PASS "<<label<<" roots="<<roots.roots().size()<<" fast="<<counts.component
             <<" whole="<<counts.whole<<" fallback="<<counts.fallback<<'\n';
}
AlgebraicRootCertificate readRoot(std::istream& in){
    AlgebraicRootCertificate c;int n;in>>n;require(n>=2&&n<=257,"saved degree invalid");
    for(int i=0;i<n;++i){std::string value;in>>value;c.coefficients.push_back(value);}
    in>>c.rootIndex>>c.lower>>c.upper>>c.rootIdentity;require(bool(in),"saved root truncated");return c;
}
P readCurve(const std::string& s){
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);
    P out(0);std::istringstream in(s);std::string term;
    while(std::getline(in,term,';')){
        auto comma=term.find(','),colon=term.find(':');require(comma!=std::string::npos&&colon!=std::string::npos,"bad support");
        auto i=std::stoi(term.substr(0,comma)),j=std::stoi(term.substr(comma+1,colon-comma-1));
        require(i>=0&&j>=0&&i+j<=4,"support degree invalid");
        out+=P(CORE::BigInt(term.substr(colon+1)))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }return out;
}
int main(int argc,char** argv){try{
    require(argc==3,"saved input and science output required");requireExactAlgebraicEnvironment();
    const P x=CGAL::shift(P(1),1,1),y=CGAL::shift(P(1),1,0);
    const P odd=-3*CGAL::ipower(x-1,3)*(x*x+1);
    evaluate("signed-right",x-y,odd,true);
    evaluate("signed-left",odd,x-y,true);
    evaluate("complex-only",x*x-y,(x*x+1)*(x*x+2),true);
    evaluate("unequal-repeated",x-y,-7*CGAL::ipower(x*x-2,2)*CGAL::ipower(x*x+1,3),true);
    evaluate("bad-prime",x-y,x*x-65537,false);
    evaluate("degree-loss",x-y,65537*x*x-1,false);
    std::ifstream in(argv[1]);std::string schema,arc,support,extra;int branch;
    in>>schema>>arc>>branch>>support;require(schema=="RANK_TWO_ACTUAL_FIBER_V1"&&arc.size()==64&&branch==0,"saved fixture mismatch");
    const auto u=readRoot(in),v=readRoot(in);require(!(in>>extra),"extra saved data");
    (void)CertifiedAlgebraicReal::reload(u);(void)CertifiedAlgebraicReal::reload(v);
    std::vector<CORE::BigInt> c;for(const auto& value:u.coefficients)c.emplace_back(value);
    evaluate("actual-U72-Y144",PT::Swap()(readCurve(support),0,1),
        PT::Swap()(P(AlgebraicPolynomial(c.begin(),c.end())),0,1),true);
    // Zero/constant elimination is NOT silently minted as a root proof.
    for(const auto& operands:std::vector<std::pair<P,P>>{{P(0),x-y},{P(1),x-y}}){
        bool rejected=false;AlgebraicPolynomial p;
        try{(void)RankTwoNativeRootSet::fromResultant(operands.first,operands.second,p);}
        catch(const spin::RankTwoSourceError&){rejected=true;}
        require(rejected,"zero/constant result minted a complete root proof");
    }
    std::ofstream out(argv[2]);for(const auto& item:records)out<<item.first<<'\t'<<item.second<<'\n';
    out.flush();require(bool(out),"science output write failed");
    std::cout<<"PASS zero-constant-rejection\nSCOPE COMPLETE_ARITHMETIC_ONLY; NO_FULL_CARRIER\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
