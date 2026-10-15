// Diagnostic harness only: clocks/records never reach the production adapter.
#include "spin/detail/RankTwoAlgebraic.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace spin::detail;
using Clock=std::chrono::steady_clock;
struct Phase {std::string name;long long ns;};
struct Observer {
    Clock::time_point previous=Clock::now();
    std::string current="native-interval-construction";
    std::vector<Phase> phases;
    void close(){const auto now=Clock::now();phases.push_back({current,
        std::chrono::duration_cast<std::chrono::nanoseconds>(now-previous).count()});previous=now;}
    static void observe(const char* phase,void* context) noexcept {
        auto& self=*static_cast<Observer*>(context);self.close();
        self.current=std::string(phase).substr(0,std::string(phase).find(' '));
    }
};
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(int argc,char** argv){try{
    require(argc==2,"profile input required");
    std::ifstream input(argv[1]);int cases=0;input>>cases;
    require(cases==9,"every captured degree class is required");
    AlgebraicKernel kernel;
    for(int k=0;k<cases;++k){
        int count=0;input>>count;require(count>=2&&count<=257,"profile polynomial degree");
        std::vector<CORE::BigInt> coefficients;std::vector<std::string> strings;
        for(int i=0;i<count;++i){std::string c;input>>c;require(!c.empty(),"truncated coefficient");
            strings.push_back(c);coefficients.emplace_back(c);}
        std::size_t ordinal=0;std::string lower,upper,identity;
        input>>ordinal>>lower>>upper>>identity;
        require(identity.size()==64,"profile root identity missing");
        const CORE::BigRat lo(lower.c_str()),hi(upper.c_str());
        const AlgebraicPolynomial original(coefficients.begin(),coefficients.end());
        require(kernel.is_square_free_1_object()(original),"captured original native polynomial is not square-free");
        const auto plainStart=Clock::now();
        const auto plainNative=kernel.construct_algebraic_real_1_object()(original,lo,hi);
        const auto plain=CertifiedAlgebraicReal::fromRoot(plainNative);
        const auto plainNs=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-plainStart).count();
        Observer observer;
        const auto measured=[&]{
            ScopedRankTwoDiagnostics observation(&Observer::observe,&observer);
            const AlgebraicPolynomial fresh(coefficients.begin(),coefficients.end());
            const auto handle=kernel.construct_algebraic_real_1_object()(fresh,lo,hi);
            return CertifiedAlgebraicReal::fromRoot(handle);
        }();
        observer.close();
        for(const auto* value:{&plain,&measured}){
            const auto& certificate=value->certificate();
            require(certificate.coefficients==strings && certificate.rootIndex==ordinal
                && certificate.rootIdentity==identity,"instrumentation changed original root/ordinal");
            require(CORE::Expr(lo)<value->expression() && value->expression()<CORE::Expr(hi),
                "instrumented/reference expression left original strict interval");
        }
        require(plain.certificate().lower==measured.certificate().lower
            && plain.certificate().upper==measured.certificate().upper,
            "instrumentation changed scientific isolator");
        requireExactAlgebraicEnvironment();
        std::cout<<"PROFILE_BEGIN "<<count-1<<' '<<ordinal<<' '<<identity<<'\n';
        std::cout<<"ORIGINAL_ISOLATOR "<<lower<<' '<<upper<<'\n'
            <<"FACTORY_ISOLATOR "<<measured.certificate().lower<<' '<<measured.certificate().upper<<'\n';
        for(const auto& phase:observer.phases)std::cout<<"PHASE "<<phase.name<<' '<<phase.ns<<'\n';
        std::cout<<"UNINSTRUMENTED_TOTAL_NS "<<plainNs<<'\n'<<"PROFILE_END "<<identity<<'\n';
        std::cout.flush();
    }
    std::string extra;require(!(input>>extra),"extra profile input");
    std::cout<<"PROFILE_SCIENCE_UNCHANGED 9\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
