#include "spin2014/SlqpTangentQp.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <random>
#include <set>
#include <algorithm>
// Link-time counting spy: every call delegates to the real production metric.
// No replacement values and no test methods added to production classes.
static std::size_t inverseCalls=0;
extern "C" Eigen::VectorXd realInverse(const spin2014::LbfgsHistory*,const Eigen::VectorXd&)
    asm("__real__ZNK8spin201412LbfgsHistory12applyInverseERKN5Eigen6MatrixIdLin1ELi1ELi0ELin1ELi1EEE");
extern "C" Eigen::VectorXd countedInverse(const spin2014::LbfgsHistory*,const Eigen::VectorXd&)
    asm("__wrap__ZNK8spin201412LbfgsHistory12applyInverseERKN5Eigen6MatrixIdLin1ELi1ELi0ELin1ELi1EEE");
extern "C" Eigen::VectorXd countedInverse(const spin2014::LbfgsHistory* h,const Eigen::VectorXd& v){
    ++inverseCalls;return realInverse(h,v);
}
void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
void dump(const spin2014::TangentQpResult& r){
    std::cout<<std::hexfloat<<r.attempted<<' '<<r.checked<<' '<<r.status<<' '<<r.equalityResidual<<' '<<r.stationarity<<'\n';
    for(Eigen::Index i=0;i<r.direction.size();++i)std::cout<<r.direction[i]<<' ';std::cout<<'\n';
    for(auto f:r.workingSet)std::cout<<f<<' ';std::cout<<'\n';
    for(Eigen::Index i=0;i<r.multipliers.size();++i)std::cout<<r.multipliers[i]<<' ';std::cout<<'\n';
    for(const auto& t:r.trace){std::cout<<t.action<<' '<<t.changedIndex<<' '<<t.equalityResidual<<' '<<t.stationarity<<' ';
        for(auto f:t.workingSet)std::cout<<f<<',';std::cout<<'\n';}
}
int main(int argc,char** argv){try{
    const bool reference=argc==2 && std::strcmp(argv[1],"--reference")==0;
    require(argc==1 || reference,"invalid option");
    for(int equalities=0;equalities<=1;++equalities)for(int sign:{1,-1}) {
        const int n=8+equalities;
        Eigen::VectorXd x=Eigen::VectorXd::Constant(n,.5),g=Eigen::VectorXd::Zero(n);
        for(int i=0;i<8;++i)g[i]=sign*(20.+i);
        Eigen::MatrixXd a=Eigen::MatrixXd::Ones(equalities,n);
        Eigen::VectorXd c=Eigen::VectorXd::Zero(equalities);
        inverseCalls=0;
        auto r=spin2014::solvePrimalBoxQp(x,g,a,c,8,{});
        require(r.checked,"analytic box not solved");
        for(int i=0;i<8;++i)require(std::abs(r.direction[i]+sign*.5)<1e-9,"wrong bound solution");
        if(equalities)require(std::abs(r.direction[8]-sign*4.)<1e-9,"wrong free coordinate");
        if(reference)dump(r);
        else {
            std::cout<<"equations="<<equalities<<" sign="<<sign<<" inverse_calls="<<inverseCalls<<'\n';
            std::size_t candidateEvaluations=0;
            for(const auto& step:r.trace)if(step.action=="BLOCKING_FACE"||
                step.action=="RELEASE_FACE"||step.action=="CHECKED")++candidateEvaluations;
            // Retain the original metric-column cache contract and add one
            // direct inverse-metric application per checked working-set
            // candidate instead of summing every cached face column over n.
            require(inverseCalls==std::size_t(1+equalities+8+candidateEvaluations),
                    "candidate inverse application count or retained-column cache changed");
        }
    }
    // Consecutive calls change dimensions, constraints and actual L-BFGS
    // histories. Include coupled metrics where faces leave and later re-enter.
    std::mt19937 generator(1729);
    auto sample=[&](){return double(generator())/double(generator.max())*2.-1.;};
    std::size_t releases=0,reentries=0;
    for(int trial=0;trial<256;++trial) {
        const int n=6+trial%3,m=trial%2;
        Eigen::VectorXd x=Eigen::VectorXd::Constant(n,.5),g(n);
        Eigen::MatrixXd a(m,n),v(n,n);
        for(int i=0;i<n;++i)g[i]=10.*sample();
        for(int j=0;j<m;++j)for(int i=0;i<n;++i)a(j,i)=sample();
        for(int j=0;j<n;++j)for(int i=0;i<n;++i)v(j,i)=sample();
        const Eigen::MatrixXd b=v.transpose()*v+Eigen::MatrixXd::Identity(n,n);
        spin2014::LbfgsHistory h;h.capacity=4;
        for(int k=0;k<4;++k) {
            Eigen::VectorXd s(n);for(int i=0;i<n;++i)s[i]=sample();
            require(h.update(s,b*s),"positive-curvature fixture rejected");
        }
        inverseCalls=0;
        const auto r=spin2014::solvePrimalBoxQp(x,g,a,Eigen::VectorXd::Zero(m),n-1,h);
        require(r.checked,"deterministic coupled QP failed");
        std::set<int> seen;
        std::size_t entries=0;
        for(const auto& t:r.trace) {
            if(t.action=="RELEASE_FACE")++releases;
            if(t.action=="BLOCKING_FACE") {
                ++entries;if(!seen.insert(t.changedIndex).second)++reentries;
            }
        }
        if(reference)dump(r);
        else {
            const auto candidateEvaluations=std::count_if(r.trace.begin(),r.trace.end(),[](const auto& step){
                return step.action=="BLOCKING_FACE"||step.action=="RELEASE_FACE"||step.action=="CHECKED";
            });
            require(inverseCalls==std::size_t(1+m+entries+candidateEvaluations),
                    "candidate inverse application count or retained-column cache changed");
        }
    }
    require(releases>0 && reentries>0,"release/re-entry paths not covered");
    if(!reference)std::cout<<"coupled trials=256 releases="<<releases<<" reentries="<<reentries<<'\n';
    // Phase I remains uncached: preserve its trace and budget accounting.
    Eigen::MatrixXd a(1,2);a<<1.,1.;
    Eigen::VectorXd c(1);c<<.1;
    const auto phase=spin2014::solvePrimalBoxQp(Eigen::Vector2d(.8,0.),Eigen::Vector2d(-3.,1.),a,c,1,{});
    require(phase.checked && std::abs(phase.direction[0]-.2)<1e-9 &&
        std::abs(phase.direction[1]+.3)<1e-9,"nonzero equality phase I changed");
    require(!phase.trace.empty() && phase.trace.front().action.find("PHASE_I_")==0,"phase I not exercised");
    if(reference)dump(phase);
    spin2014::TangentQpSettings settings;settings.maximumIterations=1;
    const auto limited=spin2014::solvePrimalBoxQp(Eigen::Vector2d(.8,0.),Eigen::Vector2d(-3.,1.),a,c,1,{},settings);
    require(!limited.checked,"phase I ignored iteration budget");
    if(reference)dump(limited);
    if(!reference)std::cout<<"column reuse checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
