#include "spin2014/SlqpTangentQp.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#ifdef BASELINE_BOX
#define solvePrimalBoxQp solveFiniteBoxQp
#endif
int main(int argc,char** argv){try{
    if(argc!=2 && argc!=3)throw std::runtime_error("fixture required");
    std::ifstream in(argv[1]);int n,m,b,k;in>>n>>m>>b>>k;
    if(n!=45||m!=5||b!=44||k!=1)throw std::runtime_error("wrong frozen failure fixture");
    Eigen::VectorXd x(n),g(n),c(m);Eigen::MatrixXd a(m,n);
    for(int i=0;i<n;++i)in>>x[i];for(int i=0;i<n;++i)in>>g[i];for(int j=0;j<m;++j)in>>c[j];
    for(int j=0;j<m;++j)for(int i=0;i<n;++i)in>>a(j,i);
    spin2014::LbfgsHistory h;h.capacity=8;
    for(int j=0;j<k;++j){Eigen::VectorXd s(n),y(n);for(int i=0;i<n;++i)in>>s[i];for(int i=0;i<n;++i)in>>y[i];
        if(!h.update(s,y))throw std::runtime_error("invalid frozen history");}
    if(!in)throw std::runtime_error("truncated fixture");
    const auto r=spin2014::solvePrimalBoxQp(x,g,a,c,b,h);
    if(!r.checked)throw std::runtime_error("frozen feasible QP still UNKNOWN: "+(r.trace.empty()?r.status:r.trace.back().action));
    if((a*r.direction+c).lpNorm<Eigen::Infinity>()>1e-9 ||
       (r.direction.head(b).array() < -x.head(b).array()-1e-10).any() ||
       (r.direction.head(b).array() > 1-x.head(b).array()+1e-10).any())throw std::runtime_error("invalid full-box direction");
    auto repeat=spin2014::solvePrimalBoxQp(x,g,a,c,b,h);
    if(r.direction!=repeat.direction||r.workingSet!=repeat.workingSet)throw std::runtime_error("repeat mismatch");
    if(argc==3){
        if(std::string(argv[2])!="--dump")throw std::runtime_error("invalid option");
        std::cout<<std::setprecision(17)<<"direction";for(Eigen::Index i=0;i<r.direction.size();++i)std::cout<<' '<<r.direction[i];std::cout<<'\n';
        std::cout<<"multipliers";for(Eigen::Index i=0;i<r.multipliers.size();++i)std::cout<<' '<<r.multipliers[i];std::cout<<'\n';
        std::cout<<"working";for(auto v:r.workingSet)std::cout<<' '<<v;std::cout<<'\n';
        return 0;
    }
    std::cout<<"frozen QP checked, trace "<<r.trace.size()<<", residual "<<r.stationarity<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
