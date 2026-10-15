// Real solver allocation regression: retained metric columns must not also
// be materialized as a dense n*(m+faces) matrix at every working-set step.
#include "spin2014/SlqpTangentQp.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
static bool measuring=false;
static std::size_t bytes=0;
extern "C" void* __real_malloc(std::size_t);
extern "C" void* __wrap_malloc(std::size_t n) {
    if(measuring)bytes+=n;
    return __real_malloc(n);
}
int main(){try{
    constexpr int n=2048,b=64;
    Eigen::VectorXd x=Eigen::VectorXd::Constant(n,.5),g=Eigen::VectorXd::Zero(n);
    for(int i=0;i<b;++i)g[i]=20.+i;
    Eigen::MatrixXd a=Eigen::MatrixXd::Zero(1,n);a(0,n-1)=1.;
    spin2014::TangentQpSettings settings;settings.maximumFaces=-1;
    measuring=true;
    const auto r=spin2014::solvePrimalBoxQp(x,g,a,Eigen::VectorXd::Zero(1),b,{},settings);
    measuring=false;
    if(!r.checked || r.workingSet.size()!=b)throw std::runtime_error("analytic QP not solved");
    for(int i=0;i<n;++i)
        if(std::abs(r.direction[i]-(i<b?-.5:0.))>1e-12)
            throw std::runtime_error("wrong analytic QP direction");
    constexpr std::size_t limit=24U*1024U*1024U;
    std::cout<<"requested_malloc_bytes="<<bytes<<" budget="<<limit<<'\n';
    if(!bytes || bytes>limit)throw std::runtime_error("redundant metric materialization allocation");
    // GNU/Linux allocator contract for this fixture, not RSS or paper speed.
    return 0;
}catch(const std::exception& e){measuring=false;std::cerr<<e.what()<<'\n';return 1;}}
