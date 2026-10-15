#include "spin2014/SlqpTangentQp.h"
#include "../src/spin2014/SchurFactorCache.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
int main(){try{
    std::cout<<std::hexfloat;
    for(double delta:{1e-4,0.}) {
        Eigen::MatrixXd a=Eigen::MatrixXd::Zero(2,3);
        a(0,0)=a(1,0)=1.;a(1,1)=delta;
        Eigen::MatrixXd s=a*a.transpose();
        Eigen::VectorXd scale=s.diagonal().array().rsqrt();
        s=scale.asDiagonal()*s*scale.asDiagonal();
        spin2014::detail::SchurFactorCache cache;
        Eigen::VectorXd out;
        if(cache.trySolve(s,Eigen::VectorXd::Zero(2),{-1,-2},out))
            throw std::runtime_error("fixture must reject fast path");
        Eigen::Vector3d x(.5,.5,0.),g(1.,0.,1.);
        const auto r=spin2014::solvePrimalBoxQp(x,g,a,Eigen::VectorXd::Zero(2),2,{});
        if(delta>0 && (!r.checked || (r.direction-Eigen::Vector3d(0.,0.,-1.)).norm()>1e-12))
            throw std::runtime_error("COD fallback did not solve near-singular fixture");
        if(delta>0 && r.multipliers.norm()<.5)
            throw std::runtime_error("fallback must solve a nonzero Schur RHS");
        if(delta==0 && (r.checked || r.trace.empty() || r.trace.back().action!="RANK_DEFICIENT"))
            throw std::runtime_error("rank-deficient reference behavior changed");
        std::cout<<delta<<' '<<r.attempted<<' '<<r.checked<<' '<<r.status<<' '
                 <<r.equalityResidual<<' '<<r.stationarity<<'\n';
        for(auto v:r.workingSet)std::cout<<v<<' ';std::cout<<'\n';
        for(Eigen::Index i=0;i<r.direction.size();++i)std::cout<<r.direction[i]<<' ';std::cout<<'\n';
        for(Eigen::Index i=0;i<r.multipliers.size();++i)std::cout<<r.multipliers[i]<<' ';std::cout<<'\n';
        for(const auto& t:r.trace){
            std::cout<<t.action<<' '<<t.changedIndex<<' '<<t.equalityResidual<<' '<<t.stationarity;
            for(auto v:t.workingSet)std::cout<<' '<<v;std::cout<<'\n';
        }
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
