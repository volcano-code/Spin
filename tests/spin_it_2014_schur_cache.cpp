#include "../src/spin2014/SchurFactorCache.h"
#include <Eigen/QR>
#include <iostream>
#include <stdexcept>
#include <limits>
void require(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}
int main(){try{
    spin2014::detail::SchurFactorCache cache;
    Eigen::VectorXd solution;
    Eigen::MatrixXd a(2,2);a<<2,1,1,2;
    Eigen::VectorXd b(2);b<<3,3;
    require(cache.trySolve(a,b,{10,20},solution),"well-conditioned solve must use guarded path");
    require((solution-Eigen::VectorXd::Ones(2)).norm()<1e-12,"wrong solution");
    // Addition in middle of sorted identity order must not reuse wrong columns.
    a.resize(3,3);a<<2,.2,1,.2,3,.4,1,.4,2;
    b=a*Eigen::Vector3d(1,2,3);
    require(cache.trySolve(a,b,{10,15,20},solution),"append rejected");
    require((solution-Eigen::Vector3d(1,2,3)).norm()<1e-12,"permuted append solve");
    Eigen::MatrixXd small(2,2);small<<3,.4,.4,2;
    require(cache.trySolve(small,Eigen::Vector2d(3.8,4.4),{15,20},solution),"release rejected");
    require((solution-Eigen::Vector2d(1,2)).norm()<1e-12,"release wrong");
    require(cache.trySolve(a,b,{10,15,20},solution),"reentry rejected");
    a(0,0)=4;b=a*Eigen::Vector3d(3,2,1);
    require(cache.trySolve(a,b,{10,15,20},solution),"changed metric rejected");
    require((solution-Eigen::Vector3d(3,2,1)).norm()<1e-12,"stale factor reused");
    a=Eigen::MatrixXd::Ones(2,2);
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,2},solution),"singular entered fast path");
    a(1,1)+=1e-10;
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,2},solution),"near singular entered fast path");
    a<<1,2,2,1;
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,2},solution),"indefinite entered fast path");
    a<<1,0,.1,1;
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,2},solution),"asymmetric entered fast path");
    require(!cache.trySolve(a,Eigen::Vector2d(1,.1),{1,2},solution),"RHS-masked asymmetry entered fast path");
    a=Eigen::MatrixXd::Identity(2,2);
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,1},solution),"duplicate identity accepted");
    a(0,0)=std::numeric_limits<double>::quiet_NaN();
    require(!cache.trySolve(a,Eigen::Vector2d(1,1),{1,2},solution),"NaN accepted");
    // Long incremental sequence with deterministic dense coupling, changing RHS.
    cache.clear();
    for(int n=1;n<=180;++n){
        Eigen::VectorXd v(n),rhs(n),want(n);
        std::vector<int> keys;for(int i=0;i<n;++i){v[i]=.2*std::sin(i+1.);want[i]=std::cos(i*.3);keys.push_back(i);}
        const Eigen::MatrixXd m=Eigen::MatrixXd::Identity(n,n)+v*v.transpose();rhs=m*want;
        require(cache.trySolve(m,rhs,keys,solution),"incremental sequence rejected");
        require((m*solution-rhs).cwiseAbs().maxCoeff()<1e-11,"sequence residual");
        require((solution-want).cwiseAbs().maxCoeff()<1e-11,"sequence solution");
    }
    std::cout<<"guarded Schur cache checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
