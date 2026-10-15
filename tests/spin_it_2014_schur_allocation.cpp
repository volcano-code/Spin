// Resource regression: explicit zero-filled box constraint matrices must not
// consume the allocation budget. The linker spy delegates every allocation to
// libc; this is the real QP and a hand-derived analytic solution, not a mock.
#include "spin2014/SlqpTangentQp.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

static bool measuring=false;
static std::size_t allocated=0;
extern "C" void* __real_malloc(std::size_t);
extern "C" void* __wrap_malloc(std::size_t size) {
    if(measuring)allocated+=size;
    return __real_malloc(size);
}
int main(){try {
    constexpr int n=2048,bounded=64;
    Eigen::VectorXd x=Eigen::VectorXd::Constant(n,.5),g=Eigen::VectorXd::Zero(n);
    for(int i=0;i<bounded;++i)g[i]=20.+i;
    const Eigen::MatrixXd a(0,n);const Eigen::VectorXd c(0);
    spin2014::TangentQpSettings settings;settings.maximumFaces=64;
    measuring=true;
    const auto r=spin2014::solvePrimalBoxQp(x,g,a,c,bounded,{},settings);
    measuring=false;
    if(!r.checked || r.workingSet.size()!=bounded)throw std::runtime_error("analytic QP not solved");
    for(int i=0;i<n;++i)
        if(std::abs(r.direction[i]-(i<bounded?-.5:0.))>1e-12)
            throw std::runtime_error("wrong analytic direction");
    std::cout<<"requested_malloc_bytes="<<allocated<<" budget="<<56*1024*1024<<'\n';
    // Fixed Linux/GNU-linker test budget, not a universal allocator/RSS claim.
    // Dense C alone costs >32 MiB over the growing 0..64 face sets.
    if(allocated==0)throw std::runtime_error("allocation instrumentation did not observe the real solver");
    if(allocated>56U*1024U*1024U)throw std::runtime_error("redundant dense working-set allocation budget exceeded");
    return 0;
}catch(const std::exception& error){measuring=false;std::cerr<<error.what()<<'\n';return 1;}}
