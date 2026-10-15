#pragma once
#include <cstdint>
#include <vector>
namespace spin::detail {
// Internal modular proof checker, not a persisted or caller-trusted proof.
// The prime is fixed; degree preservation and Bezout must BOTH hold.
inline constexpr std::uint64_t rankTwoSquareFreePrime=65537;
inline constexpr bool rankTwoSquareFreePrimeIsPrime(){
    for(std::uint64_t d=2;d*d<=rankTwoSquareFreePrime;++d)
        if(rankTwoSquareFreePrime%d==0)return false;
    return true;
}
static_assert(rankTwoSquareFreePrimeIsPrime(),"square-free certificate needs a prime field");
inline bool rankTwoVerifySquareFreeModulo(const std::vector<std::uint64_t>& q,
    const std::vector<std::uint64_t>& s,const std::vector<std::uint64_t>& t){
    // Reduction must preserve the original degree. A Bezout identity in a
    // degree-dropping image does NOT prove square-freeness over Q.
    if(q.size()<2||q.size()>1025||q.back()==0||s.size()>q.size()||t.size()>q.size())return false;
    for(const auto* coefficients:{&q,&s,&t})
        for(auto v:*coefficients)if(v>=rankTwoSquareFreePrime)return false;
    std::vector<std::uint64_t> derivative(q.size()-1),sum(2*q.size()-1,0);
    for(std::size_t i=1;i<q.size();++i)derivative[i-1]=(i*q[i])%rankTwoSquareFreePrime;
    const auto accumulate=[&](const auto& a,const auto& b){
        for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j)
            sum[i+j]=(sum[i+j]+a[i]*b[j])%rankTwoSquareFreePrime;
    };
    // This uses only bounded ordinary integer arithmetic, independently of
    // the native polynomial library that supplied S,T. Products < 2^33.
    accumulate(s,q);accumulate(t,derivative);
    if(sum[0]!=1)return false;
    for(std::size_t i=1;i<sum.size();++i)if(sum[i]!=0)return false;
    return true;
}
} // namespace spin::detail
