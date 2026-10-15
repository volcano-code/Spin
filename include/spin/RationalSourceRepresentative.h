#pragma once
#include <CGAL/Gmpq.h>
#include <stdexcept>
namespace spin {
/** Enclosed-dyadic-1024 v1; never rounds through binary64. */
inline CGAL::Gmpq rationalSourceRepresentative(
    const CGAL::Gmpq& lower,const CGAL::Gmpq& upper,bool compact)
{
    using Q=CGAL::Gmpq;
    if(upper<lower)throw std::invalid_argument("reversed rational source interval");
    const Q midpoint=(lower+upper)/2;
    if(!compact)return midpoint;
    CGAL::Gmpz denominator(1),rounded;
    mpz_mul_2exp(denominator.mpz(),denominator.mpz(),1024);
    const Q scaled=midpoint*Q(denominator)+Q(1,2);
    mpz_fdiv_q(rounded.mpz(),scaled.numerator().mpz(),scaled.denominator().mpz());
    const Q candidate(rounded,denominator);
    return lower<=candidate && candidate<=upper ? candidate : midpoint;
}
} // namespace spin
