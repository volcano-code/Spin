#pragma once
#include <CGAL/Gmpq.h>
#include <stdexcept>

namespace spin::detail {
/** Exact sign of u+v*sqrt(s) on the strictly finite positive-root branch.
 * Squaring alone loses the branch: equal nonzero signs cannot cancel.
 * No floating conversion or approximate square root is used. */
inline int positiveRadicalSumSign(const CGAL::Gmpq& u,
                                  const CGAL::Gmpq& v,
                                  const CGAL::Gmpq& s) {
    if(s<=0)throw std::invalid_argument("positive radical sign requires s>0");
    const auto sign=[](const CGAL::Gmpq& x){return x>0?1:(x<0?-1:0);};
    const int us=sign(u),vs=sign(v);
    if(us==0)return vs;
    if(vs==0 || us==vs)return us;
    return us*sign(u*u-v*v*s);
}
}
