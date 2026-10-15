#pragma once
#include <cstddef>

namespace spin::detail {
/** First true index in a pure false-then-true predicate, or count if absent.
 * The caller must establish monotonicity; this is not a feasibility oracle. */
template<class Predicate>
std::size_t firstMonotoneTrue(std::size_t count,Predicate predicate) {
    std::size_t lo=0,hi=count;
    while(lo<hi) {
        const auto mid=lo+(hi-lo)/2;
        if(predicate(mid))hi=mid;
        else lo=mid+1;
    }
    return lo;
}
}
