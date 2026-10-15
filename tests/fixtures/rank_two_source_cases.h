#pragma once
#include "spin/RankTwoSource.h"

namespace rank_two_cases {
inline spin::RankTwoExactCoefficients disk() {return {0,0,0,1,-1,-1,0,0,0,0};}
inline spin::RankTwoExactCoefficients annulus() {return {2,0,0,1,-1,-1,0,0,0,0};}
inline spin::RankTwoExactCoefficients twoLobes() {return {-2,0,0,1,-1,1,0,0,0,0};}
inline CGAL::Gmpq radius() {return CGAL::Gmpq(3,2);}
// All three edges lie strictly outside the sphere projection disk.
inline std::array<std::array<CGAL::Gmpq,2>,3> enclosingTriangle() {
    return {{{{-4,-4}},{{4,-4}},{{0,6}}}};
}
} // namespace rank_two_cases
