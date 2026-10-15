#define main frozen_regression_main
#include M1_FROZEN_REGRESSION
#undef main
#include "spin/CanonicalFullRankHollowing.h"
int main(){try {
 auto r=request();const auto es=settings();auto invalid=r;invalid.materialDensity=-1;
 throwsInvalid([&]{spin::solveFullRankQuadricHollowing(invalid,es,{fixture()});});
 throwsInvalid([&]{spin::solveCanonicalFullRankQuadricHollowing(invalid,es,{fixture()});});
 V coefficients;coefficients<<.2,0,0,-1,0,0,0,0,0,0;
 const auto plane=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
 const auto legacy=spin::evaluateFullRankMaterialMoments(r,plane,es);
 const auto canonical=spin::evaluateCanonicalFullRankMaterialMoments(r,plane,es);
 require(!legacy.coefficientJacobian && !canonical.coefficientJacobian,"plane derivative must stay unavailable");
 require(legacy.geometryHash==canonical.geometryHash && legacy.totalMomentsTarget.values()==canonical.totalMomentsTarget.values(),"canonical path changed unsupported-stratum native values");
 // With absent J there is no changed quadrature or result data to accidentally
 // separate the domains. Omitting the new policy wrapper must fail this check.
 require(legacy.scientificHash!=canonical.scientificHash,"canonical scientific policy collided with the unchanged legacy result");
 std::cout<<"[PASS] legacy and canonical APIs coexist with distinct scientific policy\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
