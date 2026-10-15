#include "spin/QuadricInverseGaussMap.h"
#include "spin/QuadricCoefficientAuthority.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        using A = spin::QuadricCoefficientAuthority10;
        A::RationalVector sphere; sphere.fill("0/1");
        sphere[0]="-1/2"; sphere[4]=sphere[5]=sphere[6]="1/2";
        A::RationalVector direction; direction.fill("0/1"); direction[0]="1/1";
        const auto changed=A::fromNormalizedGeometry(sphere).geodesicUpdate(
            direction,"1/10000000000000000000000000");
        const auto point=spin::QuadricInverseGaussMap::rationalPointFromAuthority(
            changed,Eigen::Vector3d::UnitX());
        if (!point.defined) throw std::runtime_error(point.diagnostic);
        // A positive constant-direction update decreases the sphere radius.
        // Binary64 coefficients still describe exactly the unit sphere.
        const mpq_class lo(point.lower[0]), hi(point.upper[0]), x(point.coordinates[0]);
        if (!(lo<=x && x<=hi && hi<1 && lo>mpq_class(999,1000)))
            throw std::runtime_error("source point lost the authoritative coefficient update");
        if (hi-lo>mpq_class(1,mpz_class(1)<<384))
            throw std::runtime_error("source enclosure width exceeded 2^-384");
        if (mpq_class(point.coordinates[1])!=0 || mpq_class(point.coordinates[2])!=0)
            throw std::runtime_error("axis-normal sphere point was not axis aligned");
        if (spin::QuadricInverseGaussMap::rationalPointFromAuthority(
                changed,Eigen::Vector3d::Zero()).defined)
            throw std::runtime_error("invalid normal was admitted");
        auto hyperboloid=sphere; hyperboloid[6]="-1/2";
        const auto hyper=A::fromNormalizedGeometry(hyperboloid).geodesicUpdate(
            direction,"1/10000000000000000000000000");
        const auto hp=spin::QuadricInverseGaussMap::rationalPointFromAuthority(
            hyper,Eigen::Vector3d(1,0,0));
        if (!hp.defined || mpq_class(hp.upper[0])>=1)
            throw std::runtime_error("hyperboloid source lost its coefficient authority");
        if (spin::QuadricInverseGaussMap::rationalPointFromAuthority(
                hyper,Eigen::Vector3d(0,0,1)).defined)
            throw std::runtime_error("inadmissible hyperboloid normal was admitted");
        auto rankTwo=sphere; rankTwo[6]="0/1";
        const auto unsupported=A::fromUnnormalizedGeometry(rankTwo);
        if (spin::QuadricInverseGaussMap::rationalPointFromAuthority(
                unsupported,Eigen::Vector3d(1,0,0)).defined)
            throw std::runtime_error("rank-two authority silently used a legacy fallback");
        std::cout<<"PASS: authority drives finite source point and exact enclosure\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
