#include "spin/QuadricInverseGaussMap.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

namespace {
void checkThreeScaleFiniteMap(const spin::QuadricCoefficientAuthority10& input,
                             const Eigen::Vector3d& normal) {
    const auto jacobian=spin::QuadricInverseGaussMap::rationalDifferentialFromAuthority(input,normal);
    if(!jacobian.defined)throw std::runtime_error(jacobian.diagnostic);
    mpq_class maxError=0;
    for(int scale=0;scale<3;++scale) {
        const mpq_class h(1,10000*(1<<scale));
        for(int column=0;column<10;++column) {
            auto plus=input.geometry(),minus=input.geometry();
            mpq_class positive(plus[column]),negative(minus[column]);positive+=h;negative-=h;
            plus[column]=positive.get_str();minus[column]=negative.get_str();
            // The actual source-map values independently pass through the
            // coefficient authority normalization and finite-point constructor.
            const auto p=spin::QuadricInverseGaussMap::rationalPointFromAuthority(
                spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(plus),normal);
            const auto m=spin::QuadricInverseGaussMap::rationalPointFromAuthority(
                spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(minus),normal);
            if(!p.defined || !m.defined)throw std::runtime_error("finite-map difference left admissible branch");
            for(int axis=0;axis<3;++axis) {
                const mpq_class difference=(mpq_class(p.coordinates[axis])-mpq_class(m.coordinates[axis]))/(2*h);
                const mpq_class expected(jacobian.coefficientJacobian[axis][column]);
                mpq_class error=difference-expected;if(error<0)error=-error;
                mpq_class normalization=expected;if(normalization<0)normalization=-normalization;
                if(normalization<1)normalization=1;
                error/=normalization;if(error>maxError)maxError=error;
                if(error>mpq_class(1,1000000))
                    throw std::runtime_error("non-diagonal finite-map three-scale difference exceeds 1e-6");
            }
        }
    }
    std::cout<<"finite-point difference max normalized error="<<maxError.get_d()<<'\n';
}
}

int main() {
    try {
        using Authority=spin::QuadricCoefficientAuthority10;
        Authority::RationalVector q; q.fill("0/1");
        q[0]="-1/2"; q[4]=q[5]=q[6]="1/2";
        Authority::RationalVector direction; direction.fill("0/1"); direction[0]="1/1";
        const auto changed=Authority::fromNormalizedGeometry(q).geodesicUpdate(
            direction,"1/10000000000000000000000000");
        const auto result=spin::QuadricInverseGaussMap::rationalDifferentialFromAuthority(
            changed,Eigen::Vector3d::UnitX());
        if(!result.defined)throw std::runtime_error(result.diagnostic);
        // x=sqrt(-c/a), dx/dc=-1/(2*a*x). The invisible constant
        // update increases this magnitude; a binary64 reconstruction returns -1.
        if(!(mpq_class(result.upper[0][0]) < -1))
            throw std::runtime_error("finite differential lost sub-binary64 coefficient state");
        for(int axis=0;axis<3;++axis)for(int coefficient=0;coefficient<10;++coefficient) {
            const mpq_class lo(result.lower[axis][coefficient]);
            const mpq_class hi(result.upper[axis][coefficient]);
            const mpq_class value(result.coefficientJacobian[axis][coefficient]);
            if(!(lo<=value && value<=hi) || hi-lo>mpq_class(1,mpz_class(1)<<384))
                throw std::runtime_error("finite differential enclosure is invalid");
        }
        // Independent literals for x=(1,0,0) of the unit sphere.
        const auto base=spin::QuadricInverseGaussMap::rationalDifferentialFromAuthority(
            Authority::fromNormalizedGeometry(q),Eigen::Vector3d::UnitX());
        if(!base.defined)throw std::runtime_error(base.diagnostic);
        const int selected[6][2]={{0,0},{0,1},{0,4},{1,2},{1,9},{2,8}};
        const int expected[6]={-1,-1,-1,-1,-2,-2};
        // Off-diagonal coefficients multiply 2xy/2xz. Keeping the fixed
        // normal gives dy/dAxy=-x/a=-2, not the diagonal factor -1.
        for(int i=0;i<6;++i)
            if(mpq_class(base.coefficientJacobian[selected[i][0]][selected[i][1]])!=expected[i])
                throw std::runtime_error("sphere coefficient or cross coefficient factor is wrong");
        if(spin::QuadricInverseGaussMap::rationalDifferentialFromAuthority(
            changed,Eigen::Vector3d::Zero()).defined)
            throw std::runtime_error("invalid normal admitted by authority derivative");
        auto rankTwo=q;rankTwo[6]="0/1";
        if(spin::QuadricInverseGaussMap::rationalDifferentialFromAuthority(
            Authority::fromUnnormalizedGeometry(rankTwo),Eigen::Vector3d::UnitX()).defined)
            throw std::runtime_error("authority derivative fabricated rank-two full Jacobian");
        Authority::RationalVector asymmetric{{"-2/1","1/4","-1/8","1/16",
            "3/2","1/1","3/4","1/16","-1/32","1/8"}};
        checkThreeScaleFiniteMap(Authority::fromUnnormalizedGeometry(asymmetric),Eigen::Vector3d(1,2,3));
        asymmetric[6]="-3/4";
        checkThreeScaleFiniteMap(Authority::fromUnnormalizedGeometry(asymmetric),Eigen::Vector3d(1,0.125,0.0625));
        std::cout<<"PASS: authoritative finite differential, cross factors and exact enclosures\n";
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
