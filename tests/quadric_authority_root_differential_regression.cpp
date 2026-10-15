#include "spin/TruncationRootCertificate.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

namespace {
void checkRootDifferences(const spin::QuadricCoefficientAuthority10& authority,
                          const Eigen::Vector3d& first,const Eigen::Vector3d& last) {
    const auto jacobian=spin::TruncationRootCertifier(authority,2.0).rationalDifferential(first,last);
    if(!jacobian.defined)throw std::runtime_error(jacobian.diagnostic);
    mpq_class maximum=0;
    for(int scale=0;scale<3;++scale)for(int column=0;column<10;++column) {
        const mpq_class h(1,10000*(1<<scale));
        auto plus=authority.geometry(),minus=authority.geometry();
        mpq_class p(plus[column]),m(minus[column]);p+=h;m-=h;
        plus[column]=p.get_str();minus[column]=m.get_str();
        const auto upper=spin::TruncationRootCertifier(
            spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(plus),2.0).rationalPoint(first,last);
        const auto lower=spin::TruncationRootCertifier(
            spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(minus),2.0).rationalPoint(first,last);
        if(!upper.defined || !lower.defined)
            throw std::runtime_error("root difference changed the regular unique-root branch");
        for(int axis=0;axis<3;++axis) {
            const mpq_class expected(jacobian.coefficientJacobian[axis][column]);
            const mpq_class derivative=(mpq_class(upper.coordinates[axis])-mpq_class(lower.coordinates[axis]))/(2*h);
            mpq_class error=derivative-expected;if(error<0)error=-error;
            mpq_class normalization=expected;if(normalization<0)normalization=-normalization;
            if(normalization<1)normalization=1;
            error/=normalization;if(error>maximum)maximum=error;
            if(error>mpq_class(1,1000000))
                throw std::runtime_error("root three-scale coefficient difference exceeds 1e-6");
        }
    }
    std::cout<<"root difference max normalized error="<<maximum.get_d()<<'\n';
}
}

int main() {
    try {
        using Authority=spin::QuadricCoefficientAuthority10;
        Authority::RationalVector q;q.fill("0/1");
        q[0]="-1/2";q[4]=q[5]="1/2";q[6]="-1/2";
        Authority::RationalVector direction;direction.fill("0/1");direction[0]="1/1";
        const auto changed=Authority::fromNormalizedGeometry(q).geodesicUpdate(
            direction,"1/10000000000000000000000000");
        const Eigen::Vector3d first(1,0,0),last(1,0,-0.8);
        const spin::TruncationRootCertifier certifier(changed,2.0);
        const auto result=certifier.rationalDifferential(first,last);
        if(!result.defined)throw std::runtime_error(result.diagnostic);
        const auto originalPoint=certifier.rationalPoint(first,last);
        if(result.point.coordinates!=originalPoint.coordinates ||
           result.point.lower!=originalPoint.lower || result.point.upper!=originalPoint.upper)
            throw std::runtime_error("derivative proof changed the actual frozen source point");
        // x^2=(r^2-c/a)/2, dx/dc=-1/(4*a*x). The positive constant
        // geodesic increases a enough that this magnitude becomes smaller
        // than 1/sqrt(10), even though all double coefficients are unchanged.
        const mpq_class lower(result.lower[0][0]);
        if(!(lower<0 && lower*lower<mpq_class(1,10)))
            throw std::runtime_error("truncation derivative lost authoritative sub-binary64 coefficients");
        const auto base=spin::TruncationRootCertifier(Authority::fromNormalizedGeometry(q),2.0)
            .rationalDifferential(first,last);
        if(!base.defined)throw std::runtime_error(base.diagnostic);
        // Independent constraints: x^2-z^2=1, x^2+z^2=4, gradient_y=0.
        const mpq_class linearLo(base.lower[0][1]),linearHi(base.upper[0][1]);
        if(!(linearLo<=mpq_class(-1,2) && linearHi>=mpq_class(-1,2)))
            throw std::runtime_error("root linear coefficient derivative is not -1/2");
        for(int column: {7,9}) {
            const mpq_class lo(base.lower[1][column]),hi(base.upper[1][column]);
            const mpq_class square=column==7 ? mpq_class(6) : mpq_class(10);
            if(!(lo<0 && hi<0 && hi*hi<=square && lo*lo>=square))
                throw std::runtime_error("root cross coefficient factor is not 2");
        }
        for(int axis=0;axis<3;++axis)for(int column=0;column<10;++column) {
            const mpq_class lo(result.lower[axis][column]),hi(result.upper[axis][column]);
            const mpq_class value(result.coefficientJacobian[axis][column]);
            if(!(lo<=value && value<=hi) || hi-lo>mpq_class(1,mpz_class(1)<<384))
                throw std::runtime_error("root differential lacks its certified high-precision enclosure");
        }
        if(certifier.rationalDifferential(first,first).defined)
            throw std::runtime_error("degenerate normal edge admitted by root derivative");
        if(spin::TruncationRootCertifier(changed,0.5).rationalDifferential(first,last).defined)
            throw std::runtime_error("edge without a crossing fabricated a root differential");
        checkRootDifferences(Authority::fromNormalizedGeometry(q),
            Eigen::Vector3d(1,0.1,0),Eigen::Vector3d(1,0.1,-0.8));
        std::cout<<"PASS: authority root differential, frozen point identity and exact enclosures\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
