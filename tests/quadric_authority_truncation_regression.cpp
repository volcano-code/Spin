#include "spin/TruncationRootCertificate.h"
#include "spin/QuadricCoefficientAuthority.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        using A=spin::QuadricCoefficientAuthority10;
        A::RationalVector q; q.fill("0/1");
        q[0]="-1/2"; q[4]=q[5]="1/2"; q[6]="-1/2";
        A::RationalVector direction; direction.fill("0/1"); direction[0]="1/1";
        const auto authority=A::fromNormalizedGeometry(q).geodesicUpdate(
            direction,"1/10000000000000000000000000");
        const spin::TruncationRootCertifier certifier(authority,2.0);
        const auto point=certifier.rationalPoint(Eigen::Vector3d(1,0,0),Eigen::Vector3d(1,0,-0.8));
        if (!point.defined) throw std::runtime_error(point.diagnostic);
        const mpq_class lo(point.lower[0]),hi(point.upper[0]),x(point.coordinates[0]);
        // On radius 2, x²-z²=k and x²+z²=4 imply x²=(4+k)/2.
        // Positive constant update makes k<1, despite identical double q.
        if (!(lo>1 && lo<=x && x<=hi && hi*hi<mpq_class(5,2)))
            throw std::runtime_error("truncation root lost its coefficient authority");
        for (unsigned i=0;i<3;++i)
            if (mpq_class(point.upper[i])-mpq_class(point.lower[i])>
                    mpq_class(1,mpz_class(1)<<384))
                throw std::runtime_error("truncation coordinate enclosure exceeded 2^-384");
        if (certifier.rationalPoint(Eigen::Vector3d(1,0,0),Eigen::Vector3d(1,0,0)).defined)
            throw std::runtime_error("degenerate normal edge was admitted");
        std::cout<<"PASS: authority drives exact truncation root and enclosure\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
