#include "spin/QuadricSourceMotion.h"
#include <gmpxx.h>
#include <iostream>
#include <set>
#include <stdexcept>

namespace {
void check(const spin::QuadricCoefficientAuthority10& authority,
           const Eigen::MatrixXd& design,const spin::QuadricTriangulationSettings& settings) {
    const auto source=spin::triangulateQuadric(authority,design,settings);
    const auto motion=spin::evaluateQuadricSourceMotion(authority,design,settings,source);
    if(!motion.defined)throw std::runtime_error(motion.diagnostic);
    if(motion.triangulation.rationalConstructionIdentity!=source.rationalConstructionIdentity ||
       motion.triangulation.exactVertices!=source.exactVertices)
        throw std::runtime_error("authority motion substituted a different Boolean operand");
    if(motion.coefficientJacobian.size()!=source.exactVertices.size() || motion.scientificHash.empty())
        throw std::runtime_error("authority source motion omitted actual vertices");
    for(std::size_t vertex=0;vertex<source.exactVertices.size();++vertex)for(int axis=0;axis<3;++axis) {
        mpq_class radialLow=0,radialHigh=0;
        for(int column=0;column<10;++column) {
            const mpq_class lo(motion.lower[vertex][axis][column]),hi(motion.upper[vertex][axis][column]);
            const mpq_class value(motion.coefficientJacobian[vertex][axis][column]);
            if(!(lo<=value && value<=hi) || hi-lo>mpq_class(1,mpz_class(1)<<384))
                throw std::runtime_error("source motion lacks bounded rational derivatives");
            const mpq_class q(authority.geometry()[column]);
            radialLow+=q*(q<0?hi:lo);radialHigh+=q*(q<0?lo:hi);
        }
        if(!(radialLow<=0 && radialHigh>=0))
            throw std::runtime_error("source motion broke normalization radial invariance");
    }
    // A cap is a mean of the ACTUAL loop, not a recomputed display ring.
    for(int center=source.interfaceVertexCount;center<source.closedBoundary.vertices.rows();++center) {
        std::set<int> ring;
        for(int face=source.interfaceFaceCount;face<source.closedBoundary.faces.rows();++face)
            for(int corner=0;corner<3;++corner)if(source.closedBoundary.faces(face,corner)==center)
                for(int k=0;k<3;++k)if(source.closedBoundary.faces(face,k)!=center)
                    ring.insert(source.closedBoundary.faces(face,k));
        if(ring.size()<3)throw std::runtime_error("source cap has incomplete topology");
        for(int axis=0;axis<3;++axis)for(int column=0;column<10;++column) {
            mpq_class mean=0;
            for(int vertex:ring)mean+=mpq_class(motion.coefficientJacobian[vertex][axis][column]);
            mean/=static_cast<unsigned long>(ring.size());
            if(mean!=mpq_class(motion.coefficientJacobian[center][axis][column]))
                throw std::runtime_error("authority cap differential did not retain the same exact loop");
        }
    }
    auto mutated=source;mutated.exactVertices[0][0]="0/1";
    if(spin::evaluateQuadricSourceMotion(authority,design,settings,mutated).defined)
        throw std::runtime_error("stale exact source accepted by authority motion");
    auto wrong=settings;++wrong.icosphereSubdivisions;
    if(spin::evaluateQuadricSourceMotion(authority,design,wrong,source).defined)
        throw std::runtime_error("old source used after fixed normal branch changed");
    std::cout<<"source motion vertices="<<source.exactVertices.size()<<" hash="<<motion.scientificHash<<'\n';
}
}

int main() {
    try {
        using Authority=spin::QuadricCoefficientAuthority10;
        Authority::RationalVector q;q.fill("0/1");q[0]="-1/2";q[4]=q[5]=q[6]="1/2";
        Authority::RationalVector direction;direction.fill("0/1");direction[0]="1/1";
        const auto base=Authority::fromNormalizedGeometry(q);
        const auto changed=base.geodesicUpdate(direction,"1/10000000000000000000000000");
        Eigen::MatrixXd design(8,3);
        for(int vertex=0;vertex<8;++vertex)for(int axis=0;axis<3;++axis)
            design(vertex,axis)=(vertex&(1<<axis))?.1:-.1;
        spin::QuadricTriangulationSettings settings;
        settings.icosphereSubdivisions=0;settings.maximumIcosphereSubdivisions=-1;
        settings.rationalSourceCoordinates=true;settings.compactRationalSourceCoordinates=true;
        const auto source=spin::triangulateQuadric(changed,design,settings);
        if(spin::evaluateQuadricSourceMotion(base,design,settings,source).defined)
            throw std::runtime_error("double-equal old coefficients reused authority source");
        check(changed,design,settings);
        settings.icosphereSubdivisions=1;settings.truncationRadius=8;
        for(const char* constant:{"-1/2","1/2"}) {
            q[0]=constant;q[6]="-1/2";
            check(Authority::fromNormalizedGeometry(q).geodesicUpdate(
                direction,"1/10000000000000000000000000"),design,settings);
        }
        std::cout<<"PASS: same-source authority finite/root/cap differential and stale rejection\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
