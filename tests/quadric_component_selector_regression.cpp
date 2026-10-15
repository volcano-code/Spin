#include "spin/QuadricMomentEvaluator.h"
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>

// Synthetic disconnected, closed outward components. These are selector
// contract tests, not claims that this mesh is a quadric triangulation.
int main() {
    try {
        spin::QuadricTriangulation t;
        t.paperProjectiveHyperboloid=true;
        t.closedBoundary.vertices.resize(8,3);
        t.closedBoundary.faces.resize(8,3);
        t.faceIsCap=Eigen::VectorXi::Zero(8);
        Eigen::Matrix<double,4,3> v;
        v<<1,1,1,1,-1,-1,-1,1,-1,-1,-1,1;
        Eigen::Matrix<int,4,3> f;
        f<<0,1,2,0,3,1,0,2,3,1,3,2;
        for(int i=0;i<4;++i) {
            Eigen::Vector3d a=v.row(f(i,0)),b=v.row(f(i,1)),c=v.row(f(i,2));
            if((b-a).cross(c-a).dot(a)<0)std::swap(f(i,1),f(i,2));
        }
        for(int k=0;k<2;++k) {
            t.closedBoundary.vertices.middleRows(4*k,4)=0.1*v;
            t.closedBoundary.vertices.block(4*k,2,4,1).array()+=(k?3.:-3.);
            t.closedBoundary.faces.middleRows(4*k,4)=f.array()+4*k;
        }
        spin::QuadricCoefficients10::Vector q;
        q<<-1,0,0,0,0,0,1,0,0,0;
        if(!spin::quadricSolidInteriorIsPositiveRegion(spin::QuadricCoefficients10::fromUnnormalized(q),t))
            throw std::runtime_error("unanimous component witnesses rejected");
        for(bool conflict:{false,true}) {
            q.setZero();q(conflict?3:1)=1.;
            bool rejected=false;
            try { (void)spin::quadricSolidInteriorIsPositiveRegion(spin::QuadricCoefficients10::fromUnnormalized(q),t); }
            catch(const std::runtime_error&) { rejected=true; }
            if(!rejected)throw std::runtime_error(conflict?"conflicting witnesses accepted":"zero-margin witnesses accepted");
        }
        std::cout<<"component selector positive/conflict/zero-margin checks passed\n";
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
