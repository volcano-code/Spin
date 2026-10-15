#include "../src/verification/ExactBoundaryDistance.h"
#include <CGAL/Gmpq.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace verification::detail;
void check(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
    try {
        Eigen::MatrixXd v(9,3);
        v << 0,0,0, 1,0,0, 0,1,0,
             0,0,.125, 1,0,.125, 0,1,.125,
             .25,.25,-1, .25,.25,1, .75,.25,0;
        Eigen::MatrixXi f(3,3); f << 0,1,2, 3,4,5, 6,7,8;
        auto found=exactNonincidentNearPairs(v,f,.125);
        check(found.pairs.size()==3,"inclusive distance/crossing detection");
        check(CGAL::Gmpq(found.pairs[0].squaredExact.c_str())==CGAL::Gmpq(1,64),"analytic parallel gap");
        check(found.pairs[0].parallel && !found.pairs[1].parallel,"parallel classification");
        check(CGAL::Gmpq(found.pairs[1].squaredExact.c_str())==CGAL::Gmpq(0),"edge/face crossing");
        check(exactNonincidentNearPairs(v,f,std::nextafter(.125,0.)).pairs.size()==2,"strictly outside threshold");
        Eigen::MatrixXi incident(2,3); incident << 0,1,2, 0,3,4;
        auto excluded=exactNonincidentNearPairs(v,incident,.125);
        check(excluded.pairs.empty() && excluded.incidentExcluded==1,"shared-index exclusion visible");
        bool rejected=false;
        try { exactNonincidentNearPairs(v,f,std::numeric_limits<double>::quiet_NaN()); }
        catch (const std::invalid_argument&) { rejected=true; }
        check(rejected,"NaN threshold rejection");
        incident.row(0) << 0,0,1; rejected=false;
        try { exactNonincidentNearPairs(v,incident,.125); }
        catch (const std::invalid_argument&) { rejected=true; }
        check(rejected,"degenerate triangle rejection");
        // More than four triangles exercises internal BVH nodes. Exhaustive
        // traversal oracle uses original two-boundary routine one pair at a time.
        Eigen::MatrixXd many(36,3); Eigen::MatrixXi faces(12,3);
        for (int i=0;i<12;++i) {
            many.row(3*i) << (i%4)*.2, (i/4)*.2, (i%3)*.03;
            many.row(3*i+1)=many.row(3*i)+Eigen::RowVector3d(.3,0,.01*(i%2));
            many.row(3*i+2)=many.row(3*i)+Eigen::RowVector3d(0,.3,0);
            faces.row(i) << 3*i,3*i+1,3*i+2;
        }
        auto actual=exactNonincidentNearPairs(many,faces,.075);
        std::size_t next=0;
        for (int i=0;i<12;++i) for (int j=i+1;j<12;++j) {
            Eigen::MatrixXi a=faces.middleRows(i,1), b=faces.middleRows(j,1);
            auto reference=exactBoundaryDistance(many,a,many,b,.075);
            if (CGAL::Gmpq(reference.squaredExact.c_str())<=CGAL::Gmpq(.075)*CGAL::Gmpq(.075)) {
                check(next<actual.pairs.size(),"BVH missing pair");
                check(actual.pairs[next].faces==std::array<std::size_t,2>{std::size_t(i),std::size_t(j)},"BVH order or pair mismatch");
                check(actual.pairs[next++].squaredExact==reference.squaredExact,"BVH distance mismatch");
            }
        }
        check(next==actual.pairs.size(),"BVH spurious pair");
        std::cout << "exact near-boundary regression passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
