#include "verification/ExactBoundaryAudit.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/intersections.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    using namespace verification;
    Eigen::MatrixXd v(4,3);
    v << 0,0,0, 1,0,0, 0,2,0, 0,0,3;
    Eigen::MatrixXi f(4,3);
    f << 0,2,1, 0,1,3, 0,3,2, 1,2,3;
    const auto tetra = auditExactTriangleBoundary(v,f,ExactBoundaryOrientation::OutwardPositive);
    require(tetra.facts().boundaryValid, "literal tetrahedron must pass");
    require(CGAL::Gmpq(tetra.facts().componentSignedVolumesExact[0].c_str()) == 1,
            "literal exact tetrahedron volume is 1");
    for (double delta : {1e-7,1e-8}) {
        Eigen::MatrixXd crossing(5,3);
        crossing << 0,0,0, 1,0,0, 0,1,0, 1,1,-1, -1,-1,1+delta;
        Eigen::MatrixXi faces(2,3); faces << 0,1,2, 0,3,4;
        const auto r = auditExactTriangleBoundary(crossing,faces,ExactBoundaryOrientation::AnyNonzero);
        require(r.facts().intersectionCompleted && r.facts().intersectingPairs ==
                std::vector<std::array<int,2>>{{0,1}}, "short indexed crossing lost");
    }
    // Exhaustive oracle: no shared indices, so every exact geometric contact
    // is forbidden. This independently calls CGAL do_intersect on ALL pairs;
    // it never uses the production BVH or its candidate generation.
    using K = CGAL::Simple_cartesian<CGAL::Gmpq>;
    for (int scenario=0; scenario<3; ++scenario) {
        constexpr int count=72;
        Eigen::MatrixXd points(count*3,3);
        Eigen::MatrixXi faces(count,3);
        std::vector<K::Triangle_3> tris;
        for (int i=0;i<count;++i) {
            const double x = scenario==0 ? 5*i : (i*17%23)-11;
            const double y = scenario==0 ? 0 : (i*7%19)-9;
            const double z = scenario==2 ? (i%3)-1 : 0;
            points.row(3*i) << x,y,z;
            points.row(3*i+1) << x+2,y,z;
            points.row(3*i+2) << x,y+2,z+(scenario==2?2:0);
            faces.row(i) << 3*i,3*i+1,3*i+2;
            auto p = [&](int j) { return K::Point_3(points(j,0),points(j,1),points(j,2)); };
            tris.emplace_back(p(3*i),p(3*i+1),p(3*i+2));
        }
        std::vector<std::array<int,2>> expected;
        for(int i=0;i<count;++i) for(int j=i+1;j<count;++j)
            if(CGAL::do_intersect(tris[i],tris[j])) expected.push_back({i,j});
        const auto r = auditExactTriangleBoundary(points,faces,ExactBoundaryOrientation::AnyNonzero);
        require(r.facts().intersectionCompleted, "incomplete broad phase");
        require(r.facts().intersectingPairs == expected, "BVH differs from exhaustive exact oracle");
        if(scenario==0) require(r.facts().exactPairTests==0, "separated triangles should be pruned");
        std::cout << "oracle scenario=" << scenario << " all_pairs=2556 tested="
                  << r.facts().exactPairTests << " offending=" << expected.size() << '\n';
    }
    bool rejected=false;
    try { auditExactTriangleBoundary(Eigen::MatrixXd(2,2), f, ExactBoundaryOrientation::AnyNonzero); }
    catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"native malformed shape not rejected");
    std::cout << "exact boundary native regressions passed\n";
}
