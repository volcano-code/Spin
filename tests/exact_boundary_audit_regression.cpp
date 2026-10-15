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

// Literal rational contact cases. Open two-face fixtures test intersection
// classification only; they are deliberately not valid material boundaries.
static void rationalContacts() {
    using namespace verification;
    using Row = std::array<std::string,3>;
    struct Case { const char* name; std::array<Row,3> other;
                  std::array<int,3> indices; bool forbidden; };
    const std::vector<Case> cases = {
        {"coplanar disjoint, overlapping boxes", {{{"2/3","2/3","0/1"},
            {"2/3","1/3","0/1"},{"1/3","2/3","0/1"}}}, {3,4,5}, false},
        {"unindexed point contact", {{{"2/3","0/1","0/1"},
            {"1/1","0/1","0/1"},{"2/3","-1/3","0/1"}}}, {3,4,5}, true},
        {"unindexed edge contact", {{{"0/1","0/1","0/1"},
            {"2/3","0/1","0/1"},{"0/1","-2/3","0/1"}}}, {3,4,5}, true},
        {"coplanar containment", {{{"1/9","1/9","0/1"},
            {"2/9","1/9","0/1"},{"1/9","2/9","0/1"}}}, {3,4,5}, true},
        {"transverse segment", {{{"1/3","0/1","-1/3"},
            {"1/3","0/1","1/3"},{"1/3","1/3","0/1"}}}, {3,4,5}, true},
        {"near coplanar crossing", {{{"1/9","1/9","-1/1000000000000000000000000000000"},
            {"2/9","1/9","1/1000000000000000000000000000000"},
            {"1/9","2/9","1/1000000000000000000000000000000"}}}, {3,4,5}, true},
        {"indexed point only", {{{"0/1","0/1","0/1"},
            {"0/1","-1/3","1/3"},{"-1/3","0/1","1/3"}}}, {0,4,5}, false},
        {"indexed point with area overlap", {{{"0/1","0/1","0/1"},
            {"1/3","0/1","0/1"},{"0/1","1/3","0/1"}}}, {0,4,5}, true},
        {"indexed edge only", {{{"0/1","0/1","0/1"},
            {"2/3","0/1","0/1"},{"0/1","-2/3","0/1"}}}, {0,1,5}, false},
        {"indexed edge with area overlap", {{{"0/1","0/1","0/1"},
            {"2/3","0/1","0/1"},{"1/3","1/3","0/1"}}}, {0,1,5}, true},
        {"duplicate indexed face", {{{"0/1","0/1","0/1"},
            {"2/3","0/1","0/1"},{"0/1","2/3","0/1"}}}, {0,1,2}, true}
    };
    for (const auto& c : cases) for (int reverseA=0;reverseA<2;++reverseA)
        for (int reverseB=0;reverseB<2;++reverseB) {
            std::vector<Row> vertices = {{"0/1","0/1","0/1"},
                {"2/3","0/1","0/1"},{"0/1","2/3","0/1"}};
            vertices.insert(vertices.end(),c.other.begin(),c.other.end());
            Eigen::MatrixXi faces(2,3);
            faces << 0, reverseA?2:1, reverseA?1:2,
                c.indices[0], c.indices[reverseB?2:1], c.indices[reverseB?1:2];
            const auto result=auditExactTriangleBoundaryRational(
                vertices,faces,ExactBoundaryOrientation::AnyNonzero);
            require(result.facts().degenerateFaces.empty(),c.name);
            require(result.facts().intersectionCompleted,c.name);
            const std::vector<std::array<int,2>> expected = c.forbidden ?
                std::vector<std::array<int,2>>{{0,1}} : std::vector<std::array<int,2>>{};
            require(result.facts().intersectingPairs==expected,c.name);
        }
    std::cout << "44 literal rational contact/orientation cases passed\n";
}

int main() {
    using namespace verification;
    rationalContacts();
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
    // is forbidden. Construct intersections on ALL pairs, independently of
    // the production decision-only predicate and BVH candidate generation.
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
            if(CGAL::intersection(tris[i],tris[j])) expected.push_back({i,j});
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
