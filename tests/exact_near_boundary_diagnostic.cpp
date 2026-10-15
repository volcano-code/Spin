#include "../src/verification/ExactBoundaryDistance.h"
#include <iostream>
#include <stdexcept>
int main() {
    try {
        std::size_t nv=0,nf=0; double threshold=0;
        if (!(std::cin>>nv>>nf>>threshold) || nv==0 || nf==0 || nv>10000000 || nf>10000000)
            throw std::invalid_argument("invalid diagnostic mesh header");
        Eigen::MatrixXd v(nv,3); Eigen::MatrixXi f(nf,3);
        for (std::size_t i=0;i<nv;++i) for (int j=0;j<3;++j)
            if (!(std::cin>>v(i,j))) throw std::invalid_argument("truncated vertices");
        for (std::size_t i=0;i<nf;++i) for (int j=0;j<3;++j)
            if (!(std::cin>>f(i,j))) throw std::invalid_argument("truncated faces");
        std::string extra; if (std::cin>>extra) throw std::invalid_argument("unexpected trailing input");
        const auto result=verification::detail::exactNonincidentNearPairs(v,f,threshold);
        std::cout << "{\"scope\":\"NON_SHARED_VERTEX_TRIANGLE_PAIRS_ONLY\","
                     "\"complete_wall_thickness\":\"UNKNOWN\",\"material_side\":\"UNKNOWN\","
                     "\"nodes_visited\":" << result.nodesVisited
                  << ",\"box_candidates\":" << result.boxCandidates
                  << ",\"incident_excluded\":" << result.incidentExcluded
                  << ",\"distances_evaluated\":" << result.distancesEvaluated << ",\"pairs\":[";
        bool first=true;
        for (const auto& pair:result.pairs) {
            if (!first) std::cout << ','; first=false;
            std::cout << "{\"faces\":[" << pair.faces[0] << ',' << pair.faces[1]
                      << "],\"squared_exact\":\"" << pair.squaredExact
                      << "\",\"parallel\":" << (pair.parallel ? "true" : "false") << '}';
        }
        std::cout << "]}\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
