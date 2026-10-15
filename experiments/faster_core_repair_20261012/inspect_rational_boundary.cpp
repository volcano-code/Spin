// Diagnose the actual rational boundary; no Nef, repair, or binary64 conversion.
#include "verification/ExactBoundaryAudit.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if(argc!=2)throw std::invalid_argument("rational boundary input required");
    std::ifstream input(argv[1]);std::size_t count=0,faces=0;
    if(!(input>>count>>faces) || count==0 || faces==0 || count>2000000 || faces>4000000)
        throw std::invalid_argument("invalid rational boundary dimensions");
    std::vector<std::array<std::string,3>> points(count);
    for(auto& point:points)for(auto& coordinate:point)
        if(!(input>>coordinate))throw std::invalid_argument("missing rational coordinate");
    Eigen::MatrixXi triangles(faces,3);
    for(std::size_t f=0;f<faces;++f)for(int k=0;k<3;++k)
        if(!(input>>triangles(f,k)))throw std::invalid_argument("missing rational face index");
    const auto result=verification::auditExactTriangleBoundaryRational(points,triangles,
        verification::ExactBoundaryOrientation::AnyNonzero);
    const auto& r=result.facts();
    std::cout<<"{\"boundary_valid\":"<<(r.boundaryValid?"true":"false")
        <<",\"closed\":"<<(r.closed?"true":"false")
        <<",\"edge_manifold\":"<<(r.edgeManifold?"true":"false")
        <<",\"vertex_manifold\":"<<(r.vertexManifold?"true":"false")
        <<",\"self_intersection_free\":"<<(r.intersectionFree?"true":"false")
        <<",\"exact_pair_tests\":"<<r.exactPairTests<<",\"non_manifold_vertices\":[";
    for(std::size_t i=0;i<r.nonManifoldVertices.size();++i){if(i)std::cout<<',';std::cout<<r.nonManifoldVertices[i];}
    std::cout<<"],\"degenerate_faces\":[";
    for(std::size_t i=0;i<r.degenerateFaces.size();++i){if(i)std::cout<<',';std::cout<<r.degenerateFaces[i];}
    std::cout<<"],\"intersecting_pairs\":[";
    for(std::size_t i=0;i<r.intersectingPairs.size();++i) {
        if(i)std::cout<<',';std::cout<<'['<<r.intersectingPairs[i][0]<<','<<r.intersectingPairs[i][1]<<']';
    }
    std::cout<<"],\"geometry_admission\":false,\"optimization_admission\":false}\n";
}
