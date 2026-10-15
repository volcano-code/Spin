#include "detail/RankTwoSettings.h"
#include "detail/FixedNormalSphere.h"
#include <cmath>
#include <limits>

namespace spin {
RankTwoArrangement analyzeRankTwoArrangement(
    const PreparedRankTwoSource& source,const QuadricTriangulationSettings& settings) {
    if(!settings.rankTwoArrangementSource || settings.fixedNormalParaboloid
        || settings.algebraicFiniteMap || settings.paperProjectiveHyperboloid
        || settings.maximumIcosphereSubdivisions!=-1
        || settings.icosphereSubdivisions<0
        || !std::isfinite(settings.truncationRadius) || settings.truncationRadius<=0)
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",
            "rank-two settings require explicit profile, radius and one fixed normal resolution");
    const auto& limits=settings.rankTwoLimits;
    // Icosphere counts: F=20*4^level, V=F/2+2, E=3F/2. Check before
    // allocating any grid, with division first to avoid integer overflow.
    std::size_t faces=20;
    for(int level=0;;++level) {
        if(faces>limits.maximumArrangementCells || faces/2+2>limits.maximumVertices
            || faces>static_cast<std::size_t>(std::numeric_limits<int>::max()/2))
            throw RankTwoSourceError("RESOURCE_LIMIT","requested fixed normal grid exceeds source budget");
        if(level==settings.icosphereSubdivisions)break;
        if(faces>limits.maximumArrangementCells/4)
            throw RankTwoSourceError("RESOURCE_LIMIT","requested fixed normal grid exceeds source budget");
        faces*=4;
    }
    if(faces+3*(faces/2)>limits.maximumArrangementCells)
        throw RankTwoSourceError("RESOURCE_LIMIT","normal mesh face/edge budget exceeded");
    const auto grid=detail::fixedNormalSphere(settings.icosphereSubdivisions);
    PreparedRankTwoNormalMesh::Rays rays;
    PreparedRankTwoNormalMesh::Faces indices;
    rays.reserve(grid.vertices.size());indices.reserve(grid.faces.size());
    for(const auto& v:grid.vertices)
        rays.push_back({CGAL::Gmpq(v[0]),CGAL::Gmpq(v[1]),CGAL::Gmpq(v[2])});
    for(const auto& f:grid.faces)
        indices.push_back({static_cast<std::size_t>(f[0]),static_cast<std::size_t>(f[1]),static_cast<std::size_t>(f[2])});
    const auto normals=PreparedRankTwoNormalMesh::fromExact(rays,indices,limits);
    return analyzeRankTwoArrangement(source,CGAL::Gmpq(settings.truncationRadius),normals,limits);
}
} // namespace spin
