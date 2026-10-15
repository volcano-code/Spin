#pragma once
#include "RankTwoMaterialLifts.h"
#include "RankTwoMovingMeshEdge.h"

namespace spin {
struct RankTwoStrictMaterialMeshSettings {
    std::size_t maximumNewEdgePoints=1000000;
    std::uint64_t maximumInventoryChecks=1000000;
    RankTwoMaterialLiftSettings lifts;
    RankTwoMovingMeshEdgeSettings edgePoints;
    RankTwoInteriorMeshSettings insertion;
};
// Private composition of the ACTUAL covered mesh, exact material lifts and
// complete zero-internal-edge repair. A caller cannot choose a successful
// subset of edges or manufacture a thickness grant. Required native moving
// interior nodes must already be in the input mesh. No missing region is
// deleted, and the original shared constraint edges are never resampled.
// This local strict certificate is NOT global curve isotopy, a public common
// carrier, closed 3D material, source motion, or a manufacturing wall check.
class RankTwoStrictMaterialMesh {
public:
    const std::string& parentMeshIdentity()const;
    const RankTwoInteriorMesh& mesh()const;
    const RankTwoMaterialLifts& lifts()const;
    const std::vector<RankTwoMovingMeshEdgePoint>& insertedEdgePoints()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoStrictMaterialMesh(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoStrictMaterialMesh completeRankTwoMaterialThickness(
        const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoInteriorMesh&,
        const RankTwoStrictMaterialMeshSettings&);
};
RankTwoStrictMaterialMesh completeRankTwoMaterialThickness(
    const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoInteriorMesh&,
    const RankTwoStrictMaterialMeshSettings& settings={});
} // namespace spin
