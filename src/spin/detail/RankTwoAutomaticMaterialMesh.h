#pragma once
#include "RankTwoInteriorInventory.h"
#include "RankTwoStrictMaterialMesh.h"

namespace spin {
struct RankTwoAutomaticMaterialMeshSettings {
    RankTwoInteriorInventorySettings interiors;
    RankTwoInteriorMeshSettings insertion;
    RankTwoStrictMaterialMeshSettings thickness;
};
// Consume ALL factory-owned native interiors in the SAME covered domain,
// then complete the actual strict material mesh. This local composition is
// not the public shared carrier, global isotopy, closed3D or source motion.
class RankTwoAutomaticMaterialMesh {
public:
    const std::string& baseMeshIdentity()const;
    const RankTwoInteriorInventory& interiors()const;
    const RankTwoInteriorMesh& interiorMesh()const;
    const RankTwoStrictMaterialMesh& strictMaterial()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoAutomaticMaterialMesh(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMesh(
        const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoCoveredFiniteMesh&,
        const RankTwoAutomaticMaterialMeshSettings&);
    friend RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMeshFromInventory(
        const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoCoveredFiniteMesh&,
        const RankTwoInteriorInventory&,const RankTwoInteriorMeshSettings&,const RankTwoStrictMaterialMeshSettings&);
};
RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMesh(
    const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoCoveredFiniteMesh&,
    const RankTwoAutomaticMaterialMeshSettings& settings={});
} // namespace spin
