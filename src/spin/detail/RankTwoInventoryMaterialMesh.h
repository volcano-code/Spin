#pragma once
#include "RankTwoAutomaticMaterialMesh.h"

namespace spin {
// Consume the already sealed ALL-native inventory of this exact owner. Keep
// its original policy/proofs/identity; never choose or regenerate its points.
// The actual covered mesh and all internal zero-edge repairs remain required.
RankTwoAutomaticMaterialMesh buildRankTwoAutomaticMaterialMeshFromInventory(
    const PreparedRankTwoSource&,const RankTwoArrangement&,const RankTwoCoveredFiniteMesh&,
    const RankTwoInteriorInventory&,const RankTwoInteriorMeshSettings& insertion={},
    const RankTwoStrictMaterialMeshSettings& thickness={});
} // namespace spin
