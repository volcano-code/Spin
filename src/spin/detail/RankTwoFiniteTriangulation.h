#pragma once
#include "RankTwoFinitePartition.h"
#include "RankTwoRegistryRotation.h"

namespace spin {
struct RankTwoFiniteTriangulationSettings {
    std::size_t maximumTriangles=1000000;
    std::uint64_t maximumAdjacencySteps=1000000;
};
struct RankTwoFiniteTriangle {
    std::array<std::size_t,3> points; // shared REGISTRY indices, exact CCW
    std::size_t registryRegionIndex;
};
struct RankTwoRegionTriangulation {
    std::size_t registryRegionIndex,triangleCount,vertexCount,edgeCount,connectedComponents,holes;
    std::int64_t eulerCharacteristic;
};
// Actual finite polygonal triangulation. Mandatory exact incidence and native
// rotation authorities precede No_intersection_tag constraints. Ownership is
// propagated from original LEFT/RIGHT cells, never from parity. This does NOT
// grant curve isotopy/error, positive lifted thickness, a closed source or
// Boolean/motion/manufacturing permission.
class RankTwoFiniteTriangulation {
public:
    const std::string& partitionIdentity()const;
    const std::string& registryIdentity()const;
    const std::vector<RankTwoFiniteTriangle>& triangles()const;
    const std::vector<RankTwoRegionTriangulation>& regionSummaries()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoFiniteTriangulation(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoFiniteTriangulation triangulateRankTwoFinitePartition(const RankTwoArrangement&,
        const RankTwoFinitePartition&,const RankTwoChordIncidence&,const RankTwoRegistryRotation&,
        const RankTwoFiniteTriangulationSettings&);
};
RankTwoFiniteTriangulation triangulateRankTwoFinitePartition(const RankTwoArrangement&,
    const RankTwoFinitePartition&,const RankTwoChordIncidence&,const RankTwoRegistryRotation&,
    const RankTwoFiniteTriangulationSettings& settings={});
RankTwoFiniteTriangulation reloadRankTwoFiniteTriangulation(const RankTwoArrangement&,
    const RankTwoFinitePartition&,const RankTwoChordIncidence&,const RankTwoRegistryRotation&,
    const std::string& saved,const RankTwoFiniteTriangulationSettings& settings={});
} // namespace spin
