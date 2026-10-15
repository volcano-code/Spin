#pragma once
#include "RankTwoCoveredFiniteMesh.h"
#include <cstdint>

namespace spin {
namespace detail {
// Combinatorial prerequisite only: geometry, native origin and curve error
// are supplied separately by the immutable covered-mesh authority.
struct RankTwoComplexBoundary {bool inner;std::vector<std::size_t> points;};
struct RankTwoComplexLink {std::size_t vertex;bool boundary;std::vector<std::size_t> neighbors;};
struct RankTwoComplexComponent {
    std::size_t vertexCount,edgeCount,holes;
    std::int64_t eulerCharacteristic;
    std::vector<std::size_t> triangles,boundaries;
};
struct RankTwoComplexInspection {
    std::vector<RankTwoComplexLink> links;
    std::vector<RankTwoComplexComponent> components;
    std::uint64_t workSteps;
    std::string canonicalRecord;
};
RankTwoComplexInspection inspectRankTwoTriangleComplex(
    const std::vector<std::array<std::size_t,3>>&,
    const std::vector<RankTwoComplexBoundary>&,std::uint64_t maximumSteps);
}
struct RankTwoFiniteCellTopologySettings {std::uint64_t maximumTopologySteps=2000000;};
// Factory-owned actual finite manifold/cell correspondence, not a metric
// isotopy, all-normal-owner assembly, thickness or complete carrier grant.
class RankTwoFiniteCellTopology {
public:
    const RankTwoCoveredFiniteMesh& coveredMesh()const;
    std::size_t regionCount()const;
    std::size_t nativeHoleCount()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoFiniteCellTopology(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoFiniteCellTopology certifyRankTwoFiniteCellTopology(
        const RankTwoArrangement&,const RankTwoCoveredFiniteMesh&,const RankTwoFiniteCellTopologySettings&);
};
RankTwoFiniteCellTopology certifyRankTwoFiniteCellTopology(const RankTwoArrangement&,
    const RankTwoCoveredFiniteMesh&,const RankTwoFiniteCellTopologySettings& settings={});
}
