#pragma once
#include "RankTwoInteriorMesh.h"
#include <cstdint>

namespace spin {
struct RankTwoMaterialLiftSettings { std::uint64_t maximumPredicates=1000000; };
struct RankTwoLiftedTriangle {
    std::size_t triangleIndex,registryRegionIndex;
    bool active=false;
    std::array<bool,3> materialBoundary{};
    bool strictPassed=false;
    int sumThicknessSign=0;
    std::array<bool,3> zeroInternalEdges{};
};
// Actual common mesh and original source authority are retained. Certified
// native zero-feature identities define exact coincident lifts; no tolerances,
// snapping or independently rounded owners. This is NOT a closed 3D source,
// motion proof, curve isotopy grant, or a completed general carrier.
class RankTwoMaterialLifts {
public:
    const std::string& sourceIdentity()const;
    const std::string& meshIdentity()const;
    const std::vector<RankTwoLiftedTriangle>& triangles()const;
    const CORE::Expr& lower(std::size_t regionIndex,std::size_t pointIndex)const;
    const CORE::Expr& upper(std::size_t pointIndex)const;
    bool allActiveStrictPassed()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoMaterialLifts(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoMaterialLifts buildRankTwoMaterialLifts(const PreparedRankTwoSource&,
        const RankTwoArrangement&,const RankTwoInteriorMesh&,const RankTwoMaterialLiftSettings&);
};
RankTwoMaterialLifts buildRankTwoMaterialLifts(const PreparedRankTwoSource&,
    const RankTwoArrangement&,const RankTwoInteriorMesh&,const RankTwoMaterialLiftSettings& settings={});
} // namespace spin
