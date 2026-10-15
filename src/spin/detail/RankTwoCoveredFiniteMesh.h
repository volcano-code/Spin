#pragma once
#include "RankTwoFiniteArcCoverage.h"
#include "RankTwoFiniteTriangulation.h"
#include "RankTwoChordSweep.h"

namespace spin {
// Actual common polygonal mesh of the FULL error-controlled finite domain.
// The owned final registry is the ONLY source of CT points/constraints. This
// composes local curve-error coverage with exact polygon incidence/rotation;
// it does not prove globally compatible isotopy or lifted positive thickness.
class RankTwoCoveredFiniteMesh {
public:
    const RankTwoFiniteArcCoverage& coverage()const;
    const RankTwoChordIncidence& incidence()const;
    const RankTwoRegistryRotation& rotation()const;
    const RankTwoFiniteTriangulation& triangulation()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoCoveredFiniteMesh(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoCoveredFiniteMesh triangulateRankTwoCoveredFiniteDomain(
        const RankTwoArrangement&,const RankTwoFiniteArcCoverage&,
        const RankTwoChordSweepSettings&,const RankTwoRegistryRotationSettings&,
        const RankTwoFiniteTriangulationSettings&);
};
RankTwoCoveredFiniteMesh triangulateRankTwoCoveredFiniteDomain(
    const RankTwoArrangement&,const RankTwoFiniteArcCoverage&,
    const RankTwoChordSweepSettings& sweep={},const RankTwoRegistryRotationSettings& rotation={},
    const RankTwoFiniteTriangulationSettings& triangulation={});
} // namespace spin
