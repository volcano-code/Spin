#pragma once
#include "RankTwoChordIncidence.h"

namespace spin {
struct RankTwoRegistryRotationSettings {
    std::uint64_t maximumRotationSteps=1000000;
};
// Actual finite chord directions preserve the cyclic order induced by the
// original native DCEL. Unbounded directions remain in the source permutation
// but their geometry is NOT granted. This alone proves neither global curve
// isotopy/nesting nor full strip/fan/error coverage or triangulation permission.
class RankTwoRegistryRotation {
public:
    const std::string& arrangementIdentity()const;
    const std::string& registryIdentity()const;
    const std::string& incidenceIdentity()const;
    const std::string& partitionIdentity()const;
    std::uint64_t pointCount()const;
    std::uint64_t finiteDirectionCount()const;
    std::uint64_t omittedInfiniteDirectionCount()const;
    std::uint64_t omittedOutsidePartitionDirectionCount()const;
    std::size_t maximumFiniteDegree()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoRegistryRotation(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    static RankTwoRegistryRotation certify(const RankTwoArrangement&,const RankTwoArcRegistry&,
        const RankTwoChordIncidence&,const std::string&,const std::vector<std::size_t>&,
        const std::vector<std::size_t>&,const RankTwoRegistryRotationSettings&);
    friend RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const RankTwoChordIncidence&,const RankTwoRegistryRotationSettings&);
    friend RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement&,
        const RankTwoFinitePartition&,const RankTwoChordIncidence&,const RankTwoRegistryRotationSettings&);
};
RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const RankTwoChordIncidence&,
    const RankTwoRegistryRotationSettings& settings=RankTwoRegistryRotationSettings{});
RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement&,
    const RankTwoFinitePartition&,const RankTwoChordIncidence&,
    const RankTwoRegistryRotationSettings& settings=RankTwoRegistryRotationSettings{});
} // namespace spin
