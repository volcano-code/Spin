#pragma once
#include "RankTwoArcRegistry.h"
#include <cstdint>

namespace spin {
class RankTwoFinitePartition;
struct RankTwoChordSweepSettings;
struct RankTwoChordIncidenceSettings {
    std::uint64_t maximumIncidenceTests=1000000;
};
// Factory-owned geometric incidence of ALL finite shared registry points and
// chords, or ALL of the immutable finite partition identified below. This
// is not a native-arc isotopy, complete strip/fan coverage,
// domain label, infinite clipping or triangulation permission certificate.
class RankTwoChordIncidence {
public:
    const std::string& arrangementIdentity()const;
    const std::string& registryIdentity()const;
    // Empty only for the legacy ALL-registry certificate. Nonempty scope is
    // factory-owned and cannot authorize a different partition or full graph.
    const std::string& partitionIdentity()const;
    std::uint64_t pointPairCount()const;
    std::uint64_t pointEdgeCount()const;
    std::uint64_t edgePairCount()const;
    std::uint64_t predicateGroupCount()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoChordIncidence(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    static RankTwoChordIncidence sweep(const RankTwoArrangement&,const RankTwoArcRegistry&,
        const std::string&,const std::vector<std::size_t>&,const std::vector<std::size_t>&,
        const RankTwoChordSweepSettings&);
    friend RankTwoChordIncidence certifyRankTwoChordIncidence(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const RankTwoChordIncidenceSettings&);
    friend RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement&,
        const RankTwoArcRegistry&,const RankTwoChordSweepSettings&);
    friend RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement&,
        const RankTwoFinitePartition&,const RankTwoChordSweepSettings&);
};
// Rejects distinct feature IDs at identical coordinates, a nonincident
// point lying on a chord, proper crossings and positive-length overlaps.
// Exact algebraic predicates only; no coordinate welding or approximate
// boxes. A resource error never returns a partial positive certificate.
RankTwoChordIncidence certifyRankTwoChordIncidence(const RankTwoArrangement&,
    const RankTwoArcRegistry&,const RankTwoChordIncidenceSettings& settings=RankTwoChordIncidenceSettings{});
} // namespace spin
