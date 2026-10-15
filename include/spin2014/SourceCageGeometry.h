#pragma once
#include "spin2014/SourceCagePaperFunctional.h"
#include "verification/ExactNefBudget.h"
#include <array>
#include <cstdint>
#include <memory>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceCageGeometryOutcome {ChecksPassed,RejectedGeometry,ResourceLimit,Unknown};
enum class SourceCageGeometryStage {Preflight,SourceEvaluation,OuterRegion,PrimitiveRegion,PrimitiveSubset,OuterDomainClearance,PrimitivePair};
enum class SourceCageGeometryStageOutcome {Attempted,Passed,Rejected,ResourceLimit,Unknown,BroadPhaseDisjoint};
struct SourceCageGeometrySettings {
    double minimumOuterDomainClearance=0.;
    std::size_t maximumPoints=100000,maximumTriangles=400000,maximumPrimitives=4096;
    std::size_t maximumPairs=1000000,maximumRecords=1100000;
    std::size_t maximumCheckerBytes=1024*1024*1024,maximumNefComplexity=2000000;
    std::size_t maximumConnectedShells=64,maximumExactStringBytes=8192;
    // Bounded memoization of completed exact primitive-pair overlap proofs.
    // This is an operational setting, never a relaxation of geometric checks.
    std::size_t maximumExactPairCacheEntries=8192;
    // Bounded retention of exact primitive material regions across checker
    // calls. Retained Nef complexity and represented application bytes are
    // separately capped; containment, clearance, and overlap are never cached
    // across a changed outer geometry.
    std::size_t maximumPrimitiveRegionCacheEntries=32;
    std::size_t maximumPrimitiveRegionCacheBytes=16*1024*1024;
    std::size_t maximumPrimitiveRegionCacheComplexity=250000;
};
struct SourceCageGeometryRecord {
    static constexpr std::size_t absent=std::numeric_limits<std::size_t>::max();
    SourceCageGeometryStage stage=SourceCageGeometryStage::Preflight;
    SourceCageGeometryStageOutcome outcome=SourceCageGeometryStageOutcome::Attempted;
    std::size_t first=absent,second=absent,betaOwner=absent,secondBetaOwner=absent;
    std::array<char,65> evidenceSha256{},reasonSha256{};
    std::array<char,160> reason{};
    // Diagnostic wall time for this checker event. It is deliberately excluded
    // from scientificSha256 and must never affect admission or solver status.
    std::uint64_t elapsedNanoseconds=0;
    // Absent until an applicable exact clearance completes. Bounded by settings.
    std::string squaredDistanceExact;
    std::optional<std::array<std::size_t,2>> witnessFaces;
};
struct SourceCageGeometryFacts {
    SourceCageGeometryOutcome outcome=SourceCageGeometryOutcome::Unknown;
    std::string semantics="represented-polyhedral-domain-checks-v3; memoized geometry proofs require identical deformed points and checker identity; pair proofs require identical region digests; no continuous-map Jacobian, solver feasibility, convergence, stable spin, seal, or original-paper grant";
    std::string inputSha256,sourceOutputSha256,scientificSha256;
    std::size_t pointCount=0,triangleCount=0,primitiveCount=0,pairCount=0;
    std::size_t plannedRecords=0,plannedCheckerBytes=0,pairsVisited=0,exactPairRelations=0;
    // Runtime diagnostics only; deliberately excluded from scientificSha256.
    // exactPairCacheHits also counts exact pair propositions covered by a full
    // geometry-proof memo hit (which bypasses the individual pair lookups).
    std::size_t exactPairCacheHits=0,exactPairCacheMisses=0,exactPairBooleanEvaluations=0;
    // Runtime-only counters; excluded from scientificSha256.
    std::size_t primitiveRegionCacheHits=0,primitiveRegionCacheMisses=0;
    // Reusable exact rational boundary-distance BVHs. These counters describe
    // work only and are deliberately excluded from scientificSha256.
    std::size_t outerDistanceIndexBuilds=0,boundaryDistanceIndexBuilds=0;
    std::size_t boundaryDistanceIndexCacheHits=0;
    bool fullGeometryCacheHit=false;
    std::size_t reusedGeometryProofStages=0;
    verification::ExactNefBudgetFacts nefBudget;
    std::vector<SourceCageGeometryRecord> records;
};
struct SourceCageExactPairCache;
class SourceCageGeometryAudit {
public:
    const SourceCageGeometryFacts& facts() const {return facts_;}
private:
    SourceCageGeometryAudit()=default;
    SourceCageGeometryFacts facts_;
    friend class SourceCageGeometryChecker;
};
/** Immutable checker of actual evaluateValues output. Explicit accounting bounds
 * owned mesh/evidence and known exact-audit diagnostics, not native RSS or CGAL
 * allocator internals. Run within the documented 24 GiB/60 minute process guard.
 * Invalid API settings/dimensions throw; attempted resource failure never admits.
 */
class SourceCageGeometryChecker {
public:
    SourceCageGeometryChecker(const SourceCagePaperFunctional&,const SourceCageGeometrySettings&);
    const SourceCagePaperFunctional& source() const {return source_;}
    const SourceCageGeometrySettings& settings() const {return settings_;}
    SourceCageGeometryAudit check(const Eigen::VectorXd& fullState) const;
private:
    const SourceCagePaperFunctional source_;
    const SourceCageGeometrySettings settings_;
    // Copies of this immutable checker share only bounded completed proof facts.
    std::shared_ptr<SourceCageExactPairCache> exactPairCache_;
};
const char* sourceCageGeometryOutcomeToken(SourceCageGeometryOutcome);
}
#endif
