#pragma once
#include "spin2014/SourceCagePaperFunctional.h"
#include "verification/ExactNefBudget.h"
#include <array>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceCageGeometryOutcome {ChecksPassed,RejectedGeometry,ResourceLimit,Unknown};
enum class SourceCageGeometryStage {Preflight,SourceEvaluation,OuterRegion,PrimitiveRegion,PrimitiveSubset,OuterDomainClearance,PrimitivePair};
enum class SourceCageGeometryStageOutcome {Attempted,Passed,Rejected,ResourceLimit,Unknown,BroadPhaseDisjoint};
struct SourceCageGeometrySettings {
    double minimumOuterDomainClearance=0.;
    std::size_t maximumPoints=100000,maximumTriangles=400000,maximumPrimitives=4096;
    std::size_t maximumPairs=1000000,maximumRecords=1100000;
    std::size_t maximumCheckerBytes=512*1024*1024,maximumNefComplexity=2000000;
    std::size_t maximumConnectedShells=64,maximumExactStringBytes=8192;
};
struct SourceCageGeometryRecord {
    static constexpr std::size_t absent=std::numeric_limits<std::size_t>::max();
    SourceCageGeometryStage stage=SourceCageGeometryStage::Preflight;
    SourceCageGeometryStageOutcome outcome=SourceCageGeometryStageOutcome::Attempted;
    std::size_t first=absent,second=absent,betaOwner=absent,secondBetaOwner=absent;
    std::array<char,65> evidenceSha256{},reasonSha256{};
    std::array<char,160> reason{};
    // Absent until an applicable exact clearance completes. Bounded by settings.
    std::string squaredDistanceExact;
    std::optional<std::array<std::size_t,2>> witnessFaces;
};
struct SourceCageGeometryFacts {
    SourceCageGeometryOutcome outcome=SourceCageGeometryOutcome::Unknown;
    std::string semantics="represented-polyhedral-domain-checks-v1; no continuous-map Jacobian, solver feasibility, convergence, stable spin, seal, or original-paper grant";
    std::string inputSha256,sourceOutputSha256,scientificSha256;
    std::size_t pointCount=0,triangleCount=0,primitiveCount=0,pairCount=0;
    std::size_t plannedRecords=0,plannedCheckerBytes=0,pairsVisited=0,exactPairRelations=0;
    verification::ExactNefBudgetFacts nefBudget;
    std::vector<SourceCageGeometryRecord> records;
};
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
};
const char* sourceCageGeometryOutcomeToken(SourceCageGeometryOutcome);
}
#endif
