#pragma once
#include "spin2014/SourceCageMaterialRegions.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageComponentClearanceSettings {
    double minimumComponentSeparation=0.; // finite nonnegative model units
    std::size_t maximumComponents=1000,maximumPoints=100000,maximumTriangles=400000,maximumPairs=100000;
    std::size_t maximumOwnedBytes=256*1024*1024,maximumExactDiagnosticBytes=512*1024*1024,maximumNefComplexity=1000000;
};
enum class SourceCageComponentClearanceOutcome {Validated,NotAdmitted,Rejected,ResourceLimit,Unknown};
struct SourceCageBoundaryComponentEvidence {
    std::size_t componentId=0;
    std::vector<std::size_t> sourceVertexIds,sourceFaceIds;
    int originalSign=0,parent=-1,depth=-1;
    bool checkerReversed=false;
    std::string originalInputSha256,checkerInputSha256;
    SourceCageRegionEvidence region;
};
struct SourceCageComponentPairEvidence {
    std::array<std::size_t,2> componentIds{};
    verification::ExactBoundaryClearanceFacts clearance; // local witness faces
    std::optional<std::array<std::size_t,2>> parentWitnessFaces;
};
struct SourceCageComponentClearanceFacts {
    SourceCageComponentClearanceOutcome outcome=SourceCageComponentClearanceOutcome::Unknown;
    SourceCageComponentClearanceSettings settings;
    std::string scope="BOUNDARY_COMPONENT_SEPARATION";
    std::string semantics="all distinct actual final material boundary components, including root-root air gaps; model units; no same-component feature thickness, continuous BBW Jacobian, continuous author design, seal, feasibility, physics or paper grant";
    std::string fullManufacturingThickness="NOT_ESTABLISHED",state="INCOMPLETE";
    bool admitted=false,identitiesChecked=false,coverageChecked=false;
    std::size_t componentCount=0,plannedPoints=0,plannedTriangles=0,plannedPairCount=0,plannedOwnedBytes=0,plannedExactDiagnosticBytes=0;
    std::vector<int> sourceFaceComponents;
    std::vector<SourceCageBoundaryComponentEvidence> components;
    std::vector<SourceCageComponentPairEvidence> pairs;
    std::optional<bool> applicable,atLeastMinimum,contact;
    std::optional<std::string> minimumSquaredDistanceExact;
    std::optional<std::size_t> minimumPairIndex;
    std::optional<std::array<std::size_t,2>> minimumComponentIds,minimumLocalWitnessFaces,minimumParentWitnessFaces;
    verification::ExactNefBudgetFacts nefBudget; // copy after local exact operands release
    std::string inputSha256,materialInputSha256,configSha256,scientificSha256;
    std::array<char,256> reason{};
};
class SourceCageComponentClearanceAudit;
class SourceCageComponentClearance {
public:
    const SourceCageMaterialRegionsAudit& regionsAudit() const;
    std::shared_ptr<const SourceCageMaterialRegions> regions() const;
private:
    explicit SourceCageComponentClearance(const SourceCageMaterialRegionsAudit& audit):audit_(audit){}
    SourceCageMaterialRegionsAudit audit_;
    friend SourceCageComponentClearanceAudit validateSourceCageComponentClearance(const SourceCageMaterialRegionsAudit&,const SourceCageComponentClearanceSettings&);
};
class SourceCageComponentClearanceAudit {
public:
    const SourceCageComponentClearanceFacts& facts() const;
    std::shared_ptr<const SourceCageComponentClearance> validatedClearance() const;
private:
    struct Data;
    explicit SourceCageComponentClearanceAudit(std::shared_ptr<const Data> data):data_(std::move(data)){}
    std::shared_ptr<const Data> data_;
    friend SourceCageComponentClearanceAudit validateSourceCageComponentClearance(const SourceCageMaterialRegionsAudit&,const SourceCageComponentClearanceSettings&);
};
// Genuine immutable validatedRegions only. Invalid threshold throws; all resource,
// unknown and admission failures withhold publication. Native allocator internals
// remain subject to caller's external memory/time guard; no full thickness grant.
SourceCageComponentClearanceAudit validateSourceCageComponentClearance(const SourceCageMaterialRegionsAudit&,const SourceCageComponentClearanceSettings& = {});
const char* sourceCageComponentClearanceOutcomeToken(SourceCageComponentClearanceOutcome);
}
#endif
