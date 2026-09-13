#pragma once
#include "spin2014/SourceCageBinaryMaterial.h"
#include "verification/ExactMaterialRegion.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageMaterialRegionsSettings {
    double minimumOuterWall=0.; // Finite nonnegative model units; contact always rejects.
    std::size_t maximumPrimitives=1000,maximumPoints=100000,maximumTriangles=400000;
    std::size_t maximumOwnedBytes=256*1024*1024,maximumExactDiagnosticBytes=512*1024*1024;
    std::size_t maximumNefComplexity=1000000;
};
enum class SourceCageMaterialRegionsOutcome {Validated,NotAdmitted,Rejected,ResourceLimit,Unknown};
struct SourceCageRegionEvidence {
    verification::ExactMaterialRegionFacts region;
    // Absent for explicit empty, never a fabricated empty mesh audit.
    std::optional<verification::ExactBoundaryFacts> boundary;
};
struct SourceCageSelectedRegionEvidence {
    std::size_t primitiveId=0,betaOwner=0;
    SourcePartitionKey sourceKey;
    double binaryBeta=0.;
    std::string meshSha256;
    SourceCageRegionEvidence region;
    verification::ExactMaterialRelationFacts subsetOuter;
};
struct SourceCageMaterialRegionsFacts {
    SourceCageMaterialRegionsOutcome outcome=SourceCageMaterialRegionsOutcome::Unknown;
    SourceCageMaterialRegionsSettings settings;
    std::string semantics="finite represented-polyhedral exact regions/partitions/design-membership/outer-clearance only; no seal, continuous BBW Jacobian, full wall/web thickness, physics, feasibility or paper grant";
    std::string designMembershipProof="C = union(selected immutable Full primitives) subset D; D = union(all Full primitives in terminal SourceCageBinding); continuous author design NOT_ESTABLISHED";
    std::string fullManufacturingThickness="NOT_ESTABLISHED";
    bool admitted=false,identitiesChecked=false,representedDesignMembershipEstablished=false;
    std::size_t selectedCount=0,plannedPoints=0,plannedTriangles=0,plannedOwnedBytes=0,plannedExactDiagnosticBytes=0;
    SourceCageRegionEvidence outer,carved,material;
    verification::ExactMaterialPartitionFacts materialPartition,selectedPartition;
    std::vector<SourceCageSelectedRegionEvidence> selected;
    verification::ExactBoundaryClearanceFacts outerClearance;
    verification::ExactNefBudgetFacts nefBudget; // Immutable COPY, not a live context.
    std::string inputSha256,sourceSha256,functionalSha256,partitionSha256,bindingSha256;
    std::string binaryStateSha256,binaryOutputSha256,evaluationSha256,configSha256,scientificSha256;
    std::array<char,256> reason{};
};
class SourceCageMaterialRegionsAudit;
class SourceCageMaterialRegions {
public:
    std::shared_ptr<const SourceCageBinaryMaterial> binaryMaterial() const;
    const SourceCageBinaryMaterialAudit& binaryAudit() const;
private:
    SourceCageMaterialRegions(const SourceCageBinaryMaterialAudit& audit,std::shared_ptr<const SourceCageBinaryMaterial> material):audit_(audit),material_(std::move(material)) {}
    SourceCageBinaryMaterialAudit audit_; // Shared original immutable moments/evidence.
    std::shared_ptr<const SourceCageBinaryMaterial> material_;
    friend SourceCageMaterialRegionsAudit validateSourceCageMaterialRegions(
        const SourceCageBinaryMaterialAudit&,const SourceCageMaterialRegionsSettings&);
};
class SourceCageMaterialRegionsAudit {
public:
    const SourceCageMaterialRegionsFacts& facts() const;
    std::shared_ptr<const SourceCageMaterialRegions> validatedRegions() const;
private:
    struct Data;
    explicit SourceCageMaterialRegionsAudit(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    std::shared_ptr<const Data> data_;
    friend SourceCageMaterialRegionsAudit validateSourceCageMaterialRegions(
        const SourceCageBinaryMaterialAudit&,const SourceCageMaterialRegionsSettings&);
};
// Only immutable upstream validation admits. Invalid threshold throws;
// resource and unknown never publish. Native allocator internals remain under
// the caller's external 24 GiB / 60 minute guard, not an RSS guarantee.
SourceCageMaterialRegionsAudit validateSourceCageMaterialRegions(
    const SourceCageBinaryMaterialAudit&,const SourceCageMaterialRegionsSettings& = {});
const char* sourceCageMaterialRegionsOutcomeToken(SourceCageMaterialRegionsOutcome);
}
#endif
