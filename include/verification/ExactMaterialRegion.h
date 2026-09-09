#pragma once
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactNefBudget.h"
#include <array>
#include <memory>
#include <optional>

namespace verification {
struct ExactCuboid;
struct ExactCuboidUnionSettings;
class ExactCuboidUnion;

struct ExactMaterialRegionFacts {
    std::string algorithm = "signed-shell-nef-regular-closed-v1";
    std::string backend = "CGAL-f7c3c8212b56c0d6dae63787efc99093f4383415/Nef_polyhedron_3<Simple_cartesian<Gmpq>>/GMP-6.3.0/MPFR-4.2.1";
    std::string semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string inputDigest, regionDigest, error;
    bool completed = false, valid = false, explicitEmpty = false;
    std::vector<int> parents, depths, originalSigns, checkerReversedComponents;
    // Absent until a valid region has been established.
    std::optional<int> materialComponents, boundedVoidComponents;
};

struct ExactMaterialRelationFacts {
    std::string algorithm = "nef-exact-set-relations-v1";
    std::string semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string aDigest, bDigest, relationDigest, error;
    bool completed = false, aExplicitEmpty = false, bExplicitEmpty = false;
    // Unknown after any incomplete/failed Boolean; never default clearance.
    std::optional<bool> intersectionNonempty, positiveVolumeOverlap;
    std::optional<bool> aSubsetB, bSubsetA, equal;
};

struct ExactMaterialPartitionFacts {
    std::string algorithm = "nef-inclusive-union-interior-disjoint-partition-v1";
    std::string backend, semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string emptyProtocol = "explicit-empty-allowed-nonempty-ordered-part-list-v1";
    std::string wholeDigest, partitionDigest, error;
    std::vector<std::string> partDigests;
    bool completed = false, wholeExplicitEmpty = false;
    std::vector<bool> partExplicitEmpty;
    // All computed facts absent until every exact operation succeeds.
    std::optional<bool> unionEqualsWhole, noPositiveVolumeOverlap;
    std::optional<bool> missingWholeContent, excessUnionContent, completePartition;
    std::optional<std::vector<bool>> partSubsetWhole;
    std::optional<std::vector<std::array<std::size_t,2>>> positiveVolumeOverlapPairs;
};

struct ExactBoundaryClearanceFacts {
    std::string algorithm = "exact-triangle-features-rational-dual-bvh-v1";
    std::string backend;
    std::string semantics = "indexed-triangle-boundaries-exact-squared-distance-ge-v1";
    std::string aDigest, bDigest, aInputDigest, bInputDigest, minimumClearanceExact;
    std::string clearanceDigest, error, state = "INCOMPLETE";
    bool completed = false, aExplicitEmpty = false, bExplicitEmpty = false;
    std::optional<bool> applicable, atLeastMinimum, contact;
    std::optional<std::string> squaredDistanceExact;
    std::optional<std::array<std::size_t,2>> witnessFaces;
};

class ExactMaterialRegion {
public:
    const ExactMaterialRegionFacts& facts() const { return facts_; }
    const ExactBoundaryFacts& boundaryFacts() const { return boundary_; }
private:
    ExactMaterialRegion() = default;
    struct Data;
    std::shared_ptr<const Data> data_;
    ExactMaterialRegionFacts facts_;
    ExactBoundaryFacts boundary_;
    friend ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
    friend ExactMaterialRegion emptyExactMaterialRegion();
    friend ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion&, const ExactMaterialRegion&);
    friend ExactMaterialPartitionFacts checkExactMaterialPartition(const ExactMaterialRegion&, const std::vector<ExactMaterialRegion>&);
    friend ExactBoundaryClearanceFacts checkExactBoundaryClearance(const ExactMaterialRegion&, const ExactMaterialRegion&, double);
    friend ExactCuboidUnion makeExactCuboidUnion(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings&);
    friend ExactCuboidUnion makeExactCuboidUnionWithBudget(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings&,const ExactNefBudgetPtr&);
    friend ExactMaterialRegion makeExactMaterialRegionWithBudget(const Eigen::MatrixXd&,const Eigen::MatrixXi&,const ExactNefBudgetPtr&);
    friend ExactMaterialRegion emptyExactMaterialRegionWithBudget(const ExactNefBudgetPtr&);
    friend ExactMaterialRelationFacts relateExactMaterialRegionsWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
    friend ExactMaterialPartitionFacts checkExactMaterialPartitionWithBudget(const ExactMaterialRegion&,const std::vector<ExactMaterialRegion>&,const ExactNefBudgetPtr&);
};

// Malformed input throws invalid_argument; invalid boundaries/nesting retain
// diagnostics but cannot be used as operands. No reorientation of user geometry.
ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
// Explicit empty algebra operand, never a physical outer/artifact witness.
ExactMaterialRegion emptyExactMaterialRegion();
// Invalid operands throw invalid_argument. Boolean failures return all facts
// unknown and completed=false. No thickness, provenance, seal or solve authority.
ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion&, const ExactMaterialRegion&);

// Inclusive union equality AND pairwise disjoint interiors. An empty part list
// is missing evidence and throws; explicit-empty regions are algebra operands.
ExactMaterialPartitionFacts checkExactMaterialPartition(
    const ExactMaterialRegion& whole, const std::vector<ExactMaterialRegion>& parts);
// Exact represented boundary distance, never containment or shell validity.
// Threshold is a finite nonnegative binary64, converted to rational exactly.
// Empty operands yield completed NOT_APPLICABLE with absent measured facts.
// Invalid operands/thresholds throw; operation failure publishes no clear facts.
ExactBoundaryClearanceFacts checkExactBoundaryClearance(
    const ExactMaterialRegion& a, const ExactMaterialRegion& b, double minimumClearance);

// Explicit opt-in variants. Absent-context legacy calls keep their identities
// and geometric semantics. Resource failures never publish successful facts.
ExactMaterialRegion makeExactMaterialRegionWithBudget(const Eigen::MatrixXd&,const Eigen::MatrixXi&,const ExactNefBudgetPtr&);
ExactMaterialRegion emptyExactMaterialRegionWithBudget(const ExactNefBudgetPtr&);
ExactMaterialRelationFacts relateExactMaterialRegionsWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
ExactMaterialPartitionFacts checkExactMaterialPartitionWithBudget(const ExactMaterialRegion&,const std::vector<ExactMaterialRegion>&,const ExactNefBudgetPtr&);

} // namespace verification
