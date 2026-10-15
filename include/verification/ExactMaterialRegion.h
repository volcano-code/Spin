#pragma once
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactNefBudget.h"
#include <array>
#include <memory>
#include <optional>

namespace verification {
namespace detail { class ExactBoundaryDistanceIndex; }
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

// Trial material semantics from a complete exact embedded boundary followed
// by exact component point location. This owns no Nef and grants no set-op,
// partition, clearance, final artifact or manufacturing authority.
struct ExactMaterialBoundaryFacts {
    ExactBoundaryFacts boundary;
    std::string algorithm = "exact-embedded-shell-nesting-ray-parity-v1", error;
    bool completed = false, classificationValid = false;
    std::vector<int> parents, depths, originalSigns, witnessVertices;
    std::vector<std::vector<bool>> inside;
    std::optional<int> materialComponents, boundedVoidComponents;
    std::string fixedProofIdentity;
    std::uint64_t reusedPairTests = 0;
    double boundarySeconds = 0, nestingSeconds = 0;
};

ExactMaterialBoundaryFacts auditExactRationalMaterialBoundary(
    const std::vector<std::array<std::string,3>>&, const Eigen::MatrixXi&,
    const VerifiedFixedBoundary* fixedProof = nullptr);

struct ExactMaterialRelationFacts {
    std::string algorithm = "nef-exact-set-relations-and-symmetric-difference-volume-v2";
    std::string semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string aDigest, bDigest, relationDigest, error;
    bool completed = false, aExplicitEmpty = false, bExplicitEmpty = false;
    // Unknown after any incomplete/failed Boolean; never default clearance.
    std::optional<bool> intersectionNonempty, positiveVolumeOverlap;
    std::optional<bool> aSubsetB, bSubsetA, equal;
    // Exact rational volume of the regularized symmetric difference. This is
    // a measure, not a Hausdorff bound: zero iff the regularized solids agree
    // up to lower-dimensional boundary differences.
    std::optional<std::string> symmetricDifferenceVolumeExact;
};

enum class ExactMaterialBooleanOperation { Intersection, Difference };

// Narrow exact inclusive-set containment query. It deliberately avoids the
// additional symmetric-difference volume integration performed by the full
// relation report; failures leave the subset fact unknown.
struct ExactMaterialSubsetFacts {
    std::string algorithm = "nef-exact-inclusive-material-subset-v1";
    std::string semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string aDigest, bDigest, error;
    bool completed = false;
    std::optional<bool> aSubsetB;
};

// Narrow exact query for callers that only need to reject positive-volume
// overlap. It intentionally does not compute containment or symmetric volume.
struct ExactMaterialOverlapFacts {
    std::string algorithm = "nef-exact-positive-volume-overlap-v1";
    std::string semantics = "finite-regular-closed-polyhedral-set-v1";
    std::string aDigest, bDigest, operationDigest, error;
    bool completed = false;
    std::optional<bool> positiveVolumeOverlap;
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
    // Bind a retained, immutable proof to the same indexed rational boundary.
    // No Nef reconstruction, rounding, reindexing or hash-only comparison.
    // Malformed, empty, incomplete or unequal inputs cannot borrow authority.
    bool matchesRationalBoundary(const std::vector<std::array<std::string,3>>&,
        const Eigen::MatrixXi&) const;
    // Owned copies of the authoritative indexed boundary, never rounded.
    std::vector<std::array<std::string,3>> rationalBoundaryVertices() const;
    Eigen::MatrixXi boundaryFaces() const;
    bool boundaryUsesRationals() const;
    Eigen::MatrixXd binary64BoundaryVertices() const;
    std::vector<ExactMaterialRegion> constructionOperands() const;
private:
    ExactMaterialRegion() = default;
    static ExactMaterialRegion assemble(const Eigen::MatrixXd&,
        const std::vector<std::array<std::string,3>>*, const Eigen::MatrixXi&,
        const ExactNefBudgetPtr&);
    struct Data;
    std::shared_ptr<const Data> data_;
    ExactMaterialRegionFacts facts_;
    ExactBoundaryFacts boundary_;
    friend ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
    friend ExactMaterialRegion emptyExactMaterialRegion();
    friend ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion&, const ExactMaterialRegion&);
    friend ExactMaterialRelationFacts relateExactBooleanOperationToMaterial(
        const ExactMaterialRegion&, const ExactMaterialRegion&,
        const ExactMaterialRegion&, ExactMaterialBooleanOperation);
    friend ExactMaterialSubsetFacts checkExactMaterialSubset(const ExactMaterialRegion&, const ExactMaterialRegion&);
    friend ExactMaterialOverlapFacts checkExactPositiveVolumeOverlap(const ExactMaterialRegion&,const ExactMaterialRegion&);
    friend ExactMaterialPartitionFacts checkExactMaterialPartition(const ExactMaterialRegion&, const std::vector<ExactMaterialRegion>&);
    friend ExactBoundaryClearanceFacts checkExactBoundaryClearance(const ExactMaterialRegion&, const ExactMaterialRegion&, double);
    friend ExactBoundaryClearanceFacts checkExactBoundaryClearance(const ExactMaterialRegion&,
        const ExactMaterialRegion&,double,const detail::ExactBoundaryDistanceIndex&,
        const detail::ExactBoundaryDistanceIndex&);
    friend ExactCuboidUnion makeExactCuboidUnion(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings&);
    friend ExactCuboidUnion makeExactCuboidUnionWithBudget(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings&,const ExactNefBudgetPtr&);
    friend ExactMaterialRegion makeExactMaterialRegionWithBudget(const Eigen::MatrixXd&,const Eigen::MatrixXi&,const ExactNefBudgetPtr&);
    friend ExactMaterialRegion makeExactMaterialRegionFromRationals(
        const std::vector<std::array<std::string,3>>&,const Eigen::MatrixXi&);
    friend ExactMaterialRegion makeExactInteriorDifferenceRegion(
        const ExactMaterialRegion&, const ExactMaterialRegion&);
    friend ExactMaterialRegion emptyExactMaterialRegionWithBudget(const ExactNefBudgetPtr&);
    friend ExactMaterialRelationFacts relateExactMaterialRegionsWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
    friend ExactMaterialSubsetFacts checkExactMaterialSubsetWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
    friend ExactMaterialOverlapFacts checkExactPositiveVolumeOverlapWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
    friend ExactMaterialPartitionFacts checkExactMaterialPartitionWithBudget(const ExactMaterialRegion&,const std::vector<ExactMaterialRegion>&,const ExactNefBudgetPtr&);
};

// Malformed input throws invalid_argument; invalid boundaries/nesting retain
// diagnostics but cannot be used as operands. No reorientation of user geometry.
ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
// Exact rational coordinates from construction sidecars. Every component is
// decimal numerator/positive-denominator; no conversion to binary64 occurs.
// Provides only set operations, not binary64-only boundary-clearance queries.
ExactMaterialRegion makeExactMaterialRegionFromRationals(
    const std::vector<std::array<std::string,3>>&,const Eigen::MatrixXi&);
// Narrow proof-preserving construction for hollow-only outputs. Outer must
// have one positive shell; carved must have disjoint positive root shells.
// Exact containment plus an embedded concatenated boundary proves strict
// interior removal. Reuses the actual owned Nef operands, not caller facts.
// Unsupported/contact/outside operands throw; Boolean failures remain invalid.
ExactMaterialRegion makeExactInteriorDifferenceRegion(
    const ExactMaterialRegion& outer, const ExactMaterialRegion& carved);
// Explicit empty algebra operand, never a physical outer/artifact witness.
ExactMaterialRegion emptyExactMaterialRegion();
// Invalid operands throw invalid_argument. Boolean failures return all facts
// unknown and completed=false. No thickness, provenance, seal or solve authority.
ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion&, const ExactMaterialRegion&);
// Reference operation computed exactly from the two original binary64 mesh
// operands (or rational construction operands), then compared with an
// independently supplied material. This checks represented polyhedra only;
// it neither certifies continuous quadrics nor changes a solver's output.
ExactMaterialRelationFacts relateExactBooleanOperationToMaterial(
    const ExactMaterialRegion& design, const ExactMaterialRegion& cutter,
    const ExactMaterialRegion& material, ExactMaterialBooleanOperation operation);
// Exact inclusive containment; unlike a positive-volume overlap check, this
// also rejects lower-dimensional content outside the container.
ExactMaterialSubsetFacts checkExactMaterialSubset(const ExactMaterialRegion&, const ExactMaterialRegion&);
// Computes only whether intersection has positive regularized volume. Boundary
// contact is not a positive-volume overlap; failures keep the answer unknown.
ExactMaterialOverlapFacts checkExactPositiveVolumeOverlap(const ExactMaterialRegion&,const ExactMaterialRegion&);

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
/** Exact same clearance proof using caller-retained mesh-bound rational BVHs.
 * A mismatched index is rejected before any geometric proposition is run. */
ExactBoundaryClearanceFacts checkExactBoundaryClearance(const ExactMaterialRegion& a,
    const ExactMaterialRegion& b,double minimumClearance,
    const detail::ExactBoundaryDistanceIndex& aIndex,
    const detail::ExactBoundaryDistanceIndex& bIndex);

// Explicit opt-in variants. Absent-context legacy calls keep their identities
// and geometric semantics. Resource failures never publish successful facts.
ExactMaterialRegion makeExactMaterialRegionWithBudget(const Eigen::MatrixXd&,const Eigen::MatrixXi&,const ExactNefBudgetPtr&);
ExactMaterialRegion emptyExactMaterialRegionWithBudget(const ExactNefBudgetPtr&);
ExactMaterialRelationFacts relateExactMaterialRegionsWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
ExactMaterialSubsetFacts checkExactMaterialSubsetWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
ExactMaterialOverlapFacts checkExactPositiveVolumeOverlapWithBudget(const ExactMaterialRegion&,const ExactMaterialRegion&,const ExactNefBudgetPtr&);
ExactMaterialPartitionFacts checkExactMaterialPartitionWithBudget(const ExactMaterialRegion&,const std::vector<ExactMaterialRegion>&,const ExactNefBudgetPtr&);

} // namespace verification
