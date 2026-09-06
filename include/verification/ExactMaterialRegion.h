#pragma once
#include "verification/ExactBoundaryAudit.h"
#include <memory>
#include <optional>

namespace verification {

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
};

// Malformed input throws invalid_argument; invalid boundaries/nesting retain
// diagnostics but cannot be used as operands. No reorientation of user geometry.
ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
// Explicit empty algebra operand, never a physical outer/artifact witness.
ExactMaterialRegion emptyExactMaterialRegion();
// Invalid operands throw invalid_argument. Boolean failures return all facts
// unknown and completed=false. No thickness, provenance, seal or solve authority.
ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion&, const ExactMaterialRegion&);

} // namespace verification
