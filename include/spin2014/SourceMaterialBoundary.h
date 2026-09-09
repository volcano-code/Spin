#pragma once
#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "verification/ExactCuboidUnion.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceMaterialBoundarySettings {
    verification::ExactCuboidUnionSettings unionSettings;
    // IMPLEMENTATION_DEFAULT engineering contact-clearance request, not paper
    // thickness. Measures only original boundary to newly carved boundary.
    double minimumOuterWall=0;
};
struct SourceMaterialBoundaryFacts {
    std::string algorithm="source-material-boundary-exact-partition-tetrahedra-v1";
    std::string status="NOT_ATTEMPTED_UNKNOWN",diagnostic,inputDigest,productDigest;
    std::string sourceHash,partitionHash,binaryHash,outerMeshHash,materialMeshHash;
    double density=0,minimumOuterWall=0;
    verification::ExactCuboidUnionFacts unionFacts;
    // Immutable copy at the final validation/failure point, not a live context.
    verification::ExactNefBudgetFacts resourceBudget;
    verification::ExactMaterialRegionFacts outerRegion,designRegion,materialRegion,carvedRegion;
    verification::ExactBoundaryFacts materialBoundary;
    verification::ExactMaterialPartitionFacts partition;
    verification::ExactMaterialRelationFacts carvedSubsetDesign;
    verification::ExactBoundaryClearanceFacts clearance;
    physics::RawMassMoments10 sourceSubtractionMoments,integratedMoments;
    Eigen::Matrix<double,10,1> normalizedMomentResiduals=Eigen::Matrix<double,10,1>::Constant(std::numeric_limits<double>::infinity());
    double normalizedMomentError=std::numeric_limits<double>::infinity();
    double minimumCentralInertiaEigenvalue=std::numeric_limits<double>::infinity();
};
class SourceMaterialBoundaryResult;
/** Immutable native geometry/moment product; arrays are owned copies. This is
 * NOT SEALED_SPIN_ARTIFACT, FEASIBLE_WITNESS, solver/stability/paper success,
 * cavity-web thickness, minimum feature width or manufacturing validation.
 * Original outer arrays remain byte-identical, including pre-existing holes.
 * Cavity union has its own signed outward boundary (including inward voids);
 * its complement orientation is appended to original material geometry.
 */
class SourceMaterialBoundary {
public:
    const SourceMaterialBoundaryFacts& facts() const {return facts_;}
    const Eigen::MatrixXd& outerVertices() const {return outerVertices_;}
    const Eigen::MatrixXi& outerFaces() const {return outerFaces_;}
    const Eigen::MatrixXd& cavityVertices() const {return cavityVertices_;}
    const Eigen::MatrixXi& cavityFaces() const {return cavityFaces_;}
    const Eigen::MatrixXd& materialVertices() const {return materialVertices_;}
    const Eigen::MatrixXi& materialFaces() const {return materialFaces_;}
private:
    SourceMaterialBoundary()=default;
    SourceMaterialBoundaryFacts facts_;
    Eigen::MatrixXd outerVertices_,cavityVertices_,materialVertices_;
    Eigen::MatrixXi outerFaces_,cavityFaces_,materialFaces_;
    friend SourceMaterialBoundaryResult makeSourceMaterialBoundary(const SourceAdaptiveDomain&,const Eigen::VectorXd&,double,const SourceMaterialBoundarySettings&);
};
class SourceMaterialBoundaryResult {
public:
    const SourceMaterialBoundaryFacts& facts() const {return facts_;}
    const std::shared_ptr<const SourceMaterialBoundary>& product() const {return product_;}
private:
    SourceMaterialBoundaryFacts facts_;
    std::shared_ptr<const SourceMaterialBoundary> product_;
    friend SourceMaterialBoundaryResult makeSourceMaterialBoundary(const SourceAdaptiveDomain&,const Eigen::VectorXd&,double,const SourceMaterialBoundarySettings&);
    friend class SourceAdaptiveSpinMaterialResult;
};
void validateSourceMaterialBoundarySettings(const SourceMaterialBoundarySettings&);
std::string hashSourceMaterialBoundarySettings(const SourceMaterialBoundarySettings&);
SourceMaterialBoundaryResult makeSourceMaterialBoundary(const SourceAdaptiveDomain&,
    const Eigen::VectorXd& binary,double density,const SourceMaterialBoundarySettings& = {});

/** Calls the actual native source solver and consumes a private terminal domain
 * moved from that invocation. No public mutable solve-result input is accepted.
 * Numerical and geometry outcomes are separate; UNKNOWN never grants product.
 */
class SourceAdaptiveSpinMaterialResult {
public:
    const SourceAdaptiveSpinResult& solve() const {return solve_;}
    const SourceMaterialBoundaryResult& boundary() const {return boundary_;}
    const std::string& outputDigest() const {return outputDigest_;}
private:
    SourceAdaptiveSpinResult solve_;
    SourceMaterialBoundaryResult boundary_;
    std::string outputDigest_;
    friend SourceAdaptiveSpinMaterialResult solveSourceAdaptiveSpinMaterial(const SourceAdaptiveSpinRequest&,const SourceMaterialBoundarySettings&);
};
SourceAdaptiveSpinMaterialResult solveSourceAdaptiveSpinMaterial(
    const SourceAdaptiveSpinRequest&,const SourceMaterialBoundarySettings& = {});
}
#endif
