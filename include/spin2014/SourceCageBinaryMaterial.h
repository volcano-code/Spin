#pragma once
#include "spin2014/SourceCageAdaptiveSolver.h"
#include "verification/ExactBoundaryAudit.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageBinaryMaterialSettings {
    double constraintTolerance=1e-6; // Finite positive normalized IMPLEMENTATION_DEFAULT.
    SourceCageGeometrySettings geometry;
    std::size_t maximumPoints=100000,maximumTriangles=400000;
    std::size_t maximumOwnedBytes=256*1024*1024,maximumExactDiagnosticBytes=512*1024*1024;
};
enum class SourceCageBinaryMaterialOutcome {Validated,NotAdmitted,BinaryValidationFailed,ResourceLimit,Unknown};
struct SourceCageBinaryMesh {Eigen::MatrixXd vertices;Eigen::MatrixXi triangles;};
struct SourceCageBinaryMaterialFacts {
    using Ten=Eigen::Matrix<double,10,1>;
    SourceCageBinaryMaterialOutcome outcome=SourceCageBinaryMaterialOutcome::Unknown;
    SourceCageBinaryMaterialSettings settings;
    std::string semantics="actual-deformed-binary-pre-seal-v1; no exact final region partition, continuous BBW Jacobian, full cavity-web thickness, feasible witness, spin, physics, seal or paper grant";
    std::string roundingConvention="beta>=0.5 hollowed=1, otherwise 0; IMPLEMENTATION_CONVENTION, not a recovered author tie rule";
    bool endpointAdmitted=false,emptyCavity=false,momentsAvailable=false,constraintsAvailable=false;
    bool momentGatePassed=false,constraintGatePassed=false,positiveMass=false,centralInertiaPsd=false;
    std::size_t selectedPrimitives=0,plannedOwnedBytes=0,plannedExactDiagnosticBytes=0;
    double rootSide=0.,density=0.; // Fixed INPUT-derived normalization, never output-fit.
    double normalizedMomentTolerance=1e-10; // Fixed implementation moment gate.
    double normalizedInertiaTolerance=128*std::numeric_limits<double>::epsilon(); // IMPLEMENTATION_DEFAULT PSD allowance.
    Eigen::VectorXd stateDelta;
    Ten worldIntegratedMoments=Ten::Zero(),spinIntegratedMoments=Ten::Zero();
    Ten binaryFunctionalWorldMoments=Ten::Zero(),normalizedMomentResiduals=Ten::Zero();
    Eigen::VectorXd rawConstraints,normalizedConstraints; // Degree 1 /rho L^4; degree 2 /rho L^5.
    Eigen::Matrix3d centralInertia=Eigen::Matrix3d::Zero(); // Fixed paper coordinates.
    Eigen::Vector3d centralInertiaEigenvalues=Eigen::Vector3d::Zero();
    std::shared_ptr<const SourceCageGeometryAudit> geometry;
    std::shared_ptr<const verification::ExactBoundaryAudit> exactBoundary;
    std::string inputSha256,functionalSha256,partitionSha256,configSha256,outputSha256,scientificSha256;
    std::string relaxedStateSha256,binaryStateSha256,binaryFunctionalOutputSha256;
    std::array<char,256> reason{};
};
class SourceCageBinaryMaterialAudit;
/** Immutable CANDIDATE data, not itself a validation grant. Only the audit's
 * validatedMaterial() accessor publishes validated bounded material. Neither
 * accessor returns a SealedSpinArtifact or proves final exact region partition.
 */
class SourceCageBinaryMaterial {
public:
    const SourceCagePaperFunctional& functional() const;
    const Eigen::VectorXd& relaxedState() const;
    const Eigen::VectorXd& binaryState() const;
    const SourceCageBinaryMesh& outerBoundary() const;
    const SourceCageBinaryMesh& carvedBoundary() const; // Positive cavity union; explicitly empty is legitimate.
    const SourceCageBinaryMesh& materialBoundary() const; // Outer plus reversed cavity.
private:
    struct Data;
    explicit SourceCageBinaryMaterial(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    const Data& data() const;
    std::shared_ptr<const Data> data_;
    friend SourceCageBinaryMaterialAudit validateSourceCageBinaryMaterial(
        const SourceCageAdaptiveSolverAudit&,const SourceCageBinaryMaterialSettings&);
};
/** Shared immutable evidence survives moves while an owner remains. Moved-from
 * access throws. Failed results retain diagnostic candidate data if emission
 * completed, but validatedMaterial() is always null. Copied facts cannot grant.
 */
class SourceCageBinaryMaterialAudit {
public:
    const SourceCageBinaryMaterialFacts& facts() const;
    std::shared_ptr<const SourceCageBinaryMaterial> candidate() const;
    std::shared_ptr<const SourceCageBinaryMaterial> validatedMaterial() const;
private:
    struct Data;
    explicit SourceCageBinaryMaterialAudit(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    std::shared_ptr<const Data> data_;
    friend SourceCageBinaryMaterialAudit validateSourceCageBinaryMaterial(
        const SourceCageAdaptiveSolverAudit&,const SourceCageBinaryMaterialSettings&);
};
// ONLY a privately constructed completed adaptive audit is accepted. No raw
// state, copied facts, pending checkpoint or endpoint-only overload exists.
// No optimization, pose recentering or mesh repair. Invalid tolerances throw;
// resource/unknown failures cannot publish. Known buffers are bounded separately
// from geometry checker caps and native exact-predicate scratch/RSS. Use the
// external 24 GiB / 60 minute guard for exact checking.
SourceCageBinaryMaterialAudit validateSourceCageBinaryMaterial(
    const SourceCageAdaptiveSolverAudit&,const SourceCageBinaryMaterialSettings& = {});
const char* sourceCageBinaryMaterialOutcomeToken(SourceCageBinaryMaterialOutcome);
}
#endif
