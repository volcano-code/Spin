#pragma once

#include "deformation/BoundedBiharmonicWeights.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "optimization/AlternatingBalanceOptimizer.h"
#if MAKE_IT_STAND_HAS_REACHABILITY
#include "reachability/CarvingReachability.h"
#endif
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
#include "spin/QuadricHollowingSolver.h"
#endif

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>


namespace api
{

#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER

enum class SpinItFasterEvidence
{
    PaperFormulaReproduction,
    ReconstructedAsset
};


const char* spinItFasterEvidenceToken(SpinItFasterEvidence evidence) noexcept;


/** Immutable, hash-bound result artifact for the Spin-It Faster solver. */
class SpinItFasterArtifact
{
public:
    SpinItFasterArtifact(const SpinItFasterArtifact&) noexcept = default;
    SpinItFasterArtifact(SpinItFasterArtifact&&) = delete;
    SpinItFasterArtifact& operator=(const SpinItFasterArtifact&) noexcept =
        default;
    SpinItFasterArtifact& operator=(SpinItFasterArtifact&&) = delete;
    ~SpinItFasterArtifact();

    std::string statusToken() const;
    std::string evidenceToken() const;
    bool succeeded() const;
    bool geometrySealed() const;

    const spin::QuadricHollowingResult& result() const;
    const Eigen::MatrixXd& outerVertices() const;
    const Eigen::MatrixXi& outerFaces() const;
    const Eigen::MatrixXd& designVertices() const;
    const Eigen::MatrixXi& designFaces() const;
    const Eigen::MatrixXd& materialVertices() const;
    const Eigen::MatrixXi& materialFaces() const;
    const Eigen::MatrixXd& cavityVertices() const;
    const Eigen::MatrixXi& cavityFaces() const;
    const std::string& diagnosticsJson() const;

    void saveOuterMesh(const std::string& filename) const;
    void saveDesignMesh(const std::string& filename) const;
    void saveMaterialMesh(const std::string& filename) const;
    void saveCavityMesh(const std::string& filename) const;
    void saveDiagnostics(const std::string& filename) const;

private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;

    explicit SpinItFasterArtifact(std::shared_ptr<const Storage> storage);

    friend SpinItFasterArtifact solveSpinItFaster(
        const spin::SpinDesignRequest&,
        SpinItFasterEvidence,
        const spin::QuadricMomentEvaluatorSettings&,
        const std::vector<spin::QuadricCoefficients10>&
    );
};


SpinItFasterArtifact solveSpinItFaster(
    const spin::SpinDesignRequest& request,
    SpinItFasterEvidence evidence,
    const spin::QuadricMomentEvaluatorSettings& evaluatorSettings = {},
    const std::vector<spin::QuadricCoefficients10>& initialGuesses = {}
);

#endif


/** Cell-centered, canonical x-fast material occupancy import. */
struct MaterialOccupancyData
{
    Eigen::Vector3i dimensions = Eigen::Vector3i::Zero();
    Eigen::Vector3d minimum = Eigen::Vector3d::Zero();
    Eigen::Vector3d maximum = Eigen::Vector3d::Ones();
    std::vector<unsigned char> materialMaskXFast;
    int protectedBoundaryLayers = 1;
};


/** Discretization and BBW settings that are fixed for one model. */
struct BalanceModelConfig
{
    int voxelResolution = 50;
    double minimumWallThickness = 0.0;
    deformation::WeightSolveOptions bbw;

    /** Opt-in released-beta [0,1]^3 domain with N cells per axis. */
    bool fixedUnitCubeVoxelGrid = false;
};


/** Per-run controller settings. */
struct BalanceRunConfig
{
    optimization::AlternatingBalanceOptimizerConfig optimizer;
    double initialMu = 0.2;
};


/**
 * Exact cache payload written by the archived released beta.
 *
 * Cell IDs and BBW rows retain their on-disk z-fast order. The compatibility
 * factory validates, normalizes, and remaps them by integer lattice coordinate.
 */
struct ReleasedBetaArchiveData
{
    int voxelResolution = 0;
    int hullDepth = 2;
    std::vector<int> cellIdsZFast;
    Eigen::MatrixXd betaNodeWeightsZFast;
};


/** Audited opt-in parameter and controller profile for the released beta. */
BalanceRunConfig releasedBetaRunConfig();


/**
 * Immutable, value-owned result of one facade run.
 *
 * Copying this object shares an immutable result snapshot.  No matrix or
 * report below refers to a BalanceModel, optimizer, deformation system, or
 * temporary voxel domain.
 */
class BalanceResult
{
public:
    BalanceResult(const BalanceResult&) noexcept = default;
    BalanceResult(BalanceResult&&) = delete;
    BalanceResult& operator=(const BalanceResult&) noexcept = default;
    BalanceResult& operator=(BalanceResult&&) = delete;
    ~BalanceResult();

    const optimization::AlternatingBalanceOptimizationResult& result() const;
    const manufacturing::ManufacturingGeometry& manufacturingGeometry() const;

    std::string statusToken() const;
    bool succeeded() const;

    const Eigen::MatrixXd& outerVertices() const;
    const Eigen::MatrixXi& outerFaces() const;
    const Eigen::MatrixXd& innerVertices() const;
    const Eigen::MatrixXi& innerFaces() const;

    bool manufacturingReady() const;
    const std::vector<std::string>& manufacturingFailures() const;
    std::string manufacturingReportJson() const;

    /** Stable schema, deliberately excluding paths and timestamps. */
    const std::string& diagnosticsJson() const;

    /**
     * Non-deterministic wall-clock sidecar for R10 performance studies.
     *
     * Timing is deliberately excluded from diagnosticsJson() so repeated R9
     * runs retain byte-identical scientific diagnostics.  The sidecar uses a
     * monotonic clock and reports only phase durations, never host identity,
     * paths, or timestamps.
     */
    const std::string& performanceReportJson() const;

    /** Diagnostic OBJ exports; these are not fabrication certification. */
    void saveOuterMesh(const std::string& filename) const;
    void saveInnerMesh(const std::string& filename) const;

    /** Guarded R8 export; a non-ready result is always rejected. */
    void savePrintableMesh(const std::string& filename) const;
    void savePrintableMesh(
        const std::string& filename,
        manufacturing::ManufacturingMeshFormat format
    ) const;

    void saveDiagnostics(const std::string& filename) const;

private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;

    explicit BalanceResult(std::shared_ptr<const Storage> storage);

    friend class BalanceModel;
};


/**
 * Owning, thread-safe entry point for one input mesh.
 *
 * Mutators are synchronized. optimize() takes a short immutable snapshot and
 * releases the model lock before voxelization, BBW, and optimization, making
 * it safe for a Python binding to release the GIL around the expensive call.
 */
class BalanceModel
{
public:
    explicit BalanceModel(
        const std::string& meshPath,
        const BalanceModelConfig& config = BalanceModelConfig()
    );

    BalanceModel(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& faces,
        const BalanceModelConfig& config = BalanceModelConfig()
    );

    /**
     * Construct from an immutable occupied-cell lattice. Enclosed complement
     * cells are source cavities and are never filled by optimizer refill.
     */
    static std::unique_ptr<BalanceModel> fromMaterialOccupancy(
        const MaterialOccupancyData& occupancy,
        const BalanceModelConfig& config = BalanceModelConfig()
    );

    static std::unique_ptr<BalanceModel> fromReleasedBetaArchive(
        const Eigen::MatrixXd& normalizedVertices,
        const Eigen::MatrixXi& faces,
        const ReleasedBetaArchiveData& archive,
        const BalanceModelConfig& config
    );

    BalanceModel(const BalanceModel&) = delete;
    BalanceModel& operator=(const BalanceModel&) = delete;
    BalanceModel(BalanceModel&&) = delete;
    BalanceModel& operator=(BalanceModel&&) = delete;
    ~BalanceModel();

    Eigen::MatrixXd vertices() const;
    Eigen::MatrixXi faces() const;
    Eigen::Vector3d gravity() const;
    BalanceModelConfig config() const;

    /** Gravity can change only while no standing/suspension target exists. */
    void setGravity(const Eigen::Vector3d& gravityDirection);

    /**
     * Set one actual-contact standing support.
     *
     * Every polygon point must coincide with a lowest source-mesh vertex and
     * no source vertex may lie below its plane. The polygon is canonicalized
     * by SupportRegion::fromPolygon().
     */
    void setSupportPolygon(
        const Eigen::MatrixXd& polygon,
        double gammaRadians
    );

    /** Infer the standing support from source-mesh vertices near its floor. */
    void setSupportFromMesh(
        double contactThreshold,
        double gammaRadians,
        double minimumStabilityInset = 0.0,
        bool supportPlaneAtBandTop = false
    );

    /**
     * Infer support from a provenance-bound reference surface while keeping
     * the model's mass/deformation outer surface unchanged. This is intended
     * only for explicitly labelled repairs of invalid archived surfaces.
     */
    void setSupportFromReferenceMesh(
        const Eigen::MatrixXd& referenceVertices,
        double contactThreshold,
        double gammaRadians,
        double minimumStabilityInset = 0.0,
        bool supportPlaneAtBandTop = false
    );

    /**
     * Add an independent standing pose without replacing prior constraints.
     * Each pose owns its gravity direction; this is the paper's multi-base
     * all-of problem, not multiple contact patches in one pose.
     */
    void addStandingConstraintPolygon(
        const Eigen::MatrixXd& polygon,
        const Eigen::Vector3d& gravityDirection,
        double gammaRadians
    );

    void addStandingConstraintFromMesh(
        const Eigen::Vector3d& gravityDirection,
        double contactThreshold,
        double gammaRadians,
        double minimumStabilityInset = 0.0,
        bool supportPlaneAtBandTop = false
    );

    /** One unlocked translation/uniform-scale handle with fixed rotation. */
    void addHandle(
        const Eigen::Vector3d& restPivot,
        const Eigen::Matrix3d& fixedRotation = Eigen::Matrix3d::Identity()
    );

    /** Add a fixed suspension attachment and angular target. */
    void addSuspensionConstraint(
        const Eigen::Vector3d& attachmentPoint,
        double allowedDeviationRadians
    );

    void clearHandles();
    void clearConstraints();

    BalanceResult optimize(
        const BalanceRunConfig& config = BalanceRunConfig()
    ) const;

#if MAKE_IT_STAND_HAS_REACHABILITY
    /**
     * Classify static balance under a continuous relaxation of carving only.
     * A reachable result is not a discrete or manufacturable cavity proof.
     */
    reachability::CarvingReachabilityResult analyzeCarvingReachability(
        std::size_t standingIndex = 0
    ) const;
#endif

private:
    BalanceModel();

    struct Impl;
    std::unique_ptr<Impl> impl_;
};


} // namespace api
