#pragma once

#include "deformation/BoundedBiharmonicWeights.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "optimization/AlternatingBalanceOptimizer.h"

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>


namespace api
{


/** Discretization and BBW settings that are fixed for one model. */
struct BalanceModelConfig
{
    int voxelResolution = 50;
    double minimumWallThickness = 0.0;
    deformation::WeightSolveOptions bbw;
};


/** Per-run controller settings. */
struct BalanceRunConfig
{
    optimization::AlternatingBalanceOptimizerConfig optimizer;
    double initialMu = 0.2;
};


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
        double gammaRadians
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

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};


} // namespace api
