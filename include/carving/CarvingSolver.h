#pragma once

#include <Eigen/Core>

#include "carving/CarvingDiagnostics.h"

#include <string>
#include <vector>


class VoxelGrid;


/**
 * Implements the mass-redistribution carving procedure.
 *
 * The solver owns the optimization state (mass, current center of mass,
 * iteration count and carved-voxel count), while the referenced VoxelGrid
 * owns all voxel data. It deliberately does not perform visualization or
 * cavity-mesh extraction.
 */
class CarvingSolver
{
public:
    CarvingSolver(
        VoxelGrid& voxelGrid,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& planeV,
        const Eigen::MatrixXi& planeF,
        const Eigen::Vector3d& gravity,
        const Eigen::Vector3d& targetCenterOfMass
    );

    CarvingSolver(const CarvingSolver&) = delete;
    CarvingSolver& operator=(const CarvingSolver&) = delete;

    /** Perform one carving iteration. */
    void carveStep(
        Eigen::MatrixXd& carvePlaneV,
        Eigen::MatrixXi& carvePlaneF,
        Eigen::Vector3d& centerOfMass
    );

    /** Continue carving until an original stopping condition is reached. */
    void carveAll();

    /**
     * Current projected center-of-mass error used by the original solver.
     *
     * Phase 1.6 intentionally preserves this formula so that diagnostics can
     * compare it with a direct geometric projection distance.
     */
    double getCoMDistance() const;

    /** Direct tangential distance between current and target COM on ground. */
    double getGeometricProjectedDistance() const;

    /** Return the currently active original stopping condition. */
    CarvingStopReason stopReason() const;

    /** Return true when any original stopping condition is satisfied. */
    bool isOptimized() const;

    const Eigen::Vector3d& currentCenterOfMass() const;
    const Eigen::Vector3d& targetCenterOfMass() const;

    double objectMass() const;
    int iterationCount() const;
    int carvedVoxelCount() const;

    /** Enable or disable one-line diagnostic output after every step. */
    void setDiagnosticsEnabled(bool enabled);
    bool diagnosticsEnabled() const;

    /** All step records collected since solver construction. */
    const std::vector<CarvingStepDiagnostics>& diagnosticsHistory() const;

    /** Null when no carving step has been executed yet. */
    const CarvingStepDiagnostics* lastDiagnostics() const;

    /** Write the collected diagnostics to a CSV file. */
    bool writeDiagnosticsCsv(
        const std::string& filePath
    ) const;

    int initialInsideVoxelCount() const;
    int initialCarveableVoxelCount() const;

private:
    struct ApproximateVoxelMassProperties
    {
        double mass;
        Eigen::Vector3d centerOfMass;
        int filledInsideVoxelCount;

        ApproximateVoxelMassProperties();
    };

    void ensureDiagnosticsInitialized();

    ApproximateVoxelMassProperties
    computeApproximateVoxelMassProperties() const;

    double projectedDistanceOnGround(
        const Eigen::Vector3d& point
    ) const;

    double carvedGridRatio() const;
    double carvedCarveableRatio() const;

    VoxelGrid& voxelGrid_;

    Eigen::MatrixXd planeV_;
    Eigen::MatrixXi planeF_;

    Eigen::Vector3d gravity_;
    Eigen::Vector3d targetCOM_;
    Eigen::Vector3d currCOM_;

    double mass_;

    int carvedVoxelCount_;
    int iterations_;

    bool diagnosticsEnabled_;
    bool diagnosticsInitialized_;

    int initialInsideVoxelCount_;
    int initialCarveableVoxelCount_;

    std::vector<CarvingStepDiagnostics> diagnosticsHistory_;
};
