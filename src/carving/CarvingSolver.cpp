#include "carving/CarvingSolver.h"

#include "center_of_mass.h"
#include "voxel/VoxelGrid.h"
#include "voxel/VoxelSurface.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include <igl/per_face_normals.h>
#include <igl/signed_distance.h>


namespace
{

constexpr int MAX_ITER = 20;
constexpr int MAX_CARVE = 1000;
constexpr double MAX_EMPTY = 0.80;


bool areFacingSameDirection(
    const Eigen::RowVector3d& first,
    const Eigen::RowVector3d& second
)
{
    return first.dot(second) > 0.0;
}


void flipFaces(Eigen::MatrixXi& faces)
{
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        std::swap(
            faces(row, 0),
            faces(row, 2)
        );
    }
}


/**
 * Preserve the incremental center-of-mass update used by the original
 * implementation. Phase 1.6 measures its discrepancy against a separate
 * voxel-based approximation, but does not yet change the update itself.
 */
void updateCenterOfMass(
    const Voxel& insertedVoxel,
    const Voxel& neighbor,
    const Eigen::MatrixXd& vertices,
    double mass,
    Eigen::Vector3d& globalCenterOfMass
)
{
    Eigen::MatrixXi commonSideFaces(0, 3);

    voxel_surface::appendCommonSide(
        insertedVoxel,
        neighbor,
        vertices,
        commonSideFaces,
        true
    );

    Eigen::Vector3d localCenterOfMass =
        Eigen::Vector3d::Zero();

    center_of_mass(
        vertices,
        commonSideFaces,
        mass,
        localCenterOfMass
    );

    if (voxel_surface::isVoid(neighbor))
    {
        globalCenterOfMass -= localCenterOfMass;
    }
    else
    {
        globalCenterOfMass += localCenterOfMass;
    }
}


struct CompareVoxelDistance
{
    bool operator()(
        const std::pair<Voxel*, double>& first,
        const std::pair<Voxel*, double>& second
    ) const
    {
        return first.second < second.second;
    }
};

} // namespace


CarvingSolver::ApproximateVoxelMassProperties::
ApproximateVoxelMassProperties()
    :
    mass(0.0),
    centerOfMass(Eigen::Vector3d::Zero()),
    filledInsideVoxelCount(0)
{
}


CarvingSolver::CarvingSolver(
    VoxelGrid& voxelGrid,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& planeV,
    const Eigen::MatrixXi& planeF,
    const Eigen::Vector3d& gravity,
    const Eigen::Vector3d& targetCenterOfMass
)
    :
    voxelGrid_(voxelGrid),
    planeV_(planeV),
    planeF_(planeF),
    gravity_(gravity),
    targetCOM_(targetCenterOfMass),
    currCOM_(Eigen::Vector3d::Zero()),
    mass_(0.0),
    carvedVoxelCount_(0),
    iterations_(0),
    diagnosticsEnabled_(false),
    diagnosticsInitialized_(false),
    initialInsideVoxelCount_(0),
    initialCarveableVoxelCount_(0),
    diagnosticsHistory_()
{
    mass_ = center_of_mass(
        V,
        F,
        currCOM_
    );
}


double CarvingSolver::getCoMDistance() const
{
    Eigen::MatrixXd centerOfMassPoints(2, 3);

    centerOfMassPoints <<
        targetCOM_.transpose(),
        currCOM_.transpose();

    Eigen::VectorXd planeDistances;
    Eigen::VectorXd closestFaces;
    Eigen::MatrixXd closestPoints;
    Eigen::MatrixXd closestNormals;

    igl::signed_distance(
        centerOfMassPoints,
        planeV_,
        planeF_,
        igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        planeDistances,
        closestFaces,
        closestPoints,
        closestNormals
    );

    // Preserve the original numerical formula during diagnostic collection.
    const double centerDistanceSquared =
        (
            centerOfMassPoints.row(0) -
            centerOfMassPoints.row(1)
        ).squaredNorm();

    const double centerHeightDifference =
        std::abs(
            planeDistances(1) -
            planeDistances(0)
        );

    const double firstTerm =
        std::pow(centerDistanceSquared, 2.0);

    const double secondTerm =
        std::pow(centerHeightDifference, 2.0);

    return std::sqrt(
        std::max(firstTerm, secondTerm) -
        std::min(firstTerm, secondTerm)
    );
}


double CarvingSolver::projectedDistanceOnGround(
    const Eigen::Vector3d& point
) const
{
    if (planeF_.rows() == 0 || planeV_.rows() == 0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const int firstIndex = planeF_(0, 0);
    const int secondIndex = planeF_(0, 1);
    const int thirdIndex = planeF_(0, 2);

    const Eigen::Vector3d first =
        planeV_.row(firstIndex).transpose();

    const Eigen::Vector3d second =
        planeV_.row(secondIndex).transpose();

    const Eigen::Vector3d third =
        planeV_.row(thirdIndex).transpose();

    Eigen::Vector3d normal =
        (second - first).cross(third - first);

    const double normalLength = normal.norm();

    if (normalLength <= std::numeric_limits<double>::epsilon())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    normal /= normalLength;

    const Eigen::Vector3d difference =
        point - targetCOM_;

    const Eigen::Vector3d tangentialDifference =
        difference - difference.dot(normal) * normal;

    return tangentialDifference.norm();
}


double CarvingSolver::getGeometricProjectedDistance() const
{
    return projectedDistanceOnGround(currCOM_);
}


double CarvingSolver::carvedGridRatio() const
{
    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    const int totalGridVoxelCount =
        dimensions.prod();

    if (totalGridVoxelCount <= 0)
    {
        return 0.0;
    }

    return carvedVoxelCount_ /
        static_cast<double>(totalGridVoxelCount);
}


double CarvingSolver::carvedCarveableRatio() const
{
    if (initialCarveableVoxelCount_ <= 0)
    {
        return 0.0;
    }

    return carvedVoxelCount_ /
        static_cast<double>(initialCarveableVoxelCount_);
}


CarvingStopReason CarvingSolver::stopReason() const
{
    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

    if (getCoMDistance() < voxelSize.maxCoeff())
    {
        return CarvingStopReason::TargetToleranceReached;
    }

    if (carvedGridRatio() > MAX_EMPTY)
    {
        return CarvingStopReason::MaximumEmptyRatioReached;
    }

    if (iterations_ > MAX_ITER)
    {
        return CarvingStopReason::MaximumIterationsReached;
    }

    return CarvingStopReason::None;
}


bool CarvingSolver::isOptimized() const
{
    return stopReason() != CarvingStopReason::None;
}


void CarvingSolver::ensureDiagnosticsInitialized()
{
    if (diagnosticsInitialized_)
    {
        return;
    }

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    initialInsideVoxelCount_ = 0;
    initialCarveableVoxelCount_ = 0;

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel =
                    voxelGrid_(x, y, z);

                if (voxel.distanceToMesh < 0.0)
                {
                    ++initialInsideVoxelCount_;
                }

                if (voxel.isInner)
                {
                    ++initialCarveableVoxelCount_;
                }
            }
        }
    }

    diagnosticsInitialized_ = true;
}


CarvingSolver::ApproximateVoxelMassProperties
CarvingSolver::computeApproximateVoxelMassProperties() const
{
    ApproximateVoxelMassProperties properties;

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

    const double voxelVolume =
        voxelSize(0) * voxelSize(1) * voxelSize(2);

    Eigen::Vector3d weightedCenter =
        Eigen::Vector3d::Zero();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel =
                    voxelGrid_(x, y, z);

                const bool isInsideOriginalMesh =
                    voxel.distanceToMesh < 0.0;

                if (
                    isInsideOriginalMesh &&
                    voxel.isFilled
                )
                {
                    ++properties.filledInsideVoxelCount;
                    weightedCenter += voxel.center;
                }
            }
        }
    }

    properties.mass =
        properties.filledInsideVoxelCount *
        voxelVolume;

    if (properties.filledInsideVoxelCount > 0)
    {
        properties.centerOfMass =
            weightedCenter /
            static_cast<double>(
                properties.filledInsideVoxelCount
            );
    }

    return properties;
}


void CarvingSolver::carveStep(
    Eigen::MatrixXd& carvePlaneV,
    Eigen::MatrixXi& carvePlaneF,
    Eigen::Vector3d& centerOfMass
)
{
    ensureDiagnosticsInitialized();

    CarvingStepDiagnostics diagnostics;

    diagnostics.centerOfMassBefore = currCOM_;
    diagnostics.legacyProjectedErrorBefore = getCoMDistance();
    diagnostics.geometricProjectedErrorBefore =
        getGeometricProjectedDistance();

    ++iterations_;
    diagnostics.iteration = iterations_;

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

    const Eigen::MatrixXd& corners =
        voxelGrid_.corners();


    Eigen::MatrixXd cuttingPlaneTriangleVertices(3, 3);
    Eigen::MatrixXi cuttingPlaneTriangleFaces(1, 3);

    cuttingPlaneTriangleVertices <<
        targetCOM_.transpose(),
        currCOM_.transpose(),
        (targetCOM_ + gravity_).transpose();

    cuttingPlaneTriangleFaces <<
        0, 1, 2;

    Eigen::MatrixXd cuttingPlaneTriangleNormals;

    igl::per_face_normals(
        cuttingPlaneTriangleVertices,
        cuttingPlaneTriangleFaces,
        cuttingPlaneTriangleNormals
    );


    const double gridExtent =
        dimensions.maxCoeff() *
        voxelSize.maxCoeff() *
        10.0;

    const double gravitySquaredNorm =
        gravity_.squaredNorm();

    const double planeNormalSquaredNorm =
        cuttingPlaneTriangleNormals.row(0).squaredNorm();

    assert(gravitySquaredNorm > 0.0);
    assert(planeNormalSquaredNorm > 0.0);

    const double gravityMultiplier =
        gridExtent / gravitySquaredNorm;

    const double otherVectorMultiplier =
        gridExtent / planeNormalSquaredNorm;


    carvePlaneV.resize(5, 3);

    carvePlaneV <<
        targetCOM_.transpose(),
        (
            targetCOM_ +
            gravity_ * gravityMultiplier
        ).transpose(),
        targetCOM_.transpose() +
            cuttingPlaneTriangleNormals.row(0) *
            otherVectorMultiplier,
        (
            targetCOM_ -
            gravity_ * gravityMultiplier
        ).transpose(),
        targetCOM_.transpose() -
            cuttingPlaneTriangleNormals.row(0) *
            otherVectorMultiplier;

    carvePlaneF.resize(4, 3);

    carvePlaneF <<
        0, 1, 2,
        2, 3, 0,
        0, 3, 4,
        4, 1, 0;


    Eigen::MatrixXd carvingPlaneNormals;

    igl::per_face_normals(
        carvePlaneV,
        carvePlaneF,
        carvingPlaneNormals
    );

    assert(
        areFacingSameDirection(
            carvingPlaneNormals.row(0),
            carvingPlaneNormals.row(1)
        )
    );

    assert(
        areFacingSameDirection(
            carvingPlaneNormals.row(1),
            carvingPlaneNormals.row(2)
        )
    );

    assert(
        areFacingSameDirection(
            carvingPlaneNormals.row(2),
            carvingPlaneNormals.row(3)
        )
    );


    Eigen::MatrixXd currentCenterPoint(1, 3);

    currentCenterPoint.row(0) =
        currCOM_.transpose();

    Eigen::VectorXd currentCenterDistance;
    Eigen::VectorXd currentCenterClosestFaces;
    Eigen::MatrixXd currentCenterClosestPoints;
    Eigen::MatrixXd currentCenterClosestNormals;

    igl::signed_distance(
        currentCenterPoint,
        carvePlaneV,
        carvePlaneF,
        igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        currentCenterDistance,
        currentCenterClosestFaces,
        currentCenterClosestPoints,
        currentCenterClosestNormals
    );

    if (currentCenterDistance(0) < 0.0)
    {
        flipFaces(carvePlaneF);
    }


    std::vector<std::pair<Voxel*, double> > candidateVoxels;

    for (int x = 1; x < dimensions(0) - 1; ++x)
    {
        for (int y = 1; y < dimensions(1) - 1; ++y)
        {
            for (int z = 1; z < dimensions(2) - 1; ++z)
            {
                Voxel* currentVoxel =
                    &voxelGrid_(x, y, z);

                if (
                    !currentVoxel->isInner ||
                    !currentVoxel->isFilled
                )
                {
                    continue;
                }

                Eigen::MatrixXd voxelCenterPoint(1, 3);

                voxelCenterPoint.row(0) =
                    currentVoxel->center.transpose();

                Eigen::VectorXd voxelDistance;
                Eigen::VectorXd voxelClosestFaces;
                Eigen::MatrixXd voxelClosestPoints;
                Eigen::MatrixXd voxelClosestNormals;

                igl::signed_distance(
                    voxelCenterPoint,
                    carvePlaneV,
                    carvePlaneF,
                    igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
                    std::numeric_limits<double>::min(),
                    std::numeric_limits<double>::max(),
                    voxelDistance,
                    voxelClosestFaces,
                    voxelClosestPoints,
                    voxelClosestNormals
                );

                if (voxelDistance(0) > 0.0)
                {
                    candidateVoxels.push_back(
                        std::make_pair(
                            currentVoxel,
                            voxelDistance(0)
                        )
                    );
                }
            }
        }
    }


    std::sort(
        candidateVoxels.begin(),
        candidateVoxels.end(),
        CompareVoxelDistance()
    );

    const int candidateCount =
        static_cast<int>(candidateVoxels.size());

    diagnostics.candidateVoxelCount =
        candidateCount;

    const int firstCandidateIndex =
        std::max(
            0,
            candidateCount - MAX_CARVE
        );

    const int carvedThisStep =
        candidateCount - firstCandidateIndex;

    for (
        int candidateIndex = candidateCount - 1;
        candidateIndex >= firstCandidateIndex;
        --candidateIndex
    )
    {
        Voxel* currentVoxel =
            candidateVoxels[candidateIndex].first;

        const int x = currentVoxel->xIdx;
        const int y = currentVoxel->yIdx;
        const int z = currentVoxel->zIdx;

        currentVoxel->isFilled = false;
        ++carvedVoxelCount_;

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x + 1, y, z),
            corners,
            mass_,
            currCOM_
        );

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x - 1, y, z),
            corners,
            mass_,
            currCOM_
        );

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x, y + 1, z),
            corners,
            mass_,
            currCOM_
        );

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x, y - 1, z),
            corners,
            mass_,
            currCOM_
        );

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x, y, z + 1),
            corners,
            mass_,
            currCOM_
        );

        updateCenterOfMass(
            *currentVoxel,
            voxelGrid_(x, y, z - 1),
            corners,
            mass_,
            currCOM_
        );
    }

    const ApproximateVoxelMassProperties voxelProperties =
        computeApproximateVoxelMassProperties();

    diagnostics.centerOfMassAfter = currCOM_;
    diagnostics.approximateVoxelCenterOfMassAfter =
        voxelProperties.centerOfMass;

    diagnostics.legacyProjectedErrorAfter =
        getCoMDistance();

    diagnostics.geometricProjectedErrorAfter =
        getGeometricProjectedDistance();

    diagnostics.centerOfMassDiscrepancyAfter =
        (
            currCOM_ -
            voxelProperties.centerOfMass
        ).norm();

    diagnostics.exactInitialMass = mass_;
    diagnostics.approximateRemainingVoxelMass =
        voxelProperties.mass;

    diagnostics.carvedVoxelCountThisStep =
        carvedThisStep;

    diagnostics.totalCarvedVoxelCount =
        carvedVoxelCount_;

    diagnostics.initialInsideVoxelCount =
        initialInsideVoxelCount_;

    diagnostics.initialCarveableVoxelCount =
        initialCarveableVoxelCount_;

    diagnostics.remainingInsideVoxelCount =
        voxelProperties.filledInsideVoxelCount;

    diagnostics.carvedGridRatio =
        carvedGridRatio();

    diagnostics.carvedCarveableRatio =
        carvedCarveableRatio();

    diagnostics.stopReason =
        stopReason();

    diagnosticsHistory_.push_back(diagnostics);

    if (diagnosticsEnabled_)
    {
        printCarvingStepDiagnostics(
            std::cout,
            diagnostics
        );
    }

    centerOfMass = currCOM_;
}


void CarvingSolver::carveAll()
{
    while (!isOptimized())
    {
        Eigen::MatrixXd temporaryPlaneVertices;
        Eigen::MatrixXi temporaryPlaneFaces;
        Eigen::Vector3d temporaryCenterOfMass;

        carveStep(
            temporaryPlaneVertices,
            temporaryPlaneFaces,
            temporaryCenterOfMass
        );
    }
}


const Eigen::Vector3d& CarvingSolver::currentCenterOfMass() const
{
    return currCOM_;
}


const Eigen::Vector3d& CarvingSolver::targetCenterOfMass() const
{
    return targetCOM_;
}


double CarvingSolver::objectMass() const
{
    return mass_;
}


int CarvingSolver::iterationCount() const
{
    return iterations_;
}


int CarvingSolver::carvedVoxelCount() const
{
    return carvedVoxelCount_;
}


void CarvingSolver::setDiagnosticsEnabled(
    bool enabled
)
{
    diagnosticsEnabled_ = enabled;
}


bool CarvingSolver::diagnosticsEnabled() const
{
    return diagnosticsEnabled_;
}


const std::vector<CarvingStepDiagnostics>&
CarvingSolver::diagnosticsHistory() const
{
    return diagnosticsHistory_;
}


const CarvingStepDiagnostics*
CarvingSolver::lastDiagnostics() const
{
    if (diagnosticsHistory_.empty())
    {
        return nullptr;
    }

    return &diagnosticsHistory_.back();
}


bool CarvingSolver::writeDiagnosticsCsv(
    const std::string& filePath
) const
{
    std::ofstream output(filePath.c_str());

    if (!output)
    {
        return false;
    }

    writeCarvingDiagnosticsCsvHeader(output);

    for (
        std::vector<CarvingStepDiagnostics>::const_iterator iterator =
            diagnosticsHistory_.begin();
        iterator != diagnosticsHistory_.end();
        ++iterator
    )
    {
        writeCarvingDiagnosticsCsvRow(
            output,
            *iterator
        );
    }

    return static_cast<bool>(output);
}


int CarvingSolver::initialInsideVoxelCount() const
{
    return initialInsideVoxelCount_;
}


int CarvingSolver::initialCarveableVoxelCount() const
{
    return initialCarveableVoxelCount_;
}
