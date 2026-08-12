#include "carving/CarvingSolver.h"

#include "physics/BalanceEnergy.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
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


CarvingSolver::RecomputedCarvingState::
RecomputedCarvingState()
    :
    massProperties(),
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
    gravity_(normalizedGravityDirection(gravity)),
    targetCOM_(targetCenterOfMass),
    currCOM_(Eigen::Vector3d::Zero()),
    initialMassProperties_(
        computeSurfaceMassProperties(V, F)
    ),
    currentMassProperties_(initialMassProperties_),
    carvedVoxelCount_(0),
    iterations_(0),
    candidateSearchPerformed_(false),
    lastCandidateVoxelCount_(0),
    diagnosticsEnabled_(false),
    diagnosticsInitialized_(false),
    initialInsideVoxelCount_(0),
    initialCarveableVoxelCount_(0),
    diagnosticsHistory_()
{
    static_cast<void>(planeV);
    static_cast<void>(planeF);

    if (!isClosedConsistentlyOrientedTriangleMesh(F))
    {
        throw std::invalid_argument(
            "outer mesh must be a closed, consistently oriented triangle mesh"
        );
    }

    if (
        !initialMassProperties_.isFinite() ||
        initialMassProperties_.mass <= 0.0
    )
    {
        throw std::invalid_argument(
            "outer mesh must be closed, outward-oriented and have positive mass"
        );
    }

    if (!targetCOM_.allFinite())
    {
        throw std::invalid_argument(
            "target center of mass must be finite"
        );
    }

    const Eigen::Vector3i& dimensions = voxelGrid_.dimensions();
    const Eigen::Vector3d& voxelSize = voxelGrid_.voxelSize();

    if (
        (dimensions.array() < 3).any() ||
        !voxelSize.allFinite() ||
        (voxelSize.array() <= 0.0).any()
    )
    {
        throw std::invalid_argument(
            "carving requires a built voxel grid with positive dimensions"
        );
    }

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (!voxelGrid_(x, y, z).isFilled)
                {
                    throw std::invalid_argument(
                        "carving solver requires an initially filled voxel grid"
                    );
                }
            }
        }
    }

    refreshCenterOfMass();

    ensureDiagnosticsInitialized();
}


double CarvingSolver::getCoMDistance() const
{
    return projectedCenterOfMassDistance(
        currCOM_,
        targetCOM_,
        gravity_
    );
}


double CarvingSolver::getCoMEnergy() const
{
    return projectedCenterOfMassEnergy(
        currCOM_,
        targetCOM_,
        gravity_
    );
}


double CarvingSolver::getGeometricProjectedDistance() const
{
    return getCoMDistance();
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

    if (
        candidateSearchPerformed_ &&
        lastCandidateVoxelCount_ == 0
    )
    {
        return CarvingStopReason::NoCandidatesAvailable;
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


CarvingSolver::RecomputedCarvingState
CarvingSolver::recomputeCarvingState() const
{
    RecomputedCarvingState state;
    state.massProperties = initialMassProperties_;

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel =
                    voxelGrid_(x, y, z);

                if (
                    voxel.distanceToMesh < 0.0 &&
                    voxel.isFilled
                )
                {
                    ++state.filledInsideVoxelCount;
                }

                if (
                    voxel.isInner &&
                    !voxel.isFilled
                )
                {
                    state.massProperties -=
                        computeAxisAlignedBoxMassProperties(
                            voxel.center,
                            voxelSize
                        );
                }
            }
        }
    }

    if (
        state.massProperties.mass <= 0.0 ||
        !state.massProperties.tryCenterOfMass(
            state.centerOfMass,
            0.0
        )
    )
    {
        throw std::runtime_error(
            "recomputed carving state has zero/non-finite remaining mass"
        );
    }

    return state;
}


void CarvingSolver::subtractVoxelMassProperties(
    const Voxel& voxel
)
{
    const MassProperties contribution =
        computeAxisAlignedBoxMassProperties(
            voxel.center,
            voxelGrid_.voxelSize()
        );

    const MassProperties proposed =
        currentMassProperties_ - contribution;

    Eigen::Vector3d proposedCenterOfMass =
        Eigen::Vector3d::Zero();

    if (
        proposed.mass <= 0.0 ||
        !proposed.tryCenterOfMass(
            proposedCenterOfMass,
            0.0
        )
    )
    {
        throw std::runtime_error(
            "carving would produce non-positive/non-finite remaining mass"
        );
    }

    currentMassProperties_ = proposed;
    currCOM_ = proposedCenterOfMass;
}


void CarvingSolver::refreshCenterOfMass()
{
    if (
        currentMassProperties_.mass <= 0.0 ||
        !currentMassProperties_.tryCenterOfMass(
            currCOM_,
            0.0
        )
    )
    {
        throw std::runtime_error(
            "carving produced zero/non-finite remaining mass"
        );
    }
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
    diagnostics.firstMomentBefore =
        currentMassProperties_.firstMoment;

    diagnostics.projectedDistanceBefore =
        getCoMDistance();

    diagnostics.projectedEnergyBefore =
        getCoMEnergy();

    diagnostics.initialMass =
        initialMassProperties_.mass;

    diagnostics.remainingMassBefore =
        currentMassProperties_.mass;

    ++iterations_;
    diagnostics.iteration = iterations_;

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();

    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

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

    candidateSearchPerformed_ = true;
    lastCandidateVoxelCount_ = candidateCount;

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

        subtractVoxelMassProperties(*currentVoxel);

        currentVoxel->isFilled = false;
        ++carvedVoxelCount_;
    }

    const RecomputedCarvingState recomputedState =
        recomputeCarvingState();

    diagnostics.centerOfMassAfter = currCOM_;
    diagnostics.recomputedCenterOfMassAfter =
        recomputedState.centerOfMass;

    diagnostics.firstMomentAfter =
        currentMassProperties_.firstMoment;

    diagnostics.recomputedFirstMomentAfter =
        recomputedState.massProperties.firstMoment;

    diagnostics.projectedDistanceAfter =
        getCoMDistance();

    diagnostics.projectedEnergyAfter =
        getCoMEnergy();

    diagnostics.centerOfMassDiscrepancyAfter =
        (
            currCOM_ -
            recomputedState.centerOfMass
        ).norm();

    diagnostics.firstMomentDiscrepancyAfter =
        (
            currentMassProperties_.firstMoment -
            recomputedState.massProperties.firstMoment
        ).norm();

    diagnostics.remainingMassAfter =
        currentMassProperties_.mass;

    diagnostics.recomputedRemainingMassAfter =
        recomputedState.massProperties.mass;

    diagnostics.massDiscrepancyAfter =
        std::abs(
            currentMassProperties_.mass -
            recomputedState.massProperties.mass
        );

    diagnostics.carvedVoxelCountThisStep =
        carvedThisStep;

    diagnostics.totalCarvedVoxelCount =
        carvedVoxelCount_;

    diagnostics.initialInsideVoxelCount =
        initialInsideVoxelCount_;

    diagnostics.initialCarveableVoxelCount =
        initialCarveableVoxelCount_;

    diagnostics.remainingInsideVoxelCount =
        recomputedState.filledInsideVoxelCount;

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


const MassProperties& CarvingSolver::initialMassProperties() const
{
    return initialMassProperties_;
}


const MassProperties& CarvingSolver::currentMassProperties() const
{
    return currentMassProperties_;
}


double CarvingSolver::objectMass() const
{
    return currentMassProperties_.mass;
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
