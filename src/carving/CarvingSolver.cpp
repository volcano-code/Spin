#include "carving/CarvingSolver.h"

#include "physics/BalanceEnergy.h"
#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"
#include "voxel/VoxelTopology.h"

#include <Eigen/Geometry>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>


namespace
{

constexpr int MAX_ITER = 100;


class StagedVoxelRemoval
{
public:
    StagedVoxelRemoval(
        VoxelGrid& voxelGrid,
        const PlaneCarvingCandidates& candidates,
        std::size_t prefixLength
    )
        : voxelGrid_(voxelGrid), keys_(), committed_(false)
    {
        keys_.reserve(prefixLength);

        for (std::size_t index = 0; index < prefixLength; ++index)
        {
            const VoxelKey key = candidates[index].key;
            Voxel& voxel = voxelGrid_(key.x, key.y, key.z);

            if (!voxel.isFilled)
            {
                throw std::logic_error(
                    "carving prefix contains an already empty voxel"
                );
            }

            keys_.push_back(key);
            voxel.isFilled = false;
        }
    }

    StagedVoxelRemoval(const StagedVoxelRemoval&) = delete;
    StagedVoxelRemoval& operator=(const StagedVoxelRemoval&) = delete;

    ~StagedVoxelRemoval()
    {
        if (!committed_)
        {
            for (const VoxelKey& key : keys_)
            {
                voxelGrid_(key.x, key.y, key.z).isFilled = true;
            }
        }
    }

    void commit() noexcept
    {
        committed_ = true;
    }

private:
    VoxelGrid& voxelGrid_;
    std::vector<VoxelKey> keys_;
    bool committed_;
};


BalanceObjectives makeSingleBalanceObjectives(
    const Eigen::Vector3d& targetCenterOfMass,
    const Eigen::Vector3d& gravity
)
{
    BalanceObjectives objectives;
    objectives.emplace_back(
        targetCenterOfMass,
        gravity
    );
    return objectives;
}


BalanceObjectives validateAndCopyBalanceObjectives(
    const BalanceObjectives& objectives
)
{
    if (objectives.empty())
    {
        throw std::invalid_argument(
            "at least one balance objective is required"
        );
    }

    return objectives;
}


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
    const BalanceObjectives& objectives
)
    :
    voxelGrid_(voxelGrid),
    objectives_(validateAndCopyBalanceObjectives(objectives)),
    currCOM_(Eigen::Vector3d::Zero()),
    initialMassProperties_(
        computeSurfaceMassProperties(V, F)
    ),
    currentMassProperties_(initialMassProperties_),
    carvedVoxelCount_(0),
    iterations_(0),
    minimumRemovalsToContinue_(1),
    terminalReason_(CarvingStopReason::None),
    diagnosticsEnabled_(false),
    diagnosticsInitialized_(false),
    initialInsideVoxelCount_(0),
    initialCarveableVoxelCount_(0),
    diagnosticsHistory_(),
    preCommitHookForTesting_()
{
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


CarvingSolver::CarvingSolver(
    VoxelGrid& voxelGrid,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& planeV,
    const Eigen::MatrixXi& planeF,
    const Eigen::Vector3d& gravity,
    const Eigen::Vector3d& targetCenterOfMass
)
    : CarvingSolver(
        voxelGrid,
        V,
        F,
        makeSingleBalanceObjectives(
            targetCenterOfMass,
            gravity
        )
    )
{
    static_cast<void>(planeV);
    static_cast<void>(planeF);
}


double CarvingSolver::getCoMDistance() const
{
    return aggregateProjectedCenterOfMassDistance(
        currCOM_,
        objectives_
    );
}


double CarvingSolver::getCoMEnergy() const
{
    return aggregateProjectedCenterOfMassEnergy(
        currCOM_,
        objectives_
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
    if (terminalReason_ != CarvingStopReason::None)
    {
        return terminalReason_;
    }

    if (iterations_ >= MAX_ITER)
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
    if (isOptimized())
    {
        throw std::logic_error(
            "carveStep requires a non-terminal solver"
        );
    }

    ensureDiagnosticsInitialized();

    Eigen::MatrixXd plannedPlaneVertices;
    Eigen::MatrixXi plannedPlaneFaces;

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

    const int nextIteration = iterations_ + 1;
    diagnostics.iteration = nextIteration;

    const Eigen::Vector3i& dimensions =
        voxelGrid_.dimensions();
    const Eigen::Vector3d& voxelSize =
        voxelGrid_.voxelSize();

    const CuttingPlanes planes = makeCuttingPlanes(
        currCOM_,
        objectives_
    );
    buildDebugPlaneMesh(
        planes,
        plannedPlaneVertices,
        plannedPlaneFaces
    );

    PlaneCarvingCandidates candidates;

    for (int x = 1; x < dimensions(0) - 1; ++x)
    {
        for (int y = 1; y < dimensions(1) - 1; ++y)
        {
            for (int z = 1; z < dimensions(2) - 1; ++z)
            {
                const Voxel& voxel =
                    voxelGrid_(x, y, z);

                if (!voxel.isInner || !voxel.isFilled)
                {
                    continue;
                }

                const double score =
                    aggregatePositivePlaneScore(
                        voxel.center,
                        planes
                    );

                if (score <= 0.0)
                {
                    continue;
                }

                PlaneCarvingCandidate candidate;
                candidate.key = {x, y, z};
                candidate.center = voxel.center;
                candidate.massProperties =
                    computeAxisAlignedBoxMassProperties(
                        voxel.center,
                        voxelSize
                    );
                candidate.score = score;
                candidates.push_back(candidate);
            }
        }
    }

    sortPlaneCarvingCandidates(candidates);

    if (candidates.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error(
            "carving candidate count exceeds diagnostic capacity"
        );
    }

    diagnostics.candidateVoxelCount =
        static_cast<int>(candidates.size());

    const PlaneCarvingPrefixResult prefixResult =
        selectBestCarvingPrefix(
            currentMassProperties_,
            objectives_,
            candidates
        );

    if (
        prefixResult.evaluatedPrefixCount >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        prefixResult.bestPrefixLength >
            static_cast<std::size_t>(std::numeric_limits<int>::max())
    )
    {
        throw std::overflow_error(
            "carving prefix count exceeds diagnostic capacity"
        );
    }

    diagnostics.evaluatedPrefixCount =
        static_cast<int>(prefixResult.evaluatedPrefixCount);
    diagnostics.bestPrefixLength =
        static_cast<int>(prefixResult.bestPrefixLength);
    diagnostics.bestPrefixEnergy =
        prefixResult.bestEnergy;

    const std::size_t committedPrefixLength =
        prefixResult.bestPrefixLength;

    StagedVoxelRemoval stagedRemoval(
        voxelGrid_,
        candidates,
        committedPrefixLength
    );

    const MassProperties committedMassProperties =
        prefixResult.bestMassProperties;
    const Eigen::Vector3d committedCenterOfMass =
        prefixResult.bestCenterOfMass;

    const int carvedThisStep =
        static_cast<int>(committedPrefixLength);
    const int committedCarvedVoxelCount =
        carvedVoxelCount_ + carvedThisStep;

    CarvingStopReason proposedTerminalReason =
        CarvingStopReason::None;

    if (candidates.empty())
    {
        proposedTerminalReason =
            CarvingStopReason::NoCandidatesAvailable;
    }
    else if (
        carvedThisStep < minimumRemovalsToContinue_
    )
    {
        proposedTerminalReason =
            CarvingStopReason::TooFewVoxelsRemoved;
    }

    const RecomputedCarvingState recomputedState =
        recomputeCarvingState();
    const VoxelTopologyReport topology =
        analyzeVoxelTopology(voxelGrid_);

    diagnostics.centerOfMassAfter = committedCenterOfMass;
    diagnostics.recomputedCenterOfMassAfter =
        recomputedState.centerOfMass;
    diagnostics.firstMomentAfter =
        committedMassProperties.firstMoment;
    diagnostics.recomputedFirstMomentAfter =
        recomputedState.massProperties.firstMoment;
    diagnostics.projectedDistanceAfter =
        aggregateProjectedCenterOfMassDistance(
            committedCenterOfMass,
            objectives_
        );
    diagnostics.projectedEnergyAfter =
        aggregateProjectedCenterOfMassEnergy(
            committedCenterOfMass,
            objectives_
        );
    diagnostics.centerOfMassDiscrepancyAfter =
        (committedCenterOfMass - recomputedState.centerOfMass).norm();
    diagnostics.firstMomentDiscrepancyAfter =
        (
            committedMassProperties.firstMoment -
            recomputedState.massProperties.firstMoment
        ).norm();
    diagnostics.remainingMassAfter =
        committedMassProperties.mass;
    diagnostics.recomputedRemainingMassAfter =
        recomputedState.massProperties.mass;
    diagnostics.massDiscrepancyAfter =
        std::abs(
            committedMassProperties.mass -
            recomputedState.massProperties.mass
        );
    diagnostics.carvedVoxelCountThisStep =
        carvedThisStep;
    diagnostics.totalCarvedVoxelCount =
        committedCarvedVoxelCount;
    diagnostics.initialInsideVoxelCount =
        initialInsideVoxelCount_;
    diagnostics.initialCarveableVoxelCount =
        initialCarveableVoxelCount_;
    diagnostics.remainingInsideVoxelCount =
        recomputedState.filledInsideVoxelCount;
    diagnostics.materialComponentCountAfter =
        topology.materialComponentCount;
    diagnostics.carvedVoidComponentCountAfter =
        topology.carvedVoidComponentCount;
    const double totalGridVoxelCount =
        static_cast<double>(dimensions(0)) *
        static_cast<double>(dimensions(1)) *
        static_cast<double>(dimensions(2));
    diagnostics.carvedGridRatio =
        committedCarvedVoxelCount / totalGridVoxelCount;
    diagnostics.carvedCarveableRatio =
        initialCarveableVoxelCount_ > 0 ?
            committedCarvedVoxelCount /
                static_cast<double>(initialCarveableVoxelCount_) :
            0.0;
    diagnostics.stopReason = proposedTerminalReason;

    diagnosticsHistory_.reserve(
        diagnosticsHistory_.size() + 1
    );

    if (preCommitHookForTesting_)
    {
        preCommitHookForTesting_();
    }

    /*
     * Everything below this point is non-allocating under the reserved
     * history capacity. Commit scalar state and staged flags as one step.
     */
    currentMassProperties_ = committedMassProperties;
    currCOM_ = committedCenterOfMass;
    carvedVoxelCount_ = committedCarvedVoxelCount;
    terminalReason_ = proposedTerminalReason;
    iterations_ = nextIteration;
    stagedRemoval.commit();
    diagnosticsHistory_.push_back(diagnostics);

    carvePlaneV.swap(plannedPlaneVertices);
    carvePlaneF.swap(plannedPlaneFaces);
    centerOfMass = currCOM_;

    if (diagnosticsEnabled_)
    {
        try
        {
            printCarvingStepDiagnostics(
                std::cout,
                diagnostics
            );
        }
        catch (const std::ios_base::failure&)
        {
            // Optional logging cannot roll back an already committed step.
        }
    }
}


void CarvingSolver::buildDebugPlaneMesh(
    const CuttingPlanes& planes,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
) const
{
    std::size_t activePlaneCount = 0;

    for (const CuttingPlane& plane : planes)
    {
        if (plane.active)
        {
            ++activePlaneCount;
        }
    }

    vertices.resize(
        static_cast<Eigen::Index>(4 * activePlaneCount),
        3
    );
    faces.resize(
        static_cast<Eigen::Index>(2 * activePlaneCount),
        3
    );

    if (activePlaneCount == 0)
    {
        return;
    }

    const double extent =
        voxelGrid_.dimensions().maxCoeff() *
        voxelGrid_.voxelSize().maxCoeff() *
        10.0;

    Eigen::Index planeIndex = 0;

    for (const CuttingPlane& plane : planes)
    {
        if (!plane.active)
        {
            continue;
        }

        const Eigen::Vector3d firstAxis =
            plane.gravity.normalized();
        const Eigen::Vector3d planeNormal =
            normalizedGravityDirection(plane.normal);
        const Eigen::Vector3d secondAxis =
            planeNormal.cross(firstAxis).normalized();
        const Eigen::Index vertexOffset =
            4 * planeIndex;
        const Eigen::Index faceOffset =
            2 * planeIndex;

        vertices.row(vertexOffset + 0) =
            (plane.point - extent * firstAxis - extent * secondAxis).transpose();
        vertices.row(vertexOffset + 1) =
            (plane.point + extent * firstAxis - extent * secondAxis).transpose();
        vertices.row(vertexOffset + 2) =
            (plane.point - extent * firstAxis + extent * secondAxis).transpose();
        vertices.row(vertexOffset + 3) =
            (plane.point + extent * firstAxis + extent * secondAxis).transpose();

        faces.row(faceOffset + 0) <<
            vertexOffset + 0,
            vertexOffset + 1,
            vertexOffset + 2;
        faces.row(faceOffset + 1) <<
            vertexOffset + 2,
            vertexOffset + 1,
            vertexOffset + 3;

        ++planeIndex;
    }
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
    return objectives_.front().target();
}


const BalanceObjectives& CarvingSolver::balanceObjectives() const
{
    return objectives_;
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


void CarvingSolver::setMinimumRemovalsToContinue(
    int count
)
{
    if (count < 1)
    {
        throw std::invalid_argument(
            "minimum removals to continue must be at least one"
        );
    }

    if (iterations_ != 0)
    {
        throw std::logic_error(
            "minimum removals must be configured before carving"
        );
    }

    minimumRemovalsToContinue_ = count;
}


int CarvingSolver::minimumRemovalsToContinue() const
{
    return minimumRemovalsToContinue_;
}


void CarvingSolver::setPreCommitHookForTesting(
    const std::function<void()>& hook
)
{
    if (iterations_ != 0)
    {
        throw std::logic_error(
            "pre-commit test hook must be configured before carving"
        );
    }

    preCommitHookForTesting_ = hook;
}
