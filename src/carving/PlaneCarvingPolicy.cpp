#include "carving/PlaneCarvingPolicy.h"

#include "physics/BalanceEnergy.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace
{

bool keyLess(
    const VoxelKey& first,
    const VoxelKey& second
)
{
    if (first.x != second.x)
    {
        return first.x < second.x;
    }

    if (first.y != second.y)
    {
        return first.y < second.y;
    }

    return first.z < second.z;
}


bool strictlyLowerEnergy(
    double candidate,
    double incumbent
)
{
    const double scale = std::max(
        std::numeric_limits<double>::min(),
        std::max(std::abs(candidate), std::abs(incumbent))
    );

    return candidate < incumbent -
        64.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace


PlaneCarvingCandidate::PlaneCarvingCandidate()
    :
    key{0, 0, 0},
    center(Eigen::Vector3d::Zero()),
    massProperties(),
    score(0.0)
{
}


CuttingPlane::CuttingPlane()
    :
    point(Eigen::Vector3d::Zero()),
    normal(Eigen::Vector3d::Zero()),
    gravity(Eigen::Vector3d::Zero()),
    active(false)
{
}


PlaneCarvingPrefixResult::PlaneCarvingPrefixResult()
    :
    evaluatedPrefixCount(0),
    bestPrefixLength(0),
    initialEnergy(0.0),
    bestEnergy(0.0),
    bestMassProperties(),
    bestCenterOfMass(Eigen::Vector3d::Zero())
{
}


CuttingPlanes makeCuttingPlanes(
    const Eigen::Vector3d& currentCenterOfMass,
    const BalanceObjectives& objectives
)
{
    if (!currentCenterOfMass.allFinite())
    {
        throw std::invalid_argument(
            "cutting-plane center of mass must be finite"
        );
    }

    if (objectives.empty())
    {
        throw std::invalid_argument(
            "at least one balance objective is required"
        );
    }

    CuttingPlanes planes;
    planes.reserve(objectives.size());

    for (const BalanceObjective& objective : objectives)
    {
        CuttingPlane plane;
        plane.point = objective.target();
        plane.gravity = objective.gravity();
        plane.normal = projectPerpendicularToGravity(
            currentCenterOfMass - objective.target(),
            objective.gravity()
        );

        if (!plane.normal.allFinite())
        {
            throw std::overflow_error(
                "cutting-plane normal overflowed"
            );
        }

        plane.active =
            plane.normal.cwiseAbs().maxCoeff() > 0.0;
        planes.push_back(plane);
    }

    return planes;
}


double aggregatePositivePlaneScore(
    const Eigen::Vector3d& point,
    const CuttingPlanes& planes
)
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "candidate point must be finite"
        );
    }

    double score = 0.0;

    for (const CuttingPlane& plane : planes)
    {
        if (!plane.active)
        {
            continue;
        }

        const double distance =
            (point - plane.point).dot(plane.normal);

        if (!std::isfinite(distance))
        {
            throw std::overflow_error(
                "cutting-plane score overflowed"
            );
        }

        if (distance > 0.0)
        {
            if (
                score >
                    std::numeric_limits<double>::max() -
                    distance
            )
            {
                throw std::overflow_error(
                    "aggregate cutting-plane score overflowed"
                );
            }

            score += distance;
        }
    }

    return score;
}


void sortPlaneCarvingCandidates(
    PlaneCarvingCandidates& candidates
)
{
    for (const PlaneCarvingCandidate& candidate : candidates)
    {
        if (
            !candidate.center.allFinite() ||
            !candidate.massProperties.isFinite() ||
            candidate.massProperties.mass <= 0.0 ||
            !std::isfinite(candidate.score) ||
            candidate.score <= 0.0
        )
        {
            throw std::invalid_argument(
                "carving candidate must have finite positive mass and score"
            );
        }
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](
            const PlaneCarvingCandidate& first,
            const PlaneCarvingCandidate& second
        )
        {
            if (first.score != second.score)
            {
                return first.score > second.score;
            }

            return keyLess(first.key, second.key);
        }
    );
}


PlaneCarvingPrefixResult selectBestCarvingPrefix(
    const MassProperties& initialMassProperties,
    const BalanceObjectives& objectives,
    const PlaneCarvingCandidates& orderedCandidates
)
{
    Eigen::Vector3d initialCenterOfMass =
        Eigen::Vector3d::Zero();

    if (
        !initialMassProperties.isFinite() ||
        initialMassProperties.mass <= 0.0 ||
        !initialMassProperties.tryCenterOfMass(
            initialCenterOfMass,
            0.0
        )
    )
    {
        throw std::invalid_argument(
            "initial carving mass properties must be finite and positive"
        );
    }

    PlaneCarvingPrefixResult result;
    result.evaluatedPrefixCount = 1;
    result.initialEnergy = aggregateProjectedCenterOfMassEnergy(
        initialCenterOfMass,
        objectives
    );
    result.bestEnergy = result.initialEnergy;
    result.bestMassProperties = initialMassProperties;
    result.bestCenterOfMass = initialCenterOfMass;

    MassProperties trialMassProperties =
        initialMassProperties;

    for (
        std::size_t candidateIndex = 0;
        candidateIndex < orderedCandidates.size();
        ++candidateIndex
    )
    {
        const PlaneCarvingCandidate& candidate =
            orderedCandidates[candidateIndex];

        if (
            !candidate.massProperties.isFinite() ||
            candidate.massProperties.mass <= 0.0
        )
        {
            throw std::invalid_argument(
                "candidate mass properties must be finite and positive"
            );
        }

        trialMassProperties -= candidate.massProperties;

        Eigen::Vector3d trialCenterOfMass =
            Eigen::Vector3d::Zero();

        if (
            trialMassProperties.mass <= 0.0 ||
            !trialMassProperties.tryCenterOfMass(
                trialCenterOfMass,
                0.0
            )
        )
        {
            break;
        }

        const double energy =
            aggregateProjectedCenterOfMassEnergy(
                trialCenterOfMass,
                objectives
            );

        ++result.evaluatedPrefixCount;

        if (strictlyLowerEnergy(energy, result.bestEnergy))
        {
            result.bestPrefixLength = candidateIndex + 1;
            result.bestEnergy = energy;
            result.bestMassProperties = trialMassProperties;
            result.bestCenterOfMass = trialCenterOfMass;
        }
    }

    return result;
}
