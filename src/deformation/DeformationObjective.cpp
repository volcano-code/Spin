#include "deformation/DeformationObjective.h"

#include "extraction/VoidMeshExtractor.h"
#include "physics/BalanceEnergy.h"
#include "physics/SurfaceMassDifferential.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>


namespace deformation
{
DeformationObjective::CompactCavityBinding
DeformationObjective::makeCompactCavityBinding(
    const DeformationSystem& system,
    const VoxelGrid& fixedVoxelFill
)
{
    Eigen::MatrixXi globalFaces;
    VoidMeshExtractor::extractFaces(fixedVoxelFill, globalFaces);

    CompactCavityBinding binding;
    binding.faces.resize(globalFaces.rows(), 3);

    if (globalFaces.rows() == 0)
    {
        binding.restVertices.resize(0, 3);
        binding.weights.resize(
            0,
            system.fullCornerWeights().cols()
        );
        return binding;
    }

    std::vector<int> globals(
        globalFaces.data(),
        globalFaces.data() + globalFaces.size()
    );
    std::sort(globals.begin(), globals.end());
    globals.erase(
        std::unique(globals.begin(), globals.end()),
        globals.end()
    );

    binding.restVertices.resize(
        static_cast<Eigen::Index>(globals.size()),
        3
    );
    binding.weights.resize(
        static_cast<Eigen::Index>(globals.size()),
        system.fullCornerWeights().cols()
    );

    std::vector<int> globalToCompact(
        static_cast<std::size_t>(
            system.domain().fullRestCorners().rows()
        ),
        -1
    );

    for (std::size_t compact = 0;
         compact < globals.size();
         ++compact)
    {
        const int global = globals[compact];
        if (
            global < 0 ||
            global >= system.domain().fullRestCorners().rows()
        )
        {
            throw std::logic_error(
                "cavity face references an invalid global grid corner"
            );
        }

        binding.restVertices.row(
            static_cast<Eigen::Index>(compact)
        ) = system.domain().fullRestCorners().row(global);
        binding.weights.row(
            static_cast<Eigen::Index>(compact)
        ) = system.fullCornerWeights().row(global);
        globalToCompact[static_cast<std::size_t>(global)] =
            static_cast<int>(compact);
    }

    for (Eigen::Index index = 0;
         index < globalFaces.size();
         ++index)
    {
        const int compact = globalToCompact[
            static_cast<std::size_t>(globalFaces(index))
        ];
        if (compact < 0)
        {
            throw std::logic_error(
                "cavity face compaction lost a referenced corner"
            );
        }
        binding.faces(index) = compact;
    }

    return binding;
}


DeformationObjective::DeformationObjective(
    const DeformationSystem& system,
    const VoxelGrid& fixedVoxelFill,
    const BalanceObjectives& objectives,
    double shapeScaleLambda,
    double density,
    double minimumPositiveMass
)
    : DeformationObjective(
        system,
        fixedVoxelFill,
        objectives,
        shapeScaleLambda,
        density,
        minimumPositiveMass,
        makeCompactCavityBinding(system, fixedVoxelFill)
    )
{
}


DeformationObjective::DeformationObjective(
    const DeformationSystem& system,
    const VoxelGrid& fixedVoxelFill,
    const BalanceObjectives& objectives,
    double shapeScaleLambda,
    double density,
    double minimumPositiveMass,
    const CompactCavityBinding& compactCavity
)
    :
    bindingHandles_(system.bindingHandles()),
    outerModel_(
        system.outerRestVertices(),
        system.outerWeights(),
        system.bindingHandles()
    ),
    innerModel_(
        compactCavity.restVertices,
        compactCavity.weights,
        system.bindingHandles()
    ),
    outerFaces_(system.outerFaces()),
    innerFaces_(compactCavity.faces),
    laplacianEnergy_(
        system.outerRestVertices(),
        system.outerFaces()
    ),
    objectives_(objectives),
    shapeScaleLambda_(shapeScaleLambda),
    density_(density),
    minimumPositiveMass_(minimumPositiveMass)
{
    if (objectives_.empty())
    {
        throw std::invalid_argument(
            "deformation objective requires at least one balance target"
        );
    }

    if (
        !std::isfinite(shapeScaleLambda_) ||
        shapeScaleLambda_ < 0.0 ||
        !std::isfinite(density_) ||
        density_ <= 0.0 ||
        !std::isfinite(minimumPositiveMass_) ||
        minimumPositiveMass_ < 0.0
    )
    {
        throw std::invalid_argument(
            "deformation objective scalars must be finite and physically valid"
        );
    }

    /* Reuse R5's immutable-grid contract without materializing a snapshot. */
    system.validateVoxelFill(fixedVoxelFill);
}


Eigen::VectorXd DeformationObjective::initialParameters() const
{
    return outerModel_.initialParameters();
}


HandleSet DeformationObjective::makeHandles(
    const Eigen::VectorXd& parameters
) const
{
    HandleSet handles = bindingHandles_;
    handles.setOptimizationParameters(parameters);
    return handles;
}


DeformationEnergyEvaluation DeformationObjective::evaluate(
    const Eigen::VectorXd& parameters,
    double mu
) const
{
    if (!std::isfinite(mu) || mu < 0.0 || mu > 1.0)
    {
        throw std::invalid_argument(
            "deformation tradeoff mu must lie in [0,1]"
        );
    }

    if (
        parameters.size() != bindingHandles_.optimizableParameterCount() ||
        !parameters.allFinite()
    )
    {
        throw std::invalid_argument(
            "deformation parameter vector is finite but has the wrong contract"
        );
    }
    for (Eigen::Index offset = 3;
         offset < parameters.size();
         offset += 4)
    {
        if (!(parameters(offset) > 0.0))
        {
            throw InfeasibleDeformationState(
                "a deformation trial produced a non-positive handle scale"
            );
        }
    }

    /* Validates dimensions, finiteness and strictly positive free scales. */
    static_cast<void>(makeHandles(parameters));

    const Eigen::MatrixXd outerVertices =
        outerModel_.evaluate(parameters);
    const Eigen::MatrixXd innerVertices =
        innerModel_.evaluate(parameters);

    CenterOfMassDifferential massDifferential;
    try
    {
        massDifferential = computeCenterOfMassDifferential(
            outerVertices,
            outerFaces_,
            innerVertices,
            innerFaces_,
            density_,
            minimumPositiveMass_
        );
    }
    catch (const std::domain_error& error)
    {
        throw InfeasibleDeformationState(error.what());
    }

    DeformationEnergyEvaluation result;
    result.mass = massDifferential.properties.mass;
    result.centerOfMass = massDifferential.centerOfMass;
    result.centerOfMassEnergy =
        aggregateProjectedCenterOfMassEnergy(
            result.centerOfMass,
            objectives_
        );
    result.unscaledShapeEnergy =
        laplacianEnergy_.energy(outerVertices);
    result.shapeEnergy =
        shapeScaleLambda_ * result.unscaledShapeEnergy;
    result.totalEnergy =
        (1.0 - mu) * result.centerOfMassEnergy +
        mu * result.shapeEnergy;

    Eigen::Vector3d centerEnergyGradient =
        Eigen::Vector3d::Zero();
    for (const BalanceObjective& objective : objectives_)
    {
        centerEnergyGradient += projectPerpendicularToGravity(
            result.centerOfMass - objective.target(),
            objective.gravity()
        );
    }

    const Eigen::MatrixXd centerParameterJacobian =
        massDifferential.outerJacobian * outerModel_.jacobian() +
        massDifferential.innerJacobian * innerModel_.jacobian();

    const Eigen::VectorXd centerGradient =
        centerParameterJacobian.transpose() * centerEnergyGradient;
    const Eigen::VectorXd shapeGradient =
        shapeScaleLambda_ * outerModel_.jacobian().transpose() *
        laplacianEnergy_.vertexGradient(outerVertices);

    result.gradient =
        (1.0 - mu) * centerGradient +
        mu * shapeGradient;

    if (
        !std::isfinite(result.totalEnergy) ||
        !std::isfinite(result.centerOfMassEnergy) ||
        !std::isfinite(result.shapeEnergy) ||
        !std::isfinite(result.unscaledShapeEnergy) ||
        !result.gradient.allFinite()
    )
    {
        throw std::overflow_error(
            "deformation objective evaluation overflowed"
        );
    }

    return result;
}


double DeformationObjective::shapeScaleLambda() const noexcept
{
    return shapeScaleLambda_;
}


const BalanceObjectives&
DeformationObjective::objectives() const noexcept
{
    return objectives_;
}


} // namespace deformation
