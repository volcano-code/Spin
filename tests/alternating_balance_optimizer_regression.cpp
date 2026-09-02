#include "carving/CarvingSolver.h"
#include "deformation/DeformationSystem.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "physics/MassProperties.h"
#include "optimization/AlternatingBalanceOptimizer.h"
#include "stability/BalanceTargetSchedule.h"
#include "stability/SupportRegion.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


namespace
{


constexpr double NUMERICAL_TOLERANCE = 1.0e-10;


int fail(const std::string& scenario, const std::string& reason)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


bool close(double first, double second, double tolerance)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <=
            tolerance *
            std::max({1.0, std::abs(first), std::abs(second)});
}


void makeBox(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        minimum(0), minimum(1), minimum(2),
        maximum(0), minimum(1), minimum(2),
        maximum(0), maximum(1), minimum(2),
        minimum(0), maximum(1), minimum(2),
        minimum(0), minimum(1), maximum(2),
        maximum(0), minimum(1), maximum(2),
        maximum(0), maximum(1), maximum(2),
        minimum(0), maximum(1), maximum(2);

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


/** A closed overhanging prism whose initial CoM lies beyond its fixed base. */
void makeOverhangingPrism(
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        vertices,
        faces
    );

    /* Bottom support spans x=[0,0.25], while the top spans x=[0.75,1]. */
    vertices(1, 0) = 0.25;
    vertices(2, 0) = 1.00;
    vertices(3, 0) = 0.75;
    vertices(5, 0) = 0.25;
    vertices(6, 0) = 1.00;
    vertices(7, 0) = 0.75;
}


deformation::Handle makeSupportHandle()
{
    deformation::Handle handle;
    handle.role = deformation::HandleRole::Support;
    handle.locked = true;
    return handle;
}


deformation::Handle makeUserHandle(const Eigen::Vector3d& pivot)
{
    deformation::Handle handle;
    handle.restPivot = pivot;
    handle.role = deformation::HandleRole::User;
    handle.locked = false;
    return handle;
}


struct GeometryBinding
{
    deformation::HandleSet handles;
    Eigen::MatrixXd weights;
};


GeometryBinding makeVerticalBinding(
    const deformation::GridDeformationDomain& domain
)
{
    const double minimum = domain.restNodes().col(1).minCoeff();
    const double maximum = domain.restNodes().col(1).maxCoeff();
    const double span = maximum - minimum;

    deformation::Handle support = makeSupportHandle();
    deformation::Handle user = makeUserHandle(
        Eigen::Vector3d(0.875, maximum, 0.5)
    );
    Eigen::MatrixXd weights(domain.nodeCount(), 2);

    for (Eigen::Index node = 0; node < domain.nodeCount(); ++node)
    {
        const double y = domain.restNodes()(node, 1);
        const double userWeight = (y - minimum) / span;
        weights(node, 0) = 1.0 - userWeight;
        weights(node, 1) = userWeight;

        if (close(y, minimum, 1.0e-12))
        {
            support.constrainedNodes.push_back(static_cast<int>(node));
        }
        if (close(y, maximum, 1.0e-12))
        {
            user.constrainedNodes.push_back(static_cast<int>(node));
        }
    }

    return GeometryBinding{
        deformation::HandleSet(
            deformation::HandleVector{support, user}
        ),
        weights
    };
}


std::vector<unsigned char> fillMask(const VoxelGrid& grid)
{
    std::vector<unsigned char> result;
    result.reserve(
        static_cast<std::size_t>(grid.dimensions().prod())
    );

    for (int x = 0; x < grid.dimensions()(0); ++x)
    {
        for (int y = 0; y < grid.dimensions()(1); ++y)
        {
            for (int z = 0; z < grid.dimensions()(2); ++z)
            {
                result.push_back(
                    grid(x, y, z).isFilled ? 1U : 0U
                );
            }
        }
    }
    return result;
}


bool supportNodesFixed(
    const deformation::DeformationSystem& system,
    const Eigen::MatrixXd& deformedCorners,
    double tolerance
)
{
    const deformation::GridDeformationDomain& domain = system.domain();
    const deformation::HandleSet& handles = system.bindingHandles();

    for (std::size_t handleIndex = 0;
         handleIndex < handles.size();
         ++handleIndex)
    {
        if (handles[handleIndex].role != deformation::HandleRole::Support)
        {
            continue;
        }

        for (const int local : handles[handleIndex].constrainedNodes)
        {
            const int global = domain.localToGlobalCorner()(local);
            if (!deformedCorners.row(global).isApprox(
                domain.restNodes().row(local),
                tolerance
            ))
            {
                return false;
            }
        }
    }
    return true;
}


MassProperties completeSurfaceMassOracle(
    const deformation::DeformedGeometry& geometry
)
{
    return
        computeSurfaceMassProperties(
            geometry.outerVertices,
            geometry.outerFaces,
            1.0
        ) +
        computeSurfaceMassProperties(
            geometry.innerVertices,
            geometry.innerFaces,
            1.0
        );
}


struct ProductionFixture
{
    Eigen::MatrixXd restVertices;
    Eigen::MatrixXi faces;
    VoxelGrid grid;
    deformation::GridDeformationDomain domain;
    GeometryBinding binding;
    deformation::DeformationSystem system;
    stability::SupportRegion support;

    ProductionFixture(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& inputFaces,
        int resolution
    )
        :
        restVertices(vertices),
        faces(inputFaces),
        grid(restVertices, faces, 0.0, resolution),
        domain(deformation::GridDeformationDomain::fromVoxelGrid(grid)),
        binding(makeVerticalBinding(domain)),
        system(
            grid,
            restVertices,
            faces,
            binding.handles,
            binding.weights,
            true
        ),
        support(stability::SupportRegion::fromMesh(
            restVertices,
            Eigen::Vector3d(0.0, -1.0, 0.0),
            0.0
        ))
    {
    }
};


stability::BalanceTargetSchedule makeStandingSchedule(
    const stability::SupportRegion& support,
    double gammaRadians
)
{
    stability::BalanceTargetSpecifications specifications;
    specifications.emplace_back(
        stability::StandingTargetSpecification(
            support,
            gammaRadians
        )
    );
    return stability::BalanceTargetSchedule(specifications);
}


optimization::AlternatingBalanceOptimizerConfig makeConfiguration(
    int maximumCarvingBlocks,
    double gradientTolerance
)
{
    optimization::AlternatingBalanceOptimizerConfig config;
    config.maximumCarvingBlocks = maximumCarvingBlocks;
    config.maximumDeformationIterationsPerBlock = 200;
    config.gradientTolerance = gradientTolerance;
    config.shapeScaleLambda = 0.01;
    config.lineSearch.initialStep = 1.0;
    config.lineSearch.minimumStep = 0.4;
    config.lineSearch.reductionRatio = 0.8;
    config.lineSearch.requiredRelativeDecrease = 0.01;
    config.lineSearch.muDecrement = 0.05;
    return config;
}


bool sameMassProperties(
    const MassProperties& first,
    const MassProperties& second,
    double tolerance
)
{
    return
        close(first.mass, second.mass, tolerance) &&
        first.firstMoment.isApprox(second.firstMoment, tolerance);
}


bool sameResult(
    const optimization::AlternatingBalanceOptimizationResult& first,
    const optimization::AlternatingBalanceOptimizationResult& second
)
{
    if (
        first.status != second.status ||
        first.mu != second.mu ||
        first.geometryEpoch != second.geometryEpoch ||
        first.geometryHash != second.geometryHash ||
        first.usesHypotheticalFlatteningSupport !=
            second.usesHypotheticalFlatteningSupport ||
        first.targets.generation != second.targets.generation ||
        first.targets.allReached != second.targets.allReached ||
        first.carvingBlocks.size() != second.carvingBlocks.size() ||
        first.carvingIterations.size() != second.carvingIterations.size() ||
        first.deformationIterations.size() !=
            second.deformationIterations.size() ||
        !first.parameters.isApprox(second.parameters, 0.0) ||
        fillMask(first.voxelFill) != fillMask(second.voxelFill) ||
        !first.geometry.outerVertices.isApprox(
            second.geometry.outerVertices,
            0.0
        ) ||
        !first.geometry.gridCorners.isApprox(
            second.geometry.gridCorners,
            0.0
        ) ||
        first.geometry.innerFaces != second.geometry.innerFaces ||
        !sameMassProperties(
            first.massProperties,
            second.massProperties,
            0.0
        ) ||
        !first.centerOfMass.isApprox(second.centerOfMass, 0.0)
    )
    {
        return false;
    }

    for (std::size_t index = 0;
         index < first.carvingBlocks.size();
         ++index)
    {
        const optimization::AlternatingCarvingRecord& a =
            first.carvingBlocks[index];
        const optimization::AlternatingCarvingRecord& b =
            second.carvingBlocks[index];
        if (
            a.block != b.block ||
            a.geometryEpochBeforeReset != b.geometryEpochBeforeReset ||
            a.geometryEpochAfterReset != b.geometryEpochAfterReset ||
            a.geometryEpoch != b.geometryEpoch ||
            a.refilledVoxelCount != b.refilledVoxelCount ||
            a.alphaHammingDistance != b.alphaHammingDistance ||
            a.carvingIterations != b.carvingIterations ||
            a.carvedVoxelCount != b.carvedVoxelCount ||
            a.stopReason != b.stopReason ||
            a.resetTargetEnergy != b.resetTargetEnergy ||
            a.finalTargetEnergy != b.finalTargetEnergy ||
            a.resetTargetGenerationBefore !=
                b.resetTargetGenerationBefore ||
            a.resetTargetGenerationAfter !=
                b.resetTargetGenerationAfter ||
            a.resetObjectiveContextChanged !=
                b.resetObjectiveContextChanged ||
            a.resetAllTargetsReached != b.resetAllTargetsReached ||
            a.targetGenerationBefore != b.targetGenerationBefore ||
            a.targetGenerationAfter != b.targetGenerationAfter ||
            a.objectiveContextChanged != b.objectiveContextChanged ||
            a.allTargetsReached != b.allTargetsReached ||
            a.alphaHashBeforeReset != b.alphaHashBeforeReset ||
            a.alphaHashAfterReset != b.alphaHashAfterReset ||
            a.alphaHash != b.alphaHash ||
            a.geometryHash != b.geometryHash ||
            !sameMassProperties(a.massProperties, b.massProperties, 0.0) ||
            !a.centerOfMass.isApprox(b.centerOfMass, 0.0)
        )
        {
            return false;
        }
    }

    for (std::size_t index = 0;
         index < first.carvingIterations.size();
         ++index)
    {
        const optimization::AlternatingCarvingIterationRecord& a =
            first.carvingIterations[index];
        const optimization::AlternatingCarvingIterationRecord& b =
            second.carvingIterations[index];
        if (
            a.block != b.block ||
            a.iteration != b.iteration ||
            a.geometryEpochBefore != b.geometryEpochBefore ||
            a.geometryEpochAfter != b.geometryEpochAfter ||
            a.carvedVoxelCount != b.carvedVoxelCount ||
            a.totalCarvedVoxelCount != b.totalCarvedVoxelCount ||
            a.stopReason != b.stopReason ||
            a.frozenTargetEnergyBefore != b.frozenTargetEnergyBefore ||
            a.frozenTargetEnergyAfter != b.frozenTargetEnergyAfter ||
            a.refreshedTargetEnergyAfter != b.refreshedTargetEnergyAfter ||
            a.targetGenerationBefore != b.targetGenerationBefore ||
            a.targetGenerationAfter != b.targetGenerationAfter ||
            a.objectiveContextChanged != b.objectiveContextChanged ||
            a.allTargetsReached != b.allTargetsReached ||
            a.alphaHashBefore != b.alphaHashBefore ||
            a.alphaHashAfter != b.alphaHashAfter ||
            !sameMassProperties(a.massProperties, b.massProperties, 0.0) ||
            !a.centerOfMass.isApprox(b.centerOfMass, 0.0)
        )
        {
            return false;
        }
    }

    for (std::size_t index = 0;
         index < first.deformationIterations.size();
         ++index)
    {
        const optimization::AlternatingDeformationRecord& a =
            first.deformationIterations[index];
        const optimization::AlternatingDeformationRecord& b =
            second.deformationIterations[index];
        if (
            a.block != b.block ||
            a.iteration != b.iteration ||
            a.geometryEpochBefore != b.geometryEpochBefore ||
            a.geometryEpochAfter != b.geometryEpochAfter ||
            a.decision != b.decision ||
            a.accepted != b.accepted ||
            a.rolledBack != b.rolledBack ||
            a.energyBefore != b.energyBefore ||
            a.trialEnergy != b.trialEnergy ||
            a.energyAfter != b.energyAfter ||
            a.refreshedTargetEnergyAfter !=
                b.refreshedTargetEnergyAfter ||
            a.gradientNorm != b.gradientNorm ||
            a.attemptedStep != b.attemptedStep ||
            a.muBefore != b.muBefore ||
            a.objectiveReweighted != b.objectiveReweighted ||
            a.stateAfter.step != b.stateAfter.step ||
            a.stateAfter.mu != b.stateAfter.mu ||
            a.targetGenerationBefore != b.targetGenerationBefore ||
            a.targetGenerationAfter != b.targetGenerationAfter ||
            a.objectiveContextChanged != b.objectiveContextChanged ||
            a.allTargetsReached != b.allTargetsReached ||
            a.refillCalls != b.refillCalls ||
            a.carvingCalls != b.carvingCalls ||
            a.alphaHammingDistance != b.alphaHammingDistance ||
            a.alphaHashBefore != b.alphaHashBefore ||
            a.alphaHashAfter != b.alphaHashAfter ||
            a.geometryHashBefore != b.geometryHashBefore ||
            a.geometryHashAfter != b.geometryHashAfter ||
            !a.parametersBefore.isApprox(b.parametersBefore, 0.0) ||
            !a.parametersAfter.isApprox(b.parametersAfter, 0.0) ||
            !sameMassProperties(
                a.massPropertiesAfter,
                b.massPropertiesAfter,
                0.0
            ) ||
            !a.centerOfMassAfter.isApprox(b.centerOfMassAfter, 0.0)
        )
        {
            return false;
        }
    }
    return true;
}


optimization::AlternatingBalanceOptimizationResult runOverhang(
    const ProductionFixture& fixture,
    int maximumCarvingBlocks,
    double gradientTolerance
)
{
    const stability::BalanceTargetSchedule schedule =
        makeStandingSchedule(fixture.support, 0.08);
    const optimization::AlternatingBalanceOptimizer optimizer(
        fixture.system,
        fixture.grid,
        schedule,
        makeConfiguration(maximumCarvingBlocks, gradientTolerance)
    );
    return optimizer.optimize(
        fixture.binding.handles.optimizationParameters(),
        0.2
    );
}


bool validateCarvingTrace(
    const optimization::AlternatingBalanceOptimizationResult& result,
    std::string& reason
)
{
    for (std::size_t blockIndex = 0;
         blockIndex < result.carvingBlocks.size();
         ++blockIndex)
    {
        const optimization::AlternatingCarvingRecord& block =
            result.carvingBlocks[blockIndex];
        const int resetEpochIncrement =
            block.refilledVoxelCount > 0 ? 1 : 0;
        if (
            block.block != static_cast<int>(blockIndex) ||
            block.geometryEpochAfterReset !=
                block.geometryEpochBeforeReset + resetEpochIncrement ||
            block.resetTargetGenerationAfter <
                block.resetTargetGenerationBefore ||
            block.resetObjectiveContextChanged !=
                (block.resetTargetGenerationAfter >
                 block.resetTargetGenerationBefore) ||
            block.targetGenerationBefore !=
                block.resetTargetGenerationAfter ||
            !std::isfinite(block.resetTargetEnergy) ||
            !std::isfinite(block.finalTargetEnergy)
        )
        {
            reason = "carving block reset trace is inconsistent";
            return false;
        }

        int expectedIteration = 0;
        int expectedTotalCarved = 0;
        int expectedEpoch = block.geometryEpochAfterReset;
        std::size_t expectedGeneration =
            block.resetTargetGenerationAfter;
        std::uint64_t expectedHash = block.alphaHashAfterReset;
        double expectedRefreshedEnergy = block.resetTargetEnergy;
        bool reachedInsideBlock = block.resetAllTargetsReached;

        for (const optimization::AlternatingCarvingIterationRecord& record :
             result.carvingIterations)
        {
            if (record.block != static_cast<int>(blockIndex))
            {
                continue;
            }
            if (reachedInsideBlock)
            {
                reason = "controller carved after a target-reaching state";
                return false;
            }

            expectedTotalCarved += record.carvedVoxelCount;
            const int epochIncrement =
                record.carvedVoxelCount > 0 ? 1 : 0;
            const double tolerance =
                1.0e-12 *
                std::max(1.0, record.frozenTargetEnergyBefore);
            if (
                record.iteration != expectedIteration++ ||
                record.geometryEpochBefore != expectedEpoch ||
                record.geometryEpochAfter !=
                    expectedEpoch + epochIncrement ||
                record.targetGenerationBefore != expectedGeneration ||
                !close(
                    record.frozenTargetEnergyBefore,
                    expectedRefreshedEnergy,
                    1.0e-10
                ) ||
                record.targetGenerationAfter < expectedGeneration ||
                record.objectiveContextChanged !=
                    (record.targetGenerationAfter > expectedGeneration) ||
                record.frozenTargetEnergyAfter >
                    record.frozenTargetEnergyBefore + tolerance ||
                !std::isfinite(record.refreshedTargetEnergyAfter) ||
                record.alphaHashBefore != expectedHash ||
                (
                    record.carvedVoxelCount > 0 &&
                    record.alphaHashAfter == record.alphaHashBefore
                ) ||
                (
                    record.carvedVoxelCount == 0 &&
                    (
                        record.alphaHashAfter != record.alphaHashBefore ||
                        record.stopReason == CarvingStopReason::None
                    )
                ) ||
                record.carvedVoxelCount < 0 ||
                record.totalCarvedVoxelCount != expectedTotalCarved ||
                !record.massProperties.isFinite() ||
                !record.centerOfMass.allFinite()
            )
            {
                reason =
                    "committed carving step violated energy/refresh state";
                return false;
            }

            expectedEpoch = record.geometryEpochAfter;
            expectedGeneration = record.targetGenerationAfter;
            expectedHash = record.alphaHashAfter;
            expectedRefreshedEnergy =
                record.refreshedTargetEnergyAfter;
            reachedInsideBlock = record.allTargetsReached;
        }

        if (
            expectedIteration != block.carvingIterations ||
            expectedTotalCarved != block.carvedVoxelCount ||
            expectedEpoch != block.geometryEpoch ||
            expectedGeneration != block.targetGenerationAfter ||
            expectedHash != block.alphaHash ||
            block.geometryHash == 0U ||
            !close(
                expectedRefreshedEnergy,
                block.finalTargetEnergy,
                1.0e-10
            ) ||
            block.objectiveContextChanged !=
                (block.targetGenerationAfter >
                 block.targetGenerationBefore) ||
            block.allTargetsReached != reachedInsideBlock ||
            !block.massProperties.isFinite() ||
            !block.centerOfMass.allFinite()
        )
        {
            reason = "carving block summary disagrees with its step trace";
            return false;
        }

        if (block.allTargetsReached)
        {
            const bool hasDeformationAfterReach = std::any_of(
                result.deformationIterations.begin(),
                result.deformationIterations.end(),
                [&](const optimization::AlternatingDeformationRecord& record)
                {
                    return record.block == block.block;
                }
            );
            if (
                hasDeformationAfterReach ||
                blockIndex + 1U != result.carvingBlocks.size() ||
                result.status !=
                    optimization::AlternatingOptimizationStatus::
                        StableAfterCarving
            )
            {
                reason = "controller did not stop at carving stability";
                return false;
            }
        }
    }
    return true;
}


int runRefillBlock(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeOverhangingPrism(vertices, faces);
    const ProductionFixture fixture(vertices, faces, 18);
    const optimization::AlternatingBalanceOptimizationResult result =
        runOverhang(fixture, 3, 1.0e-3);
    std::string traceReason;

    if (
        !validateCarvingTrace(result, traceReason) ||
        result.carvingBlocks.size() < 2U ||
        result.deformationIterations.empty()
    )
    {
        return fail(
            scenario,
            traceReason.empty()
                ? "fixture did not execute Carving-Deformation-Carving"
                : traceReason
        );
    }

    const optimization::AlternatingCarvingRecord& first =
        result.carvingBlocks[0];
    const optimization::AlternatingCarvingRecord& second =
        result.carvingBlocks[1];
    bool sawFirstFixedAlpha = false;
    for (const optimization::AlternatingDeformationRecord& record :
         result.deformationIterations)
    {
        if (record.block != 0)
        {
            continue;
        }
        sawFirstFixedAlpha = true;
        if (
            record.alphaHashBefore != first.alphaHash ||
            record.alphaHashAfter != first.alphaHash ||
            record.alphaHammingDistance != 0 ||
            record.refillCalls != 0 ||
            record.carvingCalls != 0
        )
        {
            return fail(
                scenario,
                "deformation block changed or rebuilt frozen alpha"
            );
        }
    }

    if (
        !sawFirstFixedAlpha ||
        first.refilledVoxelCount != 0 ||
        first.alphaHashBeforeReset != first.alphaHashAfterReset ||
        first.geometryEpochAfterReset !=
            first.geometryEpochBeforeReset ||
        first.carvedVoxelCount <= 0 ||
        second.refilledVoxelCount != first.carvedVoxelCount ||
        second.alphaHashAfterReset != first.alphaHashBeforeReset ||
        second.geometryEpochAfterReset !=
            second.geometryEpochBeforeReset + 1 ||
        second.carvedVoxelCount <= 0 ||
        second.alphaHammingDistance <= 0 ||
        second.alphaHash == first.alphaHash
    )
    {
        return fail(
            scenario,
            "next carving block did not refill A before producing pattern B"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " carved_A=" << first.carvedVoxelCount
        << " refill_B=" << second.refilledVoxelCount
        << " pattern_hamming=" << second.alphaHammingDistance
        << " deformation_trials=" << result.deformationIterations.size()
        << '\n';
    return EXIT_SUCCESS;
}


int runTransactionAndStability(const std::string& scenario)
{
    Eigen::MatrixXd boxVertices;
    Eigen::MatrixXi boxFaces;
    makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        boxVertices,
        boxFaces
    );
    const ProductionFixture stableFixture(boxVertices, boxFaces, 10);
    const stability::BalanceTargetSchedule stableSchedule =
        makeStandingSchedule(stableFixture.support, 0.08);
    const optimization::AlternatingBalanceOptimizer stableOptimizer(
        stableFixture.system,
        stableFixture.grid,
        stableSchedule,
        makeConfiguration(2, 1.0e-8)
    );
    const Eigen::VectorXd stableInitial =
        stableFixture.binding.handles.optimizationParameters();
    const optimization::AlternatingBalanceOptimizationResult stable =
        stableOptimizer.optimize(stableInitial, 0.2);
    std::string traceReason;

    if (
        !validateCarvingTrace(stable, traceReason) ||
        stable.status !=
            optimization::AlternatingOptimizationStatus::StableAfterCarving ||
        !stable.targets.allReached ||
        !stable.parameters.isApprox(stableInitial, 0.0) ||
        stable.carvingBlocks.size() != 1U ||
        !stable.deformationIterations.empty()
    )
    {
        return fail(
            scenario,
            traceReason.empty()
                ? "stable initial carving did not terminate before changing q"
                : traceReason
        );
    }

    Eigen::MatrixXd overhangVertices;
    Eigen::MatrixXi overhangFaces;
    makeOverhangingPrism(overhangVertices, overhangFaces);
    const ProductionFixture unsafeFixture(
        overhangVertices,
        overhangFaces,
        18
    );
    const optimization::AlternatingBalanceOptimizationResult unsafe =
        runOverhang(unsafeFixture, 1, 1.0e-7);
    traceReason.clear();

    if (
        !validateCarvingTrace(unsafe, traceReason) ||
        unsafe.targets.constraints.size() != 1U ||
        !unsafe.targets.constraints[0].standingState.has_value() ||
        *unsafe.targets.constraints[0].standingState !=
            stability::StandingState::FeasibleButUnsafe ||
        unsafe.targets.allReached ||
        unsafe.deformationIterations.empty()
    )
    {
        return fail(
            scenario,
            traceReason.empty()
                ? "FeasibleButUnsafe was treated as terminal stability"
                : traceReason
        );
    }

    bool sawAccepted = false;
    bool sawRejected = false;
    bool sawCommittedTargetRefresh = false;
    int expectedEpoch = unsafe.carvingBlocks.front().geometryEpoch;
    std::size_t expectedGeneration =
        unsafe.carvingBlocks.front().targetGenerationAfter;
    Eigen::VectorXd expectedParameters =
        unsafe.deformationIterations.front().parametersBefore;
    MassProperties expectedMass =
        unsafe.carvingBlocks.front().massProperties;
    Eigen::Vector3d expectedCenter =
        unsafe.carvingBlocks.front().centerOfMass;

    for (const optimization::AlternatingDeformationRecord& record :
         unsafe.deformationIterations)
    {
        if (
            record.refillCalls != 0 ||
            record.carvingCalls != 0 ||
            record.alphaHammingDistance != 0 ||
            record.alphaHashBefore != record.alphaHashAfter ||
            record.geometryEpochBefore != expectedEpoch ||
            record.targetGenerationBefore != expectedGeneration ||
            !std::isfinite(record.refreshedTargetEnergyAfter) ||
            !record.parametersBefore.isApprox(expectedParameters, 0.0)
        )
        {
            return fail(
                scenario,
                "fixed-alpha transaction pre-state was not preserved"
            );
        }

        if (record.accepted)
        {
            sawAccepted = true;
            const double maximumAcceptedEnergy =
                (1.0 - 0.01) * record.energyBefore;
            if (
                record.rolledBack ||
                record.decision !=
                    deformation::PaperTrialDecision::Accepted ||
                record.energyAfter >
                    maximumAcceptedEnergy +
                    1.0e-12 * std::max(1.0, record.energyBefore) ||
                record.geometryEpochAfter != expectedEpoch + 1 ||
                record.targetGenerationAfter < expectedGeneration
            )
            {
                return fail(
                    scenario,
                    "accepted deformation violated the 1% commit rule"
                );
            }
            if (record.targetGenerationAfter > expectedGeneration)
            {
                sawCommittedTargetRefresh = true;
            }
            expectedEpoch = record.geometryEpochAfter;
            expectedGeneration = record.targetGenerationAfter;
            expectedParameters = record.parametersAfter;
            expectedMass = record.massPropertiesAfter;
            expectedCenter = record.centerOfMassAfter;
            continue;
        }

        sawRejected = true;
        if (
            !record.rolledBack ||
            record.geometryEpochAfter != expectedEpoch ||
            record.targetGenerationAfter != expectedGeneration ||
            record.objectiveContextChanged ||
            !record.parametersAfter.isApprox(expectedParameters, 0.0) ||
            !sameMassProperties(
                record.massPropertiesAfter,
                expectedMass,
                NUMERICAL_TOLERANCE
            ) ||
            !record.centerOfMassAfter.isApprox(
                expectedCenter,
                NUMERICAL_TOLERANCE
            ) ||
            (
                record.decision ==
                    deformation::PaperTrialDecision::StepReduced &&
                record.energyAfter != record.energyBefore
            )
        )
        {
            return fail(
                scenario,
                "rejected deformation did not fully roll back geometry"
            );
        }
    }

    if (!sawAccepted || !sawRejected || !sawCommittedTargetRefresh)
    {
        return fail(
            scenario,
            "fixture did not exercise commit, rollback, and target refresh"
        );
    }

    if (!supportNodesFixed(
        unsafeFixture.system,
        unsafe.geometry.gridCorners,
        1.0e-12
    ))
    {
        return fail(scenario, "support nodes moved during deformation");
    }

    stability::BalanceTargetSpecifications multiSpecifications;
    multiSpecifications.emplace_back(
        stability::StandingTargetSpecification(
            stableFixture.support,
            0.08
        )
    );
    multiSpecifications.emplace_back(
        stability::SuspensionTargetSpecification(
            Eigen::Vector3d(0.8, 1.2, 0.5),
            Eigen::Vector3d(0.0, -1.0, 0.0),
            0.01
        )
    );
    stability::BalanceTargetSchedule multiSchedule(multiSpecifications);
    const MassProperties solidMass = computeSurfaceMassProperties(
        boxVertices,
        boxFaces,
        1.0
    );
    Eigen::Vector3d solidCenter = Eigen::Vector3d::Zero();
    if (!solidMass.tryCenterOfMass(solidCenter, 0.0))
    {
        return fail(scenario, "multi-target fixture lost positive mass");
    }
    const stability::BalanceTargetSnapshot initialMulti =
        multiSchedule.refresh(solidCenter);
    if (
        initialMulti.constraints.size() != 2U ||
        !initialMulti.constraints[0].reached ||
        initialMulti.constraints[1].reached ||
        initialMulti.allReached
    )
    {
        return fail(scenario, "multi-target fixture is not one-of-two stable");
    }

    optimization::AlternatingBalanceOptimizerConfig multiConfig =
        makeConfiguration(1, 0.0);
    multiConfig.maximumDeformationIterationsPerBlock = 1;
    const optimization::AlternatingBalanceOptimizer multiOptimizer(
        stableFixture.system,
        stableFixture.grid,
        multiSchedule,
        multiConfig
    );
    const optimization::AlternatingBalanceOptimizationResult multi =
        multiOptimizer.optimize(stableInitial, 0.2);
    traceReason.clear();
    if (
        !validateCarvingTrace(multi, traceReason) ||
        multi.status ==
            optimization::AlternatingOptimizationStatus::StableAfterCarving ||
        multi.carvingBlocks.empty() ||
        multi.carvingBlocks.front().allTargetsReached
    )
    {
        return fail(
            scenario,
            traceReason.empty()
                ? "controller used any-of instead of all-of target semantics"
                : traceReason
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " accepted=" << sawAccepted
        << " rejected=" << sawRejected
        << " terminal_state=FeasibleButUnsafe"
        << " multi_constraints=" << initialMulti.constraints.size()
        << '\n';
    return EXIT_SUCCESS;
}


int runEndToEnd(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeOverhangingPrism(vertices, faces);
    const ProductionFixture fixture(vertices, faces, 18);
    const stability::BalanceTargetSchedule schedule =
        makeStandingSchedule(fixture.support, 0.08);
    const optimization::AlternatingBalanceOptimizerConfig config =
        makeConfiguration(3, 1.0e-3);
    const optimization::AlternatingBalanceOptimizer optimizer(
        fixture.system,
        fixture.grid,
        schedule,
        config
    );
    const Eigen::VectorXd initial =
        fixture.binding.handles.optimizationParameters();
    const optimization::AlternatingBalanceOptimizationResult first =
        optimizer.optimize(initial, 0.2);
    const optimization::AlternatingBalanceOptimizationResult repeated =
        optimizer.optimize(initial, 0.2);
    std::string traceReason;

    if (
        !validateCarvingTrace(first, traceReason) ||
        first.carvingBlocks.size() < 2U ||
        first.deformationIterations.empty() ||
        first.carvingBlocks[0].carvedVoxelCount <= 0 ||
        first.carvingBlocks[1].refilledVoxelCount !=
            first.carvingBlocks[0].carvedVoxelCount ||
        first.carvingBlocks[1].alphaHash ==
            first.carvingBlocks[0].alphaHash ||
        first.carvingBlocks[1].alphaHammingDistance <= 0
    )
    {
        return fail(
            scenario,
            traceReason.empty()
                ? "end-to-end run did not complete distinct carving patterns"
                : traceReason
        );
    }

    const MassProperties oracle = completeSurfaceMassOracle(first.geometry);
    Eigen::Vector3d oracleCenter = Eigen::Vector3d::Zero();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> inertiaSolver(
        first.finalInertiaTensorAboutCenter
    );
    if (
        !oracle.tryCenterOfMass(oracleCenter, 0.0) ||
        !sameMassProperties(
            first.massProperties,
            oracle,
            1.0e-9
        ) ||
        !first.centerOfMass.isApprox(oracleCenter, 1.0e-9) ||
        !close(
            first.finalMassMoments.mass,
            first.massProperties.mass,
            1.0e-9
        ) ||
        !first.finalMassMoments.firstMoment.isApprox(
            first.massProperties.firstMoment,
            1.0e-9
        ) ||
        !first.finalInertiaTensorAboutCenter.allFinite() ||
        !first.finalInertiaTensorAboutCenter.isApprox(
            first.finalInertiaTensorAboutCenter.transpose(),
            1.0e-12
        ) ||
        inertiaSolver.info() != Eigen::Success ||
        !(inertiaSolver.eigenvalues().array() > 0.0).all() ||
        !supportNodesFixed(
            fixture.system,
            first.geometry.gridCorners,
            1.0e-12
        )
    )
    {
        return fail(
            scenario,
            "final geometry, mass oracle, CoM, or fixed support disagrees"
        );
    }

    if (!sameResult(first, repeated))
    {
        return fail(
            scenario,
            "identical production runs were not deterministic/idempotent"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " blocks=" << first.carvingBlocks.size()
        << " deformation_trials=" << first.deformationIterations.size()
        << " final_mass=" << first.massProperties.mass
        << " final_com=" << first.centerOfMass.transpose()
        << '\n';
    return EXIT_SUCCESS;
}


int runReleasedBetaOuterLoop(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeOverhangingPrism(vertices, faces);
    const ProductionFixture fixture(vertices, faces, 14);
    const stability::BalanceTargetSchedule schedule =
        makeStandingSchedule(fixture.support, 0.08);
    optimization::AlternatingBalanceOptimizerConfig config =
        makeConfiguration(1, 0.0);
    config.maximumDeformationIterationsPerBlock = 8;
    config.shapeScaleLambda = 20.0;
    config.lineSearch.requiredRelativeDecrease = 0.03;
    config.lineSearch.releasedBetaSemantics = true;
    config.freezeTargetsAfterInitialRefresh = true;
    config.releasedBetaOuterLoop = true;

    const optimization::AlternatingBalanceOptimizer optimizer(
        fixture.system,
        fixture.grid,
        schedule,
        config
    );
    const optimization::AlternatingBalanceOptimizationResult result =
        optimizer.optimize(
            fixture.binding.handles.optimizationParameters(),
            0.75
        );

    if (
        !result.releasedBetaOuterLoopEnabled ||
        result.deformationIterations.empty() ||
        result.carvingBlocks.size() != 1U ||
        result.carvingBlocks.front().carvingIterations != 1
    )
    {
        return fail(
            scenario,
            "released-beta outer loop did not begin with exactly one carve"
        );
    }

    bool sawAccepted = false;
    bool sawRollback = false;
    for (const optimization::AlternatingDeformationRecord& record :
         result.deformationIterations)
    {
        if (
            !record.releasedBetaOuterLoop ||
            record.trialBalanceCalls != 1 ||
            record.trialCarvingIterations != 1 ||
            record.refillCalls != record.carvingCalls ||
            record.refillCalls !=
                record.trialBalanceCalls + record.rollbackBalanceCalls ||
            record.objectiveReweighted ||
            record.targetGenerationAfter != record.targetGenerationBefore
        )
        {
            return fail(
                scenario,
                "a beta trial skipped recarving or refreshed its frozen target"
            );
        }

        if (record.accepted)
        {
            sawAccepted = true;
            if (
                record.rolledBack ||
                record.rollbackBalanceCalls != 0 ||
                record.refillCalls != 1 ||
                record.alphaHashAfter != record.trialAlphaHash
            )
            {
                return fail(
                    scenario,
                    "accepted beta trial did not retain its recarved state"
                );
            }
        }
        else
        {
            sawRollback = true;
            if (
                !record.rolledBack ||
                record.rollbackBalanceCalls != 1 ||
                record.rollbackCarvingIterations != 1 ||
                record.refillCalls != 2 ||
                record.alphaHashAfter != record.rollbackAlphaHash ||
                record.alphaHashAfter != record.alphaHashBefore ||
                record.geometryHashAfter != record.geometryHashBefore ||
                !record.parametersAfter.isApprox(
                    record.parametersBefore,
                    0.0
                )
            )
            {
                return fail(
                    scenario,
                    "uphill beta trial was not fully rebalanced after rollback"
                );
            }
        }
    }

    const MassProperties oracle = completeSurfaceMassOracle(result.geometry);
    if (
        !sawAccepted ||
        !sawRollback ||
        !sameMassProperties(result.massProperties, oracle, 1.0e-9)
    )
    {
        return fail(
            scenario,
            "fixture did not prove accepted and rejected beta balance calls"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " trials=" << result.deformationIterations.size()
        << " accepted=" << sawAccepted
        << " rollback=" << sawRollback
        << " final_status=" << static_cast<int>(result.status)
        << '\n';
    return EXIT_SUCCESS;
}


int runInvalidOuterBoundaryGate(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        vertices,
        faces
    );
    Eigen::MatrixXd crossingVertices;
    Eigen::MatrixXi crossingFaces;
    makeBox(
        Eigen::Vector3d(0.5, 0.5, 0.5),
        Eigen::Vector3d(1.5, 1.5, 1.5),
        crossingVertices,
        crossingFaces
    );
    Eigen::MatrixXd joinedVertices(vertices.rows() + crossingVertices.rows(), 3);
    joinedVertices << vertices, crossingVertices;
    Eigen::MatrixXi joinedFaces(faces.rows() + crossingFaces.rows(), 3);
    joinedFaces.topRows(faces.rows()) = faces;
    joinedFaces.bottomRows(crossingFaces.rows()) =
        crossingFaces.array() + vertices.rows();
    const VoxelGrid grid(joinedVertices, joinedFaces, 0.0, 8);
    const deformation::GridDeformationDomain domain =
        deformation::GridDeformationDomain::fromVoxelGrid(grid);
    const GeometryBinding binding = makeVerticalBinding(domain);
    const deformation::DeformationSystem invalidSystem(
        grid,
        joinedVertices,
        joinedFaces,
        binding.handles,
        binding.weights,
        true
    );
    const stability::SupportRegion support =
        stability::SupportRegion::fromMesh(
            joinedVertices,
            Eigen::Vector3d(0.0, -1.0, 0.0),
            0.0
        );
    const optimization::AlternatingBalanceOptimizer optimizer(
        invalidSystem,
        grid,
        makeStandingSchedule(support, 0.08),
        makeConfiguration(1, 1.0e-8)
    );
    bool rejected = false;
    try
    {
        static_cast<void>(optimizer.optimize(
            binding.handles.optimizationParameters(),
            0.2
        ));
    }
    catch (const deformation::InfeasibleDeformationState& error)
    {
        rejected = std::string(error.what()).find("outer boundary") !=
            std::string::npos;
    }
    if (!rejected)
    {
        return fail(
            scenario,
            "invalid outer boundary was not rejected by deformation evaluation"
        );
    }
    std::cout << "[PASS] scenario=" << scenario << " rejected=1\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail(
            "unknown",
            "expected refill-block, transaction, e2e, beta-outer, or invalid-outer"
        );
    }

    const std::string scenario = argv[1];
    if (scenario == "refill-block")
    {
        return runRefillBlock(scenario);
    }
    if (scenario == "transaction")
    {
        return runTransactionAndStability(scenario);
    }
    if (scenario == "e2e")
    {
        return runEndToEnd(scenario);
    }
    if (scenario == "beta-outer")
    {
        return runReleasedBetaOuterLoop(scenario);
    }
    if (scenario == "invalid-outer")
    {
        return runInvalidOuterBoundaryGate(scenario);
    }
    return fail(scenario, "unknown scenario");
}
