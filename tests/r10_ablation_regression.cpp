#include "deformation/DeformationSystem.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "optimization/AlternatingBalanceOptimizer.h"
#include "stability/BalanceTargetSchedule.h"
#include "stability/SupportRegion.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& reason)
{
    std::cerr << "[FAIL] scenario=r10-ablations reason=" << reason << '\n';
    return EXIT_FAILURE;
}


bool close(double first, double second, double tolerance = 1.0e-12)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <= tolerance *
            std::max({1.0, std::abs(first), std::abs(second)});
}


void makeOverhangingPrism(
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        0.00, 0.0, 0.0,
        0.25, 0.0, 0.0,
        1.00, 1.0, 0.0,
        0.75, 1.0, 0.0,
        0.00, 0.0, 1.0,
        0.25, 0.0, 1.0,
        1.00, 1.0, 1.0,
        0.75, 1.0, 1.0;

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


deformation::Handle makeSupportHandle()
{
    deformation::Handle result;
    result.role = deformation::HandleRole::Support;
    result.locked = true;
    return result;
}


deformation::Handle makeUserHandle(const Eigen::Vector3d& pivot)
{
    deformation::Handle result;
    result.restPivot = pivot;
    result.role = deformation::HandleRole::User;
    result.locked = false;
    return result;
}


struct Fixture
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    VoxelGrid grid;
    deformation::GridDeformationDomain domain;
    deformation::HandleSet handles;
    Eigen::MatrixXd weights;
    deformation::DeformationSystem system;
    stability::SupportRegion support;

    explicit Fixture(int resolution)
        :
        vertices(),
        faces(),
        grid(makeGrid(vertices, faces, resolution)),
        domain(deformation::GridDeformationDomain::fromVoxelGrid(grid)),
        handles(makeHandles(domain)),
        weights(makeWeights(domain)),
        system(grid, vertices, faces, handles, weights, true),
        support(stability::SupportRegion::fromMesh(
            vertices,
            Eigen::Vector3d(0.0, -1.0, 0.0),
            0.0
        ))
    {
    }

private:
    static VoxelGrid makeGrid(
        Eigen::MatrixXd& vertices,
        Eigen::MatrixXi& faces,
        int resolution
    )
    {
        makeOverhangingPrism(vertices, faces);
        return VoxelGrid(vertices, faces, 0.0, resolution);
    }

    static deformation::HandleSet makeHandles(
        const deformation::GridDeformationDomain& domain
    )
    {
        const double minimum = domain.restNodes().col(1).minCoeff();
        const double maximum = domain.restNodes().col(1).maxCoeff();
        deformation::Handle supportHandle = makeSupportHandle();
        deformation::Handle userHandle = makeUserHandle(
            Eigen::Vector3d(0.875, maximum, 0.5)
        );

        for (Eigen::Index node = 0; node < domain.nodeCount(); ++node)
        {
            const double y = domain.restNodes()(node, 1);
            if (close(y, minimum))
            {
                supportHandle.constrainedNodes.push_back(
                    static_cast<int>(node)
                );
            }
            if (close(y, maximum))
            {
                userHandle.constrainedNodes.push_back(
                    static_cast<int>(node)
                );
            }
        }
        return deformation::HandleSet(
            deformation::HandleVector{supportHandle, userHandle}
        );
    }

    static Eigen::MatrixXd makeWeights(
        const deformation::GridDeformationDomain& domain
    )
    {
        const double minimum = domain.restNodes().col(1).minCoeff();
        const double maximum = domain.restNodes().col(1).maxCoeff();
        const double span = maximum - minimum;
        Eigen::MatrixXd result(domain.nodeCount(), 2);
        for (Eigen::Index node = 0; node < domain.nodeCount(); ++node)
        {
            const double user =
                (domain.restNodes()(node, 1) - minimum) / span;
            result(node, 0) = 1.0 - user;
            result(node, 1) = user;
        }
        return result;
    }
};


stability::BalanceTargetSchedule makeSchedule(
    const stability::SupportRegion& support
)
{
    stability::BalanceTargetSpecifications specifications;
    specifications.emplace_back(
        stability::StandingTargetSpecification(support, 0.08)
    );
    return stability::BalanceTargetSchedule(specifications);
}


optimization::AlternatingBalanceOptimizerConfig makeConfig()
{
    optimization::AlternatingBalanceOptimizerConfig config;
    config.maximumCarvingBlocks = 2;
    config.maximumDeformationIterationsPerBlock = 200;
    config.gradientTolerance = 1.0e-7;
    config.shapeScaleLambda = 0.01;
    config.lineSearch.initialStep = 1.0;
    config.lineSearch.minimumStep = 0.4;
    config.lineSearch.reductionRatio = 0.8;
    config.lineSearch.requiredRelativeDecrease = 0.01;
    config.lineSearch.muDecrement = 0.05;
    return config;
}


optimization::AlternatingBalanceOptimizationResult run(
    const Fixture& fixture,
    const optimization::AlternatingBalanceOptimizerConfig& config
)
{
    const optimization::AlternatingBalanceOptimizer optimizer(
        fixture.system,
        fixture.grid,
        makeSchedule(fixture.support),
        config
    );
    return optimizer.optimize(
        fixture.handles.optimizationParameters(),
        0.2
    );
}


std::vector<unsigned char> fillMask(const VoxelGrid& grid)
{
    std::vector<unsigned char> result;
    const Eigen::Vector3i dimensions = grid.dimensions();
    result.reserve(static_cast<std::size_t>(dimensions.prod()));
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                result.push_back(grid(x, y, z).isFilled ? 1U : 0U);
            }
        }
    }
    return result;
}


bool validateEventLedger(
    const optimization::AlternatingBalanceOptimizationResult& result,
    std::string& reason
)
{
    if (
        result.events.empty() ||
        result.events.front().kind !=
            optimization::AlternatingEventKind::Initial ||
        result.events.front().ordinal != 0 ||
        result.events.front().block != -1 ||
        result.events.front().iteration != -1 ||
        result.events.front().geometryChanged
    )
    {
        reason = "event ledger does not begin with the immutable input state";
        return false;
    }

    int earliestFeasible = -1;
    int earliestSafe = -1;
    int previousEpoch = -1;
    for (std::size_t index = 0; index < result.events.size(); ++index)
    {
        const optimization::AlternatingEventRecord& event =
            result.events[index];
        if (
            event.ordinal != static_cast<int>(index) ||
            event.geometryEpoch < previousEpoch ||
            !event.massProperties.isFinite() ||
            !event.centerOfMass.allFinite() ||
            !std::isfinite(event.projectedTargetEnergy)
        )
        {
            reason = "event ordinals or committed-state payload are invalid";
            return false;
        }
        if (event.geometryChanged && event.geometryEpoch <= previousEpoch)
        {
            reason = "a geometry-changing event did not advance its epoch";
            return false;
        }
        if (
            index > 0U &&
            event.kind ==
                optimization::AlternatingEventKind::CarvingIteration &&
            event.geometryChanged &&
            event.geometryHash == result.events[index - 1U].geometryHash
        )
        {
            reason = "a committed carve retained the preceding geometry hash";
            return false;
        }
        previousEpoch = event.geometryEpoch;
        if (earliestFeasible < 0 && event.allStandingFeasible)
        {
            earliestFeasible = event.ordinal;
        }
        if (earliestSafe < 0 && event.allTargetsReached)
        {
            earliestSafe = event.ordinal;
        }
    }

    if (
        result.firstStandingFeasibleEvent != earliestFeasible ||
        result.firstSafeEvent != earliestSafe
    )
    {
        reason = "first feasible/safe event indices disagree with the ledger";
        return false;
    }
    const optimization::AlternatingEventRecord& terminal =
        result.events.back();
    if (
        terminal.geometryHash != result.geometryHash ||
        terminal.geometryEpoch != result.geometryEpoch
    )
    {
        reason = "terminal event does not bind the returned geometry state";
        return false;
    }
    return true;
}


bool containsEvent(
    const optimization::AlternatingBalanceOptimizationResult& result,
    optimization::AlternatingEventKind kind
)
{
    return std::any_of(
        result.events.begin(),
        result.events.end(),
        [kind](const optimization::AlternatingEventRecord& event)
        {
            return event.kind == kind;
        }
    );
}


bool equivalentZeroLaplacianResults(
    const optimization::AlternatingBalanceOptimizationResult& first,
    const optimization::AlternatingBalanceOptimizationResult& second
)
{
    if (
        first.status != second.status ||
        first.geometryEpoch != second.geometryEpoch ||
        first.geometryHash != second.geometryHash ||
        first.firstStandingFeasibleEvent !=
            second.firstStandingFeasibleEvent ||
        first.firstSafeEvent != second.firstSafeEvent ||
        first.events.size() != second.events.size() ||
        first.carvingBlocks.size() != second.carvingBlocks.size() ||
        first.carvingIterations.size() != second.carvingIterations.size() ||
        first.deformationIterations.size() !=
            second.deformationIterations.size() ||
        !first.parameters.isApprox(second.parameters, 0.0) ||
        fillMask(first.voxelFill) != fillMask(second.voxelFill) ||
        !close(first.finalProjectedTargetEnergy,
               second.finalProjectedTargetEnergy, 0.0) ||
        !close(first.finalUnscaledShapeEnergy,
               second.finalUnscaledShapeEnergy, 0.0) ||
        !close(first.finalShapeEnergy, second.finalShapeEnergy, 0.0) ||
        !close(first.finalTotalEnergy, second.finalTotalEnergy, 0.0)
    )
    {
        return false;
    }

    for (std::size_t index = 0; index < first.events.size(); ++index)
    {
        const auto& a = first.events[index];
        const auto& b = second.events[index];
        if (
            a.kind != b.kind ||
            a.ordinal != b.ordinal ||
            a.geometryEpoch != b.geometryEpoch ||
            a.geometryChanged != b.geometryChanged ||
            a.allStandingFeasible != b.allStandingFeasible ||
            a.allTargetsReached != b.allTargetsReached ||
            !close(a.projectedTargetEnergy, b.projectedTargetEnergy, 0.0) ||
            !a.centerOfMass.isApprox(b.centerOfMass, 0.0)
        )
        {
            return false;
        }
    }
    return true;
}


int runAblations()
{
    const Fixture fixture(18);
    const optimization::AlternatingBalanceOptimizerConfig base = makeConfig();

    const optimization::AlternatingBalanceOptimizationResult full =
        run(fixture, base);
    std::string ledgerReason;
    if (!validateEventLedger(full, ledgerReason))
    {
        return fail("full: " + ledgerReason);
    }
    if (
        full.carvingBlocks.empty() ||
        full.carvingIterations.empty() ||
        !containsEvent(
            full,
            optimization::AlternatingEventKind::CarvingIteration
        ) ||
        !close(full.effectiveShapeScaleLambda, base.shapeScaleLambda)
    )
    {
        return fail("full profile did not exercise carving and Laplacian");
    }

    optimization::AlternatingBalanceOptimizerConfig carvingOnly = base;
    carvingOnly.enableDeformation = false;
    const optimization::AlternatingBalanceOptimizationResult carved =
        run(fixture, carvingOnly);
    ledgerReason.clear();
    if (
        !validateEventLedger(carved, ledgerReason) ||
        carved.status !=
            optimization::AlternatingOptimizationStatus::
                DeformationDisabledUnstable ||
        carved.carvingIterations.empty() ||
        !carved.deformationIterations.empty() ||
        containsEvent(
            carved,
            optimization::AlternatingEventKind::AcceptedDeformation
        )
    )
    {
        return fail(
            ledgerReason.empty()
                ? "carving-only profile changed q or lacked its terminal status"
                : "carving-only: " + ledgerReason
        );
    }

    optimization::AlternatingBalanceOptimizerConfig deformationOnly = base;
    deformationOnly.enableCarving = false;
    const optimization::AlternatingBalanceOptimizationResult deformed =
        run(fixture, deformationOnly);
    ledgerReason.clear();
    if (
        !validateEventLedger(deformed, ledgerReason) ||
        !deformed.carvingBlocks.empty() ||
        !deformed.carvingIterations.empty() ||
        containsEvent(deformed, optimization::AlternatingEventKind::SolidRefill) ||
        containsEvent(
            deformed,
            optimization::AlternatingEventKind::CarvingIteration
        ) ||
        deformed.deformationIterations.empty() ||
        fillMask(deformed.voxelFill) != fillMask(fixture.grid)
    )
    {
        return fail(
            ledgerReason.empty()
                ? "deformation-only profile changed alpha or skipped q"
                : "deformation-only: " + ledgerReason
        );
    }

    optimization::AlternatingBalanceOptimizerConfig noLaplacian = base;
    noLaplacian.enableLaplacian = false;
    noLaplacian.shapeScaleLambda = 123.0;
    const optimization::AlternatingBalanceOptimizationResult noLap =
        run(fixture, noLaplacian);
    ledgerReason.clear();
    if (
        !validateEventLedger(noLap, ledgerReason) ||
        noLap.effectiveShapeScaleLambda != 0.0 ||
        noLap.finalShapeEnergy != 0.0 ||
        !std::isfinite(noLap.finalUnscaledShapeEnergy)
    )
    {
        return fail(
            ledgerReason.empty()
                ? "no-Laplacian profile retained a scaled shape penalty"
                : "no-Laplacian: " + ledgerReason
        );
    }

    optimization::AlternatingBalanceOptimizerConfig noCarvingNoLap =
        deformationOnly;
    noCarvingNoLap.enableLaplacian = false;
    const optimization::AlternatingBalanceOptimizationResult deformNoLap =
        run(fixture, noCarvingNoLap);
    if (
        !deformNoLap.carvingBlocks.empty() ||
        !deformNoLap.carvingIterations.empty() ||
        deformNoLap.effectiveShapeScaleLambda != 0.0 ||
        deformNoLap.finalShapeEnergy != 0.0
    )
    {
        return fail("no-carving/no-Laplacian did not compose both switches");
    }

    optimization::AlternatingBalanceOptimizerConfig lambdaZero = base;
    lambdaZero.shapeScaleLambda = 0.0;
    const optimization::AlternatingBalanceOptimizationResult zero =
        run(fixture, lambdaZero);
    optimization::AlternatingBalanceOptimizerConfig disabledLaplacian = base;
    disabledLaplacian.enableLaplacian = false;
    disabledLaplacian.shapeScaleLambda = 0.01;
    const optimization::AlternatingBalanceOptimizationResult switchedOff =
        run(fixture, disabledLaplacian);
    if (!equivalentZeroLaplacianResults(zero, switchedOff))
    {
        return fail("lambda=0 and disabled-Laplacian semantics diverged");
    }

    optimization::AlternatingBalanceOptimizerConfig disabled = base;
    disabled.enableCarving = false;
    disabled.enableDeformation = false;
    bool rejected = false;
    try
    {
        static_cast<void>(optimization::AlternatingBalanceOptimizer(
            fixture.system,
            fixture.grid,
            makeSchedule(fixture.support),
            disabled
        ));
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    if (!rejected)
    {
        return fail("configuration accepted both algorithms disabled");
    }

    std::cout
        << "[PASS] scenario=r10-ablations"
        << " full_events=" << full.events.size()
        << " carving_only_events=" << carved.events.size()
        << " deformation_only_events=" << deformed.events.size()
        << " first_feasible=" << full.firstStandingFeasibleEvent
        << " first_safe=" << full.firstSafeEvent
        << '\n';
    return EXIT_SUCCESS;
}


} // namespace


int main()
{
    return runAblations();
}
