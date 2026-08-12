#include "carving/CarvingSolver.h"
#include "deformation/DeformedVoxelMassField.h"
#include "deformation/GridDeformationDomain.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


struct CarvingSolverTestAccess
{
    static void setPreCommitHook(
        CarvingSolver& solver,
        const std::function<void()>& hook
    )
    {
        solver.setPreCommitHookForTesting(hook);
    }
};


namespace
{


using deformation::DeformedVoxelMassField;
using deformation::GridDeformationDomain;


constexpr double TOLERANCE = 2.0e-10;


int fail(
    const std::string& scenario,
    const std::string& reason
)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


bool approximatelyEqual(
    double first,
    double second,
    double tolerance = TOLERANCE
)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <= tolerance *
            std::max({1.0, std::abs(first), std::abs(second)});
}


bool approximatelyEqual(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& second,
    double tolerance = TOLERANCE
)
{
    return
        first.allFinite() &&
        second.allFinite() &&
        (first - second).norm() <= tolerance *
            std::max({1.0, first.norm(), second.norm()});
}


bool approximatelyEqual(
    const MassProperties& first,
    const MassProperties& second,
    double tolerance = TOLERANCE
)
{
    return
        approximatelyEqual(first.mass, second.mass, tolerance) &&
        approximatelyEqual(
            first.firstMoment,
            second.firstMoment,
            tolerance
        );
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


Eigen::Vector3d warpPoint(const Eigen::Vector3d& point)
{
    Eigen::Vector3d result;
    result(0) =
        1.18 * point(0) +
        0.10 * point(1) +
        0.025 * point(2) * point(2) +
        0.31;
    result(1) =
        0.92 * point(1) +
        0.06 * point(2) +
        0.015 * point(0) * point(0) -
        0.17;
    result(2) =
        1.07 * point(2) +
        0.04 * point(0) +
        0.012 * point(1) * point(1) +
        0.23;
    return result;
}


Eigen::MatrixXd warpPoints(const Eigen::MatrixXd& points)
{
    Eigen::MatrixXd result(points.rows(), 3);
    for (Eigen::Index row = 0; row < points.rows(); ++row)
    {
        result.row(row) =
            warpPoint(points.row(row).transpose()).transpose();
    }
    return result;
}


struct Fixture
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    VoxelGrid grid;
    GridDeformationDomain domain;

    Fixture()
        : outerVertices(), outerFaces(), grid(), domain()
    {
        makeBox(
            Eigen::Vector3d(-1.0, -0.8, -0.65),
            Eigen::Vector3d(1.0, 0.8, 0.65),
            outerVertices,
            outerFaces
        );
        grid.build(outerVertices, outerFaces, 0.0, 8);
        domain = GridDeformationDomain::fromVoxelGrid(grid);
    }
};


int linearCellIndex(
    const Eigen::Vector3i& dimensions,
    int x,
    int y,
    int z
)
{
    return (x * dimensions(1) + y) * dimensions(2) + z;
}


std::vector<unsigned char> fillSnapshot(const VoxelGrid& grid)
{
    const Eigen::Vector3i dimensions = grid.dimensions();
    std::vector<unsigned char> result(
        static_cast<std::size_t>(dimensions.prod()),
        0U
    );

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                result[static_cast<std::size_t>(
                    linearCellIndex(dimensions, x, y, z)
                )] = grid(x, y, z).isFilled ? 1U : 0U;
            }
        }
    }
    return result;
}


MassProperties recomputeFromField(
    const VoxelGrid& grid,
    const MassProperties& outer,
    const DeformedVoxelMassField& field
)
{
    MassProperties result = outer;
    const Eigen::Vector3i dimensions = grid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel = grid(x, y, z);
                if (voxel.isInner && !voxel.isFilled)
                {
                    result -= field.massProperties(
                        linearCellIndex(dimensions, x, y, z)
                    );
                }
            }
        }
    }
    return result;
}


BalanceObjectives makeObjectives(const Eigen::Vector3d& fullCenter)
{
    BalanceObjectives objectives;
    objectives.emplace_back(
        fullCenter - Eigen::Vector3d(0.42, 0.0, 0.0),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );
    return objectives;
}


bool equalCommittedResult(
    const CarvingSolver& first,
    const CarvingSolver& second,
    const VoxelGrid& firstGrid,
    const VoxelGrid& secondGrid
)
{
    const CarvingStepDiagnostics* firstDiagnostics =
        first.lastDiagnostics();
    const CarvingStepDiagnostics* secondDiagnostics =
        second.lastDiagnostics();

    return
        firstDiagnostics != nullptr &&
        secondDiagnostics != nullptr &&
        fillSnapshot(firstGrid) == fillSnapshot(secondGrid) &&
        approximatelyEqual(
            first.currentMassProperties(),
            second.currentMassProperties(),
            3.0e-11
        ) &&
        approximatelyEqual(
            first.currentCenterOfMass(),
            second.currentCenterOfMass(),
            3.0e-11
        ) &&
        first.carvedVoxelCount() == second.carvedVoxelCount() &&
        firstDiagnostics->candidateVoxelCount ==
            secondDiagnostics->candidateVoxelCount &&
        firstDiagnostics->bestPrefixLength ==
            secondDiagnostics->bestPrefixLength &&
        approximatelyEqual(
            firstDiagnostics->bestPrefixEnergy,
            secondDiagnostics->bestPrefixEnergy,
            3.0e-11
        );
}


int runIdentity()
{
    const std::string scenario = "identity";
    Fixture fixture;
    VoxelGrid legacyGrid = fixture.grid;
    VoxelGrid sidecarGrid = fixture.grid;

    const MassProperties outer = computeSurfaceMassProperties(
        fixture.outerVertices,
        fixture.outerFaces
    );
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    if (!outer.tryCenterOfMass(center))
    {
        return fail(scenario, "outer center of mass is invalid");
    }

    const DeformedVoxelMassField identityField(
        fixture.domain,
        fixture.domain.fullRestCorners()
    );
    const BalanceObjectives objectives = makeObjectives(center);

    CarvingSolver legacy(
        legacyGrid,
        fixture.outerVertices,
        fixture.outerFaces,
        objectives
    );
    CarvingSolver sidecar(
        sidecarGrid,
        outer,
        identityField,
        objectives
    );

    Eigen::MatrixXd legacyPlane;
    Eigen::MatrixXi legacyFaces;
    Eigen::Vector3d legacyCenter;
    legacy.carveStep(legacyPlane, legacyFaces, legacyCenter);

    Eigen::MatrixXd sidecarPlane;
    Eigen::MatrixXi sidecarFaces;
    Eigen::Vector3d sidecarCenter;
    sidecar.carveStep(sidecarPlane, sidecarFaces, sidecarCenter);

    if (
        sidecar.carvedVoxelCount() <= 0 ||
        !equalCommittedResult(
            legacy,
            sidecar,
            legacyGrid,
            sidecarGrid
        ) ||
        !legacyPlane.isApprox(sidecarPlane, 3.0e-11) ||
        legacyFaces.rows() != sidecarFaces.rows() ||
        legacyFaces.cols() != sidecarFaces.cols() ||
        (legacyFaces.array() != sidecarFaces.array()).any() ||
        !approximatelyEqual(legacyCenter, sidecarCenter, 3.0e-11)
    )
    {
        return fail(
            scenario,
            "identity sidecar changed the legacy axis-aligned result"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " carved=" << sidecar.carvedVoxelCount() << '\n';
    return EXIT_SUCCESS;
}


struct DeformedFixture
{
    Fixture rest;
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXd corners;
    MassProperties outerMass;
    DeformedVoxelMassField field;
    BalanceObjectives objectives;

    DeformedFixture()
        :
        rest(),
        outerVertices(warpPoints(rest.outerVertices)),
        corners(warpPoints(rest.domain.fullRestCorners())),
        outerMass(computeSurfaceMassProperties(
            outerVertices,
            rest.outerFaces
        )),
        field(rest.domain, corners),
        objectives()
    {
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        if (!outerMass.tryCenterOfMass(center))
        {
            throw std::runtime_error(
                "deformed outer fixture has no positive center of mass"
            );
        }
        objectives = makeObjectives(center);
    }
};


int runDeformed()
{
    const std::string scenario = "deformed";
    DeformedFixture fixture;
    VoxelGrid grid = fixture.rest.grid;
    CarvingSolver solver(
        grid,
        fixture.outerMass,
        fixture.field,
        fixture.objectives
    );

    Eigen::MatrixXd plane;
    Eigen::MatrixXi planeFaces;
    Eigen::Vector3d center;
    solver.carveStep(plane, planeFaces, center);

    const MassProperties independent = recomputeFromField(
        grid,
        fixture.outerMass,
        fixture.field
    );
    Eigen::Vector3d independentCenter = Eigen::Vector3d::Zero();
    const CarvingStepDiagnostics* diagnostics = solver.lastDiagnostics();

    if (
        solver.carvedVoxelCount() <= 0 ||
        diagnostics == nullptr ||
        diagnostics->bestPrefixLength <= 0 ||
        !independent.tryCenterOfMass(independentCenter) ||
        !approximatelyEqual(
            solver.currentMassProperties(),
            independent
        ) ||
        !approximatelyEqual(
            solver.currentCenterOfMass(),
            independentCenter
        ) ||
        !approximatelyEqual(center, independentCenter) ||
        diagnostics->massDiscrepancyAfter > TOLERANCE ||
        diagnostics->firstMomentDiscrepancyAfter > TOLERANCE ||
        diagnostics->centerOfMassDiscrepancyAfter > TOLERANCE
    )
    {
        return fail(
            scenario,
            "candidate selection or full recompute did not use deformed exact mass"
        );
    }

    VoxelGrid wrongGrid;
    Eigen::MatrixXd wrongVertices;
    Eigen::MatrixXi wrongFaces;
    makeBox(
        Eigen::Vector3d(-1.0, -0.8, -0.65),
        Eigen::Vector3d(1.0, 0.8, 0.65),
        wrongVertices,
        wrongFaces
    );
    wrongGrid.build(wrongVertices, wrongFaces, 0.0, 7);

    bool sizeRejected = false;
    try
    {
        CarvingSolver invalid(
            wrongGrid,
            fixture.outerMass,
            fixture.field,
            fixture.objectives
        );
        static_cast<void>(invalid);
    }
    catch (const std::invalid_argument&)
    {
        sizeRejected = true;
    }

    if (!sizeRejected)
    {
        return fail(scenario, "a mismatched mass-field cell count was accepted");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " carved=" << solver.carvedVoxelCount()
        << " mass=" << solver.objectMass() << '\n';
    return EXIT_SUCCESS;
}


int runTransaction()
{
    const std::string scenario = "transaction";
    DeformedFixture fixture;
    VoxelGrid grid = fixture.rest.grid;
    CarvingSolver solver(
        grid,
        fixture.outerMass,
        fixture.field,
        fixture.objectives
    );

    const std::vector<unsigned char> fillsBefore = fillSnapshot(grid);
    const MassProperties massBefore = solver.currentMassProperties();
    const Eigen::Vector3d centerBefore = solver.currentCenterOfMass();
    bool observedStagedRemoval = false;

    CarvingSolverTestAccess::setPreCommitHook(
        solver,
        [&]()
        {
            observedStagedRemoval = fillSnapshot(grid) != fillsBefore;
            throw std::runtime_error("intentional deformed-carving failure");
        }
    );

    Eigen::MatrixXd plane = Eigen::MatrixXd::Constant(1, 3, 17.0);
    Eigen::MatrixXi planeFaces = Eigen::MatrixXi::Constant(1, 3, 19);
    Eigen::Vector3d center = Eigen::Vector3d::Constant(23.0);
    bool threw = false;

    try
    {
        solver.carveStep(plane, planeFaces, center);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    if (
        !threw ||
        !observedStagedRemoval ||
        fillSnapshot(grid) != fillsBefore ||
        solver.iterationCount() != 0 ||
        solver.carvedVoxelCount() != 0 ||
        solver.lastDiagnostics() != nullptr ||
        !approximatelyEqual(solver.currentMassProperties(), massBefore, 0.0) ||
        !approximatelyEqual(solver.currentCenterOfMass(), centerBefore, 0.0) ||
        plane.rows() != 1 || plane.cols() != 3 ||
        !plane.isConstant(17.0) ||
        planeFaces.rows() != 1 || planeFaces.cols() != 3 ||
        !planeFaces.isConstant(19) ||
        (center.array() != 23.0).any()
    )
    {
        return fail(
            scenario,
            "an exception did not roll back fill and all committed state"
        );
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runObjectiveUpdate()
{
    const std::string scenario = "objective-update";
    Fixture fixture;
    VoxelGrid refreshedGrid = fixture.grid;
    VoxelGrid staleGrid = fixture.grid;

    const MassProperties outer = computeSurfaceMassProperties(
        fixture.outerVertices,
        fixture.outerFaces
    );
    Eigen::Vector3d outerCenter = Eigen::Vector3d::Zero();
    if (!outer.tryCenterOfMass(outerCenter))
    {
        return fail(scenario, "outer center of mass is invalid");
    }

    const DeformedVoxelMassField field(
        fixture.domain,
        fixture.domain.fullRestCorners()
    );
    const BalanceObjectives initialObjectives =
        makeObjectives(outerCenter);

    CarvingSolver refreshed(
        refreshedGrid,
        outer,
        field,
        initialObjectives
    );
    CarvingSolver stale(
        staleGrid,
        outer,
        field,
        initialObjectives
    );

    Eigen::MatrixXd temporaryPlane;
    Eigen::MatrixXi temporaryFaces;
    Eigen::Vector3d temporaryCenter;
    refreshed.carveStep(
        temporaryPlane,
        temporaryFaces,
        temporaryCenter
    );
    stale.carveStep(
        temporaryPlane,
        temporaryFaces,
        temporaryCenter
    );

    if (
        refreshed.isOptimized() ||
        stale.isOptimized() ||
        !equalCommittedResult(
            refreshed,
            stale,
            refreshedGrid,
            staleGrid
        )
    )
    {
        return fail(
            scenario,
            "fixture did not provide a committed non-terminal first step"
        );
    }

    const std::vector<unsigned char> fillsBefore =
        fillSnapshot(refreshedGrid);
    const MassProperties massBefore =
        refreshed.currentMassProperties();
    const Eigen::Vector3d centerBefore =
        refreshed.currentCenterOfMass();
    const int iterationBefore = refreshed.iterationCount();
    const int carvedBefore = refreshed.carvedVoxelCount();
    const CarvingStopReason terminalBefore = refreshed.stopReason();
    const std::size_t historySizeBefore =
        refreshed.diagnosticsHistory().size();
    const CarvingStepDiagnostics historicalRecord =
        *refreshed.lastDiagnostics();

    BalanceObjectives refreshedObjectives;
    refreshedObjectives.emplace_back(
        centerBefore + Eigen::Vector3d(0.42, 0.0, 0.0),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );
    const double expectedRefreshedEnergy =
        aggregateProjectedCenterOfMassEnergy(
            centerBefore,
            refreshedObjectives
        );
    const double staleEnergy = refreshed.getCoMEnergy();

    refreshed.updateBalanceObjectives(refreshedObjectives);

    const CarvingStepDiagnostics* historyAfterUpdate =
        refreshed.lastDiagnostics();
    if (
        fillSnapshot(refreshedGrid) != fillsBefore ||
        !approximatelyEqual(
            refreshed.currentMassProperties(),
            massBefore,
            0.0
        ) ||
        !approximatelyEqual(
            refreshed.currentCenterOfMass(),
            centerBefore,
            0.0
        ) ||
        refreshed.iterationCount() != iterationBefore ||
        refreshed.carvedVoxelCount() != carvedBefore ||
        refreshed.stopReason() != terminalBefore ||
        refreshed.diagnosticsHistory().size() != historySizeBefore ||
        historyAfterUpdate == nullptr ||
        historyAfterUpdate->iteration != historicalRecord.iteration ||
        historyAfterUpdate->projectedEnergyAfter !=
            historicalRecord.projectedEnergyAfter ||
        refreshed.balanceObjectives().size() != 1 ||
        !refreshed.balanceObjectives().front().target().isApprox(
            refreshedObjectives.front().target(),
            0.0
        ) ||
        !approximatelyEqual(
            refreshed.getCoMEnergy(),
            expectedRefreshedEnergy,
            3.0e-11
        ) ||
        approximatelyEqual(staleEnergy, expectedRefreshedEnergy, 1.0e-6)
    )
    {
        return fail(
            scenario,
            "target refresh mutated committed state or rewrote history"
        );
    }

    const BalanceObjectives validBeforeEmpty =
        refreshed.balanceObjectives();
    bool emptyRejected = false;
    try
    {
        refreshed.updateBalanceObjectives(BalanceObjectives());
    }
    catch (const std::invalid_argument&)
    {
        emptyRejected = true;
    }

    if (
        !emptyRejected ||
        fillSnapshot(refreshedGrid) != fillsBefore ||
        !approximatelyEqual(
            refreshed.currentMassProperties(),
            massBefore,
            0.0
        ) ||
        !approximatelyEqual(
            refreshed.currentCenterOfMass(),
            centerBefore,
            0.0
        ) ||
        refreshed.iterationCount() != iterationBefore ||
        refreshed.carvedVoxelCount() != carvedBefore ||
        refreshed.stopReason() != terminalBefore ||
        refreshed.diagnosticsHistory().size() != historySizeBefore ||
        refreshed.balanceObjectives().size() != validBeforeEmpty.size() ||
        !refreshed.balanceObjectives().front().target().isApprox(
            validBeforeEmpty.front().target(),
            0.0
        )
    )
    {
        return fail(
            scenario,
            "empty target update violated the strong exception guarantee"
        );
    }

    Eigen::MatrixXd refreshedPlane;
    Eigen::MatrixXi refreshedFaces;
    Eigen::Vector3d refreshedCenter;
    refreshed.carveStep(
        refreshedPlane,
        refreshedFaces,
        refreshedCenter
    );

    Eigen::MatrixXd stalePlane;
    Eigen::MatrixXi staleFaces;
    Eigen::Vector3d staleCenter;
    stale.carveStep(stalePlane, staleFaces, staleCenter);

    const CarvingStepDiagnostics* refreshedStep =
        refreshed.lastDiagnostics();
    const CarvingStepDiagnostics* staleStep =
        stale.lastDiagnostics();
    if (refreshedStep == nullptr || staleStep == nullptr)
    {
        return fail(scenario, "second carving step has no diagnostics");
    }
    const bool resultsDiffer =
        fillSnapshot(refreshedGrid) != fillSnapshot(staleGrid) ||
        refreshedStep->candidateVoxelCount !=
            staleStep->candidateVoxelCount ||
        refreshedStep->bestPrefixLength !=
            staleStep->bestPrefixLength ||
        !approximatelyEqual(
            refreshed.currentCenterOfMass(),
            stale.currentCenterOfMass(),
            1.0e-8
        );

    if (
        refreshed.iterationCount() != iterationBefore + 1 ||
        !approximatelyEqual(
            refreshedStep->projectedEnergyBefore,
            expectedRefreshedEnergy,
            3.0e-11
        ) ||
        !resultsDiffer
    )
    {
        return fail(
            scenario,
            "the next committed step did not use the refreshed objectives"
        );
    }
    const int refreshedNextPrefix = refreshedStep->bestPrefixLength;

    VoxelGrid reentrantGrid = fixture.grid;
    CarvingSolver reentrant(
        reentrantGrid,
        outer,
        field,
        initialObjectives
    );
    const std::vector<unsigned char> reentrantFills =
        fillSnapshot(reentrantGrid);
    CarvingSolverTestAccess::setPreCommitHook(
        reentrant,
        [&]()
        {
            reentrant.updateBalanceObjectives(refreshedObjectives);
        }
    );
    bool reentrantRejected = false;
    try
    {
        reentrant.carveStep(
            temporaryPlane,
            temporaryFaces,
            temporaryCenter
        );
    }
    catch (const std::logic_error&)
    {
        reentrantRejected = true;
    }

    if (
        !reentrantRejected ||
        fillSnapshot(reentrantGrid) != reentrantFills ||
        reentrant.iterationCount() != 0 ||
        reentrant.carvedVoxelCount() != 0 ||
        reentrant.lastDiagnostics() != nullptr ||
        reentrant.balanceObjectives().size() !=
            initialObjectives.size() ||
        !reentrant.balanceObjectives().front().target().isApprox(
            initialObjectives.front().target(),
            0.0
        )
    )
    {
        return fail(
            scenario,
            "an in-step target update was not transactionally rejected"
        );
    }

    refreshed.carveAll();
    const CarvingStopReason finalTerminal = refreshed.stopReason();
    const int finalIterations = refreshed.iterationCount();
    const int finalCarved = refreshed.carvedVoxelCount();
    const MassProperties finalMass = refreshed.currentMassProperties();
    const std::vector<unsigned char> finalFills =
        fillSnapshot(refreshedGrid);
    const std::size_t finalHistorySize =
        refreshed.diagnosticsHistory().size();

    BalanceObjectives terminalObjectives;
    terminalObjectives.emplace_back(
        refreshed.currentCenterOfMass() -
            Eigen::Vector3d(0.17, 0.0, 0.0),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );
    refreshed.updateBalanceObjectives(terminalObjectives);

    if (
        finalTerminal == CarvingStopReason::None ||
        !refreshed.isOptimized() ||
        refreshed.stopReason() != finalTerminal ||
        refreshed.iterationCount() != finalIterations ||
        refreshed.carvedVoxelCount() != finalCarved ||
        !approximatelyEqual(
            refreshed.currentMassProperties(),
            finalMass,
            0.0
        ) ||
        fillSnapshot(refreshedGrid) != finalFills ||
        refreshed.diagnosticsHistory().size() != finalHistorySize
    )
    {
        return fail(
            scenario,
            "updating a terminal solver unexpectedly reopened or mutated it"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " old_energy=" << staleEnergy
        << " refreshed_energy=" << expectedRefreshedEnergy
        << " next_prefix=" << refreshedNextPrefix << '\n';
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: deformed-carving-solver-regression "
            << "<identity|deformed|transaction|objective-update>\n";
        return EXIT_FAILURE;
    }

    const std::string scenario(argv[1]);
    try
    {
        if (scenario == "identity")
        {
            return runIdentity();
        }
        if (scenario == "deformed")
        {
            return runDeformed();
        }
        if (scenario == "transaction")
        {
            return runTransaction();
        }
        if (scenario == "objective-update")
        {
            return runObjectiveUpdate();
        }
        return fail(scenario, "unknown scenario");
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
}
