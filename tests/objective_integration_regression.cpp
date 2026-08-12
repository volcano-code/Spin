#include <cstdint>

#include "carving/CarvingSolver.h"
#include "grid_util.h"
#include "physics/BalanceObjective.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <igl/read_triangle_mesh.h>


namespace
{

bool approximatelyEqual(
    double actual,
    double expected
)
{
    const double scale = std::max(
        std::abs(actual),
        std::abs(expected)
    );

    return
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            1e-12 + 1e-10 * scale;
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=solver-boundary "
              << message << '\n';
    return EXIT_FAILURE;
}

} // namespace


int main(
    int argc,
    char** argv
)
{
    if (argc != 2)
    {
        return fail("usage: objective-integration-regression <mesh>");
    }

    Eigen::MatrixXd inputVertices;
    Eigen::MatrixXi faces;

    if (
        !igl::read_triangle_mesh(
            argv[1],
            inputVertices,
            faces
        ) ||
        inputVertices.rows() == 0 ||
        faces.rows() == 0
    )
    {
        return fail("could not load the closed integration fixture");
    }

    Eigen::MatrixXd vertices;
    alignToAxis(inputVertices, vertices);

    VoxelGrid voxelGrid;
    voxelGrid.build(vertices, faces);

    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -3.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d(0.0, 0.0, 1.0),
        Eigen::Vector3d(-2.0, 0.0, 0.0)
    );

    CarvingSolver solver(
        voxelGrid,
        vertices,
        faces,
        objectives
    );

    if (solver.minimumRemovalsToContinue() != 1)
    {
        return fail("the paper's unspecified progress threshold has the wrong default");
    }

    solver.setMinimumRemovalsToContinue(1);

    bool rejectedMinimum = false;

    try
    {
        solver.setMinimumRemovalsToContinue(0);
    }
    catch (const std::invalid_argument&)
    {
        rejectedMinimum = true;
    }

    if (!rejectedMinimum)
    {
        return fail("the solver accepted a non-positive progress threshold");
    }

    const double expectedDistance =
        aggregateProjectedCenterOfMassDistance(
            solver.currentCenterOfMass(),
            objectives
        );

    const double expectedEnergy =
        aggregateProjectedCenterOfMassEnergy(
            solver.currentCenterOfMass(),
            objectives
        );

    if (
        solver.balanceObjectives().size() != 3 ||
        !approximatelyEqual(
            solver.getCoMDistance(),
            expectedDistance
        ) ||
        !approximatelyEqual(
            solver.getCoMEnergy(),
            expectedEnergy
        )
    )
    {
        return fail("the production solver did not consume all objectives");
    }

    const int iterationsBefore = solver.iterationCount();
    const int carvedBefore = solver.carvedVoxelCount();
    const double massBefore = solver.objectMass();
    const Eigen::Vector3d firstMomentBefore =
        solver.currentMassProperties().firstMoment;
    const std::size_t historyBefore =
        solver.diagnosticsHistory().size();

    if (
        solver.stopReason() != CarvingStopReason::None ||
        solver.isOptimized()
    )
    {
        return fail("multi-target solver started in a terminal state");
    }

    Eigen::MatrixXd planeVertices;
    Eigen::MatrixXi planeFaces;
    Eigen::Vector3d returnedCenterOfMass;

    solver.carveStep(
        planeVertices,
        planeFaces,
        returnedCenterOfMass
    );

    bool rejectedLateMinimum = false;

    try
    {
        solver.setMinimumRemovalsToContinue(2);
    }
    catch (const std::logic_error&)
    {
        rejectedLateMinimum = true;
    }

    if (!rejectedLateMinimum)
    {
        return fail("the solver changed its progress threshold after carving began");
    }

    const CarvingStepDiagnostics* diagnostics =
        solver.lastDiagnostics();

    const CarvingStopReason firstStepReason =
        diagnostics == nullptr ?
            CarvingStopReason::MaximumIterationsReached :
            diagnostics->stopReason;

    if (
        diagnostics == nullptr ||
        solver.iterationCount() != iterationsBefore + 1 ||
        solver.carvedVoxelCount() < carvedBefore ||
        solver.objectMass() > massBefore ||
        solver.diagnosticsHistory().size() != historyBefore + 1 ||
        !approximatelyEqual(
            diagnostics->projectedEnergyAfter,
            solver.getCoMEnergy()
        ) ||
        diagnostics->projectedEnergyAfter >
            diagnostics->projectedEnergyBefore + 1e-12 ||
        (returnedCenterOfMass - solver.currentCenterOfMass()).norm() > 1e-12 ||
        planeVertices.rows() != 12 ||
        planeFaces.rows() != 6
    )
    {
        return fail("multi-target production carving contract failed");
    }

    if (diagnostics->carvedVoxelCountThisStep <= 0)
    {
        return fail("multi-target fixture did not produce a positive first prefix");
    }

    VoxelGrid thresholdGrid;
    thresholdGrid.build(vertices, faces);
    CarvingSolver thresholdSolver(
        thresholdGrid,
        vertices,
        faces,
        objectives
    );
    thresholdSolver.setMinimumRemovalsToContinue(
        diagnostics->carvedVoxelCountThisStep + 1
    );

    Eigen::MatrixXd thresholdPlaneVertices;
    Eigen::MatrixXi thresholdPlaneFaces;
    Eigen::Vector3d thresholdCenterOfMass;
    thresholdSolver.carveStep(
        thresholdPlaneVertices,
        thresholdPlaneFaces,
        thresholdCenterOfMass
    );

    if (
        thresholdSolver.carvedVoxelCount() !=
            diagnostics->carvedVoxelCountThisStep ||
        thresholdSolver.stopReason() !=
            CarvingStopReason::TooFewVoxelsRemoved ||
        thresholdSolver.lastDiagnostics() == nullptr ||
        thresholdSolver.lastDiagnostics()->stopReason !=
            CarvingStopReason::TooFewVoxelsRemoved
    )
    {
        return fail("positive improving prefix was not committed before too-few termination");
    }

    VoxelGrid boundaryGrid;
    boundaryGrid.build(vertices, faces);
    CarvingSolver boundarySolver(
        boundaryGrid,
        vertices,
        faces,
        objectives
    );
    boundarySolver.setMinimumRemovalsToContinue(
        diagnostics->carvedVoxelCountThisStep
    );
    boundarySolver.carveStep(
        thresholdPlaneVertices,
        thresholdPlaneFaces,
        thresholdCenterOfMass
    );

    if (
        boundarySolver.carvedVoxelCount() !=
            diagnostics->carvedVoxelCountThisStep ||
        boundarySolver.stopReason() ==
            CarvingStopReason::TooFewVoxelsRemoved
    )
    {
        return fail("too-few boundary used <= instead of the documented strict comparison");
    }

    if (firstStepReason == CarvingStopReason::None)
    {
        solver.carveAll();
    }

    if (!solver.isOptimized())
    {
        return fail("multi-target carveAll did not terminate");
    }

    static_cast<void>(firstMomentBefore);

    bool rejectedEmpty = false;

    try
    {
        CarvingSolver emptySolver(
            voxelGrid,
            vertices,
            faces,
            BalanceObjectives()
        );
        static_cast<void>(emptySolver);
    }
    catch (const std::invalid_argument&)
    {
        rejectedEmpty = true;
    }
    catch (...)
    {
        return fail("empty objectives threw the wrong exception type");
    }

    if (!rejectedEmpty)
    {
        return fail("the production solver accepted an empty objective set");
    }

    std::cout
        << "[PASS] scenario=solver-boundary objectives=3 "
        << "aggregate=consumed carving=paper-r4\n";

    return EXIT_SUCCESS;
}
