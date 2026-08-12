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

    bool rejectedStopQuery = false;

    try
    {
        static_cast<void>(solver.stopReason());
    }
    catch (const std::logic_error&)
    {
        rejectedStopQuery = true;
    }
    catch (...)
    {
        return fail("multi-target stopReason threw the wrong exception type");
    }

    bool rejectedOptimizedQuery = false;

    try
    {
        static_cast<void>(solver.isOptimized());
    }
    catch (const std::logic_error&)
    {
        rejectedOptimizedQuery = true;
    }
    catch (...)
    {
        return fail("multi-target isOptimized threw the wrong exception type");
    }

    if (!rejectedStopQuery || !rejectedOptimizedQuery)
    {
        return fail("multi-target solver reported incomplete carving as optimized");
    }

    Eigen::MatrixXd planeVertices;
    Eigen::MatrixXi planeFaces;
    Eigen::Vector3d returnedCenterOfMass;

    bool rejectedStep = false;

    try
    {
        solver.carveStep(
            planeVertices,
            planeFaces,
            returnedCenterOfMass
        );
    }
    catch (const std::logic_error&)
    {
        rejectedStep = true;
    }
    catch (...)
    {
        return fail("multi-target carveStep threw the wrong exception type");
    }

    if (!rejectedStep)
    {
        return fail("multi-target carveStep silently used one target");
    }

    bool rejectedAll = false;

    try
    {
        solver.carveAll();
    }
    catch (const std::logic_error&)
    {
        rejectedAll = true;
    }
    catch (...)
    {
        return fail("multi-target carveAll threw the wrong exception type");
    }

    if (
        !rejectedAll ||
        solver.iterationCount() != iterationsBefore ||
        solver.carvedVoxelCount() != carvedBefore ||
        solver.objectMass() != massBefore ||
        solver.currentMassProperties().firstMoment != firstMomentBefore ||
        solver.diagnosticsHistory().size() != historyBefore
    )
    {
        return fail("the R4 boundary mutated production carving state");
    }

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
        << "aggregate=consumed carving=deferred-to-R4\n";

    return EXIT_SUCCESS;
}
