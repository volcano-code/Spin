#include <cstdint>

#include "carving/CarvingSolver.h"
#include "grid_util.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include <igl/read_triangle_mesh.h>


namespace
{

bool isFiniteVector(
    const Eigen::Vector3d& value
)
{
    return
        std::isfinite(value(0)) &&
        std::isfinite(value(1)) &&
        std::isfinite(value(2));
}


int fail(
    const std::string& message
)
{
    std::cerr
        << "[FAIL] "
        << message
        << '\n';

    return EXIT_FAILURE;
}


Eigen::RowVector3d createGroundPlane(
    const Eigen::MatrixXd& vertices,
    Eigen::MatrixXd& planeVertices,
    Eigen::MatrixXi& planeFaces
)
{
    Eigen::AlignedBox3d boundingBox;
    createAlignedBox(
        vertices,
        boundingBox
    );

    planeVertices.resize(4, 3);

    planeVertices <<
        boundingBox.corner(
            boundingBox.BottomLeftCeil
        ).transpose(),
        boundingBox.corner(
            boundingBox.BottomRightCeil
        ).transpose(),
        boundingBox.corner(
            boundingBox.BottomLeftFloor
        ).transpose(),
        boundingBox.corner(
            boundingBox.BottomRightFloor
        ).transpose();

    planeFaces.resize(2, 3);

    planeFaces <<
        0, 1, 2,
        2, 1, 3;

    return planeVertices.colwise().sum() / 4.0;
}

} // namespace


int main(
    int argc,
    char* argv[]
)
{
    const std::string meshPath =
        argc > 1
            ? argv[1]
            : "../data/bunny.off";

    const std::string csvPath =
        argc > 2
            ? argv[2]
            : "carving_regression.csv";

    const int maximumSteps =
        argc > 3
            ? std::max(1, std::atoi(argv[3]))
            : 1;


    Eigen::MatrixXd inputVertices;
    Eigen::MatrixXi faces;

    if (
        !igl::read_triangle_mesh(
            meshPath,
            inputVertices,
            faces
        )
    )
    {
        return fail(
            "could not load mesh: " + meshPath
        );
    }

    if (
        inputVertices.rows() == 0 ||
        faces.rows() == 0
    )
    {
        return fail(
            "loaded mesh is empty"
        );
    }


    Eigen::MatrixXd vertices;

    alignToAxis(
        inputVertices,
        vertices
    );


    Eigen::MatrixXd planeVertices;
    Eigen::MatrixXi planeFaces;

    const Eigen::RowVector3d planeCenter =
        createGroundPlane(
            vertices,
            planeVertices,
            planeFaces
        );


    Eigen::AlignedBox3d boundingBox;

    createAlignedBox(
        vertices,
        boundingBox
    );

    /*
     * Use a deterministic target on the ground plane with a horizontal
     * offset. The offset avoids a degenerate cutting-plane triangle when the
     * original COM happens to lie exactly above the plane center.
     */
    Eigen::Vector3d targetCenterOfMass =
        planeCenter.transpose();

    targetCenterOfMass(0) +=
        0.15 * boundingBox.sizes()(0);


    VoxelGrid voxelGrid;

    voxelGrid.build(
        vertices,
        faces
    );


    CarvingSolver solver(
        voxelGrid,
        vertices,
        faces,
        planeVertices,
        planeFaces,
        Eigen::Vector3d(0.0, -1.0, 0.0),
        targetCenterOfMass
    );

    solver.setDiagnosticsEnabled(true);


    if (
        !std::isfinite(solver.objectMass()) ||
        std::abs(solver.objectMass()) <=
            std::numeric_limits<double>::epsilon()
    )
    {
        return fail(
            "initial mass is zero or non-finite"
        );
    }

    if (
        !isFiniteVector(
            solver.currentCenterOfMass()
        )
    )
    {
        return fail(
            "initial center of mass is non-finite"
        );
    }

    if (
        !std::isfinite(
            solver.getCoMDistance()
        ) ||
        !std::isfinite(
            solver.getGeometricProjectedDistance()
        )
    )
    {
        return fail(
            "initial projected error is non-finite"
        );
    }


    int previousCarvedVoxelCount =
        solver.carvedVoxelCount();

    int executedSteps = 0;

    while (
        !solver.isOptimized() &&
        executedSteps < maximumSteps
    )
    {
        Eigen::MatrixXd carvePlaneVertices;
        Eigen::MatrixXi carvePlaneFaces;
        Eigen::Vector3d centerOfMass;

        solver.carveStep(
            carvePlaneVertices,
            carvePlaneFaces,
            centerOfMass
        );

        ++executedSteps;

        const CarvingStepDiagnostics* diagnostics =
            solver.lastDiagnostics();

        if (diagnostics == nullptr)
        {
            return fail(
                "carveStep did not create a diagnostic record"
            );
        }

        if (
            diagnostics->iteration !=
            solver.iterationCount()
        )
        {
            return fail(
                "diagnostic iteration does not match solver state"
            );
        }

        if (
            diagnostics->candidateVoxelCount <
            diagnostics->carvedVoxelCountThisStep
        )
        {
            return fail(
                "carved more voxels than the candidate set contains"
            );
        }

        if (
            diagnostics->totalCarvedVoxelCount <
            previousCarvedVoxelCount
        )
        {
            return fail(
                "total carved-voxel count decreased"
            );
        }

        if (
            diagnostics->totalCarvedVoxelCount !=
            solver.carvedVoxelCount()
        )
        {
            return fail(
                "diagnostic carved count does not match solver state"
            );
        }

        if (
            !isFiniteVector(
                diagnostics->centerOfMassAfter
            ) ||
            !isFiniteVector(
                diagnostics->approximateVoxelCenterOfMassAfter
            ) ||
            !std::isfinite(
                diagnostics->legacyProjectedErrorAfter
            ) ||
            !std::isfinite(
                diagnostics->geometricProjectedErrorAfter
            ) ||
            !std::isfinite(
                diagnostics->centerOfMassDiscrepancyAfter
            ) ||
            !std::isfinite(
                diagnostics->approximateRemainingVoxelMass
            )
        )
        {
            return fail(
                "a diagnostic metric is non-finite"
            );
        }

        if (
            diagnostics->approximateRemainingVoxelMass < 0.0
        )
        {
            return fail(
                "approximate remaining voxel mass is negative"
            );
        }

        previousCarvedVoxelCount =
            diagnostics->totalCarvedVoxelCount;
    }


    if (executedSteps == 0)
    {
        return fail(
            "solver was already optimized; regression step was not exercised"
        );
    }

    if (
        solver.diagnosticsHistory().size() !=
        static_cast<std::size_t>(executedSteps)
    )
    {
        return fail(
            "diagnostic-history size does not match executed steps"
        );
    }

    if (!solver.writeDiagnosticsCsv(csvPath))
    {
        return fail(
            "could not write diagnostics CSV: " + csvPath
        );
    }


    std::cout
        << "[PASS] mesh=" << meshPath
        << " steps=" << executedSteps
        << " carved_voxels="
        << solver.carvedVoxelCount()
        << " legacy_error="
        << solver.getCoMDistance()
        << " geometric_error="
        << solver.getGeometricProjectedDistance()
        << " stop="
        << carvingStopReasonName(
            solver.stopReason()
        )
        << " csv=" << csvPath
        << '\n';

    return EXIT_SUCCESS;
}
