#include <cstdint>

#include "carving/CarvingSolver.h"
#include "grid_util.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <igl/read_triangle_mesh.h>


namespace
{

enum class RegressionScenario
{
    InitialState,
    OneStep,
    FullKnownBaseline,
    NoCandidates
};


const char* const CSV_HEADER =
    "iteration,"
    "com_before_x,com_before_y,com_before_z,"
    "com_after_x,com_after_y,com_after_z,"
    "voxel_com_after_x,voxel_com_after_y,voxel_com_after_z,"
    "legacy_error_before,legacy_error_after,"
    "geometric_error_before,geometric_error_after,"
    "com_discrepancy_after,"
    "exact_initial_mass,approximate_remaining_voxel_mass,"
    "candidate_voxel_count,carved_voxel_count_step,"
    "total_carved_voxel_count,"
    "initial_inside_voxel_count,initial_carveable_voxel_count,"
    "remaining_inside_voxel_count,"
    "carved_grid_ratio,carved_carveable_ratio,stop_reason";

constexpr std::size_t CSV_COLUMN_COUNT = 26;
constexpr double CONTINUITY_TOLERANCE = 1e-12;
constexpr double BASELINE_TOLERANCE = 1e-6;


struct VoxelStateSnapshot
{
    int filledInsideVoxelCount;
    double mass;
    Eigen::Vector3d centerOfMass;
};


bool parseScenario(
    const std::string& value,
    RegressionScenario& scenario
)
{
    if (value == "initial")
    {
        scenario = RegressionScenario::InitialState;
        return true;
    }

    if (value == "one-step")
    {
        scenario = RegressionScenario::OneStep;
        return true;
    }

    if (value == "full-known-baseline")
    {
        scenario = RegressionScenario::FullKnownBaseline;
        return true;
    }

    if (value == "no-candidates")
    {
        scenario = RegressionScenario::NoCandidates;
        return true;
    }

    return false;
}


const char* scenarioName(
    RegressionScenario scenario
)
{
    switch (scenario)
    {
        case RegressionScenario::InitialState:
            return "initial";

        case RegressionScenario::OneStep:
            return "one-step";

        case RegressionScenario::FullKnownBaseline:
            return "full-known-baseline";

        case RegressionScenario::NoCandidates:
            return "no-candidates";
    }

    return "unknown";
}


bool isFiniteVector(
    const Eigen::Vector3d& value
)
{
    return
        std::isfinite(value(0)) &&
        std::isfinite(value(1)) &&
        std::isfinite(value(2));
}


bool approximatelyEqual(
    double first,
    double second,
    double tolerance
)
{
    return std::abs(first - second) <= tolerance;
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


void removeAllCarveableVoxels(
    VoxelGrid& voxelGrid
)
{
    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                voxelGrid(x, y, z).isInner = false;
            }
        }
    }
}


VoxelStateSnapshot captureVoxelState(
    const VoxelGrid& voxelGrid
)
{
    VoxelStateSnapshot snapshot;
    snapshot.filledInsideVoxelCount = 0;
    snapshot.mass = 0.0;
    snapshot.centerOfMass = Eigen::Vector3d::Zero();

    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    Eigen::Vector3d weightedCenter =
        Eigen::Vector3d::Zero();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel =
                    voxelGrid(x, y, z);

                if (
                    voxel.distanceToMesh < 0.0 &&
                    voxel.isFilled
                )
                {
                    ++snapshot.filledInsideVoxelCount;
                    weightedCenter += voxel.center;
                }
            }
        }
    }

    const Eigen::Vector3d& voxelSize =
        voxelGrid.voxelSize();

    snapshot.mass =
        snapshot.filledInsideVoxelCount *
        voxelSize(0) * voxelSize(1) * voxelSize(2);

    if (snapshot.filledInsideVoxelCount > 0)
    {
        snapshot.centerOfMass =
            weightedCenter /
            static_cast<double>(
                snapshot.filledInsideVoxelCount
            );
    }

    return snapshot;
}


std::vector<std::string> splitCsvLine(
    const std::string& line
)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;

    while (std::getline(stream, field, ','))
    {
        fields.push_back(field);
    }

    return fields;
}


bool parseFiniteNumber(
    const std::string& value,
    double& parsed
)
{
    errno = 0;

    char* end = nullptr;
    parsed = std::strtod(value.c_str(), &end);

    return
        end != value.c_str() &&
        end != nullptr &&
        *end == '\0' &&
        errno != ERANGE &&
        std::isfinite(parsed);
}


bool parseInteger(
    const std::string& value,
    int& parsed
)
{
    errno = 0;

    char* end = nullptr;
    const long result =
        std::strtol(value.c_str(), &end, 10);

    if (
        end == value.c_str() ||
        end == nullptr ||
        *end != '\0' ||
        errno == ERANGE ||
        result < std::numeric_limits<int>::min() ||
        result > std::numeric_limits<int>::max()
    )
    {
        return false;
    }

    parsed = static_cast<int>(result);
    return true;
}


bool csvIntegerMatches(
    const std::vector<std::string>& fields,
    std::size_t column,
    int expected,
    std::string& error
)
{
    int parsed = 0;

    if (
        !parseInteger(fields[column], parsed) ||
        parsed != expected
    )
    {
        error =
            "diagnostics CSV integer column " +
            std::to_string(column) +
            " disagrees with solver history";
        return false;
    }

    return true;
}


bool csvDoubleMatches(
    const std::vector<std::string>& fields,
    std::size_t column,
    double expected,
    std::string& error
)
{
    double parsed = 0.0;

    const double tolerance =
        4.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(expected));

    if (
        !parseFiniteNumber(fields[column], parsed) ||
        !approximatelyEqual(parsed, expected, tolerance)
    )
    {
        error =
            "diagnostics CSV floating-point column " +
            std::to_string(column) +
            " disagrees with solver history";
        return false;
    }

    return true;
}


bool validateCsv(
    const std::string& csvPath,
    const std::vector<CarvingStepDiagnostics>& history,
    std::string& error
)
{
    std::ifstream input(csvPath.c_str());

    if (!input)
    {
        error = "could not reopen diagnostics CSV: " + csvPath;
        return false;
    }

    std::string line;

    if (!std::getline(input, line))
    {
        error = "diagnostics CSV is empty";
        return false;
    }

    if (!line.empty() && line.back() == '\r')
    {
        line.pop_back();
    }

    if (line != CSV_HEADER)
    {
        error = "diagnostics CSV header does not match the 26-column schema";
        return false;
    }

    std::size_t rowIndex = 0;

    while (std::getline(input, line))
    {
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }

        if (rowIndex >= history.size())
        {
            error = "diagnostics CSV contains more rows than solver history";
            return false;
        }

        const std::vector<std::string> fields =
            splitCsvLine(line);

        if (fields.size() != CSV_COLUMN_COUNT)
        {
            error = "diagnostics CSV data row does not contain 26 columns";
            return false;
        }

        const CarvingStepDiagnostics& diagnostics =
            history[rowIndex];

        if (
            !csvIntegerMatches(fields, 0, diagnostics.iteration, error) ||
            !csvDoubleMatches(fields, 1, diagnostics.centerOfMassBefore(0), error) ||
            !csvDoubleMatches(fields, 2, diagnostics.centerOfMassBefore(1), error) ||
            !csvDoubleMatches(fields, 3, diagnostics.centerOfMassBefore(2), error) ||
            !csvDoubleMatches(fields, 4, diagnostics.centerOfMassAfter(0), error) ||
            !csvDoubleMatches(fields, 5, diagnostics.centerOfMassAfter(1), error) ||
            !csvDoubleMatches(fields, 6, diagnostics.centerOfMassAfter(2), error) ||
            !csvDoubleMatches(fields, 7, diagnostics.approximateVoxelCenterOfMassAfter(0), error) ||
            !csvDoubleMatches(fields, 8, diagnostics.approximateVoxelCenterOfMassAfter(1), error) ||
            !csvDoubleMatches(fields, 9, diagnostics.approximateVoxelCenterOfMassAfter(2), error) ||
            !csvDoubleMatches(fields, 10, diagnostics.legacyProjectedErrorBefore, error) ||
            !csvDoubleMatches(fields, 11, diagnostics.legacyProjectedErrorAfter, error) ||
            !csvDoubleMatches(fields, 12, diagnostics.geometricProjectedErrorBefore, error) ||
            !csvDoubleMatches(fields, 13, diagnostics.geometricProjectedErrorAfter, error) ||
            !csvDoubleMatches(fields, 14, diagnostics.centerOfMassDiscrepancyAfter, error) ||
            !csvDoubleMatches(fields, 15, diagnostics.exactInitialMass, error) ||
            !csvDoubleMatches(fields, 16, diagnostics.approximateRemainingVoxelMass, error) ||
            !csvIntegerMatches(fields, 17, diagnostics.candidateVoxelCount, error) ||
            !csvIntegerMatches(fields, 18, diagnostics.carvedVoxelCountThisStep, error) ||
            !csvIntegerMatches(fields, 19, diagnostics.totalCarvedVoxelCount, error) ||
            !csvIntegerMatches(fields, 20, diagnostics.initialInsideVoxelCount, error) ||
            !csvIntegerMatches(fields, 21, diagnostics.initialCarveableVoxelCount, error) ||
            !csvIntegerMatches(fields, 22, diagnostics.remainingInsideVoxelCount, error) ||
            !csvDoubleMatches(fields, 23, diagnostics.carvedGridRatio, error) ||
            !csvDoubleMatches(fields, 24, diagnostics.carvedCarveableRatio, error)
        )
        {
            return false;
        }

        if (
            fields.back() !=
            carvingStopReasonName(history[rowIndex].stopReason)
        )
        {
            error = "diagnostics CSV stop reason disagrees with solver history";
            return false;
        }

        ++rowIndex;
    }

    if (rowIndex != history.size())
    {
        error = "diagnostics CSV contains fewer rows than solver history";
        return false;
    }

    return true;
}


bool validateHistory(
    const CarvingSolver& solver,
    std::string& error
)
{
    const std::vector<CarvingStepDiagnostics>& history =
        solver.diagnosticsHistory();

    int previousTotalCarved = 0;

    for (std::size_t index = 0; index < history.size(); ++index)
    {
        const CarvingStepDiagnostics& diagnostics =
            history[index];

        if (
            diagnostics.iteration !=
            static_cast<int>(index + 1)
        )
        {
            error = "diagnostic iterations are not consecutive";
            return false;
        }

        if (
            diagnostics.candidateVoxelCount < 0 ||
            diagnostics.carvedVoxelCountThisStep < 0 ||
            diagnostics.candidateVoxelCount <
                diagnostics.carvedVoxelCountThisStep
        )
        {
            error = "candidate/carved count invariant failed";
            return false;
        }

        if (
            diagnostics.totalCarvedVoxelCount !=
            previousTotalCarved +
                diagnostics.carvedVoxelCountThisStep
        )
        {
            error = "total carved count is not the sum of step removals";
            return false;
        }

        if (
            diagnostics.initialInsideVoxelCount < 0 ||
            diagnostics.initialCarveableVoxelCount < 0 ||
            diagnostics.remainingInsideVoxelCount < 0 ||
            diagnostics.remainingInsideVoxelCount >
                diagnostics.initialInsideVoxelCount
        )
        {
            error = "inside-voxel count invariant failed";
            return false;
        }

        if (
            !isFiniteVector(diagnostics.centerOfMassBefore) ||
            !isFiniteVector(diagnostics.centerOfMassAfter) ||
            !isFiniteVector(
                diagnostics.approximateVoxelCenterOfMassAfter
            ) ||
            !std::isfinite(diagnostics.legacyProjectedErrorBefore) ||
            !std::isfinite(diagnostics.legacyProjectedErrorAfter) ||
            !std::isfinite(diagnostics.geometricProjectedErrorBefore) ||
            !std::isfinite(diagnostics.geometricProjectedErrorAfter) ||
            !std::isfinite(diagnostics.centerOfMassDiscrepancyAfter) ||
            !std::isfinite(diagnostics.exactInitialMass) ||
            !std::isfinite(
                diagnostics.approximateRemainingVoxelMass
            ) ||
            !std::isfinite(diagnostics.carvedGridRatio) ||
            !std::isfinite(diagnostics.carvedCarveableRatio)
        )
        {
            error = "diagnostic history contains a non-finite value";
            return false;
        }

        if (
            diagnostics.approximateRemainingVoxelMass < 0.0 ||
            diagnostics.carvedGridRatio < 0.0 ||
            diagnostics.carvedGridRatio > 1.0 + CONTINUITY_TOLERANCE ||
            diagnostics.carvedCarveableRatio < 0.0 ||
            diagnostics.carvedCarveableRatio >
                1.0 + CONTINUITY_TOLERANCE
        )
        {
            error = "diagnostic mass or ratio is outside its valid range";
            return false;
        }

        if (index > 0)
        {
            const CarvingStepDiagnostics& previous =
                history[index - 1];

            if (
                (
                    diagnostics.centerOfMassBefore -
                    previous.centerOfMassAfter
                ).norm() > CONTINUITY_TOLERANCE ||
                !approximatelyEqual(
                    diagnostics.legacyProjectedErrorBefore,
                    previous.legacyProjectedErrorAfter,
                    CONTINUITY_TOLERANCE
                ) ||
                !approximatelyEqual(
                    diagnostics.geometricProjectedErrorBefore,
                    previous.geometricProjectedErrorAfter,
                    CONTINUITY_TOLERANCE
                )
            )
            {
                error = "diagnostic before/after continuity invariant failed";
                return false;
            }
        }

        if (
            index + 1 < history.size() &&
            diagnostics.stopReason != CarvingStopReason::None
        )
        {
            error = "a non-final diagnostic row reports a stop reason";
            return false;
        }

        previousTotalCarved =
            diagnostics.totalCarvedVoxelCount;
    }

    if (
        solver.iterationCount() !=
        static_cast<int>(history.size())
    )
    {
        error = "solver iteration count disagrees with diagnostic history";
        return false;
    }

    if (
        solver.carvedVoxelCount() !=
        previousTotalCarved
    )
    {
        error = "solver carved count disagrees with diagnostic history";
        return false;
    }

    if (history.empty())
    {
        if (solver.lastDiagnostics() != nullptr)
        {
            error = "zero-step solver unexpectedly exposes last diagnostics";
            return false;
        }
    }
    else
    {
        if (solver.lastDiagnostics() != &history.back())
        {
            error = "lastDiagnostics does not reference the final record";
            return false;
        }

        if (history.back().stopReason != solver.stopReason())
        {
            error = "final diagnostic stop reason disagrees with solver";
            return false;
        }

        if (
            (
                history.back().centerOfMassAfter -
                solver.currentCenterOfMass()
            ).norm() > CONTINUITY_TOLERANCE ||
            !approximatelyEqual(
                history.back().legacyProjectedErrorAfter,
                solver.getCoMDistance(),
                CONTINUITY_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.back().geometricProjectedErrorAfter,
                solver.getGeometricProjectedDistance(),
                CONTINUITY_TOLERANCE
            )
        )
        {
            error = "final diagnostics disagree with solver getters";
            return false;
        }
    }

    return true;
}

} // namespace


int main(
    int argc,
    char* argv[]
)
{
    if (argc != 4)
    {
        return fail(
            "usage: carving-regression <scenario> <mesh> <csv>"
        );
    }

    RegressionScenario scenario;

    if (!parseScenario(argv[1], scenario))
    {
        return fail(
            "unknown regression scenario: " +
            std::string(argv[1])
        );
    }

    const std::string meshPath = argv[2];
    const std::string csvPath = argv[3];

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
     * The deterministic horizontal offset avoids a degenerate cutting-plane
     * triangle if the initial CoM projects exactly to the plane center.
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

    if (scenario == RegressionScenario::NoCandidates)
    {
        removeAllCarveableVoxels(voxelGrid);
    }

    const VoxelStateSnapshot initialVoxelState =
        captureVoxelState(voxelGrid);

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
        !isFiniteVector(solver.currentCenterOfMass()) ||
        !std::isfinite(solver.getCoMDistance()) ||
        !std::isfinite(
            solver.getGeometricProjectedDistance()
        )
    )
    {
        return fail(
            "initial CoM or projected error is non-finite"
        );
    }

    const Eigen::Vector3d initialCenterOfMass =
        solver.currentCenterOfMass();

    const double initialLegacyError =
        solver.getCoMDistance();

    const double initialGeometricError =
        solver.getGeometricProjectedDistance();

    Eigen::Vector3d oneStepReturnedCenterOfMass =
        Eigen::Vector3d::Zero();

    if (scenario == RegressionScenario::OneStep)
    {
        Eigen::MatrixXd carvePlaneVertices;
        Eigen::MatrixXi carvePlaneFaces;
        solver.carveStep(
            carvePlaneVertices,
            carvePlaneFaces,
            oneStepReturnedCenterOfMass
        );
    }
    else if (
        scenario == RegressionScenario::FullKnownBaseline ||
        scenario == RegressionScenario::NoCandidates
    )
    {
        solver.carveAll();
    }

    std::string validationError;

    if (!validateHistory(solver, validationError))
    {
        return fail(validationError);
    }

    const std::vector<CarvingStepDiagnostics>& history =
        solver.diagnosticsHistory();

    const VoxelStateSnapshot finalVoxelState =
        captureVoxelState(voxelGrid);

    if (scenario == RegressionScenario::InitialState)
    {
        if (
            solver.iterationCount() != 0 ||
            solver.carvedVoxelCount() != 0 ||
            !history.empty() ||
            solver.isOptimized()
        )
        {
            return fail(
                "zero-step state changed or reported optimization"
            );
        }
    }
    else if (scenario == RegressionScenario::OneStep)
    {
        if (
            history.size() != 1 ||
            history.front().candidateVoxelCount <= 0 ||
            history.front().carvedVoxelCountThisStep <= 0 ||
            history.front().stopReason != CarvingStopReason::None ||
            (
                history.front().centerOfMassBefore -
                initialCenterOfMass
            ).norm() > CONTINUITY_TOLERANCE ||
            !approximatelyEqual(
                history.front().legacyProjectedErrorBefore,
                initialLegacyError,
                CONTINUITY_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.front().geometricProjectedErrorBefore,
                initialGeometricError,
                CONTINUITY_TOLERANCE
            ) ||
            (
                oneStepReturnedCenterOfMass -
                solver.currentCenterOfMass()
            ).norm() > CONTINUITY_TOLERANCE
        )
        {
            return fail(
                "one-step scenario did not exercise ordinary carving"
            );
        }
    }
    else if (scenario == RegressionScenario::NoCandidates)
    {
        if (
            history.size() != 1 ||
            history.front().candidateVoxelCount != 0 ||
            history.front().carvedVoxelCountThisStep != 0 ||
            history.front().totalCarvedVoxelCount != 0 ||
            history.front().initialCarveableVoxelCount != 0 ||
            history.front().remainingInsideVoxelCount !=
                history.front().initialInsideVoxelCount ||
            history.front().stopReason !=
                CarvingStopReason::NoCandidatesAvailable ||
            !solver.isOptimized() ||
            (
                history.front().centerOfMassAfter -
                initialCenterOfMass
            ).norm() > CONTINUITY_TOLERANCE ||
            !approximatelyEqual(
                history.front().legacyProjectedErrorAfter,
                initialLegacyError,
                CONTINUITY_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.front().geometricProjectedErrorAfter,
                initialGeometricError,
                CONTINUITY_TOLERANCE
            ) ||
            finalVoxelState.filledInsideVoxelCount !=
                initialVoxelState.filledInsideVoxelCount ||
            !approximatelyEqual(
                finalVoxelState.mass,
                initialVoxelState.mass,
                CONTINUITY_TOLERANCE
            ) ||
            (
                finalVoxelState.centerOfMass -
                initialVoxelState.centerOfMass
            ).norm() > CONTINUITY_TOLERANCE ||
            history.front().remainingInsideVoxelCount !=
                finalVoxelState.filledInsideVoxelCount ||
            !approximatelyEqual(
                history.front().approximateRemainingVoxelMass,
                finalVoxelState.mass,
                CONTINUITY_TOLERANCE
            ) ||
            (
                history.front().approximateVoxelCenterOfMassAfter -
                finalVoxelState.centerOfMass
            ).norm() > CONTINUITY_TOLERANCE ||
            history.front().carvedGridRatio != 0.0 ||
            history.front().carvedCarveableRatio != 0.0
        )
        {
            return fail(
                "empty candidate set did not terminate without mutation"
            );
        }
    }
    else
    {
        if (
            history.size() != 17 ||
            history.back().iteration != 17 ||
            history.back().candidateVoxelCount != 0 ||
            history.back().carvedVoxelCountThisStep != 0 ||
            history.back().totalCarvedVoxelCount != 14028 ||
            history.back().stopReason !=
                CarvingStopReason::NoCandidatesAvailable ||
            !solver.isOptimized()
        )
        {
            return fail(
                "full run did not stop at the frozen candidate-exhaustion state"
            );
        }

        for (
            std::size_t index = 0;
            index + 1 < history.size();
            ++index
        )
        {
            if (history[index].candidateVoxelCount <= 0)
            {
                return fail(
                    "candidate exhaustion occurred before iteration 17"
                );
            }
        }

        if (
            !approximatelyEqual(
                history.front().geometricProjectedErrorBefore,
                0.028996164626118211,
                BASELINE_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.back().geometricProjectedErrorAfter,
                0.81139676456917131,
                BASELINE_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.back().centerOfMassDiscrepancyAfter,
                0.80074263484909125,
                BASELINE_TOLERANCE
            ) ||
            history.back().geometricProjectedErrorAfter <=
                history.front().geometricProjectedErrorBefore
        )
        {
            return fail(
                "known divergent R0 numerical baseline changed"
            );
        }

        if (
            history.size() < 2 ||
            (
                history.back().centerOfMassAfter -
                history.back().centerOfMassBefore
            ).norm() > CONTINUITY_TOLERANCE ||
            !approximatelyEqual(
                history.back().legacyProjectedErrorAfter,
                history.back().legacyProjectedErrorBefore,
                CONTINUITY_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.back().geometricProjectedErrorAfter,
                history.back().geometricProjectedErrorBefore,
                CONTINUITY_TOLERANCE
            ) ||
            history.back().remainingInsideVoxelCount !=
                history[history.size() - 2].remainingInsideVoxelCount ||
            !approximatelyEqual(
                history.back().approximateRemainingVoxelMass,
                history[history.size() - 2].approximateRemainingVoxelMass,
                CONTINUITY_TOLERANCE
            ) ||
            (
                history.back().approximateVoxelCenterOfMassAfter -
                history[history.size() - 2].approximateVoxelCenterOfMassAfter
            ).norm() > CONTINUITY_TOLERANCE ||
            !approximatelyEqual(
                history.back().carvedGridRatio,
                history[history.size() - 2].carvedGridRatio,
                CONTINUITY_TOLERANCE
            ) ||
            !approximatelyEqual(
                history.back().carvedCarveableRatio,
                history[history.size() - 2].carvedCarveableRatio,
                CONTINUITY_TOLERANCE
            ) ||
            finalVoxelState.filledInsideVoxelCount !=
                history.back().remainingInsideVoxelCount ||
            !approximatelyEqual(
                finalVoxelState.mass,
                history.back().approximateRemainingVoxelMass,
                CONTINUITY_TOLERANCE
            ) ||
            (
                finalVoxelState.centerOfMass -
                history.back().approximateVoxelCenterOfMassAfter
            ).norm() > CONTINUITY_TOLERANCE
        )
        {
            return fail(
                "candidate-exhaustion probe changed the carving state"
            );
        }
    }

    if (
        scenario == RegressionScenario::FullKnownBaseline ||
        scenario == RegressionScenario::NoCandidates
    )
    {
        const int terminalIterationCount =
            solver.iterationCount();

        const int terminalCarvedVoxelCount =
            solver.carvedVoxelCount();

        const std::size_t terminalHistorySize =
            history.size();

        const Eigen::Vector3d terminalCenterOfMass =
            solver.currentCenterOfMass();

        solver.carveAll();

        if (
            solver.iterationCount() != terminalIterationCount ||
            solver.carvedVoxelCount() != terminalCarvedVoxelCount ||
            history.size() != terminalHistorySize ||
            (
                solver.currentCenterOfMass() -
                terminalCenterOfMass
            ).norm() > CONTINUITY_TOLERANCE
        )
        {
            return fail(
                "carveAll is not idempotent after a terminal state"
            );
        }
    }

    if (!solver.writeDiagnosticsCsv(csvPath))
    {
        return fail(
            "could not write diagnostics CSV: " + csvPath
        );
    }

    if (!validateCsv(csvPath, history, validationError))
    {
        return fail(validationError);
    }

    std::cout
        << "[PASS] scenario=" << scenarioName(scenario)
        << " mesh=" << meshPath
        << " steps=" << solver.iterationCount()
        << " carved_voxels=" << solver.carvedVoxelCount()
        << " legacy_error=" << solver.getCoMDistance()
        << " geometric_error="
        << solver.getGeometricProjectedDistance()
        << " stop="
        << carvingStopReasonName(solver.stopReason())
        << " csv=" << csvPath
        << '\n';

    return EXIT_SUCCESS;
}
