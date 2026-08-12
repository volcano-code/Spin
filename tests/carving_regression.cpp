#include <cstdint>

#include "carving/CarvingSolver.h"
#include "grid_util.h"
#include "physics/BalanceEnergy.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
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
    FullMassConsistency,
    NoCandidates
};


const char* const CSV_HEADER =
    "iteration,"
    "com_before_x,com_before_y,com_before_z,"
    "com_after_x,com_after_y,com_after_z,"
    "recomputed_com_after_x,recomputed_com_after_y,"
    "recomputed_com_after_z,"
    "first_moment_before_x,first_moment_before_y,"
    "first_moment_before_z,"
    "first_moment_after_x,first_moment_after_y,"
    "first_moment_after_z,"
    "recomputed_first_moment_after_x,"
    "recomputed_first_moment_after_y,"
    "recomputed_first_moment_after_z,"
    "projected_distance_before,projected_distance_after,"
    "projected_energy_before,projected_energy_after,"
    "com_discrepancy_after,first_moment_discrepancy_after,"
    "initial_mass,remaining_mass_before,remaining_mass_after,"
    "recomputed_remaining_mass_after,mass_discrepancy_after,"
    "candidate_voxel_count,carved_voxel_count_step,"
    "total_carved_voxel_count,"
    "initial_inside_voxel_count,initial_carveable_voxel_count,"
    "remaining_inside_voxel_count,"
    "carved_grid_ratio,carved_carveable_ratio,stop_reason";

constexpr std::size_t CSV_COLUMN_COUNT = 39;
constexpr int MAXIMUM_TEST_STEPS = 100;

constexpr double MASS_MOMENT_ABSOLUTE_TOLERANCE = 1e-12;
constexpr double MASS_MOMENT_RELATIVE_TOLERANCE = 1e-10;
constexpr double COM_ABSOLUTE_TOLERANCE = 1e-10;
constexpr double COM_RELATIVE_TOLERANCE = 1e-10;
constexpr double ENERGY_ABSOLUTE_TOLERANCE = 1e-12;
constexpr double ENERGY_RELATIVE_TOLERANCE = 1e-10;


struct LongDoubleMassProperties
{
    long double mass;
    std::array<long double, 3> firstMoment;
    int filledInsideVoxelCount;
    std::vector<unsigned char> filledFlags;

    LongDoubleMassProperties()
        :
        mass(0.0L),
        firstMoment{{0.0L, 0.0L, 0.0L}},
        filledInsideVoxelCount(0),
        filledFlags()
    {
    }
};


struct FixtureGeometry
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    Eigen::MatrixXd planeVertices;
    Eigen::MatrixXi planeFaces;
    Eigen::Vector3d targetCenterOfMass;
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

    if (value == "full-mass-consistency")
    {
        scenario = RegressionScenario::FullMassConsistency;
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

        case RegressionScenario::FullMassConsistency:
            return "full-mass-consistency";

        case RegressionScenario::NoCandidates:
            return "no-candidates";
    }

    return "unknown";
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


bool isFiniteVector(
    const Eigen::Vector3d& value
)
{
    return value.allFinite();
}


bool approximatelyEqual(
    double actual,
    long double expected,
    double absoluteTolerance,
    double relativeTolerance
)
{
    if (!std::isfinite(actual) || !std::isfinite(expected))
    {
        return false;
    }

    const long double difference =
        std::abs(
            static_cast<long double>(actual) -
            expected
        );

    const long double scale =
        std::max(
            std::abs(static_cast<long double>(actual)),
            std::abs(expected)
        );

    return difference <=
        static_cast<long double>(absoluteTolerance) +
        static_cast<long double>(relativeTolerance) * scale;
}


bool massMomentEqual(
    double actual,
    long double expected
)
{
    return approximatelyEqual(
        actual,
        expected,
        MASS_MOMENT_ABSOLUTE_TOLERANCE,
        MASS_MOMENT_RELATIVE_TOLERANCE
    );
}


bool centerEqual(
    double actual,
    long double expected
)
{
    return approximatelyEqual(
        actual,
        expected,
        COM_ABSOLUTE_TOLERANCE,
        COM_RELATIVE_TOLERANCE
    );
}


bool centerEqual(
    const Eigen::Vector3d& actual,
    const std::array<long double, 3>& expected
)
{
    return
        centerEqual(actual(0), expected[0]) &&
        centerEqual(actual(1), expected[1]) &&
        centerEqual(actual(2), expected[2]);
}


bool vectorMassMomentEqual(
    const Eigen::Vector3d& actual,
    const std::array<long double, 3>& expected
)
{
    return
        massMomentEqual(actual(0), expected[0]) &&
        massMomentEqual(actual(1), expected[1]) &&
        massMomentEqual(actual(2), expected[2]);
}


bool energyEqual(
    double actual,
    long double expected
)
{
    return approximatelyEqual(
        actual,
        expected,
        ENERGY_ABSOLUTE_TOLERANCE,
        ENERGY_RELATIVE_TOLERANCE
    );
}


std::array<long double, 3> oracleCenterOfMass(
    const LongDoubleMassProperties& state
)
{
    return
    {{
        state.firstMoment[0] / state.mass,
        state.firstMoment[1] / state.mass,
        state.firstMoment[2] / state.mass
    }};
}


Eigen::RowVector3d createGroundPlane(
    const Eigen::MatrixXd& vertices,
    Eigen::MatrixXd& planeVertices,
    Eigen::MatrixXi& planeFaces
)
{
    Eigen::AlignedBox3d boundingBox;
    createAlignedBox(vertices, boundingBox);

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


bool prepareGeometry(
    const std::string& meshPath,
    FixtureGeometry& geometry,
    std::string& error
)
{
    Eigen::MatrixXd inputVertices;

    if (
        !igl::read_triangle_mesh(
            meshPath,
            inputVertices,
            geometry.faces
        )
    )
    {
        error = "could not load mesh: " + meshPath;
        return false;
    }

    if (
        inputVertices.rows() == 0 ||
        geometry.faces.rows() == 0
    )
    {
        error = "loaded mesh is empty";
        return false;
    }

    alignToAxis(
        inputVertices,
        geometry.vertices
    );

    const Eigen::RowVector3d planeCenter =
        createGroundPlane(
            geometry.vertices,
            geometry.planeVertices,
            geometry.planeFaces
        );

    Eigen::AlignedBox3d boundingBox;
    createAlignedBox(
        geometry.vertices,
        boundingBox
    );

    geometry.targetCenterOfMass =
        planeCenter.transpose();

    geometry.targetCenterOfMass(0) +=
        0.15 * boundingBox.sizes()(0);

    return true;
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


int countInsideVoxels(
    const VoxelGrid& voxelGrid
)
{
    int count = 0;
    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (voxelGrid(x, y, z).distanceToMesh < 0.0)
                {
                    ++count;
                }
            }
        }
    }

    return count;
}


int countCarveableVoxels(
    const VoxelGrid& voxelGrid
)
{
    int count = 0;
    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (voxelGrid(x, y, z).isInner)
                {
                    ++count;
                }
            }
        }
    }

    return count;
}


LongDoubleMassProperties captureOracle(
    const MassProperties& initialOuter,
    const VoxelGrid& voxelGrid
)
{
    LongDoubleMassProperties result;

    result.mass =
        static_cast<long double>(initialOuter.mass);

    for (int component = 0; component < 3; ++component)
    {
        result.firstMoment[component] =
            static_cast<long double>(
                initialOuter.firstMoment(component)
            );
    }

    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    const Eigen::Vector3d& voxelSize =
        voxelGrid.voxelSize();

    const long double voxelVolume =
        static_cast<long double>(voxelSize(0)) *
        static_cast<long double>(voxelSize(1)) *
        static_cast<long double>(voxelSize(2));

    result.filledFlags.reserve(
        static_cast<std::size_t>(dimensions.prod())
    );

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel =
                    voxelGrid(x, y, z);

                result.filledFlags.push_back(
                    voxel.isFilled ? 1U : 0U
                );

                if (
                    voxel.distanceToMesh < 0.0 &&
                    voxel.isFilled
                )
                {
                    ++result.filledInsideVoxelCount;
                }

                /*
                 * Only carveable inner voxels can be removed by this solver.
                 * The continuous outer surface remains the immutable positive
                 * contribution; each removed inner voxel is subtracted.
                 */
                if (
                    voxel.isInner &&
                    !voxel.isFilled
                )
                {
                    result.mass -= voxelVolume;

                    for (int component = 0; component < 3; ++component)
                    {
                        result.firstMoment[component] -=
                            voxelVolume *
                            static_cast<long double>(
                                voxel.center(component)
                            );
                    }
                }
            }
        }
    }

    return result;
}


bool validateMassProperties(
    const MassProperties& actual,
    const Eigen::Vector3d& actualCenterOfMass,
    const LongDoubleMassProperties& oracle,
    const std::string& label,
    std::string& error
)
{
    if (
        !actual.isFinite() ||
        actual.mass <= 0.0 ||
        !isFiniteVector(actualCenterOfMass) ||
        oracle.mass <= 0.0L
    )
    {
        error = label + " contains non-finite/non-positive mass properties";
        return false;
    }

    if (!massMomentEqual(actual.mass, oracle.mass))
    {
        error = label + " mass disagrees with the long-double oracle";
        return false;
    }

    if (
        !vectorMassMomentEqual(
            actual.firstMoment,
            oracle.firstMoment
        )
    )
    {
        error = label + " first moment disagrees with the long-double oracle";
        return false;
    }

    const std::array<long double, 3> expectedCenter =
        oracleCenterOfMass(oracle);

    if (!centerEqual(actualCenterOfMass, expectedCenter))
    {
        error = label + " center of mass disagrees with the oracle";
        return false;
    }

    Eigen::Vector3d derivedCenter =
        Eigen::Vector3d::Zero();

    if (
        !actual.tryCenterOfMass(derivedCenter) ||
        !centerEqual(derivedCenter, expectedCenter)
    )
    {
        error = label + " does not satisfy center=firstMoment/mass";
        return false;
    }

    return true;
}


bool validateEnergyPair(
    double distance,
    double energy,
    const std::string& label,
    std::string& error
)
{
    if (
        !std::isfinite(distance) ||
        distance < 0.0 ||
        !std::isfinite(energy) ||
        energy < 0.0 ||
        !energyEqual(
            energy,
            0.5L *
                static_cast<long double>(distance) *
                static_cast<long double>(distance)
        )
    )
    {
        error = label + " violates E=0.5*d^2";
        return false;
    }

    return true;
}


bool validateDiagnosticStep(
    const CarvingSolver& solver,
    const CarvingStepDiagnostics& diagnostics,
    const LongDoubleMassProperties& before,
    const LongDoubleMassProperties& after,
    const MassProperties& initialOuter,
    int expectedIteration,
    int previousTotalCarved,
    int initialInsideCount,
    int initialCarveableCount,
    std::string& error
)
{
    if (diagnostics.iteration != expectedIteration)
    {
        error = "diagnostic iteration is not consecutive";
        return false;
    }

    if (
        diagnostics.candidateVoxelCount < 0 ||
        diagnostics.carvedVoxelCountThisStep < 0 ||
        diagnostics.carvedVoxelCountThisStep >
            diagnostics.candidateVoxelCount ||
        diagnostics.totalCarvedVoxelCount !=
            previousTotalCarved +
            diagnostics.carvedVoxelCountThisStep ||
        diagnostics.totalCarvedVoxelCount !=
            solver.carvedVoxelCount()
    )
    {
        error = "candidate/carved count invariant failed";
        return false;
    }

    if (
        diagnostics.initialInsideVoxelCount != initialInsideCount ||
        diagnostics.initialCarveableVoxelCount != initialCarveableCount ||
        diagnostics.remainingInsideVoxelCount !=
            after.filledInsideVoxelCount ||
        before.filledInsideVoxelCount -
            after.filledInsideVoxelCount !=
            diagnostics.carvedVoxelCountThisStep
    )
    {
        error = "inside-voxel count disagrees with material state";
        return false;
    }

    if (
        !massMomentEqual(
            diagnostics.initialMass,
            static_cast<long double>(initialOuter.mass)
        ) ||
        !massMomentEqual(
            diagnostics.remainingMassBefore,
            before.mass
        ) ||
        !massMomentEqual(
            diagnostics.remainingMassAfter,
            after.mass
        ) ||
        !massMomentEqual(
            diagnostics.recomputedRemainingMassAfter,
            after.mass
        ) ||
        !vectorMassMomentEqual(
            diagnostics.firstMomentBefore,
            before.firstMoment
        ) ||
        !vectorMassMomentEqual(
            diagnostics.firstMomentAfter,
            after.firstMoment
        ) ||
        !vectorMassMomentEqual(
            diagnostics.recomputedFirstMomentAfter,
            after.firstMoment
        ) ||
        !centerEqual(
            diagnostics.centerOfMassBefore,
            oracleCenterOfMass(before)
        ) ||
        !centerEqual(
            diagnostics.centerOfMassAfter,
            oracleCenterOfMass(after)
        ) ||
        !centerEqual(
            diagnostics.recomputedCenterOfMassAfter,
            oracleCenterOfMass(after)
        )
    )
    {
        error = "diagnostic mass/first-moment/CoM disagrees with oracle";
        return false;
    }

    if (
        !validateMassProperties(
            solver.currentMassProperties(),
            solver.currentCenterOfMass(),
            after,
            "solver incremental state",
            error
        )
    )
    {
        return false;
    }

    if (
        !validateEnergyPair(
            diagnostics.projectedDistanceBefore,
            diagnostics.projectedEnergyBefore,
            "diagnostic before energy",
            error
        ) ||
        !validateEnergyPair(
            diagnostics.projectedDistanceAfter,
            diagnostics.projectedEnergyAfter,
            "diagnostic after energy",
            error
        ) ||
        !validateEnergyPair(
            solver.getCoMDistance(),
            solver.getCoMEnergy(),
            "solver energy",
            error
        )
    )
    {
        return false;
    }

    if (
        !centerEqual(
            solver.getCoMDistance(),
            static_cast<long double>(
                diagnostics.projectedDistanceAfter
            )
        ) ||
        !energyEqual(
            solver.getCoMEnergy(),
            static_cast<long double>(
                diagnostics.projectedEnergyAfter
            )
        )
    )
    {
        error = "final diagnostic energy disagrees with solver getters";
        return false;
    }

    const double expectedMassDiscrepancy =
        std::abs(
            diagnostics.remainingMassAfter -
            diagnostics.recomputedRemainingMassAfter
        );

    const double expectedFirstMomentDiscrepancy =
        (
            diagnostics.firstMomentAfter -
            diagnostics.recomputedFirstMomentAfter
        ).norm();

    const double expectedCenterDiscrepancy =
        (
            diagnostics.centerOfMassAfter -
            diagnostics.recomputedCenterOfMassAfter
        ).norm();

    if (
        !massMomentEqual(
            diagnostics.massDiscrepancyAfter,
            static_cast<long double>(expectedMassDiscrepancy)
        ) ||
        !massMomentEqual(
            diagnostics.firstMomentDiscrepancyAfter,
            static_cast<long double>(expectedFirstMomentDiscrepancy)
        ) ||
        !centerEqual(
            diagnostics.centerOfMassDiscrepancyAfter,
            static_cast<long double>(expectedCenterDiscrepancy)
        )
    )
    {
        error = "reported recompute discrepancy is internally inconsistent";
        return false;
    }

    if (
        !std::isfinite(diagnostics.carvedGridRatio) ||
        diagnostics.carvedGridRatio < 0.0 ||
        diagnostics.carvedGridRatio > 1.0 + COM_ABSOLUTE_TOLERANCE ||
        !std::isfinite(diagnostics.carvedCarveableRatio) ||
        diagnostics.carvedCarveableRatio < 0.0 ||
        diagnostics.carvedCarveableRatio >
            1.0 + COM_ABSOLUTE_TOLERANCE
    )
    {
        error = "carved ratio is outside its valid range";
        return false;
    }

    if (diagnostics.candidateVoxelCount == 0)
    {
        if (
            diagnostics.carvedVoxelCountThisStep != 0 ||
            before.filledFlags != after.filledFlags ||
            !massMomentEqual(
                diagnostics.remainingMassAfter,
                before.mass
            ) ||
            !vectorMassMomentEqual(
                diagnostics.firstMomentAfter,
                before.firstMoment
            ) ||
            !centerEqual(
                diagnostics.centerOfMassAfter,
                oracleCenterOfMass(before)
            ) ||
            diagnostics.stopReason !=
                CarvingStopReason::NoCandidatesAvailable
        )
        {
            error = "zero-candidate step changed material state";
            return false;
        }
    }

    return true;
}


std::vector<std::string> splitCsvLine(
    const std::string& line
)
{
    std::vector<std::string> fields;
    std::stringstream input(line);
    std::string field;

    while (std::getline(input, field, ','))
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
            "CSV integer column " +
            std::to_string(column) +
            " disagrees with diagnostic history";
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

    if (
        !parseFiniteNumber(fields[column], parsed) ||
        parsed != expected
    )
    {
        error =
            "CSV floating-point column " +
            std::to_string(column) +
            " disagrees with diagnostic history";
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
        error = "diagnostics CSV header does not match the 39-column schema";
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
            error = "diagnostics CSV contains excess rows";
            return false;
        }

        const std::vector<std::string> fields =
            splitCsvLine(line);

        if (fields.size() != CSV_COLUMN_COUNT)
        {
            error = "diagnostics CSV row does not contain 39 columns";
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
            !csvDoubleMatches(fields, 7, diagnostics.recomputedCenterOfMassAfter(0), error) ||
            !csvDoubleMatches(fields, 8, diagnostics.recomputedCenterOfMassAfter(1), error) ||
            !csvDoubleMatches(fields, 9, diagnostics.recomputedCenterOfMassAfter(2), error) ||
            !csvDoubleMatches(fields, 10, diagnostics.firstMomentBefore(0), error) ||
            !csvDoubleMatches(fields, 11, diagnostics.firstMomentBefore(1), error) ||
            !csvDoubleMatches(fields, 12, diagnostics.firstMomentBefore(2), error) ||
            !csvDoubleMatches(fields, 13, diagnostics.firstMomentAfter(0), error) ||
            !csvDoubleMatches(fields, 14, diagnostics.firstMomentAfter(1), error) ||
            !csvDoubleMatches(fields, 15, diagnostics.firstMomentAfter(2), error) ||
            !csvDoubleMatches(fields, 16, diagnostics.recomputedFirstMomentAfter(0), error) ||
            !csvDoubleMatches(fields, 17, diagnostics.recomputedFirstMomentAfter(1), error) ||
            !csvDoubleMatches(fields, 18, diagnostics.recomputedFirstMomentAfter(2), error) ||
            !csvDoubleMatches(fields, 19, diagnostics.projectedDistanceBefore, error) ||
            !csvDoubleMatches(fields, 20, diagnostics.projectedDistanceAfter, error) ||
            !csvDoubleMatches(fields, 21, diagnostics.projectedEnergyBefore, error) ||
            !csvDoubleMatches(fields, 22, diagnostics.projectedEnergyAfter, error) ||
            !csvDoubleMatches(fields, 23, diagnostics.centerOfMassDiscrepancyAfter, error) ||
            !csvDoubleMatches(fields, 24, diagnostics.firstMomentDiscrepancyAfter, error) ||
            !csvDoubleMatches(fields, 25, diagnostics.initialMass, error) ||
            !csvDoubleMatches(fields, 26, diagnostics.remainingMassBefore, error) ||
            !csvDoubleMatches(fields, 27, diagnostics.remainingMassAfter, error) ||
            !csvDoubleMatches(fields, 28, diagnostics.recomputedRemainingMassAfter, error) ||
            !csvDoubleMatches(fields, 29, diagnostics.massDiscrepancyAfter, error) ||
            !csvIntegerMatches(fields, 30, diagnostics.candidateVoxelCount, error) ||
            !csvIntegerMatches(fields, 31, diagnostics.carvedVoxelCountThisStep, error) ||
            !csvIntegerMatches(fields, 32, diagnostics.totalCarvedVoxelCount, error) ||
            !csvIntegerMatches(fields, 33, diagnostics.initialInsideVoxelCount, error) ||
            !csvIntegerMatches(fields, 34, diagnostics.initialCarveableVoxelCount, error) ||
            !csvIntegerMatches(fields, 35, diagnostics.remainingInsideVoxelCount, error) ||
            !csvDoubleMatches(fields, 36, diagnostics.carvedGridRatio, error) ||
            !csvDoubleMatches(fields, 37, diagnostics.carvedCarveableRatio, error)
        )
        {
            return false;
        }

        if (
            fields[38] !=
            carvingStopReasonName(diagnostics.stopReason)
        )
        {
            error = "CSV stop reason disagrees with diagnostic history";
            return false;
        }

        ++rowIndex;
    }

    if (rowIndex != history.size())
    {
        error = "diagnostics CSV contains too few rows";
        return false;
    }

    return true;
}


std::string serializeHistory(
    const std::vector<CarvingStepDiagnostics>& history
)
{
    std::ostringstream output;
    writeCarvingDiagnosticsCsvHeader(output);

    for (const CarvingStepDiagnostics& diagnostics : history)
    {
        writeCarvingDiagnosticsCsvRow(output, diagnostics);
    }

    return output.str();
}


bool isNormalTerminalReason(
    CarvingStopReason reason
)
{
    return
        reason == CarvingStopReason::TargetToleranceReached ||
        reason == CarvingStopReason::MaximumEmptyRatioReached ||
        reason == CarvingStopReason::NoCandidatesAvailable ||
        reason == CarvingStopReason::MaximumIterationsReached;
}


bool runToTerminal(
    CarvingSolver& solver,
    VoxelGrid& voxelGrid,
    const MassProperties& initialOuter,
    int initialInsideCount,
    int initialCarveableCount,
    std::string& error
)
{
    int previousTotalCarved = 0;
    int expectedIteration = 1;

    while (!solver.isOptimized())
    {
        if (expectedIteration > MAXIMUM_TEST_STEPS)
        {
            error = "full run exceeded the 100-step safety limit";
            return false;
        }

        const LongDoubleMassProperties before =
            captureOracle(initialOuter, voxelGrid);

        Eigen::MatrixXd carvePlaneVertices;
        Eigen::MatrixXi carvePlaneFaces;
        Eigen::Vector3d returnedCenterOfMass;

        solver.carveStep(
            carvePlaneVertices,
            carvePlaneFaces,
            returnedCenterOfMass
        );

        const LongDoubleMassProperties after =
            captureOracle(initialOuter, voxelGrid);

        const CarvingStepDiagnostics* diagnostics =
            solver.lastDiagnostics();

        if (diagnostics == nullptr)
        {
            error = "carveStep did not produce diagnostics";
            return false;
        }

        if (
            !validateDiagnosticStep(
                solver,
                *diagnostics,
                before,
                after,
                initialOuter,
                expectedIteration,
                previousTotalCarved,
                initialInsideCount,
                initialCarveableCount,
                error
            )
        )
        {
            return false;
        }

        if (
            !centerEqual(
                returnedCenterOfMass,
                oracleCenterOfMass(after)
            )
        )
        {
            error = "carveStep returned CoM disagrees with oracle";
            return false;
        }

        if (
            !solver.isOptimized() &&
            diagnostics->stopReason != CarvingStopReason::None
        )
        {
            error = "non-terminal step reports a stop reason";
            return false;
        }

        previousTotalCarved =
            diagnostics->totalCarvedVoxelCount;

        ++expectedIteration;
    }

    if (
        solver.diagnosticsHistory().empty() ||
        solver.iterationCount() > MAXIMUM_TEST_STEPS ||
        !isNormalTerminalReason(solver.stopReason()) ||
        solver.diagnosticsHistory().back().stopReason !=
            solver.stopReason()
    )
    {
        error = "full run did not reach a normal terminal state";
        return false;
    }

    return true;
}


bool validateTerminalIdempotence(
    CarvingSolver& solver,
    const VoxelGrid& voxelGrid,
    const MassProperties& initialOuter,
    std::string& error
)
{
    const int iterationCount = solver.iterationCount();
    const int carvedCount = solver.carvedVoxelCount();
    const std::size_t historySize =
        solver.diagnosticsHistory().size();

    const LongDoubleMassProperties before =
        captureOracle(initialOuter, voxelGrid);

    solver.carveAll();

    const LongDoubleMassProperties after =
        captureOracle(initialOuter, voxelGrid);

    if (
        solver.iterationCount() != iterationCount ||
        solver.carvedVoxelCount() != carvedCount ||
        solver.diagnosticsHistory().size() != historySize ||
        before.filledFlags != after.filledFlags ||
        !massMomentEqual(solver.objectMass(), before.mass) ||
        !vectorMassMomentEqual(
            solver.currentMassProperties().firstMoment,
            before.firstMoment
        ) ||
        !centerEqual(
            solver.currentCenterOfMass(),
            oracleCenterOfMass(before)
        )
    )
    {
        error = "carveAll is not idempotent after termination";
        return false;
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

    FixtureGeometry geometry;
    std::string error;

    if (!prepareGeometry(meshPath, geometry, error))
    {
        return fail(error);
    }

    VoxelGrid voxelGrid;
    voxelGrid.build(geometry.vertices, geometry.faces);

    if (scenario == RegressionScenario::NoCandidates)
    {
        removeAllCarveableVoxels(voxelGrid);
    }

    const int initialInsideCount =
        countInsideVoxels(voxelGrid);

    const int initialCarveableCount =
        countCarveableVoxels(voxelGrid);

    const MassProperties initialOuter =
        computeSurfaceMassProperties(
            geometry.vertices,
            geometry.faces
        );

    if (scenario == RegressionScenario::InitialState)
    {
        CarvingSolver unitGravitySolver(
            voxelGrid,
            geometry.vertices,
            geometry.faces,
            geometry.planeVertices,
            geometry.planeFaces,
            Eigen::Vector3d(0.0, -1.0, 0.0),
            geometry.targetCenterOfMass
        );

        CarvingSolver tinyGravitySolver(
            voxelGrid,
            geometry.vertices,
            geometry.faces,
            geometry.planeVertices,
            geometry.planeFaces,
            Eigen::Vector3d(0.0, -1e-300, 0.0),
            geometry.targetCenterOfMass
        );

        CarvingSolver hugeGravitySolver(
            voxelGrid,
            geometry.vertices,
            geometry.faces,
            geometry.planeVertices,
            geometry.planeFaces,
            Eigen::Vector3d(0.0, -1e308, 0.0),
            geometry.targetCenterOfMass
        );

        if (
            !energyEqual(
                tinyGravitySolver.getCoMDistance(),
                unitGravitySolver.getCoMDistance()
            ) ||
            !energyEqual(
                hugeGravitySolver.getCoMDistance(),
                unitGravitySolver.getCoMDistance()
            ) ||
            !energyEqual(
                tinyGravitySolver.getCoMEnergy(),
                unitGravitySolver.getCoMEnergy()
            ) ||
            !energyEqual(
                hugeGravitySolver.getCoMEnergy(),
                unitGravitySolver.getCoMEnergy()
            )
        )
        {
            return fail(
                "solver depends on finite non-zero gravity magnitude"
            );
        }
    }

    CarvingSolver solver(
        voxelGrid,
        geometry.vertices,
        geometry.faces,
        geometry.planeVertices,
        geometry.planeFaces,
        Eigen::Vector3d(0.0, -1.0, 0.0),
        geometry.targetCenterOfMass
    );

    solver.setDiagnosticsEnabled(false);

    const LongDoubleMassProperties initialOracle =
        captureOracle(initialOuter, voxelGrid);

    if (
        !validateMassProperties(
            solver.initialMassProperties(),
            solver.currentCenterOfMass(),
            initialOracle,
            "initial solver state",
            error
        ) ||
        !validateMassProperties(
            solver.currentMassProperties(),
            solver.currentCenterOfMass(),
            initialOracle,
            "current zero-step state",
            error
        ) ||
        !validateEnergyPair(
            solver.getCoMDistance(),
            solver.getCoMEnergy(),
            "initial solver energy",
            error
        )
    )
    {
        return fail(error);
    }

    if (scenario == RegressionScenario::InitialState)
    {
        if (
            solver.iterationCount() != 0 ||
            solver.carvedVoxelCount() != 0 ||
            !solver.diagnosticsHistory().empty() ||
            solver.lastDiagnostics() != nullptr ||
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
        const LongDoubleMassProperties before =
            captureOracle(initialOuter, voxelGrid);

        Eigen::MatrixXd carvePlaneVertices;
        Eigen::MatrixXi carvePlaneFaces;
        Eigen::Vector3d returnedCenterOfMass;

        solver.carveStep(
            carvePlaneVertices,
            carvePlaneFaces,
            returnedCenterOfMass
        );

        const LongDoubleMassProperties after =
            captureOracle(initialOuter, voxelGrid);

        const CarvingStepDiagnostics* diagnostics =
            solver.lastDiagnostics();

        if (
            diagnostics == nullptr ||
            !validateDiagnosticStep(
                solver,
                *diagnostics,
                before,
                after,
                initialOuter,
                1,
                0,
                initialInsideCount,
                initialCarveableCount,
                error
            ) ||
            diagnostics->candidateVoxelCount <= 0 ||
            diagnostics->carvedVoxelCountThisStep <= 0 ||
            diagnostics->stopReason != CarvingStopReason::None ||
            !centerEqual(
                returnedCenterOfMass,
                oracleCenterOfMass(after)
            )
        )
        {
            return fail(
                error.empty()
                    ? "one-step scenario did not exercise ordinary carving"
                    : error
            );
        }
    }
    else
    {
        if (
            !runToTerminal(
                solver,
                voxelGrid,
                initialOuter,
                initialInsideCount,
                initialCarveableCount,
                error
            ) ||
            !validateTerminalIdempotence(
                solver,
                voxelGrid,
                initialOuter,
                error
            )
        )
        {
            return fail(error);
        }

        if (scenario == RegressionScenario::NoCandidates)
        {
            if (
                solver.iterationCount() != 1 ||
                solver.carvedVoxelCount() != 0 ||
                solver.stopReason() !=
                    CarvingStopReason::NoCandidatesAvailable ||
                initialCarveableCount != 0 ||
                captureOracle(
                    initialOuter,
                    voxelGrid
                ).filledFlags != initialOracle.filledFlags
            )
            {
                return fail(
                    "initially empty candidate set did not stop without mutation"
                );
            }
        }
        else
        {
            /*
             * Repeat the complete live run. Exact serialized diagnostics and
             * voxel flags form the deterministic-state contract without
             * freezing any obsolete R0 trajectory value.
             */
            VoxelGrid repeatedGrid;
            repeatedGrid.build(
                geometry.vertices,
                geometry.faces
            );

            CarvingSolver repeatedSolver(
                repeatedGrid,
                geometry.vertices,
                geometry.faces,
                geometry.planeVertices,
                geometry.planeFaces,
                Eigen::Vector3d(0.0, -1.0, 0.0),
                geometry.targetCenterOfMass
            );

            if (
                !runToTerminal(
                    repeatedSolver,
                    repeatedGrid,
                    initialOuter,
                    countInsideVoxels(repeatedGrid),
                    countCarveableVoxels(repeatedGrid),
                    error
                ) ||
                serializeHistory(
                    repeatedSolver.diagnosticsHistory()
                ) !=
                    serializeHistory(
                        solver.diagnosticsHistory()
                    ) ||
                captureOracle(
                    initialOuter,
                    repeatedGrid
                ).filledFlags !=
                    captureOracle(
                        initialOuter,
                        voxelGrid
                    ).filledFlags
            )
            {
                return fail(
                    error.empty()
                        ? "two identical full runs produced different states"
                        : error
                );
            }
        }
    }

    if (!solver.writeDiagnosticsCsv(csvPath))
    {
        return fail(
            "could not write diagnostics CSV: " + csvPath
        );
    }

    if (
        !validateCsv(
            csvPath,
            solver.diagnosticsHistory(),
            error
        )
    )
    {
        return fail(error);
    }

    std::cout
        << "[PASS] scenario=" << scenarioName(scenario)
        << " mesh=" << meshPath
        << " steps=" << solver.iterationCount()
        << " carved_voxels=" << solver.carvedVoxelCount()
        << " mass=" << solver.objectMass()
        << " distance=" << solver.getCoMDistance()
        << " energy=" << solver.getCoMEnergy()
        << " stop=" << carvingStopReasonName(solver.stopReason())
        << " csv=" << csvPath
        << '\n';

    return EXIT_SUCCESS;
}
