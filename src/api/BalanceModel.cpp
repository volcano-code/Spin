#include "api/BalanceModel.h"

#include "deformation/DeformationSystem.h"
#include "deformation/DeformedVoxelMassField.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "physics/BalanceEnergy.h"
#include "stability/BalanceTargetSchedule.h"
#include "stability/SupportRegion.h"
#include "voxel/MaterialOccupancy.h"
#include "voxel/VoxelDomain.h"
#include "voxel/VoxelGrid.h"

#include <igl/read_triangle_mesh.h>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <fcntl.h>
#include <unistd.h>


namespace api
{
namespace
{


constexpr double HALF_PI =
    1.57079632679489661923132169163975144;


struct DiagnosticsStanding
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::MatrixXd polygon;
    Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -1.0, 0.0);
    double gammaRadians = 0.0;
};


struct DiagnosticsSuspension
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::Vector3d attachment = Eigen::Vector3d::Zero();
    double angleRadians = 0.0;
};


struct DiagnosticsHandle
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::Vector3d pivot = Eigen::Vector3d::Zero();
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
};


using DiagnosticsSuspensions = std::vector<
    DiagnosticsSuspension,
    Eigen::aligned_allocator<DiagnosticsSuspension>
>;


using DiagnosticsStandings = std::vector<
    DiagnosticsStanding,
    Eigen::aligned_allocator<DiagnosticsStanding>
>;


using DiagnosticsHandles = std::vector<
    DiagnosticsHandle,
    Eigen::aligned_allocator<DiagnosticsHandle>
>;


void appendNumber(std::ostream& output, double value);


struct PerformanceDurations
{
    double inputSnapshot = 0.0;
    double voxelization = 0.0;
    double domainAndHandles = 0.0;
    double deformationSystemAndBbw = 0.0;
    double targetAndOptimizerSetup = 0.0;
    double optimization = 0.0;
    double manufacturingAudit = 0.0;
    double diagnosticsPayloadSetup = 0.0;
    double diagnosticsSerialization = 0.0;
    double total = 0.0;
};


using PerformanceClock = std::chrono::steady_clock;


double elapsedSeconds(
    const PerformanceClock::time_point& start,
    const PerformanceClock::time_point& end
)
{
    return std::chrono::duration<double>(end - start).count();
}


std::string buildPerformanceReportJson(
    const PerformanceDurations& durations
)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "{\n"
           << "  \"schema\": \"make-it-stand-r10-performance-report-v1\",\n"
           << "  \"schema_version\": 1,\n"
           << "  \"clock\": \"steady_clock\",\n"
           << "  \"scope\": \"BalanceModel.optimize_after_config_validation_"
              "through_diagnostics_serialization\",\n"
           << "  \"phases\": {\n"
           << "    \"input_snapshot_seconds\": ";
    appendNumber(output, durations.inputSnapshot);
    output << ",\n    \"voxelization_seconds\": ";
    appendNumber(output, durations.voxelization);
    output << ",\n    \"domain_and_handles_seconds\": ";
    appendNumber(output, durations.domainAndHandles);
    output << ",\n    \"deformation_system_and_bbw_seconds\": ";
    appendNumber(output, durations.deformationSystemAndBbw);
    output << ",\n    \"target_and_optimizer_setup_seconds\": ";
    appendNumber(output, durations.targetAndOptimizerSetup);
    output << ",\n    \"optimization_seconds\": ";
    appendNumber(output, durations.optimization);
    output << ",\n    \"manufacturing_audit_seconds\": ";
    appendNumber(output, durations.manufacturingAudit);
    output << ",\n    \"diagnostics_payload_setup_seconds\": ";
    appendNumber(output, durations.diagnosticsPayloadSetup);
    output << ",\n    \"diagnostics_serialization_seconds\": ";
    appendNumber(output, durations.diagnosticsSerialization);
    output << "\n  },\n  \"total_seconds\": ";
    appendNumber(output, durations.total);
    output << "\n}\n";
    return output.str();
}


void validateModelConfig(const BalanceModelConfig& config)
{
    if (
        config.voxelResolution < 3 ||
        !std::isfinite(config.minimumWallThickness) ||
        config.minimumWallThickness < 0.0 ||
        config.bbw.maxIterations <= 0 ||
        !std::isfinite(config.bbw.tolerance) ||
        config.bbw.tolerance <= 0.0
    )
    {
        throw std::invalid_argument(
            "balance model configuration is invalid"
        );
    }
}


void validateRunConfig(const BalanceRunConfig& run)
{
    const optimization::AlternatingBalanceOptimizerConfig& config =
        run.optimizer;
    const deformation::PaperLineSearchConfig& line = config.lineSearch;
    if (
        !std::isfinite(run.initialMu) ||
        run.initialMu < 0.0 ||
        run.initialMu > 1.0 ||
        config.maximumCarvingBlocks <= 0 ||
        config.maximumDeformationIterationsPerBlock <= 0 ||
        config.minimumCarvingRemovalsToContinue < 0 ||
        !std::isfinite(config.gradientTolerance) ||
        config.gradientTolerance < 0.0 ||
        !std::isfinite(config.shapeScaleLambda) ||
        config.shapeScaleLambda < 0.0 ||
        !std::isfinite(config.density) ||
        config.density <= 0.0 ||
        !std::isfinite(config.minimumPositiveMass) ||
        config.minimumPositiveMass < 0.0 ||
        !std::isfinite(config.minimumPositiveCellMass) ||
        config.minimumPositiveCellMass < 0.0 ||
        !std::isfinite(line.initialStep) ||
        !std::isfinite(line.minimumStep) ||
        !std::isfinite(line.reductionRatio) ||
        !std::isfinite(line.requiredRelativeDecrease) ||
        !std::isfinite(line.muDecrement) ||
        line.initialStep <= 0.0 ||
        line.minimumStep <= 0.0 ||
        line.minimumStep > line.initialStep ||
        line.reductionRatio <= 0.0 ||
        line.reductionRatio >= 1.0 ||
        line.requiredRelativeDecrease < 0.0 ||
        line.requiredRelativeDecrease >= 1.0 ||
        line.muDecrement <= 0.0 ||
        line.muDecrement > 1.0 ||
        (!config.enableCarving && !config.enableDeformation)
    )
    {
        throw std::invalid_argument(
            "balance run configuration is invalid"
        );
    }
}


void validateSourceMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (
        vertices.rows() == 0 ||
        vertices.cols() != 3 ||
        !vertices.allFinite() ||
        faces.rows() == 0 ||
        faces.cols() != 3
    )
    {
        throw std::invalid_argument(
            "source mesh must contain finite N-by-3 vertices and M-by-3 faces"
        );
    }
    if (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows())
    {
        throw std::invalid_argument(
            "source mesh face indices are out of range"
        );
    }

    const manufacturing::BoundaryAudit audit =
        manufacturing::auditTriangleBoundary(
            vertices,
            faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (
        !audit.strictBoundaryValid() ||
        audit.faceComponentCount != 1 ||
        audit.unreferencedVertexCount != 0
    )
    {
        throw std::domain_error(
            "source mesh must be one fully referenced, outward, closed, "
            "manifold, non-self-intersecting triangle boundary"
        );
    }
}


void validateReleasedBetaSourceMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (
        vertices.rows() == 0 ||
        vertices.cols() != 3 ||
        !vertices.allFinite() ||
        faces.rows() == 0 ||
        faces.cols() != 3
    )
    {
        throw std::invalid_argument(
            "released-beta source must contain finite normalized vertices and triangles"
        );
    }
    if (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows())
    {
        throw std::invalid_argument(
            "released-beta source face indices are out of range"
        );
    }
    constexpr double unitCubeTolerance = 2.0e-10;
    if (
        (vertices.array() < -unitCubeTolerance).any() ||
        (vertices.array() > 1.0 + unitCubeTolerance).any()
    )
    {
        throw std::invalid_argument(
            "released-beta source vertices must use normalized [0,1]^3 coordinates"
        );
    }

    /*
     * Boundary predicates now treat indexed shared-edge/shared-vertex contact
     * explicitly, so no enlarged compatibility tolerance is necessary.  The
     * normal scale-derived tolerance accepts the archived decimal sphere
     * facets while retaining the horse preset's genuine non-adjacent crossing.
     */
    const manufacturing::BoundaryAudit audit =
        manufacturing::auditTriangleBoundary(
            vertices,
            faces,
            manufacturing::BoundaryOrientation::OutwardPositive,
            false,
            0.0
        );
    if (
        !audit.strictBoundaryValid() ||
        audit.faceComponentCount != 1 ||
        audit.unreferencedVertexCount != 0
    )
    {
        throw std::domain_error(
            "released-beta source failed the compatibility boundary gate; "
            "true or non-adjacent intersections are never waived"
        );
    }
}


Eigen::MatrixXd normalizeAndCoordinateMapReleasedBetaWeights(
    const ReleasedBetaArchiveData& archive
)
{
    const int resolution = archive.voxelResolution;
    if (resolution < 3 || archive.hullDepth < 0)
    {
        throw std::invalid_argument(
            "released-beta archive has an invalid resolution or hull depth"
        );
    }
    const long long maximumIndex =
        static_cast<long long>(std::numeric_limits<int>::max());
    const auto checkedCubeCount = [maximumIndex](
        long long side,
        const char* label
    )
    {
        if (
            side <= 0LL ||
            side > maximumIndex / side ||
            side * side > maximumIndex / side
        )
        {
            throw std::overflow_error(
                std::string("released-beta ") + label +
                " exceeds integer indexing"
            );
        }
        return side * side * side;
    };
    const long long cellCount64 = checkedCubeCount(
        static_cast<long long>(resolution),
        "cell grid"
    );
    const long long nodeSide64 = static_cast<long long>(resolution) + 1LL;
    const long long nodeCount64 = checkedCubeCount(nodeSide64, "node grid");
    if (
        archive.cellIdsZFast.size() !=
            static_cast<std::size_t>(cellCount64)
    )
    {
        throw std::invalid_argument(
            "released-beta .vox cache must contain exactly N^3 entries"
        );
    }

    std::size_t occupiedCount = 0U;
    for (const int id : archive.cellIdsZFast)
    {
        if (id < -1)
        {
            throw std::invalid_argument(
                "released-beta .vox entries must be -1 or nonnegative IDs"
            );
        }
        occupiedCount += id >= 0 ? 1U : 0U;
    }
    if (occupiedCount == 0U)
    {
        throw std::invalid_argument(
            "released-beta .vox cache contains no occupied cells"
        );
    }
    std::vector<unsigned char> seen(occupiedCount, 0U);
    for (const int id : archive.cellIdsZFast)
    {
        if (id < 0)
        {
            continue;
        }
        if (
            static_cast<std::size_t>(id) >= occupiedCount ||
            seen[static_cast<std::size_t>(id)] != 0U
        )
        {
            throw std::invalid_argument(
                "released-beta occupied IDs must be unique and contiguous"
            );
        }
        seen[static_cast<std::size_t>(id)] = 1U;
    }

    const int nodeSide = resolution + 1;
    const auto betaCellIndex = [resolution](int x, int y, int z)
    {
        return static_cast<std::size_t>(
            (x * resolution + y) * resolution + z
        );
    };
    const auto betaNodeCoordinateIndex = [nodeSide](int x, int y, int z)
    {
        return static_cast<std::size_t>(
            (x * nodeSide + y) * nodeSide + z
        );
    };
    const auto currentGlobalCorner = [nodeSide](int x, int y, int z)
    {
        return x + nodeSide * (y + nodeSide * z);
    };

    std::vector<unsigned char> active(
        static_cast<std::size_t>(nodeCount64),
        0U
    );
    for (int x = 0; x < resolution; ++x)
    {
        for (int y = 0; y < resolution; ++y)
        {
            for (int z = 0; z < resolution; ++z)
            {
                if (archive.cellIdsZFast[betaCellIndex(x, y, z)] < 0)
                {
                    continue;
                }
                for (int dx = 0; dx <= 1; ++dx)
                {
                    for (int dy = 0; dy <= 1; ++dy)
                    {
                        for (int dz = 0; dz <= 1; ++dz)
                        {
                            active[betaNodeCoordinateIndex(
                                x + dx,
                                y + dy,
                                z + dz
                            )] = 1U;
                        }
                    }
                }
            }
        }
    }

    const Eigen::Index activeNodeCount = static_cast<Eigen::Index>(
        std::count(active.begin(), active.end(), 1U)
    );
    const Eigen::MatrixXd& raw = archive.betaNodeWeightsZFast;
    if (
        raw.rows() != activeNodeCount ||
        raw.cols() <= 0 ||
        !raw.allFinite()
    )
    {
        throw std::invalid_argument(
            "released-beta .bbw shape must be active-node by handle and finite"
        );
    }

    std::vector<int> betaRowByCoordinate(
        static_cast<std::size_t>(nodeCount64),
        -1
    );
    int betaRow = 0;
    for (int x = 0; x < nodeSide; ++x)
    {
        for (int y = 0; y < nodeSide; ++y)
        {
            for (int z = 0; z < nodeSide; ++z)
            {
                const std::size_t coordinate =
                    betaNodeCoordinateIndex(x, y, z);
                if (active[coordinate] != 0U)
                {
                    betaRowByCoordinate[coordinate] = betaRow++;
                }
            }
        }
    }

    Eigen::MatrixXd betaNormalized(raw.rows(), raw.cols());
    for (Eigen::Index row = 0; row < raw.rows(); ++row)
    {
        double sum = 0.0;
        for (Eigen::Index column = 0; column < raw.cols(); ++column)
        {
            const double value = raw(row, column);
            if (value < 0.0)
            {
                throw std::invalid_argument(
                    "released-beta .bbw weights must be nonnegative"
                );
            }
            sum += value;
        }
        if (!std::isfinite(sum) || !(sum > 0.0))
        {
            throw std::invalid_argument(
                "released-beta .bbw rows must have a positive finite sum"
            );
        }
        betaNormalized.row(row) = raw.row(row) / sum;
    }

    Eigen::MatrixXd mapped(activeNodeCount, raw.cols());
    Eigen::Index destination = 0;
    for (int z = 0; z < nodeSide; ++z)
    {
        for (int y = 0; y < nodeSide; ++y)
        {
            for (int x = 0; x < nodeSide; ++x)
            {
                const std::size_t coordinate =
                    betaNodeCoordinateIndex(x, y, z);
                if (active[coordinate] == 0U)
                {
                    continue;
                }
                const int source = betaRowByCoordinate[coordinate];
                if (source < 0 || currentGlobalCorner(x, y, z) < 0)
                {
                    throw std::logic_error(
                        "released-beta active-node coordinate map is inconsistent"
                    );
                }
                mapped.row(destination++) = betaNormalized.row(source);
            }
        }
    }
    if (destination != activeNodeCount)
    {
        throw std::logic_error(
            "released-beta active-node coordinate map has the wrong size"
        );
    }
    return mapped;
}


int releasedBetaClosestActiveNode(
    const deformation::GridDeformationDomain& domain,
    const Eigen::Vector3d& point
)
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "released-beta handle point must be finite"
        );
    }
    const Eigen::Vector3i dimensions = domain.dimensions();
    if (
        dimensions(0) != dimensions(1) ||
        dimensions(1) != dimensions(2)
    )
    {
        throw std::logic_error(
            "released-beta handle binding requires a cubic lattice"
        );
    }
    const int resolution = dimensions(0);
    Eigen::Vector3i coordinate;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double scaled = point(axis) * static_cast<double>(resolution);
        const int cell = static_cast<int>(std::floor(scaled));
        if (cell < 0 || cell >= resolution)
        {
            throw std::domain_error(
                "released-beta handle point lies outside the voxel domain"
            );
        }
        const double fraction = scaled - static_cast<double>(cell);
        coordinate(axis) = cell + (fraction <= 0.5 ? 0 : 1);
    }

    const int side = resolution + 1;
    const int global =
        coordinate(0) + side * (coordinate(1) + side * coordinate(2));
    const int local = domain.globalToLocalCorner(global);
    if (local < 0)
    {
        throw std::domain_error(
            "released-beta handle maps to an inactive lattice node"
        );
    }
    return local;
}


bool isProperRotation(const Eigen::Matrix3d& rotation)
{
    constexpr double tolerance = 1.0e-10;
    return
        rotation.allFinite() &&
        (rotation.transpose() * rotation).isApprox(
            Eigen::Matrix3d::Identity(),
            tolerance
        ) &&
        std::abs(rotation.determinant() - 1.0) <= tolerance;
}


double meshScale(const Eigen::MatrixXd& vertices)
{
    return (
        vertices.colwise().maxCoeff() -
        vertices.colwise().minCoeff()
    ).norm();
}


double matchingTolerance(const Eigen::MatrixXd& vertices)
{
    const double extent = meshScale(vertices);
    const double coordinate = vertices.cwiseAbs().maxCoeff();
    return 4096.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, extent, coordinate});
}


struct StandingConstruction
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    stability::SupportRegion support;
    double gammaRadians = 0.0;
    Eigen::MatrixXd fixedPoints;
};


StandingConstruction constructStandingFromPolygon(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& polygon,
    const Eigen::Vector3d& gravityDirection,
    double gammaRadians
)
{
    const Eigen::Vector3d gravity =
        normalizedGravityDirection(gravityDirection);
    const stability::SupportRegion support =
        stability::SupportRegion::fromPolygon(polygon, gravity);
    static_cast<void>(stability::StandingTargetSpecification(
        support,
        gammaRadians
    ));

    const double tolerance = matchingTolerance(vertices);
    Eigen::VectorXd sourceHeights(vertices.rows());
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        sourceHeights(row) = -gravity.dot(vertices.row(row).transpose());
    }
    const double trueMinimum = sourceHeights.minCoeff();
    if (std::abs(support.supportHeight() - trueMinimum) > tolerance)
    {
        throw std::domain_error(
            "support polygon plane does not coincide with the true mesh floor"
        );
    }
    if ((sourceHeights.array() < support.supportHeight() - tolerance).any())
    {
        throw std::domain_error(
            "a source-mesh vertex lies below the support polygon plane"
        );
    }
    for (Eigen::Index point = 0; point < polygon.rows(); ++point)
    {
        bool matchedLowestVertex = false;
        for (Eigen::Index vertex = 0; vertex < vertices.rows(); ++vertex)
        {
            if (
                std::abs(sourceHeights(vertex) - trueMinimum) <= tolerance &&
                (
                    vertices.row(vertex).transpose() -
                    polygon.row(point).transpose()
                ).norm() <= tolerance
            )
            {
                matchedLowestVertex = true;
                break;
            }
        }
        if (!matchedLowestVertex)
        {
            throw std::domain_error(
                "every support polygon point must match a lowest source-mesh vertex"
            );
        }
    }

    Eigen::MatrixXd fixedPoints(support.supportPolygon2D().rows(), 3);
    for (Eigen::Index row = 0; row < fixedPoints.rows(); ++row)
    {
        fixedPoints.row(row) = support.fromPlaneCoordinates(
            support.supportPolygon2D().row(row).transpose()
        ).transpose();
    }
    return StandingConstruction{
        support,
        gammaRadians,
        std::move(fixedPoints)
    };
}


StandingConstruction constructStandingFromMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::Vector3d& gravityDirection,
    double contactThreshold,
    double gammaRadians,
    double minimumStabilityInset = 0.0,
    bool supportPlaneAtBandTop = false
)
{
    const Eigen::Vector3d gravity =
        normalizedGravityDirection(gravityDirection);
    const stability::SupportRegion support = stability::SupportRegion::fromMesh(
        vertices,
        gravity,
        contactThreshold,
        minimumStabilityInset,
        supportPlaneAtBandTop
    );
    static_cast<void>(stability::StandingTargetSpecification(
        support,
        gammaRadians
    ));
    Eigen::MatrixXd fixedPoints(
        static_cast<Eigen::Index>(support.contactVertexIndices().size()),
        3
    );
    for (std::size_t index = 0;
         index < support.contactVertexIndices().size();
         ++index)
    {
        fixedPoints.row(static_cast<Eigen::Index>(index)) =
            vertices.row(support.contactVertexIndices()[index]);
    }
    return StandingConstruction{
        support,
        gammaRadians,
        std::move(fixedPoints)
    };
}


std::string statusToken(
    optimization::AlternatingOptimizationStatus status
)
{
    using Status = optimization::AlternatingOptimizationStatus;
    switch (status)
    {
        case Status::StableAfterCarving:
            return "stable_after_carving";
        case Status::StableAfterDeformation:
            return "stable_after_deformation";
        case Status::StableWithoutCarving:
            return "stable_without_carving";
        case Status::NoDeformationDegreesOfFreedom:
            return "no_deformation_degrees_of_freedom";
        case Status::DeformationDisabledUnstable:
            return "deformation_disabled_unstable";
        case Status::DeformationConvergedUnstable:
            return "deformation_converged_unstable";
        case Status::MuExhaustedUnstable:
            return "mu_exhausted_unstable";
        case Status::AlternatingStalled:
            return "alternating_stalled";
        case Status::MaximumCarvingBlocks:
            return "maximum_carving_blocks";
        case Status::MaximumCarvingIterations:
            return "maximum_carving_iterations";
        case Status::MaximumDeformationIterations:
            return "maximum_deformation_iterations";
    }
    return "unknown";
}


bool optimizationSucceeded(
    const optimization::AlternatingBalanceOptimizationResult& result
)
{
    using Status = optimization::AlternatingOptimizationStatus;
    return
        result.targets.allReached &&
        (
            result.status == Status::StableAfterCarving ||
            result.status == Status::StableAfterDeformation ||
            result.status == Status::StableWithoutCarving
        );
}


const char* eventKindToken(
    optimization::AlternatingEventKind kind
)
{
    using Kind = optimization::AlternatingEventKind;
    switch (kind)
    {
        case Kind::Initial:
            return "initial";
        case Kind::SolidRefill:
            return "solid_refill";
        case Kind::CarvingIteration:
            return "carving_iteration";
        case Kind::AcceptedDeformation:
            return "accepted_deformation";
    }
    return "unknown";
}


const char* standingStateToken(stability::StandingState state)
{
    switch (state)
    {
        case stability::StandingState::OutsideSupport:
            return "outside_support";
        case stability::StandingState::FeasibleButUnsafe:
            return "feasible_but_unsafe";
        case stability::StandingState::Stable:
            return "stable";
    }
    return "unknown";
}


const char* trialDecisionToken(deformation::PaperTrialDecision decision)
{
    switch (decision)
    {
        case deformation::PaperTrialDecision::Accepted:
            return "accepted";
        case deformation::PaperTrialDecision::StepReduced:
            return "step_reduced";
        case deformation::PaperTrialDecision::MuReduced:
            return "mu_reduced";
        case deformation::PaperTrialDecision::MuExhausted:
            return "mu_exhausted";
    }
    return "unknown";
}


std::string jsonBoolean(bool value)
{
    return value ? "true" : "false";
}


void appendNumber(std::ostream& output, double value)
{
    if (std::isfinite(value))
    {
        output << std::setprecision(17) << value;
    }
    else
    {
        output << "null";
    }
}


void appendVector3(
    std::ostream& output,
    const Eigen::Vector3d& value
)
{
    output << '[';
    appendNumber(output, value.x());
    output << ',';
    appendNumber(output, value.y());
    output << ',';
    appendNumber(output, value.z());
    output << ']';
}


void appendVector(
    std::ostream& output,
    const Eigen::VectorXd& value
)
{
    output << '[';
    for (Eigen::Index index = 0; index < value.size(); ++index)
    {
        if (index != 0)
        {
            output << ',';
        }
        appendNumber(output, value(index));
    }
    output << ']';
}


void appendPointMatrix(
    std::ostream& output,
    const Eigen::MatrixXd& points
)
{
    output << '[';
    for (Eigen::Index row = 0; row < points.rows(); ++row)
    {
        if (row != 0)
        {
            output << ',';
        }
        appendVector3(output, points.row(row).transpose());
    }
    output << ']';
}


void appendMatrix3(
    std::ostream& output,
    const Eigen::Matrix3d& value
)
{
    output << '[';
    for (Eigen::Index row = 0; row < 3; ++row)
    {
        if (row != 0)
        {
            output << ',';
        }
        appendVector3(output, value.row(row).transpose());
    }
    output << ']';
}


void appendMassProperties(
    std::ostream& output,
    const MassProperties& properties
)
{
    output << "{\"mass\":";
    appendNumber(output, properties.mass);
    output << ",\"first_moment\":";
    appendVector3(output, properties.firstMoment);
    output << '}';
}


void appendOptionalBool(
    std::ostream& output,
    const std::optional<bool>& value
)
{
    output << (value.has_value() ? jsonBoolean(*value) : "null");
}


void appendOptionalNumber(
    std::ostream& output,
    const std::optional<double>& value
)
{
    if (value.has_value())
    {
        appendNumber(output, *value);
    }
    else
    {
        output << "null";
    }
}


std::string buildDiagnosticsJson(
    const BalanceModelConfig& modelConfig,
    const BalanceRunConfig& runConfig,
    bool releasedBetaArchiveImport,
    int sourceVertexCount,
    int sourceFaceCount,
    const Eigen::Vector3d& gravity,
    const DiagnosticsStandings& standings,
    const DiagnosticsSuspensions& suspensions,
    const DiagnosticsHandles& handles,
    const optimization::AlternatingBalanceOptimizationResult& result,
    const manufacturing::ManufacturingGeometry& manufacturing
)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::boolalpha;
    const auto& config = runConfig.optimizer;
    const auto& line = config.lineSearch;

    output << "{\n";
    output << "  \"schema\": \"make-it-stand-r9-balance-result-v1\",\n";
    output << "  \"schema_version\": 1,\n";
    output << "  \"status\": \"" << statusToken(result.status) << "\",\n";
    output << "  \"succeeded\": "
           << jsonBoolean(optimizationSucceeded(result)) << ",\n";
    output << "  \"source_mesh\": {\"vertices\":" << sourceVertexCount
           << ",\"faces\":" << sourceFaceCount << "},\n";
    output << "  \"configuration\": {\n";
    output << "    \"voxel_resolution\": " << modelConfig.voxelResolution << ",\n";
    output << "    \"fixed_unit_cube_voxel_grid\": "
           << jsonBoolean(modelConfig.fixedUnitCubeVoxelGrid) << ",\n";
    if (releasedBetaArchiveImport)
    {
        output << "    \"released_beta_archive_import\": true,\n";
    }
    output << "    \"minimum_wall_thickness\": ";
    appendNumber(output, modelConfig.minimumWallThickness);
    output << ",\n    \"bbw_max_iterations\": "
           << modelConfig.bbw.maxIterations << ",\n";
    output << "    \"bbw_tolerance\": ";
    appendNumber(output, modelConfig.bbw.tolerance);
    output << ",\n    \"initial_mu\": ";
    appendNumber(output, runConfig.initialMu);
    output << ",\n    \"maximum_carving_blocks\": "
           << config.maximumCarvingBlocks << ",\n";
    output << "    \"maximum_deformation_iterations_per_block\": "
           << config.maximumDeformationIterationsPerBlock << ",\n";
    output << "    \"gradient_tolerance\": ";
    appendNumber(output, config.gradientTolerance);
    output << ",\n    \"shape_scale_lambda\": ";
    appendNumber(output, config.shapeScaleLambda);
    output << ",\n    \"effective_shape_scale_lambda\": ";
    appendNumber(output, result.effectiveShapeScaleLambda);
    output << ",\n    \"enable_carving\": "
           << jsonBoolean(config.enableCarving);
    output << ",\n    \"enable_deformation\": "
           << jsonBoolean(config.enableDeformation);
    output << ",\n    \"enable_laplacian\": "
           << jsonBoolean(config.enableLaplacian);
    output << ",\n    \"freeze_targets_after_initial_refresh\": "
           << jsonBoolean(config.freezeTargetsAfterInitialRefresh);
    output << ",\n    \"released_beta_outer_loop\": "
           << jsonBoolean(config.releasedBetaOuterLoop);
    output << ",\n    \"density\": ";
    appendNumber(output, config.density);
    output << ",\n    \"minimum_positive_mass\": ";
    appendNumber(output, config.minimumPositiveMass);
    output << ",\n    \"minimum_positive_cell_mass\": ";
    appendNumber(output, config.minimumPositiveCellMass);
    output << ",\n    \"minimum_carving_removals_to_continue\": "
           << config.minimumCarvingRemovalsToContinue << ",\n";
    output << "    \"line_search\": {\"initial_step\":";
    appendNumber(output, line.initialStep);
    output << ",\"minimum_step\":";
    appendNumber(output, line.minimumStep);
    output << ",\"reduction_ratio\":";
    appendNumber(output, line.reductionRatio);
    output << ",\"required_relative_decrease\":";
    appendNumber(output, line.requiredRelativeDecrease);
    output << ",\"mu_decrement\":";
    appendNumber(output, line.muDecrement);
    output << ",\"released_beta_semantics\":"
           << jsonBoolean(line.releasedBetaSemantics);
    output << "}\n  },\n";
    output << "  \"problem\": {\"gravity\":";
    appendVector3(output, gravity);
    output << ",\"gravity_scope\":\"legacy_default_for_support_and_suspension_constructors\"";
    output << ",\"standing_constraint\":";
    if (standings.size() == 1U)
    {
        output << "{\"gamma_radians\":";
        appendNumber(output, standings.front().gammaRadians);
        output << ",\"gravity\":";
        appendVector3(output, standings.front().gravity);
        output << ",\"support_polygon\":";
        appendPointMatrix(output, standings.front().polygon);
        output << '}';
    }
    else
    {
        output << "null";
    }
    output << ",\"standing_constraints\":[";
    for (std::size_t index = 0; index < standings.size(); ++index)
    {
        if (index != 0U) output << ',';
        output << "{\"gravity\":";
        appendVector3(output, standings[index].gravity);
        output << ",\"gamma_radians\":";
        appendNumber(output, standings[index].gammaRadians);
        output << ",\"support_polygon\":";
        appendPointMatrix(output, standings[index].polygon);
        output << '}';
    }
    output << "],\"suspension_constraints\":[";
    for (std::size_t index = 0; index < suspensions.size(); ++index)
    {
        if (index != 0U) output << ',';
        output << "{\"attachment\":";
        appendVector3(output, suspensions[index].attachment);
        output << ",\"allowed_deviation_radians\":";
        appendNumber(output, suspensions[index].angleRadians);
        output << '}';
    }
    output << "],\"user_handles\":[";
    for (std::size_t index = 0; index < handles.size(); ++index)
    {
        if (index != 0U) output << ',';
        output << "{\"rest_pivot\":";
        appendVector3(output, handles[index].pivot);
        output << ",\"fixed_rotation\":";
        appendMatrix3(output, handles[index].rotation);
        output << '}';
    }
    output << "]},\n";
    output << "  \"r10_metrics\": {\n";
    output << "    \"semantics\": "
           << "\"paper_profile_committed_events_not_beta2013_iterations\",\n";
    output << "    \"grid_dimensions\": ["
           << result.gridDimensions(0) << ','
           << result.gridDimensions(1) << ','
           << result.gridDimensions(2) << "],\n";
    output << "    \"voxel_size\": ";
    appendVector3(output, result.voxelSize);
    output << ",\n    \"initial_material_voxels\": "
           << result.initialMaterialVoxelCount;
    output << ",\n    \"initial_carveable_voxels\": "
           << result.initialCarveableVoxelCount;
    output << ",\n    \"active_deformation_nodes\": "
           << result.activeDeformationNodeCount;
    output << ",\n    \"fixed_deformation_nodes\": "
           << result.fixedDeformationNodeCount;
    output << ",\n    \"fixed_deformation_node_fraction\": ";
    if (result.activeDeformationNodeCount > 0)
    {
        appendNumber(
            output,
            static_cast<double>(result.fixedDeformationNodeCount) /
                static_cast<double>(result.activeDeformationNodeCount)
        );
    }
    else
    {
        output << "null";
    }
    output << ",\n    \"user_handles\": " << result.userHandleCount;
    output << ",\n    \"initial_mass_properties\": ";
    appendMassProperties(output, result.initialMassProperties);
    output << ",\n    \"initial_center_of_mass\": ";
    appendVector3(output, result.initialCenterOfMass);
    output << ",\n    \"final_projected_target_energy\": ";
    appendNumber(output, result.finalProjectedTargetEnergy);
    output << ",\n    \"final_unscaled_shape_energy\": ";
    appendNumber(output, result.finalUnscaledShapeEnergy);
    output << ",\n    \"final_shape_energy\": ";
    appendNumber(output, result.finalShapeEnergy);
    output << ",\n    \"final_total_energy\": ";
    appendNumber(output, result.finalTotalEnergy);
    output << ",\n    \"mass_moment_source\": "
              "\"outer_surface_minus_deformed_carved_cells\"";
    output << ",\n    \"final_mass_moment_mass\": ";
    appendNumber(output, result.finalMassMoments.mass);
    output << ",\n    \"final_mass_moment_first\": ";
    appendVector3(output, result.finalMassMoments.firstMoment);
    output << ",\n    \"final_raw_second_moment\": ";
    appendMatrix3(output, result.finalMassMoments.secondMoment);
    output << ",\n    \"final_inertia_tensor_about_center\": ";
    appendMatrix3(output, result.finalInertiaTensorAboutCenter);
    output << ",\n    \"normalized_rms_vertex_displacement\": ";
    appendNumber(output, result.normalizedRmsVertexDisplacement);
    output << ",\n    \"normalized_maximum_vertex_displacement\": ";
    appendNumber(output, result.normalizedMaximumVertexDisplacement);
    output << ",\n    \"maximum_user_handle_translation\": ";
    appendNumber(output, result.maximumUserHandleTranslation);
    output << ",\n    \"maximum_user_handle_scale_deviation\": ";
    appendNumber(output, result.maximumUserHandleScaleDeviation);
    output << ",\n    \"initial_outer_volume\": ";
    appendNumber(output, result.initialOuterVolume);
    output << ",\n    \"final_outer_volume\": ";
    appendNumber(output, result.finalOuterVolume);
    output << ",\n    \"final_cavity_volume\": ";
    appendNumber(output, result.finalCavityVolume);
    output << ",\n    \"final_cavity_volume_fraction\": ";
    appendNumber(output, result.finalCavityVolumeFraction);
    output << ",\n    \"material_mass_change_fraction\": ";
    appendNumber(output, result.finalMaterialMassChangeFraction);
    output << ",\n    \"first_standing_feasible_event\": "
           << result.firstStandingFeasibleEvent;
    output << ",\n    \"first_safe_event\": " << result.firstSafeEvent;
    output << ",\n    \"committed_event_count\": "
           << result.events.size() << "\n  },\n";
    output << "  \"geometry_hash\": \"" << result.geometryHash << "\",\n";
    output << "  \"geometry_epoch\": " << result.geometryEpoch << ",\n";
    output << "  \"mass_properties\": ";
    appendMassProperties(output, result.massProperties);
    output << ",\n  \"center_of_mass\": ";
    appendVector3(output, result.centerOfMass);
    output << ",\n  \"final_mu\": ";
    appendNumber(output, result.mu);
    output << ",\n  \"parameters\": ";
    appendVector(output, result.parameters);
    output << ",\n  \"targets\": {\"generation\":"
           << result.targets.generation << ",\"all_reached\":"
           << jsonBoolean(result.targets.allReached)
           << ",\"constraints\":[";
    for (std::size_t index = 0;
         index < result.targets.constraints.size();
         ++index)
    {
        if (index != 0U)
        {
            output << ',';
        }
        const stability::BalanceConstraintSnapshot& target =
            result.targets.constraints[index];
        output << "{\"type\":\""
               << (target.type == stability::BalanceConstraintType::Standing
                       ? "standing" : "suspension")
               << "\",\"reached\":" << jsonBoolean(target.reached)
               << ",\"standing_state\":";
        if (target.standingState.has_value())
        {
            output << '\"' << standingStateToken(*target.standingState)
                   << '\"';
        }
        else
        {
            output << "null";
        }
        output << ",\"requested_margin_attainable\":";
        appendOptionalBool(output, target.requestedMarginAttainable);
        output << ",\"suspension_deviation_radians\":";
        appendOptionalNumber(output, target.suspensionDeviationRadians);
        output << ",\"suspension_downstream\":";
        appendOptionalBool(output, target.suspensionDownstream);
        output << '}';
    }
    output << "]},\n";

    output << "  \"carving_blocks\": [";
    for (std::size_t index = 0; index < result.carvingBlocks.size(); ++index)
    {
        if (index != 0U) output << ',';
        const auto& record = result.carvingBlocks[index];
        output << "{\"block\":" << record.block
               << ",\"geometry_epoch_before_reset\":"
               << record.geometryEpochBeforeReset
               << ",\"geometry_epoch_after_reset\":"
               << record.geometryEpochAfterReset
               << ",\"geometry_epoch\":" << record.geometryEpoch
               << ",\"refilled_voxels\":" << record.refilledVoxelCount
               << ",\"alpha_hamming_distance\":"
               << record.alphaHammingDistance
               << ",\"iterations\":" << record.carvingIterations
               << ",\"carved_voxels\":" << record.carvedVoxelCount
               << ",\"stop_reason\":\""
               << carvingStopReasonName(record.stopReason) << "\""
               << ",\"reset_target_energy\":";
        appendNumber(output, record.resetTargetEnergy);
        output << ",\"final_target_energy\":";
        appendNumber(output, record.finalTargetEnergy);
        output << ",\"target_generation_before\":"
               << record.targetGenerationBefore
               << ",\"target_generation_after\":"
               << record.targetGenerationAfter
               << ",\"objective_context_changed\":"
               << jsonBoolean(record.objectiveContextChanged)
               << ",\"all_targets_reached\":"
               << jsonBoolean(record.allTargetsReached)
               << ",\"alpha_hash\":\"" << record.alphaHash
               << "\",\"geometry_hash\":\"" << record.geometryHash
               << "\",\"mass_properties\":";
        appendMassProperties(output, record.massProperties);
        output << ",\"center_of_mass\":";
        appendVector3(output, record.centerOfMass);
        output << '}';
    }
    output << "],\n";

    output << "  \"carving_iterations\": [";
    for (std::size_t index = 0; index < result.carvingIterations.size(); ++index)
    {
        if (index != 0U) output << ',';
        const auto& record = result.carvingIterations[index];
        output << "{\"block\":" << record.block
               << ",\"iteration\":" << record.iteration
               << ",\"geometry_epoch_before\":"
               << record.geometryEpochBefore
               << ",\"geometry_epoch_after\":"
               << record.geometryEpochAfter
               << ",\"carved_voxels\":" << record.carvedVoxelCount
               << ",\"total_carved_voxels\":"
               << record.totalCarvedVoxelCount
               << ",\"stop_reason\":\""
               << carvingStopReasonName(record.stopReason) << "\""
               << ",\"frozen_energy_before\":";
        appendNumber(output, record.frozenTargetEnergyBefore);
        output << ",\"frozen_energy_after\":";
        appendNumber(output, record.frozenTargetEnergyAfter);
        output << ",\"refreshed_energy_after\":";
        appendNumber(output, record.refreshedTargetEnergyAfter);
        output << ",\"target_generation_before\":"
               << record.targetGenerationBefore
               << ",\"target_generation_after\":"
               << record.targetGenerationAfter
               << ",\"objective_context_changed\":"
               << jsonBoolean(record.objectiveContextChanged)
               << ",\"all_targets_reached\":"
               << jsonBoolean(record.allTargetsReached)
               << ",\"alpha_hash_before\":\""
               << record.alphaHashBefore << "\",\"alpha_hash_after\":\""
               << record.alphaHashAfter << "\",\"mass_properties\":";
        appendMassProperties(output, record.massProperties);
        output << ",\"center_of_mass\":";
        appendVector3(output, record.centerOfMass);
        output << '}';
    }
    output << "],\n";

    output << "  \"deformation_iterations\": [";
    for (std::size_t index = 0;
         index < result.deformationIterations.size();
         ++index)
    {
        if (index != 0U) output << ',';
        const auto& record = result.deformationIterations[index];
        output << "{\"block\":" << record.block
               << ",\"iteration\":" << record.iteration
               << ",\"decision\":\""
               << trialDecisionToken(record.decision) << "\""
               << ",\"accepted\":" << jsonBoolean(record.accepted)
               << ",\"rolled_back\":" << jsonBoolean(record.rolledBack)
               << ",\"energy_before\":";
        appendNumber(output, record.energyBefore);
        output << ",\"trial_energy\":";
        appendNumber(output, record.trialEnergy);
        output << ",\"energy_after\":";
        appendNumber(output, record.energyAfter);
        output << ",\"refreshed_energy_after\":";
        appendNumber(output, record.refreshedTargetEnergyAfter);
        output << ",\"gradient_norm\":";
        appendNumber(output, record.gradientNorm);
        output << ",\"attempted_step\":";
        appendNumber(output, record.attemptedStep);
        output << ",\"mu_before\":";
        appendNumber(output, record.muBefore);
        output << ",\"mu_after\":";
        appendNumber(output, record.stateAfter.mu);
        output << ",\"step_after\":";
        appendNumber(output, record.stateAfter.step);
        output << ",\"objective_reweighted\":"
               << jsonBoolean(record.objectiveReweighted)
               << ",\"released_beta_outer_loop\":"
               << jsonBoolean(record.releasedBetaOuterLoop)
               << ",\"refill_calls\":" << record.refillCalls
               << ",\"carving_calls\":" << record.carvingCalls
               << ",\"trial_balance_calls\":"
               << record.trialBalanceCalls
               << ",\"rollback_balance_calls\":"
               << record.rollbackBalanceCalls
               << ",\"trial_carving_iterations\":"
               << record.trialCarvingIterations
               << ",\"trial_carved_voxels\":"
               << record.trialCarvedVoxelCount
               << ",\"rollback_carving_iterations\":"
               << record.rollbackCarvingIterations
               << ",\"rollback_carved_voxels\":"
               << record.rollbackCarvedVoxelCount
               << ",\"geometry_epoch_before\":"
               << record.geometryEpochBefore
               << ",\"geometry_epoch_after\":"
               << record.geometryEpochAfter
               << ",\"target_generation_before\":"
               << record.targetGenerationBefore
               << ",\"target_generation_after\":"
               << record.targetGenerationAfter
               << ",\"objective_context_changed\":"
               << jsonBoolean(record.objectiveContextChanged)
               << ",\"all_targets_reached\":"
               << jsonBoolean(record.allTargetsReached)
               << ",\"alpha_hash_before\":\""
               << record.alphaHashBefore << "\",\"alpha_hash_after\":\""
               << record.alphaHashAfter << "\",\"trial_alpha_hash\":\""
               << record.trialAlphaHash << "\",\"rollback_alpha_hash\":\""
               << record.rollbackAlphaHash
               << "\",\"geometry_hash_before\":\""
               << record.geometryHashBefore
               << "\",\"geometry_hash_after\":\""
               << record.geometryHashAfter << "\",\"parameters_before\":";
        appendVector(output, record.parametersBefore);
        output << ",\"parameters_after\":";
        appendVector(output, record.parametersAfter);
        output << ",\"mass_properties_after\":";
        appendMassProperties(output, record.massPropertiesAfter);
        output << ",\"center_of_mass_after\":";
        appendVector3(output, record.centerOfMassAfter);
        output << '}';
    }
    output << "],\n";
    output << "  \"r10_events\": [";
    for (std::size_t index = 0; index < result.events.size(); ++index)
    {
        if (index != 0U) output << ',';
        const auto& record = result.events[index];
        output << "{\"ordinal\":" << record.ordinal
               << ",\"kind\":\"" << eventKindToken(record.kind) << "\""
               << ",\"block\":" << record.block
               << ",\"iteration\":" << record.iteration
               << ",\"geometry_epoch\":" << record.geometryEpoch
               << ",\"geometry_changed\":"
               << jsonBoolean(record.geometryChanged)
               << ",\"has_standing_constraints\":"
               << jsonBoolean(record.hasStandingConstraints)
               << ",\"all_standing_feasible\":"
               << jsonBoolean(record.allStandingFeasible)
               << ",\"all_targets_reached\":"
               << jsonBoolean(record.allTargetsReached)
               << ",\"target_generation\":" << record.targetGeneration
               << ",\"target_snapshot_hash\":\""
               << record.targetSnapshotHash << "\""
               << ",\"alpha_hash\":\"" << record.alphaHash << "\""
               << ",\"geometry_hash\":\"" << record.geometryHash << "\""
               << ",\"parameters\":";
        appendVector(output, record.parameters);
        output
               << ",\"projected_target_energy\":";
        appendNumber(output, record.projectedTargetEnergy);
        output << ",\"mass_properties\":";
        appendMassProperties(output, record.massProperties);
        output << ",\"center_of_mass\":";
        appendVector3(output, record.centerOfMass);
        output << '}';
    }
    output << "],\n";
    output << "  \"manufacturing\": "
           << manufacturing.report.toJson() << '\n';
    output << "}\n";
    return output.str();
}


std::string systemErrorMessage(
    const std::string& operation,
    int error
)
{
    return operation + ": " + std::strerror(error);
}


void writeAtomicBytes(
    const std::string& filename,
    const std::string& bytes
)
{
    if (filename.empty())
    {
        throw std::invalid_argument("output filename must not be empty");
    }

    const std::filesystem::path destination(filename);
    const std::filesystem::path directory = destination.has_parent_path()
        ? destination.parent_path()
        : std::filesystem::path(".");
    const std::string destinationName = destination.filename().string();
    if (destinationName.empty() || destinationName == "." ||
        destinationName == "..")
    {
        throw std::invalid_argument(
            "output filename must identify a regular directory entry"
        );
    }

    const int directoryDescriptor = ::open(
        directory.c_str(),
        O_RDONLY | O_DIRECTORY | O_CLOEXEC
    );
    if (directoryDescriptor < 0)
    {
        const int error = errno;
        throw std::runtime_error(systemErrorMessage(
            "failed to open output directory",
            error
        ));
    }

    if (::fsync(directoryDescriptor) != 0)
    {
        const int error = errno;
        static_cast<void>(::close(directoryDescriptor));
        throw std::runtime_error(systemErrorMessage(
            "failed to pre-sync output directory",
            error
        ));
    }

    static std::atomic<std::uint64_t> sequence{0U};
    std::string temporaryName;
    int descriptor = -1;
    for (int attempt = 0; attempt < 128; ++attempt)
    {
        const std::uint64_t id = sequence.fetch_add(
            1U,
            std::memory_order_relaxed
        );
        temporaryName =
            ".make_it_stand_r9." +
            std::to_string(static_cast<long long>(::getpid())) + "." +
            std::to_string(id) + ".tmp";
        descriptor = ::openat(
            directoryDescriptor,
            temporaryName.c_str(),
            O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC,
            0666
        );
        if (descriptor >= 0)
        {
            break;
        }
        if (errno != EEXIST)
        {
            const int error = errno;
            static_cast<void>(::close(directoryDescriptor));
            throw std::runtime_error(systemErrorMessage(
                "failed to create output temporary file",
                error
            ));
        }
    }
    if (descriptor < 0)
    {
        static_cast<void>(::close(directoryDescriptor));
        throw std::runtime_error(
            "output temporary filename retries were exhausted"
        );
    }

    bool committed = false;
    try
    {
        std::size_t offset = 0U;
        while (offset < bytes.size())
        {
            const std::size_t chunk = std::min(
                bytes.size() - offset,
                static_cast<std::size_t>(
                    std::numeric_limits<ssize_t>::max()
                )
            );
            const ssize_t written = ::write(
                descriptor,
                bytes.data() + offset,
                chunk
            );
            if (written < 0)
            {
                if (errno == EINTR)
                {
                    continue;
                }
                const int error = errno;
                throw std::runtime_error(systemErrorMessage(
                    "failed to write output temporary file",
                    error
                ));
            }
            if (written == 0)
            {
                throw std::runtime_error(
                    "output temporary write made no progress"
                );
            }
            offset += static_cast<std::size_t>(written);
        }

        if (::fsync(descriptor) != 0)
        {
            const int error = errno;
            throw std::runtime_error(systemErrorMessage(
                "failed to fsync output temporary file",
                error
            ));
        }
        if (::close(descriptor) != 0)
        {
            const int error = errno;
            descriptor = -1;
            throw std::runtime_error(systemErrorMessage(
                "failed to close output temporary file",
                error
            ));
        }
        descriptor = -1;

        if (::renameat(
                directoryDescriptor,
                temporaryName.c_str(),
                directoryDescriptor,
                destinationName.c_str()
            ) != 0)
        {
            const int error = errno;
            throw std::runtime_error(systemErrorMessage(
                "failed to atomically replace output file",
                error
            ));
        }
        committed = true;

        // The rename is already committed; directory fsync is best effort so
        // a durability limitation cannot be misreported as target preservation.
        static_cast<void>(::fsync(directoryDescriptor));
        static_cast<void>(::close(directoryDescriptor));
    }
    catch (...)
    {
        if (descriptor >= 0)
        {
            static_cast<void>(::close(descriptor));
        }
        if (!committed)
        {
            static_cast<void>(::unlinkat(
                directoryDescriptor,
                temporaryName.c_str(),
                0
            ));
        }
        static_cast<void>(::close(directoryDescriptor));
        throw;
    }
}


std::string serializeDiagnosticObj(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (
        vertices.cols() != 3 ||
        faces.cols() != 3 ||
        !vertices.allFinite() ||
        (faces.size() > 0 &&
         (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows()))
    )
    {
        throw std::logic_error("diagnostic mesh snapshot is invalid");
    }
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(17);
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        output << "v " << vertices(row, 0) << ' '
               << vertices(row, 1) << ' ' << vertices(row, 2) << '\n';
    }
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        output << "f " << faces(row, 0) + 1 << ' '
               << faces(row, 1) + 1 << ' ' << faces(row, 2) + 1 << '\n';
    }
    if (!output)
    {
        throw std::runtime_error("failed to serialize diagnostic OBJ");
    }
    return output.str();
}


manufacturing::ManufacturingMeshFormat formatFromFilename(
    const std::string& filename
)
{
    std::string extension = std::filesystem::path(filename).extension().string();
    std::transform(
        extension.begin(),
        extension.end(),
        extension.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        }
    );
    if (extension == ".obj")
    {
        return manufacturing::ManufacturingMeshFormat::Obj;
    }
    if (extension == ".stl")
    {
        return manufacturing::ManufacturingMeshFormat::BinaryStl;
    }
    throw std::invalid_argument(
        "printable mesh filename must end in .obj or .stl"
    );
}


} // namespace


BalanceRunConfig releasedBetaRunConfig()
{
    BalanceRunConfig config;
    config.initialMu = 0.75;
    config.optimizer.shapeScaleLambda = 20.0;
    config.optimizer.lineSearch.initialStep = 1.0;
    config.optimizer.lineSearch.reductionRatio = 0.8;
    config.optimizer.lineSearch.minimumStep = 0.4;
    config.optimizer.lineSearch.requiredRelativeDecrease = 0.03;
    config.optimizer.lineSearch.muDecrement = 0.05;
    config.optimizer.lineSearch.releasedBetaSemantics = true;
    config.optimizer.freezeTargetsAfterInitialRefresh = true;
    config.optimizer.releasedBetaOuterLoop = true;
    return config;
}


#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER

const char* spinItFasterEvidenceToken(SpinItFasterEvidence evidence) noexcept
{
    switch (evidence)
    {
    case SpinItFasterEvidence::PaperFormulaReproduction:
        return "PAPER_FORMULA_REPRODUCTION";
    case SpinItFasterEvidence::ReconstructedAsset:
        return "RECONSTRUCTED_ASSET";
    }
    return "UNKNOWN";
}


struct SpinItFasterArtifact::Storage
{
    spin::QuadricHollowingResult result;
    spin::SpinDesignRequest request;
    SpinItFasterEvidence evidence =
        SpinItFasterEvidence::PaperFormulaReproduction;
    bool geometrySealed = false;
    std::string diagnostics;
};


namespace
{


void appendSpinVector(std::ostream& output, const Eigen::VectorXd& values)
{
    output << '[';
    for (Eigen::Index index = 0; index < values.size(); ++index)
    {
        if (index != 0) output << ',';
        appendNumber(output, values(index));
    }
    output << ']';
}


std::string spinDiagnosticsJson(
    const spin::QuadricHollowingResult& result,
    const spin::SpinDesignRequest& request,
    SpinItFasterEvidence evidence,
    bool geometrySealed,
    const manufacturing::BoundaryAudit& materialAudit,
    const manufacturing::BoundaryAudit& cavityAudit
)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "{\n"
           << "  \"schema\": \"spin-it-faster-artifact-v1\",\n"
           << "  \"schema_version\": 1,\n"
           << "  \"evidence\": \""
           << spinItFasterEvidenceToken(evidence) << "\",\n"
           << "  \"status\": \""
           << spin::spinSolverStatusToken(result.status) << "\",\n"
           << "  \"mode\": \"" << spin::spinModeToken(request.mode)
           << "\",\n"
           << "  \"static_objective\": \""
           << spin::staticObjectiveToken(request.staticObjective) << "\",\n"
           << "  \"geometry_sealed\": "
           << jsonBoolean(geometrySealed) << ",\n"
           << "  \"input_hash\": \"" << result.inputHash << "\",\n"
           << "  \"geometry_hash\": \"" << result.geometryHash << "\",\n"
           << "  \"trace_hash\": \"" << result.traceHash << "\",\n"
           << "  \"output_hash\": \"" << result.outputHash << "\",\n"
           << "  \"normalized_equation_residual\": ";
    appendNumber(output, result.normalizedEquationResidual);
    output << ",\n  \"physical_relative_residual\": ";
    appendNumber(output, result.physicalRelativeResidual);
    if (request.mode == spin::SpinMode::BuoyantStability)
    {
        output << ",\n  \"buoyant_metacenter_height\": ";
        appendNumber(output, result.buoyantMetacenterHeight);
        output << ",\n  \"buoyant_stability_margin\": ";
        appendNumber(output, result.buoyantStabilityMargin);
    }
    else
    {
        output << ",\n  \"buoyant_metacenter_height\": null"
               << ",\n  \"buoyant_stability_margin\": null";
    }
    output << ",\n  \"quadric_coefficients\": ";
    appendSpinVector(output, result.quadricCoefficients);
    output << ",\n  \"final_raw_moments\": ";
    appendSpinVector(output, result.finalMoments.values());
    const bool validMoments = result.finalMoments.isFinite() &&
        result.finalMoments.mass() > 0.0;
    if (validMoments)
    {
        const Eigen::Vector3d centerTarget =
            result.finalMoments.centerOfMass();
        const Eigen::Vector3d centerModel =
            request.targetFrame.axes * centerTarget +
            request.targetFrame.origin;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(
            result.finalMoments.inertiaAboutCenter()
        );
        const Eigen::Matrix3d principalAxesModel =
            request.targetFrame.axes * eigen.eigenvectors();
        output << ",\n  \"center_of_mass_model\": ";
        appendVector3(output, centerModel);
        output << ",\n  \"target_origin_model\": ";
        appendVector3(output, request.targetFrame.origin);
        output << ",\n  \"target_axes_model\": ";
        appendMatrix3(output, request.targetFrame.axes);
        output << ",\n  \"principal_axes_model\": ";
        appendMatrix3(output, principalAxesModel);
    }
    else
    {
        output << ",\n  \"center_of_mass_model\": null"
               << ",\n  \"target_origin_model\": null"
               << ",\n  \"target_axes_model\": null"
               << ",\n  \"principal_axes_model\": null";
    }
    output << ",\n  \"geometry\": {"
           << "\"outer_vertices\":" << result.outerSurface.vertices.rows()
           << ",\"outer_faces\":" << result.outerSurface.faces.rows()
           << ",\"design_vertices\":" << request.designDomain.vertices.rows()
           << ",\"design_faces\":" << request.designDomain.faces.rows()
           << ",\"material_vertices\":"
           << result.materialBoundary.vertices.rows()
           << ",\"material_faces\":"
           << result.materialBoundary.faces.rows()
           << ",\"material_boundary_valid\":"
           << jsonBoolean(materialAudit.strictBoundaryValid())
           << ",\"cavity_vertices\":"
           << result.cavityBoundary.vertices.rows()
           << ",\"cavity_faces\":" << result.cavityBoundary.faces.rows()
           << ",\"cavity_boundary_valid\":"
           << jsonBoolean(cavityAudit.strictBoundaryValid()) << "},\n"
           << "  \"trace\": [";
    for (std::size_t index = 0; index < result.trace.size(); ++index)
    {
        if (index != 0U) output << ',';
        const spin::SpinIterationTrace& entry = result.trace[index];
        output << "{\"restart\":" << entry.restart
               << ",\"iteration\":" << entry.iteration
               << ",\"function_evaluations\":"
               << entry.functionEvaluations
               << ",\"normalized_equation_residual\":";
        appendNumber(output, entry.normalizedEquationResidual);
        output << ",\"physical_relative_residual\":";
        appendNumber(output, entry.physicalRelativeResidual);
        output << ",\"objective\":";
        appendNumber(output, entry.objective);
        output << ",\"multiplier_mu\":";
        appendNumber(output, entry.multiplierMu);
        output << ",\"accepted_geodesic_step\":";
        appendNumber(output, entry.acceptedGeodesicStep);
        output << ",\"quadric_hash\":\"" << entry.quadricHash
               << "\",\"material_geometry_hash\":\""
               << entry.materialGeometryHash
               << "\",\"material_topology_hash\":\""
               << entry.materialTopologyHash << "\"}";
    }
    output << "]\n}\n";
    return output.str();
}


} // namespace


SpinItFasterArtifact::SpinItFasterArtifact(
    std::shared_ptr<const Storage> storage
)
    : storage_(std::move(storage))
{
    if (!storage_)
    {
        throw std::invalid_argument("Spin-It Faster artifact storage must exist");
    }
}


SpinItFasterArtifact::~SpinItFasterArtifact() = default;


std::string SpinItFasterArtifact::statusToken() const
{
    return spin::spinSolverStatusToken(storage_->result.status);
}


std::string SpinItFasterArtifact::evidenceToken() const
{
    return spinItFasterEvidenceToken(storage_->evidence);
}


bool SpinItFasterArtifact::succeeded() const
{
    return storage_->result.status == spin::SpinSolverStatus::Converged;
}


bool SpinItFasterArtifact::geometrySealed() const
{
    return storage_->geometrySealed;
}


const spin::QuadricHollowingResult& SpinItFasterArtifact::result() const
{
    return storage_->result;
}


const Eigen::MatrixXd& SpinItFasterArtifact::outerVertices() const
{
    return storage_->result.outerSurface.vertices;
}


const Eigen::MatrixXi& SpinItFasterArtifact::outerFaces() const
{
    return storage_->result.outerSurface.faces;
}


const Eigen::MatrixXd& SpinItFasterArtifact::designVertices() const
{
    return storage_->request.designDomain.vertices;
}


const Eigen::MatrixXi& SpinItFasterArtifact::designFaces() const
{
    return storage_->request.designDomain.faces;
}


const Eigen::MatrixXd& SpinItFasterArtifact::materialVertices() const
{
    return storage_->result.materialBoundary.vertices;
}


const Eigen::MatrixXi& SpinItFasterArtifact::materialFaces() const
{
    return storage_->result.materialBoundary.faces;
}


const Eigen::MatrixXd& SpinItFasterArtifact::cavityVertices() const
{
    return storage_->result.cavityBoundary.vertices;
}


const Eigen::MatrixXi& SpinItFasterArtifact::cavityFaces() const
{
    return storage_->result.cavityBoundary.faces;
}


const std::string& SpinItFasterArtifact::diagnosticsJson() const
{
    return storage_->diagnostics;
}


void SpinItFasterArtifact::saveOuterMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(outerVertices(), outerFaces())
    );
}


void SpinItFasterArtifact::saveDesignMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(designVertices(), designFaces())
    );
}


void SpinItFasterArtifact::saveMaterialMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(materialVertices(), materialFaces())
    );
}


void SpinItFasterArtifact::saveCavityMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(cavityVertices(), cavityFaces())
    );
}


void SpinItFasterArtifact::saveDiagnostics(const std::string& filename) const
{
    writeAtomicBytes(filename, diagnosticsJson());
}


SpinItFasterArtifact solveSpinItFaster(
    const spin::SpinDesignRequest& request,
    SpinItFasterEvidence evidence,
    const spin::QuadricMomentEvaluatorSettings& evaluatorSettings,
    const std::vector<spin::QuadricCoefficients10>& initialGuesses
)
{
    spin::QuadricHollowingResult result = spin::solveQuadricHollowing(
        request,
        evaluatorSettings,
        initialGuesses
    );
    const manufacturing::BoundaryAudit materialAudit =
        manufacturing::auditTriangleBoundary(
            result.materialBoundary.vertices,
            result.materialBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive,
            result.status != spin::SpinSolverStatus::Converged
        );
    const manufacturing::BoundaryAudit cavityAudit =
        manufacturing::auditTriangleBoundary(
            result.cavityBoundary.vertices,
            result.cavityBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive,
            result.status != spin::SpinSolverStatus::Converged
        );
    const bool geometrySealed =
        result.status == spin::SpinSolverStatus::Converged &&
        materialAudit.strictBoundaryValid() &&
        cavityAudit.strictBoundaryValid() &&
        result.inputHash.size() == 64U && result.geometryHash.size() == 64U &&
        result.traceHash.size() == 64U && result.outputHash.size() == 64U;
    auto storage = std::make_shared<SpinItFasterArtifact::Storage>();
    storage->result = std::move(result);
    storage->request = request;
    storage->evidence = evidence;
    storage->geometrySealed = geometrySealed;
    storage->diagnostics = spinDiagnosticsJson(
        storage->result,
        storage->request,
        evidence,
        geometrySealed,
        materialAudit,
        cavityAudit
    );
    return SpinItFasterArtifact(std::move(storage));
}

#endif


struct BalanceResult::Storage
{
    optimization::AlternatingBalanceOptimizationResult result;
    manufacturing::ManufacturingGeometry manufacturing;
    std::string diagnostics;
    std::string performanceReport;

    Storage(
        optimization::AlternatingBalanceOptimizationResult resultValue,
        manufacturing::ManufacturingGeometry manufacturingValue,
        std::string diagnosticsValue,
        std::string performanceReportValue
    )
        :
        result(std::move(resultValue)),
        manufacturing(std::move(manufacturingValue)),
        diagnostics(std::move(diagnosticsValue)),
        performanceReport(std::move(performanceReportValue))
    {
    }
};


struct BalanceModel::Impl
{
    struct ReleasedBetaCache
    {
        int voxelResolution = 0;
        int hullDepth = 2;
        std::vector<int> cellIdsZFast;
        Eigen::MatrixXd coordinateMappedNodeWeights;
    };

    struct MaterialOccupancyCache
    {
        std::shared_ptr<const voxel::VoxelDomain> domain;
        int protectedBoundaryLayers = 1;
    };

    struct Standing
    {
        stability::SupportRegion support;
        double gammaRadians = 0.0;
        Eigen::MatrixXd fixedPoints;
    };

    using Standings = std::vector<
        Standing,
        Eigen::aligned_allocator<Standing>
    >;

    struct Suspension
    {
        Eigen::Vector3d attachment = Eigen::Vector3d::Zero();
        double angleRadians = 0.0;
    };

    struct UserHandle
    {
        Eigen::Vector3d pivot = Eigen::Vector3d::Zero();
        Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    };

    struct Snapshot
    {
        Eigen::MatrixXd vertices;
        Eigen::MatrixXi faces;
        BalanceModelConfig config;
        Eigen::Vector3d gravity;
        Standings standings;
        std::vector<Suspension, Eigen::aligned_allocator<Suspension>>
            suspensions;
        std::vector<UserHandle, Eigen::aligned_allocator<UserHandle>>
            userHandles;
        std::shared_ptr<const ReleasedBetaCache> releasedBetaCache;
        std::shared_ptr<const MaterialOccupancyCache> materialOccupancyCache;
    };

    mutable std::mutex mutex;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    BalanceModelConfig config;
    Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -1.0, 0.0);
    Standings standings;
    std::vector<Suspension, Eigen::aligned_allocator<Suspension>> suspensions;
    std::vector<UserHandle, Eigen::aligned_allocator<UserHandle>> userHandles;
    std::shared_ptr<const ReleasedBetaCache> releasedBetaCache;
    std::shared_ptr<const MaterialOccupancyCache> materialOccupancyCache;

    Snapshot captureSnapshot() const
    {
        std::lock_guard<std::mutex> lock(mutex);
        Snapshot snapshot;
        snapshot.vertices = vertices;
        snapshot.faces = faces;
        snapshot.config = config;
        snapshot.gravity = gravity;
        snapshot.standings = standings;
        snapshot.suspensions = suspensions;
        snapshot.userHandles = userHandles;
        snapshot.releasedBetaCache = releasedBetaCache;
        snapshot.materialOccupancyCache = materialOccupancyCache;
        return snapshot;
    }

    static VoxelGrid makeVoxelGrid(const Snapshot& snapshot)
    {
        if (snapshot.releasedBetaCache)
        {
            return VoxelGrid::fromReleasedBetaOccupancy(
                snapshot.releasedBetaCache->voxelResolution,
                snapshot.releasedBetaCache->cellIdsZFast,
                snapshot.releasedBetaCache->hullDepth
            );
        }
        if (snapshot.materialOccupancyCache)
        {
            return VoxelGrid::fromMaterialOccupancy(
                snapshot.materialOccupancyCache->domain,
                snapshot.config.minimumWallThickness,
                snapshot.materialOccupancyCache->protectedBoundaryLayers
            );
        }
        return VoxelGrid(
            snapshot.vertices,
            snapshot.faces,
            snapshot.config.minimumWallThickness,
            snapshot.config.voxelResolution,
            snapshot.config.fixedUnitCubeVoxelGrid
        );
    }
};


BalanceResult::BalanceResult(std::shared_ptr<const Storage> storage)
    : storage_(std::move(storage))
{
    if (!storage_)
    {
        throw std::invalid_argument("balance result storage must exist");
    }
}


BalanceResult::~BalanceResult() = default;


const optimization::AlternatingBalanceOptimizationResult&
BalanceResult::result() const
{
    return storage_->result;
}


const manufacturing::ManufacturingGeometry&
BalanceResult::manufacturingGeometry() const
{
    return storage_->manufacturing;
}


std::string BalanceResult::statusToken() const
{
    return api::statusToken(storage_->result.status);
}


bool BalanceResult::succeeded() const
{
    return optimizationSucceeded(storage_->result);
}


const Eigen::MatrixXd& BalanceResult::outerVertices() const
{
    return storage_->result.geometry.outerVertices;
}


const Eigen::MatrixXi& BalanceResult::outerFaces() const
{
    return storage_->result.geometry.outerFaces;
}


const Eigen::MatrixXd& BalanceResult::innerVertices() const
{
    return storage_->result.geometry.innerVertices;
}


const Eigen::MatrixXi& BalanceResult::innerFaces() const
{
    return storage_->result.geometry.innerFaces;
}


bool BalanceResult::manufacturingReady() const
{
    return storage_->manufacturing.report.ready;
}


const std::vector<std::string>&
BalanceResult::manufacturingFailures() const
{
    return storage_->manufacturing.report.failures;
}


std::string BalanceResult::manufacturingReportJson() const
{
    return storage_->manufacturing.report.toJson();
}


const std::string& BalanceResult::diagnosticsJson() const
{
    return storage_->diagnostics;
}


const std::string& BalanceResult::performanceReportJson() const
{
    return storage_->performanceReport;
}


void BalanceResult::saveOuterMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(outerVertices(), outerFaces())
    );
}


void BalanceResult::saveInnerMesh(const std::string& filename) const
{
    writeAtomicBytes(
        filename,
        serializeDiagnosticObj(innerVertices(), innerFaces())
    );
}


void BalanceResult::savePrintableMesh(const std::string& filename) const
{
    savePrintableMesh(filename, formatFromFilename(filename));
}


void BalanceResult::savePrintableMesh(
    const std::string& filename,
    manufacturing::ManufacturingMeshFormat format
) const
{
    // The sealed R8 exporter already performs its own fd-based transactional
    // write.  Do not wrap it in a second, weaker filename transaction.
    manufacturing::exportReadyManufacturingMesh(
        storage_->manufacturing,
        filename,
        format
    );
}


void BalanceResult::saveDiagnostics(const std::string& filename) const
{
    writeAtomicBytes(filename, storage_->diagnostics);
}


BalanceModel::BalanceModel()
    : impl_(std::make_unique<Impl>())
{
}


BalanceModel::BalanceModel(
    const std::string& meshPath,
    const BalanceModelConfig& config
)
    : impl_(std::make_unique<Impl>())
{
    validateModelConfig(config);
    if (meshPath.empty())
    {
        throw std::invalid_argument("source mesh path must not be empty");
    }
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    if (!igl::read_triangle_mesh(meshPath, vertices, faces))
    {
        throw std::runtime_error("failed to read source triangle mesh");
    }
    validateSourceMesh(vertices, faces);
    impl_->vertices = std::move(vertices);
    impl_->faces = std::move(faces);
    impl_->config = config;
}


BalanceModel::BalanceModel(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const BalanceModelConfig& config
)
    : impl_(std::make_unique<Impl>())
{
    validateModelConfig(config);
    validateSourceMesh(vertices, faces);
    impl_->vertices = vertices;
    impl_->faces = faces;
    impl_->config = config;
}


std::unique_ptr<BalanceModel> BalanceModel::fromMaterialOccupancy(
    const MaterialOccupancyData& occupancy,
    const BalanceModelConfig& config
)
{
    validateModelConfig(config);
    if (occupancy.protectedBoundaryLayers < 0)
    {
        throw std::invalid_argument(
            "material occupancy protectedBoundaryLayers must be nonnegative"
        );
    }
    if (config.voxelResolution != occupancy.dimensions.maxCoeff())
    {
        throw std::invalid_argument(
            "material occupancy dimensions must match config.voxelResolution"
        );
    }

    voxel::GridSpec specification;
    specification.dimensions = occupancy.dimensions;
    specification.minimum = occupancy.minimum;
    specification.maximum = occupancy.maximum;
    const std::shared_ptr<const voxel::VoxelDomain> domain =
        voxel::VoxelDomain::create(
            specification,
            occupancy.materialMaskXFast
        );
    const voxel::MaterialBoundaryMesh boundaries =
        voxel::extractMaterialBoundaries(domain);
    if (boundaries.materialComponentCount != 1)
    {
        throw std::domain_error(
            "material occupancy must contain one six-connected component"
        );
    }
    const manufacturing::BoundaryAudit outer =
        manufacturing::auditTriangleBoundary(
            boundaries.outerVertices,
            boundaries.outerFaces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    const manufacturing::BoundaryAudit inner =
        manufacturing::auditTriangleBoundary(
            boundaries.innerVertices,
            boundaries.innerFaces,
            manufacturing::BoundaryOrientation::InwardNegative,
            true
        );
    if (
        !outer.strictBoundaryValid() ||
        outer.faceComponentCount != 1 ||
        (!inner.empty && !inner.strictBoundaryValid())
    )
    {
        throw std::domain_error(
            "material occupancy boundary reconstruction failed strict audit"
        );
    }

    auto cache = std::make_shared<Impl::MaterialOccupancyCache>();
    cache->domain = domain;
    cache->protectedBoundaryLayers = occupancy.protectedBoundaryLayers;

    std::unique_ptr<BalanceModel> model(new BalanceModel());
    model->impl_->vertices = boundaries.outerVertices;
    model->impl_->faces = boundaries.outerFaces;
    model->impl_->config = config;
    model->impl_->materialOccupancyCache = std::move(cache);
    return model;
}


std::unique_ptr<BalanceModel> BalanceModel::fromReleasedBetaArchive(
    const Eigen::MatrixXd& normalizedVertices,
    const Eigen::MatrixXi& faces,
    const ReleasedBetaArchiveData& archive,
    const BalanceModelConfig& config
)
{
    validateModelConfig(config);
    if (
        !config.fixedUnitCubeVoxelGrid ||
        config.minimumWallThickness != 0.0 ||
        config.voxelResolution != archive.voxelResolution
    )
    {
        throw std::invalid_argument(
            "released-beta archive import requires a matching fixed unit-cube "
            "configuration with zero modern wall thickness"
        );
    }
    validateReleasedBetaSourceMesh(normalizedVertices, faces);
    Eigen::MatrixXd mappedWeights =
        normalizeAndCoordinateMapReleasedBetaWeights(archive);

    auto cache = std::make_shared<Impl::ReleasedBetaCache>();
    cache->voxelResolution = archive.voxelResolution;
    cache->hullDepth = archive.hullDepth;
    cache->cellIdsZFast = archive.cellIdsZFast;
    cache->coordinateMappedNodeWeights = std::move(mappedWeights);

    std::unique_ptr<BalanceModel> model(new BalanceModel());
    model->impl_->vertices = normalizedVertices;
    model->impl_->faces = faces;
    model->impl_->config = config;
    model->impl_->releasedBetaCache = std::move(cache);
    return model;
}


BalanceModel::~BalanceModel() = default;


Eigen::MatrixXd BalanceModel::vertices() const
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    return impl_->vertices;
}


Eigen::MatrixXi BalanceModel::faces() const
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    return impl_->faces;
}


Eigen::Vector3d BalanceModel::gravity() const
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    return impl_->gravity;
}


BalanceModelConfig BalanceModel::config() const
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    return impl_->config;
}


void BalanceModel::setGravity(const Eigen::Vector3d& gravityDirection)
{
    const Eigen::Vector3d normalized =
        normalizedGravityDirection(gravityDirection);
    std::lock_guard<std::mutex> lock(impl_->mutex);
    if (!impl_->standings.empty() || !impl_->suspensions.empty())
    {
        throw std::logic_error(
            "gravity cannot change after balance targets have been added"
        );
    }
    impl_->gravity = normalized;
}


void BalanceModel::setSupportPolygon(
    const Eigen::MatrixXd& polygon,
    double gammaRadians
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromPolygon(
        impl_->vertices,
        polygon,
        impl_->gravity,
        gammaRadians
    );
    impl_->standings.clear();
    impl_->standings.push_back(Impl::Standing{
        std::move(standing.support),
        standing.gammaRadians,
        std::move(standing.fixedPoints)
    });
}


void BalanceModel::setSupportFromMesh(
    double contactThreshold,
    double gammaRadians,
    double minimumStabilityInset,
    bool supportPlaneAtBandTop
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromMesh(
        impl_->vertices,
        impl_->gravity,
        contactThreshold,
        gammaRadians,
        minimumStabilityInset,
        supportPlaneAtBandTop
    );
    impl_->standings.clear();
    impl_->standings.push_back(Impl::Standing{
        std::move(standing.support),
        standing.gammaRadians,
        std::move(standing.fixedPoints)
    });
}


void BalanceModel::setSupportFromReferenceMesh(
    const Eigen::MatrixXd& referenceVertices,
    double contactThreshold,
    double gammaRadians,
    double minimumStabilityInset,
    bool supportPlaneAtBandTop
)
{
    if (
        referenceVertices.rows() == 0 ||
        referenceVertices.cols() != 3 ||
        !referenceVertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "reference support vertices must be a nonempty finite N-by-3 matrix"
        );
    }
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromMesh(
        referenceVertices,
        impl_->gravity,
        contactThreshold,
        gammaRadians,
        minimumStabilityInset,
        supportPlaneAtBandTop
    );
    impl_->standings.clear();
    impl_->standings.push_back(Impl::Standing{
        std::move(standing.support),
        standing.gammaRadians,
        std::move(standing.fixedPoints)
    });
}


void BalanceModel::addStandingConstraintPolygon(
    const Eigen::MatrixXd& polygon,
    const Eigen::Vector3d& gravityDirection,
    double gammaRadians
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromPolygon(
        impl_->vertices,
        polygon,
        gravityDirection,
        gammaRadians
    );
    impl_->standings.push_back(Impl::Standing{
        std::move(standing.support),
        standing.gammaRadians,
        std::move(standing.fixedPoints)
    });
}


void BalanceModel::addStandingConstraintFromMesh(
    const Eigen::Vector3d& gravityDirection,
    double contactThreshold,
    double gammaRadians,
    double minimumStabilityInset,
    bool supportPlaneAtBandTop
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromMesh(
        impl_->vertices,
        gravityDirection,
        contactThreshold,
        gammaRadians,
        minimumStabilityInset,
        supportPlaneAtBandTop
    );
    impl_->standings.push_back(Impl::Standing{
        std::move(standing.support),
        standing.gammaRadians,
        std::move(standing.fixedPoints)
    });
}


void BalanceModel::addHandle(
    const Eigen::Vector3d& restPivot,
    const Eigen::Matrix3d& fixedRotation
)
{
    if (!restPivot.allFinite())
    {
        throw std::invalid_argument("handle pivot must be finite");
    }
    if (!isProperRotation(fixedRotation))
    {
        throw std::invalid_argument(
            "handle rotation must be finite, orthonormal, and proper"
        );
    }
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->userHandles.push_back(Impl::UserHandle{
        restPivot,
        fixedRotation
    });
}


void BalanceModel::addSuspensionConstraint(
    const Eigen::Vector3d& attachmentPoint,
    double allowedDeviationRadians
)
{
    if (!attachmentPoint.allFinite())
    {
        throw std::invalid_argument(
            "suspension attachment must be finite"
        );
    }
    if (
        !std::isfinite(allowedDeviationRadians) ||
        allowedDeviationRadians < 0.0 ||
        allowedDeviationRadians >= HALF_PI
    )
    {
        throw std::invalid_argument(
            "suspension angle must be finite and in [0, pi/2)"
        );
    }
    std::lock_guard<std::mutex> lock(impl_->mutex);
    static_cast<void>(stability::SuspensionTargetSpecification(
        attachmentPoint,
        impl_->gravity,
        allowedDeviationRadians
    ));
    impl_->suspensions.push_back(Impl::Suspension{
        attachmentPoint,
        allowedDeviationRadians
    });
}


void BalanceModel::clearHandles()
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->userHandles.clear();
}


void BalanceModel::clearConstraints()
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->standings.clear();
    impl_->suspensions.clear();
}


#if MAKE_IT_STAND_HAS_REACHABILITY
reachability::CarvingReachabilityResult
BalanceModel::analyzeCarvingReachability(std::size_t standingIndex) const
{
    const Impl::Snapshot snapshot = impl_->captureSnapshot();
    if (standingIndex >= snapshot.standings.size())
    {
        throw std::out_of_range(
            "standing constraint index is outside the model"
        );
    }

    const VoxelGrid grid = Impl::makeVoxelGrid(snapshot);
    const deformation::GridDeformationDomain domain =
        deformation::GridDeformationDomain::fromVoxelGrid(grid);
    const deformation::DeformedVoxelMassField massField(
        domain,
        domain.fullRestCorners(),
        1.0,
        0.0
    );
    const Impl::Standing& standing = snapshot.standings[standingIndex];

    MassMoments totalMoments;
    double fixedMass = 0.0;
    Eigen::Vector2d fixedFirstMoment2D = Eigen::Vector2d::Zero();
    std::vector<reachability::ReachabilityCell> variableCells;
    variableCells.reserve(static_cast<std::size_t>(domain.activeCellCount()));

    for (Eigen::Index cell = 0; cell < domain.cellCount(); ++cell)
    {
        if (!domain.isActiveCell(static_cast<int>(cell))) continue;
        const Eigen::Vector3i coordinate =
            domain.cellCoordinates(static_cast<int>(cell));
        const Voxel& voxelState = grid(
            coordinate.x(),
            coordinate.y(),
            coordinate.z()
        );
        if (!voxelState.isFilled) continue;

        const MassMoments moments = massField.massMoments(
            static_cast<int>(cell)
        );
        totalMoments += moments;
        const Eigen::Vector3d center =
            moments.firstMoment / moments.mass;
        const Eigen::Vector2d projectedFirstMoment = moments.mass *
            standing.support.toPlaneCoordinates(center);
        if (voxelState.isInner)
        {
            variableCells.push_back(reachability::ReachabilityCell{
                static_cast<int>(cell),
                moments.mass,
                projectedFirstMoment
            });
        }
        else
        {
            fixedMass += moments.mass;
            fixedFirstMoment2D += projectedFirstMoment;
        }
    }
    if (!(totalMoments.mass > 0.0) || !(fixedMass > 0.0))
    {
        throw std::domain_error(
            "reachability requires positive total and fixed shell mass"
        );
    }

    const stability::StandingTarget target = stability::makeStandingTarget(
        standing.support,
        totalMoments.firstMoment / totalMoments.mass,
        standing.gammaRadians
    );

    reachability::CarvingReachabilityInput input;
    input.fixedMass = fixedMass;
    input.fixedFirstMoment2D = fixedFirstMoment2D;
    input.variableCells = variableCells;
    input.safePolygon2D = target.safeRegion2D();
    if (!target.requestedMarginAttainable() || target.safeRegion2D().rows() < 3)
    {
        reachability::CarvingReachabilityResult result;
        result.inputHash = reachability::computeInputHash(input);
        result.status = reachability::ReachabilityStatus::Indeterminate;
        result.inputCellCount = static_cast<int>(variableCells.size());
        result.variableCellCount = result.inputCellCount;
        result.fixedMass = fixedMass;
        result.primalSolverStatus = "DEGENERATE_OR_UNATTAINABLE_SAFE_REGION";
        return result;
    }

    input.variableCells = std::move(variableCells);
    return reachability::CarvingReachabilityAnalyzer().analyze(input);
}
#endif


BalanceResult BalanceModel::optimize(const BalanceRunConfig& runConfig) const
{
    validateRunConfig(runConfig);
    const PerformanceClock::time_point totalStart = PerformanceClock::now();
    PerformanceDurations durations;
    PerformanceClock::time_point phaseStart = totalStart;
    const Impl::Snapshot snapshot = impl_->captureSnapshot();
    durations.inputSnapshot = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );
    if (snapshot.standings.empty() && snapshot.suspensions.empty())
    {
        throw std::logic_error(
            "at least one standing or suspension target is required"
        );
    }

    phaseStart = PerformanceClock::now();
    VoxelGrid grid = Impl::makeVoxelGrid(snapshot);
    durations.voxelization = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );

    phaseStart = PerformanceClock::now();
    const deformation::GridDeformationDomain domain =
        deformation::GridDeformationDomain::fromVoxelGrid(grid);

    Eigen::Index fixedPointCount =
        static_cast<Eigen::Index>(snapshot.suspensions.size());
    for (const auto& standing : snapshot.standings)
    {
        fixedPointCount += standing.fixedPoints.rows();
    }
    Eigen::MatrixXd fixedPoints(fixedPointCount, 3);
    Eigen::Index destination = 0;
    for (const auto& standing : snapshot.standings)
    {
        fixedPoints.middleRows(
            destination,
            standing.fixedPoints.rows()
        ) = standing.fixedPoints;
        destination += standing.fixedPoints.rows();
    }
    for (const auto& suspension : snapshot.suspensions)
    {
        fixedPoints.row(destination++) = suspension.attachment.transpose();
    }

    deformation::Handle fixed;
    if (snapshot.releasedBetaCache)
    {
        fixed.restPivot = fixedPoints.colwise().mean().transpose();
        fixed.role = deformation::HandleRole::Support;
        fixed.locked = true;
        fixed.constrainedNodes.reserve(
            static_cast<std::size_t>(fixedPoints.rows())
        );
        for (Eigen::Index row = 0; row < fixedPoints.rows(); ++row)
        {
            fixed.constrainedNodes.push_back(releasedBetaClosestActiveNode(
                domain,
                fixedPoints.row(row).transpose()
            ));
        }
    }
    else
    {
        fixed = deformation::makeFixedSupportHandle(domain, fixedPoints);
    }
    std::set<int> occupiedNodes(
        fixed.constrainedNodes.begin(),
        fixed.constrainedNodes.end()
    );
    deformation::HandleVector handles;
    handles.reserve(snapshot.userHandles.size() + 1U);
    handles.push_back(std::move(fixed));
    for (const auto& specification : snapshot.userHandles)
    {
        const int node = snapshot.releasedBetaCache
            ? releasedBetaClosestActiveNode(domain, specification.pivot)
            : domain.nearestActiveNode(specification.pivot);
        if (!snapshot.releasedBetaCache)
        {
            static_cast<void>(domain.embed(specification.pivot));
        }
        if (!occupiedNodes.insert(node).second)
        {
            throw std::domain_error(
                "user handle conflicts with a fixed or previous handle node"
            );
        }
        deformation::Handle handle;
        handle.restPivot = snapshot.releasedBetaCache
            ? domain.restNodes().row(node).transpose()
            : specification.pivot;
        handle.fixedRotation = specification.rotation;
        handle.role = deformation::HandleRole::User;
        handle.locked = false;
        handle.constrainedNodes = {node};
        handles.push_back(std::move(handle));
    }
    const deformation::HandleSet binding(std::move(handles));
    durations.domainAndHandles = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );

    phaseStart = PerformanceClock::now();
    std::unique_ptr<deformation::DeformationSystem> deformationSystem;
    if (snapshot.releasedBetaCache)
    {
        if (
            snapshot.releasedBetaCache->coordinateMappedNodeWeights.cols() !=
            static_cast<Eigen::Index>(binding.size())
        )
        {
            throw std::invalid_argument(
                "released-beta .bbw column count does not match support/user handles"
            );
        }
        deformationSystem = std::make_unique<deformation::DeformationSystem>(
            grid,
            snapshot.vertices,
            snapshot.faces,
            binding,
            snapshot.releasedBetaCache->coordinateMappedNodeWeights,
            true
        );
    }
    else
    {
        deformationSystem = std::make_unique<deformation::DeformationSystem>(
            grid,
            snapshot.vertices,
            snapshot.faces,
            binding,
            snapshot.config.bbw
        );
    }
    durations.deformationSystemAndBbw = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );

    phaseStart = PerformanceClock::now();
    stability::BalanceTargetSpecifications targets;
    targets.reserve(
        snapshot.suspensions.size() +
        snapshot.standings.size()
    );
    for (const auto& standing : snapshot.standings)
    {
        targets.emplace_back(stability::StandingTargetSpecification(
            standing.support,
            standing.gammaRadians
        ));
    }
    for (const auto& suspension : snapshot.suspensions)
    {
        targets.emplace_back(stability::SuspensionTargetSpecification(
            suspension.attachment,
            snapshot.gravity,
            suspension.angleRadians
        ));
    }
    const stability::BalanceTargetSchedule schedule(
        targets,
        runConfig.optimizer.freezeTargetsAfterInitialRefresh
    );
    const optimization::AlternatingBalanceOptimizer optimizer(
        *deformationSystem,
        grid,
        schedule,
        runConfig.optimizer
    );
    durations.targetAndOptimizerSetup = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );
    phaseStart = PerformanceClock::now();
    optimization::AlternatingBalanceOptimizationResult result =
        optimizer.optimize(
            binding.optimizationParameters(),
            runConfig.initialMu
        );
    durations.optimization = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );

    phaseStart = PerformanceClock::now();
    manufacturing::ManufacturingOptions manufacturingOptions;
    manufacturingOptions.density = runConfig.optimizer.density;
    manufacturing::ManufacturingGeometry manufacturing =
        manufacturing::prepareManufacturingGeometry(
            result,
            manufacturingOptions
        );
    durations.manufacturingAudit = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );
    phaseStart = PerformanceClock::now();
    DiagnosticsStandings diagnosticsStandings;
    diagnosticsStandings.reserve(snapshot.standings.size());
    for (const auto& standing : snapshot.standings)
    {
        diagnosticsStandings.push_back(DiagnosticsStanding{
            standing.fixedPoints,
            standing.support.gravity(),
            standing.gammaRadians
        });
    }
    DiagnosticsSuspensions diagnosticsSuspensions;
    diagnosticsSuspensions.reserve(snapshot.suspensions.size());
    for (const auto& suspension : snapshot.suspensions)
    {
        diagnosticsSuspensions.push_back(DiagnosticsSuspension{
            suspension.attachment,
            suspension.angleRadians
        });
    }
    DiagnosticsHandles diagnosticsHandles;
    diagnosticsHandles.reserve(snapshot.userHandles.size());
    for (const auto& handle : snapshot.userHandles)
    {
        diagnosticsHandles.push_back(DiagnosticsHandle{
            handle.pivot,
            handle.rotation
        });
    }
    durations.diagnosticsPayloadSetup = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );
    phaseStart = PerformanceClock::now();
    const std::string diagnostics = buildDiagnosticsJson(
        snapshot.config,
        runConfig,
        static_cast<bool>(snapshot.releasedBetaCache),
        static_cast<int>(snapshot.vertices.rows()),
        static_cast<int>(snapshot.faces.rows()),
        snapshot.gravity,
        diagnosticsStandings,
        diagnosticsSuspensions,
        diagnosticsHandles,
        result,
        manufacturing
    );
    durations.diagnosticsSerialization = elapsedSeconds(
        phaseStart,
        PerformanceClock::now()
    );
    durations.total = elapsedSeconds(totalStart, PerformanceClock::now());
    const std::string performanceReport =
        buildPerformanceReportJson(durations);
    return BalanceResult(std::make_shared<const BalanceResult::Storage>(
        std::move(result),
        std::move(manufacturing),
        diagnostics,
        performanceReport
    ));
}


} // namespace api
