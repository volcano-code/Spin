#include "api/BalanceModel.h"

#include "deformation/DeformationSystem.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "physics/BalanceEnergy.h"
#include "stability/BalanceTargetSchedule.h"
#include "stability/SupportRegion.h"
#include "voxel/VoxelGrid.h"

#include <igl/read_triangle_mesh.h>

#include <Eigen/Geometry>

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
    double gammaRadians
)
{
    const Eigen::Vector3d gravity =
        normalizedGravityDirection(gravityDirection);
    const stability::SupportRegion support = stability::SupportRegion::fromMesh(
        vertices,
        gravity,
        contactThreshold
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
               << record.alphaHashAfter << "\",\"geometry_hash_before\":\""
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
    };

    mutable std::mutex mutex;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    BalanceModelConfig config;
    Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -1.0, 0.0);
    Standings standings;
    std::vector<Suspension, Eigen::aligned_allocator<Suspension>> suspensions;
    std::vector<UserHandle, Eigen::aligned_allocator<UserHandle>> userHandles;
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
    double gammaRadians
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromMesh(
        impl_->vertices,
        impl_->gravity,
        contactThreshold,
        gammaRadians
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
    double gammaRadians
)
{
    std::lock_guard<std::mutex> lock(impl_->mutex);
    StandingConstruction standing = constructStandingFromMesh(
        impl_->vertices,
        gravityDirection,
        contactThreshold,
        gammaRadians
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


BalanceResult BalanceModel::optimize(const BalanceRunConfig& runConfig) const
{
    validateRunConfig(runConfig);
    const PerformanceClock::time_point totalStart = PerformanceClock::now();
    PerformanceDurations durations;
    PerformanceClock::time_point phaseStart = totalStart;
    Impl::Snapshot snapshot;
    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        snapshot.vertices = impl_->vertices;
        snapshot.faces = impl_->faces;
        snapshot.config = impl_->config;
        snapshot.gravity = impl_->gravity;
        snapshot.standings = impl_->standings;
        snapshot.suspensions = impl_->suspensions;
        snapshot.userHandles = impl_->userHandles;
    }
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
    VoxelGrid grid(
        snapshot.vertices,
        snapshot.faces,
        snapshot.config.minimumWallThickness,
        snapshot.config.voxelResolution
    );
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

    deformation::Handle fixed =
        deformation::makeFixedSupportHandle(domain, fixedPoints);
    std::set<int> occupiedNodes(
        fixed.constrainedNodes.begin(),
        fixed.constrainedNodes.end()
    );
    deformation::HandleVector handles;
    handles.reserve(snapshot.userHandles.size() + 1U);
    handles.push_back(std::move(fixed));
    for (const auto& specification : snapshot.userHandles)
    {
        static_cast<void>(domain.embed(specification.pivot));
        const int node = domain.nearestActiveNode(specification.pivot);
        if (!occupiedNodes.insert(node).second)
        {
            throw std::domain_error(
                "user handle conflicts with a fixed or previous handle node"
            );
        }
        deformation::Handle handle;
        handle.restPivot = specification.pivot;
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
    const deformation::DeformationSystem deformationSystem(
        grid,
        snapshot.vertices,
        snapshot.faces,
        binding,
        snapshot.config.bbw
    );
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
    const stability::BalanceTargetSchedule schedule(targets);
    const optimization::AlternatingBalanceOptimizer optimizer(
        deformationSystem,
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
