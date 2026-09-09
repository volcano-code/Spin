#include "api/SpinIt2014.h"
#include "spin2014/SlqpAuditTrace.h"

#include "manufacturing/ManufacturingGeometry.h"

#include <Eigen/Eigenvalues>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>


namespace
{


std::string objText(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(17);
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        output << "v " << vertices(row, 0) << ' ' << vertices(row, 1) << ' '
               << vertices(row, 2) << '\n';
    }
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        output << "f " << faces(row, 0) + 1 << ' ' << faces(row, 1) + 1 << ' '
               << faces(row, 2) + 1 << '\n';
    }
    return output.str();
}


void writeText(const std::string& filename, const std::string& text)
{
    const std::filesystem::path target(filename);
    if (target.has_parent_path())
    {
        std::filesystem::create_directories(target.parent_path());
    }
    const std::filesystem::path temporary = target.string() + ".tmp";
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output || !(output << text))
        {
            throw std::runtime_error("failed to write Spin-It 2014 artifact");
        }
    }
    std::error_code error;
    std::filesystem::rename(temporary, target, error);
    if (error)
    {
        std::filesystem::remove(temporary);
        throw std::runtime_error("failed to publish Spin-It 2014 artifact");
    }
}


void vectorJson(std::ostream& output, const Eigen::VectorXd& value)
{
    output << '[' << std::setprecision(17);
    for (Eigen::Index index = 0; index < value.size(); ++index)
    {
        if (index != 0) output << ',';
        output << value(index);
    }
    output << ']';
}


void matrixJson(std::ostream& output, const Eigen::Matrix3d& value)
{
    output << '[';
    for (int row = 0; row < 3; ++row)
    {
        if (row != 0) output << ',';
        vectorJson(output, value.row(row).transpose());
    }
    output << ']';
}


std::string diagnostics(
    const spin2014::AdaptiveVoxelSpinResult& result,
    const spin2014::AdaptiveVoxelSpinRequest& request,
    api::SpinIt2014Evidence evidence,
    bool sealed
)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    const bool massPropertiesAvailable =
        result.finalMoments.isFinite() && result.finalMoments.mass() > 0.0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    Eigen::Matrix3d principalAxes = Eigen::Matrix3d::Identity();
    if (massPropertiesAvailable)
    {
        center = result.finalMoments.centerOfMass();
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(
            result.finalMoments.inertiaAboutCenter()
        );
        if (eigen.info() == Eigen::Success)
        {
            principalAxes = eigen.eigenvectors();
        }
    }
    output << "{\n  \"schema\": \"spin-it-2014-artifact-v1\",\n"
           << "  \"schema_version\": 1,\n"
           << "  \"evidence\": \""
           << api::spinIt2014EvidenceToken(evidence) << "\",\n"
           << "  \"status\": \""
           << spin2014::adaptiveSpinStatusToken(result.status) << "\",\n"
           << "  \"mode\": \"" << spin::spinModeToken(request.mode)
           << "\",\n  \"objective\": \""
           << spin2014::spin2014ObjectiveToken(request.objective.value_or(
                  spin2014::defaultSpin2014Objective(request.mode)
              ))
           << "\",\n  \"objective_resolution\": \""
           << (request.objective.has_value() ? "EXPLICIT" : "MODE_DEFAULT_LEGACY")
           << "\",\n  \"objective_semantics\": \""
           << (request.formulation == spin2014::Spin2014Formulation::PaperGivens ?
               "PAPER_GIVENS_RAW_MOMENTS" : "ANGLE_ELIMINATED_TRANSVERSE_INERTIA")
           << "\",\n"
           << "  \"original_phi_trajectory\": false,\n"
           << "  \"inertia_weight\": " << std::setprecision(17)
           << request.inertiaWeight
           << ",\n  \"center_of_mass_weight\": "
           << request.centerOfMassWeight
           << ",\n  \"geometry_sealed\": "
           << (sealed ? "true" : "false") << ",\n"
           << "  \"input_hash\": \"" << result.inputHash << "\",\n"
           << "  \"trace_hash\": \"" << result.adaptiveTraceHash << "\",\n"
           << "  \"geometry_hash\": \"" << result.geometryHash << "\",\n"
           << "  \"output_hash\": \"" << result.outputHash << "\",\n"
           << "  \"binary_constraint_infinity_norm\": "
           << std::setprecision(17) << result.binaryConstraintInfinityNorm
           << ",\n  \"binary_projection_changed_leaves\": "
           << result.binaryProjectionChangedLeaves
           << ",\n  \"binary_projection_max_delta\": "
           << result.binaryProjectionMaxDelta
           << ",\n  \"output_grid_dimensions\": ["
           << result.outputGridDimensions(0) << ','
           << result.outputGridDimensions(1) << ','
           << result.outputGridDimensions(2) << ']'
           << ",\n  \"mass_properties_available\": "
           << (massPropertiesAvailable ? "true" : "false")
           << ",\n  \"center_of_mass_model\": ";
    if (massPropertiesAvailable) vectorJson(output, center);
    else output << "null";
    output << ",\n  \"target_origin_model\": [0,0,0],\n"
           << "  \"target_axes_model\": ";
    matrixJson(output, Eigen::Matrix3d::Identity());
    output << ",\n  \"principal_axes_model\": ";
    if (massPropertiesAvailable) matrixJson(output, principalAxes);
    else output << "null";
    output << ",\n  \"final_raw_moments\": ";
    vectorJson(output, result.finalMoments.values());
    if(request.formulation == spin2014::Spin2014Formulation::PaperGivens) {
        output << ",\n  \"formulation\": \"PAPER_GIVENS\",\n  \"initial_phi\": " << result.initialPhi
               << ",\n  \"final_phi\": " << result.finalPhi
               << ",\n  \"joint_angle_optimization\": " << (request.mode == spin::SpinMode::StaticBalance ? "false" : "true")
               << ",\n  \"binary_constraints\": ";
        vectorJson(output,result.binaryConstraints);
    }
    if(request.slqp.algorithm==spin2014::SlqpAlgorithm::AuditedKktRestorationV1) {
        output<<",\n  \"slqp_algorithm\":\"AUDITED_KKT_RESTORATION_V1\",\n  \"optimality_scope\":\"FIRST_ORDER_KKT_ONLY\",\n  \"slqp_audit_controls\":";
        spin2014::writeSlqpAuditControlsJson(output,request.slqp);
    }
    output << ",\n  \"adaptive_passes\": [";
    for (std::size_t index = 0; index < result.adaptiveTrace.size(); ++index)
    {
        if (index != 0U) output << ',';
        const spin2014::AdaptiveSpinPass& pass = result.adaptiveTrace[index];
        output << "{\"pass\":" << pass.pass
               << ",\"leaves\":" << pass.leafCountBefore
               << ",\"fractional\":" << pass.fractionalLeafCount
               << ",\"split\":" << pass.splitLeafCount
               << ",\"merged\":" << pass.mergedParentCount
               << ",\"slqp_status\":\""
               << spin2014::slqpStatusToken(pass.slqpStatus) << "\""
               << ",\"final_equality_infinity_norm\":"
               << pass.finalEqualityInfinityNorm
               << ",\"final_projected_gradient_infinity_norm\":"
               << pass.finalProjectedGradientInfinityNorm
               << ",\"line_search_evaluations\":"
               << pass.totalLineSearchEvaluations
               << ",\"active_set_transitions\":"
               << pass.activeSetTransitions
               << ",\"partition_hash\":\"" << pass.partitionHash
               << "\",\"slqp_trace_hash\":\"" << pass.slqpTraceHash << "\"";
        if(request.formulation == spin2014::Spin2014Formulation::PaperGivens) {
            output << ",\"initial_phi\":" << pass.initialPhi << ",\"final_phi\":" << pass.finalPhi << ",\"slqp_trace\":[";
            for(std::size_t j=0;j<pass.slqpTrace.size();++j) {
                if(j)output<<',';const auto& t=pass.slqpTrace[j];
                output << "{\"iteration\":"<<t.iteration<<",\"phi\":"<<t.phi
                       <<",\"objective\":"<<t.objective<<",\"equality_infinity_norm\":"<<t.equalityInfinityNorm
                       <<",\"accepted_step\":"<<t.acceptedStep;
                if(request.slqp.algorithm==spin2014::SlqpAlgorithm::AuditedKktRestorationV1)spin2014::writeSlqpAuditJson(output,t);
                output<<"}";
            }output<<']';
        }
        output << "}";
    }
    output << "]\n}\n";
    return output.str();
}


} // namespace


namespace api
{


const char* spinIt2014EvidenceToken(SpinIt2014Evidence evidence) noexcept
{
    switch (evidence)
    {
    case SpinIt2014Evidence::PaperFormulaReproduction:
        return "PAPER_FORMULA_REPRODUCTION";
    case SpinIt2014Evidence::ReconstructedAsset:
        return "RECONSTRUCTED_ASSET";
    }
    return "UNKNOWN";
}


struct SpinIt2014Artifact::Storage
{
    spin2014::AdaptiveVoxelSpinResult result;
    spin2014::AdaptiveVoxelSpinRequest request;
    SpinIt2014Evidence evidence = SpinIt2014Evidence::PaperFormulaReproduction;
    bool sealed = false;
    std::string diagnostics;
};


SpinIt2014Artifact::SpinIt2014Artifact(std::shared_ptr<const Storage> storage)
    : storage_(std::move(storage))
{
    if (!storage_) throw std::invalid_argument("Spin-It 2014 storage must exist");
}


SpinIt2014Artifact::~SpinIt2014Artifact() = default;
std::string SpinIt2014Artifact::statusToken() const
{
    return spin2014::adaptiveSpinStatusToken(storage_->result.status);
}
std::string SpinIt2014Artifact::evidenceToken() const
{
    return spinIt2014EvidenceToken(storage_->evidence);
}
bool SpinIt2014Artifact::succeeded() const
{
    return storage_->result.status == spin2014::AdaptiveSpinStatus::Converged;
}
bool SpinIt2014Artifact::geometrySealed() const { return storage_->sealed; }
const spin2014::AdaptiveVoxelSpinResult& SpinIt2014Artifact::result() const
{
    return storage_->result;
}
const Eigen::MatrixXd& SpinIt2014Artifact::outerVertices() const
{
    return storage_->result.sourceBoundary.outerVertices;
}
const Eigen::MatrixXi& SpinIt2014Artifact::outerFaces() const
{
    return storage_->result.sourceBoundary.outerFaces;
}
const Eigen::MatrixXd& SpinIt2014Artifact::designVertices() const
{
    return storage_->result.designBoundary.outerVertices;
}
const Eigen::MatrixXi& SpinIt2014Artifact::designFaces() const
{
    return storage_->result.designBoundary.outerFaces;
}
const Eigen::MatrixXd& SpinIt2014Artifact::materialVertices() const
{
    return storage_->result.materialBoundary.outerVertices;
}
const Eigen::MatrixXi& SpinIt2014Artifact::materialFaces() const
{
    return storage_->result.materialBoundary.outerFaces;
}
const Eigen::MatrixXd& SpinIt2014Artifact::cavityVertices() const
{
    return storage_->result.materialBoundary.innerVertices;
}
const Eigen::MatrixXi& SpinIt2014Artifact::cavityFaces() const
{
    return storage_->result.materialBoundary.innerFaces;
}
const std::string& SpinIt2014Artifact::diagnosticsJson() const
{
    return storage_->diagnostics;
}


void SpinIt2014Artifact::saveOuterMesh(const std::string& filename) const
{
    writeText(filename, objText(outerVertices(), outerFaces()));
}
void SpinIt2014Artifact::saveDesignMesh(const std::string& filename) const
{
    writeText(filename, objText(designVertices(), designFaces()));
}
void SpinIt2014Artifact::saveMaterialMesh(const std::string& filename) const
{
    writeText(filename, objText(materialVertices(), materialFaces()));
}
void SpinIt2014Artifact::saveCavityMesh(const std::string& filename) const
{
    writeText(filename, objText(cavityVertices(), cavityFaces()));
}
void SpinIt2014Artifact::saveDiagnostics(const std::string& filename) const
{
    writeText(filename, diagnosticsJson());
}


SpinIt2014Artifact solveSpinIt2014(
    const spin2014::AdaptiveVoxelSpinRequest& request,
    SpinIt2014Evidence evidence
)
{
    auto storage = std::make_shared<SpinIt2014Artifact::Storage>();
    storage->request = request;
    storage->evidence = evidence;
    storage->result = spin2014::solveAdaptiveVoxelSpin(request);
    if (storage->result.status == spin2014::AdaptiveSpinStatus::Converged)
    {
        const manufacturing::BoundaryAudit outer =
            manufacturing::auditTriangleBoundary(
                storage->result.materialBoundary.outerVertices,
                storage->result.materialBoundary.outerFaces,
                manufacturing::BoundaryOrientation::OutwardPositive
            );
        const manufacturing::BoundaryAudit cavity =
            manufacturing::auditTriangleBoundary(
                storage->result.materialBoundary.innerVertices,
                storage->result.materialBoundary.innerFaces,
                manufacturing::BoundaryOrientation::InwardNegative,
                true
            );
        storage->sealed = outer.strictBoundaryValid() &&
            (cavity.empty || cavity.strictBoundaryValid()) &&
            storage->result.inputHash.size() == 64U &&
            storage->result.adaptiveTraceHash.size() == 64U &&
            storage->result.geometryHash.size() == 64U &&
            storage->result.outputHash.size() == 64U;
    }
    storage->diagnostics = diagnostics(
        storage->result,
        storage->request,
        evidence,
        storage->sealed
    );
    return SpinIt2014Artifact(std::move(storage));
}


} // namespace api
