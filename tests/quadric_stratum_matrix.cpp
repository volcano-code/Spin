#include "spin/CanonicalFullRankHollowing.h"
#include "spin/QuadricSourceMotion.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Coefficients = spin::QuadricCoefficients10::Vector;

struct Fixture {
    const char* name;
    const char* role;
    Coefficients coefficients;
    bool expectedDifferential;
    const char* expectedSubtype;
};

struct CallRecord {
    bool returned = false;
    std::string exceptionType;
    std::string exceptionMessage;
};

struct CaseRecord {
    Fixture fixture;
    Coefficients normalizedCoefficients = Coefficients::Zero();
    CallRecord classificationCall;
    std::string classificationToken;
    std::string classificationSubtype;
    CallRecord sourceMotionCall;
    bool sourceMotionDefined = false;
    std::string sourceMotionDiagnostic;
    std::string sourceMotionScientificHash;
    int sourceVertexCount = 0;
    int sourceFaceCount = 0;
    CallRecord materialCall;
    bool materialMomentsFinite = false;
    double materialMass = 0.0;
    double materialMomentNorm = 0.0;
    std::string materialGeometryHash;
    std::string materialScientificHash;
    std::string materialDiagnostic;
    bool jacobianAvailable = false;
    bool jacobianFinite = false;
    double jacobianNorm = 0.0;
    CallRecord solverCall;
    std::string solverStatus;
    std::string solverDiagnostic;
    int solverFunctionEvaluations = 0;
    bool solverFinalGeometryVerified = false;
    bool solverFinalEvaluationAvailable = false;
    std::string solverVerificationDiagnostic;
    std::vector<std::string> violations;
};

spin::TriangleMesh box(double lower, double upper)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.faces.resize(12, 3);
    mesh.vertices <<
        lower, lower, lower,
        upper, lower, lower,
        upper, upper, lower,
        lower, upper, lower,
        lower, lower, upper,
        upper, lower, upper,
        upper, upper, upper,
        lower, upper, upper;
    mesh.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
    return mesh;
}

spin::SpinDesignRequest request()
{
    spin::SpinDesignRequest value;
    value.designDomain = box(-1.0, 1.0);
    value.outerSurface = box(-2.0, 2.0);
    value.materialDensity = 1.7;
    value.minimumWallThickness = 0.2;
    value.settings.equationProfile = spin::SpinEquationProfile::DualMonomialRawV1;
    value.settings.maximumRestarts = 1;
    value.settings.maximumFunctionEvaluationsPerRestart = 2;
    value.settings.equationTolerance = 1e-12;
    value.settings.physicalRelativeTolerance = 1e-6;
    return value;
}

spin::QuadricMomentEvaluatorSettings evaluatorSettings()
{
    spin::QuadricMomentEvaluatorSettings value;
    value.triangulation.icosphereSubdivisions = 2;
    value.triangulation.truncationRadius = 8.0;
    return value;
}

Coefficients coefficients(
    double a, double x, double y, double z,
    double x2, double y2, double z2, double yz, double xz, double xy)
{
    Coefficients value;
    value << a, x, y, z, x2, y2, z2, yz, xz, xy;
    return value;
}

std::vector<Fixture> fixtures()
{
    return {
        {"plane", "unsupported", coefficients(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), false, ""},
        {"elliptic-cylinder", "unsupported", coefficients(-0.25, 0, 0, 0, 1, 1, 0, 0, 0, 0), false, ""},
        {"elliptic-paraboloid", "unsupported", coefficients(0, 0, 0, -1, 1, 1, 0, 0, 0, 0), false, ""},
        {"hyperbolic-paraboloid", "unsupported", coefficients(0, 0, 0, -1, 1, -1, 0, 0, 0, 0), false, ""},
        {"parabolic-cylinder", "unsupported", coefficients(0, 0, -1, 0, 1, 0, 0, 0, 0, 0), false, ""},
        {"cone", "unsupported", coefficients(0, 0, 0, 0, 1, 1, -1, 0, 0, 0), false, ""},
        {"parallel-planes", "unsupported", coefficients(-0.25, 0, 0, 0, 1, 0, 0, 0, 0, 0), false, ""},
        {"intersecting-planes", "unsupported", coefficients(0, 0, 0, 0, 1, -1, 0, 0, 0, 0), false, ""},
        {"double-plane", "unsupported", coefficients(0, 0, 0, 0, 1, 0, 0, 0, 0, 0), false, ""},
        {"isolated-point", "unsupported", coefficients(0, 0, 0, 0, 1, 1, 1, 0, 0, 0), false, ""},
        {"imaginary-ellipsoid", "unsupported", coefficients(1, 0, 0, 0, 1, 1, 1, 0, 0, 0), false, ""},
        {"ellipsoid-control", "positive-control", coefficients(-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, 0.72, 0.06, -0.04, 0.03), true, "ellipsoid"},
        {"one-sheet-hyperboloid-control", "positive-control", coefficients(-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03), true, "one-sheet-hyperboloid"},
        {"two-sheet-hyperboloid-control", "positive-control", coefficients(0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03), true, "two-sheet-hyperboloid"}
    };
}

const char* semanticSubtype(spin::QuadricSurfaceSubtype subtype)
{
    switch (subtype) {
    case spin::QuadricSurfaceSubtype::Plane: return "plane";
    case spin::QuadricSurfaceSubtype::EllipticCylinder: return "elliptic-cylinder";
    case spin::QuadricSurfaceSubtype::Ellipsoid: return "ellipsoid";
    case spin::QuadricSurfaceSubtype::OneSheetHyperboloid: return "one-sheet-hyperboloid";
    case spin::QuadricSurfaceSubtype::TwoSheetHyperboloid: return "two-sheet-hyperboloid";
    case spin::QuadricSurfaceSubtype::Paraboloid: return "paraboloid";
    }
    throw std::runtime_error("unhandled quadric surface subtype");
}

template <class Function>
void capture(CallRecord& call, Function&& function)
{
    try {
        function();
        call.returned = true;
    } catch (const std::invalid_argument& error) {
        call.exceptionType = "std::invalid_argument";
        call.exceptionMessage = error.what();
    } catch (const std::runtime_error& error) {
        call.exceptionType = "std::runtime_error";
        call.exceptionMessage = error.what();
    }
}

std::string escapeJson(const std::string& value)
{
    std::ostringstream out;
    for (const unsigned char character : value) {
        switch (character) {
        case '\"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (character < 0x20) {
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<unsigned int>(character) << std::dec << std::setfill(' ');
            } else {
                out << character;
            }
        }
    }
    return out.str();
}

void writeString(std::ostream& out, const std::string& value)
{
    out << '\"' << escapeJson(value) << '\"';
}

void writeCall(std::ostream& out, const CallRecord& call)
{
    out << "{\"returned\":" << (call.returned ? "true" : "false")
        << ",\"exception_type\":";
    writeString(out, call.exceptionType);
    out << ",\"exception_message\":";
    writeString(out, call.exceptionMessage);
    out << '}';
}

void writeVector(std::ostream& out, const Coefficients& value)
{
    out << '[';
    for (int index = 0; index < value.size(); ++index) {
        if (index != 0) out << ',';
        out << value(index);
    }
    out << ']';
}

CaseRecord evaluateCase(
    const Fixture& fixture,
    const spin::SpinDesignRequest& design,
    const spin::QuadricMomentEvaluatorSettings& settings)
{
    CaseRecord record;
    record.fixture = fixture;

    if (!fixture.coefficients.allFinite() || fixture.coefficients.norm() == 0.0)
        throw std::runtime_error(std::string("malformed fixture: ") + fixture.name);
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(fixture.coefficients);
    record.normalizedCoefficients = quadric.values();

    spin::QuadricSurfaceSubtype classifiedSubtype = spin::QuadricSurfaceSubtype::Plane;
    capture(record.classificationCall, [&] {
        classifiedSubtype = spin::classifyQuadricSurfaceSubtype(quadric);
    });
    if (record.classificationCall.returned) {
        const auto subtype = classifiedSubtype;
        record.classificationToken = spin::quadricSurfaceSubtypeToken(subtype);
        record.classificationSubtype = semanticSubtype(subtype);
    }

    Eigen::MatrixXd domainBounds(2, 3);
    domainBounds << -1, -1, -1, 1, 1, 1;
    capture(record.sourceMotionCall, [&] {
        const auto source = spin::evaluateQuadricSourceMotion(
            quadric, domainBounds, settings.triangulation);
        record.sourceMotionDefined = source.defined;
        record.sourceMotionDiagnostic = source.diagnostic;
        record.sourceMotionScientificHash = source.scientificHash;
        record.sourceVertexCount = static_cast<int>(source.triangulation.closedBoundary.vertices.rows());
        record.sourceFaceCount = static_cast<int>(source.triangulation.closedBoundary.faces.rows());
    });

    capture(record.materialCall, [&] {
        const auto material = spin::evaluateCanonicalFullRankMaterialMoments(
            design, quadric, settings);
        record.materialMomentsFinite = material.totalMomentsTarget.values().allFinite();
        record.materialMass = material.totalMomentsTarget.values()(0);
        record.materialMomentNorm = material.totalMomentsTarget.values().norm();
        record.materialGeometryHash = material.geometryHash;
        record.materialScientificHash = material.scientificHash;
        record.materialDiagnostic = material.diagnostic;
        record.jacobianAvailable = material.coefficientJacobian.has_value();
        if (material.coefficientJacobian) {
            record.jacobianFinite = material.coefficientJacobian->allFinite();
            record.jacobianNorm = material.coefficientJacobian->norm();
        }
    });

    capture(record.solverCall, [&] {
        const auto solver = spin::solveCanonicalFullRankQuadricHollowing(
            design, settings, {quadric});
        record.solverStatus = spin::spinSolverStatusToken(solver.numerical.status);
        record.solverDiagnostic = solver.numerical.diagnostic;
        record.solverFunctionEvaluations = solver.numerical.functionEvaluations;
        record.solverFinalGeometryVerified = solver.finalGeometryVerified;
        record.solverFinalEvaluationAvailable = solver.finalEvaluation.has_value();
        record.solverVerificationDiagnostic = solver.verificationDiagnostic;
    });

    if (record.solverCall.returned && record.solverFunctionEvaluations > 2)
        record.violations.push_back("solver exceeded the preregistered two-evaluation limit");

    if (!fixture.expectedDifferential) {
        if (record.sourceMotionCall.returned && record.sourceMotionDefined)
            record.violations.push_back("unsupported stratum exposed defined source motion");
        if (record.materialCall.returned && record.jacobianAvailable)
            record.violations.push_back("unsupported differential became available");
        if (record.solverCall.returned &&
            (record.solverStatus == "CONVERGED" || record.solverStatus == "INFEASIBLE"))
            record.violations.push_back("unsupported stratum received a mathematical solver verdict");
        if (record.solverCall.returned &&
            (record.solverFinalGeometryVerified || record.solverFinalEvaluationAvailable))
            record.violations.push_back("unsupported stratum exposed a certified/final solver artifact");
    } else {
        if (!record.classificationCall.returned)
            record.violations.push_back("positive control classification threw");
        else if (record.classificationSubtype != fixture.expectedSubtype)
            record.violations.push_back("positive control classification subtype changed");
        if (!record.sourceMotionCall.returned || !record.sourceMotionDefined)
            record.violations.push_back("positive control source motion unavailable");
        if (!record.materialCall.returned)
            record.violations.push_back("positive control material evaluation threw");
        else {
            if (!record.materialMomentsFinite || !std::isfinite(record.materialMass) || record.materialMass <= 0.0)
                record.violations.push_back("positive control native material values are not finite positive material");
            if (!record.jacobianAvailable)
                record.violations.push_back("positive control native Jacobian unavailable");
            else if (!record.jacobianFinite || !std::isfinite(record.jacobianNorm) || record.jacobianNorm <= 0.0)
                record.violations.push_back("positive control native Jacobian is not finite nonzero");
        }
        if (!record.solverCall.returned)
            record.violations.push_back("positive control solver call threw instead of returning a bounded outcome");
    }
    return record;
}

void writeResult(const std::string& path, const std::vector<CaseRecord>& records)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open matrix JSON output: " + path);
    out << std::setprecision(17);
    bool passed = true;
    for (const auto& record : records) passed = passed && record.violations.empty();
    out << "{\n  \"schema\":\"m1-stratum-matrix-result-v1\",\n"
        << "  \"preregistration\":\"stratum/preregistration.json\",\n"
        << "  \"scope_limitation\":\"finite readiness characterization; not a proof over every rank transition\",\n"
        << "  \"negative_control_interpretation\":\"harness sensitivity control only; not production mutation coverage\",\n"
        << "  \"overall_pass\":" << (passed ? "true" : "false") << ",\n"
        << "  \"settings\":{\"design_box\":[-1,1],\"outer_box\":[-2,2],\"density\":1.7,\"wall\":0.2,"
        << "\"subdivisions\":2,\"truncation\":8,\"profile\":\"dual-monomial-raw-v1\","
        << "\"maximum_restarts\":1,\"maximum_numerical_evaluations_per_restart\":2},\n"
        << "  \"cases\":[\n";
    for (std::size_t index = 0; index < records.size(); ++index) {
        const auto& record = records[index];
        out << "    {\"name\":";
        writeString(out, record.fixture.name);
        out << ",\"role\":";
        writeString(out, record.fixture.role);
        out << ",\"literal_coefficients\":";
        writeVector(out, record.fixture.coefficients);
        out << ",\"normalized_coefficients\":";
        writeVector(out, record.normalizedCoefficients);
        out << ",\"expected_differential_admissible\":"
            << (record.fixture.expectedDifferential ? "true" : "false")
            << ",\"classification\":{\"call\":";
        writeCall(out, record.classificationCall);
        out << ",\"token\":";
        writeString(out, record.classificationToken);
        out << ",\"semantic_subtype\":";
        writeString(out, record.classificationSubtype);
        out << "},\"source_motion\":{\"call\":";
        writeCall(out, record.sourceMotionCall);
        out << ",\"defined\":" << (record.sourceMotionDefined ? "true" : "false")
            << ",\"diagnostic\":";
        writeString(out, record.sourceMotionDiagnostic);
        out << ",\"scientific_hash\":";
        writeString(out, record.sourceMotionScientificHash);
        out << ",\"vertices\":" << record.sourceVertexCount
            << ",\"faces\":" << record.sourceFaceCount
            << "},\"material_value_and_j\":{\"call\":";
        writeCall(out, record.materialCall);
        out << ",\"moments_finite\":" << (record.materialMomentsFinite ? "true" : "false")
            << ",\"mass\":" << record.materialMass
            << ",\"moment_norm\":" << record.materialMomentNorm
            << ",\"geometry_hash\":";
        writeString(out, record.materialGeometryHash);
        out << ",\"scientific_hash\":";
        writeString(out, record.materialScientificHash);
        out << ",\"diagnostic\":";
        writeString(out, record.materialDiagnostic);
        out << ",\"jacobian_available\":" << (record.jacobianAvailable ? "true" : "false")
            << ",\"jacobian_finite\":" << (record.jacobianFinite ? "true" : "false")
            << ",\"jacobian_norm\":" << record.jacobianNorm
            << "},\"solver\":{\"call\":";
        writeCall(out, record.solverCall);
        out << ",\"status\":";
        writeString(out, record.solverStatus);
        out << ",\"diagnostic\":";
        writeString(out, record.solverDiagnostic);
        out << ",\"function_evaluations\":" << record.solverFunctionEvaluations
            << ",\"final_geometry_verified\":" << (record.solverFinalGeometryVerified ? "true" : "false")
            << ",\"final_evaluation_available\":" << (record.solverFinalEvaluationAvailable ? "true" : "false")
            << ",\"verification_diagnostic\":";
        writeString(out, record.solverVerificationDiagnostic);
        out << "},\"violations\":[";
        for (std::size_t violation = 0; violation < record.violations.size(); ++violation) {
            if (violation != 0) out << ',';
            writeString(out, record.violations[violation]);
        }
        out << "]}" << (index + 1 == records.size() ? "\n" : ",\n");
    }
    out << "  ]\n}\n";
    if (!out) throw std::runtime_error("failed writing matrix JSON output: " + path);
}

int demandUnsupportedJacobian()
{
    const auto design = request();
    const auto settings = evaluatorSettings();
    spin::validateSpinDesignRequest(design);
    spin::validateQuadricMomentEvaluatorSettings(settings);
    const auto fixture = fixtures().front();
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(fixture.coefficients);
    const auto material = spin::evaluateCanonicalFullRankMaterialMoments(design, quadric, settings);
    if (!material.coefficientJacobian)
        throw std::runtime_error("demand-unsupported-j: plane Jacobian unavailable: " + material.diagnostic);
    std::cout << "[NEGATIVE-CONTROL-UNEXPECTED-PASS] plane Jacobian norm="
              << material.coefficientJacobian->norm() << '\n';
    return 0;
}

int runMatrix(const std::string& outputPath)
{
    const auto design = request();
    const auto settings = evaluatorSettings();
    spin::validateSpinDesignRequest(design);
    spin::validateQuadricMomentEvaluatorSettings(settings);

    std::vector<CaseRecord> records;
    for (const auto& fixture : fixtures()) {
        records.push_back(evaluateCase(fixture, design, settings));
        std::cout << "case=" << fixture.name
                  << " material_returned=" << records.back().materialCall.returned
                  << " J=" << records.back().jacobianAvailable
                  << " solver_returned=" << records.back().solverCall.returned
                  << " violations=" << records.back().violations.size() << '\n';
    }
    writeResult(outputPath, records);

    bool passed = true;
    for (const auto& record : records) {
        for (const auto& violation : record.violations) {
            passed = false;
            std::cerr << "[VIOLATION] " << record.fixture.name << ": " << violation << '\n';
        }
    }
    if (!passed) return 1;
    std::cout << "[PASS] finite stratum matrix; result=" << outputPath << '\n';
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::cout << std::setprecision(17);
        if (argc == 3 && std::string(argv[1]) == "--mode" &&
            std::string(argv[2]) == "demand-unsupported-j")
            return demandUnsupportedJacobian();
        if (argc != 2)
            throw std::invalid_argument(
                "usage: quadric_stratum_matrix OUTPUT.json | --mode demand-unsupported-j");
        return runMatrix(argv[1]);
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
