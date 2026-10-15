#include "spin/FullRankHollowingSolver.h"
#include "spin/QuadricInverseGaussMap.h"
#include "spin/QuadricSourceMotion.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

constexpr int kPreRegisteredMaximumRestarts = 2;
constexpr int kMaximumEvaluationsPerRestart = 60;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string hexEncode(const std::string& value) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string encoded;
    encoded.reserve(value.size() * 2);
    for (const unsigned char byte : value) {
        encoded.push_back(digits[byte >> 4]);
        encoded.push_back(digits[byte & 0x0f]);
    }
    return encoded;
}

spin::TriangleMesh box(double lo, double hi) {
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices << lo, lo, lo, hi, lo, lo, hi, hi, lo, lo, hi, lo,
                    lo, lo, hi, hi, lo, hi, hi, hi, hi, lo, hi, hi;
    mesh.faces.resize(12, 3);
    mesh.faces << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
                  0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2,
                  0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5;
    return mesh;
}

spin::SpinDesignRequest fixedRequest(unsigned int seed) {
    spin::SpinDesignRequest request;
    request.outerSurface = box(-2.0, 2.0);
    request.designDomain = box(-1.0, 1.0);
    request.mode = spin::SpinMode::YoYo;
    request.materialDensity = 1.7;
    request.minimumWallThickness = 0.0;
    request.settings.equationProfile = spin::SpinEquationProfile::DualMonomialRawV1;
    request.settings.maximumIterations = 60;
    request.settings.maximumRestarts = kPreRegisteredMaximumRestarts;
    request.settings.maximumFunctionEvaluationsPerRestart =
        kMaximumEvaluationsPerRestart;
    request.settings.deterministicSeed = seed;

    const auto outer = physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(request.outerSurface.vertices,
                                  request.outerSurface.faces,
                                  request.materialDensity));
    const auto design = physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(request.designDomain.vertices,
                                  request.designDomain.faces,
                                  request.materialDensity));
    request.fixedShellMoments = outer - design;
    return request;
}

spin::QuadricMomentEvaluatorSettings evaluatorSettings() {
    spin::QuadricMomentEvaluatorSettings settings;
    settings.triangulation.icosphereSubdivisions = 2;
    settings.triangulation.truncationRadius = 8.0;
    // The strict solver profile must override the source-compatible legacy
    // derivative setting; the output records the effective path per seed.
    settings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
    return settings;
}

void runUnknownSeedSweep(int selectedSeed = -1) {
    std::cout << std::setprecision(17);
    const auto settings = evaluatorSettings();
    for (unsigned int seed = 0; seed < 5; ++seed) {
        if (selectedSeed >= 0 && seed != static_cast<unsigned int>(selectedSeed))
            continue;
        const auto request = fixedRequest(seed);
        const auto generated = spin::deterministicInitialQuadrics(
            request.mode, request.settings);
        require(generated.size() == static_cast<std::size_t>(kPreRegisteredMaximumRestarts),
                "each seed must generate its preregistered deterministic restart starts");
        for (std::size_t index = 0; index < generated.size(); ++index) {
            std::cout << "initial-start seed=" << seed
                      << " restart=" << index
                      << " hash=" << generated[index].hash()
                      << " coefficients=";
            for (int coefficient = 0; coefficient < 10; ++coefficient) {
                if (coefficient) std::cout << ',';
                std::cout << generated[index].values()(coefficient);
            }
            std::cout << '\n';
        }
        const auto subtype = spin::classifyQuadricSurfaceSubtype(generated.front());
        const auto initialSourceMotion = spin::evaluateQuadricSourceMotion(
            generated.front(), request.designDomain.vertices, settings.triangulation);

        // Deliberately omit initialGuesses. The solver must generate its own
        // start only from the frozen request and registered seed.
        const auto result = spin::solveFullRankQuadricHollowing(request, settings);
        const spin::SpinSolverEvent* firstDerivativeFailure = nullptr;
        const spin::SpinSolverEvent* booleanAuditFailure = nullptr;
        for (const auto& event : result.numerical.events) {
            if (event.type == "evaluation" && event.outcome == "failed" &&
                event.detail.find(
                    "quadric material Boolean is not a positive closed manifold") !=
                    std::string::npos) {
                booleanAuditFailure = &event;
            }
            if (event.type == "evaluation" && event.outcome == "failed" &&
                event.detail.find("moment evaluator derivative is undefined") !=
                    std::string::npos) {
                firstDerivativeFailure = &event;
                break;
            }
        }
        std::string firstFailureSubtype = "NONE";
        std::string firstFailureDetail;
        std::string firstFailureCoefficients = "none";
        int witnessVertex = -1;
        int witnessNormal = -1;
        double witnessNearestDistance = std::numeric_limits<double>::infinity();
        double witnessRadiusDelta = std::numeric_limits<double>::infinity();
        double witnessImplicitResidual = std::numeric_limits<double>::infinity();
        std::string witnessPoint = "none";
        double firstFailureStep = -1.0;
        int firstFailureEvaluation = -1;
        int firstFailureIteration = -1;
        if (firstDerivativeFailure != nullptr) {
            firstFailureEvaluation = firstDerivativeFailure->evaluationIndex;
            firstFailureIteration = firstDerivativeFailure->iteration;
            firstFailureStep = firstDerivativeFailure->step.value_or(-1.0);
            firstFailureDetail = firstDerivativeFailure->detail;
            if (firstDerivativeFailure->quadricCoefficients) {
                spin::QuadricCoefficients10::Vector coefficients;
                for (int index = 0; index < coefficients.size(); ++index)
                    coefficients(index) = (*firstDerivativeFailure->quadricCoefficients)[
                        static_cast<std::size_t>(index)];
                std::ostringstream coefficientText;
                coefficientText << std::setprecision(17);
                for (int index = 0; index < coefficients.size(); ++index) {
                    if (index != 0) coefficientText << ',';
                    coefficientText << coefficients(index);
                }
                firstFailureCoefficients = coefficientText.str();
                const auto candidate = spin::QuadricCoefficients10::fromUnnormalized(
                    coefficients);
                firstFailureSubtype = spin::quadricSurfaceSubtypeToken(
                    spin::classifyQuadricSurfaceSubtype(candidate));

                auto triangulationSettings = settings.triangulation;
                triangulationSettings.retainFixedNormalProvenance = true;
                const auto source = spin::triangulateQuadric(
                    candidate, request.designDomain.vertices, triangulationSettings);
                auto referenceSettings = triangulationSettings;
                referenceSettings.truncationRadius = source.truncationRadius;
                referenceSettings.maximumRadiusDoublings = 0;
                spin::QuadricCoefficients10::Vector sphereParameters;
                sphereParameters << 1, 0, 0, 0, -1, -1, -1, 0, 0, 0;
                const auto referenceSphere = spin::triangulateQuadric(
                    spin::QuadricCoefficients10::fromUnnormalized(sphereParameters),
                    request.designDomain.vertices, referenceSettings);
                const spin::QuadricInverseGaussMap inverseGauss(candidate);
                for (int vertex = 0; vertex < source.interfaceVertexCount; ++vertex) {
                    const Eigen::Vector3d point =
                        source.closedBoundary.vertices.row(vertex).transpose();
                    double nearest = std::numeric_limits<double>::infinity();
                    int nearestIndex = -1;
                    for (int normal = 0; normal < referenceSphere.fixedNormals.rows(); ++normal) {
                        const auto mapped = inverseGauss.point(
                            referenceSphere.fixedNormals.row(normal));
                        if (mapped.valid()) {
                            const double distance = (point - *mapped.point).norm();
                            if (distance < nearest) {
                                nearest = distance;
                                nearestIndex = normal;
                            }
                        }
                    }
                    const double identityTolerance = 1e-12 * std::max(1.0, point.norm());
                    const double radiusDelta = std::abs(
                        point.norm() - source.truncationRadius);
                    const double rootTolerance =
                        1e-10 * std::max(1.0, source.truncationRadius);
                    if (nearest <= identityTolerance ||
                        (source.type == spin::QuadricSurfaceType::Hyperboloid &&
                         radiusDelta <= rootTolerance))
                        continue;
                    witnessVertex = vertex;
                    witnessNormal = nearestIndex;
                    witnessNearestDistance = nearest;
                    witnessRadiusDelta = radiusDelta;
                    witnessImplicitResidual = std::abs(candidate.evaluate(point));
                    std::ostringstream pointText;
                    pointText << std::setprecision(17) << point.x() << ','
                              << point.y() << ',' << point.z();
                    witnessPoint = pointText.str();
                    break;
                }
            }
        }
        require(result.inputHash.size() == 64 && result.outputHash.size() == 64,
                "seed run lacks bound input/output scientific identities");
        require(result.numerical.status != spin::SpinSolverStatus::InvalidRequest,
                "frozen unknown-answer request was unexpectedly invalid");
        require(result.numerical.restartCount <= kPreRegisteredMaximumRestarts &&
                    result.numerical.functionEvaluations <=
                        kPreRegisteredMaximumRestarts * kMaximumEvaluationsPerRestart &&
                    result.numerical.restartFunctionEvaluations.size() ==
                        static_cast<std::size_t>(result.numerical.restartCount) &&
                    std::all_of(result.numerical.restartFunctionEvaluations.begin(),
                                result.numerical.restartFunctionEvaluations.end(),
                                [](int count) {
                                    return count >= 0 &&
                                           count <= kMaximumEvaluationsPerRestart;
                                }) &&
                    std::accumulate(
                        result.numerical.restartFunctionEvaluations.begin(),
                        result.numerical.restartFunctionEvaluations.end(), 0) ==
                        result.numerical.functionEvaluations,
                "seed run exceeded its preregistered per-restart budget");
        require(result.numerical.equationProfile ==
                    spin::SpinEquationProfile::DualMonomialRawV1,
                "seed run escaped the strict equation profile");

        std::cout << "seed=" << seed
                  << " generated_start=" << generated.front().hash()
                  << " status=" << spin::spinSolverStatusToken(result.numerical.status)
                  << " evaluations=" << result.numerical.functionEvaluations
                  << " trace=" << result.numerical.trace.size()
                  << " restarts=" << result.numerical.restartCount
                  << " restart_evaluations=";
        for (std::size_t index = 0;
             index < result.numerical.restartFunctionEvaluations.size(); ++index) {
            if (index) std::cout << ',';
            std::cout << result.numerical.restartFunctionEvaluations[index];
        }
        std::cout
                  << " residual=" << result.numerical.normalizedEquationResidual
                  << " final_verified=" << result.finalGeometryVerified
                  << " input_hash=" << result.inputHash
                  << " event_hash=" << result.numerical.eventHash
                  << " trace_hash=" << result.numerical.traceHash
                  << " output_hash=" << result.outputHash
                  << " initial_subtype=" << spin::quadricSurfaceSubtypeToken(subtype)
                  << " initial_source_motion_defined=" << initialSourceMotion.defined
                  << " initial_source_diagnostic_hex=" << hexEncode(initialSourceMotion.diagnostic)
                  << " first_derivative_failure_evaluation=" << firstFailureEvaluation
                  << " first_derivative_failure_iteration=" << firstFailureIteration
                  << " first_derivative_failure_subtype=" << firstFailureSubtype
                  << " first_derivative_failure_step=" << firstFailureStep
                  << " first_derivative_failure_hex=" << hexEncode(firstFailureDetail)
                  << " first_derivative_failure_q=" << firstFailureCoefficients
                  << " source_witness_vertex=" << witnessVertex
                  << " source_witness_nearest_normal=" << witnessNormal
                  << " source_witness_nearest_distance=" << witnessNearestDistance
                  << " source_witness_radius_delta=" << witnessRadiusDelta
                  << " source_witness_q_residual=" << witnessImplicitResidual
                  << " source_witness_point=" << witnessPoint
                  << " solver_diagnostic_hex=" << hexEncode(result.numerical.diagnostic)
                  << " verification_diagnostic_hex=" << hexEncode(result.verificationDiagnostic)
                  << '\n';
        for (const auto& event : result.numerical.events) {
            if (event.type != "attempt_start" && event.type != "attempt_terminal" &&
                event.outcome != "failed" && event.outcome != "rejected")
                continue;
            std::cout << "solver-event seed=" << seed
                      << " restart=" << event.restart
                      << " iteration=" << event.iteration
                      << " evaluation=" << event.evaluationIndex
                      << " type=" << event.type
                      << " outcome=" << event.outcome
                      << " reason=" << event.reason
                      << " detail_hex=" << hexEncode(event.detail)
                      << " coefficients=";
            if (event.quadricCoefficients) {
                for (int coefficient = 0; coefficient < 10; ++coefficient) {
                    if (coefficient) std::cout << ',';
                    std::cout << (*event.quadricCoefficients)[
                        static_cast<std::size_t>(coefficient)];
                }
            } else {
                std::cout << "none";
            }
            std::cout << '\n';
        }
        if (booleanAuditFailure && booleanAuditFailure->quadricCoefficients) {
            std::cout << "boolean-audit-witness seed=" << seed << " eval="
                      << booleanAuditFailure->evaluationIndex << " q=";
            for (int index = 0; index < 10; ++index) {
                if (index != 0) std::cout << ',';
                std::cout << (*booleanAuditFailure->quadricCoefficients)[
                    static_cast<std::size_t>(index)];
            }
            std::cout << " detail_hex=" << hexEncode(booleanAuditFailure->detail)
                      << '\n';
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--diagnose-seed") {
            const int seed = std::stoi(argv[2]);
            require(seed >= 0 && seed < 5, "diagnostic seed is out of range");
            runUnknownSeedSweep(seed);
            return 0;
        }
        runUnknownSeedSweep();
        std::cout << "[PASS] unknown-seed-sweep\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] unknown-seed-sweep: " << error.what() << '\n';
        return 1;
    }
}
