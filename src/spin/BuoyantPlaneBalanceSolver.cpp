#include "spin/BuoyantPlaneBalanceSolver.h"

#include "spin/SpinDesignRequest.h"
#include "spin/SpinMomentProblem.h"

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace spin {
namespace {

void validateSettings(const BuoyantPlaneBalanceSettings& settings) {
    if (settings.maximumIterations < 0 || settings.maximumLineSearchSteps <= 0 ||
        !std::isfinite(settings.constraintTolerance) || settings.constraintTolerance <= 0.0 ||
        !std::isfinite(settings.backtrackingFactor) || settings.backtrackingFactor <= 0.0 ||
        settings.backtrackingFactor >= 1.0 || !std::isfinite(settings.minimumStep) ||
        settings.minimumStep <= 0.0)
        throw std::invalid_argument("invalid buoyant plane solver settings");
}

Eigen::Vector3d fixedConstraintScales(const SpinDesignRequest& request) {
    double length = 0.0;
    Eigen::Vector3d lower = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d upper = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
    for (Eigen::Index i = 0; i < request.designDomain.vertices.rows(); ++i) {
        const Eigen::Vector3d world = request.designDomain.vertices.row(i).transpose();
        const Eigen::Vector3d target = request.targetFrame.axes.transpose() *
            (world - request.targetFrame.origin);
        lower = lower.cwiseMin(target);
        upper = upper.cwiseMax(target);
    }
    length = (upper - lower).maxCoeff();
    if (!std::isfinite(length) || length <= 0.0)
        throw std::invalid_argument("buoyant plane fixed-input length scale is invalid");
    const double density = request.materialDensity;
    Eigen::Vector3d scales(density * length * length * length,
                           density * std::pow(length, 4),
                           density * std::pow(length, 4));
    if (!scales.allFinite() || (scales.array() <= 0.0).any())
        throw std::invalid_argument("buoyant plane fixed-input moment scales are not representable");
    return scales;
}

struct EvaluatedPlane {
    Eigen::Vector3d parameters = Eigen::Vector3d::Zero();
    Eigen::Vector4d coefficients = Eigen::Vector4d(0.0, 0.0, 0.0, -1.0);
    PlaneMomentEvaluation evaluation;
    Eigen::Vector3d constraints = Eigen::Vector3d::Zero();
    Eigen::Vector3d normalizedConstraints = Eigen::Vector3d::Zero();
    double residual = std::numeric_limits<double>::infinity();
};

EvaluatedPlane evaluate(const SpinDesignRequest& request,
                        const Eigen::Vector3d& parameters,
                        const Eigen::Vector3d& scales,
                        double targetMass,
                        const BuoyantPlaneBalanceSettings& settings) {
    EvaluatedPlane result;
    result.parameters = parameters;
    result.coefficients << parameters.x(), parameters.y(), parameters.z(), -1.0;
    result.evaluation = evaluatePlaneMaterialMoments(
        request, result.coefficients, settings.planeEvaluation);
    const auto& moments = result.evaluation.totalMomentsTarget.values();
    const double mass = moments[0];
    const Eigen::Vector2d center = request.buoyancy.centerOfBuoyancyTarget.head<2>();
    result.constraints << mass - targetMass,
        moments[1] - center.x() * mass,
        moments[2] - center.y() * mass;
    result.normalizedConstraints = result.constraints.cwiseQuotient(scales);
    result.residual = result.normalizedConstraints.norm();
    if (!result.constraints.allFinite() || !result.normalizedConstraints.allFinite() ||
        !std::isfinite(result.residual))
        throw std::domain_error("buoyant plane constraints are nonfinite");
    return result;
}

Eigen::Matrix3d reducedJacobian(const SpinDesignRequest& request,
                                const PlaneMomentEvaluation& evaluation,
                                const Eigen::Vector3d& scales,
                                double lengthScale) {
    if (!evaluation.coefficientJacobian)
        throw std::runtime_error("buoyant plane-family differential unavailable");
    const auto& jac = *evaluation.coefficientJacobian;
    const double cbx = request.buoyancy.centerOfBuoyancyTarget.x();
    const double cby = request.buoyancy.centerOfBuoyancyTarget.y();
    Eigen::Matrix3d result;
    for (int column = 0; column < 3; ++column) {
        // Internal variable zero is a/L, so its chain-rule factor is L.
        const double parameterScale = column == 0 ? lengthScale : 1.0;
        result(0, column) = jac(0, column) * parameterScale / scales(0);
        result(1, column) = (jac(1, column) - cbx * jac(0, column)) *
            parameterScale / scales(1);
        result(2, column) = (jac(2, column) - cby * jac(0, column)) *
            parameterScale / scales(2);
    }
    return result;
}

void append(BuoyantPlaneBalanceResult& result, const EvaluatedPlane& current,
            int iteration, double step, bool accepted, int evaluationIndex) {
    BuoyantPlaneBalanceTraceRow row;
    row.evaluation = evaluationIndex;
    row.iteration = iteration;
    row.planeCoefficients = current.coefficients;
    row.constraints = current.constraints;
    row.normalizedConstraints = current.normalizedConstraints;
    row.residualNorm = std::isfinite(current.residual) ? current.residual : -1.0;
    row.step = step;
    row.accepted = accepted;
    row.evaluationSucceeded = std::isfinite(current.residual);
    row.derivativeAvailable = row.evaluationSucceeded &&
        current.evaluation.coefficientJacobian.has_value();
    row.diagnostic = current.evaluation.diagnostic;
    result.trace.push_back(std::move(row));
}

void retain(BuoyantPlaneBalanceResult& result, const EvaluatedPlane& current,
            const SpinDesignRequest& request) {
    result.planeCoefficients = current.coefficients;
    result.constraints = current.constraints;
    result.normalizedConstraints = current.normalizedConstraints;
    result.totalMoments = current.evaluation.totalMomentsTarget;
    result.variableMoments = current.evaluation.variableMomentsTarget;
    result.evaluation = current.evaluation;
    if (result.totalMoments.mass() > 0.0) {
        const auto& moments = result.totalMoments.values();
        result.centerOfMassHeight = moments[3] / moments[0];
        result.metacenterHeight = buoyantMetacenterHeight(request.buoyancy);
        result.stabilityMargin = result.metacenterHeight - result.centerOfMassHeight;
    }
}

bool verifyConvergedPlaneFamily(BuoyantPlaneBalanceResult& result,
                                const EvaluatedPlane& current,
                                const SpinDesignRequest& request,
                                const BuoyantPlaneBalanceSettings& settings,
                                double lengthScale) {
    if (!current.evaluation.coefficientJacobian) {
        result.status = BuoyantPlaneBalanceStatus::DerivativeUnavailable;
        result.diagnostic = "equilibrium constraints reached without a usable plane-family differential";
        return false;
    }
    const Eigen::Matrix3d jacobian = reducedJacobian(
        request, current.evaluation, result.constraintScales, lengthScale);
    if (!jacobian.allFinite()) {
        result.status = BuoyantPlaneBalanceStatus::DerivativeUnavailable;
        result.diagnostic = "equilibrium constraints reached with a nonfinite plane-family differential";
        return false;
    }
    Eigen::FullPivLU<Eigen::Matrix3d> factor(jacobian);
    if (factor.rank() != 3 || factor.rcond() < 1e-13) {
        result.status = BuoyantPlaneBalanceStatus::SingularReducedSystem;
        result.diagnostic = "equilibrium constraints reached but the reduced plane system is singular";
        return false;
    }
    const auto& moments = result.totalMoments.values();
    const double mass = moments[0];
    const double centerHeight = moments[3] / mass;
    const double cbx = request.buoyancy.centerOfBuoyancyTarget.x();
    const double cby = request.buoyancy.centerOfBuoyancyTarget.y();
    const double bx = current.coefficients.y();
    const double by = current.coefficients.z();
    result.equalityMultipliers.resize(3);
    result.equalityMultipliers << bx / mass, by / mass,
        (current.coefficients.x() - centerHeight + cbx * bx + cby * by) / mass;
    result.multiplierMu = current.coefficients.norm() / mass;

    QuadricCoefficients10::Vector rawQuadric = QuadricCoefficients10::Vector::Zero();
    rawQuadric.head<4>() = current.coefficients;
    const auto normalizedQuadric = QuadricCoefficients10::fromUnnormalized(rawQuadric);
    const auto momentProblem = SpinMomentProblem::buoyantStability(
        result.targetTotalMass, request.buoyancy.centerOfBuoyancyTarget.head<2>());
    const auto momentEvaluation = momentProblem.evaluate(moments);
    const double auditTolerance = std::max(1e-12,
        settings.constraintTolerance * result.constraintScales.maxCoeff());
    result.theoremOneNumericalAudit = auditFasterTheoremOneNumericalConditions(
        normalizedQuadric, momentEvaluation, result.multiplierMu,
        result.equalityMultipliers, auditTolerance);
    if (result.theoremOneNumericalAudit->status !=
        FasterTheoremConditionStatus::SatisfiedNumerically) {
        result.status = BuoyantPlaneBalanceStatus::OptimalityAuditFailed;
        result.diagnostic = "reduced buoyant plane root failed the independent Eq. 7 numerical audit: " +
            result.theoremOneNumericalAudit->diagnostic;
        return false;
    }
    result.status = BuoyantPlaneBalanceStatus::Converged;
    result.equilibriumConstraintsVerified = true;
    result.eq7NumericalConditionsVerified = true;
    result.stableAtPrescribedWaterline = result.stabilityMargin > 0.0;
    result.diagnostic = "fixed-waterline equilibrium constraints and finite-precision Eq. 7 audit passed; this is not a theorem proof";
    return true;
}

} // namespace

const char* buoyantPlaneBalanceStatusToken(BuoyantPlaneBalanceStatus status) noexcept {
    switch (status) {
    case BuoyantPlaneBalanceStatus::Converged: return "CONVERGED";
    case BuoyantPlaneBalanceStatus::IterationLimit: return "ITERATION_LIMIT";
    case BuoyantPlaneBalanceStatus::DerivativeUnavailable: return "DERIVATIVE_UNAVAILABLE";
    case BuoyantPlaneBalanceStatus::SingularReducedSystem: return "SINGULAR_REDUCED_SYSTEM";
    case BuoyantPlaneBalanceStatus::LineSearchExhausted: return "LINE_SEARCH_EXHAUSTED";
    case BuoyantPlaneBalanceStatus::GeometryEvaluationFailed: return "GEOMETRY_EVALUATION_FAILED";
    case BuoyantPlaneBalanceStatus::OptimalityAuditFailed: return "OPTIMALITY_AUDIT_FAILED";
    case BuoyantPlaneBalanceStatus::InvalidRequest: return "INVALID_REQUEST";
    }
    return "UNKNOWN";
}

BuoyantPlaneBalanceResult solveBuoyantStabilityByPlane(
    const SpinDesignRequest& request,
    const Eigen::Vector3d& initialParameters,
    const BuoyantPlaneBalanceSettings& settings) {
    BuoyantPlaneBalanceResult result;
    try {
        validateSpinDesignRequest(request);
        validateSettings(settings);
        if (request.mode != SpinMode::BuoyantStability)
            throw std::invalid_argument("buoyant plane solver requires buoyant-stability mode");
        if (!initialParameters.allFinite())
            throw std::invalid_argument("initial buoyant plane parameters must be finite");
        result.constraintScales = fixedConstraintScales(request);
        result.targetTotalMass = request.buoyancy.liquidDensity *
            request.buoyancy.displacedLiquidVolume;
        if (!std::isfinite(result.targetTotalMass) || result.targetTotalMass <= 0.0)
            throw std::invalid_argument("buoyant total target mass is not representable");
    } catch (const std::exception& error) {
        result.status = BuoyantPlaneBalanceStatus::InvalidRequest;
        result.diagnostic = error.what();
        return result;
    }

    const double lengthScale = std::cbrt(result.constraintScales(0) /
                                          request.materialDensity);
    const double initialMassResidual = std::abs(
        request.fixedShellMoments.mass() - result.targetTotalMass);
    if (initialMassResidual > 1e-14 * std::max(1.0, result.targetTotalMass) &&
        request.fixedShellMoments.mass() >= result.targetTotalMass) {
        result.status = BuoyantPlaneBalanceStatus::InvalidRequest;
        result.diagnostic = "fixed shell alone exceeds the prescribed displaced-liquid mass";
        return result;
    }

    Eigen::Vector3d currentParameters = initialParameters;
    currentParameters.x() /= lengthScale;
    EvaluatedPlane current;
    try {
        Eigen::Vector3d actual = currentParameters;
        actual.x() *= lengthScale;
        current = evaluate(request, actual, result.constraintScales,
                           result.targetTotalMass, settings);
    } catch (const std::exception& error) {
        result.status = BuoyantPlaneBalanceStatus::GeometryEvaluationFailed;
        result.diagnostic = error.what();
        EvaluatedPlane failed;
        failed.coefficients << initialParameters.x(), initialParameters.y(),
            initialParameters.z(), -1.0;
        append(result, failed, 0, 0.0, false, 0);
        result.trace.back().diagnostic = error.what();
        return result;
    }
    retain(result, current, request);
    append(result, current, 0, 0.0, false, 0);
    if (current.residual <= settings.constraintTolerance) {
        verifyConvergedPlaneFamily(result, current, request, settings, lengthScale);
        return result;
    }

    int evaluationCount = 1;
    for (int iteration = 0; iteration < settings.maximumIterations; ++iteration) {
        Eigen::Matrix3d jacobian;
        try {
            jacobian = reducedJacobian(request, current.evaluation,
                                      result.constraintScales, lengthScale);
        } catch (const std::exception& error) {
            result.status = BuoyantPlaneBalanceStatus::DerivativeUnavailable;
            result.diagnostic = error.what();
            return result;
        }
        if (!jacobian.allFinite()) {
            result.status = BuoyantPlaneBalanceStatus::DerivativeUnavailable;
            result.diagnostic = "nonfinite reduced buoyant plane Jacobian";
            return result;
        }
        Eigen::FullPivLU<Eigen::Matrix3d> factor(jacobian);
        if (factor.rank() != 3 || factor.rcond() < 1e-13) {
            result.status = BuoyantPlaneBalanceStatus::SingularReducedSystem;
            result.diagnostic = "three-constraint buoyant plane Jacobian is singular or ill-conditioned";
            return result;
        }
        const Eigen::Vector3d direction = factor.solve(-current.normalizedConstraints);
        if (!direction.allFinite() || direction.norm() == 0.0) {
            result.status = BuoyantPlaneBalanceStatus::SingularReducedSystem;
            result.diagnostic = "reduced buoyant Newton direction is zero or nonfinite";
            return result;
        }

        bool accepted = false;
        double step = 1.0;
        for (int line = 0; line < settings.maximumLineSearchSteps &&
             step >= settings.minimumStep; ++line) {
            Eigen::Vector3d candidateScaled(current.parameters.x() / lengthScale,
                                            current.parameters.y(),
                                            current.parameters.z());
            candidateScaled += step * direction;
            Eigen::Vector3d candidateActual = candidateScaled;
            candidateActual.x() *= lengthScale;
            try {
                auto candidate = evaluate(request, candidateActual,
                    result.constraintScales, result.targetTotalMass, settings);
                ++evaluationCount;
                const bool reaches = candidate.residual <= settings.constraintTolerance;
                if (reaches && !candidate.evaluation.coefficientJacobian) {
                    append(result, candidate, iteration + 1, step, false, evaluationCount - 1);
                    retain(result, candidate, request);
                    result.status = BuoyantPlaneBalanceStatus::DerivativeUnavailable;
                    result.diagnostic = "equilibrium constraints reached at a plane-family differential singularity";
                    return result;
                }
                if (candidate.residual < current.residual &&
                    candidate.evaluation.coefficientJacobian.has_value()) {
                    current = std::move(candidate);
                    ++result.acceptedSteps;
                    result.iterations = iteration + 1;
                    append(result, current, iteration + 1, step, true, evaluationCount - 1);
                    retain(result, current, request);
                    accepted = true;
                    if (reaches) {
                        verifyConvergedPlaneFamily(result, current, request, settings, lengthScale);
                        return result;
                    }
                    break;
                }
                append(result, candidate, iteration + 1, step, false, evaluationCount - 1);
            } catch (const std::exception& error) {
                EvaluatedPlane failed;
                failed.coefficients << candidateActual.x(), candidateActual.y(),
                    candidateActual.z(), -1.0;
                ++evaluationCount;
                append(result, failed, iteration + 1, step, false, evaluationCount - 1);
                result.trace.back().diagnostic = error.what();
                result.diagnostic = error.what();
            }
            step *= settings.backtrackingFactor;
        }
        if (!accepted) {
            result.status = BuoyantPlaneBalanceStatus::LineSearchExhausted;
            if (result.diagnostic.empty())
                result.diagnostic = "backtracking could not reduce fixed-waterline constraint residual";
            return result;
        }
    }
    result.status = BuoyantPlaneBalanceStatus::IterationLimit;
    result.diagnostic = "buoyant plane iteration limit reached; no infeasibility claim";
    return result;
}

} // namespace spin
