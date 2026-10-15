#include "spin/QuadricHollowingSolver.h"
#include "QuadricMaterialExportGuard.h"
#include "spin/DualMonomialCoordinates.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/LU>
#include <Eigen/QR>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


constexpr int momentDimension = 10;


struct NewtonState
{
    std::optional<spin::SpinMomentSubspaceLinearization> subspace;
    spin::QuadricCoefficients10::Vector solverCoordinates;
    double rawResidual = std::numeric_limits<double>::infinity();
    spin::SpinMomentMapEvaluation map;
    spin::SpinMomentEvaluation problem;
    Eigen::VectorXd lambda;
    double mu = 0.0;
    Eigen::Matrix<double, momentDimension, 1> optimalityResidual =
        Eigen::Matrix<double, momentDimension, 1>::Zero();
    double normalizedResidual = std::numeric_limits<double>::infinity();
    double physicalResidual = std::numeric_limits<double>::infinity();
};


class InvalidNumericalState : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


void validateSettings(const spin::SpinSolverSettings& settings)
{
    if ((settings.equationProfile != spin::SpinEquationProfile::LegacyNativeNormalizedV1 &&
         settings.equationProfile != spin::SpinEquationProfile::DualMonomialRawV1) ||
        (settings.equationProfile == spin::SpinEquationProfile::DualMonomialRawV1 &&
         settings.maximumFunctionEvaluationsPerRestart > 60))
        throw std::invalid_argument("invalid equation profile or strict evaluation budget");
    if (
        settings.maximumIterations <= 0 || settings.maximumRestarts <= 0 ||
        settings.maximumFunctionEvaluationsPerRestart <= 0 ||
        !std::isfinite(settings.equationTolerance) ||
        settings.equationTolerance <= 0.0 ||
        !std::isfinite(settings.physicalRelativeTolerance) ||
        settings.physicalRelativeTolerance <= 0.0 ||
        settings.lineSearchReduction != 0.5 ||
        settings.maximumGeodesicStep != 0.5
    )
    {
        throw std::invalid_argument("invalid spherical Newton settings");
    }
}


NewtonState evaluateState(
    const spin::MassMomentProblem& problem,
    const spin::QuadricCoefficients10& quadric,
    spin::SpinMomentMapEvaluation map,
    spin::SpinEquationProfile profile
)
{
    NewtonState state;
    state.map = std::move(map);
    state.solverCoordinates = quadric.values();
    if (
        !state.map.rawMoments.allFinite() ||
        !state.map.momentJacobian.allFinite()
    )
    {
        throw InvalidNumericalState(
            "moment evaluator returned non-finite data"
        );
    }
    const double rawMomentNorm = state.map.rawMoments.norm();
    const double momentJacobianNorm = state.map.momentJacobian.norm();
    if (
        !std::isfinite(rawMomentNorm) ||
        !std::isfinite(momentJacobianNorm)
    )
    {
        throw InvalidNumericalState(
            "moment evaluator data overflowed its numerical scale"
        );
    }
    if (profile == spin::SpinEquationProfile::DualMonomialRawV1)
    {
        // Validate the native map before projecting: a radial projection must
        // not conceal invalid evaluator data. Check the transformed map too.
        const auto coordinates =
            spin::DualMonomialCoordinates10::fromGeometry(quadric);
        state.solverCoordinates = coordinates.values();
        state.map.momentJacobian =
            (state.map.momentJacobian * coordinates.geometryJacobian()).eval();
        if (!state.map.momentJacobian.allFinite() ||
            !std::isfinite(state.map.momentJacobian.norm()))
        {
            throw InvalidNumericalState("dual moment Jacobian overflowed");
        }
    }
    if (!state.map.derivativeDefined)
    {
        throw std::runtime_error("moment evaluator derivative is undefined" +
            (state.map.derivativeDiagnostic.empty() ? std::string() :
                ": " + state.map.derivativeDiagnostic));
    }
    state.problem = problem.evaluate(state.map.rawMoments);
    const Eigen::Index constraintCount = state.problem.constraints.size();
    if (
        state.problem.constraintJacobian.rows() != constraintCount ||
        state.problem.constraintJacobian.cols() != momentDimension ||
        state.problem.constraintHessians.size() !=
            static_cast<std::size_t>(constraintCount)
    )
    {
        throw std::runtime_error("moment problem derivative shape is invalid");
    }
    bool finiteConstraintHessians = true;
    for (const auto& hessian : state.problem.constraintHessians)
    {
        finiteConstraintHessians =
            finiteConstraintHessians && hessian.allFinite();
    }
    if (
        !std::isfinite(state.problem.objective) ||
        !state.problem.objectiveGradient.allFinite() ||
        !state.problem.objectiveHessian.allFinite() ||
        !state.problem.constraints.allFinite() ||
        !state.problem.constraintJacobian.allFinite() ||
        !finiteConstraintHessians
    )
    {
        throw InvalidNumericalState(
            "moment problem produced a non-finite objective or derivative"
        );
    }

    Eigen::MatrixXd projection(momentDimension, constraintCount + 1);
    if (constraintCount > 0)
    {
        projection.leftCols(constraintCount) =
            -state.problem.constraintJacobian.transpose();
    }
    projection.col(constraintCount) = state.solverCoordinates;
    const Eigen::VectorXd coefficients =
        projection.completeOrthogonalDecomposition().solve(
            -state.problem.objectiveGradient
        );
    if (!coefficients.allFinite())
    {
        throw InvalidNumericalState(
            "moment multiplier solve produced non-finite coefficients"
        );
    }
    state.lambda = coefficients.head(constraintCount);
    state.mu = coefficients(constraintCount);
    state.optimalityResidual =
        state.problem.objectiveGradient + state.mu * state.solverCoordinates;
    if (constraintCount > 0)
    {
        state.optimalityResidual -=
            state.problem.constraintJacobian.transpose() * state.lambda;
    }

    if (!state.optimalityResidual.allFinite())
    {
        throw InvalidNumericalState(
            "optimality residual produced non-finite coefficients"
        );
    }
    const double optimalitySquared =
        state.optimalityResidual.squaredNorm();
    const double constraintSquared = state.problem.constraints.squaredNorm();
    const double combinedSquared = optimalitySquared + constraintSquared;
    const double objectiveGradientNorm =
        state.problem.objectiveGradient.norm();
    const double constraintJacobianNorm =
        state.problem.constraintJacobian.norm();
    const double equationScale = std::max({
        1.0,
        objectiveGradientNorm,
        constraintJacobianNorm
    });
    const double physicalScale = std::max(1.0, rawMomentNorm);
    if (
        !std::isfinite(optimalitySquared) ||
        !std::isfinite(constraintSquared) ||
        !std::isfinite(combinedSquared) ||
        !std::isfinite(objectiveGradientNorm) ||
        !std::isfinite(constraintJacobianNorm) ||
        !std::isfinite(equationScale) ||
        !std::isfinite(physicalScale)
    )
    {
        throw InvalidNumericalState(
            "Eq. 14 residual normalization overflowed"
        );
    }
    state.rawResidual = std::sqrt(combinedSquared);
    state.normalizedResidual = state.rawResidual / equationScale;
    state.physicalResidual = std::sqrt(constraintSquared) / physicalScale;
    if (
        !std::isfinite(state.normalizedResidual) ||
        !std::isfinite(state.physicalResidual) || !std::isfinite(state.mu)
    )
    {
        throw InvalidNumericalState("Eq. 14 evaluation overflowed");
    }
    return state;
}


spin::SpinIterationTrace traceEntry(
    int restart,
    int iteration,
    int evaluations,
    double step,
    const NewtonState& state
)
{
    spin::SpinIterationTrace entry;
    entry.restart = restart;
    entry.iteration = iteration;
    entry.functionEvaluations = evaluations;
    entry.normalizedEquationResidual = state.normalizedResidual;
    entry.rawEquationResidual = state.rawResidual;
    entry.physicalRelativeResidual = state.physicalResidual;
    entry.objective = state.problem.objective;
    entry.multiplierMu = state.mu;
    entry.acceptedGeodesicStep = step;
    entry.quadricHash.clear();
    entry.materialGeometryHash = state.map.geometryHash;
    entry.materialTopologyHash = state.map.topologyHash;
    return entry;
}


std::string traceHash(const std::vector<spin::SpinIterationTrace>& trace,
                      spin::SpinEquationProfile profile)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    const bool strict=profile==spin::SpinEquationProfile::DualMonomialRawV1;
    stream << (strict ? "SpinNewtonTrace-v2;" : "SpinNewtonTrace-v1;") << std::hexfloat;
    if(strict) stream << spin::spinEquationProfileToken(profile) << ';'
                      << spin::spinSolverCoordinateIdentity(profile) << ';';
    for (const spin::SpinIterationTrace& entry : trace)
    {
        if(strict) {
            stream << *entry.rawEquationResidual << ';';
            for(double value:*entry.quadricCoefficients) stream << value << ';';
        }
        stream << entry.restart << ';' << entry.iteration << ';'
               << entry.functionEvaluations << ';'
               << entry.normalizedEquationResidual << ';'
               << entry.physicalRelativeResidual << ';'
               << entry.objective << ';' << entry.multiplierMu << ';'
               << entry.acceptedGeodesicStep << ';'
               << entry.quadricHash << ';'
               << entry.materialGeometryHash << ';'
               << entry.materialTopologyHash << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


std::array<double, momentDimension> coefficientArray(
    const spin::QuadricCoefficients10& quadric
)
{
    std::array<double, momentDimension> result{};
    for (int index = 0; index < momentDimension; ++index)
    {
        result[static_cast<std::size_t>(index)] = quadric.values()(index);
    }
    return result;
}


std::optional<double> finiteMetric(double value)
{
    return std::isfinite(value) ? std::optional<double>(value) : std::nullopt;
}


spin::SpinSolverEvent restartEvent(
    const char* type,
    int restart,
    int iteration,
    int evaluationIndex,
    int restartEvaluationIndex,
    const char* outcome,
    const char* reason
)
{
    spin::SpinSolverEvent event;
    event.type = type;
    event.restart = restart;
    event.iteration = iteration;
    event.evaluationIndex = evaluationIndex;
    event.restartEvaluationIndex = restartEvaluationIndex;
    event.outcome = outcome;
    event.reason = reason;
    return event;
}


std::string solverEventHash(const std::vector<spin::SpinSolverEvent>& events,
                            spin::SpinEquationProfile profile,
                            bool captureEvaluationPayload)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    const bool strict=profile==spin::SpinEquationProfile::DualMonomialRawV1;
    if (captureEvaluationPayload)
        stream << (strict ? "SpinSolverEvents-v6-capture-v1;" :
                             "SpinSolverEvents-v5-capture-v1;");
    else
        stream << (strict ? "SpinSolverEvents-v2;" : "SpinSolverEvents-v1;");
    stream << std::hexfloat;
    if(strict) stream << spin::spinEquationProfileToken(profile) << ';'
                      << spin::spinSolverCoordinateIdentity(profile) << ';';
    const auto optionalDouble = [&stream](const std::optional<double>& value) {
        if (value) stream << "some:" << *value;
        else stream << "none";
        stream << ';';
    };
    const auto optionalString = [&stream](
        const std::optional<std::string>& value
    ) {
        if (value) stream << "some:" << value->size() << ':' << *value;
        else stream << "none";
        stream << ';';
    };
    for (const spin::SpinSolverEvent& event : events)
    {
        if(strict) optionalDouble(event.rawEquationResidual);
        stream << event.type.size() << ':' << event.type << ';'
               << event.restart << ';' << event.iteration << ';'
               << event.evaluationIndex << ';'
               << event.restartEvaluationIndex << ';';
        optionalDouble(event.step);
        stream << event.outcome.size() << ':' << event.outcome << ';'
               << event.reason.size() << ':' << event.reason << ';'
               << event.detail.size() << ':' << event.detail << ';';
        if (event.quadricCoefficients)
        {
            stream << "some;";
            for (double coefficient : *event.quadricCoefficients)
            {
                stream << coefficient << ';';
            }
        }
        else
        {
            stream << "none;";
        }
        stream << event.quadricHash.size() << ':' << event.quadricHash << ';';
        optionalDouble(event.normalizedEquationResidual);
        optionalDouble(event.physicalRelativeResidual);
        optionalDouble(event.objective);
        optionalDouble(event.multiplierMu);
        optionalString(event.materialGeometryHash);
        optionalString(event.materialTopologyHash);
        if (captureEvaluationPayload)
        {
            if (!event.evaluationPayload)
            {
                stream << "payload:none;";
            }
            else
            {
                const auto& payload = *event.evaluationPayload;
                stream << "payload:spin-evaluation-payload-v1;";
                for (double coordinate : payload.solverCoordinates)
                    optionalDouble(coordinate);
                for (double moment : payload.solverRawMoments)
                    optionalDouble(moment);
                stream << (payload.variableRawMoments ? "variable:some;" :
                                                          "variable:none;");
                if (payload.variableRawMoments)
                    for (double moment : *payload.variableRawMoments)
                        optionalDouble(moment);
                stream << (payload.solverMomentJacobian ? "jacobian:some;" :
                                                           "jacobian:none;");
                if (payload.solverMomentJacobian)
                {
                    for (Eigen::Index row = 0; row < momentDimension; ++row)
                        for (Eigen::Index column = 0; column < momentDimension; ++column)
                            optionalDouble((*payload.solverMomentJacobian)(row, column));
                }
                stream << (payload.derivativeDefined ? "defined;" : "undefined;");
                optionalString(payload.jacobianCoordinateIdentity);
                optionalString(payload.evaluationMapIdentity);
                optionalString(payload.booleanBackendIdentity);
                optionalString(payload.derivativeIdentity);
                optionalString(payload.derivativeDiagnostic);
                optionalString(payload.geometryHash);
                optionalString(payload.topologyHash);
                optionalDouble(payload.evaluatorWallSeconds);
                optionalDouble(payload.residualEvaluationSeconds);
                optionalDouble(payload.derivativeEvaluationSeconds);
            }
            if (!event.newtonSystemPayload)
            {
                stream << "newton-system:none;";
            }
            else
            {
                const auto& payload = *event.newtonSystemPayload;
                stream << "newton-system:spin-newton-system-payload-v2;"
                       << payload.matrix.rows() << ';' << payload.matrix.cols() << ';'
                       << payload.computedRank << ';'
                       << payload.equilibratedRank << ';';
                for (Eigen::Index row = 0; row < payload.matrix.rows(); ++row)
                    for (Eigen::Index column = 0; column < payload.matrix.cols(); ++column)
                        optionalDouble(payload.matrix(row, column));
                stream << payload.rightHandSide.size() << ';';
                for (Eigen::Index i = 0; i < payload.rightHandSide.size(); ++i)
                    optionalDouble(payload.rightHandSide(i));
                optionalDouble(payload.multiplierMu);
                stream << payload.constraintMultipliers.size() << ';';
                for (Eigen::Index i = 0; i < payload.constraintMultipliers.size(); ++i)
                    optionalDouble(payload.constraintMultipliers(i));
                stream << payload.rowScales.size() << ';';
                for (Eigen::Index i = 0; i < payload.rowScales.size(); ++i)
                    optionalDouble(payload.rowScales(i));
                stream << payload.columnScales.size() << ';';
                for (Eigen::Index i = 0; i < payload.columnScales.size(); ++i)
                    optionalDouble(payload.columnScales(i));
                optionalDouble(payload.linearResidualInfinityNorm);
                optionalDouble(payload.linearBackwardError);
            }
        }
    }
    return manufacturing::computeContentSha256(stream.str());
}


double stoppingResidual(const NewtonState& state, const spin::SpinSolverSettings& settings)
{
    return settings.equationProfile==spin::SpinEquationProfile::DualMonomialRawV1
        ? state.rawResidual : state.normalizedResidual;
}

bool converged(
    const NewtonState& state,
    const spin::SpinSolverSettings& settings
)
{
    return stoppingResidual(state, settings) <= settings.equationTolerance &&
        state.physicalResidual <= settings.physicalRelativeTolerance &&
        state.mu > 0.0;
}


std::uint64_t splitMix64(std::uint64_t* state)
{
    *state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t value = *state;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}


double deterministicSignedUnit(std::uint64_t* state)
{
    const std::uint64_t bits = splitMix64(state) >> 11U;
    const double unit = static_cast<double>(bits) / 9007199254740992.0;
    return 2.0 * unit - 1.0;
}


spin::SpinMomentProblem problemForRequest(
    const spin::SpinDesignRequest& request
)
{
    switch (request.mode)
    {
    case spin::SpinMode::StaticBalance:
        return request.staticObjective == spin::StaticObjective::ParetoTradeoff
            ? spin::SpinMomentProblem::staticBalancePareto(
                request.staticParetoWeight, spin::staticObjectiveMassScale(request)
            )
            : spin::SpinMomentProblem::staticBalance(
                request.staticObjective, spin::staticObjectiveMassScale(request)
            );
    case spin::SpinMode::BuoyantStability:
        return spin::SpinMomentProblem::buoyantStability(
            request.buoyancy.liquidDensity *
                request.buoyancy.displacedLiquidVolume,
            request.buoyancy.centerOfBuoyancyTarget.head<2>()
        );
    case spin::SpinMode::FidgetSpinner:
        return spin::SpinMomentProblem::fidgetSpinner();
    case spin::SpinMode::YoYo:
        return spin::SpinMomentProblem::yoYo();
    case spin::SpinMode::SpinningTop:
        return spin::SpinMomentProblem::spinningTop(
            request.centerOfMassWeight, request.spinningTopObjective
        );
    }
    throw std::invalid_argument("unsupported Spin-It Faster mode");
}


std::string hollowingOutputHash(
    const spin::QuadricHollowingResult& result
)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    const bool strict=result.equationProfile==spin::SpinEquationProfile::DualMonomialRawV1;
    if(strict) {
        stream << "QuadricHollowingResult-v3;" << spin::spinEquationProfileToken(result.equationProfile)
               << ';' << spin::spinSolverCoordinateIdentity(result.equationProfile) << ';' << std::hexfloat;
        if(result.rawEquationResidual) stream << *result.rawEquationResidual;
        else stream << "none";
        stream << ';';
    } else stream << "QuadricHollowingResult-v2;";
    stream << result.inputHash << ';'
           << result.geometryHash << ';' << result.traceHash << ';'
           << result.eventHash << ';'
           << spin::spinSolverStatusToken(result.status) << ';'
           << std::hexfloat << result.normalizedEquationResidual << ';'
           << result.physicalRelativeResidual << ';';
    if (std::isfinite(result.buoyantMetacenterHeight))
    {
        stream << result.buoyantMetacenterHeight << ';'
               << result.buoyantStabilityMargin << ';';
    }
    else
    {
        stream << "not-buoyant;";
    }
    for (int index = 0; index < result.quadricCoefficients.size(); ++index)
    {
        stream << result.quadricCoefficients(index) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


} // namespace


namespace spin
{

SpinMomentSubspaceLinearization linearizeSpinMomentSubspace(
    const SpinMomentProblem& problem,
    const QuadricCoefficients10& quadric,
    const Eigen::Matrix<double, 10, 1>& rawMoments,
    const Eigen::MatrixXd& supportBasis,
    const Eigen::MatrixXd& nativeTangentBasis,
    const Eigen::MatrixXd& nativeDirectionalJacobian
)
{
    using Matrix10 = Eigen::Matrix<double, 10, 10>;
    const Eigen::Index k = supportBasis.cols(), d = k - 1;
    const auto require = [](bool condition, const char* message) {
        if (!condition) throw std::invalid_argument(message);
    };
    require(k >= 2 && k <= 10 && supportBasis.rows() == 10 &&
        nativeTangentBasis.rows() == 10 && nativeTangentBasis.cols() == d &&
        nativeDirectionalJacobian.rows() == 10 && nativeDirectionalJacobian.cols() == d,
        "supported Eq.14 basis/derivative dimensions are invalid");
    require(rawMoments.allFinite() && supportBasis.allFinite() &&
        nativeTangentBasis.allFinite() && nativeDirectionalJacobian.allFinite() &&
        std::isfinite(nativeDirectionalJacobian.norm()),
        "supported Eq.14 inputs are non-finite");
    const auto coordinates = DualMonomialCoordinates10::fromGeometry(quadric);
    const auto s = coordinates.values();
    const Matrix10 complement = Matrix10::Identity() - supportBasis * supportBasis.transpose();
    // These are coordinate-contract checks, not a relaxed scientific stop gate.
    constexpr double contractTolerance = 1e-12;
    require((supportBasis.transpose()*supportBasis-Eigen::MatrixXd::Identity(k,k)).norm() <= contractTolerance &&
        (nativeTangentBasis.transpose()*nativeTangentBasis-Eigen::MatrixXd::Identity(d,d)).norm() <= contractTolerance &&
        (complement*s).norm() <= contractTolerance &&
        (nativeTangentBasis.transpose()*quadric.values()).norm() <= contractTolerance,
        "supported Eq.14 basis is not orthonormal/tangent or excludes current coordinates");

    Eigen::VectorXd householder = supportBasis.transpose()*s;
    householder(0) += householder(0) >= 0.0 ? 1.0 : -1.0;
    householder.normalize();
    SpinMomentSubspaceLinearization result;
    result.solverTangent = supportBasis *
        (Eigen::MatrixXd::Identity(k,k)-2.0*householder*householder.transpose()).rightCols(d);
    const Eigen::MatrixXd nativeMotion = coordinates.geometryJacobian()*result.solverTangent;
    const Eigen::MatrixXd change = nativeTangentBasis.transpose()*nativeMotion;
    require((nativeTangentBasis*change-nativeMotion).norm() <= contractTolerance,
        "material derivative does not span supported dual tangent");
    result.momentTangentJacobian = nativeDirectionalJacobian*change;
    const auto value = problem.evaluate(rawMoments);
    const Eigen::Index n = value.constraints.size();
    require(value.constraintJacobian.rows()==n && value.constraintJacobian.cols()==10 &&
        value.constraintHessians.size()==static_cast<std::size_t>(n) &&
        std::isfinite(value.objective) && value.objectiveGradient.allFinite() &&
        value.objectiveHessian.allFinite() && value.constraints.allFinite() &&
        value.constraintJacobian.allFinite(), "supported moment problem is invalid");
    require((complement*value.objectiveGradient).norm() <= contractTolerance &&
        (complement*value.constraintJacobian.transpose()).norm() <= contractTolerance,
        "support excludes objective or constraint gradient equations");
    // Check individual derivatives before combining them: cancellation for one
    // multiplier is not evidence that omitted equations are identically zero.
    require((complement*value.objectiveHessian*result.momentTangentJacobian).norm() <= contractTolerance,
        "support excludes objective gradient differential");
    for (const auto& h : value.constraintHessians)
        require(h.allFinite() &&
            (complement*h*result.momentTangentJacobian).norm() <= contractTolerance,
            "support excludes constraint gradient differential");

    Eigen::MatrixXd projection(10,n+1);
    projection.leftCols(n) = -value.constraintJacobian.transpose();
    projection.col(n) = s;
    const Eigen::VectorXd multipliers = projection.completeOrthogonalDecomposition().solve(-value.objectiveGradient);
    result.multipliers = multipliers.head(n);
    result.multiplierMu = multipliers(n);
    result.fullResidual.resize(10+n);
    result.fullResidual.head(10) = value.objectiveGradient + projection*multipliers;
    result.fullResidual.tail(n) = value.constraints;
    result.rawResidual = result.fullResidual.norm();
    Matrix10 hessian = value.objectiveHessian;
    for (Eigen::Index i=0;i<n;++i)hessian-=multipliers(i)*value.constraintHessians[i];
    result.matrix = Eigen::MatrixXd::Zero(k+n,k+n);
    result.matrix.topLeftCorner(k,d) = supportBasis.transpose()*
        (hessian*result.momentTangentJacobian+result.multiplierMu*result.solverTangent);
    result.matrix.block(0,d,k,n) = -supportBasis.transpose()*value.constraintJacobian.transpose();
    result.matrix.block(0,d+n,k,1) = supportBasis.transpose()*s;
    result.matrix.block(k,0,n,d) = value.constraintJacobian*result.momentTangentJacobian;
    result.rightHandSide.resize(k+n);
    result.rightHandSide.head(k) = -supportBasis.transpose()*result.fullResidual.head(10);
    result.rightHandSide.tail(n) = -value.constraints;
    require(result.matrix.allFinite() && result.rightHandSide.allFinite() &&
        multipliers.allFinite() && std::isfinite(result.rawResidual),
        "supported Eq.14 linearization overflowed");
    return result;
}


SpinNewtonLinearSolveResult solveSpinNewtonLinearSystem(
    const Eigen::MatrixXd& matrix,
    const Eigen::VectorXd& rightHandSide
)
{
    SpinNewtonLinearSolveResult result;
    if (matrix.rows() == 0 || matrix.rows() != matrix.cols() ||
        rightHandSide.size() != matrix.rows() || !matrix.allFinite() ||
        !rightHandSide.allFinite())
    {
        result.diagnostic = "Newton linear system must be nonempty, square, and finite";
        return result;
    }

    const Eigen::Index systemSize = matrix.rows();
    const double backwardErrorLimit = 128.0 *
        std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max<Eigen::Index>(systemSize, 1));
    const auto assessOriginalSystem = [&](const Eigen::VectorXd& solution) {
        const Eigen::VectorXd residual = matrix * solution - rightHandSide;
        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        const double matrixInfinityNorm = matrix.cwiseAbs()
            .rowwise().sum().maxCoeff();
        const double denominator = matrixInfinityNorm *
            solution.lpNorm<Eigen::Infinity>() +
            rightHandSide.lpNorm<Eigen::Infinity>();
        result.backwardError = denominator > 0.0
            ? result.residualInfinityNorm / denominator
            : std::numeric_limits<double>::infinity();
        return residual.allFinite() && std::isfinite(result.backwardError) &&
            result.backwardError <= backwardErrorLimit;
    };

    Eigen::FullPivLU<Eigen::MatrixXd> originalFactorization(matrix);
    result.originalRank = originalFactorization.rank();
    result.equilibratedRank = result.originalRank;
    if (result.originalRank == systemSize)
    {
        result.solution = originalFactorization.solve(rightHandSide);
        if (result.solution.allFinite() && assessOriginalSystem(result.solution))
        {
            result.solved = true;
            result.diagnostic = "unscaled solve passed original-system backward-error check";
            return result;
        }
    }

    // Max-norm diagonal scaling is exactly reversible: A x = b becomes
    // (D_r^-1 A D_c^-1) y = D_r^-1 b with x = D_c^-1 y.
    // No regularizer, pseudoinverse, or altered equation is introduced.
    result.usedEquilibration = true;
    result.rowScales.resize(systemSize);
    Eigen::MatrixXd rowScaled = matrix;
    Eigen::VectorXd scaledRightHandSide = rightHandSide;
    for (Eigen::Index row = 0; row < systemSize; ++row)
    {
        const double scale = matrix.row(row).cwiseAbs().maxCoeff();
        if (!std::isfinite(scale) || scale == 0.0)
        {
            result.diagnostic = "row equilibration found a zero or nonfinite row";
            return result;
        }
        result.rowScales(row) = scale;
        rowScaled.row(row) /= scale;
        scaledRightHandSide(row) /= scale;
    }
    result.columnScales.resize(systemSize);
    Eigen::MatrixXd equilibrated = rowScaled;
    for (Eigen::Index column = 0; column < systemSize; ++column)
    {
        const double scale = rowScaled.col(column).cwiseAbs().maxCoeff();
        if (!std::isfinite(scale) || scale == 0.0)
        {
            result.diagnostic = "column equilibration found a zero or nonfinite column";
            return result;
        }
        result.columnScales(column) = scale;
        equilibrated.col(column) /= scale;
    }
    if (!equilibrated.allFinite() || !scaledRightHandSide.allFinite())
    {
        result.diagnostic = "reversibly equilibrated Newton system is nonfinite";
        return result;
    }

    Eigen::FullPivLU<Eigen::MatrixXd> equilibratedFactorization(equilibrated);
    result.equilibratedRank = equilibratedFactorization.rank();
    if (result.equilibratedRank < systemSize)
    {
        result.diagnostic = "Newton system remains rank deficient after reversible equilibration";
        return result;
    }
    const Eigen::VectorXd equilibratedSolution =
        equilibratedFactorization.solve(scaledRightHandSide);
    result.solution = equilibratedSolution.cwiseQuotient(result.columnScales);
    if (!result.solution.allFinite() || !assessOriginalSystem(result.solution))
    {
        result.diagnostic = "equilibrated solve failed the original-system backward-error check";
        return result;
    }
    result.solved = true;
    result.diagnostic = "reversible diagonal equilibration recovered full rank and passed the original-system backward-error check";
    return result;
}


static SpinNewtonResult solveSpinMomentSystemImpl(
    const MassMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink,
    const std::optional<SpinMomentSubspaceLinearization>* subspaceState = nullptr,
    const std::string& subspaceIdentity = {}
)
{
    validateSettings(settings);
    if (!evaluator)
    {
        throw std::invalid_argument("spherical Newton evaluator is empty");
    }
    if (initialGuesses.empty())
    {
        throw std::invalid_argument(
            "spherical Newton requires deterministic initial guesses"
        );
    }

    SpinNewtonResult result;
    result.equationProfile=settings.equationProfile;
    const auto emitEvent = [&eventSink, &result](const SpinSolverEvent& event) {
        if (eventSink)
        {
            SpinSolverEvent snapshot = event;
            snapshot.equationProfile = result.equationProfile;
            eventSink(snapshot);
        }
    };
    const bool strict=settings.equationProfile==SpinEquationProfile::DualMonomialRawV1;
    result.eventStreamVersion = settings.captureEvaluationPayload
        ? (strict ? 6 : 5)
        : (strict ? 2 : 1);
    bool sawGeometryFailure = false;
    const auto finish = [&result, &settings]() {
        for(auto& event:result.events) event.equationProfile=result.equationProfile;
        for(auto& entry:result.trace) entry.equationProfile=result.equationProfile;
        result.traceHash = ::traceHash(result.trace, result.equationProfile);
        result.eventHash = ::solverEventHash(
            result.events, result.equationProfile,
            settings.captureEvaluationPayload
        );
    };
    const int restartLimit = std::min(
        settings.maximumRestarts,
        static_cast<int>(initialGuesses.size())
    );
    for (int restart = 0; restart < restartLimit; ++restart)
    {
        result.restartCount = restart + 1;
        result.diagnostic.clear();
        QuadricCoefficients10 current = initialGuesses[restart];
        int evaluationsThisRestart = 0;
        result.restartFunctionEvaluations.push_back(0);
        SpinSolverEvent start = restartEvent(
            "attempt_start", restart, 0, result.functionEvaluations, 0,
            "started", "restart_started"
        );
        start.quadricCoefficients = coefficientArray(current);
        start.quadricHash = current.hash();
        if (subspaceState) start.detail = "supported-space-v1;support=" + subspaceIdentity;
        result.events.push_back(std::move(start));
        emitEvent(result.events.back());

        struct EvaluationResult
        {
            std::optional<NewtonState> state;
            std::string error;
            bool invalidNumericalState = false;
        };
        const auto evaluateAttempt = [&](const QuadricCoefficients10& quadric,
                                         int iteration,
                                         std::optional<double> step) {
            ++evaluationsThisRestart;
            ++result.functionEvaluations;
            result.restartFunctionEvaluations.back() =
                evaluationsThisRestart;
            SpinSolverEvent event = restartEvent(
                "evaluation", restart, iteration,
                result.functionEvaluations, evaluationsThisRestart,
                "failed", "evaluation_exception"
            );
            event.step = step;
            event.quadricCoefficients = coefficientArray(quadric);
            event.quadricHash = quadric.hash();
            SpinEvaluationTiming timing;
            timing.evaluationIndex = result.functionEvaluations;
            timing.restart = restart;
            timing.restartEvaluationIndex = evaluationsThisRestart;
            const auto started = std::chrono::steady_clock::now();
            const auto recordFailure = [&](const std::exception& error,
                                           const char* reason,
                                           bool invalidNumericalState) {
                const auto failed = std::chrono::steady_clock::now();
                if (!timing.evaluatorReturned)
                {
                    timing.wallSeconds = std::chrono::duration<double>(
                        failed - started
                    ).count();
                    result.evaluatorWallSeconds += timing.wallSeconds;
                }
                event.reason = reason;
                event.detail = error.what();
                result.evaluationTimings.push_back(std::move(timing));
                result.events.push_back(std::move(event));
                emitEvent(result.events.back());
                return EvaluationResult{
                    std::nullopt,
                    error.what(),
                    invalidNumericalState
                };
            };
            try
            {
                SpinMomentMapEvaluation map = evaluator(quadric);
                const auto returned = std::chrono::steady_clock::now();
                timing.evaluatorReturned = true;
                timing.wallSeconds =
                    std::chrono::duration<double>(returned - started).count();
                timing.reportedResidualSeconds =
                    map.residualEvaluationSeconds;
                timing.reportedDerivativeSeconds =
                    map.derivativeEvaluationSeconds;
                result.evaluatorWallSeconds += timing.wallSeconds;
                result.residualEvaluationSeconds +=
                    map.residualEvaluationSeconds;
                result.derivativeEvaluationSeconds +=
                    map.derivativeEvaluationSeconds;
                if (!map.geometryHash.empty())
                {
                    event.materialGeometryHash = map.geometryHash;
                }
                if (!map.topologyHash.empty())
                {
                    event.materialTopologyHash = map.topologyHash;
                }
                if (settings.captureEvaluationPayload &&
                    map.rawMoments.allFinite())
                {
                    SpinSolverEvent::EvaluationPayload payload;
                    Eigen::Matrix<double, momentDimension, 1> solverCoordinates =
                        quadric.values();
                    Eigen::Matrix<double, momentDimension, momentDimension>
                        solverJacobian = map.momentJacobian;
                    if (settings.equationProfile ==
                        SpinEquationProfile::DualMonomialRawV1)
                    {
                        const auto coordinates =
                            DualMonomialCoordinates10::fromGeometry(quadric);
                        solverCoordinates = coordinates.values();
                        if (map.derivativeDefined)
                        {
                            solverJacobian = (solverJacobian *
                                coordinates.geometryJacobian()).eval();
                        }
                    }
                    bool payloadFinite = solverCoordinates.allFinite();
                    if (map.derivativeDefined)
                        payloadFinite = payloadFinite && solverJacobian.allFinite();
                    if (payloadFinite)
                    {
                        for (Eigen::Index i = 0; i < momentDimension; ++i)
                        {
                            payload.solverCoordinates[static_cast<std::size_t>(i)] =
                                solverCoordinates(i);
                            payload.solverRawMoments[static_cast<std::size_t>(i)] =
                                map.rawMoments(i);
                        }
                        if (map.variableRawMoments &&
                            map.variableRawMoments->allFinite())
                        {
                            std::array<double, momentDimension> variable{};
                            for (Eigen::Index i = 0; i < momentDimension; ++i)
                                variable[static_cast<std::size_t>(i)] =
                                    (*map.variableRawMoments)(i);
                            payload.variableRawMoments = variable;
                        }
                        if (map.derivativeDefined && !subspaceState)
                            payload.solverMomentJacobian = solverJacobian;
                        payload.derivativeDefined = map.derivativeDefined;
                        payload.jacobianCoordinateIdentity =
                            spin::spinSolverCoordinateIdentity(
                                settings.equationProfile);
                        if (subspaceState)
                            payload.jacobianCoordinateIdentity =
                                "NATIVE_DIRECTIONAL_IN_SUBSPACE_ATTACHMENT_V1";
                        payload.evaluationMapIdentity = map.evaluationMapIdentity;
                        payload.booleanBackendIdentity = map.booleanBackendIdentity;
                        payload.derivativeIdentity = map.derivativeIdentity;
                        payload.derivativeDiagnostic = map.derivativeDiagnostic;
                        payload.geometryHash = map.geometryHash;
                        payload.topologyHash = map.topologyHash;
                        payload.evaluatorWallSeconds = timing.wallSeconds;
                        payload.residualEvaluationSeconds =
                            map.residualEvaluationSeconds;
                        payload.derivativeEvaluationSeconds =
                            map.derivativeEvaluationSeconds;
                        event.evaluationPayload = std::move(payload);
                    }
                }
                NewtonState state = evaluateState(
                    problem, quadric, std::move(map), settings.equationProfile
                );
                if (subspaceState)
                {
                    if (!*subspaceState)
                        throw InvalidNumericalState("supported evaluation has no current linearization");
                    state.subspace = **subspaceState;
                }
                event.outcome = "evaluated";
                event.rawEquationResidual = finiteMetric(state.rawResidual);
                event.reason = "evaluation_succeeded";
                event.detail = "evaluation_succeeded";
                if (!state.map.evaluationMapIdentity.empty())
                    event.detail += ";evaluation_map=" + state.map.evaluationMapIdentity;
                if (!state.map.booleanBackendIdentity.empty())
                    event.detail += ";backend=" + state.map.booleanBackendIdentity;
                if (!state.map.derivativeIdentity.empty())
                    event.detail += ";derivative=" + state.map.derivativeIdentity;
                event.normalizedEquationResidual = finiteMetric(
                    state.normalizedResidual
                );
                event.physicalRelativeResidual = finiteMetric(
                    state.physicalResidual
                );
                event.objective = finiteMetric(state.problem.objective);
                event.multiplierMu = finiteMetric(state.mu);
                result.evaluationTimings.push_back(std::move(timing));
                result.events.push_back(std::move(event));
                return EvaluationResult{std::move(state), "", false};
            }
            catch (const InvalidNumericalState& error)
            {
                return recordFailure(
                    error, "invalid_numerical_state", true
                );
            }
            catch (const std::exception& error)
            {
                return recordFailure(error, "evaluation_exception", false);
            }
            catch (...)
            {
                const auto failed = std::chrono::steady_clock::now();
                timing.wallSeconds = std::chrono::duration<double>(
                    failed - started
                ).count();
                result.evaluatorWallSeconds += timing.wallSeconds;
                result.evaluationTimings.push_back(std::move(timing));
                throw;
            }
        };
        const auto terminate = [&](int iteration, const char* reason) {
            result.events.push_back(restartEvent(
                "attempt_terminal", restart, iteration,
                result.functionEvaluations, evaluationsThisRestart,
                "terminated", reason
            ));
            emitEvent(result.events.back());
        };

        EvaluationResult initial = evaluateAttempt(current, 0, 0.0);
        if (!initial.state)
        {
            sawGeometryFailure = true;
            result.diagnostic = initial.error;
            terminate(
                0,
                initial.invalidNumericalState
                    ? "initial_numerical_state_invalid"
                    : "initial_evaluation_failed"
            );
            continue;
        }
        NewtonState state = std::move(*initial.state);
        if (strict && state.mu <= 0.0)
        {
            result.events.back().outcome = "rejected";
            result.events.back().reason = "initial_multiplier_nonpositive";
            result.events.back().detail +=
                ";strict Eq. 14 initial guess requires mu > 0";
            emitEvent(result.events.back());
            result.diagnostic =
                "strict Eq. 14 initial guess has nonpositive multiplier";
            terminate(0, "initial_multiplier_nonpositive");
            continue;
        }
        result.events.back().outcome = "accepted";
        result.events.back().reason = "initial_evaluation_succeeded";
        emitEvent(result.events.back());
        result.trace.push_back(
            traceEntry(restart, 0, evaluationsThisRestart, 0.0, state)
        );
        result.trace.back().quadricHash = current.hash();
        if(strict) result.trace.back().quadricCoefficients=coefficientArray(current);
        if ((state.subspace ? state.subspace->momentTangentJacobian.norm()
                            : state.map.momentJacobian.norm()) <= 1e-15)
        {
            result.status = SpinSolverStatus::OptimizationDidNotConverge;
            result.diagnostic = "initial moment Jacobian vanishes";
            terminate(0, "initial_sensitivity_vanished");
            continue;
        }
        if (converged(state, settings))
        {
            result.status = SpinSolverStatus::Converged;
            result.quadric = current;
            result.normalizedEquationResidual = state.normalizedResidual;
            result.rawEquationResidual = state.rawResidual;
            result.physicalRelativeResidual = state.physicalResidual;
            result.multiplierMu = state.mu;
            result.multipliers = state.lambda;
            result.diagnostic.clear();
            terminate(0, "converged");
            finish();
            return result;
        }

        for (int iteration = 1; iteration <= settings.maximumIterations;
             ++iteration)
        {
            const Eigen::Index constraintCount =
                state.problem.constraints.size();
            Eigen::Matrix<double, momentDimension, momentDimension> hessianL =
                state.problem.objectiveHessian;
            for (Eigen::Index constraint = 0;
                 constraint < constraintCount;
                 ++constraint)
            {
                hessianL -= state.lambda(constraint) *
                    state.problem.constraintHessians[
                        static_cast<std::size_t>(constraint)
                    ];
            }
            if (!hessianL.allFinite())
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic =
                    "Eq. 14 Lagrangian Hessian overflowed";
                terminate(iteration, "newton_system_invalid");
                break;
            }

            const Eigen::Index systemSize = state.subspace
                ? state.subspace->matrix.rows() : momentDimension + constraintCount + 1;
            Eigen::MatrixXd system;
            Eigen::VectorXd rightHandSide;
            if (state.subspace)
            {
                system = state.subspace->matrix;
                rightHandSide = state.subspace->rightHandSide;
            }
            else
            {
            system = Eigen::MatrixXd::Zero(systemSize, systemSize);
            rightHandSide = Eigen::VectorXd::Zero(systemSize);
            system.topLeftCorner(momentDimension, momentDimension) =
                hessianL * state.map.momentJacobian +
                state.mu *
                    Eigen::Matrix<double, momentDimension, momentDimension>::Identity();
            if (constraintCount > 0)
            {
                system.block(
                    0,
                    momentDimension,
                    momentDimension,
                    constraintCount
                ) = -state.problem.constraintJacobian.transpose();
                system.block(
                    momentDimension,
                    0,
                    constraintCount,
                    momentDimension
                ) = state.problem.constraintJacobian *
                    state.map.momentJacobian;
            }
            system.block(
                0,
                momentDimension + constraintCount,
                momentDimension,
                1
            ) = state.solverCoordinates;
            system.block(
                momentDimension + constraintCount,
                0,
                1,
                momentDimension
            ) = 2.0 * state.solverCoordinates.transpose();
            rightHandSide.head(momentDimension) =
                -state.optimalityResidual;
            if (constraintCount > 0)
            {
                rightHandSide.segment(momentDimension, constraintCount) =
                    -state.problem.constraints;
            }
            }
            if (!system.allFinite() || !rightHandSide.allFinite())
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Eq. 14 bordered Newton system overflowed";
                terminate(iteration, "newton_system_invalid");
                break;
            }

            const SpinNewtonLinearSolveResult linearSolve =
                solveSpinNewtonLinearSystem(system, rightHandSide);
            const auto captureNewtonSystem = [&](SpinSolverEvent& event) {
                if (!settings.captureEvaluationPayload)
                    return;
                SpinSolverEvent::NewtonSystemPayload payload;
                payload.matrix = system;
                payload.rightHandSide = rightHandSide;
                payload.computedRank = linearSolve.originalRank;
                payload.equilibratedRank = linearSolve.equilibratedRank;
                payload.multiplierMu = state.mu;
                payload.constraintMultipliers = state.lambda;
                payload.rowScales = linearSolve.rowScales;
                payload.columnScales = linearSolve.columnScales;
                payload.linearResidualInfinityNorm =
                    linearSolve.residualInfinityNorm;
                payload.linearBackwardError = linearSolve.backwardError;
                event.newtonSystemPayload = std::move(payload);
            };
            if (!linearSolve.solved)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = linearSolve.diagnostic;
                const char* terminalReason =
                    linearSolve.equilibratedRank < systemSize
                    ? "bordered_system_rank_deficient"
                    : "newton_linear_solve_backward_error_failed";
                SpinSolverEvent event = restartEvent(
                    "attempt_terminal", restart, iteration,
                    result.functionEvaluations, evaluationsThisRestart,
                    "terminated", terminalReason
                );
                event.detail = linearSolve.diagnostic;
                if (strict)
                {
                    event.quadricCoefficients = coefficientArray(current);
                    event.quadricHash = current.hash();
                    event.multiplierMu = state.mu;
                    event.normalizedEquationResidual = state.normalizedResidual;
                    event.physicalRelativeResidual = state.physicalResidual;
                    event.rawEquationResidual = state.rawResidual;
                }
                captureNewtonSystem(event);
                result.events.push_back(std::move(event));
                emitEvent(result.events.back());
                break;
            }
            if (linearSolve.usedEquilibration || settings.captureEvaluationPayload)
            {
                std::ostringstream detail;
                detail.imbue(std::locale::classic());
                detail << std::hexfloat
                       << "original_rank=" << linearSolve.originalRank
                       << ";equilibrated_rank=" << linearSolve.equilibratedRank
                       << ";residual_inf=" << linearSolve.residualInfinityNorm
                       << ";backward_error=" << linearSolve.backwardError;
                SpinSolverEvent event = restartEvent(
                    "newton_linear_solve", restart, iteration,
                    result.functionEvaluations, evaluationsThisRestart,
                    "accepted", linearSolve.usedEquilibration
                        ? "rank_recovered_by_reversible_equilibration"
                        : "regular_system_payload_captured"
                );
                event.detail = detail.str();
                if (strict)
                {
                    event.quadricCoefficients = coefficientArray(current);
                    event.quadricHash = current.hash();
                    event.multiplierMu = state.mu;
                    event.normalizedEquationResidual = state.normalizedResidual;
                    event.physicalRelativeResidual = state.physicalResidual;
                    event.rawEquationResidual = state.rawResidual;
                }
                captureNewtonSystem(event);
                result.events.push_back(std::move(event));
                emitEvent(result.events.back());
            }
            Eigen::Matrix<double, momentDimension, 1> direction;
            if (state.subspace)
                direction = state.subspace->solverTangent *
                    linearSolve.solution.head(state.subspace->solverTangent.cols());
            else direction = linearSolve.solution.head(momentDimension);
            direction = strict
                ? DualMonomialCoordinates10::fromGeometry(current).projectToTangent(direction)
                : current.projectToTangent(direction);
            const double directionNorm = direction.norm();
            if (!direction.allFinite() || !std::isfinite(directionNorm))
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Newton direction is non-finite";
                terminate(iteration, "newton_direction_invalid");
                break;
            }
            if (directionNorm <= 1e-15)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Newton direction vanished before convergence";
                terminate(iteration, "newton_direction_vanished");
                break;
            }

            double step = std::min(
                directionNorm,
                settings.maximumGeodesicStep
            );
            bool accepted = false;
            bool stationaryNonpositiveMu = false;
            bool topologyTransitionRejected = false;
            NewtonState candidateState;
            QuadricCoefficients10 candidate = current;
            while (
                evaluationsThisRestart <
                    settings.maximumFunctionEvaluationsPerRestart &&
                step > std::numeric_limits<double>::epsilon()
            )
            {
                candidate = strict
                    ? DualMonomialCoordinates10::fromGeometry(current).geodesicUpdate(direction, step).geometry()
                    : current.geodesicUpdate(direction, step);
                EvaluationResult trial = evaluateAttempt(
                    candidate, iteration, step
                );
                if (trial.state)
                {
                    candidateState = std::move(*trial.state);
                    if (strict && candidateState.mu <= 0.0 &&
                        stoppingResidual(candidateState, settings) <=
                            settings.equationTolerance &&
                        candidateState.physicalResidual <=
                            settings.physicalRelativeTolerance)
                    {
                        result.events.back().outcome = "rejected";
                        result.events.back().reason =
                            "stationary_multiplier_nonpositive";
                        result.events.back().detail +=
                            ";Eq. 14 residual converged but mu is nonpositive; restart required";
                        emitEvent(result.events.back());
                        stationaryNonpositiveMu = true;
                        break;
                    }
                    const bool topologyChanged =
                        !state.map.topologyHash.empty() &&
                        !candidateState.map.topologyHash.empty() &&
                        state.map.topologyHash != candidateState.map.topologyHash;
                    if (topologyChanged)
                    {
                        // A derivative is valid only within its current
                        // Boolean topology stratum.  The candidate was
                        // evaluated from its actual Boolean geometry, so it
                        // has its own fresh map/Jacobian.  Accept it only
                        // after an independently measured residual decrease;
                        // the next Newton iteration then relinearizes at the
                        // new stratum and never reuses the old derivative.
                        if (stoppingResidual(candidateState, settings) <
                            stoppingResidual(state, settings))
                        {
                            result.events.back().outcome = "accepted";
                            result.events.back().reason =
                                "topology_transition_relinearized";
                            result.events.back().detail +=
                                ";accepted residual-decreasing Boolean topology transition; next iteration uses fresh derivative";
                            emitEvent(result.events.back());
                            accepted = true;
                            break;
                        }
                        topologyTransitionRejected = true;
                        result.events.back().outcome = "rejected";
                        result.events.back().reason = "topology_transition";
                        result.events.back().detail +=
                            ";topology transition did not decrease independently evaluated residual";
                        emitEvent(result.events.back());
                    }
                    else if (stoppingResidual(candidateState, settings) <
                        stoppingResidual(state, settings))
                    {
                        result.events.back().outcome = "accepted";
                        result.events.back().reason = "residual_decreased";
                        emitEvent(result.events.back());
                        accepted = true;
                        break;
                    }
                    else
                    {
                        result.events.back().outcome = "rejected";
                        result.events.back().reason = "no_residual_decrease";
                        emitEvent(result.events.back());
                    }
                }
                else
                {
                    sawGeometryFailure = true;
                    result.diagnostic = trial.error;
                }
                step *= settings.lineSearchReduction;
            }
            if (stationaryNonpositiveMu)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic =
                    "strict Eq. 14 stationary candidate has nonpositive multiplier; restarting";
                terminate(
                    iteration,
                    "stationary_multiplier_nonpositive_restart"
                );
                break;
            }
            if (!accepted)
            {
                result.status =
                    evaluationsThisRestart >=
                        settings.maximumFunctionEvaluationsPerRestart
                    ? SpinSolverStatus::ResourceLimitReached
                    : SpinSolverStatus::LineSearchExhausted;
                if (result.diagnostic.empty() && topologyTransitionRejected)
                {
                    result.diagnostic =
                        "fixed-topology derivative guard rejected Boolean topology transition";
                }
                else if (result.diagnostic.empty())
                {
                    result.diagnostic = "Eq. 14 line search found no decrease";
                }
                terminate(
                    iteration,
                    evaluationsThisRestart >=
                            settings.maximumFunctionEvaluationsPerRestart
                        ? "maximum_evaluations_reached"
                        : "line_search_exhausted"
                );
                break;
            }

            current = candidate;
            state = candidateState;
            result.trace.push_back(traceEntry(
                restart,
                iteration,
                evaluationsThisRestart,
                step,
                state
            ));
            result.trace.back().quadricHash = current.hash();
            if(strict) result.trace.back().quadricCoefficients=coefficientArray(current);
            if (converged(state, settings))
            {
                result.status = SpinSolverStatus::Converged;
                result.quadric = current;
                result.normalizedEquationResidual = state.normalizedResidual;
                result.rawEquationResidual = state.rawResidual;
                result.physicalRelativeResidual = state.physicalResidual;
                result.multiplierMu = state.mu;
                result.multipliers = state.lambda;
                result.diagnostic.clear();
                terminate(iteration, "converged");
                finish();
                return result;
            }
            if (iteration == settings.maximumIterations)
            {
                result.status = SpinSolverStatus::ResourceLimitReached;
                result.diagnostic = "maximum Newton iterations reached";
                terminate(iteration, "maximum_iterations_reached");
            }
        }
    }

    if (result.status == SpinSolverStatus::Indeterminate)
    {
        result.status = sawGeometryFailure
            ? SpinSolverStatus::GeometryEvaluationFailed
            : SpinSolverStatus::OptimizationDidNotConverge;
    }
    if (!result.trace.empty())
    {
        const SpinIterationTrace& final = result.trace.back();
        result.normalizedEquationResidual =
            final.normalizedEquationResidual;
        result.rawEquationResidual = final.rawEquationResidual;
        result.physicalRelativeResidual = final.physicalRelativeResidual;
        result.multiplierMu = final.multiplierMu;
    }
    finish();
    return result;
}


SpinNewtonResult solveSpinMomentSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink
)
{
    return solveSpinMomentSystemImpl(problem, settings, evaluator, initialGuesses, eventSink);
}

SpinNewtonResult solveMassMomentSystem(
    const MassMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink
)
{
    return solveSpinMomentSystemImpl(problem, settings, evaluator, initialGuesses, eventSink);
}

LegacyTransactionResult solveSpinMomentSystemLegacyTransactional(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& trueEvaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const std::string& contextIdentity,
    const Vector10& fixedMomentsTarget,
    const LegacyPredictor& predictor,
    const LegacyTransactionSettings& transactionSettings
)
{
    validateLegacyTransactionSettings(transactionSettings);
    if (settings.equationProfile != SpinEquationProfile::LegacyNativeNormalizedV1)
        throw std::invalid_argument("prediction transactions require LEGACY_NATIVE_NORMALIZED_V1");
    if (contextIdentity.empty() || !fixedMomentsTarget.allFinite())
        throw std::invalid_argument("prediction transactions require bound input and finite fixed moments");
    // The disabled path is deliberately direct: not an external-evaluator shim
    // and not a second Newton implementation. No extra real evaluations occur.
    if (predictor)
        throw std::invalid_argument("prediction transaction certification is not installed yet");
    LegacyTransactionResult result;
    result.attachment.contextIdentity = contextIdentity;
    result.solver = solveSpinMomentSystemImpl(problem, settings, trueEvaluator,
                                             initialGuesses, {});
    return result;
}

std::string spinSolverScientificEventHash(const SpinNewtonResult& result)
{
    auto events = result.events;
    for (auto& event : events)
    {
        if (!event.evaluationPayload) continue;
        event.evaluationPayload->evaluatorWallSeconds = 0;
        event.evaluationPayload->residualEvaluationSeconds = 0;
        event.evaluationPayload->derivativeEvaluationSeconds = 0;
    }
    return manufacturing::computeContentSha256(
        "SpinSolverScientificEvents-v1;" + solverEventHash(
            events, result.equationProfile, result.eventStreamVersion >= 5));
}

SpinNewtonSubspaceResult solveSpinMomentSubspaceSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const Eigen::MatrixXd& supportBasis,
    const SpinMomentSubspaceEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink
)
{
    validateSettings(settings);
    if (settings.equationProfile != SpinEquationProfile::DualMonomialRawV1 || !evaluator)
        throw std::invalid_argument("supported-space Newton requires strict dual profile and evaluator");
    const Eigen::Index k = supportBasis.cols();
    if (supportBasis.rows()!=10 || k<2 || k>10 || !supportBasis.allFinite() ||
        (supportBasis.transpose()*supportBasis-Eigen::MatrixXd::Identity(k,k)).norm()>1e-12)
        throw std::invalid_argument("invalid fixed supported-space basis");
    SpinNewtonSubspaceResult result;
    result.supportBasis = supportBasis; // Freeze against caller mutation.
    std::ostringstream science;
    science.imbue(std::locale::classic());
    science << "SpinNewtonSubspace-v1;" << std::hexfloat;
    const auto matrix = [&science](const auto& value) {
        science << value.rows() << ',' << value.cols() << ';';
        for(Eigen::Index r=0;r<value.rows();++r)
            for(Eigen::Index c=0;c<value.cols();++c)science << value(r,c) << ';';
    };
    const auto string = [&science](const std::string& value) {
        science << value.size() << ':' << value << ';';
    };
    matrix(result.supportBasis);
    const std::string supportIdentity = manufacturing::computeContentSha256(science.str());
    std::optional<SpinMomentSubspaceLinearization> current;
    const SpinMomentEvaluator proxy = [&](const QuadricCoefficients10& q) {
        current.reset(); // A failed trial can never expose a previous Jacobian.
        SpinMomentSubspaceEvaluationRecord record;
        record.evaluationIndex = static_cast<int>(result.evaluations.size())+1;
        record.quadricCoefficients = coefficientArray(q);
        try
        {
            record.evaluation = evaluator(q);
            const auto& e = *record.evaluation;
            if (!e.material.momentJacobian.isZero(0))
                throw std::invalid_argument("directional evaluator must not also supply an ambient Jacobian");
            if (!e.material.derivativeDefined)
                throw std::runtime_error("directional derivative undefined: " + e.material.derivativeDiagnostic);
            current = linearizeSpinMomentSubspace(problem,q,e.material.rawMoments,
                result.supportBasis,e.nativeTangentBasis,e.nativeDirectionalJacobian);
            record.linearization = current;
            const auto material = e.material;
            result.evaluations.push_back(std::move(record));
            return material;
        }
        catch(const std::exception& error)
        {
            record.error = error.what();
            record.errorType = typeid(error).name();
            result.evaluations.push_back(std::move(record));
            throw;
        }
        catch(...)
        {
            record.error = "non-standard evaluator exception";
            record.errorType = "non-standard";
            result.evaluations.push_back(std::move(record));
            // Keep the original category in the directional attachment, then
            // enter the shared loop's normal failed-evaluation/backtrack path.
            // The legacy full-coordinate evaluator's exception contract is
            // deliberately unchanged.
            throw std::runtime_error("non-standard evaluator exception");
        }
    };
    result.solver = solveSpinMomentSystemImpl(problem,settings,proxy,initialGuesses,
        eventSink,&current,supportIdentity);
    // The ordinary capture hash includes historical timing fields. The new
    // attachment explicitly hashes scientific decisions without capture data,
    // then binds genuine directional maps/linearizations below (no timings).
    string(solverEventHash(result.solver.events,result.solver.equationProfile,false));
    string(result.solver.traceHash);
    science << static_cast<int>(result.solver.status) << ';';
    for(const auto& record:result.evaluations)
    {
        science << record.evaluationIndex << ';';
        for(double p:record.quadricCoefficients)science << p << ';';
        string(record.error);string(record.errorType);
        science << bool(record.evaluation) << ';';
        if(record.evaluation)
        {
            const auto& e=*record.evaluation;
            matrix(e.material.rawMoments);
            // Invalid returns must also have distinct evidence identities.
            // This legacy slot is required to be zero on valid directional
            // calls and is never interpreted as their derivative.
            matrix(e.material.momentJacobian);
            science << bool(e.material.variableRawMoments) << ';';
            if(e.material.variableRawMoments)matrix(*e.material.variableRawMoments);
            matrix(e.nativeTangentBasis);matrix(e.nativeDirectionalJacobian);
            science << e.material.derivativeDefined << ';';
            string(e.material.geometryHash);string(e.material.topologyHash);
            string(e.material.evaluationMapIdentity);string(e.material.booleanBackendIdentity);
            string(e.material.derivativeIdentity);string(e.material.derivativeDiagnostic);
        }
        science << bool(record.linearization) << ';';
        if(record.linearization)
        {
            const auto& l=*record.linearization;
            matrix(l.solverTangent);matrix(l.momentTangentJacobian);
            matrix(l.multipliers);science << l.multiplierMu << ';';
            matrix(l.fullResidual);science << l.rawResidual << ';';
            matrix(l.matrix);matrix(l.rightHandSide);
        }
    }
    result.scientificHash = manufacturing::computeContentSha256(science.str());
    return result;
}

std::vector<QuadricCoefficients10> deterministicInitialQuadrics(
    SpinMode mode,
    const SpinSolverSettings& settings
)
{
    validateSettings(settings);
    std::vector<QuadricCoefficients10> guesses;
    guesses.reserve(static_cast<std::size_t>(settings.maximumRestarts));
    std::uint64_t state =
        static_cast<std::uint64_t>(settings.deterministicSeed) ^
        0x5350494e49544641ULL;
    for (int restart = 0; restart < settings.maximumRestarts; ++restart)
    {
        QuadricCoefficients10::Vector values;
        for (int coordinate = 0; coordinate < values.size(); ++coordinate)
        {
            values(coordinate) = deterministicSignedUnit(&state);
        }
        switch (mode)
        {
        case SpinMode::StaticBalance:
            values(0) = std::abs(values(0)) + 0.25;
            break;
        case SpinMode::BuoyantStability:
            values.tail<6>().setZero();
            values(0) = std::abs(values(0)) + 0.25;
            break;
        case SpinMode::FidgetSpinner:
            values(3) = 0.0;
            values(6) = 0.0;
            values(7) = 0.0;
            values(8) = 0.0;
            break;
        case SpinMode::YoYo:
            values(0) = 0.0;
            values(6) = -std::abs(values(6)) - 0.25;
            break;
        case SpinMode::SpinningTop:
            values(6) = -std::abs(values(6)) - 0.1;
            break;
        }
        guesses.push_back(
            QuadricCoefficients10::fromUnnormalized(values)
        );
    }
    return guesses;
}


QuadricHollowingResult solveQuadricHollowing(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& evaluatorSettings,
    const std::vector<QuadricCoefficients10>& initialGuesses
)
{
    validateSpinDesignRequest(request);
    validateQuadricMomentEvaluatorSettings(evaluatorSettings);
    if (evaluatorSettings.derivativeMethod ==
        QuadricMomentDerivativeMethod::DiscreteReducedStratum)
    {
        throw std::invalid_argument(
            "reduced-stratum differential is evaluation-only; the Eq. 14 solver requires a full 10-D Jacobian"
        );
    }
    const std::vector<QuadricCoefficients10> guesses =
        initialGuesses.empty()
        ? deterministicInitialQuadrics(request.mode, request.settings)
        : initialGuesses;
    std::string lastDerivativeDiagnostic;
    std::string lastDerivativeScientificHash;
    const SpinMomentEvaluator evaluator =
        [&](const QuadricCoefficients10& quadric) {
            const QuadricMomentEvaluation evaluation =
                evaluateQuadricMaterialMoments(
                    request,
                    quadric,
                    evaluatorSettings
                );
            SpinMomentMapEvaluation map;
            lastDerivativeDiagnostic = evaluation.derivativeDiagnostic;
            lastDerivativeScientificHash = evaluation.derivativeScientificHash;
            map.rawMoments = evaluation.totalMomentsTarget.values();
            map.variableRawMoments = evaluation.variableMomentsTarget.values();
            map.momentJacobian = evaluation.momentJacobian;
            map.derivativeDefined = evaluation.derivativeDefined;
            map.derivativeDiagnostic = evaluation.derivativeDiagnostic;
            map.residualEvaluationSeconds =
                evaluation.residualEvaluationSeconds;
            map.derivativeEvaluationSeconds =
                evaluation.derivativeEvaluationSeconds;
            map.geometryHash = evaluation.geometryHash;
            map.topologyHash = evaluation.topologyHash;
            map.evaluationMapIdentity = evaluation.derivativeScientificHash;
            map.booleanBackendIdentity = evaluation.booleanBackendIdentity;
            map.derivativeIdentity = quadricMomentDerivativeMethodToken(
                evaluation.derivativeMethod);
            return map;
        };
    const SpinNewtonResult solved = solveSpinMomentSystem(
        problemForRequest(request),
        request.settings,
        evaluator,
        guesses
    );

    QuadricHollowingResult result;
    result.equationProfile=solved.equationProfile;
    result.rawEquationResidual=solved.rawEquationResidual;
    result.derivativeMethod = quadricMomentDerivativeMethodToken(evaluatorSettings.derivativeMethod);
    result.derivativeDiagnostic = lastDerivativeDiagnostic;
    result.derivativeScientificHash = lastDerivativeScientificHash;
    result.status = solved.status;
    result.diagnostic = solved.diagnostic;
    result.inputHash = computeQuadricEvaluatorRequestHash(request, evaluatorSettings);
    result.outerSurface = request.outerSurface;
    result.normalizedEquationResidual = solved.normalizedEquationResidual;
    result.physicalRelativeResidual = solved.physicalRelativeResidual;
    result.trace = solved.trace;
    result.traceHash = solved.traceHash;
    result.solveFunctionEvaluations = solved.functionEvaluations;
    result.solveRestartFunctionEvaluations =
        solved.restartFunctionEvaluations;
    result.solveEvaluatorWallSeconds = solved.evaluatorWallSeconds;
    result.eventStreamVersion = solved.eventStreamVersion;
    result.events = solved.events;
    result.evaluationTimings = solved.evaluationTimings;
    result.eventHash = solved.eventHash;
    result.residualEvaluationSeconds = solved.residualEvaluationSeconds;
    result.derivativeEvaluationSeconds = solved.derivativeEvaluationSeconds;
    if (solved.quadric.has_value())
    {
        const auto evaluateForMaterialExport =
            [&](const QuadricCoefficients10& quadric) {
                ++result.postSolveMaterialEvaluationCount;
                const auto started = std::chrono::steady_clock::now();
                try
                {
                    QuadricMomentEvaluation evaluation =
                        evaluateQuadricMaterialMoments(
                            request, quadric, evaluatorSettings
                        );
                    result.postSolveMaterialEvaluationWallSeconds +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - started
                        ).count();
                    result.postSolveMaterialResidualEvaluationSeconds +=
                        evaluation.residualEvaluationSeconds;
                    result.postSolveMaterialDerivativeEvaluationSeconds +=
                        evaluation.derivativeEvaluationSeconds;
                    return evaluation;
                }
                catch (...)
                {
                    result.postSolveMaterialEvaluationWallSeconds +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - started
                        ).count();
                    throw;
                }
            };
        result.quadricCoefficients = solved.quadric->values();
        try
        {
            const QuadricMomentEvaluation finalEvaluation =
                evaluateForMaterialExport(*solved.quadric);
            detail::requireNonemptyQuadricExportRegion(
                finalEvaluation, detail::QuadricExportRole::VariableMaterial);
            result.finalMoments = finalEvaluation.totalMomentsTarget;
            result.derivativeDiagnostic = finalEvaluation.derivativeDiagnostic;
            result.derivativeScientificHash = finalEvaluation.derivativeScientificHash;
            if (request.mode == SpinMode::BuoyantStability)
            {
                result.buoyantMetacenterHeight =
                    buoyantMetacenterHeight(request.buoyancy);
                result.buoyantStabilityMargin = buoyantStabilityMargin(
                    result.finalMoments.values(),
                    request.buoyancy
                );
                if (result.buoyantStabilityMargin <= 0.0)
                {
                    result.status = SpinSolverStatus::Infeasible;
                    result.diagnostic =
                        "equilibrium center of mass is not below the "
                        "prescribed-waterline metacenter";
                }
            }
            result.materialBoundary = finalEvaluation.materialBoundary;
            const QuadricCoefficients10 complement =
                QuadricCoefficients10::fromUnnormalized(
                    -solved.quadric->values()
                );
            const QuadricMomentEvaluation cavityEvaluation =
                evaluateForMaterialExport(complement);
            detail::requireNonemptyQuadricExportRegion(
                cavityEvaluation, detail::QuadricExportRole::Cavity);
            result.cavityBoundary = cavityEvaluation.materialBoundary;
            result.geometryHash = manufacturing::computeContentSha256(
                std::string("QuadricHollowingGeometry-v2;") +
                result.inputHash + ';' + finalEvaluation.geometryHash + ';' +
                cavityEvaluation.geometryHash
            );
        }
        catch (const std::exception& error)
        {
            result.status = SpinSolverStatus::InvalidMaterialArtifact;
            result.diagnostic = error.what();
        }
    }
    result.outputHash = hollowingOutputHash(result);
    return result;
}


} // namespace spin
