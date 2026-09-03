#include "spin2014/SlqpLbfgsSolver.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{


struct LbfgsHistory
{
    int capacity = 0;
    std::deque<Eigen::VectorXd> s;
    std::deque<Eigen::VectorXd> y;
    std::deque<double> inverseCurvature;

    void update(const Eigen::VectorXd& step, const Eigen::VectorXd& gradientChange)
    {
        const double curvature = step.dot(gradientChange);
        if (
            !std::isfinite(curvature) ||
            curvature <= 1e-12 * step.norm() * gradientChange.norm()
        )
        {
            return;
        }
        if (static_cast<int>(s.size()) == capacity)
        {
            s.pop_front();
            y.pop_front();
            inverseCurvature.pop_front();
        }
        s.push_back(step);
        y.push_back(gradientChange);
        inverseCurvature.push_back(1.0 / curvature);
    }

    Eigen::VectorXd applyInverse(const Eigen::VectorXd& vector) const
    {
        if (s.empty())
        {
            return vector;
        }
        Eigen::VectorXd result = vector;
        std::vector<double> alpha(s.size(), 0.0);
        for (std::size_t offset = 0; offset < s.size(); ++offset)
        {
            const std::size_t index = s.size() - 1U - offset;
            alpha[index] = inverseCurvature[index] * s[index].dot(result);
            result -= alpha[index] * y[index];
        }
        const double yy = y.back().squaredNorm();
        const double scale = yy > 0.0
            ? s.back().dot(y.back()) / yy
            : 1.0;
        result *= std::clamp(scale, 1e-6, 1e6);
        for (std::size_t index = 0; index < s.size(); ++index)
        {
            const double beta =
                inverseCurvature[index] * y[index].dot(result);
            result += s[index] * (alpha[index] - beta);
        }
        return result;
    }
};


Eigen::VectorXd solveSymmetric(
    const Eigen::MatrixXd& matrix,
    const Eigen::VectorXd& rightHandSide
)
{
    if (matrix.rows() == 0)
    {
        return Eigen::VectorXd();
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(matrix);
    decomposition.setThreshold(1e-12);
    return decomposition.solve(rightHandSide);
}


Eigen::VectorXd leastSquaresMultipliers(
    const spin2014::AdaptiveFillEvaluation& evaluation
)
{
    const Eigen::MatrixXd& a = evaluation.constraintJacobian;
    if (a.rows() == 0)
    {
        return Eigen::VectorXd();
    }
    return solveSymmetric(
        a * a.transpose(),
        -a * evaluation.gradient
    );
}


Eigen::VectorXd lagrangianGradient(
    const spin2014::AdaptiveFillEvaluation& evaluation,
    const Eigen::VectorXd& multipliers
)
{
    return evaluation.gradient +
        evaluation.constraintJacobian.transpose() * multipliers;
}


double merit(
    const spin2014::AdaptiveFillEvaluation& evaluation,
    double penalty
)
{
    return evaluation.objective + penalty * evaluation.constraints.lpNorm<1>();
}


double maximumBoundStep(
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& direction
)
{
    double result = 1.0;
    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        if (direction(index) > 0.0)
        {
            result = std::min(result, (1.0 - beta(index)) / direction(index));
        }
        else if (direction(index) < 0.0)
        {
            result = std::min(result, -beta(index) / direction(index));
        }
    }
    return std::clamp(result, 0.0, 1.0);
}


Eigen::VectorXd projectedTangentGradient(
    const spin2014::AdaptiveFillEvaluation& evaluation,
    const Eigen::VectorXd& beta,
    double activeTolerance,
    std::vector<unsigned char>* active
)
{
    const Eigen::MatrixXd& a = evaluation.constraintJacobian;
    active->assign(static_cast<std::size_t>(beta.size()), 0U);
    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        if (
            beta(index) <= activeTolerance ||
            beta(index) >= 1.0 - activeTolerance
        )
        {
            (*active)[static_cast<std::size_t>(index)] = 1U;
        }
    }
    Eigen::VectorXd projected = evaluation.gradient;
    for (Eigen::Index pass = 0; pass < 4 * beta.size() + 1; ++pass)
    {
        Eigen::MatrixXd freeA = a;
        for (Eigen::Index index = 0; index < beta.size(); ++index)
        {
            if ((*active)[static_cast<std::size_t>(index)] != 0U)
            {
                freeA.col(index).setZero();
            }
        }
        projected = evaluation.gradient;
        if (freeA.rows() > 0)
        {
            projected += a.transpose() * solveSymmetric(
                freeA * freeA.transpose(), -freeA * projected
            );
        }
        Eigen::Index change = -1;
        double largestViolation = activeTolerance;
        bool release = false;
        for (Eigen::Index index = 0; index < beta.size(); ++index)
        {
            const bool lower = beta(index) <= activeTolerance;
            const bool upper = beta(index) >= 1.0 - activeTolerance;
            if (!lower && !upper) continue;
            const double direction = -projected(index);
            const bool fixed = (*active)[static_cast<std::size_t>(index)] != 0U;
            const bool pointsInward =
                (lower && direction > activeTolerance) ||
                (upper && direction < -activeTolerance);
            const bool pointsOutward =
                (lower && direction < -activeTolerance) ||
                (upper && direction > activeTolerance);
            const bool violation = fixed ? pointsInward : pointsOutward;
            if (violation && std::abs(direction) > largestViolation)
            {
                largestViolation = std::abs(direction);
                change = index;
                release = fixed;
            }
        }
        if (change < 0) break;
        (*active)[static_cast<std::size_t>(change)] = release ? 0U : 1U;
    }
    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        if ((*active)[static_cast<std::size_t>(index)] != 0U)
        {
            projected(index) = 0.0;
        }
    }
    return projected;
}


void appendDouble(std::string& bytes, double value)
{
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int shift = 0; shift < 64; shift += 8)
    {
        bytes.push_back(static_cast<char>((bits >> shift) & 0xffU));
    }
}


std::string hashTrace(const std::vector<spin2014::SlqpIteration>& trace)
{
    std::string bytes("spin-it-2014-slqp-trace-v3");
    for (const spin2014::SlqpIteration& entry : trace)
    {
        bytes.append(reinterpret_cast<const char*>(&entry.iteration), sizeof(int));
        appendDouble(bytes, entry.objective);
        appendDouble(bytes, entry.equalityInfinityNorm);
        appendDouble(bytes, entry.projectedGradientInfinityNorm);
        appendDouble(bytes, entry.kktInfinityNorm);
        appendDouble(bytes, entry.meritValue);
        appendDouble(bytes, entry.acceptedStep);
        bytes.append(
            reinterpret_cast<const char*>(&entry.lineSearchEvaluations),
            sizeof(int)
        );
        bytes.append(
            reinterpret_cast<const char*>(&entry.activeVariables), sizeof(int)
        );
        bytes.append(
            reinterpret_cast<const char*>(&entry.lbfgsHistorySize), sizeof(int)
        );
        bytes.push_back(entry.usedFirstOrderFallback ? '\1' : '\0');
    }
    return manufacturing::computeContentSha256(bytes);
}


void validateSettings(const spin2014::SlqpLbfgsSettings& settings)
{
    if (
        settings.maximumIterations <= 0 ||
        settings.maximumLineSearchSteps <= 0 ||
        settings.lbfgsMemory <= 0 ||
        !std::isfinite(settings.equalityTolerance) ||
        settings.equalityTolerance <= 0.0 ||
        !std::isfinite(settings.projectedGradientTolerance) ||
        settings.projectedGradientTolerance <= 0.0 ||
        !std::isfinite(settings.activeSetTolerance) ||
        settings.activeSetTolerance < 0.0 ||
        !std::isfinite(settings.meritPenalty) ||
        settings.meritPenalty <= 0.0 ||
        !std::isfinite(settings.lineSearchReduction) ||
        settings.lineSearchReduction <= 0.0 ||
        settings.lineSearchReduction >= 1.0 ||
        !std::isfinite(settings.armijo) ||
        settings.armijo <= 0.0 || settings.armijo >= 1.0
    )
    {
        throw std::invalid_argument("invalid Spin-It 2014 SLQP settings");
    }
}


} // namespace


namespace spin2014
{


const char* slqpStatusToken(SlqpStatus status) noexcept
{
    switch (status)
    {
    case SlqpStatus::Converged: return "CONVERGED";
    case SlqpStatus::MaximumIterations: return "MAXIMUM_ITERATIONS";
    case SlqpStatus::LineSearchExhausted: return "LINE_SEARCH_EXHAUSTED";
    case SlqpStatus::InvalidProblem: return "INVALID_PROBLEM";
    }
    return "INVALID_PROBLEM";
}


SlqpLbfgsResult solveSlqpLbfgs(
    const AdaptiveFillFunctional& functional,
    const Eigen::VectorXd& initialBeta,
    const SlqpLbfgsSettings& settings
)
{
    validateSettings(settings);
    if (
        initialBeta.size() == 0 || !initialBeta.allFinite() ||
        (initialBeta.array() < 0.0).any() ||
        (initialBeta.array() > 1.0).any()
    )
    {
        throw std::invalid_argument("initial fill field must lie in [0,1]");
    }

    SlqpLbfgsResult result;
    result.beta = initialBeta;
    result.evaluation = functional.evaluate(result.beta);
    LbfgsHistory history;
    history.capacity = settings.lbfgsMemory;

    for (int iteration = 0; iteration < settings.maximumIterations; ++iteration)
    {
        const AdaptiveFillEvaluation current = result.evaluation;
        const Eigen::VectorXd multipliers = leastSquaresMultipliers(current);
        const Eigen::VectorXd lagrangian =
            lagrangianGradient(current, multipliers);
        std::vector<unsigned char> active;
        const Eigen::VectorXd projected = projectedTangentGradient(
            current,
            result.beta,
            settings.activeSetTolerance,
            &active
        );
        const double equalityNorm = current.constraints.size() > 0
            ? current.constraints.lpNorm<Eigen::Infinity>()
            : 0.0;
        const double projectedNorm = projected.lpNorm<Eigen::Infinity>();

        SlqpIteration trace;
        trace.iteration = iteration;
        trace.objective = current.objective;
        trace.equalityInfinityNorm = equalityNorm;
        trace.projectedGradientInfinityNorm = projectedNorm;
        trace.kktInfinityNorm = projectedNorm;
        trace.meritValue = merit(current, settings.meritPenalty);
        trace.lbfgsHistorySize = static_cast<int>(history.s.size());
        for (unsigned char fixed : active)
        {
            trace.activeVariables += fixed != 0U;
        }

        if (
            equalityNorm <= settings.equalityTolerance &&
            projectedNorm <= settings.projectedGradientTolerance
        )
        {
            result.trace.push_back(trace);
            result.status = SlqpStatus::Converged;
            result.diagnostic = "SLQP equality and projected-gradient tolerances met";
            result.traceHash = hashTrace(result.trace);
            return result;
        }

        const auto applyFreeInverse = [&history, &active](
            const Eigen::VectorXd& vector
        )
        {
            Eigen::VectorXd freeVector = vector;
            for (Eigen::Index index = 0; index < freeVector.size(); ++index)
            {
                if (active[static_cast<std::size_t>(index)] != 0U)
                {
                    freeVector(index) = 0.0;
                }
            }
            Eigen::VectorXd result = history.applyInverse(freeVector);
            for (Eigen::Index index = 0; index < result.size(); ++index)
            {
                if (active[static_cast<std::size_t>(index)] != 0U)
                {
                    result(index) = 0.0;
                }
            }
            return result;
        };
        Eigen::VectorXd inverseGradient = applyFreeInverse(current.gradient);
        Eigen::MatrixXd inverseConstraint(
            result.beta.size(), current.constraints.size()
        );
        for (Eigen::Index row = 0; row < current.constraints.size(); ++row)
        {
            inverseConstraint.col(row) = applyFreeInverse(
                current.constraintJacobian.row(row).transpose()
            );
        }
        const Eigen::MatrixXd schur =
            current.constraintJacobian * inverseConstraint;
        const Eigen::VectorXd qpMultipliers = solveSymmetric(
            schur,
            current.constraints -
                current.constraintJacobian * inverseGradient
        );
        Eigen::VectorXd direction =
            -inverseGradient - inverseConstraint * qpMultipliers;

        if (!direction.allFinite() || direction.norm() <= 1e-15)
        {
            result.trace.push_back(trace);
            result.status = SlqpStatus::LineSearchExhausted;
            result.diagnostic = "SLQP subproblem produced no usable direction";
            result.traceHash = hashTrace(result.trace);
            return result;
        }

        double step = maximumBoundStep(result.beta, direction);
        const double oldMerit = merit(current, settings.meritPenalty);
        const double predicted = std::max(
            1e-14,
            -current.gradient.dot(direction) +
                settings.meritPenalty * current.constraints.lpNorm<1>()
        );
        bool accepted = false;
        Eigen::VectorXd candidateBeta;
        AdaptiveFillEvaluation candidate;
        for (int lineSearch = 0;
             lineSearch < settings.maximumLineSearchSteps && step > 1e-14;
             ++lineSearch)
        {
            ++trace.lineSearchEvaluations;
            candidateBeta = (result.beta + step * direction)
                .cwiseMax(0.0).cwiseMin(1.0);
            candidate = functional.evaluate(candidateBeta);
            if (
                merit(candidate, settings.meritPenalty) <=
                oldMerit - settings.armijo * step * predicted
            )
            {
                accepted = true;
                break;
            }
            step *= settings.lineSearchReduction;
        }
        if (!accepted)
        {
            // A stale L-BFGS metric can cease to be a descent metric after a
            // large active-set change.  Retry the same linearized equality
            // subproblem with the identity metric.  This is still an SLQP
            // step; it merely discards the unreliable quasi-Newton scaling.
            Eigen::MatrixXd freeA = current.constraintJacobian;
            for (Eigen::Index index = 0; index < result.beta.size(); ++index)
            {
                if (active[static_cast<std::size_t>(index)] != 0U)
                {
                    freeA.col(index).setZero();
                }
            }
            Eigen::VectorXd normal = Eigen::VectorXd::Zero(result.beta.size());
            if (freeA.rows() > 0)
            {
                normal = -freeA.transpose() * solveSymmetric(
                    freeA * freeA.transpose(), current.constraints
                );
            }
            direction = normal - projected;
            for (Eigen::Index index = 0; index < result.beta.size(); ++index)
            {
                if (active[static_cast<std::size_t>(index)] != 0U)
                {
                    direction(index) = 0.0;
                }
            }
            step = maximumBoundStep(result.beta, direction);
            const double fallbackPredicted = std::max(
                1e-14,
                -current.gradient.dot(direction) +
                    settings.meritPenalty * current.constraints.lpNorm<1>()
            );
            for (int lineSearch = 0;
                 lineSearch < settings.maximumLineSearchSteps && step > 1e-14;
                 ++lineSearch)
            {
                ++trace.lineSearchEvaluations;
                candidateBeta = (result.beta + step * direction)
                    .cwiseMax(0.0).cwiseMin(1.0);
                candidate = functional.evaluate(candidateBeta);
                if (
                    merit(candidate, settings.meritPenalty) <=
                    oldMerit - settings.armijo * step * fallbackPredicted
                )
                {
                    accepted = true;
                    trace.usedFirstOrderFallback = true;
                    history = LbfgsHistory();
                    history.capacity = settings.lbfgsMemory;
                    break;
                }
                step *= settings.lineSearchReduction;
            }
        }
        trace.acceptedStep = accepted ? step : 0.0;
        result.trace.push_back(trace);
        if (!accepted)
        {
            result.status = SlqpStatus::LineSearchExhausted;
            result.diagnostic = "SLQP merit-function line search was exhausted";
            result.traceHash = hashTrace(result.trace);
            return result;
        }

        const Eigen::VectorXd nextMultipliers = leastSquaresMultipliers(candidate);
        const Eigen::VectorXd stepVector = candidateBeta - result.beta;
        const Eigen::VectorXd gradientChange =
            lagrangianGradient(candidate, nextMultipliers) - lagrangian;
        history.update(stepVector, gradientChange);
        result.beta = std::move(candidateBeta);
        result.evaluation = std::move(candidate);
    }

    result.status = SlqpStatus::MaximumIterations;
    result.diagnostic = "SLQP reached its iteration limit";
    result.traceHash = hashTrace(result.trace);
    return result;
}


} // namespace spin2014
