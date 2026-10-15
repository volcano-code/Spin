#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SlqpIterationStore.h"
#include "spin2014/SlqpAuditTrace.h"
#include "spin2014/SlqpHistory.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{


using spin2014::LbfgsHistory;


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
    const Eigen::VectorXd& direction,
    Eigen::Index betaCount
)
{
    double result = 1.0;
    for (Eigen::Index index = 0; index < betaCount; ++index)
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
    std::vector<unsigned char>* active,
    Eigen::Index betaCount
)
{
    const Eigen::MatrixXd& a = evaluation.constraintJacobian;
    active->assign(static_cast<std::size_t>(beta.size()), 0U);
    for (Eigen::Index index = 0; index < betaCount; ++index)
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
            const bool lower = index < betaCount && beta(index) <= activeTolerance;
            const bool upper = index < betaCount && beta(index) >= 1.0 - activeTolerance;
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


std::string hashTrace(const std::vector<spin2014::SlqpIteration>& trace, bool paper)
{
    std::string bytes(paper ? "spin-it-2014-paper-givens-slqp-trace-v1" : "spin-it-2014-slqp-trace-v3");
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
        if(paper) appendDouble(bytes,entry.phi);
    }
    return manufacturing::computeContentSha256(bytes);
}


void validateSettings(const spin2014::SlqpLbfgsSettings& settings)
{
    spin2014::validateSlqpDirectionSettings(settings);
    if(settings.algorithm != spin2014::SlqpAlgorithm::LegacyV3 &&
       settings.algorithm != spin2014::SlqpAlgorithm::AuditedKktRestorationV1)
        throw std::invalid_argument("unsupported SLQP algorithm");
    if(settings.algorithm == spin2014::SlqpAlgorithm::AuditedKktRestorationV1) {
        spin2014::validateSlqpAuditSettings(settings.audit);
        if(settings.activeSetTolerance>=.5)throw std::invalid_argument("overlapping active faces");
    }
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

void validateSlqpDirectionSettings(const SlqpLbfgsSettings& s) {
    if(s.directionPolicy==SlqpDirectionPolicy::ProjectedGradientV1)return;
    if((s.directionPolicy!=SlqpDirectionPolicy::MetricTangentRescueV1 &&
        s.directionPolicy!=SlqpDirectionPolicy::FiniteBoxQpV1 &&
        s.directionPolicy!=SlqpDirectionPolicy::FiniteBoxPrimalV2) ||
       s.algorithm!=SlqpAlgorithm::AuditedKktRestorationV1 ||
       ((s.metricTangent.maximumFaces<=0 || s.metricTangent.maximumFaces>100000) &&
        !(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2 && s.metricTangent.maximumFaces==-1)) ||
       s.metricTangent.maximumIterations<=0 ||
       !std::isfinite(s.metricTangent.residualTolerance) || s.metricTangent.residualTolerance<=0)
        throw std::invalid_argument("invalid metric tangent rescue controls");
}


const char* slqpStatusToken(SlqpStatus status) noexcept
{
    switch (status)
    {
    case SlqpStatus::Converged: return "CONVERGED";
    case SlqpStatus::MaximumIterations: return "MAXIMUM_ITERATIONS";
    case SlqpStatus::LineSearchExhausted: return "LINE_SEARCH_EXHAUSTED";
    case SlqpStatus::InvalidProblem: return "INVALID_PROBLEM";
    case SlqpStatus::RestorationStalled: return "RESTORATION_STALLED";
    case SlqpStatus::NumericalUnknown: return "NUMERICAL_UNKNOWN";
    }
    return "INVALID_PROBLEM";
}

const char* slqpAlgorithmToken(SlqpAlgorithm algorithm) noexcept {
    switch(algorithm) {
    case SlqpAlgorithm::LegacyV3:return "LEGACY_V3";
    case SlqpAlgorithm::AuditedKktRestorationV1:return "AUDITED_KKT_RESTORATION_V1";
    }return "UNKNOWN";
}


static SlqpLbfgsResult solveInternal(
    const EqualityBoxProblem& problem,
    const Eigen::VectorXd& initialBeta,
    const SlqpLbfgsSettings& settings,
    bool paper
)
{
    if(!problem.evaluate || problem.variableCount!=initialBeta.size())
        throw std::invalid_argument("invalid equality-box problem dimensions/callback");
    const auto& evaluate=problem.evaluate;
    const Eigen::Index betaCount=problem.betaCount;
    validateSettings(settings);
    const bool audited=settings.algorithm==SlqpAlgorithm::AuditedKktRestorationV1;
    if(!settings.iterationStorePath.empty()&&!audited)
        throw std::invalid_argument("disk iteration storage requires audited SLQP");
    if (
        initialBeta.size() == 0 || !initialBeta.allFinite() ||
        betaCount < 0 || betaCount > initialBeta.size() ||
        (initialBeta.head(betaCount).array() < 0.0).any() ||
        (initialBeta.head(betaCount).array() > 1.0).any()
    )
    {
        throw std::invalid_argument("initial fill field must lie in [0,1]");
    }

    SlqpLbfgsResult result;
    const auto last=[&]() -> SlqpIteration& {return result.iterationStore?result.iterationStore->back():result.trace.back();};
    const auto append=[&](SlqpIteration t){
        if(result.iterationStore)result.iterationStore->append(std::move(t));
        else result.trace.push_back(std::move(t));
        if(settings.iterationObserver)settings.iterationObserver(last());
    };
    result.algorithm=settings.algorithm;
    result.boundedVariableCount=betaCount;
    result.beta = initialBeta;
    Eigen::Index equalityCount=-1;
    const auto checkedEvaluate=[&](const Eigen::VectorXd& x) {
        auto e=evaluate(x);
        if(audited && (!validSlqpEvaluation(e,x.size()) ||
            (equalityCount>=0 && e.constraints.size()!=equalityCount)))
            throw std::invalid_argument("nonfinite or dimension-invalid callback evaluation");
        equalityCount=e.constraints.size();
        return e;
    };
    result.evaluation = checkedEvaluate(result.beta);
    if(!settings.iterationStorePath.empty())
        result.iterationStore=std::make_shared<SlqpIterationStore>(settings.iterationStorePath,settings.iterationStoreMaximumBytes);
    LbfgsHistory history;
    history.capacity = settings.lbfgsMemory;
    std::vector<unsigned char> previousMembership;
    const auto finish=[&]() {
        if(audited) {
            result.finalState=result.beta;
            result.finalKkt=computeKktWitness(result.evaluation,result.beta,betaCount,settings.activeSetTolerance,settings.audit);
            if(slqpTraceCount(result)==0 || !last().state.isApprox(result.beta,0.)) {
                SlqpIteration t;t.iteration=static_cast<int>(slqpTraceCount(result));t.finalState=true;
                t.objective=result.evaluation.objective;t.kkt=result.finalKkt;t.kktInfinityNorm=t.kkt.full;
                t.equalityInfinityNorm=t.kkt.equalityPrimal;
                t.meritValue=merit(result.evaluation,settings.meritPenalty);
                std::vector<unsigned char> active;
                auto p=projectedTangentGradient(result.evaluation,result.beta,settings.activeSetTolerance,&active,betaCount);
                t.projectedGradientInfinityNorm=p.lpNorm<Eigen::Infinity>();
                t.activeMembership=active;for(int i=0;i<betaCount;++i)if(active[i])t.activeMembership[i]=result.beta(i)<=settings.activeSetTolerance?1:2;
                for(auto fixed:active)t.activeVariables+=fixed!=0;
                t.membershipHash=slqpMembershipHash(t.activeMembership);
                t.phi=result.beta.size()>betaCount?result.beta(betaCount):0.;
                t.historyBefore=t.historyAfter=history.s.size();
                if(!previousMembership.empty()&&previousMembership!=t.activeMembership) {
                    t.historyResetReason="ACTIVE_MEMBERSHIP_CHANGED";history.resetForMembership(previousMembership,t.activeMembership);t.historyAfter=0;
                }
                t.lbfgsHistorySize=history.s.size();
                captureSlqpState(t,result.evaluation,result.beta);append(std::move(t));
            } else last().finalState=true;
            if(result.iterationStore)result.iterationStore->seal();
            result.traceHash=hashAuditedSlqp(result,settings);
        } else result.traceHash=hashTrace(result.trace,paper);
        return result;
    };

    for (int iteration = 0; iteration < settings.maximumIterations; ++iteration)
    {
        const AdaptiveFillEvaluation current = result.evaluation;
        const KktResidual witness=audited?computeKktWitness(current,result.beta,betaCount,settings.activeSetTolerance,settings.audit):KktResidual();
        const Eigen::VectorXd multipliers = audited ?
            (witness.valid?witness.lambda:Eigen::VectorXd::Zero(current.constraints.size()).eval()) : leastSquaresMultipliers(current);
        const Eigen::VectorXd lagrangian =
            lagrangianGradient(current, multipliers);
        std::vector<unsigned char> active;
        const Eigen::VectorXd projected = projectedTangentGradient(
            current,
            result.beta,
            settings.activeSetTolerance,
            &active,
            betaCount
        );
        const double equalityNorm = current.constraints.size() > 0
            ? current.constraints.lpNorm<Eigen::Infinity>()
            : 0.0;
        const double projectedNorm = projected.lpNorm<Eigen::Infinity>();

        SlqpIteration trace;
        trace.phi = result.beta.size() > betaCount ? result.beta(betaCount) : 0.;
        trace.iteration = iteration;
        trace.objective = current.objective;
        trace.equalityInfinityNorm = equalityNorm;
        trace.projectedGradientInfinityNorm = projectedNorm;
        trace.kktInfinityNorm = projectedNorm;
        trace.meritValue = merit(current, settings.meritPenalty);
        trace.lbfgsHistorySize = static_cast<int>(history.s.size());
        if(audited) {
            trace.kkt=witness;trace.kktInfinityNorm=witness.full;
            captureSlqpState(trace,current,result.beta);trace.activeMembership=active;
            for(int i=0;i<betaCount;++i)if(active[i])trace.activeMembership[i]=result.beta(i)<=settings.activeSetTolerance?1:2;
            trace.membershipHash=slqpMembershipHash(trace.activeMembership);
            trace.historyBefore=history.s.size();
            if(iteration==0 || previousMembership!=trace.activeMembership) {
                trace.historyResetReason=iteration==0?"NEW_SOLVE_OR_PARTITION":"ACTIVE_MEMBERSHIP_CHANGED";
                history.resetForMembership(previousMembership,trace.activeMembership);
            }
            trace.historyAfter=history.s.size();previousMembership=trace.activeMembership;
            trace.lbfgsHistorySize=history.s.size();
            if(!witness.valid) {
                append(std::move(trace));result.status=SlqpStatus::NumericalUnknown;
                result.diagnostic=witness.status;return finish();
            }
        }
        for (unsigned char fixed : active)
        {
            trace.activeVariables += fixed != 0U;
        }

        if (
            equalityNorm <= settings.equalityTolerance &&
            (audited ? (witness.valid && witness.boundPrimal<=settings.projectedGradientTolerance &&
                witness.stationarity<=settings.projectedGradientTolerance && witness.dualSign<=settings.projectedGradientTolerance &&
                witness.complementarity<=settings.projectedGradientTolerance) : projectedNorm <= settings.projectedGradientTolerance)
        )
        {
            append(trace);
            result.status = SlqpStatus::Converged;
            result.diagnostic = audited?"FIRST_ORDER_KKT_TOLERANCES_MET; no minimum or physical certificate":"SLQP equality and projected-gradient tolerances met";
            return finish();
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

        trace.linearizedConstraintResidual=(current.constraintJacobian*direction+current.constraints).lpNorm<1>();
        const bool infeasible=equalityNorm>settings.equalityTolerance;
        bool needRestoration=audited && infeasible && (!direction.allFinite()||direction.norm()<=1e-15 ||
            !std::isfinite(trace.linearizedConstraintResidual)||trace.linearizedConstraintResidual>settings.equalityTolerance);

        if ((!direction.allFinite() || direction.norm() <= 1e-15) && !needRestoration &&
            settings.directionPolicy!=SlqpDirectionPolicy::FiniteBoxQpV1 &&
            settings.directionPolicy!=SlqpDirectionPolicy::FiniteBoxPrimalV2)
        {
            append(trace);
            result.status = SlqpStatus::LineSearchExhausted;
            result.diagnostic = "SLQP subproblem produced no usable direction";
            return finish();
        }

        double step = needRestoration?0.:maximumBoundStep(result.beta, direction, betaCount);
        // Independent opt-in primary QP, including finite displacement bounds.
        // A failed QP stops UNKNOWN rather than silently selecting another QP.
        if(audited && (settings.directionPolicy==SlqpDirectionPolicy::FiniteBoxQpV1 ||
                       settings.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2)) {
            auto controls=settings.metricTangent;controls.activeTolerance=settings.activeSetTolerance;
            const auto solve=settings.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2 ? solvePrimalBoxQp : solveFiniteBoxQp;
            trace.metricTangent=solve(result.beta,current.gradient,
                current.constraintJacobian,current.constraints,betaCount,history,controls);
            if(!trace.metricTangent.checked) {
                append(std::move(trace));result.status=SlqpStatus::NumericalUnknown;
                result.diagnostic="finite-box QP was not verified; UNKNOWN, no infeasibility certificate";
                return finish();
            }
            direction=trace.metricTangent.direction;
            step=maximumBoundStep(result.beta,direction,betaCount);
            needRestoration=false;
            trace.metricTangentSelected=true;
            trace.linearizedConstraintResidual=(current.constraintJacobian*direction+current.constraints).lpNorm<1>();
        }
        // Opt-in rescue of the measured zero-box-step failure. The original
        // identity-selected membership remains separately auditable. The QP
        // uses the full SPD metric and its own add/release working set.
        if(audited && settings.directionPolicy==SlqpDirectionPolicy::MetricTangentRescueV1 &&
           !infeasible && !needRestoration && step<=1e-14) {
            auto controls=settings.metricTangent;controls.activeTolerance=settings.activeSetTolerance;
            trace.metricTangent=solveTangentQp(result.beta,current.gradient,
                current.constraintJacobian,current.constraints,betaCount,history,controls);
            if(trace.metricTangent.checked && trace.metricTangent.direction.size()==result.beta.size() &&
               trace.metricTangent.direction.allFinite()) {
                const double trialStep=maximumBoundStep(result.beta,trace.metricTangent.direction,betaCount);
                if(trialStep>1e-14) {
                    direction=trace.metricTangent.direction;step=trialStep;
                    trace.metricTangentSelected=true;
                    trace.linearizedConstraintResidual=(current.constraintJacobian*direction+current.constraints).lpNorm<1>();
                }
            }
        }
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
            candidateBeta = result.beta + step * direction;
            candidateBeta.head(betaCount) = candidateBeta.head(betaCount).cwiseMax(0.0).cwiseMin(1.0);
            try {candidate = checkedEvaluate(candidateBeta);} catch(...) {if(!audited)throw;trace.objectiveTrials.push_back({step,INFINITY,false,false,"INVALID_EVALUATION_OR_CALLBACK_EXCEPTION"});step*=settings.lineSearchReduction;continue;}
            if(audited)trace.objectiveTrials.push_back({step,candidate.constraints.lpNorm<1>(),true,false,"INSUFFICIENT_MERIT_DECREASE"});
            if (
                merit(candidate, settings.meritPenalty) <=
                oldMerit - settings.armijo * step * predicted
            )
            {
                accepted = true;
                if(audited){trace.objectiveTrials.back().accepted=true;trace.objectiveTrials.back().reason="ACCEPTED";}
                break;
            }
            step *= settings.lineSearchReduction;
        }
        if(trace.metricTangentSelected)
            trace.metricTangentTrialCount=static_cast<int>(trace.objectiveTrials.size());
        if (!accepted && !(audited&&infeasible))
        {
            trace.metricTangentFallbackAttempted=trace.metricTangent.attempted;
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
            step = maximumBoundStep(result.beta, direction, betaCount);
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
                candidateBeta = result.beta + step * direction;
                candidateBeta.head(betaCount) = candidateBeta.head(betaCount).cwiseMax(0.0).cwiseMin(1.0);
                try {candidate = checkedEvaluate(candidateBeta);} catch(...) {if(!audited)throw;trace.objectiveTrials.push_back({step,INFINITY,false,false,"INVALID_EVALUATION_OR_CALLBACK_EXCEPTION"});step*=settings.lineSearchReduction;continue;}
                if(audited)trace.objectiveTrials.push_back({step,candidate.constraints.lpNorm<1>(),true,false,"INSUFFICIENT_MERIT_DECREASE"});
                if (
                    merit(candidate, settings.meritPenalty) <=
                    oldMerit - settings.armijo * step * fallbackPredicted
                )
                {
                    accepted = true;
                    if(audited) {trace.objectiveTrials.back().accepted=true;trace.objectiveTrials.back().reason="ACCEPTED";trace.historyResetReason+=";FIRST_ORDER_FALLBACK";}
                    trace.usedFirstOrderFallback = true;
                    history = LbfgsHistory();
                    history.capacity = settings.lbfgsMemory;
                    break;
                }
                step *= settings.lineSearchReduction;
            }
        }
        if(audited && !accepted && infeasible) {
            trace.restoration=solveRestoration(current,result.beta,betaCount,settings.audit);
            trace.restoration.reason=needRestoration?"UNUSABLE_OR_INACCURATE_EQUALITY_DIRECTION":"MERIT_SEARCH_FAILED_WHILE_INFEASIBLE";
            const double oldResidual=current.constraints.lpNorm<1>();
            const double improvement=oldResidual-trace.restoration.predictedResidual;
            if(trace.restoration.direction.size()==result.beta.size() && trace.restoration.direction.allFinite() &&
               std::isfinite(improvement)&&improvement>0.) {
                direction=trace.restoration.direction;step=1.;
                for(int lineSearch=0;lineSearch<settings.maximumLineSearchSteps&&step>1e-14;++lineSearch) {
                    ++trace.lineSearchEvaluations;
                    candidateBeta=result.beta+step*direction;
                    bool valid=true;
                    try {candidate=checkedEvaluate(candidateBeta);}catch(...) {valid=false;}
                    const double residual=valid?candidate.constraints.lpNorm<1>():INFINITY;
                    const bool sufficient=valid&&std::isfinite(residual)&&residual<oldResidual &&
                        residual<=oldResidual-settings.armijo*step*improvement;
                    trace.restoration.trials.push_back({step,residual,valid,sufficient,
                        sufficient?"ACCEPTED":valid?"INSUFFICIENT_ACTUAL_CONSTRAINT_DECREASE":"INVALID_EVALUATION_OR_CALLBACK_EXCEPTION"});
                    if(sufficient) {accepted=true;trace.restoration.accepted=true;break;}
                    step*=settings.lineSearchReduction;
                }
            }
            if(accepted) {
                history=LbfgsHistory();history.capacity=settings.lbfgsMemory;
                trace.historyResetReason+=";RESTORATION_ACCEPTED";
            } else {
                result.status=trace.restoration.status=="Optimal"?SlqpStatus::RestorationStalled:SlqpStatus::NumericalUnknown;
                result.diagnostic="local restoration made no accepted progress; UNKNOWN, no infeasibility certificate";
            }
        }
        trace.acceptedStep = accepted ? step : 0.0;
        trace.metricTangentAccepted=trace.metricTangentSelected && accepted &&
            !trace.usedFirstOrderFallback && !trace.restoration.accepted;
        if (!accepted)
        {
            append(trace);
            if(!trace.restoration.attempted) {
                result.status = SlqpStatus::LineSearchExhausted;
                result.diagnostic = "SLQP merit-function line search was exhausted";
            }return finish();
        }

        const Eigen::VectorXd nextMultipliers = audited?multipliers:leastSquaresMultipliers(candidate);
        const Eigen::VectorXd stepVector = candidateBeta - result.beta;
        const Eigen::VectorXd gradientChange =
            lagrangianGradient(candidate, nextMultipliers) - lagrangian;
        if(!trace.restoration.accepted) {
            trace.curvaturePairAccepted=history.update(stepVector, gradientChange);
            trace.curvaturePairRejected=!trace.curvaturePairAccepted;
        }
        trace.historyAfter=history.s.size();append(std::move(trace));
        result.beta = std::move(candidateBeta);
        result.evaluation = std::move(candidate);
    }

    result.status = SlqpStatus::MaximumIterations;
    result.diagnostic = "SLQP reached its iteration limit";
    return finish();
}

SlqpLbfgsResult solveSlqpLbfgs(const AdaptiveFillFunctional& functional,
    const Eigen::VectorXd& initialBeta,const SlqpLbfgsSettings& settings)
{
    if(settings.algorithm!=SlqpAlgorithm::LegacyV3)
        throw std::invalid_argument("audited algorithm requires PAPER_GIVENS at spin interfaces");
    return solveInternal({initialBeta.size(),initialBeta.size(),[&functional](const Eigen::VectorXd& beta){return functional.evaluate(beta);}},initialBeta,settings,false);
}

SlqpLbfgsResult solveSlqpLbfgs(const EqualityBoxProblem& problem,
    const Eigen::VectorXd& initial,const SlqpLbfgsSettings& settings)
{
    return solveInternal(problem, initial, settings, false);
}

SlqpLbfgsResult solveSlqpLbfgsPaper(const AdaptiveFillFunctional& functional,
    const Eigen::VectorXd& initialBeta,double initialPhi,const SlqpLbfgsSettings& settings)
{
    const Eigen::Index n=functional.betaCount();const bool angle=functional.paperHasAngle();
    if(initialBeta.size()!=n||!std::isfinite(initialPhi)||(!angle&&initialPhi!=0.))
        throw std::invalid_argument("invalid joint Givens dimensions or initial phi (static requires zero)");
    Eigen::VectorXd state(n+(angle?1:0));state.head(n)=initialBeta;if(angle)state(n)=initialPhi;
    auto result=solveInternal({state.size(),n,[&](const Eigen::VectorXd& x) {
        if(!x.allFinite())throw std::invalid_argument("nonfinite joint Givens state");
        return functional.evaluatePaper(x.head(n),angle?x(n):0.);
    }},state,settings,true);
    result.formulation=Spin2014Formulation::PaperGivens;
    result.initialPhi=initialPhi;result.finalPhi=angle?result.beta(n):0.;
    result.beta=result.beta.head(n).eval();
    // Bind final accepted state as well: a max-iteration trace ends before its
    // last accepted step. Legacy serialization is deliberately unchanged.
    std::string bytes=result.traceHash;appendDouble(bytes,result.finalPhi);
    for(Eigen::Index i=0;i<n;++i)appendDouble(bytes,result.beta(i));
    result.traceHash=manufacturing::computeContentSha256(bytes);
    return result;
}


} // namespace spin2014
