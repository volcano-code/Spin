#include "spin/LegacyPredictionTransaction.h"
#include "spin/QuadricHollowingSolver.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
spin::LegacyPredictionRequest request() {
    spin::LegacyPredictionRequest r;
    r.contextIdentity = "actual-input-and-evaluator-v1";
    r.checkpointIdentity = "true-checkpoint-v1";
    r.predictorIdentity = "pure-predictor-v1";
    r.stateVersion = 1;
    r.transactionId = 1;
    r.tentativeStep = 1;
    r.q.setZero(); r.q(0) = 1;
    r.stateIdentity = spin::legacyPredictionStateIdentity(r);
    return r;
}
spin::LegacyPredictionResponse response(const spin::LegacyPredictionRequest& r) {
    spin::LegacyPredictionResponse v;
    v.valid = true;
    v.variableMoments = r.q;
    v.operationalDirectionMatrix.setIdentity();
    v.contextIdentity = r.contextIdentity;
    v.stateIdentity = r.stateIdentity;
    v.predictorIdentity = r.predictorIdentity;
    return v;
}
void contract() {
    for (int size : {0, 4}) {
        spin::LegacyTransactionSettings s; s.maximumTentativeSteps = size;
        bool rejected = false;
        try { spin::validateLegacyTransactionSettings(s); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid tentative block length admitted");
    }
    for (int size : {1, 2, 3}) {
        spin::LegacyTransactionSettings s; s.maximumTentativeSteps = size;
        spin::validateLegacyTransactionSettings(s);
    }
    const auto r = request();
    const auto v = response(r);
    require(spin::validateLegacyPrediction(r, v).empty(), "valid response rejected");
    auto changed = r; ++changed.stateVersion;
    require(spin::legacyPredictionStateIdentity(changed) != r.stateIdentity,
            "same q at another state version reused checkpoint identity");
    changed = r; changed.contextIdentity += "-other-input";
    require(spin::legacyPredictionStateIdentity(changed) != r.stateIdentity,
            "different input reused output-only identity");
    changed = r; changed.checkpointIdentity += "-other-checkpoint";
    require(spin::legacyPredictionStateIdentity(changed) != r.stateIdentity,
            "checkpoint omitted from prediction identity");
    auto wrong = v; wrong.stateIdentity += "-stale";
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "stale state accepted");
    wrong = v; wrong.variableMoments(2) = std::numeric_limits<double>::quiet_NaN();
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "nonfinite moments accepted");
    wrong = v; wrong.operationalDirectionMatrix(1, 2) = std::numeric_limits<double>::infinity();
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "nonfinite operation matrix accepted");
    wrong = v; wrong.reportedSeconds = -1;
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "invalid reported time accepted");
    changed = r; changed.equationProfile = spin::SpinEquationProfile::DualMonomialRawV1;
    require(!spin::validateLegacyPrediction(changed, v).empty(), "DUAL admitted operational matrix");
    changed = r; changed.q(0) = std::numeric_limits<double>::quiet_NaN();
    require(!spin::validateLegacyPrediction(changed, v).empty(), "nonfinite request admitted");
    changed = r; changed.contextIdentity.clear();
    require(!spin::validateLegacyPrediction(changed, v).empty(), "missing input identity admitted");
    changed = r; changed.q *= 2;
    require(!spin::validateLegacyPrediction(changed, v).empty(), "unnormalized coefficients admitted");
    changed = r; changed.tentativeStep = 2;
    require(!spin::validateLegacyPrediction(changed, v).empty(), "missing tentative predecessor admitted");
    changed.previousTentativeIdentity = "first-tentative-state";
    changed.stateIdentity = spin::legacyPredictionStateIdentity(changed);
    require(spin::validateLegacyPrediction(changed, response(changed)).empty(),
            "valid tentative predecessor rejected");
    wrong = v; wrong.equationProfile = spin::SpinEquationProfile::DualMonomialRawV1;
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "DUAL response admitted");
    wrong = v; wrong.contextIdentity += "-other";
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "wrong input response admitted");
    wrong = v; wrong.predictorIdentity.clear();
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "anonymous predictor admitted");
    wrong = v; wrong.reportedSeconds = std::numeric_limits<double>::quiet_NaN();
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "nonfinite reported time admitted");
    wrong = v; wrong.schemaVersion = 2;
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "unknown schema admitted");
    changed = r; changed.stateIdentity = r.checkpointIdentity;
    require(!spin::validateLegacyPrediction(changed, v).empty(), "checkpoint reused as tentative identity");
    wrong = v; wrong.valid = false;
    require(!spin::validateLegacyPrediction(r, wrong).empty(), "invalid payload admitted");
}
spin::QuadricCoefficients10 initialPoint() {
    spin::Vector10 p;
    p << 1.0, 0.12, -0.08, 0.09, -0.07, 0.05, -0.04, 0.03, -0.02, 0.01;
    return spin::QuadricCoefficients10::fromUnnormalized(p);
}
spin::SpinMomentMapEvaluation identityMap(const spin::QuadricCoefficients10& q) {
    spin::SpinMomentMapEvaluation value;
    value.rawMoments = q.values();
    value.momentJacobian.setIdentity();
    value.derivativeDefined = true;
    value.geometryHash = q.hash();
    value.topologyHash = "synthetic-identity";
    return value;
}
void disabled() {
    const auto problem = spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    spin::SpinSolverSettings settings;
    for (bool payload : {false, true}) {
        settings.captureEvaluationPayload = payload;
        std::vector<std::string> baseCalls, wrappedCalls;
        const auto base = spin::solveSpinMomentSystem(problem, settings,
            [&](const auto& q) { baseCalls.push_back(q.hash()); return identityMap(q); },
            {initialPoint()});
        const auto wrapped = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
            [&](const auto& q) { wrappedCalls.push_back(q.hash()); return identityMap(q); },
            {initialPoint()}, "kernel-fixture-v1", spin::Vector10::Zero(), {},
            spin::LegacyTransactionSettings{});
        require(base.status == spin::SpinSolverStatus::Converged && base.trace.size() > 1,
                "fixture did not exercise noninitial convergence");
        require(base.status == wrapped.solver.status, "disabled status changed");
        require(base.quadric->values() == wrapped.solver.quadric->values(), "disabled q changed");
        require(base.multiplierMu == wrapped.solver.multiplierMu &&
                base.multipliers == wrapped.solver.multipliers, "disabled multipliers changed");
        require(base.normalizedEquationResidual == wrapped.solver.normalizedEquationResidual &&
                base.rawEquationResidual == wrapped.solver.rawEquationResidual &&
                base.physicalRelativeResidual == wrapped.solver.physicalRelativeResidual,
                "disabled residual changed");
        require(base.traceHash == wrapped.solver.traceHash, "disabled trace changed");
        if (!payload)
            require(base.eventHash == wrapped.solver.eventHash, "disabled events changed");
        require(spin::spinSolverScientificEventHash(base) ==
                spin::spinSolverScientificEventHash(wrapped.solver),
                "disabled scientific events changed (timings excluded)");
        if (payload) {
            auto changed = base;
            auto& details = *changed.events.at(1).evaluationPayload;
            details.evaluatorWallSeconds = 123.0;
            details.residualEvaluationSeconds = 456.0;
            details.derivativeEvaluationSeconds = 789.0;
            require(spin::spinSolverScientificEventHash(base) ==
                    spin::spinSolverScientificEventHash(changed),
                    "diagnostic timing contaminated scientific identity");
            details.solverRawMoments[0] += 0.25;
            require(spin::spinSolverScientificEventHash(base) !=
                    spin::spinSolverScientificEventHash(changed),
                    "changed moments retained scientific identity");
        }
        require(base.functionEvaluations == wrapped.solver.functionEvaluations &&
                baseCalls == wrappedCalls, "disabled true call count or order changed");
        require(base.restartCount == wrapped.solver.restartCount &&
                base.restartFunctionEvaluations == wrapped.solver.restartFunctionEvaluations,
                "disabled restart accounting changed");
        require(wrapped.attachment.predictorCalls == 0, "disabled hook called predictor");
        std::cout << "DISABLED_BASELINE trace=" << base.traceHash << " events="
                  << base.eventHash << " wrapped_events=" << wrapped.solver.eventHash
                  << " scientific_events=" << spin::spinSolverScientificEventHash(base)
                  << " true_calls=" << baseCalls.size() << '\n';
    }
    settings.equationProfile = spin::SpinEquationProfile::DualMonomialRawV1;
    bool rejected = false;
    try {
        spin::solveSpinMomentSystemLegacyTransactional(problem, settings, identityMap,
            {initialPoint()}, "kernel-fixture-v1", spin::Vector10::Zero(), {}, {});
    } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "transaction entry admitted DUAL even with disabled hook");
}
void requireSameAcceptedStates(const spin::SpinNewtonResult& base,
                               const spin::SpinNewtonResult& actual) {
    require(base.trace.size() == actual.trace.size(), "certified prefix lost or tentative state leaked");
    for (std::size_t i = 0; i < base.trace.size(); ++i) {
        require(base.trace[i].quadricHash == actual.trace[i].quadricHash,
                "rollback changed certified coefficients");
        require(base.trace[i].objective == actual.trace[i].objective &&
                base.trace[i].multiplierMu == actual.trace[i].multiplierMu,
                "tentative moments or multipliers contaminated certified state");
    }
}
void checkLedger(const spin::LegacyTransactionResult& r, int observedCalls) {
    require(r.attachment.trueCalls == static_cast<std::uint64_t>(observedCalls),
            "host ledger omitted a real call");
    require(r.solver.functionEvaluations == observedCalls, "solver budget omitted a real call");
    std::uint64_t calls = 0, trials = 0, sequence = 0;
    for (const auto& event : r.attachment.events) {
        require(event.sequence == ++sequence, "cost event was deleted or reordered");
        require(event.trueCalls >= calls && event.logicalTrials >= trials,
                "rollback refunded consumed cost or logical budget");
        require(std::isfinite(event.exclusiveSeconds) && event.exclusiveSeconds >= 0,
                "invalid host exclusive time");
        calls = event.trueCalls; trials = event.logicalTrials;
    }
    require(!r.attachment.scientificHash.empty(), "transaction scientific identity absent");
}
void singleStep() {
    const auto problem = spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    spin::SpinSolverSettings settings;
    settings.captureEvaluationPayload = true;
    spin::Vector10 fixed = spin::Vector10::Zero(); fixed(0) = 0.25;
    auto map = [&](const auto& q) {
        auto result = identityMap(q); result.rawMoments += fixed;
        result.variableRawMoments = q.values(); return result;
    };
    const auto base = spin::solveSpinMomentSystem(problem, settings, map, {initialPoint()});
    int calls = 0;
    std::function<spin::LegacyPredictionResponse()> expiredQuery;
    const auto solved = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
        [&](const auto& q) { ++calls; return map(q); }, {initialPoint()}, "kernel-fixed-v1", fixed,
        [&](const spin::LegacyPredictionRequest& r) {
            expiredQuery = r.queryTrueResponse;
            auto value = r.queryTrueResponse();
            require(value.variableMoments == r.q, "query exposed total rather than variable moments");
            return value;
        }, {});
    require(solved.solver.status == spin::SpinSolverStatus::Converged,
            "transparent single-step failed despite sufficient original budget");
    requireSameAcceptedStates(base, solved.solver);
    require(solved.attachment.committedBlocks > 0 && solved.attachment.rolledBackBlocks == 0,
            "valid transparent proposals were not truly certified");
    require(calls > base.functionEvaluations, "transparent real queries were hidden in predictor timing");
    require(solved.attachment.predictorCalls > 0 && solved.attachment.logicalTrials > 0,
            "transaction did not execute actual proposals");
    bool sawMoments = false;
    for (const auto& event : solved.attachment.events) {
        if (!event.tentativeMoments) continue;
        require((*event.tentativeMoments)(0) == event.coefficients(0) + 0.25,
                "fixed shell added twice or density repeated");
        require(!event.geometryCertified, "prediction self-certified geometry");
        sawMoments = true;
    }
    require(sawMoments, "no tentative material evidence");
    bool rejected = false;
    try { expiredQuery(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected && calls == solved.solver.functionEvaluations,
            "escaped predictor query outlived budget and checkpoint scope");
    checkLedger(solved, calls);
}
void rollback() {
    const auto problem = spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    spin::SpinSolverSettings settings;
    const auto base = spin::solveSpinMomentSystem(problem, settings, identityMap, {initialPoint()});
    for (const std::string fault : {"stale", "nonfinite", "exception", "certification"}) {
        int calls = 0, predictions = 0;
        bool failCertification = false;
        const auto solved = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
            [&](const auto& q) {
                ++calls;
                if (failCertification) { failCertification = false;
                    throw std::runtime_error("frozen true certification failure"); }
                return identityMap(q);
            }, {initialPoint()}, "kernel-fault-v1", spin::Vector10::Zero(),
            [&](const spin::LegacyPredictionRequest& r) {
                ++predictions;
                auto value = response(r);
                if (fault == "stale") value.stateIdentity += "-stale";
                else if (fault == "nonfinite") value.variableMoments(0) =
                    std::numeric_limits<double>::quiet_NaN();
                else if (fault == "exception") throw std::runtime_error("frozen predictor error");
                else failCertification = true;
                return value;
            }, {});
        require(solved.solver.status == spin::SpinSolverStatus::Converged,
                "true fallback did not recover certified baseline");
        requireSameAcceptedStates(base, solved.solver);
        require(solved.attachment.rolledBackBlocks == 2 && predictions == 2 &&
                solved.attachment.predictorDisabled,
                "two failed blocks did not disable prediction for this solve");
        require(solved.attachment.committedBlocks == 0, "fault granted tentative convergence");
        checkLedger(solved, calls);
    }
    bool rejected = false;
    int calls = 0;
    const auto nested = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
        [&](const auto& q) { ++calls; return identityMap(q); }, {initialPoint()},
        "kernel-reentry-v1", spin::Vector10::Zero(),
        [&](const spin::LegacyPredictionRequest& r) {
            try { spin::solveSpinMomentSystem(problem, settings, identityMap, {initialPoint()}); }
            catch (const std::invalid_argument&) { rejected = true; }
            return r.queryTrueResponse();
        }, {});
    require(rejected && nested.solver.status == spin::SpinSolverStatus::Converged,
            "prediction callback solver reentry was not rejected");
    checkLedger(nested, calls);
}
void budget() {
    const auto problem = spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    for (bool failCertificate : {false, true}) {
        spin::SpinSolverSettings settings;
        settings.maximumFunctionEvaluationsPerRestart = 2;
        int calls = 0;
        const auto solved = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
            [&](const auto& q) {
                ++calls;
                if (failCertificate && calls == 2) throw std::runtime_error("last budget failed");
                return identityMap(q);
            }, {initialPoint()}, "kernel-budget-v1", spin::Vector10::Zero(),
            [](const spin::LegacyPredictionRequest& r) { return response(r); }, {});
        require(calls == 2 && solved.solver.status == spin::SpinSolverStatus::ResourceLimitReached,
                "certification spent last query but fallback was free or tentative converged");
        require(solved.attachment.lastCertifiedCoefficients.has_value(), "lost last certified state");
        if (failCertificate) {
            require(solved.solver.trace.size() == 1 &&
                    *solved.attachment.lastCertifiedCoefficients == initialPoint().values(),
                    "failed last certification leaked candidate or partial trace");
            require(solved.attachment.rolledBackBlocks == 1, "budget failure omitted rollback");
        }
        checkLedger(solved, calls);
    }
    spin::SpinSolverSettings settings;
    settings.maximumFunctionEvaluationsPerRestart = 1;
    int calls = 0, predictions = 0;
    const auto solved = spin::solveSpinMomentSystemLegacyTransactional(problem, settings,
        [&](const auto& q) { ++calls; return identityMap(q); }, {initialPoint()},
        "kernel-initial-only-v1", spin::Vector10::Zero(),
        [&](const auto& r) { ++predictions; return response(r); }, {});
    require(calls == 1 && predictions == 0 && solved.solver.trace.size() == 1,
            "initial-only budget permitted an unbudgeted proposal");
    checkLedger(solved, calls);
}
}
int main(int argc, char** argv) {
    try {
        if (argc == 1 || std::string(argv[1]) == "--contract") contract();
        else if (std::string(argv[1]) == "--disabled") disabled();
        else if (std::string(argv[1]) == "--single-step") singleStep();
        else if (std::string(argv[1]) == "--rollback") rollback();
        else if (std::string(argv[1]) == "--budget") budget();
        else throw std::invalid_argument("unknown regression scenario");
        std::cout << "LEGACY_TRANSACTION_CONTRACT_PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "LEGACY_TRANSACTION_FAIL " << e.what() << '\n';
        return 1;
    }
}
