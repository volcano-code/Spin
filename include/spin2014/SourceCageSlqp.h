#pragma once
#include "spin2014/SourceCagePaperFunctional.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include <array>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
// All defaults and accounting factors are IMPLEMENTATION_DEFAULT.
struct SourceCageSlqpBudget {
    // Cannot exceed INT_MAX/6: generic trace matrices use int coefficient indices.
    std::size_t maximumStateValues=100000;
    std::size_t maximumSourceEvaluations=10000; // Includes the final direct evaluation.
    std::size_t maximumTraceRows=1000;
    std::size_t maximumCallbackRecords=100000;
    std::size_t maximumAdapterBytes=512U*1024U*1024U;
};
struct SourceCageSlqpSettings {
    SlqpLbfgsSettings solver=[] {
        SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;return s;
    }();
    SourceCageSlqpBudget budget;
};
struct SourceCageSlqpResources {
    std::size_t stateValues=0,constraints=0,traceRows=0,callbackRecords=0;
    std::size_t tapeBytes=0,sourceOutputBytes=0,sourceWorkspaceBytes=0;
    std::size_t traceBytes=0,historyBytes=0,solverWorkspaceBytes=0;
    std::size_t explicitLpBytes=0,hashWorkspaceBytes=0,requiredBytes=0;
};
enum class SourceCageSlqpOutcome { Completed, ResourceLimit, InitialEvaluationFailure };
enum class SourceCageSlqpGeometryStatus { NotEvaluated };
enum class SourceCageSlqpFeasibility { Unknown };
enum class SourceCageSlqpEvaluationPurpose { SolverCallback, FinalAcceptedState };
enum class SourceCageSlqpEvaluationStatus { Unused, Evaluated, InvalidNumericalTrial, ResourceLimit };
struct SourceCageSlqpEvaluationRecord {
    std::size_t sequence=0;
    SourceCageSlqpEvaluationPurpose purpose=SourceCageSlqpEvaluationPurpose::SolverCallback;
    SourceCageSlqpEvaluationStatus status=SourceCageSlqpEvaluationStatus::Unused;
    bool sourceInvoked=false;
    std::array<char,65> stateSha256{},outputSha256{},reasonSha256{};
    std::array<char,256> reason{}; // Truncated display only; digest binds full message.
};
struct SourceCageSlqpResult {
    SourceCageSlqpOutcome outcome=SourceCageSlqpOutcome::Completed;
    SourceCageSlqpGeometryStatus geometryStatus=SourceCageSlqpGeometryStatus::NotEvaluated;
    SourceCageSlqpFeasibility feasibility=SourceCageSlqpFeasibility::Unknown;
    // False when the first actual source callback failed before the generic
    // engine returned anything. Then numerical is only default storage, not a
    // manufactured termination or KKT witness; the attempted-state tape remains.
    bool numericalResultAvailable=false;
    SlqpLbfgsResult numerical; // Underlying termination preserved whenever available.
    SourceCageSlqpSettings settings;
    SourceCageSlqpResources resources;
    Eigen::VectorXd initialState;
    // Only the first recordCount columns/records are events. Remaining capacity
    // stays preallocated and zero/Unused; no tape growth after exhaustion.
    Eigen::MatrixXd callbackStates;
    std::vector<SourceCageSlqpEvaluationRecord> records;
    std::size_t callbackCount=0,sourceEvaluationCount=0,recordCount=0;
    bool finalEvaluationAvailable=false;
    SourceCagePaperEvaluation finalEvaluation;
    std::string configSha256,callbackTapeSha256,scientificSha256;
};
// Checked conservative simultaneous known-buffer envelope, not a process RSS or
// opaque HiGHS allocator guarantee. Immutable upstream storage has separate caps.
SourceCageSlqpResources planSourceCageSlqpResources(const SourceCagePaperFunctional&,
    const SourceCageSlqpSettings& = {});
// No user geometry callback or acceptance boolean. Rebind by constructing a new
// immutable functional and a fresh solve. Numeric convergence grants first-order
// tolerances only, never geometry admission, feasibility, or sealed material.
SourceCageSlqpResult solveSourceCageSlqp(const SourceCagePaperFunctional&,
    const Eigen::VectorXd& initialState,const SourceCageSlqpSettings& = {});
}
#endif
