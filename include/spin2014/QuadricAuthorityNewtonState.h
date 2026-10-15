#pragma once
#include "spin2014/QuadricAuthorityMomentProblem.h"
#include "spin/QuadricCoefficientAuthority.h"
#include "spin/QuadricAuthorityLinearSolve.h"
#include <optional>
#include <functional>

namespace spin2014 {
/** Fixed-state Eq.14 assembly, NOT convergence or material admission. The
 * caller supplies the actual same-state native coefficient Jacobian and total
 * density-weighted moments; a nonempty map identity binds the attachment.
 * Only full-rank multiplier projections are accepted, without regularization.
 */
struct QuadricAuthorityNewtonState {
    unsigned precisionBits=1024;
    std::string coefficientPolicy;
    std::string evaluationMapIdentity;
    std::string multiplierMethod;
    spin::QuadricCoefficientAuthority10::RationalVector geometryAuthority;
    spin::QuadricCoefficientAuthority10::RationalVector dualAuthority;
    QuadricAuthorityMomentEvaluation objective;
    AuthorityMomentMatrix10 nativeMomentJacobian;
    AuthorityMomentMatrix10 dualMomentJacobian;
    std::vector<std::string> multipliers;
    std::string multiplierMu;
    AuthorityMomentVector10 optimalityResidual;
    std::string rawEquationResidual;
    std::string multiplierProjectionResidual;
    std::vector<std::vector<std::string>> matrix;
    std::vector<std::string> rightHandSide;
    std::string scientificHash;
};
QuadricAuthorityNewtonState linearizeSpin2014AuthorityNewtonState(
    const QuadricObjectiveSettings& settings,
    const spin::QuadricCoefficientAuthority10& coefficientAuthority,
    const AuthorityMomentVector10& totalRawMoments,
    const AuthorityMomentMatrix10& nativeMomentJacobian,
    const std::string& evaluationMapIdentity);

/** Numerical screening only; does not authenticate caller state, prove
 * material legality or grant final convergence. The enclosing solver must
 * reassemble owned state and independently validate any stationary candidate.
 * Comparisons use canonical rational authority, never its binary64 view. */
std::string classifySpin2014AuthorityNumericalState(
    const QuadricAuthorityNewtonState&,bool initialState);

/** Protected original-equation solve followed by MPFR1024 dual tangent
 * projection and the existing min(norm(direction),1/2) initial arc. The state
 * is reassembled from its owned inputs to reject stale/modified attachments.
 * Neither a defined direction nor its geodesic update proves nonlinear
 * decrease, convergence, derivative applicability or material legality.
 */
struct QuadricAuthorityNewtonDirection {
    bool defined=false;
    spin::SpinNewtonAuthorityLinearSolveResult linearSolve;
    spin::QuadricCoefficientAuthority10::RationalVector dualDirection;
    std::string directionNorm,initialArcLength,stateIdentity,scientificHash,diagnostic;
};
QuadricAuthorityNewtonDirection solveSpin2014AuthorityNewtonDirection(
    const QuadricObjectiveSettings&,const QuadricAuthorityNewtonState&,
    const std::string& absoluteLinearResidualLimit);

/** Actual same-operation Boolean material and its authority Eq.14 state.
 * This is a numerical point, not a final material seal or convergence claim.
 */
struct QuadricAuthorityNewtonPoint {
    spin::QuadricAuthorityTotalMomentEvaluation material;
    QuadricAuthorityNewtonState state;
    /** Frozen effective evaluator recipe; never inferred from a display or
     * borrowed from a differently discretized base point. */
    std::string evaluatorRequestIdentity;
};
QuadricAuthorityNewtonPoint evaluateSpin2014AuthorityNewtonPoint(
    const spin::QuadricAuthorityMaterialContext&,const QuadricObjectiveSettings&,
    const spin::QuadricCoefficientAuthority10&,const spin::QuadricMomentEvaluatorSettings&);

struct QuadricAuthorityNewtonTrial {
    QuadricAuthorityNewtonPoint point;
    bool accepted=false,topologyChanged=false;
    std::string arcLength,baseStateIdentity,reason,scientificHash;
};
/** Replay a specified direction/arc against ACTUAL new material. Supports
 * frozen diagnostic directions as well as protected Newton directions.
 * It never substitutes a linear prediction for the residual-decrease check;
 * geometry/derivative failures propagate for the enclosing attempt ledger.
 */
QuadricAuthorityNewtonTrial evaluateSpin2014AuthorityNewtonTrial(
    const spin::QuadricAuthorityMaterialContext&,const QuadricObjectiveSettings&,
    const QuadricAuthorityNewtonPoint&,const spin::QuadricCoefficientAuthority10::RationalVector&,
    const std::string& arcLength,const spin::QuadricMomentEvaluatorSettings&);

struct QuadricAuthorityNewtonTrialRecord {
    std::string arcLength,reason,diagnostic,rawEquationResidual,trialIdentity;
    spin::QuadricCoefficientAuthority10::RationalVector geometryAuthority,dualAuthority;
    bool accepted=false;
};
/** Observation only. Emitted outside geometry exception handlers so observer
 * interruption propagates; no convergence or checkpoint semantics implied. */
struct QuadricAuthorityNewtonEvent {
    unsigned sequence=0,restart=0,evaluationIndex=0;
    std::string type,status,diagnostic,stateIdentity,materialIdentity,trialIdentity,rawResidual,arcLength;
    spin::QuadricCoefficientAuthority10::RationalVector geometryAuthority{},dualAuthority{};
    bool accepted=false;
};
using QuadricAuthorityEventSink=std::function<void(const QuadricAuthorityNewtonEvent&)>;
struct QuadricAuthorityNewtonAdvance {
    bool accepted=false;
    unsigned productionEvaluations=0;
    QuadricAuthorityNewtonDirection direction;
    std::vector<QuadricAuthorityNewtonTrialRecord> trials;
    /** Present only for an accepted update. Rejected numerical trial evidence
     * stays in trials and must not escape as a committable material point. */
    std::optional<QuadricAuthorityNewtonTrial> finalTrial;
    std::string status,scientificHash;
};
/** One protected Newton step with existing 1/2 cap and 1/2 backtracking,
 * freshly constructed actual trial material, and all failures retained.
 * Caller supplies the REMAINING per-attempt budget (1..60). No extra
 * budget, cap growth, backend switch or convergence claim is introduced.
 */
QuadricAuthorityNewtonAdvance advanceSpin2014AuthorityNewtonPoint(
    const spin::QuadricAuthorityMaterialContext&,const QuadricObjectiveSettings&,
    const QuadricAuthorityNewtonPoint&,const spin::QuadricMomentEvaluatorSettings&,
    const std::string& absoluteLinearResidualLimit,unsigned remainingEvaluationBudget,
    const QuadricAuthorityEventSink& eventSink={});

struct QuadricAuthorityNewtonAttempt {
    unsigned productionEvaluations=0,acceptedUpdates=0;
    spin::QuadricCoefficientAuthority10::RationalVector initialGeometryAuthority,initialDualAuthority;
    std::string initialStateIdentity,initialRawResidual,initialMu;
    /** Actual first construction, absent on initial geometry/derivative failure.
     * Moment/identity evidence is not an independently audited geometry seal. */
    std::optional<AuthorityMomentVector10> initialTotalMoments;
    std::string initialObjective,initialMaterialIdentity;
    std::vector<QuadricAuthorityNewtonTrialRecord> trials;
    std::vector<std::string> stepIdentities,acceptedStateIdentities;
    std::optional<QuadricAuthorityNewtonPoint> stationaryCandidate;
    std::string finalStateIdentity,finalRawResidual,status,diagnostic,scientificHash;
};
/** Reuse the registered platform-stable seed sequence, lift each represented
 * binary64 start exactly to canonical rational input, then normalize at
 * MPFR1024. No solution/checkpoint/old endpoint input is accepted. This is
 * initialization only, not a claim that a start is geometrically admissible. */
std::vector<spin::QuadricCoefficientAuthority10> makeSpin2014AuthorityInitialSequence(
    spin::SpinMode,const spin::SpinSolverSettings&);
/** Actual one-attempt loop. The initial Boolean counts toward the 1..60
 * budget. Positive mu is an INITIAL screen, not an intermediate constraint.
 * A stationary candidate is numerical only: no seal, independent check or
 * restart/answer-free initialization is implied by this function. */
QuadricAuthorityNewtonAttempt runSpin2014AuthorityNewtonAttempt(
    const spin::QuadricAuthorityMaterialContext&,const QuadricObjectiveSettings&,
    const spin::QuadricCoefficientAuthority10&,const spin::QuadricMomentEvaluatorSettings&,
    const std::string& absoluteLinearResidualLimit,unsigned evaluationBudget=60,
    const QuadricAuthorityEventSink& eventSink={});
struct QuadricAuthorityNewtonRestarts {
    std::vector<QuadricAuthorityNewtonAttempt> attempts;
    std::optional<unsigned> stationaryAttemptIndex;
    unsigned productionEvaluations=0,seed=0;
    std::string status,protocolScope,scientificHash;
    std::vector<QuadricAuthorityNewtonEvent> events;
};
/** Starts derive ONLY from the context request seed/mode/settings. All
 * executed attempts, including initial failures, are retained. A reduced
 * budget is explicitly diagnostic, never a frozen full-run substitute.
 * Numerical candidate status cannot grant a final material seal. */
QuadricAuthorityNewtonRestarts runSpin2014AuthorityNewtonRestarts(
    const spin::QuadricAuthorityMaterialContext&,const QuadricObjectiveSettings&,
    const spin::QuadricMomentEvaluatorSettings&,const std::string& absoluteLinearResidualLimit,
    unsigned evaluationBudgetPerAttempt=60,const QuadricAuthorityEventSink& eventSink={});
} // namespace spin2014
