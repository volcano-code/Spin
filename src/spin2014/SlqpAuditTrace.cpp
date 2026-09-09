#include "spin2014/SlqpAuditTrace.h"
#include "manufacturing/SlicerValidation.h"
#include <iomanip>
#include <sstream>
namespace spin2014 {
std::string slqpMembershipHash(const std::vector<unsigned char>& m) {
    std::string b="spin2014-active-membership-v1";for(auto v:m)b.push_back(static_cast<char>(v));
    return manufacturing::computeContentSha256(b);
}
void captureSlqpState(SlqpIteration& t,const AdaptiveFillEvaluation& e,const Eigen::VectorXd& x) {
    t.state=x;t.gradient=e.gradient;t.constraints=e.constraints;t.constraintJacobian=e.constraintJacobian;
}
std::string hashAuditedSlqp(const SlqpLbfgsResult& r,const SlqpLbfgsSettings& s) {
    std::ostringstream o;o.imbue(std::locale::classic());o<<std::hexfloat;
    auto v=[&](const auto& a) {o<<'['<<a.rows()<<','<<a.cols()<<':';for(int i=0;i<a.size();++i)o<<a.data()[i]<<',';o<<']';};
    auto k=[&](const KktResidual& a) {o<<a.valid<<a.status<<a.equalityPrimal<<','<<a.boundPrimal<<','<<a.stationarity<<','<<a.dualSign<<','<<a.complementarity<<','<<a.full;v(a.lambda);v(a.alpha);v(a.eta);};
    auto trials=[&](const auto& list) {for(const auto& a:list)o<<a.step<<','<<a.constraintL1<<','<<a.valid<<a.accepted<<a.reason<<';';};
    o<<"AUDITED_KKT_RESTORATION_V1 IMPLEMENTATION_DEFAULT FIRST_ORDER_ONLY "
     <<s.maximumIterations<<','<<s.maximumLineSearchSteps<<','<<s.lbfgsMemory<<','
     <<s.equalityTolerance<<','<<s.projectedGradientTolerance<<','<<s.activeSetTolerance<<','
     <<s.meritPenalty<<','<<s.lineSearchReduction<<','<<s.armijo<<','
     <<s.audit.restorationTrustBound<<','<<s.audit.lpIterationLimit<<','<<s.audit.lpTimeLimit<<','<<s.audit.lpFeasibilityTolerance
     <<" HIGHS:threads=1,parallel=off,seed=0,solver=simplex,presolve=on;";
    for(const auto& t:r.trace) {
        o<<t.iteration<<','<<t.objective<<','<<t.equalityInfinityNorm<<','<<t.projectedGradientInfinityNorm<<','<<t.kktInfinityNorm<<','<<t.meritValue<<','<<t.acceptedStep<<','<<t.lineSearchEvaluations<<','<<t.usedFirstOrderFallback<<','<<t.phi<<',';
        v(t.state);v(t.gradient);v(t.constraints);v(t.constraintJacobian);k(t.kkt);
        o<<t.membershipHash<<t.historyResetReason<<t.historyBefore<<','<<t.historyAfter<<','<<t.curvaturePairAccepted<<t.curvaturePairRejected<<','<<t.linearizedConstraintResidual<<','<<t.finalState;
        o<<t.restoration.attempted<<t.restoration.accepted<<t.restoration.status<<t.restoration.reason<<t.restoration.predictedResidual;
        v(t.restoration.direction);trials(t.restoration.trials);trials(t.objectiveTrials);
    }
    o<<slqpStatusToken(r.status)<<r.diagnostic<<" bounded_prefix="<<r.boundedVariableCount;v(r.finalState);k(r.finalKkt);
    return manufacturing::computeContentSha256(o.str());
}
}
