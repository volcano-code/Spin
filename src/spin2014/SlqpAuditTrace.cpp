#include "spin2014/SlqpAuditTrace.h"
#include "spin2014/SlqpIterationStore.h"
#include "manufacturing/SlicerValidation.h"
#include "verification/Sha256Stream.h"
#include <iomanip>
#include <stdexcept>
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
    const SlqpIterationTraversal traverse=[&](const auto& visit) {
        forEachSlqpIteration(r,visit);
    };
    return hashAuditedSlqp(r,s,traverse);
}
std::string hashAuditedSlqp(const SlqpLbfgsResult& r,const SlqpLbfgsSettings& s,
    const SlqpIterationTraversal& traverse) {
    if(!traverse)throw std::invalid_argument("audited SLQP hash requires a trace traversal");
    verification::Sha256Stream hash;auto& o=hash.stream();o.imbue(std::locale::classic());o<<std::hexfloat;
    auto v=[&](const auto& a) {o<<'['<<a.rows()<<','<<a.cols()<<':';for(int i=0;i<a.size();++i)o<<a.data()[i]<<',';o<<']';};
    auto k=[&](const KktResidual& a) {o<<a.valid<<a.status<<a.equalityPrimal<<','<<a.boundPrimal<<','<<a.stationarity<<','<<a.dualSign<<','<<a.complementarity<<','<<a.full;v(a.lambda);v(a.alpha);v(a.eta);};
    auto trials=[&](const auto& list) {for(const auto& a:list)o<<a.step<<','<<a.constraintL1<<','<<a.valid<<a.accepted<<a.reason<<';';};
    o<<"AUDITED_KKT_RESTORATION_V1 IMPLEMENTATION_DEFAULT FIRST_ORDER_ONLY "
     <<s.maximumIterations<<','<<s.maximumLineSearchSteps<<','<<s.lbfgsMemory<<','
     <<s.equalityTolerance<<','<<s.projectedGradientTolerance<<','<<s.activeSetTolerance<<','
     <<s.meritPenalty<<','<<s.lineSearchReduction<<','<<s.armijo<<','
     <<s.audit.restorationTrustBound<<','<<s.audit.lpIterationLimit<<','<<s.audit.lpTimeLimit<<','<<s.audit.lpFeasibilityTolerance
     <<" HIGHS:threads=1,parallel=off,seed=0,solver=simplex,presolve=on;";
    if(s.directionPolicy==SlqpDirectionPolicy::MetricTangentRescueV1)
        o<<" METRIC_TANGENT_RESCUE_V1 "<<s.metricTangent.maximumFaces<<','
         <<s.metricTangent.maximumIterations<<','<<s.metricTangent.residualTolerance<<';';
    if(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxQpV1 || s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2)
        o<<(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2?" FINITE_BOX_PRIMAL_V2 ":" FINITE_BOX_QP_V1 ")<<s.metricTangent.maximumFaces<<','
         <<s.metricTangent.maximumIterations<<','<<s.metricTangent.residualTolerance<<';';
    traverse([&](const SlqpIteration& t) {
        o<<t.iteration<<','<<t.objective<<','<<t.equalityInfinityNorm<<','<<t.projectedGradientInfinityNorm<<','<<t.kktInfinityNorm<<','<<t.meritValue<<','<<t.acceptedStep<<','<<t.lineSearchEvaluations<<','<<t.usedFirstOrderFallback<<','<<t.phi<<',';
        v(t.state);v(t.gradient);v(t.constraints);v(t.constraintJacobian);k(t.kkt);
        o<<t.membershipHash<<t.historyResetReason<<t.historyBefore<<','<<t.historyAfter<<','<<t.curvaturePairAccepted<<t.curvaturePairRejected<<','<<t.linearizedConstraintResidual<<','<<t.finalState;
        o<<t.restoration.attempted<<t.restoration.accepted<<t.restoration.status<<t.restoration.reason<<t.restoration.predictedResidual;
        v(t.restoration.direction);trials(t.restoration.trials);trials(t.objectiveTrials);
        if(t.metricTangent.attempted) {
            const auto& q=t.metricTangent;
            o<<(q.status.rfind("BOX_QP_",0)==0 ? " finite_box " : " metric_tangent ")<<q.attempted<<q.checked<<std::quoted(q.status)
             <<t.metricTangentSelected<<t.metricTangentAccepted<<','<<t.metricTangentTrialCount
             <<','<<t.metricTangentFallbackAttempted<<','<<q.equalityResidual<<','<<q.stationarity;
            v(q.direction);v(q.multipliers);
            for(auto i:q.workingSet)o<<i<<',';o<<';';
            for(const auto& w:q.trace) {
                o<<std::quoted(w.action)<<w.changedIndex<<','<<w.equalityResidual<<','<<w.stationarity<<'[';
                for(auto i:w.workingSet)o<<i<<',';o<<']';
            }
        }
    });
    o<<slqpStatusToken(r.status)<<r.diagnostic<<" bounded_prefix="<<r.boundedVariableCount;v(r.finalState);k(r.finalKkt);
    return hash.digest();
}
}
