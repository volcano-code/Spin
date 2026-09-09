#include "spin2014/SlqpAuditTrace.h"
#include <cmath>
#include <iomanip>
namespace spin2014 {
namespace {
void number(std::ostream& o,double v) {if(std::isfinite(v))o<<v;else o<<"null";}
void vector(std::ostream& o,const Eigen::VectorXd& v) {o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';number(o,v(i));}o<<']';}
void trials(std::ostream& o,const std::vector<RestorationTrial>& list) {
    o<<'[';bool first=true;for(const auto& t:list) {
        if(!first)o<<',';first=false;o<<"{\"step\":";number(o,t.step);
        o<<",\"constraint_l1\":";number(o,t.constraintL1);
        o<<",\"valid\":"<<(t.valid?"true":"false")<<",\"accepted\":"<<(t.accepted?"true":"false")<<",\"reason\":"<<std::quoted(t.reason)<<'}';
    }o<<']';
}
}
void writeSlqpAuditControlsJson(std::ostream& o,const SlqpLbfgsSettings& s) {
    o<<"{\"provenance\":\"IMPLEMENTATION_DEFAULT\",\"restoration_trust_bound\":"<<s.audit.restorationTrustBound
     <<",\"lp_iteration_limit\":"<<s.audit.lpIterationLimit<<",\"lp_time_limit\":"<<s.audit.lpTimeLimit
     <<",\"lp_feasibility_tolerance\":"<<s.audit.lpFeasibilityTolerance
     <<",\"threads\":1,\"parallel\":\"off\",\"random_seed\":0,\"solver\":\"simplex\",\"presolve\":\"on\""
     <<",\"maximum_iterations\":"<<s.maximumIterations<<",\"maximum_line_search_steps\":"<<s.maximumLineSearchSteps
     <<",\"lbfgs_memory\":"<<s.lbfgsMemory<<",\"equality_tolerance\":"<<s.equalityTolerance
     <<",\"projected_gradient_tolerance\":"<<s.projectedGradientTolerance<<",\"active_set_tolerance\":"<<s.activeSetTolerance
     <<",\"merit_penalty\":"<<s.meritPenalty<<",\"line_search_reduction\":"<<s.lineSearchReduction<<",\"armijo\":"<<s.armijo<<'}';
}
void writeSlqpAuditJson(std::ostream& o,const SlqpIteration& t) {
    o<<",\"state\":";vector(o,t.state);o<<",\"gradient\":";vector(o,t.gradient);
    o<<",\"constraints\":";vector(o,t.constraints);o<<",\"constraint_jacobian\":[";
    for(int j=0;j<t.constraintJacobian.rows();++j) {if(j)o<<',';vector(o,t.constraintJacobian.row(j).transpose());}o<<']';
    const auto& k=t.kkt;
    o<<",\"kkt\":{\"valid\":"<<(k.valid?"true":"false")<<",\"status\":"<<std::quoted(k.status)
     <<",\"lambda\":";vector(o,k.lambda);o<<",\"alpha\":";vector(o,k.alpha);o<<",\"eta\":";vector(o,k.eta);
    o<<",\"equality_primal\":";number(o,k.equalityPrimal);o<<",\"bound_primal\":";number(o,k.boundPrimal);
    o<<",\"stationarity\":";number(o,k.stationarity);o<<",\"dual_sign\":";number(o,k.dualSign);
    o<<",\"complementarity\":";number(o,k.complementarity);o<<",\"full\":";number(o,k.full);o<<'}';
    o<<",\"active_membership\":[";for(std::size_t i=0;i<t.activeMembership.size();++i) {if(i)o<<',';o<<int(t.activeMembership[i]);}o<<']';
    o<<",\"membership_hash\":"<<std::quoted(t.membershipHash)<<",\"history_reset_reason\":"<<std::quoted(t.historyResetReason)
     <<",\"history_before\":"<<t.historyBefore<<",\"history_after\":"<<t.historyAfter
     <<",\"curvature_pair_accepted\":"<<(t.curvaturePairAccepted?"true":"false")
     <<",\"curvature_pair_rejected\":"<<(t.curvaturePairRejected?"true":"false")
     <<",\"used_first_order_fallback\":"<<(t.usedFirstOrderFallback?"true":"false")
     <<",\"linearized_constraint_residual\":";number(o,t.linearizedConstraintResidual);
    o<<",\"final_state\":"<<(t.finalState?"true":"false")<<",\"objective_trials\":";trials(o,t.objectiveTrials);
    const auto& r=t.restoration;
    o<<",\"restoration\":{\"attempted\":"<<(r.attempted?"true":"false")<<",\"accepted\":"<<(r.accepted?"true":"false")
     <<",\"lp_status\":"<<std::quoted(r.status)<<",\"reason\":"<<std::quoted(r.reason)<<",\"direction\":";vector(o,r.direction);
    o<<",\"predicted_residual\":";number(o,r.predictedResidual);o<<",\"trials\":";trials(o,r.trials);o<<'}';
}
}
