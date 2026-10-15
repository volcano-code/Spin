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
     <<",\"merit_penalty\":"<<s.meritPenalty<<",\"line_search_reduction\":"<<s.lineSearchReduction<<",\"armijo\":"<<s.armijo;
    if(s.directionPolicy==SlqpDirectionPolicy::MetricTangentRescueV1)
        o<<",\"direction_policy\":\"METRIC_TANGENT_RESCUE_V1\",\"metric_tangent_maximum_faces\":"<<s.metricTangent.maximumFaces
         <<",\"metric_tangent_maximum_iterations\":"<<s.metricTangent.maximumIterations
         <<",\"metric_tangent_residual_tolerance\":"<<s.metricTangent.residualTolerance;
    if(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxQpV1 || s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2) {
        o<<",\"direction_policy\":"<<std::quoted(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2?"FINITE_BOX_PRIMAL_V2":"FINITE_BOX_QP_V1")
         <<",\"finite_box_maximum_faces\":";
        if(s.directionPolicy==SlqpDirectionPolicy::FiniteBoxPrimalV2 && s.metricTangent.maximumFaces==-1)
            o<<"\"BOUNDED_VARIABLES\"";
        else o<<s.metricTangent.maximumFaces;
        o<<",\"finite_box_maximum_iterations\":"<<s.metricTangent.maximumIterations
         <<",\"finite_box_residual_tolerance\":"<<s.metricTangent.residualTolerance;
    }
    o<<'}';
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
    if(t.metricTangent.attempted) {
        const auto& q=t.metricTangent;
        auto indices=[&](const auto& list){o<<'[';bool first=true;for(auto i:list){if(!first)o<<',';first=false;o<<i;}o<<']';};
        o<< (q.status.rfind("BOX_QP_",0)==0 ? ",\"finite_box\":" : ",\"metric_tangent\":")
         <<"{\"attempted\":true,\"checked\":"<<(q.checked?"true":"false")
         <<",\"status\":"<<std::quoted(q.status)<<",\"selected\":"<<(t.metricTangentSelected?"true":"false")
         <<",\"accepted\":"<<(t.metricTangentAccepted?"true":"false")
         <<",\"trial_count\":"<<t.metricTangentTrialCount
         <<",\"fallback_attempted\":"<<(t.metricTangentFallbackAttempted?"true":"false")<<",\"direction\":";vector(o,q.direction);
        o<<",\"multipliers\":";vector(o,q.multipliers);o<<",\"working_set\":";indices(q.workingSet);
        o<<",\"equality_residual\":";number(o,q.equalityResidual);o<<",\"stationarity\":";number(o,q.stationarity);
        o<<",\"trace\":[";bool first=true;
        for(const auto& w:q.trace){if(!first)o<<',';first=false;o<<"{\"working_set\":";indices(w.workingSet);
            o<<",\"action\":"<<std::quoted(w.action)<<",\"changed_index\":"<<w.changedIndex<<",\"equality_residual\":";
            number(o,w.equalityResidual);o<<",\"stationarity\":";number(o,w.stationarity);o<<'}';}
        o<<"]}";
    }
}
}
