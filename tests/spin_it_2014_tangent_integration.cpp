#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SlqpAuditTrace.h"
#include <Eigen/Cholesky>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace spin2014;
static void require(bool c,const char* why){if(!c)throw std::runtime_error(why);}

// Independent finite enumeration of all 3^4 active faces of a convex box QP.
static double oracle(const Eigen::Matrix4d& q,const Eigen::Vector4d& target) {
    double best=INFINITY;
    for(int code=0;code<81;++code){
        int value=code;Eigen::Vector4d x=Eigen::Vector4d::Zero();std::vector<int> free;
        for(int i=0;i<4;++i){int face=value%3;value/=3;if(face==0)free.push_back(i);else x(i)=face==2?1.:0.;}
        Eigen::MatrixXd a(free.size(),free.size());Eigen::VectorXd rhs(free.size());
        auto b=q*target-q*x;
        for(int i=0;i<int(free.size());++i){rhs(i)=b(free[i]);for(int j=0;j<int(free.size());++j)a(i,j)=q(free[i],free[j]);}
        if(!free.empty()){Eigen::VectorXd solved=a.ldlt().solve(rhs);for(int i=0;i<int(free.size());++i)x(free[i])=solved(i);}
        if((x.array()<-1e-12).any()||(x.array()>1+1e-12).any())continue;
        best=std::min(best,.5*(x-target).dot(q*(x-target)));
    }return best;
}
int main(){try{
    Eigen::Matrix4d a;a<<1,2,-3,4,-2,4,1,3,3,-1,4,2,4,3,2,-1;
    const Eigen::Matrix4d q=a.transpose()*a+.25*Eigen::Matrix4d::Identity();
    const std::vector<Eigen::Vector4d> targets={{1.2,-.3,.8,1.4},{-.2,.8,1.1,.2},
        {.2,1.3,-.4,.6},{1.1,.2,1.5,-.3},{-.1,-.2,.7,1.3},{.9,1.2,-.6,1.4}};
    int attempts=0,accepted=0,unknownFallbacks=0,rejectedFallbacks=0;
    for(const auto& target:targets){
        EqualityBoxProblem p;p.variableCount=p.betaCount=4;
        p.evaluate=[q,target](const Eigen::VectorXd& x){AdaptiveFillEvaluation e;
            e.objective=.5*(x-target).dot(q*(x-target));e.gradient=q*(x-target);
            e.constraints.resize(0);e.constraintJacobian.resize(0,4);return e;};
        SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;s.maximumIterations=240;
        const Eigen::Vector4d initial(.8,.2,.8,.2);
        auto old=solveSlqpLbfgs(p,initial,s);auto oldAgain=solveSlqpLbfgs(p,initial,s);
        require(old.traceHash==oldAgain.traceHash,"legacy repeat mismatch");
        s.directionPolicy=SlqpDirectionPolicy::MetricTangentRescueV1;
        auto r=solveSlqpLbfgs(p,initial,s);
        require(r.status==SlqpStatus::Converged,"rescue convex box solve incomplete");
        require(r.finalKkt.valid&&r.finalKkt.full<=1e-7,"new profile did not meet KKT");
        require(std::abs(r.evaluation.objective-oracle(q,target))<1e-6,"independent box oracle mismatch");
        require(r.traceHash!=old.traceHash,"policy missing from trace hash");
        for(const auto& t:r.trace)if(t.metricTangent.attempted){
            ++attempts;accepted+=t.metricTangentAccepted;
            if(t.metricTangentAccepted)require(t.metricTangent.checked&&t.metricTangentSelected&&!t.usedFirstOrderFallback,"rescued acceptance bypassed gate");
            if(t.metricTangentAccepted) {
                auto limited=s;limited.metricTangent.maximumFaces=1;
                // Add an independent permanently lower-active coordinate so
                // the one-face budget is genuinely exceeded, not merely small.
                EqualityBoxProblem padded;padded.variableCount=padded.betaCount=5;
                padded.evaluate=[q,target](const Eigen::VectorXd& x){AdaptiveFillEvaluation e;
                    const Eigen::Vector4d delta=x.head<4>()-target;
                    e.objective=.5*delta.dot(q*delta)+x(4)+.5*x(4)*x(4);
                    e.gradient.resize(5);e.gradient.head<4>()=q*delta;e.gradient(4)=1+x(4);
                    e.constraints.resize(0);e.constraintJacobian.resize(0,5);return e;};
                Eigen::VectorXd paddedInitial(5);paddedInitial.head<4>()=initial;paddedInitial(4)=0;
                auto paddedSettings=s;paddedSettings.directionPolicy=SlqpDirectionPolicy::ProjectedGradientV1;
                auto paddedOld=solveSlqpLbfgs(padded,paddedInitial,paddedSettings);
                auto unknown=solveSlqpLbfgs(padded,paddedInitial,limited);
                for(const auto& u:unknown.trace)if(u.metricTangent.attempted&&!u.metricTangent.checked){
                    require(!u.metricTangentSelected&&!u.metricTangentAccepted,"UNKNOWN QP was admitted");
                    require(u.metricTangentFallbackAttempted&&u.usedFirstOrderFallback,"UNKNOWN lost old fallback");
                    ++unknownFallbacks;
                }
                require(unknown.status==paddedOld.status&&unknown.beta.isApprox(paddedOld.beta,0.),"UNKNOWN changed legacy terminal state");
                // Deliberately invalid callback ONLY along the witnessed rescue ray.
                // This is a failure-path fixture, not an optimization/paper case.
                EqualityBoxProblem guarded=p;
                const Eigen::VectorXd origin=t.state, ray=t.metricTangent.direction;
                guarded.evaluate=[q,target,origin,ray](const Eigen::VectorXd& x){
                    const Eigen::VectorXd offset=x-origin;
                    const double length=offset.norm(),alpha=offset.dot(ray)/ray.squaredNorm();
                    if(length>1e-14&&alpha>0&&(offset-alpha*ray).norm()<1e-8*length)
                        throw std::runtime_error("TEST_REJECT_RESCUE_RAY");
                    AdaptiveFillEvaluation e;e.objective=.5*(x-target).dot(q*(x-target));e.gradient=q*(x-target);
                    e.constraints.resize(0);e.constraintJacobian.resize(0,4);return e;
                };
                auto rejected=solveSlqpLbfgs(guarded,initial,s);
                for(const auto& u:rejected.trace)if(u.metricTangentSelected&&!u.metricTangentAccepted){
                    require(u.metricTangent.checked&&u.metricTangentFallbackAttempted&&u.usedFirstOrderFallback,
                        "rejected rescue did not return to accepted fallback");
                    require(u.metricTangentTrialCount==s.maximumLineSearchSteps,"rescue ray not fully tested");
                    for(int i=0;i<u.metricTangentTrialCount;++i)
                        require(!u.objectiveTrials[i].valid&&!u.objectiveTrials[i].accepted,"invalid rescue trial accepted");
                    require(u.historyResetReason.find("FIRST_ORDER_FALLBACK")!=std::string::npos,"fallback did not reset history");
                    ++rejectedFallbacks;
                }
                require(rejected.status==SlqpStatus::Converged,"fallback after rejected rescue incomplete");
            }
            std::ostringstream json;writeSlqpAuditJson(json,t);
            require(json.str().find("metric_tangent")!=std::string::npos,"missing tangent evidence JSON");
        }
        if(!r.trace.empty()){
            auto changed=r;changed.trace.front().metricTangent.attempted=true;
            changed.trace.front().metricTangent.status="TAMPERED";
            require(hashAuditedSlqp(changed,s)!=hashAuditedSlqp(r,s),"QP evidence omitted from digest");
        }
        for(std::size_t i=0;i<r.trace.size();++i)if(r.trace[i].metricTangentAccepted){
            const auto bound=[&](const auto& mutate){auto copy=r;mutate(copy.trace[i]);
                require(hashAuditedSlqp(copy,s)!=hashAuditedSlqp(r,s),"individual QP field omitted from digest");};
            bound([](auto& t){t.metricTangent.checked=false;});
            bound([](auto& t){t.metricTangent.direction(0)+=.125;});
            bound([](auto& t){t.metricTangent.multipliers(0)+=.125;});
            bound([](auto& t){t.metricTangent.workingSet.push_back(99);});
            bound([](auto& t){t.metricTangent.equalityResidual+=.125;});
            bound([](auto& t){t.metricTangent.stationarity+=.125;});
            bound([](auto& t){t.metricTangentSelected=false;});
            bound([](auto& t){t.metricTangentAccepted=false;});
            bound([](auto& t){++t.metricTangentTrialCount;});
            bound([](auto& t){t.metricTangentFallbackAttempted=true;});
            bound([](auto& t){t.metricTangent.trace.front().action="TAMPERED";});
            bound([](auto& t){t.metricTangent.trace.front().changedIndex+=1;});
            bound([](auto& t){t.metricTangent.trace.front().workingSet.push_back(99);});
            bound([](auto& t){t.metricTangent.trace.front().equalityResidual+=.125;});
            bound([](auto& t){t.metricTangent.trace.front().stationarity+=.125;});
        }
        std::cout<<"legacy="<<slqpStatusToken(old.status)<<" rescue="<<slqpStatusToken(r.status)
                 <<" old_records="<<old.trace.size()<<" new_records="<<r.trace.size()<<'\n';
    }
    std::cout<<"tangent_attempts="<<attempts<<" nonlinear_accepted="<<accepted
             <<" unknown_fallbacks="<<unknownFallbacks<<" rejected_fallbacks="<<rejectedFallbacks<<'\n';
    require(attempts>0&&accepted>0,"native integration never exercised actual rescue");
    require(unknownFallbacks>0&&rejectedFallbacks>0,"native integration never exercised rescue refusal");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
