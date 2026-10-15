#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SlqpAuditTrace.h"
#include "spin2014/SlqpIterationStore.h"
#include <Eigen/Cholesky>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace spin2014;
static void require(bool c,const char* m){if(!c)throw std::runtime_error(m);}

// Independent enumeration of every face of the four-dimensional final QP.
static double oracle(const Eigen::Matrix4d& q,const Eigen::Vector4d& target) {
    double best=INFINITY;
    for(int code=0;code<81;++code) {
        int v=code;Eigen::Vector4d x=Eigen::Vector4d::Zero();std::vector<int> free;
        for(int i=0;i<4;++i){int k=v%3;v/=3;if(!k)free.push_back(i);else x[i]=k==2?1.:0.;}
        Eigen::MatrixXd a(free.size(),free.size());Eigen::VectorXd b(free.size());
        const Eigen::Vector4d rhs=q*(target-x);
        for(int i=0;i<int(free.size());++i){b[i]=rhs[free[i]];for(int j=0;j<int(free.size());++j)a(i,j)=q(free[i],free[j]);}
        if(!free.empty()){Eigen::VectorXd z=a.ldlt().solve(b);for(int i=0;i<int(free.size());++i)x[free[i]]=z[i];}
        if((x.array()<-1e-12).any()||(x.array()>1+1e-12).any())continue;
        best=std::min(best,.5*(x-target).dot(q*(x-target)));
    }return best;
}
int main(){try{
    Eigen::Matrix4d a;a<<1,2,-3,4,-2,4,1,3,3,-1,4,2,4,3,2,-1;
    const Eigen::Matrix4d q=a.transpose()*a+.25*Eigen::Matrix4d::Identity();
    int accepted=0;
    for(const Eigen::Vector4d target:std::vector<Eigen::Vector4d>{{1.2,-.3,.8,1.4},{-.2,.8,1.1,.2},
            {.2,1.3,-.4,.6},{1.1,.2,1.5,-.3},{-.1,-.2,.7,1.3},{.9,1.2,-.6,1.4}}) {
        EqualityBoxProblem p;p.variableCount=p.betaCount=4;
        p.evaluate=[q,target](const Eigen::VectorXd& x){AdaptiveFillEvaluation e;
            e.objective=.5*(x-target).dot(q*(x-target));e.gradient=q*(x-target);
            e.constraints.resize(0);e.constraintJacobian.resize(0,4);return e;};
        SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;
        s.directionPolicy=SlqpDirectionPolicy::FiniteBoxQpV1;s.maximumIterations=240;
#ifdef TEST_PRIMAL
        s.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;
#endif
        const Eigen::Vector4d initial(.8,.2,.8,.2);
        auto r=solveSlqpLbfgs(p,initial,s);
        require(r.status==SlqpStatus::Converged,"full box integration incomplete");
        require(r.finalKkt.valid&&r.finalKkt.full<=1e-7,"nonlinear KKT gate failed");
        require(std::abs(r.evaluation.objective-oracle(q,target))<1e-6,"enumerated oracle disagrees");
        require(r.traceHash==solveSlqpLbfgs(p,initial,s).traceHash,"repeat trace differs");
        for(const auto& t:r.trace)if(t.metricTangent.attempted) {
            require(t.metricTangent.checked,"unknown QP used");
            if(t.metricTangentAccepted){++accepted;require(t.metricTangentSelected&&!t.usedFirstOrderFallback,"false accepted flag");}
            std::ostringstream json;writeSlqpAuditJson(json,t);
            require(json.str().find("finite_box")!=std::string::npos&&json.str().find("metric_tangent")==std::string::npos,"wrong JSON discriminator");
        }
        auto changed=r;changed.trace.front().metricTangent.direction[0]+=.125;
        require(hashAuditedSlqp(changed,s)!=r.traceHash,"direction not hash bound");
        auto limited=s;limited.metricTangent.maximumIterations=1;
        auto stopped=solveSlqpLbfgs(p,initial,limited);
        require(stopped.status==SlqpStatus::NumericalUnknown,"QP limit did not stop UNKNOWN");
        require(stopped.finalState.isApprox(initial,0.),"unverified QP changed state");
    }
    require(accepted>6,"primary QP never accepted");
    std::cout<<"finite box integration passed; accepted "<<accepted<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
