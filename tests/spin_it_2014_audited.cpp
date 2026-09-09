#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SlqpKkt.h"
#include "spin2014/SlqpHistory.h"
#include "spin2014/AdaptiveVoxelSpinSolver.h"
#include "spin2014/DualMaterialSpinSolver.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <iomanip>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> static void rejects(F f) {bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}require(caught,"invalid input accepted");}
using namespace spin2014;
static SlqpLbfgsSettings audited() {SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;return s;}
static void independentlyCheck(const SlqpLbfgsResult& result, int betaCount) {
    for(const auto& t:result.trace) {
        if(!t.kkt.valid)continue;
        double equality=0,bound=0,stationarity=0,sign=0,comp=0;
        for(int j=0;j<t.constraints.size();++j)equality=std::max(equality,std::abs(t.constraints(j)));
        for(int i=0;i<t.state.size();++i) {
            double r=t.gradient(i)-t.kkt.alpha(i)+t.kkt.eta(i);
            for(int j=0;j<t.constraints.size();++j)r+=t.constraintJacobian(j,i)*t.kkt.lambda(j);
            stationarity=std::max(stationarity,std::abs(r));
            sign=std::max({sign,-t.kkt.alpha(i),-t.kkt.eta(i)});
            if(i<betaCount) {
                bound=std::max({bound,-t.state(i),t.state(i)-1.});
                comp=std::max({comp,std::abs(t.kkt.alpha(i)*t.state(i)),std::abs(t.kkt.eta(i)*(1.-t.state(i)))});
            } else require(t.kkt.alpha(i)==0 && t.kkt.eta(i)==0,"phi bound dual");
        }
        require(std::abs(equality-t.kkt.equalityPrimal)<1e-12 && std::abs(bound-t.kkt.boundPrimal)<1e-12 &&
            std::abs(stationarity-t.kkt.stationarity)<1e-10 && std::abs(sign-t.kkt.dualSign)<1e-12 &&
            std::abs(comp-t.kkt.complementarity)<1e-12,"independent literal-loop KKT recomputation mismatch");
    }
}
static void suppliedWitnesses() {
    AdaptiveFillEvaluation e;e.objective=0;e.gradient=Eigen::VectorXd::Zero(1);e.constraints.resize(0);e.constraintJacobian.resize(0,1);
    Eigen::VectorXd x=Eigen::VectorXd::Zero(1),zero=x,none(0);
    auto check=[&](double g,double state,double alpha,double eta) {e.gradient(0)=g;x(0)=state;
        return evaluateKkt(e,x,1,none,Eigen::VectorXd::Constant(1,alpha),Eigen::VectorXd::Constant(1,eta));};
    require(check(0,.5,0,0).full==0,"unconstrained stationary");
    require(check(1,0,1,0).full==0,"lower bound sign");
    require(check(-1,1,0,1).full==0,"upper bound sign");
    require(check(0,-.2,0,0).boundPrimal==.2,"outside beta");
    require(check(-1,0,-1,0).dualSign==1,"negative alpha");
    require(check(1,1,0,-1).dualSign==1,"negative eta");
    require(check(1,.5,1,0).complementarity==.5,"nonactive multiplier");
    e.gradient(0)=1;x(0)=0;
    require(!evaluateKkt(e,x,0,none,Eigen::VectorXd::Ones(1),zero).valid,"angle multiplier rejected");
    require(!evaluateKkt(e,x,1,none,none,zero).valid,"bad multiplier dimensions");
    e.constraintJacobian.resize(0,2);require(!evaluateKkt(e,x,1,none,zero,zero).valid,"bad Jacobian dimensions");e.constraintJacobian.resize(0,1);
    for(double bad:{NAN,INFINITY,-INFINITY}) {e.objective=bad;require(!evaluateKkt(e,x,1,none,zero,zero).valid,"nonfinite objective");}
    e.objective=0;e.constraints=Eigen::VectorXd::Zero(1);e.constraintJacobian=Eigen::MatrixXd::Constant(1,1,1e308);
    auto overflow=evaluateKkt(e,x,1,Eigen::VectorXd::Constant(1,1e308),zero,zero);
    require(!overflow.valid&&!std::isfinite(overflow.full),"overflow cannot return zero KKT");
    e.constraintJacobian.setOnes();x(0)=.75;
    require(evaluateKkt(e,x,1,Eigen::VectorXd::Constant(1,-1),zero,zero).full==0,"interior equality witness");
    // Near-face multiplier must pay complementarity in the minimax objective.
    e.constraints.resize(0);e.constraintJacobian.resize(0,1);e.gradient(0)=1e6;x(0)=1e-5;
    auto near=computeKktWitness(e,x,1,1e-4,{});
    require(near.valid&&near.complementarity==10.&&near.full>=10.,"near-bound multiplier hid complementarity");
}
static EqualityBoxProblem linearProblem(int variant) {
    EqualityBoxProblem p;p.variableCount=p.betaCount=1;
    p.evaluate=[variant](const Eigen::VectorXd& x) {
        AdaptiveFillEvaluation e;e.objective=x(0);e.gradient=Eigen::VectorXd::Ones(1);
        e.constraints=Eigen::VectorXd::Constant(variant==0?1:2,x(0)-.75);
        e.constraintJacobian=Eigen::MatrixXd::Ones(e.constraints.size(),1);
        if(variant==1){e.constraints(1)=2*x(0)-1.5;e.constraintJacobian(1,0)=2;}
        if(variant==2)e.constraints(1)=x(0)-.25;return e;
    };return p;
}
static void analyticProblems() {
    // A subnormal positive curvature can have an infinite reciprocal. A
    // rejected update must preserve the already-full metric, including its
    // oldest entry, rather than evicting it before validation completes.
    LbfgsHistory fullHistory;fullHistory.capacity=1;
    require(fullHistory.update(Eigen::VectorXd::Constant(1,2.),Eigen::VectorXd::Constant(1,4.)),
        "valid full-history setup failed");
    const Eigen::VectorXd tiny=Eigen::VectorXd::Constant(1,1e-160);
    const double tinyCurvature=tiny.dot(tiny);
    require(std::isfinite(tinyCurvature) && tinyCurvature>0. && !std::isfinite(1./tinyCurvature),
        "fixed subnormal fixture must exercise reciprocal overflow, not zero-curvature rejection");
    const bool tinyAccepted=fullHistory.update(tiny,tiny);
    const auto retainedMetric=fullHistory.applyInverse(Eigen::VectorXd::Constant(1,2.));
    std::cout<<"subnormal_curvature="<<tinyCurvature<<" accepted="<<tinyAccepted
        <<" retained_inverse="<<fullHistory.inverseCurvature.front()
        <<" retained_metric="<<retainedMetric(0)<<'\n';
    require(!tinyAccepted,"nonfinite reciprocal curvature accepted");
    require(fullHistory.s.size()==1 && fullHistory.y.size()==1 && fullHistory.inverseCurvature.size()==1 &&
        fullHistory.s.front()(0)==2. && fullHistory.y.front()(0)==4. && fullHistory.inverseCurvature.front()==.125 &&
        retainedMetric.allFinite() && retainedMetric(0)==1.,
        "rejected reciprocal overflow mutated full history or its inverse metric");
    LbfgsHistory history;history.capacity=3;
    require(history.update(Eigen::Vector2d(1,0),Eigen::Vector2d(2,0)),"positive curvature rejected");
    require(!history.update(Eigen::Vector2d(1,0),Eigen::Vector2d(-2,0)),"negative curvature accepted");
    require(!history.update(Eigen::Vector2d::Zero(),Eigen::Vector2d::Zero()),"zero curvature accepted");
    require(history.resetForMembership({1,0},{0,1}) && history.s.empty(),"same-count different membership did not reset real history");
    auto s=audited();
    const auto duplicate=solveSlqpLbfgs(linearProblem(1),Eigen::VectorXd::Zero(1),s);
    require(duplicate.status==SlqpStatus::Converged&&duplicate.finalKkt.full<=1e-7&&std::abs(duplicate.beta(0)-.75)<1e-8,"duplicate equality rank failure");
    require(std::abs(duplicate.finalKkt.lambda(0)+2*duplicate.finalKkt.lambda(1)+1)<1e-7,"rank-deficient witness");
    independentlyCheck(duplicate,1);
    auto inconsistent=solveSlqpLbfgs(linearProblem(2),Eigen::VectorXd::Zero(1),s);
    require(inconsistent.status!=SlqpStatus::Converged&&inconsistent.finalKkt.full>=.25,"inconsistent pair certified");
    require(std::abs(inconsistent.evaluation.constraints.lpNorm<1>()-.5)<1e-10,"fixed inconsistent min L1");
    auto singular=linearProblem(0);singular.evaluate=[](const Eigen::VectorXd& x) {
        auto e=linearProblem(0).evaluate(x);e.constraints(0)=x(0)*x(0)-.25;e.constraintJacobian(0,0)=2*x(0);return e;};
    const auto stalled=solveSlqpLbfgs(singular,Eigen::VectorXd::Zero(1),s);
    require(stalled.status==SlqpStatus::RestorationStalled && stalled.beta(0)==0 && !stalled.trace.front().restoration.accepted,"singular linearization not stalled");
    EqualityBoxProblem quadratic;quadratic.variableCount=quadratic.betaCount=2;
    quadratic.evaluate=[](const Eigen::VectorXd& x) {AdaptiveFillEvaluation e;Eigen::Vector2d optimum(.8,-.2);
        e.gradient=x-optimum;e.objective=.5*e.gradient.squaredNorm();e.constraints.resize(0);e.constraintJacobian.resize(0,2);return e;};
    Eigen::Vector2d initial(0,1);auto q=solveSlqpLbfgs(quadratic,initial,s);
    require(q.status==SlqpStatus::Converged && (q.beta-Eigen::Vector2d(.8,0)).norm()<1e-7,"fixed box quadratic bound release");
    require(std::abs(q.finalKkt.alpha(1)-.2)<1e-7,"final active multiplier");
    bool reset=false;for(const auto& t:q.trace)reset|=t.historyResetReason.find("ACTIVE_MEMBERSHIP_CHANGED")!=std::string::npos;
    require(reset,"real active transition history reset absent");independentlyCheck(q,2);
    auto oneStep=audited();oneStep.maximumIterations=1;
    const auto qLimit=solveSlqpLbfgs(quadratic,initial,oneStep);
    require(qLimit.status==SlqpStatus::MaximumIterations && qLimit.trace.back().activeMembership==std::vector<unsigned char>({0,1}) &&
        qLimit.trace.back().activeVariables==1 && qLimit.trace.back().historyBefore==1 &&
        qLimit.trace.back().historyAfter==0 && qLimit.trace.back().lbfgsHistorySize==0,
        "final snapshot active count/history must match actual accepted state");
    quadratic.betaCount=1;quadratic.evaluate=[](const Eigen::VectorXd& x) {AdaptiveFillEvaluation e;e.gradient=x-Eigen::Vector2d(.4,3.);
        e.objective=.5*e.gradient.squaredNorm();e.constraints.resize(0);e.constraintJacobian.resize(0,2);return e;};
    auto angle=solveSlqpLbfgs(quadratic,Eigen::Vector2d(.2,0),s);
    require(angle.status==SlqpStatus::Converged && (angle.beta-Eigen::Vector2d(.4,3)).norm()<1e-7,"angle clamped");independentlyCheck(angle,1);
    s.maximumIterations=1;auto last=solveSlqpLbfgs(linearProblem(0),Eigen::VectorXd::Zero(1),s);
    require(last.status==SlqpStatus::MaximumIterations && last.trace.back().finalState && last.trace.back().state(0)==.75 && last.finalKkt.full<1e-7,"max iteration final accepted state audit");
    require(last.trace.back().meritValue==.75,"final accepted-state merit must be evaluated, not default zero");
    auto bad=linearProblem(0);bad.evaluate=[](const Eigen::VectorXd& x) {auto e=linearProblem(0).evaluate(x);if(x(0)>0)e.objective=NAN;return e;};
    auto invalid=solveSlqpLbfgs(bad,Eigen::VectorXd::Zero(1),audited());
    require(invalid.status!=SlqpStatus::Converged && invalid.trace.front().restoration.trials.size()==24,"invalid restoration trials lost");
    for(const auto& t:invalid.trace.front().restoration.trials)require(!t.valid&&!t.accepted&&!t.reason.empty(),"invalid trial accepted or failure reason missing");
    s=audited();s.audit.lpTimeLimit=NAN;rejects([&]{solveSlqpLbfgs(linearProblem(0),Eigen::VectorXd::Zero(1),s);});
    s=audited();s.algorithm=static_cast<SlqpAlgorithm>(9);rejects([&]{solveSlqpLbfgs(linearProblem(0),Eigen::VectorXd::Zero(1),s);});
    auto curved=linearProblem(0);curved.evaluate=[](const Eigen::VectorXd& x) {
        auto e=linearProblem(0).evaluate(x);e.constraints(0)=.25-x(0)+1e30*x(0)*x(0);e.constraintJacobian(0,0)=-1+2e30*x(0);return e;};
    auto noDecrease=solveSlqpLbfgs(curved,Eigen::VectorXd::Zero(1),audited());
    require(noDecrease.status==SlqpStatus::RestorationStalled && noDecrease.trace.front().restoration.trials.size()==24,"finite restoration line search failure not retained");
    for(auto t:noDecrease.trace.front().restoration.trials)require(t.valid&&!t.accepted&&t.constraintL1>=.25,"finite nondecrease accepted");
    AdaptiveFillEvaluation cap;cap.objective=0;cap.gradient.resize(12);cap.constraints=Eigen::VectorXd::Zero(5);cap.constraintJacobian.resize(5,12);
    for(int i=0;i<12;++i){cap.gradient(i)=((i*17+3)%19)-9.;for(int j=0;j<5;++j)cap.constraintJacobian(j,i)=((i*7+j*11+i*j*3)%23)-11.;}
    SlqpAuditSettings limit;limit.lpIterationLimit=1;
    auto capped=computeKktWitness(cap,Eigen::VectorXd::Constant(12,.5),12,1e-10,limit);
    require(!capped.valid && capped.status.find("Iteration limit")!=std::string::npos,"actual LP cap not reported UNKNOWN");
    std::cout<<"duplicate="<<slqpStatusToken(duplicate.status)<<" inconsistent_L1="<<inconsistent.evaluation.constraints.lpNorm<1>()
        <<" singular="<<slqpStatusToken(stalled.status)<<" finite_failed_trials="<<noDecrease.trace.front().restoration.trials.size()
        <<" invalid_trials="<<invalid.trace.front().restoration.trials.size()<<" cap="<<capped.status<<'\n';
    auto changing=linearProblem(0);changing.evaluate=[](const Eigen::VectorXd& x) {auto e=linearProblem(0).evaluate(x);if(x(0)>0){e.constraints.resize(0);e.constraintJacobian.resize(0,1);}return e;};
    auto changedDimensions=solveSlqpLbfgs(changing,Eigen::VectorXd::Zero(1),audited());
    require(changedDimensions.status!=SlqpStatus::Converged && !changedDimensions.trace.front().restoration.trials.front().valid,"callback equality dimension changed during solve");
    EqualityBoxProblem constant;constant.variableCount=2;constant.betaCount=0;
    constant.evaluate=[](const Eigen::VectorXd&) {AdaptiveFillEvaluation e;e.gradient=Eigen::VectorXd::Zero(2);e.constraints.resize(0);e.constraintJacobian.resize(0,2);return e;};
    const auto freeHash=solveSlqpLbfgs(constant,Eigen::VectorXd::Constant(2,.5),audited()).traceHash;
    constant.betaCount=1;
    require(solveSlqpLbfgs(constant,Eigen::VectorXd::Constant(2,.5),audited()).traceHash!=freeHash,"bounded prefix policy omitted from trace identity");
}
static std::shared_ptr<const voxel::VoxelDomain> givensFixture() {
    voxel::GridSpec grid;grid.dimensions=Eigen::Vector3i(4,4,4);grid.minimum=Eigen::Vector3d(-1.5,-3,-1.5);grid.maximum=Eigen::Vector3d(1.5,3,1.5);
    std::vector<unsigned char> material(64,1),carveable(64,0);for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x)carveable[x+4*(y+4*z)]=1;
    return voxel::VoxelDomain::create(grid,material,{},carveable);
}
static void givensIntegration() {
    AdaptiveVoxelDomain domain(givensFixture(),2,2);
    AdaptiveFillFunctional f(domain,SpinIt2014MomentProblem::fromObjective(Spin2014Objective::TopInertiaOnly),.03);
    auto s=audited();s.maximumIterations=180;s.equalityTolerance=1e-9;
    Eigen::VectorXd x=Eigen::VectorXd::Constant(domain.leaves().size(),.2);
    const auto a=solveSlqpLbfgsPaper(f,x,1.7,s),b=solveSlqpLbfgsPaper(f,x,1.7,s);
    std::cout<<"givens="<<slqpStatusToken(a.status)<<" kkt="<<a.finalKkt.full<<" phi="<<a.finalPhi<<'\n';
    require(a.status==SlqpStatus::Converged && a.finalKkt.full<=1e-7 && a.finalPhi>1,"audited real Givens did not converge");
    require(a.traceHash==b.traceHash && a.finalState.isApprox(b.finalState,0),"audited replay mismatch");independentlyCheck(a,x.size());
    rejects([&]{solveSlqpLbfgs(f,x,s);});
    AdaptiveVoxelSpinRequest r;r.domain=givensFixture();r.objective=Spin2014Objective::TopInertiaOnly;r.initialLevel=2;r.maximumLevel=2;r.initialFill=.2;
    r.slqp=s;rejects([&]{validateAdaptiveVoxelSpinRequest(r);});
    r.formulation=Spin2014Formulation::PaperGivens;r.initialPhi=1.7;
    auto adaptive=solveAdaptiveVoxelSpin(r);require(!adaptive.binaryFill.empty(),"audited binary validation not reached");
    const auto h=computeAdaptiveVoxelSpinRequestHash(r);r.slqp.audit.restorationTrustBound=.5;require(computeAdaptiveVoxelSpinRequestHash(r)!=h,"audit controls omitted from identity");r.slqp=s;
    r.binaryConstraintTolerance=1e-14;require(solveAdaptiveVoxelSpin(r).status==AdaptiveSpinStatus::BinaryValidationFailed,"binary failure promoted");
    r.binaryConstraintTolerance=1e-6;r.initialLevel=1;auto refined=solveAdaptiveVoxelSpin(r);require(refined.adaptiveTrace.size()>1,"audited adaptive carry missing");
    for(std::size_t i=1;i<refined.adaptiveTrace.size();++i)require(refined.adaptiveTrace[i].initialPhi==refined.adaptiveTrace[i-1].finalPhi,"phi not carried");
}
int main() {try {
    std::cout<<std::setprecision(17);
    spin2014::EqualityBoxProblem problem;
    problem.variableCount = problem.betaCount = 1;
    problem.evaluate = [](const Eigen::VectorXd& x) {
        spin2014::AdaptiveFillEvaluation e;
        e.objective=x(0);e.gradient=Eigen::VectorXd::Ones(1);
        e.constraints=Eigen::VectorXd::Constant(1,x(0)-.75);
        e.constraintJacobian=Eigen::MatrixXd::Ones(1,1);return e;
    };
    spin2014::SlqpLbfgsSettings settings;
    settings.algorithm=spin2014::SlqpAlgorithm::AuditedKktRestorationV1;
    const auto result=spin2014::solveSlqpLbfgs(problem,Eigen::VectorXd::Zero(1),settings);
    std::cout << "status=" << spin2014::slqpStatusToken(result.status)
              << " x=" << result.beta.transpose() << " residual="
              << result.evaluation.constraints.norm() << std::endl;
    require(result.status==spin2014::SlqpStatus::Converged && std::abs(result.beta(0)-.75)<1e-8,
            "fixed equality restoration must release x=0 for min f=x, c=x-.75");
    require(result.finalKkt.valid && result.finalKkt.full<=1e-7 && std::abs(result.finalKkt.lambda(0)+1)<1e-7,
            "independent first-order KKT witness must certify the equality optimum");
    require(result.trace.front().restoration.attempted && result.trace.front().restoration.accepted,
            "bound release must use actual restoration");
    auto e=problem.evaluate(Eigen::VectorXd::Constant(1,1.));
    auto k=spin2014::evaluateKkt(e,Eigen::VectorXd::Constant(1,1.),1,
        Eigen::VectorXd::Constant(1,-1),Eigen::VectorXd::Zero(1),Eigen::VectorXd::Zero(1));
    require(k.valid && k.stationarity==0 && k.equalityPrimal==.25 && k.full>=.25,
            "stationarity alone cannot conceal equality infeasibility");
    independentlyCheck(result,1);suppliedWitnesses();analyticProblems();givensIntegration();
    std::cout<<"[PASS] audited actual KKT/restoration/analytic/Givens tests\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
