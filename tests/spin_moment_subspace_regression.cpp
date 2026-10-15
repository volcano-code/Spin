#include "spin/QuadricHollowingSolver.h"
#include "spin/DualMonomialCoordinates.h"
#include <Eigen/QR>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using V=Eigen::Matrix<double,10,1>;
using M=Eigen::Matrix<double,10,10>;
using Q=spin::QuadricCoefficients10;
bool dumpLedger=false;
void require(bool b,const char* m) {if(!b)throw std::runtime_error(m);}
Eigen::MatrixXd support() {
    Eigen::MatrixXd u=Eigen::MatrixXd::Zero(10,6);
    const int ids[]={0,1,2,3,7,8};
    for(int i=0;i<6;++i)u(ids[i],i)=1;
    return u;
}
Eigen::MatrixXd tangent(const Eigen::MatrixXd& u,const V& p) {
    Eigen::VectorXd a=u.transpose()*p;
    a(0)+=a(0)>=0?1:-1; a.normalize();
    return u*(Eigen::MatrixXd::Identity(a.size(),a.size())-2*a*a.transpose()).rightCols(a.size()-1);
}
// Independent literal full residual, not SpinMomentProblem derivatives.
Eigen::VectorXd residual(const V& r,const V& s,const Eigen::VectorXd& lambda,double mu) {
    double m=r(0),x=r(1),y=r(2),z=r(3);
    V g=V::Zero();g(0)=-z/(m*m);g(3)=1/m;
    Eigen::MatrixXd c=Eigen::MatrixXd::Zero(4,10);
    c(0,1)=1;c(1,2)=1;
    c(2,0)=-x*z/(m*m);c(2,1)=z/m;c(2,3)=x/m;c(2,8)=-1;
    c(3,0)=-y*z/(m*m);c(3,2)=z/m;c(3,3)=y/m;c(3,7)=-1;
    Eigen::VectorXd f(14);
    f.head<10>()=g-c.transpose()*lambda+mu*s;
    f.tail<4>()<<x,y,-r(8)+x*z/m,-r(7)+y*z/m;
    return f;
}
void regression() {
    V p;p<<.125,.3,.7,-1,0,0,0,.2,.4,0;p.normalize();
    V r;r<<2,.1,-.05,.8,3,4,5,.03,-.02,.1;
    const auto q=Q::fromUnnormalized(p);
    const auto c=spin::DualMonomialCoordinates10::fromGeometry(q);
    auto u=support();auto b=tangent(u,p);
    M map=.02*M::Identity();map(3,7)=.07;map(8,3)=-.04;
    const auto problem=spin::SpinMomentProblem::spinningTop(0,spin::SpinningTopObjective::LowerCenterOfMass);
    auto e=spin::linearizeSpinMomentSubspace(problem,q,r,u,b,map*b);
    require(e.matrix.rows()==10 && e.matrix.cols()==10,"supported system dimensions");
    require(e.momentTangentJacobian.rows()==10 && e.momentTangentJacobian.cols()==5,"derivative must retain five columns");
    auto f=residual(r,c.values(),e.multipliers,e.multiplierMu);
    require((e.fullResidual-f).norm()<1e-14,"full residual differs from literal four-constraint equations");
    require(std::abs(e.rawResidual-f.norm())<1e-14,"raw residual was projected or normalized");
    require((e.solverTangent.transpose()*e.solverTangent-Eigen::MatrixXd::Identity(5,5)).norm()<1e-13,"orthonormal tangent");
    require((e.solverTangent.transpose()*c.values()).norm()<1e-13,"not tangent to dual sphere");
    for(double h:{1e-4,3e-5,1e-5}) {
        Eigen::MatrixXd fd(10,10);
        for(int j=0;j<10;++j) {
            Eigen::VectorXd lplus=e.multipliers,lminus=e.multipliers;
            double muplus=e.multiplierMu,muminus=e.multiplierMu;
            V sp=c.values(),sm=sp,rp=r,rm=r;
            if(j<5) {
                sp=(sp+h*e.solverTangent.col(j)).normalized();
                sm=(sm-h*e.solverTangent.col(j)).normalized();
                V pp=sp,pm=sm;pp.tail<3>()*=.5;pm.tail<3>()*=.5;
                rp+=map*(pp.normalized()-p);rm+=map*(pm.normalized()-p);
            } else if(j<9) {lplus(j-5)+=h;lminus(j-5)-=h;}
            else {muplus+=h;muminus-=h;}
            const Eigen::VectorXd diff=(residual(rp,sp,lplus,muplus)-residual(rm,sm,lminus,muminus))/(2*h);
            fd.col(j).head(6)=u.transpose()*diff.head(10);fd.col(j).tail(4)=diff.tail(4);
        }
        double error=(e.matrix-fd).norm();
        if(error>=1e-6)std::cerr<<"difference at "<<h<<"\n"<<(e.matrix-fd)<<'\n';
        require(error<1e-6,"supported Eq14 matrix fails literal residual finite differences");
        std::cout<<"residual FD h="<<h<<" error="<<error<<'\n';
    }
    const auto solved=spin::solveSpinNewtonLinearSystem(e.matrix,e.rightHandSide);
    require(solved.solved && (e.matrix*solved.solution-e.rightHandSide).norm()<1e-10,"supported Newton equation unsolved");
    auto reject=[&](const Eigen::MatrixXd& badU,const Eigen::MatrixXd& badB,const Eigen::MatrixXd& badJ,const spin::SpinMomentProblem& pr) {
        bool caught=false;try{(void)spin::linearizeSpinMomentSubspace(pr,q,r,badU,badB,badJ);}catch(const std::invalid_argument&){caught=true;}
        require(caught,"invalid supported derivative accepted");
    };
    reject(u,b,Eigen::MatrixXd::Zero(10,10),problem);
    auto bad=b;bad.col(0)=bad.col(1);reject(u,bad,map*b,problem);
    bad=b;bad(4,0)=.1;reject(u,bad,map*b,problem);
    auto badU=u;badU.col(0)=badU.col(1);reject(badU,b,map*b,problem);
    reject(u,b,map*b,spin::SpinMomentProblem::yoYo());
    auto nan=(map*b).eval();nan(0,0)=std::numeric_limits<double>::quiet_NaN();reject(u,b,nan,problem);
    // Isolate the two distinct guards previously left without branch tests.
    bad=b;bad.col(0)=std::cos(.2)*b.col(0)+std::sin(.2)*V::Unit(4);
    try{(void)spin::linearizeSpinMomentSubspace(problem,q,r,u,bad,map*b);throw std::runtime_error("span guard missing");}
    catch(const std::invalid_argument& e){require(std::string(e.what()).find("span")!=std::string::npos,"span case rejected for wrong reason");}
    Eigen::MatrixXd wide=Eigen::MatrixXd::Zero(10,9);
    for(int i=0;i<4;++i)wide(i,i)=1;
    wide(4,4)=wide(5,4)=std::sqrt(.5);
    for(int i=6;i<10;++i)wide(i,i-1)=1;
    const auto wb=tangent(wide,p);
    V yr;yr<<2,0,0,0,1,1,1,0,0,0;
    Eigen::MatrixXd wj=Eigen::MatrixXd::Zero(10,8);wj(4,0)=.1;wj(5,0)=-.1;
    try{(void)spin::linearizeSpinMomentSubspace(spin::SpinMomentProblem::yoYo(),q,yr,wide,wb,wj);throw std::runtime_error("Hessian support guard missing");}
    catch(const std::invalid_argument& e){require(std::string(e.what()).find("gradient differential")!=std::string::npos,"Hessian case rejected before differential guard");}
    std::cout<<"[PASS] supported Eq14 assembly and rejection contracts\n";
}
void stream() {
    V p,r;for(int i=0;i<10;++i)std::cin>>p(i);for(int i=0;i<10;++i)std::cin>>r(i);
    Eigen::MatrixXd b(10,5),j(10,5);
    for(int a=0;a<10;++a)for(int k=0;k<5;++k)std::cin>>b(a,k);
    for(int a=0;a<10;++a)for(int k=0;k<5;++k)std::cin>>j(a,k);
    if(!std::cin)throw std::runtime_error("invalid replay input");
    auto e=spin::linearizeSpinMomentSubspace(spin::SpinMomentProblem::spinningTop(0,spin::SpinningTopObjective::LowerCenterOfMass),Q::fromUnnormalized(p),r,support(),b,j);
    std::cout<<std::setprecision(17);
    for(auto const* mat:{&e.solverTangent,&e.momentTangentJacobian,&e.matrix}) {
        for(int a=0;a<mat->rows();++a){for(int k=0;k<mat->cols();++k)std::cout<<(*mat)(a,k)<<' ';std::cout<<'\n';}
    }
    std::cout<<e.multipliers.transpose()<<'\n'<<e.multiplierMu<<'\n'<<e.rawResidual<<'\n';
}

void dump(const spin::SpinNewtonSubspaceResult& result,const char* name) {
    if(!dumpLedger)return;
    const auto matrix=[](const char* label,const auto& value) {
        std::cout<<label<<' '<<value.rows()<<' '<<value.cols();
        for(int row=0;row<value.rows();++row)for(int col=0;col<value.cols();++col)
            std::cout<<' '<<value(row,col);
        std::cout<<'\n';
    };
    std::cout<<std::hexfloat<<"BEGIN "<<name<<' '<<result.scientificHash<<' '
        <<spin::spinSolverStatusToken(result.solver.status)<<' '
        <<result.solver.functionEvaluations<<' '<<result.solver.restartCount<<'\n';
    matrix("support",result.supportBasis);
    for(const auto& e:result.solver.events) {
        std::cout<<"event "<<e.type<<' '<<e.restart<<' '<<e.iteration<<' '
            <<e.evaluationIndex<<' '<<e.restartEvaluationIndex<<' '
            <<e.outcome<<' '<<e.reason<<' '<<std::quoted(e.detail)<<'\n';
        if(e.rawEquationResidual)std::cout<<"event_raw "<<*e.rawEquationResidual<<'\n';
        if(e.step)std::cout<<"event_step "<<*e.step<<'\n';
    }
    for(const auto& record:result.evaluations) {
        std::cout<<"evaluation "<<record.evaluationIndex;
        for(double x:record.quadricCoefficients)std::cout<<' '<<x;
        std::cout<<' '<<std::quoted(record.errorType)<<' '<<std::quoted(record.error)<<'\n';
        if(record.evaluation) {
            const auto& e=*record.evaluation;
            matrix("raw_moments",e.material.rawMoments);
            matrix("unused_ambient",e.material.momentJacobian);
            matrix("native_basis",e.nativeTangentBasis);
            matrix("native_jacobian",e.nativeDirectionalJacobian);
            std::cout<<"map "<<e.material.derivativeDefined<<' '<<e.material.geometryHash<<' '<<e.material.topologyHash<<'\n';
        }
        if(record.linearization) {
            const auto& l=*record.linearization;
            matrix("dual_tangent",l.solverTangent);matrix("dual_jacobian",l.momentTangentJacobian);
            matrix("lambda",l.multipliers);matrix("residual",l.fullResidual);
            std::cout<<"mu_raw "<<l.multiplierMu<<' '<<l.rawResidual<<'\n';
            matrix("system",l.matrix);matrix("rhs",l.rightHandSide);
        }
    }
    std::cout<<"END "<<name<<'\n'<<std::defaultfloat;
}

void nonlinear() {
    const auto u=support();
    V start;start<<.125,.3,.7,-1,0,0,0,.2,.4,0;
    const auto q=Q::fromUnnormalized(start);
    M a=M::Zero();a(0,0)=.1;a(1,1)=a(2,2)=.2;a(3,3)=.3;
    a(7,7)=a(8,8)=.1;
    V offset;offset<<2,0,0,.8,3,4,5,0,0,0;
    auto map=[&](const Q& p) {
        spin::SpinMomentSubspaceMapEvaluation e;
        e.material.rawMoments=offset+a*p.values();
        e.material.derivativeDefined=true;
        e.material.geometryHash=p.hash();
        e.material.topologyHash=p.values()(8)>.1?"chart-a":"chart-b";
        e.nativeTangentBasis=tangent(u,p.values());
        e.nativeDirectionalJacobian=a*e.nativeTangentBasis;
        return e;
    };
    const auto problem=spin::SpinMomentProblem::spinningTop(0,spin::SpinningTopObjective::LowerCenterOfMass);
    spin::SpinSolverSettings settings;
    settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    settings.maximumRestarts=2;
    const auto run=spin::solveSpinMomentSubspaceSystem(problem,settings,u,map,{q});
    const auto repeat=spin::solveSpinMomentSubspaceSystem(problem,settings,u,map,{q});
    require(run.solver.status==spin::SpinSolverStatus::Converged,"nonlinear supported solver did not converge");
    require(run.solver.rawEquationResidual.value()<1e-12 && run.solver.trace.size()>2,"initial feasibility mistaken for optimization");
    require(run.evaluations.size()==static_cast<std::size_t>(run.solver.functionEvaluations),"directional evaluation ledger incomplete");
    require(run.scientificHash==repeat.scientificHash,"supported scientific repeat mismatch");
    require(run.solver.trace.front().rawEquationResidual.value()>1e-3,"fixture must need real optimization");
    // Verify the first evaluated trial against the separately FD-checked
    // restricted assembly, lifting five coordinates before the geodesic.
    const auto initial=map(q);
    const auto assembly=spin::linearizeSpinMomentSubspace(problem,q,initial.material.rawMoments,u,
        initial.nativeTangentBasis,initial.nativeDirectionalJacobian);
    const auto linear=spin::solveSpinNewtonLinearSystem(assembly.matrix,assembly.rightHandSide);
    require(linear.solved,"independent first-direction solve failed");
    const V dual=spin::DualMonomialCoordinates10::fromGeometry(q).values();
    V direction=assembly.solverTangent*linear.solution.head(5);
    direction-=dual*dual.dot(direction);
    const double step=std::min(direction.norm(),settings.maximumGeodesicStep);
    V expected=std::cos(step)*dual+std::sin(step)*direction.normalized();
    expected.tail<3>()*=.5;expected.normalize();
    V actual;for(int i=0;i<10;++i)actual(i)=run.evaluations[1].quadricCoefficients[i];
    require((actual-expected).norm()<1e-13,"first trial used wrong supported-space direction");
    bool transition=false;
    for(const auto& e:run.solver.events) {
        if(e.reason=="topology_transition_relinearized")transition=true;
        if(e.evaluationPayload)require(!e.evaluationPayload->solverMomentJacobian,"directional Jacobian padded in old payload");
    }
    require(transition,"accepted topology transition must use own derivative");
    for(const auto& record:run.evaluations) {
        require(record.evaluation && record.linearization,"valid evaluated state missing attachment");
        V p;for(int i=0;i<10;++i)p(i)=record.quadricCoefficients[i];
        for(int i:{4,5,6,9})require(p(i)==0,"geodesic left coefficient support");
        const auto& l=*record.linearization;
        const auto f=residual(offset+a*p,spin::DualMonomialCoordinates10::fromGeometry(Q::fromUnnormalized(p)).values(),l.multipliers,l.multiplierMu);
        require((l.fullResidual-f).norm()<1e-13,"stale residual/multiplier/geometry state");
    }
    // A real callback exception on the first trial must be retained and cause
    // backtracking, not erase a cost or abort all future attempts.
    int calls=0;
    auto rejecting=[&](const Q& p){if(++calls==2)throw std::domain_error("intentional invalid trial");return map(p);};
    const auto failedTrial=spin::solveSpinMomentSubspaceSystem(problem,settings,u,rejecting,{q});
    require(failedTrial.solver.status==spin::SpinSolverStatus::Converged,"defined trial rejection did not backtrack");
    require(!failedTrial.evaluations[1].evaluation && failedTrial.evaluations[1].error=="intentional invalid trial","failed callback lost");
    require(calls==failedTrial.solver.functionEvaluations,"exception call uncounted");
    calls=0;
    auto nonstandard=[&](const Q& p){if(++calls==2)throw 7;return map(p);};
    spin::SpinNewtonSubspaceResult nonstandardTrial;
    try{nonstandardTrial=spin::solveSpinMomentSubspaceSystem(problem,settings,u,nonstandard,{q});}
    catch(...){throw std::runtime_error("non-standard callback escaped and lost ledger");}
    require(nonstandardTrial.solver.status==spin::SpinSolverStatus::Converged,"non-standard trial did not backtrack");
    require(nonstandardTrial.evaluations[1].errorType=="non-standard" &&
        nonstandardTrial.evaluations[1].error=="non-standard evaluator exception","non-standard failure classification lost");
    require(calls==nonstandardTrial.solver.functionEvaluations &&
        nonstandardTrial.evaluations.size()==static_cast<std::size_t>(calls),"non-standard failure cost missing");
    auto budget=settings;budget.maximumFunctionEvaluationsPerRestart=2;
    const auto limited=spin::solveSpinMomentSubspaceSystem(problem,budget,u,map,{q});
    require(limited.solver.status!=spin::SpinSolverStatus::Converged && limited.solver.functionEvaluations==2,"budget bypass");
    V bad=start;bad(4)=.2;
    const auto restart=spin::solveSpinMomentSubspaceSystem(problem,settings,u,map,{Q::fromUnnormalized(bad),q});
    require(restart.solver.status==spin::SpinSolverStatus::Converged && restart.solver.restartCount==2 && !restart.evaluations.front().error.empty(),"unsupported initial state must fail and restart");
    auto zero=[&](const Q& p){auto e=map(p);e.nativeDirectionalJacobian.setZero();return e;};
    const auto insensitive=spin::solveSpinMomentSubspaceSystem(problem,settings,u,zero,{q});
    require(insensitive.solver.status!=spin::SpinSolverStatus::Converged && insensitive.solver.events.back().reason=="initial_sensitivity_vanished","zero directional sensitivity accepted");
    const auto negative=spin::solveSpinMomentSubspaceSystem(problem,settings,u,map,{Q::fromUnnormalized(-start),q});
    require(negative.solver.status==spin::SpinSolverStatus::Converged && negative.solver.restartCount==2 && negative.solver.events[1].reason=="initial_multiplier_nonpositive","negative initial mu not screened");
    settings.captureEvaluationPayload=true;
    const auto captured=spin::solveSpinMomentSubspaceSystem(problem,settings,u,map,{q});
    require(captured.scientificHash==run.scientificHash,"diagnostic timing/capture changed science hash");
    for(const auto& event:captured.solver.events)if(event.evaluationPayload)
        require(!event.evaluationPayload->solverMomentJacobian,"capture invented ambient derivative");
    auto timed=[&](const Q& p){auto e=map(p);e.material.residualEvaluationSeconds=123;e.material.derivativeEvaluationSeconds=456;return e;};
    require(spin::solveSpinMomentSubspaceSystem(problem,settings,u,timed,{q}).scientificHash==run.scientificHash,"reported timings entered science hash");
    auto wrong1=[&](const Q& p){auto e=map(p);e.material.momentJacobian(0,0)=1;return e;};
    auto wrong2=[&](const Q& p){auto e=map(p);e.material.momentJacobian(0,0)=2;return e;};
    const auto w1=spin::solveSpinMomentSubspaceSystem(problem,settings,u,wrong1,{q});
    const auto w2=spin::solveSpinMomentSubspaceSystem(problem,settings,u,wrong2,{q});
    require(w1.solver.status!=spin::SpinSolverStatus::Converged && w1.scientificHash!=w2.scientificHash,"different invalid evaluator returns alias scientific identity");
    dump(run,"success");dump(failedTrial,"trial_exception");dump(limited,"budget");
    dump(nonstandardTrial,"nonstandard_exception");
    dump(restart,"unsupported_start");dump(insensitive,"zero_sensitivity");
    dump(negative,"negative_mu");dump(captured,"capture");dump(w1,"invalid1");dump(w2,"invalid2");
    std::cout<<"[PASS] supported nonlinear iterations="<<run.solver.trace.size()-1<<" evaluations="<<run.solver.functionEvaluations<<" raw="<<*run.solver.rawEquationResidual<<" hash="<<run.scientificHash<<'\n';
}
}
int main(int argc,char** argv) {
    try{dumpLedger=argc==2 && std::string(argv[1])=="--ledger";if(argc>1 && !dumpLedger)stream();else {regression();nonlinear();}return 0;}
    catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}
}
