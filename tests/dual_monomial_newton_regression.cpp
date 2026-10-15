#include "spin/DualMonomialCoordinates.h"
#include "spin/QuadricHollowingSolver.h"
#include <Eigen/LU>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using V = Eigen::Matrix<double,10,1>;
using M = Eigen::Matrix<double,10,10>;
using Q = spin::QuadricCoefficients10;
using C = spin::DualMonomialCoordinates10;
void require(bool ok, const char* why) { if(!ok) throw std::runtime_error(why); }
V literal(V p) { p.tail<3>()*=2; return p.normalized(); }
Q geometry(V s) { s.tail<3>()*=.5; return Q::fromUnnormalized(s); }

void adapter() {
    const spin::SpinSolverSettings positional{7};
    require(positional.maximumIterations==7 && positional.equationProfile==spin::SpinEquationProfile::LegacyNativeNormalizedV1,
        "legacy aggregate settings source compatibility");
    V v; v << -.7,.2,-.3,.4,1.1,.8,1.4,.12,-.08,.17;
    const auto p=Q::fromUnnormalized(v);
    const auto c=C::fromGeometry(p);
    require((c.values()-literal(p.values())).norm()<1e-15,"literal monomial conversion");
    require((c.geometry().values()-p.values()).norm()<1e-15,"geometry round trip");
    require((C::fromUnnormalized(-c.values()).geometry().values()+p.values()).norm()<1e-15,"sign lost");
    const Eigen::Vector3d x(.4,-.8,1.3);
    V monomials; monomials << 1,x.x(),x.y(),x.z(),x.x()*x.x(),x.y()*x.y(),x.z()*x.z(),
        x.y()*x.z(),x.x()*x.z(),x.x()*x.y();
    V weighted=p.values(); weighted.tail<3>()*=2;
    require(std::abs(p.evaluate(x)-weighted.dot(monomials))<1e-15,"all cross-term polynomial duality");
    require(std::abs(c.values().dot(monomials)-p.evaluate(x)/weighted.norm())<1e-15,"dual polynomial normalization");
    std::cout << "[PASS] literal coefficient/raw-moment cross-term duality\n";
    for(double h : {3e-4,1e-4,3e-5}) {
        M fd;
        for(int j=0;j<10;++j) {
            V a=c.values(),b=a; a(j)+=h; b(j)-=h;
            fd.col(j)=(C::fromUnnormalized(a).geometry().values()-C::fromUnnormalized(b).geometry().values())/(2*h);
        }
        const double error=(fd-c.geometryJacobian()).norm();
        require(error<5*h*h,"adapter normalization differential FD");
        std::cout << "adapter h=" << h << " error=" << error << '\n';
    }
    require((c.geometryJacobian()*c.values()).norm()<1e-15,"adapter radial null");
    std::cout << "[PASS] literal coefficient-to-geometry differential\n";
    const V d=c.projectToTangent(V::LinSpaced(10,-.3,.8));
    require(std::abs(c.values().dot(d))<1e-15,"dual tangent projection");
    const auto moved=c.geodesicUpdate(d,.25);
    const V expected=std::cos(.25)*c.values()+std::sin(.25)*d.normalized();
    require((moved.values()-expected).norm()<1e-15 && std::abs(moved.values().norm()-1)<1e-15,"dual geodesic");
    for(V bad : {V::Zero().eval(), (V::Constant(std::numeric_limits<double>::infinity())).eval()}) {
        bool rejected=false; try { (void)C::fromUnnormalized(bad); } catch(const std::invalid_argument&) {rejected=true;}
        require(rejected,"invalid coordinates accepted");
    }
}

void newton() {
    // The published-mode pairing diagnostic is only an algebraic map fixture.
    V r; r << 1,0,0,0,2,3,4,0,0,.5;
    V g=V::Zero(); g(4)=-111./125; g(5)=-101./125; g(6)=26./25; g(9)=2./25;
    const auto problem=spin::SpinMomentProblem::yoYo();
    const auto evaluation=problem.evaluate(r);
    require((evaluation.objectiveGradient-g).norm()<1e-14,"independent yo-yo rational gradient");
    const auto stationary=geometry(-g);
    const Eigen::Vector3d x(1,1,std::sqrt(20./13));
    require(stationary.evaluate(x)>0 && Q::fromUnnormalized(-g).evaluate(x)<0,"yo-yo pairing sign diagnostic");
    spin::SpinSolverSettings strict; strict.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    strict.maximumRestarts=1; strict.maximumFunctionEvaluationsPerRestart=1;
    const auto mapAt=[&](const Q& center) {
        return [&,center](const Q& p) {
            spin::SpinMomentMapEvaluation e;
            e.rawMoments=r+.01*(p.values()-center.values());
            e.momentJacobian=.01*M::Identity(); e.derivativeDefined=true;
            return e;
        };
    };
    const auto solved=spin::solveSpinMomentSystem(problem,strict,mapAt(stationary),{stationary});
    require(solved.status==spin::SpinSolverStatus::Converged && solved.rawEquationResidual && *solved.rawEquationResidual<1e-12,"strict paired stationary state");
    auto legacy=strict; legacy.equationProfile=spin::SpinEquationProfile::LegacyNativeNormalizedV1;
    const auto old=spin::solveSpinMomentSystem(problem,legacy,mapAt(stationary),{stationary});
    require(old.status!=spin::SpinSolverStatus::Converged && old.events[1].normalizedEquationResidual.value()>1e-3,"legacy pairing distinction through evaluator events");
    require(solved.eventHash!=old.eventHash && solved.traceHash!=old.traceHash,"profile hash alias");

    V start=-g; start(4)+=.16; start(9)-=.1;
    const V s=start.normalized(); const Q p=geometry(s);
    strict.maximumFunctionEvaluationsPerRestart=60;
    const auto run=spin::solveSpinMomentSystem(problem,strict,mapAt(p),{p});
    const auto repeat=spin::solveSpinMomentSystem(problem,strict,mapAt(p),{p});
    require(run.events.size()>3 && run.events[2].step && run.events[2].quadricCoefficients,"Newton trial absent");
    // Independent reduced 4-coordinate bordered system. Constraints force the
    // other active coordinates to zero at L=Qxz=Qyz=0.
    Eigen::Vector4d n(12,14,26,2), w(1,1,0,0), ss(s(4),s(5),s(6),s(9));
    Eigen::Matrix4d hn=Eigen::Matrix4d::Zero(); hn.diagonal()<<2,2,4,4;
    hn(0,2)=hn(2,0)=hn(1,2)=hn(2,1)=2;
    const Eigen::Matrix4d h=hn/25-2*(n*w.transpose()+w*n.transpose())/125+6*85.5*w*w.transpose()/625;
    V ts=s; ts.tail<3>()*=.5;
    M t=M::Identity(); t.diagonal().tail<3>().setConstant(.5);
    const M chain=.01*(M::Identity()-p.values()*p.values().transpose())*t/ts.norm();
    const int indices[4]={4,5,6,9}; Eigen::Matrix4d j;
    for(int a=0;a<4;++a) for(int b=0;b<4;++b) j(a,b)=chain(indices[a],indices[b]);
    const double mu=-g.dot(s);
    Eigen::Matrix<double,5,5> system=Eigen::Matrix<double,5,5>::Zero();
    system.topLeftCorner<4,4>()=h*j+mu*Eigen::Matrix4d::Identity();
    system.topRightCorner<4,1>()=ss; system.bottomLeftCorner<1,4>()=2*ss.transpose();
    Eigen::Matrix<double,5,1> rhs; rhs << -(Eigen::Vector4d(g(4),g(5),g(6),g(9))+mu*ss),0;
    const Eigen::Vector4d reduced=system.fullPivLu().solve(rhs).head<4>();
    V d=V::Zero(); for(int a=0;a<4;++a) d(indices[a])=reduced(a);
    d-=s*s.dot(d); const double step=std::min(.5,d.norm());
    const V expected=std::cos(step)*s+std::sin(step)*d.normalized();
    V actual; for(int i=0;i<10;++i) actual(i)=run.events[2].quadricCoefficients.value()[i];
    require(std::abs(run.events[2].step.value()-step)<1e-13 && (literal(actual)-expected).norm()<1e-12,"actual Newton direction differs from independent reduced system");
    double previous=std::numeric_limits<double>::infinity();
    for(const auto& e:run.events) {
        require(e.equationProfile==strict.equationProfile,"event profile lost");
        if(e.type=="evaluation" && e.rawEquationResidual) {
            V q; for(int i=0;i<10;++i) q(i)=e.quadricCoefficients.value()[i];
            const auto state=problem.evaluate(r+.01*(q-p.values()));
            const double freshMu=-state.objectiveGradient.dot(literal(q));
            require(std::abs(e.multiplierMu.value()-freshMu)<1e-12,"multipliers not recomputed");
            if(e.outcome=="accepted") {require(*e.rawEquationResidual<previous,"strict raw line-search decrease");previous=*e.rawEquationResidual;}
        }
    }
    require(run.eventHash==repeat.eventHash && run.traceHash==repeat.traceHash && run.eventStreamVersion==2,"strict scientific determinism");
    require(run.status==spin::SpinSolverStatus::Converged && run.rawEquationResidual.value()<1e-12,"strict nonlinear algebraic Newton solve");
    for(const auto& e:run.trace) require(e.quadricCoefficients && e.equationProfile==strict.equationProfile && e.rawEquationResidual,"strict trace coordinate identity");

    const auto staticProblem=spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    const auto failed=spin::solveSpinMomentSystem(staticProblem,strict,[](const Q&) -> spin::SpinMomentMapEvaluation {throw std::runtime_error("fixture failure");},{p});
    require(!failed.rawEquationResidual && !failed.events[1].rawEquationResidual && failed.equationProfile==strict.equationProfile,"failed residual null/profile");
    const auto nonfinite=spin::solveSpinMomentSystem(staticProblem,strict,[](const Q&) {spin::SpinMomentMapEvaluation e;e.rawMoments.setConstant(std::numeric_limits<double>::quiet_NaN());return e;},{p});
    require(!nonfinite.rawEquationResidual && nonfinite.events[1].reason=="invalid_numerical_state","strict nonfinite guard");
    strict.equationProfile=static_cast<spin::SpinEquationProfile>(99);
    bool rejected=false;try{(void)spin::solveSpinMomentSystem(problem,strict,mapAt(p),{p});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid equation profile accepted");
    std::cout << "strict nonlinear evaluations=" << run.functionEvaluations << " raw=" << run.rawEquationResidual.value() << " event=" << run.eventHash << '\n';
}

void stoppingAndSensitivity() {
    spin::SpinSolverSettings strict; strict.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    strict.maximumRestarts=1; strict.maximumFunctionEvaluationsPerRestart=2;
    V start=V::Zero();start(0)=-.96209441;start(3)=-.22572537;start(4)=-.15304383;
    const Q p=Q::fromUnnormalized(start);
    V r=V::Zero();r(0)=.54149167;r(3)=-1.13932953;
    M a=.1*M::Identity();a(0,0)=-2.59280884;a(3,3)=-.93656388;a(4,4)=-3.68855061;
    const auto map=[&](const Q& q) {spin::SpinMomentMapEvaluation e;
        e.rawMoments=r+a*(q.values()-p.values());e.momentJacobian=a;e.derivativeDefined=true;return e;};
    const auto problem=spin::SpinMomentProblem::staticBalance(spin::StaticObjective::LowerCenterOfMass);
    const auto run=spin::solveSpinMomentSystem(problem,strict,map,{p});
    require(run.events.size()==4 && run.events[2].reason=="no_residual_decrease" &&
        run.events[2].outcome=="rejected","strict must reject raw increase despite normalized decrease");
    require(*run.events[2].rawEquationResidual>*run.events[1].rawEquationResidual &&
        *run.events[2].normalizedEquationResidual<*run.events[1].normalizedEquationResidual,
        "line-search adversarial fixture must separate raw and normalized decrease");
    auto legacy=strict;legacy.equationProfile=spin::SpinEquationProfile::LegacyNativeNormalizedV1;
    const auto old=spin::solveSpinMomentSystem(problem,legacy,map,{p});
    require(old.events[2].outcome=="accepted","legacy normalized line-search changed");
    std::cout << "raw line-search before=" << *run.events[1].rawEquationResidual << " trial=" << *run.events[2].rawEquationResidual
              << " normalized before=" << *run.events[1].normalizedEquationResidual << " trial=" << *run.events[2].normalizedEquationResidual << '\n';

    V unit=V::Zero();unit(0)=1;const auto q=Q::fromUnnormalized(unit);
    const auto radial=[](const Q& p) {spin::SpinMomentMapEvaluation e;
        e.rawMoments(0)=1+.5*(p.values()(0)*p.values()(0)-1);
        e.momentJacobian(0,0)=p.values()(0);e.derivativeDefined=true;return e;};
    const auto mass=spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    const auto insensitive=spin::solveSpinMomentSystem(mass,strict,radial,{q});
    require(insensitive.status!=spin::SpinSolverStatus::Converged &&
        insensitive.events.back().reason=="initial_sensitivity_vanished","native radial-only Jacobian must not satisfy dual sensitivity");
    const auto overflow=spin::solveSpinMomentSystem(mass,strict,[](const Q&) {
        spin::SpinMomentMapEvaluation e;e.rawMoments(0)=1;e.momentJacobian(0,0)=1e308;
        e.derivativeDefined=true;return e;},{q});
    require(overflow.events[1].reason=="invalid_numerical_state" && !overflow.events[1].rawEquationResidual,
        "coordinate projection must not conceal overflowing native Jacobian norm");
    strict.maximumFunctionEvaluationsPerRestart=61;
    bool rejected=false;try{(void)spin::solveSpinMomentSystem(mass,strict,radial,{q});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"strict evaluation budget above 60 accepted");
}

void supportedModes() {
    V r; r << 1,0,0,0,2,3,4,0,0,.5;
    V g=V::Zero();g(4)=-111./125;g(5)=-101./125;g(6)=26./25;g(9)=2./25;
    const std::vector<spin::SpinMomentProblem> problems={
        spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass),
        spin::SpinMomentProblem::fidgetSpinner(),spin::SpinMomentProblem::yoYo(),
        spin::SpinMomentProblem::spinningTop(0),spin::SpinMomentProblem::buoyantStability(1,Eigen::Vector2d::Zero())};
    for(std::size_t i=0;i<problems.size();++i) {
        V s=-g;
        if(i==0){s.setZero();s(0)=1;}
        if(i==1){s.setZero();s(4)=s(5)=1;}
        if(i==4){s.setZero();s(3)=-1;}
        const Q p=geometry(s);
        spin::SpinSolverSettings settings;settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;settings.maximumRestarts=1;
        const auto result=spin::solveSpinMomentSystem(problems[i],settings,[&](const Q& q){
            spin::SpinMomentMapEvaluation e;e.rawMoments=r+.01*(q.values()-p.values());
            e.momentJacobian=.01*M::Identity();e.derivativeDefined=true;return e;},{p});
        require(result.status==spin::SpinSolverStatus::Converged && result.rawEquationResidual.value()<1e-12,
            "supported-mode strict algebraic known solution");
    }
}
}

int runDualMonomialRegression() {
    try {adapter();newton();stoppingAndSensitivity();supportedModes();std::cout << "[PASS] dual monomial coordinate Newton\n";return 0;}
    catch(const std::exception& e) {std::cerr << "[FAIL] " << e.what() << '\n';return 1;}
}
