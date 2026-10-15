#include "spin/MassMomentProblem.h"
#include "spin/QuadricHollowingSolver.h"
#include "spin2014/QuadricMomentProblem.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using V = Eigen::Matrix<double, 10, 1>;
using M = Eigen::Matrix<double, 10, 10>;
using O = spin2014::Spin2014Objective;
void require(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
V input() { V s; s << 2,.2,-.4,1,3,4,5,.15,-.1,.25; return s; }

void formulas() {
    // Independent hand integration coordinates: m=2, c=(.1,-.2,.5).
    // Central I=(8.42,7.48,6.90); Ixy=-.29, Ixz=.20, Iyz=-.35.
    const double ratio = 127.015 / 47.61;
    const O objectives[] = {O::YoYo,O::TopInertiaOnly,O::TopLoweringOnly,
                            O::TopMassAndLowering,O::StaticLowering};
    const double expected[] = {ratio,ratio,ratio+.25,ratio+1,.25};
    const int constraintCounts[] = {5,4,4,4,2};
    for (int k=0; k<5; ++k) {
        spin2014::QuadricObjectiveSettings settings;
        settings.objective=objectives[k];
        settings.declaredLengthUnit="model_length";
        settings.declaredMassUnit="model_mass";
        const auto problem=spin2014::makeSpin2014QuadricMomentProblem(settings);
        const auto e=problem.evaluate(input());
        require(std::abs(e.objective-expected[k])<1e-13,"2014 objective replaced by an approximate Faster objective");
        require(e.constraints.size()==constraintCounts[k],"physical mode constraints aliased");
        require(std::abs(e.constraints(0)-.2)<1e-14 && std::abs(e.constraints(1)+.4)<1e-14,
                "target first moment constraints lost");
        if (k<4) {
            const int offset=k==0?3:2;
            require(std::abs(e.constraints(offset)-.2)<1e-14 &&
                    std::abs(e.constraints(offset+1)+.35)<1e-14,"non-diagonal central inertia incorrect");
        }
        for(double h : {1e-4,3e-5,1e-5}) {
            V gradient; M hessian;
            for(int j=0;j<10;++j) {
                V plus=input(),minus=input(); plus(j)+=h;minus(j)-=h;
                const auto a=problem.evaluate(plus),b=problem.evaluate(minus);
                gradient(j)=(a.objective-b.objective)/(2*h);
                hessian.col(j)=(a.objectiveGradient-b.objectiveGradient)/(2*h);
            }
            require((gradient-e.objectiveGradient).norm()/std::max(1.,e.objectiveGradient.norm())<1e-6,
                    "objective gradient does not differentiate actual 2014 objective");
            require((hessian-e.objectiveHessian).norm()/std::max(1.,e.objectiveHessian.norm())<1e-6,
                    "objective Hessian does not differentiate actual gradient");
        }
        settings.solverMassUnitInDeclaredUnits=2;
        settings.solverLengthUnitInDeclaredUnits=3;
        V d; d << 2,6,6,6,18,18,18,18,18,18;
        const auto converted=spin2014::makeSpin2014QuadricMomentProblem(settings);
        const auto a=converted.evaluate(input().cwiseQuotient(d));
        require(std::abs(a.objective-e.objective)<1e-13,"normalized moments changed declared-unit objective");
        require((a.objectiveGradient-d.asDiagonal()*e.objectiveGradient).norm()<1e-12,
                "declared-unit gradient chain rule missing");
        require((a.objectiveHessian-d.asDiagonal()*e.objectiveHessian*d.asDiagonal()).norm()<1e-11,
                "declared-unit Hessian chain rule missing");
        require((a.constraints-e.constraints).norm()<1e-13 &&
                (a.constraintJacobian-e.constraintJacobian*d.asDiagonal()).norm()<1e-12,
                "constraint coordinate transform missing");
        for(Eigen::Index j=0;j<a.constraints.size();++j)
            require((a.constraintHessians[j]-d.asDiagonal()*e.constraintHessians[j]*d.asDiagonal()).norm()<1e-12,
                    "constraint Hessian coordinate transform missing");
        for(double h : {1e-5,3e-6,1e-6}) {
            V gradient; M hessian;
            const V normalized=input().cwiseQuotient(d);
            for(int j=0;j<10;++j) {
                V plus=normalized,minus=normalized;plus(j)+=h;minus(j)-=h;
                const auto a=converted.evaluate(plus),b=converted.evaluate(minus);
                gradient(j)=(a.objective-b.objective)/(2*h);
                hessian.col(j)=(a.objectiveGradient-b.objectiveGradient)/(2*h);
            }
            require((gradient-a.objectiveGradient).norm()/std::max(1.,a.objectiveGradient.norm())<1e-6,
                    "normalized full gradient failed finite differences");
            require((hessian-a.objectiveHessian).norm()/std::max(1.,a.objectiveHessian.norm())<1e-6,
                    "normalized full Hessian failed finite differences");
        }
    }
    spin2014::QuadricObjectiveSettings settings;
    settings.objective=O::TopMassAndLowering;
    const auto originalHash=spin2014::spin2014QuadricObjectiveHash(settings);
    settings.objective=O::TopLoweringOnly;
    require(originalHash!=spin2014::spin2014QuadricObjectiveHash(settings),"strategies share request identity");
    settings.objective=O::TopMassAndLowering;
    settings.declaredMassUnit="kg";
    require(originalHash!=spin2014::spin2014QuadricObjectiveHash(settings),"units absent from objective identity");
    for(double bad : {0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        settings.solverLengthUnitInDeclaredUnits=bad;
        bool rejected=false;try { (void)spin2014::makeSpin2014QuadricMomentProblem(settings); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"invalid unit scale accepted");
    }
    std::cout << "[PASS] five original objectives, non-diagonal inertia, declared-unit chain rule and three-step differences\n";
}

void kernel() {
    V p; p << 1,.12,-.08,.09,-.07,.05,-.04,.03,-.02,.01;
    const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
    const auto identity=[](const auto& q) {
        spin::SpinMomentMapEvaluation map;map.rawMoments=q.values();map.momentJacobian.setIdentity();
        map.derivativeDefined=true;map.geometryHash=q.hash();map.topologyHash="synthetic-identity";return map;
    };
    const auto oldProblem=spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass);
    spin::SpinSolverSettings settings;
    const auto old=spin::solveSpinMomentSystem(oldProblem,settings,identity,{q});
    const auto shared=spin::solveMassMomentSystem(spin::MassMomentProblem(oldProblem),settings,identity,{q});
    require(old.status==shared.status && old.traceHash==shared.traceHash && old.eventHash==shared.eventHash &&
            old.functionEvaluations==shared.functionEvaluations,"common interface changed historical Newton path");

    settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    settings.maximumRestarts=1;settings.maximumFunctionEvaluationsPerRestart=1;
    for(O objective : {O::YoYo,O::TopInertiaOnly,O::TopLoweringOnly,O::TopMassAndLowering,O::StaticLowering}) {
        spin2014::QuadricObjectiveSettings objectiveSettings;objectiveSettings.objective=objective;
        const auto problem=spin2014::makeSpin2014QuadricMomentProblem(objectiveSettings);
        const auto map=[](const auto& q) {
            spin::SpinMomentMapEvaluation value;value.rawMoments=input();
            value.momentJacobian=.01*M::Identity();value.derivativeDefined=true;value.geometryHash=q.hash();return value;
        };
        const auto result=spin::solveMassMomentSystem(problem,settings,map,{q});
        bool actual=false;
        for(const auto& event:result.events)
            if(event.objective) { require(std::abs(*event.objective-problem.evaluate(input()).objective)<1e-13,
                "Newton used the wrong objective factory");actual=true; }
        require(actual && result.functionEvaluations==1,"2014 problem not consumed by actual Newton evaluation");
    }
    bool rejected=false;
    try { spin::MassMomentProblem problem{spin::MassMomentProblem::Evaluator{}}; }
    catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"empty mass moment problem accepted");
    std::cout << "[PASS] unchanged historical kernel and actual five-objective Newton calls (synthetic, not material success)\n";
}
}

int main(int argc,char** argv) {
    try {
        const std::string scenario=argc>1?argv[1]:"--all";
        if(scenario=="--formulas" || scenario=="--all") formulas();
        if(scenario=="--kernel" || scenario=="--all") kernel();
        return EXIT_SUCCESS;
    } catch(const std::exception& e) {std::cerr << "[FAIL] " << e.what() << '\n';return EXIT_FAILURE;}
}
