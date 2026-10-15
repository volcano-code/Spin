#include "spin2014/Spin2014GivensProblem.h"
#include "spin2014/AdaptiveVoxelSpinSolver.h"
#include "spin2014/DualMaterialSpinSolver.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

using namespace spin2014;
using V = Spin2014GivensProblem::Vector11;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
template<class F> void rejects(F f) { bool caught=false; try { f(); } catch(const std::invalid_argument&) { caught=true; } check(caught,"invalid input accepted"); }

std::shared_ptr<const voxel::VoxelDomain> fixture() {
    voxel::GridSpec grid;grid.dimensions=Eigen::Vector3i(4,4,4);
    grid.minimum=Eigen::Vector3d(-1.5,-3.,-1.5);grid.maximum=Eigen::Vector3d(1.5,3.,1.5);
    std::vector<unsigned char> material(64,1),carveable(64,0);
    for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x)carveable[x+4*(y+4*z)]=1;
    return voxel::VoxelDomain::create(grid,material,{},carveable);
}

void integration() {
    AdaptiveVoxelDomain domain(fixture(),2,2);
    auto problem=SpinIt2014MomentProblem::fromObjective(Spin2014Objective::TopInertiaOnly);
    AdaptiveFillFunctional f(domain,problem,0.03),unregularized(domain,problem,0.);
    Eigen::VectorXd beta=Eigen::VectorXd::Constant(domain.leaves().size(),0.2);
    for(int j=0;j<beta.size();++j)beta(j)=.2+.01*j;
    const double phi=1.7;const auto exact=f.evaluatePaper(beta,phi);
    check(exact.gradient.size()==beta.size()+1,"phi must be a joint unknown");
    for(int j=0;j<exact.gradient.size();++j) {
        auto plus=beta,minus=beta;double p=phi,m=phi,h=1e-5;
        if(j<beta.size()){plus(j)+=h;minus(j)-=h;}else{p+=h;m-=h;}
        auto a=f.evaluatePaper(plus,p),b=f.evaluatePaper(minus,m);
        check(std::abs((a.objective-b.objective)/(2*h)-exact.gradient(j))<1e-7,"functional objective chain rule");
        check(((a.constraints-b.constraints)/(2*h)-exact.constraintJacobian.col(j)).norm()<1e-7,"functional constraint chain rule");
    }
    const auto plain=unregularized.evaluatePaper(beta,phi);
    check(exact.regularizationEnergy>0.,"regularization fixture must be nonconstant");
    check((exact.gradient.head(beta.size())-plain.gradient.head(beta.size())-.03*f.laplacian()*beta).norm()<1e-13,"beta-only Laplacian gradient");
    check(std::abs(exact.gradient(beta.size())-plain.gradient(beta.size()))<1e-15,"phi received regularization");
    AdaptiveFillFunctional weighted(domain,problem,0.,2.,3.);
    const auto w=weighted.evaluatePaper(beta,phi);
    const auto expectedMoments=2.*domain.solidMoments().values()+domain.variableMoments().values()-3.*f.momentColumns()*beta;
    check((w.finalMoments.values()-expectedMoments).norm()<1e-12,"direct paper weighted adaptive moments");
    beta.setConstant(.2);
    SlqpLbfgsSettings settings;settings.maximumIterations=180;settings.equalityTolerance=1e-9;
    const auto result=solveSlqpLbfgsPaper(f,beta,phi,settings);
    const auto repeat=solveSlqpLbfgsPaper(f,beta,phi,settings);
    std::cout<<"joint status="<<slqpStatusToken(result.status)<<" residual="<<result.evaluation.constraints.lpNorm<Eigen::Infinity>()<<" final_phi="<<result.finalPhi<<" beta="<<result.beta.transpose()<<'\n';
    int accepted=0;for(const auto& t:result.trace){if(t.acceptedStep>0.)++accepted;std::cout<<"iteration="<<t.iteration<<" phi="<<t.phi<<" accepted_step="<<t.acceptedStep<<" residual="<<t.equalityInfinityNorm<<'\n';}
    check(accepted>0&&std::abs(result.finalPhi-phi)>1e-4,"no accepted angle optimization");
    check(result.status==SlqpStatus::Converged&&result.evaluation.constraints.lpNorm<Eigen::Infinity>()<=settings.equalityTolerance,"joint solver equality convergence");
    check(result.beta.size()==beta.size()&&result.finalPhi>1.,"angle clamped or returned as fill");
    check(result.traceHash==repeat.traceHash&&result.beta.isApprox(repeat.beta,0.)&&result.finalPhi==repeat.finalPhi,"joint deterministic replay");
    settings.maximumIterations=1;check(solveSlqpLbfgsPaper(f,beta,phi,settings).status!=SlqpStatus::Converged,"maxiter feasible certification");
    AdaptiveVoxelSpinRequest request;request.domain=fixture();request.objective=Spin2014Objective::TopInertiaOnly;
    request.initialLevel=2;request.maximumLevel=2;request.initialFill=.2;
    const auto legacy=computeAdaptiveVoxelSpinRequestHash(request);
    request.formulation=Spin2014Formulation::PaperGivens;request.initialPhi=phi;
    const auto identity=computeAdaptiveVoxelSpinRequestHash(request);check(identity!=legacy,"formulation identity absent");
    request.initialPhi=-.3;check(computeAdaptiveVoxelSpinRequestHash(request)!=identity,"phi identity absent");request.initialPhi=phi;
    const auto solved=solveAdaptiveVoxelSpin(request),again=solveAdaptiveVoxelSpin(request);
    check(solved.outputHash==again.outputHash&&solved.adaptiveTraceHash==again.adaptiveTraceHash,"adaptive replay");
    check(!solved.binaryFill.empty(),"joint binary validation not reached");
    V raw;raw.head<10>()=solved.finalMoments.values();raw(10)=solved.finalPhi;
    auto binary=Spin2014GivensProblem(*request.objective).evaluate(raw);
    check((binary.constraints-solved.binaryConstraints).norm()<1e-13,"binary selected evaluator/final phi");
    request.binaryConstraintTolerance=1e-14;
    const auto roundedFailure=solveAdaptiveVoxelSpin(request);
    check(roundedFailure.status==AdaptiveSpinStatus::BinaryValidationFailed,"binary strict tolerance must fail without repair");
    request.binaryConstraintTolerance=1e-6;
    request.slqp.maximumIterations=1;
    const auto nonconverged=solveAdaptiveVoxelSpin(request);
    check(nonconverged.status==AdaptiveSpinStatus::RelaxationDidNotConverge&&nonconverged.finalPhi!=phi,"failed relaxation angle provenance");
    request.slqp.maximumIterations=120;
    request.maximumLeafCount=1;
    const auto initialFailure=solveAdaptiveVoxelSpin(request);
    check(initialFailure.status==AdaptiveSpinStatus::ResourceLimitReached&&initialFailure.finalPhi==phi,"resource failure initial angle provenance");
    request.maximumLeafCount=1000000;
    // Same deterministic iteration allowance used by the converged joint
    // functional fixture above; this is not a paper-case parameter search.
    request.slqp.maximumIterations=180;
    DualMaterialSpinRequest dual;dual.hollowing=request;
    dual.boundaryDensity=1.;dual.interiorDensity=1.;
    const auto directPaper=solveAdaptiveVoxelSpin(request);
    const auto dualResult=solveDualMaterialSpin(dual);
    std::cout<<"dual paper status="<<adaptiveSpinStatusToken(dualResult.hollowing.status)
             <<" diagnostic="<<dualResult.hollowing.diagnostic
             <<" final_phi="<<dualResult.hollowing.finalPhi<<'\n';
    check(dualResult.hollowing.status==AdaptiveSpinStatus::Converged&&
          dualResult.hollowing.formulation==Spin2014Formulation::PaperGivens,
          "dual-material product must preserve a converged PAPER_GIVENS solve");
    check(dualResult.hollowing.outputHash==directPaper.outputHash&&
          dualResult.hollowing.finalPhi==directPaper.finalPhi,
          "equal-density dual-material product changed the paper Givens result");
    check(dualResult.hollowing.finalMoments.values().isApprox(
              (dualResult.boundaryContribution+dualResult.interiorContribution).values(),1e-11),
          "dual-material PAPER_GIVENS contributions do not reproduce final moments");
    request.initialLevel=1;
    const auto refined=solveAdaptiveVoxelSpin(request);
    check(refined.adaptiveTrace.size()>1,"fixture must exercise adaptive angle carry");
    for(std::size_t i=1;i<refined.adaptiveTrace.size();++i) {
        check(refined.adaptiveTrace[i].initialPhi==refined.adaptiveTrace[i-1].finalPhi,"optimized angle was reset between passes");
        check(refined.adaptiveTrace[i].slqpTrace.front().phi==refined.adaptiveTrace[i].initialPhi,"carried angle missing from actual solver trace");
    }
    std::cout<<"adaptive angle carry passes="<<refined.adaptiveTrace.size()<<'\n';
    request.mode=spin::SpinMode::StaticBalance;request.objective=Spin2014Objective::StaticLowering;
    rejects([&]{validateAdaptiveVoxelSpinRequest(request);});request.initialPhi=0.;
    AdaptiveFillFunctional stat(domain,SpinIt2014MomentProblem::staticBalance(),.03);
    check(stat.evaluatePaper(beta,0.).gradient.size()==beta.size(),"static has angle unknown");
    rejects([&]{solveSlqpLbfgsPaper(stat,beta,1.);});
    request.maximumLeafCount=1;auto failure=solveAdaptiveVoxelSpin(request);
    check(failure.status==AdaptiveSpinStatus::ResourceLimitReached,"resource fail closed");
    voxel::GridSpec badGrid;badGrid.dimensions=Eigen::Vector3i(3,3,3);
    badGrid.minimum=Eigen::Vector3d(-1.5,-3.,-1.5);badGrid.maximum=Eigen::Vector3d(1.5,3.,1.5);
    std::vector<unsigned char> mask(27,1),carveable(27,0);carveable[13]=1;
    AdaptiveVoxelDomain badDomain(voxel::VoxelDomain::create(badGrid,mask,{},carveable),2,2);
    rejects([&]{badDomain.maximumResolutionMaterialDomain(Eigen::VectorXd::Zero(badDomain.leaves().size()));});
}

int main() { try {
    std::cout << std::setprecision(17);
    V x; x << 2.,1./3.,-1./4.,1./2.,3.,2.,4.,1./7.,-1./5.,3./4.,std::atan2(4.,3.);
    const Spin2014Objective objectives[] = {Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering};
    const double expected[] = {107474./15625.,1655209./250000.,425521./62500.,1842709./250000.,3./16.};
    for(int k=0;k<5;++k) {
        Spin2014GivensProblem p(objectives[k],2.,3.);
        auto e=p.evaluate(x);
        check(std::abs(e.objective-expected[k])<1e-12,"independent rational objective");
        Eigen::VectorXd constraints(k==0?6:(k==4?2:5));
        if(k==0) constraints << 1./3.,-1./4.,1./2.,-1./5.,1./7.,27./100.;
        else if(k==4) constraints << 1./3.,-1./4.;
        else constraints << 1./3.,-1./4.,-1./5.,1./7.,27./100.;
        check((e.constraints-constraints).norm()<1e-13,"raw constraint ordering");
        if(k!=4) {
            Eigen::Matrix2d r,b;r<<.6,-.8,.8,.6;b<<6.,-.75,-.75,7.;
            check(std::abs((r*b*r.transpose())(0,1)+e.constraints.tail(1)(0))<1e-13,"R convention rotated offdiagonal=-g");
        }
        if(k!=4) {
            check(std::abs(e.objectiveGradient(10)-2322./15625.)<1e-12,"angle gradient oracle");
            check(std::abs(e.objectiveHessian(10,10)+13334./15625.)<1e-12,"angle Hessian oracle");
            auto old=SpinIt2014MomentProblem::fromObjective(objectives[k],2.,3.).evaluate(x.head<10>());
            check(std::abs(old.objective-e.objective)>1e-3,"off-manifold legacy difference");
        }
        double previous=1.;
        for(double h:{1e-3,3e-4,1e-4}) {
            double error=0.,secondOrderError=0.;
            for(int j=0;j<11;++j) {
                V a=x,b=x; a(j)+=h; b(j)-=h;
                const auto ep=p.evaluate(a),em=p.evaluate(b);
                a=x;b=x;a(j)+=2*h;b(j)-=2*h;
                const auto epp=p.evaluate(a),emm=p.evaluate(b);
                secondOrderError=std::max(secondOrderError,((ep.objectiveGradient-em.objectiveGradient)/(2*h)-e.objectiveHessian.col(j)).cwiseAbs().maxCoeff());
                // Fourth-order centered differences at the fixed three steps:
                // suppress known sin(4 phi) truncation without relaxing tolerance.
                auto compare=[&](double approx,double exact) {
                    double er=std::abs(approx-exact)/(std::abs(exact)>1e-7?std::abs(exact):1.);
                    if(er>1e-6)std::cerr<<"coordinate="<<j<<" exact="<<exact<<" approximate="<<approx<<" error="<<er<<'\n';
                    check(er<=1e-6,"11-coordinate finite difference tolerance"); error=std::max(error,er);
                };
                compare((8*(ep.objective-em.objective)-(epp.objective-emm.objective))/(12*h),e.objectiveGradient(j));
                for(int i=0;i<11;++i) compare((8*(ep.objectiveGradient(i)-em.objectiveGradient(i))-(epp.objectiveGradient(i)-emm.objectiveGradient(i)))/(12*h),e.objectiveHessian(i,j));
                for(int r=0;r<e.constraints.size();++r) {
                    compare((8*(ep.constraints(r)-em.constraints(r))-(epp.constraints(r)-emm.constraints(r)))/(12*h),e.constraintJacobian(r,j));
                    for(int i=0;i<11;++i) compare((8*(ep.constraintJacobian(r,i)-em.constraintJacobian(r,i))-(epp.constraintJacobian(r,i)-emm.constraintJacobian(r,i)))/(12*h),e.constraintHessians[r](i,j));
                }
            }
            check(error<=previous*1.05+1e-8,"finite difference trend above roundoff floor"); previous=error;
            std::cout << "objective="<<k<<" h="<<h<<" fourth_order_max_relative_or_zero_absolute_error="<<error<<" second_order_hessian_absolute_diagnostic="<<secondOrderError<<'\n';
        }
        for(double phi:{-2.3,1.7}) { V a=x;a(10)=phi;V b=a;b(10)+=std::acos(-1.);check(std::abs(p.evaluate(a).objective-p.evaluate(b).objective)<1e-12,"angle periodicity"); }
        V a=x;a(4)=a(5)=3.;a(9)=0.; auto repeated=p.evaluate(a);a(10)=-4.;check(std::abs(p.evaluate(a).objective-repeated.objective)<1e-12,"repeated transverse eigenvalues");
        a=x;a.head<10>()*=3.; double scaleExpected=expected[k]+(k==3?6.:0.);check(std::abs(p.evaluate(a).objective-scaleExpected)<1e-12,"density covariance");
        a=x;a(0)=0;rejects([&]{p.evaluate(a);});a=x;a(10)=std::numeric_limits<double>::infinity();rejects([&]{p.evaluate(a);});
        a=x;a(4)=-a(5);if(k!=4) rejects([&]{p.evaluate(a);});
        a=x;a(1)=a(2)=a(7)=a(8)=0.;if(k==0)a(3)=0.;
        a(10)=.5*std::atan2(-2*a(9),a(4)-a(5));
        const auto diagonal=p.evaluate(a);
        check(std::abs(diagonal.objective-SpinIt2014MomentProblem::fromObjective(objectives[k],2.,3.).evaluate(a.head<10>()).objective)<1e-12,"feasible diagonalization agreement");
        if(k!=4)check(std::abs(diagonal.constraints.tail(1)(0))<1e-12,"feasible g");
        for(int j=0;j<11;++j){a=x;a(j)=std::numeric_limits<double>::quiet_NaN();rejects([&]{p.evaluate(a);});}
    }
    integration();
    std::cout<<"[PASS] spin-it-2014-givens analytic and integration contract\n";
    return 0;
} catch(const std::exception& e) { std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;} }
