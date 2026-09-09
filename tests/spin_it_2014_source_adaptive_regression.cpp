#if __has_include("spin2014/SourceAdaptiveSpinSolver.h")
#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "spin2014/SlqpAuditTrace.h"
#include "source_surface_test_utils.h"
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
#include <set>
#include <sys/resource.h>
using namespace spin2014;
using namespace surface_test;
std::shared_ptr<const SourceSurfaceOctree> source(const Mesh& outer,const Mesh& design,int depth=5,
    Eigen::Vector3d origin={0,0,0},double side=1) {
    SourceSurfaceOctreeConfig c;c.maximumDepth=depth;c.rootOrigin=origin;c.rootSide=side;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
}
void graph(const SourceAdaptiveDomain& d) {
    std::set<std::pair<std::size_t,std::size_t>> expected;
    for(std::size_t i=0;i<d.leaves().size();++i) for(std::size_t j=i+1;j<d.leaves().size();++j) {
        const auto& a=d.leaves()[i];const auto& b=d.leaves()[j];int touching=0,overlap=0;
        for(int k=0;k<3;++k) {
            const auto al=a.coordinates[k]<<(20-a.level), ah=(a.coordinates[k]+1)<<(20-a.level);
            const auto bl=b.coordinates[k]<<(20-b.level), bh=(b.coordinates[k]+1)<<(20-b.level);
            touching+=(ah==bl || bh==al);overlap+=std::min(ah,bh)>std::max(al,bl);
        }
        if(touching==1 && overlap==2) expected.emplace(i,j);
    }
    require(std::set<std::pair<std::size_t,std::size_t>>(d.edges().begin(),d.edges().end())==expected,"sparse face graph disagrees with independent integer boxes");
    Eigen::MatrixXd l=d.uniformLaplacian();require(l.isApprox(l.transpose()) && l.rowwise().sum().norm()<1e-12,"symmetric zero row sum Laplacian");
    if(l.rows()) require(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(l).eigenvalues().minCoeff()>-1e-10,"PSD Laplacian");
    Eigen::VectorXd b(l.rows());for(int i=0;i<b.size();++i)b(i)=.01*(i%19);
    double energy=0;for(auto e:expected)energy+=std::pow(b(e.first)-b(e.second),2);
    require(std::abs(b.dot(l*b)-energy)<1e-10,"edge energy");
}
void partitionTests() {
    auto s=source(cube(0,1),cube(.125,.875));SourceAdaptiveDomain d(s,2,5);
    graph(d);Eigen::VectorXd b=Eigen::VectorXd::Constant(d.leaves().size(),.37);
    const auto before=d.removedMoments(b);const auto identity=d.partitionHash();
    auto split=d.splitFractionalLeaves(b,.05);require(split.changes.size()>0 && d.partitionHash()!=identity,"actual fractional split");
    require((b.array()==.37).all(),"exact beta inheritance");close(d.removedMoments(b).values(),before.values(),"split conservation all ten");graph(d);d.validateCoverage();
    b.setConstant(.01);auto merge=d.mergeBinarySiblings(b,.05);
    require(merge.changes.size()>1 && d.leaves().size()==1 && d.leaves()[0].level==0 && b(0)==.01,"recursive merge below initial level without snapping");
    close(d.removedMoments(b).values(),s->carveableMoments().scaledDensity(.01).values(),"near bound mass preserved");
    SourceAdaptiveDomain unequal(s,2,5);b=Eigen::VectorXd::Constant(unequal.leaves().size(),.01);
    for(int i=0;i<b.size();++i)b(i)+=i*1e-8;
    const auto n=unequal.leaves().size();unequal.mergeBinarySiblings(b,.05);require(unequal.leaves().size()==n,"unequal near-bound values must not merge");
    // Mixed levels, holes and face/corner contacts exercised against exhaustive boxes.
    SourceAdaptiveDomain mixed(s,2,5);b=Eigen::VectorXd::Zero(mixed.leaves().size());b(0)=.37;
    const auto mixedBefore=mixed.removedMoments(b);mixed.splitFractionalLeaves(b,.05);graph(mixed);mixed.validateCoverage();
    close(mixed.removedMoments(b).values(),mixedBefore.values(),"nonuniform beta conservation");
    auto cavity=source(cube(0,1),join(cube(.125,.875),cube(.375,.625),true));
    SourceAdaptiveDomain holes(cavity,3,5);graph(holes);holes.validateCoverage();
    SourceAdaptiveDomain cap(s,2,5,n);b=Eigen::VectorXd::Constant(cap.leaves().size(),.5);
    auto h=cap.partitionHash();const Eigen::VectorXd saved=b;
    rejects([&]{cap.splitFractionalLeaves(b,.05);},"leaf limit must fail");
    require(h==cap.partitionHash() && b==saved,"strong leaf budget atomicity");
    b(0)=std::numeric_limits<double>::quiet_NaN();rejects([&]{cap.splitFractionalLeaves(b,.05);},"NaN beta rejected");require(h==cap.partitionHash() && std::isnan(b(0)),"invalid beta atomicity");
    b=saved;rejects([&]{cap.mergeBinarySiblings(b,.5);},"epsilon rejected");require(h==cap.partitionHash() && b==saved,"invalid epsilon atomicity");
    // Non-default represented root: every queried direct-child sum and stored-parent identity.
    auto shifted=source(box({-.2,-.1,0},{1,1,1}),cube(.25,.75),4,{-.3,-.2,-.1},1.4);
    for(const auto& cell:shifted->cells()) if(cell.level<4) {
        Ten total=Ten::Zero();for(int k=0;k<8;++k) {auto xyz=cell.coordinates;for(int a=0;a<3;++a)xyz[a]=xyz[a]*2+((k>>a)&1);auto c=shifted->cell(cell.level+1,xyz);total+=c.moments.values();require(c.parent==cell.index,"shifted query parent");}
        if(cell.moments.mass())close(total,cell.moments.values(),"shifted parent descendant conservation");else require(total.isZero(),"zero descendants");
    }
}
SourceAdaptiveSpinRequest fixedRequest(bool shifted=true) {
    const double offset=shifted?1./8192:0.;SourceAdaptiveSpinRequest r;
    r.source=source(box({-.5+offset,-.5,0},{.5+offset,.5,1}),box({.125,-.125,.625},{.375,.125,.875}),5,{-1,-1,0},2);
    r.objective=Spin2014Objective::StaticLowering;r.initialLevel=2;r.maximumOptimizationLevel=5;
    r.laplacianWeight=1;r.slqp.maximumIterations=300;r.slqp.equalityTolerance=1e-10;r.slqp.projectedGradientTolerance=1e-8;return r;
}
void functionalTests() {
    auto r=fixedRequest();SourceAdaptiveDomain d(r.source,2,5);require(d.leaves().size()==2,"fixed two source leaves");
    close(r.source->carveableMoments().values(),boxMoments({3./16,-1./16,11./16},{5./16,1./16,13./16}),"fixed conservative accepted box");
    const Ten outer=boxMoments({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    Eigen::VectorXd b(2);b<<.2,.4;
    Ten removed=.2*boxMoments({3./16,-1./16,11./16},{5./16,0,13./16})+.4*boxMoments({3./16,0,11./16},{5./16,1./16,13./16});
    for(auto objective:{Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering}) {
        AdaptiveFillFunctional f(d,SpinIt2014MomentProblem::fromObjective(objective,1.3,.7),.8,2.);
        double phi=objective==Spin2014Objective::StaticLowering?0:.23;auto e=f.evaluatePaper(b,phi);
        close(e.finalMoments.values(),2*(outer-removed),"original outer density and analytic source columns");
        close(f.evaluatePaper(Eigen::VectorXd::Zero(2),phi).finalMoments.values(),2*outer,"beta0 original outer");
        close(f.evaluatePaper(Eigen::VectorXd::Ones(2),phi).finalMoments.values(),2*(outer-boxMoments({3./16,-1./16,11./16},{5./16,1./16,13./16})),"beta1 retains unresolved boundary");
        const auto raw=e.finalMoments.values();double c=std::cos(phi),sn=std::sin(phi),delta=objective==Spin2014Objective::YoYo?0:raw(3)*raw(3)/raw(0);
        double ia=c*c*(raw(5)+raw(6))+sn*sn*(raw(4)+raw(6))+2*c*sn*raw(9)-delta;
        double ib=sn*sn*(raw(5)+raw(6))+c*c*(raw(4)+raw(6))-2*c*sn*raw(9)-delta;
        double expected=objective==Spin2014Objective::StaticLowering?0:1.3*(ia*ia+ib*ib)/std::pow(raw(4)+raw(5),2);
        if(objective==Spin2014Objective::StaticLowering || objective==Spin2014Objective::TopLoweringOnly)expected+=.7*std::pow(raw(3)/raw(0),2);
        if(objective==Spin2014Objective::TopMassAndLowering)expected+=.7*raw(3)*raw(3);
        require(std::abs(e.paperObjective-expected)<1e-12 && std::abs(e.regularizationEnergy-.016)<1e-12,"all five hand-derived paper formulas");
        for(double step:{1e-4,1e-5,1e-6})for(int i=0;i<e.gradient.size();++i) {
            auto plus=b,minus=b;double pp=phi,pm=phi;if(i<2){plus(i)+=step;minus(i)-=step;}else{pp+=step;pm-=step;}
            auto ep=f.evaluatePaper(plus,pp),em=f.evaluatePaper(minus,pm);
            require(std::abs((ep.objective-em.objective)/(2*step)-e.gradient(i))<2e-7,"three-step functional finite differences");
            require(((ep.constraints-em.constraints)/(2*step)-e.constraintJacobian.col(i)).norm()<1e-8,"constraint finite differences");
        }
    }
}
void solverTests() {
    auto r=fixedRequest();auto a=solveSourceAdaptiveSpin(r);
    std::cout<<std::setprecision(17)<<"FIXED SYNTHETIC_REGRESSION outcome="<<sourceAdaptiveSpinStatusToken(a.status)<<" passes="<<a.passes.size()<<" beta="<<a.relaxedFill.transpose()<<" input_hash="<<a.inputHash<<" trace_hash="<<a.traceHash<<" output_hash="<<a.outputHash<<"\n";
    for(const auto& p:a.passes) {std::cout<<"pass="<<p.index<<" before="<<p.before.size()<<" after="<<p.after.size()<<" solve="<<slqpStatusToken(p.solve.status)<<" beta="<<p.solve.beta.transpose()<<" phi="<<p.solve.finalPhi<<" before_hash="<<p.beforeHash<<" after_hash="<<p.afterHash<<" slqp_hash="<<p.solve.traceHash<<"\n";for(const auto& t:p.solve.trace){std::cout<<"{\"iteration\":"<<t.iteration<<",\"objective\":"<<t.objective<<",\"accepted_step\":"<<t.acceptedStep;writeSlqpAuditJson(std::cout,t);std::cout<<"}\n";}}
    require(a.status==SourceAdaptiveSpinStatus::BinaryValidationFailed,"STOP: fixed actual trajectory differs from required binary failure; inspect complete trace");
    require(a.passes.size()>=4 && a.passes.front().solve.trace.size()>1,"actual accepted step and source refinement");
    require((a.passes.front().solve.beta.array()-.25).abs().maxCoeff()<1e-6,"fixed first solve beta quarter");
    require(a.binaryCheck && a.binaryCheck->normalizedMomentError<=1e-10 && a.binaryCheck->normalizedConstraintResidual>r.binaryConstraintTolerance,"rounded failure remains failed with independent moments");
    require(std::abs(a.binaryCheck->normalizedConstraintResidual-1./131072)<1e-15,"fixed root-degree scaled binary residual");
    require(std::abs(a.binaryCheck->rawConstraintResidual-1./8192)<1e-15,"fixed raw first-x residual");
    std::cout<<"fixed_binary raw_residual="<<a.binaryCheck->rawConstraintResidual<<" normalized_residual="<<a.binaryCheck->normalizedConstraintResidual<<" moment_error="<<a.binaryCheck->normalizedMomentError<<" cuboids="<<a.binaryCuboids->cuboids.size()<<'\n';
    auto repeat=solveSourceAdaptiveSpin(r);require(a.traceHash==repeat.traceHash && a.outputHash==repeat.outputHash,"repeat trace and binary identity");
    for(std::size_t i=1;i<a.passes.size();++i)require(a.passes[i].inputPhi==a.passes[i-1].solve.finalPhi && a.passes[i].solve.trace.front().historyBefore==0,"phi carried and history resets");
    auto sym=fixedRequest(false);auto good=solveSourceAdaptiveSpin(sym);
    require(good.status==SourceAdaptiveSpinStatus::BinaryMomentValidated && good.relaxedFill.cwiseAbs().maxCoeff()<1e-6,"symmetric outer with positive x carving forces beta zero");
    require(good.passes.size()>=2 && good.passes.back().before.size()==1,"post recursive merge actual solve");
    auto limited=r;limited.maximumAdaptivePasses=1;auto pass=solveSourceAdaptiveSpin(limited);require(pass.status==SourceAdaptiveSpinStatus::AdaptivePassLimit && !pass.binaryCheck,"pass exhaustion UNKNOWN");
    limited=r;limited.maximumLeafCount=2;auto leaf=solveSourceAdaptiveSpin(limited);require(leaf.status==SourceAdaptiveSpinStatus::ResourceLimitReached && leaf.finalPartition.size()==2 && !leaf.binaryCheck,"leaf limit last valid state");
    auto hash=computeSourceAdaptiveSpinRequestHash(r);limited=r;limited.density=2;require(computeSourceAdaptiveSpinRequestHash(limited)!=hash,"density hash");limited=r;limited.objective=Spin2014Objective::TopLoweringOnly;require(computeSourceAdaptiveSpinRequestHash(limited)!=hash,"objective hash");limited=r;limited.maximumOutputCuboids+=1;require(computeSourceAdaptiveSpinRequestHash(limited)!=hash,"guard hash");limited=r;limited.binaryConstraintTolerance*=2;require(computeSourceAdaptiveSpinRequestHash(limited)!=hash,"tolerance hash");
    const std::vector<std::function<void(SourceAdaptiveSpinRequest&)>> mutations={
        [](auto& q){q.initialLevel=3;},[](auto& q){q.maximumOptimizationLevel=4;},[](auto& q){++q.maximumAdaptivePasses;},[](auto& q){++q.maximumLeafCount;},
        [](auto& q){q.initialFill=.4;},[](auto& q){q.fractionalEpsilon=.1;},[](auto& q){q.laplacianWeight=.7;},[](auto& q){q.inertiaWeight=.8;},[](auto& q){q.centerWeight=.9;},
        [](auto& q){++q.slqp.maximumIterations;},[](auto& q){++q.slqp.maximumLineSearchSteps;},[](auto& q){++q.slqp.lbfgsMemory;},
        [](auto& q){q.slqp.equalityTolerance*=2;},[](auto& q){q.slqp.projectedGradientTolerance*=2;},[](auto& q){q.slqp.activeSetTolerance*=2;},
        [](auto& q){q.slqp.meritPenalty*=2;},[](auto& q){q.slqp.lineSearchReduction=.4;},[](auto& q){q.slqp.armijo*=2;},
        [](auto& q){q.slqp.audit.restorationTrustBound*=2;},[](auto& q){++q.slqp.audit.lpIterationLimit;},[](auto& q){q.slqp.audit.lpTimeLimit*=2;},[](auto& q){q.slqp.audit.lpFeasibilityTolerance*=2;}};
    for(const auto& mutate:mutations){limited=r;mutate(limited);require(computeSourceAdaptiveSpinRequestHash(limited)!=hash,"every implementation control binds request hash");}
    limited=r;limited.slqp.algorithm=SlqpAlgorithm::LegacyV3;rejects([&]{solveSourceAdaptiveSpin(limited);},"legacy algorithm rejected");
    limited=r;limited.initialPhi=.1;rejects([&]{solveSourceAdaptiveSpin(limited);},"static phi rejected");
    limited=r;limited.slqp.maximumIterations=1;auto unknown=solveSourceAdaptiveSpin(limited);require(unknown.status==SourceAdaptiveSpinStatus::RelaxationUnknown && !unknown.binaryCheck,"nonconvergence UNKNOWN not infeasible");
    SourceAdaptiveSpinRequest angle;angle.source=source(cube(-.5,.5),cube(-.25,.25),5,{-1,-1,-1},2);angle.objective=Spin2014Objective::YoYo;angle.initialPhi=.37;angle.initialLevel=2;angle.maximumOptimizationLevel=3;
    auto carry=solveSourceAdaptiveSpin(angle);require(carry.status==SourceAdaptiveSpinStatus::BinaryMomentValidated && carry.passes.size()==2 && carry.finalPhi==.37,"nonzero phi carried across actual solves");
    for(const auto& p:carry.passes)require(p.inputPhi==.37 && p.solve.trace.front().historyBefore==0,"nonzero phi history reset");
    auto changedPhi=angle;changedPhi.initialPhi=.38;require(computeSourceAdaptiveSpinRequestHash(angle)!=computeSourceAdaptiveSpinRequestHash(changedPhi),"carried input phi binds hash");
    // Every removable point lies above current COM. Thus removal decreases the
    // static objective; symmetric full removal is the known global box optimum.
    auto top=r;top.source=source(box({-.5,-.5,0},{.5,.5,1}),box({-.25,-.25,.625},{.25,.25,.875}),5,{-1,-1,0},2);
    auto carved=solveSourceAdaptiveSpin(top);require(carved.status==SourceAdaptiveSpinStatus::BinaryMomentValidated && carved.binaryCuboids->cuboids.size()>1 && (carved.binaryFill.array()==1).all(),"known symmetric upper carve optimum");
    auto capped=top;capped.maximumOutputCuboids=1;auto fail=solveSourceAdaptiveSpin(capped);require(fail.status==SourceAdaptiveSpinStatus::ResourceLimitReached && !fail.binaryCuboids && !fail.binaryCheck && fail.finalPartitionHash==carved.finalPartitionHash && fail.relaxedFill==carved.relaxedFill,"output cap preserves valid continuous state and no partial artifact");
}
void binaryTests() {
    auto r=fixedRequest();SourceAdaptiveDomain d(r.source,2,5);Eigen::VectorXd binary=Eigen::VectorXd::Ones(d.leaves().size());
    auto output=sourceBinaryCuboids(d,binary,1000);const auto identity=d.partitionHash();
    const Ten outer=boxMoments({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    const Ten removed=boxMoments({3./16,-1./16,11./16},{5./16,1./16,13./16});
    auto expected=r.source->outerSurfaceMoments()-d.removedMoments(binary);
    auto check=checkSourceBinaryMoments(*r.source,output,expected,r.objective,0,1,1e-6);
    require(check.independentMoments.has_value() && check.normalizedMomentError<1e-10,"independent tetra cuboid moments");close(check.independentMoments->values(),outer-removed,"hand-derived selected intersections only");
    for(const auto& box:output.cuboids)require(box.lower.x()>=3./16 && box.upper.x()<=5./16 && box.lower.y()>=-1./16 && box.upper.y()<=1./16 && box.lower.z()>=11./16 && box.upper.z()<=13./16,"do not delete bounding-cell volume beyond accepted geometry");
    rejects([&]{sourceBinaryCuboids(d,binary,1);},"direct output cap");require(d.partitionHash()==identity && (binary.array()==1).all(),"output cap atomicity");
    auto outerShell=join(cube(0,1),cube(.375,.625),true);
    auto shell=source(outerShell,join(cube(.125,.875),cube(.375,.625),true));SourceAdaptiveDomain cavity(shell,2,5);
    Eigen::VectorXd all=Eigen::VectorXd::Ones(cavity.leaves().size());auto cubes=sourceBinaryCuboids(cavity,all,10000);
    const Ten exactOuter=boxMoments({0,0,0},{1,1,1})-boxMoments({.375,.375,.375},{.625,.625,.625});
    const Ten accepted=boxMoments(Eigen::Vector3d::Constant(5./32),Eigen::Vector3d::Constant(27./32))-boxMoments(Eigen::Vector3d::Constant(11./32),Eigen::Vector3d::Constant(21./32));
    AdaptiveFillFunctional functional(cavity,SpinIt2014MomentProblem::staticBalance(),1,2);
    close(functional.evaluatePaper(Eigen::VectorXd::Zero(all.size()),0).finalMoments.values(),2*exactOuter,"signed original cavity base");
    const auto evaluation=functional.evaluatePaper(all,0);close(evaluation.finalMoments.values(),2*(exactOuter-accepted),"signed conservative cavity analytic moments");
    auto checked=checkSourceBinaryMoments(*shell,cubes,evaluation.finalMoments,Spin2014Objective::StaticLowering,0,2,1e-6);require(checked.normalizedMomentError<1e-10,"signed cavity independent tetra moment gate");close(checked.independentMoments->values(),2*(exactOuter-accepted),"streamed signed outer and accepted cavity boxes");
    auto wrong=checkSourceBinaryMoments(*shell,cubes,shell->outerSurfaceMoments(),Spin2014Objective::StaticLowering,0,2,1e-6);require(!wrong.passed && wrong.normalizedMomentError>1e-10,"independent checker rejects incorrect expected vector");
}
void numericalFailureTests() {
    SourceAdaptiveSpinRequest r;r.source=source(cube(0,1),cube(.39,.40),3);r.initialLevel=2;r.maximumOptimizationLevel=3;
    auto empty=solveSourceAdaptiveSpin(r);require(empty.status==SourceAdaptiveSpinStatus::RelaxationUnknown && !empty.binaryCheck && empty.relaxedFill.size()==0 && empty.passes.size()==1,"zero-variable audited solve is explicit UNKNOWN without fake KKT");
    r=fixedRequest();r.objective=Spin2014Objective::YoYo;r.inertiaWeight=std::numeric_limits<double>::max();
    auto overflow=solveSourceAdaptiveSpin(r);require(overflow.status==SourceAdaptiveSpinStatus::RelaxationUnknown && !overflow.binaryCheck && overflow.finalPartition.size()==2 && (overflow.relaxedFill.array()==.5).all(),"nonfinite objective retains initial state and explicit UNKNOWN");
}
int main(int argc,char**) {try {
    if(argc>1){auto s=source(cube(0,1),cube(.25,.75),9);SourceAdaptiveDomain d(s,2,9);require(s->stats().maximumVisitedLevel==9 && d.leaves().size()==8,"genuine geometry9 coarse partition");graph(d);rusage u{};getrusage(RUSAGE_SELF,&u);std::cout<<"geometry_depth9 nodes="<<s->cells().size()<<" optimization_leaves="<<d.leaves().size()<<" face_keys="<<d.faceKeyCount()<<" rss_KiB="<<u.ru_maxrss<<'\n';return 0;}
    numericalFailureTests();partitionTests();functionalTests();binaryTests();solverTests();std::cout<<"source adaptive tests passed; faster="<<MAKE_IT_STAND_HAS_SPIN_IT_FASTER<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#else
#include <iostream>
int main(){std::cerr<<"RED: SourceAdaptiveSpinSolver source-native adaptive integration is unavailable\n";return 1;}
#endif
