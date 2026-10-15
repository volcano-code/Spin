#include "spin2014/SourceCageGuardedSlqp.h"
#include "source_surface_test_utils.h"

#include <Eigen/LU>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace spin2014;
using namespace surface_test;

namespace {
SourceConformingBoundary asymmetricBoundary() {
    // Frozen, deliberately offset and anisotropic geometry. These dimensions
    // are a method fixture, not parameters attributed to Spin-It.
    const auto outer=box({.15,.18,.12},{.85,.82,.88});
    const auto design=box({.31,.26,.29},{.60,.57,.67});
    SourceSurfaceOctreeConfig sourceConfig;
    sourceConfig.maximumDepth=3;
    sourceConfig.rootOrigin=Eigen::Vector3d::Zero();
    sourceConfig.rootSide=1.;
    auto source=std::make_shared<const SourceSurfaceOctree>(
        outer.v,outer.f,design.v,design.f,sourceConfig);
    SourceAdaptiveDomain domain(source,2,3);

    SourceSurfaceCageConfig cageConfig;
    cageConfig.cageLevel=2;
    cageConfig.rootOrigin=Eigen::Vector3d::Zero();
    cageConfig.rootSide=1.;
    SourceCageGridConfig gridConfig;
    gridConfig.fineLevel=3;
    // Match the vendored libigl BBW protocol: solve each handle's bounded
    // biharmonic QP independently, then normalize rows. The alternative
    // coupled-partition protocol is retained as a separate diagnostic result.
    SourceCageBbwConfig bbwConfig;
    bbwConfig.protocol=SourceCageBbwProtocol::INDEPENDENT_NORMALIZED;
    const SourceCageFem fem(SourceCageGrid(
        SourceSurfaceCage(outer.v,outer.f,cageConfig),gridConfig),{});
    SourceCageBbw bbw(fem,bbwConfig);
    require(bbw.success(),"asymmetric fixture solves actual full-handle BBW");
    std::size_t refinedHandles=0;
    double worstBackendKkt=0.,worstRefinedKkt=0.;
    for(const auto& attempt:bbw.attempts()) {
        require(attempt.assessment.optimalityPassed,"every independent handle meets original FEM KKT gate");
        worstBackendKkt=std::max(worstBackendKkt,attempt.backendStationarityResidual);
        worstRefinedKkt=std::max(worstRefinedKkt,attempt.refinedStationarityResidual);
        if(attempt.activeSetRefinementApplied) {
            ++refinedHandles;
            require(attempt.backendStationarityResidual>bbwConfig.optimalityTolerance,
                    "refinement only repairs a candidate that originally misses KKT");
            require(attempt.refinedStationarityResidual<=bbwConfig.optimalityTolerance,
                    "refined candidate meets the unchanged KKT threshold");
        }
    }
    require(refinedHandles>0,"asymmetric fixture exercises active-set numerical refinement");
    std::cout<<"[BBW] independent_handles="<<bbw.attempts().size()
             <<" polished_handles="<<refinedHandles
             <<" worst_backend_stationarity="<<worstBackendKkt
             <<" worst_post_polish_stationarity="<<worstRefinedKkt
             <<" output_digest="<<bbw.outputSha256()<<'\n';
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}

Ten orientedTetraMoments(const Eigen::MatrixXd& points,
                         SourceConformingBoundary::Triangle face) {
    Eigen::Matrix3d p;
    for(int i=0;i<3;++i) p.row(i)=points.row(static_cast<Eigen::Index>(face[i]));
    const double volume=p.determinant()/6.;
    const Eigen::Vector3d sum=p.colwise().sum().transpose();
    const Eigen::Matrix3d q=volume*(sum*sum.transpose()+p.transpose()*p)/20.;
    Ten result;
    result<<volume,volume*sum/4.,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);
    return result;
}

Ten independentlyIntegrate(const SourceConformingBoundary& boundary,
                           const SourceCagePaperEvaluation& evaluation,
                           const Eigen::VectorXd& beta,double density) {
    Ten result=Ten::Zero();
    const auto& points=evaluation.world.deformedPoints;
    const auto& outerFaces=boundary.binding().source().outerFaces();
    for(Eigen::Index i=0;i<outerFaces.rows();++i)
        result+=orientedTetraMoments(points,{{
            static_cast<std::size_t>(outerFaces(i,0)),
            static_cast<std::size_t>(outerFaces(i,1)),
            static_cast<std::size_t>(outerFaces(i,2))}});
    for(std::size_t owner=0;owner<boundary.primitiveBoundaries().size();++owner)
        for(const auto& ref:boundary.primitiveBoundaries()[owner])
            result-=beta[boundary.binding().primitives()[owner].betaIndex]
                *static_cast<double>(ref.orientation)
                *orientedTetraMoments(points,boundary.triangles()[ref.triangle]);
    return density*result;
}
}

#include "../experiments/spin_it_reproduction/strict/source_cage_joint_batch.h"

int main(int argc,char** argv) {
    try {
        if(argc==3&&(std::string(argv[1])=="--batch-preflight"||std::string(argv[1])=="--batch-run"))
            return joint_batch::run(argv[2],std::string(argv[1])=="--batch-preflight");
        if(argc!=1)throw std::invalid_argument("expected --batch-preflight DIR or --batch-run DIR");
        const auto boundary=asymmetricBoundary();
        const SourceCageMomentFunctional moments(boundary,1.3,{});
        require(moments.handleCount()==98,"all level-2 cage boundary handles are retained");

        SourceCagePaperFunctionalConfig config;
        config.objective=Spin2014Objective::TopMassAndLowering;
        config.inertiaWeight=.8;
        config.centerWeight=.6;
        config.fillWeight=.03;
        config.shapeWeight=.04;
        config.displacementWeight=.05;
        SourceCagePaperFrame frame;
        frame.origin=Eigen::Vector3d(.5,.5,.5);
        const SourceCagePaperFunctional functional(moments,config,frame);

        Eigen::VectorXd beta=Eigen::VectorXd::LinSpaced(
            static_cast<Eigen::Index>(functional.betaCount()),.2,.65);
        const Eigen::VectorXd initial=functional.initialState(beta,.17);
        const auto initialEvaluation=functional.evaluate(initial);
        require(initialEvaluation.hasJacobians,"actual joint derivatives are available");
        require(initialEvaluation.constraints.norm()>1e-5,
                "the frozen asymmetric input begins outside the target constraints");

        SourceCageGuardedSlqpSettings settings;
        settings.source.solver.maximumIterations=4;
        settings.source.solver.maximumLineSearchSteps=10;
        settings.geometry.minimumOuterDomainClearance=0.;
        const auto audit=solveSourceCageGuardedSlqp(functional,initial,settings);
        const auto& facts=audit.facts();
        require(facts.source.numericalResultAvailable,"guarded SLQP returns an actual numerical result");
        require(facts.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed,
                "final non-symmetric joint state passes fresh exact geometry checks");
        require(facts.finalRecheckCount==1,"final accepted state receives one fresh geometry recheck");

        const Eigen::VectorXd& final=facts.source.numerical.finalState;
        require(final.size()==static_cast<Eigen::Index>(functional.stateCount()),
                "complete beta/angle/all-handle state is retained");
        const double fillChange=(final.head(static_cast<Eigen::Index>(functional.betaCount()))-beta).norm();
        const double handleChange=(final.segment(static_cast<Eigen::Index>(functional.handleOffset()),
            static_cast<Eigen::Index>(functional.handleCount()*4))-initial.segment(
            static_cast<Eigen::Index>(functional.handleOffset()),
            static_cast<Eigen::Index>(functional.handleCount()*4))).norm();
        require(fillChange>1e-9,"joint solve changes internal fill variables");
        require(handleChange>1e-9,"joint solve changes cage handle parameters");

        const auto finalEvaluation=functional.evaluateValues(final);
        const Ten independent=independentlyIntegrate(boundary,finalEvaluation,
            final.head(static_cast<Eigen::Index>(functional.betaCount())),moments.density());
        const Ten native=finalEvaluation.world.finalMoments.values();
        const double momentError=(independent-native).cwiseAbs().maxCoeff()
            /std::max(1.,native.cwiseAbs().maxCoeff());
        require(momentError<=1e-10,"independent signed-tetrahedron moments match the final deformed material");
        require(finalEvaluation.shapeEnergy>0.,"accepted deformation has nonzero shape regularization energy");

        std::cout<<"[PASS] reconstructed asymmetric source-cage joint case"
            <<" beta_count="<<functional.betaCount()
            <<" handle_count="<<functional.handleCount()
            <<" iterations="<<facts.source.numerical.trace.size()
            <<" initial_constraint_norm="<<initialEvaluation.constraints.norm()
            <<" final_constraint_norm="<<finalEvaluation.constraints.norm()
            <<" fill_change="<<fillChange
            <<" handle_change="<<handleChange
            <<" independent_moment_error="<<momentError
            <<" source_digest="<<facts.scientificSha256<<'\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"[FAIL] reconstructed asymmetric source-cage joint case: "
                 <<error.what()<<'\n';
        return 1;
    }
}
