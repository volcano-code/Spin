#include "spin2014/SourceCageBinaryMaterial.h"
#include "source_surface_test_utils.h"
#include "task-5w-controller-export-helpers.h"

void meshJson(const spin2014::SourceCageBinaryMesh& mesh) {
    std::cout<<"{\"vertices\":";matrixJson(mesh.vertices);
    std::cout<<",\"faces\":";matrixJson(mesh.triangles);std::cout<<'}';
}

int main() {try {
    using namespace spin2014;
    std::cout<<std::setprecision(17)<<std::boolalpha;
    const auto outer=surface_test::cube(.3125,.6875),design=surface_test::cube(.425,.575);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    const SourceAdaptiveDomain domain(source,4,4);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    surface_test::require(bbw.success(),"full source BBW required");
    const SourceCageMomentFunctional material(SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{}),1.17,{});
    const char* names[]={"empty_affine","full_affine","tie_affine","mixed_rounding_failure"};
    for(int test=0;test<4;++test) {
        SourceCagePaperFunctionalConfig pc;
        pc.objective=test==3?Spin2014Objective::StaticLowering:Spin2014Objective::TopInertiaOnly;
        pc.inertiaWeight=0.;pc.centerWeight=0.;
        // KNOWN_SOLUTION_VALIDATION: zero objective, fixed feasible relaxed axis.
        SourceCagePaperFrame frame;
        const Eigen::Vector3d translation(.011,-.007,.019);
        frame.origin=1.02*Eigen::Vector3d(.5,.5,.3)+translation;
        const SourceCagePaperFunctional original(material,pc,frame);
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(original.betaCount(),test==0?.25:(test==1?.75:.5));
        if(test==3)for(std::size_t i=0;i<domain.leaves().size();++i)
            beta[i]=domain.leaves()[i].coordinates[0]<8?.4:.6;
        Eigen::VectorXd x=original.initialState(beta,test==3?0.:.23);
        const auto& cage=bbw.fem().grid().cage();
        for(std::size_t h=0;h<original.handleCount();++h) {
            const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
            const auto offset=original.handleOffset()+4*h;
            x[offset]=1.02;
            x.segment<3>(offset+1)=(1.02-1.)*center+translation;
        }
        if(test==3) {
            const auto m=original.evaluateValues(x).world.finalMoments.values();
            // Explicitly known-solution setup, NOT optimization or paper input.
            frame.origin.x()=m[1]/m[0];frame.origin.y()=m[2]/m[0];
        }
        const SourceCagePaperFunctional f(material,pc,frame);
        SourceCageAdaptiveSolverSettings settings;settings.maximumPasses=2;
        settings.solve.source.solver.maximumIterations=1;
        settings.solve.source.solver.maximumLineSearchSteps=2;
        settings.solve.budget.maximumGuardCalls=8;
        const auto adaptive=solveSourceCageAdaptive(f,x,settings);
        surface_test::require(adaptive.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete,
            "fixed stationary protocol must reach actual relaxed endpoint");
        const SourceCageBinaryMaterialSettings bs;
        const auto binary=validateSourceCageBinaryMaterial(adaptive,bs);
        const auto& facts=binary.facts();const auto candidate=binary.candidate();
        surface_test::require(bool(candidate),"actual emitted candidate required even for rounding failure");
        std::cout<<"{\"id\":"<<std::quoted(names[test])
            <<",\"scope\":\"KNOWN_SOLUTION_VALIDATION\",\"sealed\":false,\"paper_reproduction\":false"
            <<",\"adaptive_outcome\":\"RELAXED_ADAPTIVE_COMPLETE\",\"adaptive_hash\":"<<std::quoted(adaptive.facts().scientificSha256)
            <<",\"outcome\":"<<std::quoted(sourceCageBinaryMaterialOutcomeToken(facts.outcome))
            <<",\"scientific_hash\":"<<std::quoted(facts.scientificSha256)
            <<",\"validated\":"<<bool(binary.validatedMaterial())
            <<",\"empty_cavity\":"<<facts.emptyCavity
            <<",\"constraint_tolerance\":"<<bs.constraintTolerance
            <<",\"constraint_residual\":"<<facts.normalizedConstraints.cwiseAbs().maxCoeff()
            <<",\"geometry_passed\":"<<(facts.geometry&&facts.geometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed)
            <<",\"material_boundary_valid\":"<<(facts.exactBoundary&&facts.exactBoundary->facts().boundaryValid)
            <<",\"relaxed_state\":";arrayJson(candidate->relaxedState());
        std::cout<<",\"delta\":";arrayJson(facts.stateDelta);
        std::cout<<",\"world_integrated\":";arrayJson(facts.worldIntegratedMoments);
        std::cout<<",\"spin_integrated\":";arrayJson(facts.spinIntegratedMoments);
        std::cout<<",\"binary_functional_moments\":";arrayJson(facts.binaryFunctionalWorldMoments);
        std::cout<<",\"constraints\":";arrayJson(facts.rawConstraints);
        std::cout<<",\"normalized_constraints\":";arrayJson(facts.normalizedConstraints);
        std::cout<<",\"inertia_eigenvalues\":";arrayJson(facts.centralInertiaEigenvalues);
        std::cout<<",\"inertia_psd\":"<<facts.centralInertiaPsd;
        std::cout<<",\"central_inertia\":";matrixJson(facts.centralInertia);
        std::cout<<",\"origin\":";arrayJson(candidate->functional().frame().origin);
        std::cout<<",\"rotation\":";matrixJson(candidate->functional().frame().rotation);
        std::cout<<",\"outer\":";meshJson(candidate->outerBoundary());
        std::cout<<",\"cavity\":";meshJson(candidate->carvedBoundary());
        std::cout<<",\"material\":";meshJson(candidate->materialBoundary());
        std::cout<<",\"source\":";dumpSource(candidate->functional(),candidate->binaryState());
        std::cout<<"}\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
