// Preserved failed synthetic attempt1. Source function as compiled into
// task-5aa-green-test-1.log's executable. Used with the includes and helpers of
// tests/spin_it_2014_source_cage_component_clearance_regression.cpp.
// Failure text: actual full BBW succeeds. No terminal reason was exported in
// that run, so the observed approximately300s delay suggests but does not prove
// a BBW timeout. No downstream component result came from this wide fixture.
SourceCagePaperFunctional fixture(bool wide=false) {
    const auto outer=wide?cube(.25,.75):cube(.3125,.6875);
    const auto design=wide?box({.3625,.425,.425},{.6375,.575,.575}):cube(.425,.575);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,4,4,1000);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;pc.inertiaWeight=0.;pc.centerWeight=0.;
    return SourceCagePaperFunctional(SourceCageMomentFunctional(SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.17,mc),pc);
}
// Wide-case state used all full handles: scale1.02, translation(.011,-.007,.019),
// phi.23. Beta1 iff depth4 source leaf's x coordinate equals6 or8; else beta0.
// Intended independently specified component threshold.08, never reached.
