#include "spin2014/SourceCageBbw.h"
#include "source_surface_test_utils.h"
#include <functional>
#include <iostream>
#include <iomanip>
using namespace spin2014;
using namespace surface_test;
using Protocol=SourceCageBbwProtocol;
SourceCageGrid makeGrid(const Mesh& mesh,Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1) {
    SourceSurfaceCageConfig c; c.cageLevel=2; c.rootOrigin=origin; c.rootSide=side;
    SourceCageGridConfig g; g.fineLevel=3;
    return {SourceSurfaceCage(mesh.v,mesh.f,c),g};
}
void near(double a,double b,double tolerance=1e-9) {
    require(std::isfinite(a) && std::abs(a-b)<=tolerance*std::max(1.,std::abs(b)),"literal or covariance mismatch");
}
std::vector<int> freeNodes(const SourceCageFem& f) {
    std::vector<bool> fixed(f.grid().nodeKeys().size(),false);
    for(const auto& b:f.boundaryConstraints()) fixed[b.nodeIndex]=true;
    std::vector<int> out;
    for(int n=0;n<int(fixed.size());++n) if(!fixed[n]) out.push_back(n);
    return out;
}
void validated(const SourceCageBbw& run,Protocol p) {
    if(!run.success()) {
        for(const auto& a:run.attempts()) {
            std::cerr<<a.diagnostics<<"ASSESS boundary="<<a.assessment.boundaryResidual
                <<" partition="<<a.assessment.partitionResidual<<" stationary="<<a.assessment.stationarityResidual
                <<" projected="<<a.assessment.projectedGradientResidual<<'\n';
        }
        throw std::runtime_error("actual solve rejected: "+run.terminalReason());
    }
    require(run.rawAssessment().available && run.rawAssessment().primalPassed && run.rawAssessment().optimalityPassed,"raw original FEM KKT");
    require(run.finalAssessment().available && run.finalAssessment().primalPassed,"final primal validity");
    if(p==Protocol::COUPLED_PARTITION) require(run.finalAssessment().optimalityPassed,"coupled final KKT");
    require(run.finalWeights().allFinite(),"finite actual final field");
    require(!run.attempts().empty(),"actual backend invocation required");
    for(const auto& a:run.attempts()) {
        require(a.runInvoked && a.candidateAvailable && a.assessment.optimalityPassed,"every backend candidate independently assessed");
        require(!a.diagnostics.empty(),"actual backend diagnostic text retained");
        require(a.qpIterations>=0 && a.endUnixSeconds>=a.startUnixSeconds && a.elapsedSeconds>=0,"actual trace metadata");
    }
}
int main() {
    int failures=0;
    auto test=[&](const char* name,const std::function<void()>& body) {
        try {body(); std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& e) {++failures; std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';}
    };
    std::cout<<std::setprecision(17);
    for(auto mass:{SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,SourceCageFemMassRule::BARYCENTRIC_TET}) {
        SourceCageFemConfig fc; fc.massRule=mass;
        for(auto protocol:{Protocol::COUPLED_PARTITION,Protocol::INDEPENDENT_NORMALIZED}) {
            SourceCageBbwConfig c; c.protocol=protocol;
            test("actual_small_analytic_solve_all_boundary_rows",[&] {
                const SourceCageFem f(makeGrid(cube(.3,.4)),fc);
                require(f.grid().nodeKeys().size()==27 && f.tetrahedra().size()==48 && f.boundaryConstraints().size()==26,"fixed analytic fixture");
                const SourceCageBbw run(f,c); validated(run,protocol);
                const int n=freeNodes(f).at(0);
                for(int h=0;h<8;++h) {
                    const auto& key=f.grid().nodeKeys()[f.grid().handleNodeIndices()[h]];
                    const bool diagonal=(key==SourceCageGrid::Key{2,2,2} || key==SourceCageGrid::Key{4,4,4});
                    near(run.rawWeights()(n,h),diagonal?5./16:1./16);
                    near(run.finalWeights()(n,h),diagonal?5./16:1./16);
                }
                // Independent cube-face barycentric formula: order the two
                // varying coordinates. Triangle weights are 1-max, |u-v|, min.
                for(int i=0;i<27;++i) if(i!=n) {
                    const auto key=f.grid().nodeKeys()[i];
                    int axis=0; while(key[axis]==3) ++axis;
                    int a=(axis+1)%3,b=(axis+2)%3;
                    for(int h=0;h<8;++h) {
                        const auto hk=f.grid().nodeKeys()[f.grid().handleNodeIndices()[h]];
                        double expected=0;
                        if(hk[axis]==key[axis]) {
                            const double u=(key[a]-2)/2.,v=(key[b]-2)/2.;
                            if(hk[a]==2 && hk[b]==2) expected=1-std::max(u,v);
                            if(hk[a]==4 && hk[b]==4) expected=std::min(u,v);
                            if(hk[a]==4 && hk[b]==2) expected=std::max(u-v,0.);
                            if(hk[a]==2 && hk[b]==4) expected=std::max(v-u,0.);
                        }
                        near(run.finalWeights()(i,h),expected);
                    }
                }
                require(run.attempts().size()==(protocol==Protocol::COUPLED_PARTITION?1:8),"actual complete handle calls");
                for(const auto& attempt:run.attempts()) require(attempt.qpIterations>0,"backend actually iterates from its own non-solution start");
            });
            test("actual_702_variable_and_26_column_box_qps",[&] {
                const SourceCageFem f(makeGrid(cube(.3,.7)),fc);
                require(f.grid().cells().size()==64 && f.grid().nodeKeys().size()==125 && freeNodes(f).size()==27 && f.grid().handleNodeIndices().size()==26,"fixed nontrivial fixture");
                const auto inputHash=f.outputSha256();
                const auto inputK=f.stiffness(); const auto inputMass=f.lumpedMass();
                const auto inputNodes=f.grid().nodePositions();
                const SourceCageBbw run(f,c),again(f,c);
                validated(run,protocol); validated(again,protocol);
                require(run.outputSha256()==again.outputSha256() && run.rawWeights()==again.rawWeights(),"two actual solves deterministic");
                require(f.outputSha256()==inputHash,"immutable FEM input");
                require((f.stiffness()-inputK).norm()==0 && f.lumpedMass()==inputMass && f.grid().nodePositions()==inputNodes,"actual source/FEM arrays unchanged");
                require(run.attempts().size()==(protocol==Protocol::COUPLED_PARTITION?1:26),"no omitted column QPs");
                int active=0;
                for(int n:freeNodes(f)) for(int h=0;h<26;++h) if(std::abs(run.rawWeights()(n,h))<1e-9) ++active;
                require(active>0,"nontrivial actual active bounds retained");
                if(protocol==Protocol::INDEPENDENT_NORMALIZED) {
                    require(run.rawAssessment().partitionResidual>1e-4 && run.maximumNormalizationChange()>1e-4,"raw independent partition differs; normalization recorded");
                    require(!run.finalAssessment().optimalityPassed && run.finalAssessment().projectedGradientResidual>1e-4,"raw optimality is not transferred to normalized final field");
                    SourceCageBbwConfig coupledConfig;
                    require(!assessSourceCageBbw(f,run.rawWeights(),coupledConfig).primalPassed,"raw independent partition violation rejected by coupled assessment");
                    require(std::abs(run.finalAssessment().energy-run.rawAssessment().energy)>1e-4,"normalization energy difference retained");
                }
                Eigen::MatrixXd points(9,3); points.topRows(8)=cube(.3,.7).v; points.row(8)<<.4375,.5625,.4375;
                const auto transferred=f.grid().interpolate(run.finalWeights(),points);
                require(transferred.allFinite() && transferred.minCoeff()>=-1e-9 && transferred.maxCoeff()<=1+1e-9,"surface/new dyadic vertex transfer bounds");
                near((transferred.rowwise().sum().array()-1).abs().maxCoeff(),0);
                std::cout<<"FACT mass="<<int(mass)<<" protocol="<<int(protocol)<<" rawEnergy="<<run.rawAssessment().energy<<" finalEnergy="<<run.finalAssessment().energy
                    <<" rawPartition="<<run.rawAssessment().partitionResidual<<" normalizationChange="<<run.maximumNormalizationChange()
                    <<" rawKkt="<<run.rawAssessment().projectedGradientResidual<<" finalCoupledKkt="<<run.finalAssessment().projectedGradientResidual<<" active="<<active<<" hash="<<run.outputSha256()<<'\n';
            });
            test("disconnected_actual_solve_preserves_zero_influence",[&] {
                const SourceCageFem f(makeGrid(join(cube(.03,.09),cube(.91,.97))),fc);
                const SourceCageBbw run(f,c); validated(run,protocol);
                require(f.grid().handleNodeIndices().size()==16 && freeNodes(f).size()==2,"both disconnected domains retained");
                for(int n=0;n<run.rawWeights().rows();++n) for(int h=0;h<16;++h) {
                    const bool nodeLeft=f.grid().nodePositions()(n,0)<.5;
                    const bool handleLeft=f.grid().nodePositions()(f.grid().handleNodeIndices()[h],0)<.5;
                    if(nodeLeft!=handleLeft) near(run.rawWeights()(n,h),0);
                }
                auto limited=c; limited.maximumQPairContributions=1;
                const SourceCageBbw failed(f,limited);
                require(!failed.success() && failed.terminal()==SourceCageBbwTerminal::RESOURCE_LIMIT && failed.finalWeights().size()==0,"resource failure never validated");
            });
            test("translation_scale_weights_and_unscaled_energy",[&] {
                const SourceCageFem base(makeGrid(cube(.3,.7)),fc); const SourceCageBbw a(base,c); validated(a,protocol);
                for(double s:{2.,1.7}) {
                    const Eigen::Vector3d origin=s==2?Eigen::Vector3d(-.5,-.25,-.125):Eigen::Vector3d(-.3,-.2,-.1);
                    auto mesh=cube(.3,.7); mesh.v=(s*mesh.v).rowwise()+origin.transpose();
                    const SourceCageFem shifted(makeGrid(mesh,origin,s),fc);
                    const SourceCageBbw b(shifted,c),repeat(shifted,c);
                    validated(b,protocol); validated(repeat,protocol);
                    require(b.outputSha256()==repeat.outputSha256() && b.finalWeights()==repeat.finalWeights(),"transformed actual solves deterministic");
                    const double difference=(a.finalWeights()-b.finalWeights()).cwiseAbs().maxCoeff();
                    near(difference,0);
                    near(b.rawAssessment().energy,a.rawAssessment().energy/s);
                    std::cout<<"COVARIANCE mass="<<int(mass)<<" protocol="<<int(protocol)<<" scale="<<s
                             <<" weightDifference="<<difference<<" rawKkt="<<b.rawAssessment().projectedGradientResidual
                             <<" energy="<<b.rawAssessment().energy<<" hash="<<b.outputSha256()<<'\n';
                }
            });
        }
    }
    test("independent_assessor_rejects_false_acceptance",[&] {
        const SourceCageFem f(makeGrid(cube(.3,.4)),{}); const SourceCageBbw run(f); validated(run,Protocol::COUPLED_PARTITION);
        const int n=freeNodes(f)[0]; const auto good=run.finalWeights();
        auto bad=good; bad.row(n).setConstant(1./8);
        auto a=assessSourceCageBbw(f,bad,{}); require(a.primalPassed && !a.optimalityPassed,"feasible nonoptimal interior rejected");
        SourceCageBbwConfig independentConfig; independentConfig.protocol=Protocol::INDEPENDENT_NORMALIZED;
        a=assessSourceCageBbw(f,bad,independentConfig);
        require(a.primalPassed && !a.optimalityPassed,"independent feasible nonoptimal interior rejected");
        bad=good; bad.row(n).setZero(); bad(n,0)=1;
        a=assessSourceCageBbw(f,bad,{}); require(a.primalPassed && !a.optimalityPassed && a.projectedGradientResidual>1e-8,"wrong pure active candidate rejected");
        for(int which=0;which<6;++which) {
            bad=good;
            if(which==0) bad(0,0)+=.1;
            if(which==1) bad.row(n)*=.5;
            if(which==2) {bad(n,0)=-.1; bad(n,1)+=.1+good(n,0);}
            if(which==3) bad(n,0)=std::numeric_limits<double>::quiet_NaN();
            if(which==4) bad(n,0)=std::numeric_limits<double>::infinity();
            if(which==5) bad.conservativeResize(26,8);
            a=assessSourceCageBbw(f,bad,{}); require(!a.primalPassed,"bad primal candidate rejected");
        }
    });
    test("iteration_time_resource_and_invalid_configuration_fail_closed",[&] {
        const SourceCageFem f(makeGrid(cube(.3,.7)),{});
        for(auto p:{Protocol::COUPLED_PARTITION,Protocol::INDEPENDENT_NORMALIZED}) {
            SourceCageBbwConfig c; c.protocol=p; c.iterationLimit=0;
            const SourceCageBbw run(f,c);
            require(!run.success() && run.terminal()==SourceCageBbwTerminal::ITERATION_LIMIT && !run.attempts().empty() && run.finalWeights().size()==0,"actual zero-iteration limit is not convergence");
            c.iterationLimit=100000; c.timeLimitSeconds=1e-12;
            const SourceCageBbw timed(f,c);
            require(!timed.success() && timed.terminal()==SourceCageBbwTerminal::TIME_LIMIT && timed.finalWeights().size()==0,"actual time exhaustion fails closed");
        }
        for(int which=0;which<5;++which) for(auto limit:{std::size_t(1),std::numeric_limits<std::size_t>::max()}) {
            SourceCageBbwConfig c;
            if(which==0)c.maximumQPairContributions=limit;
            if(which==1)c.maximumQNonzeros=limit;
            if(which==2)c.maximumVariables=limit;
            if(which==3)c.maximumDenseValues=limit;
            if(which==4)c.maximumBytes=limit;
            const SourceCageBbw run(f,c);
            require(!run.success() && run.terminal()==SourceCageBbwTerminal::RESOURCE_LIMIT,"preflight resource/index/byte overflow rejected");
        }
        for(int which=0;which<5;++which) {
            SourceCageBbwConfig c;
            if(which==0)c.primalTolerance=std::numeric_limits<double>::quiet_NaN();
            if(which==1)c.optimalityTolerance=-1;
            if(which==2)c.timeLimitSeconds=std::numeric_limits<double>::infinity();
            if(which==3)c.iterationLimit=-1;
            if(which==4)c.protocol=static_cast<Protocol>(123);
            const SourceCageBbw run(f,c);
            require(!run.success() && run.terminal()==SourceCageBbwTerminal::INVALID_CONFIG,"invalid configs cannot grant status");
        }
    });
    std::cout<<"source cage BBW SYNTHETIC_REGRESSION: "<<failures<<" failed groups\n";
    return failures?1:0;
}
