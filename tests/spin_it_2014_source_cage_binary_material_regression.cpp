#include <iostream>
#include <algorithm>
#if !__has_include("spin2014/SourceCageBinaryMaterial.h")
int main(){std::cerr<<"FAIL: actual deformed binary material admission and validation layer is absent\n";return 1;}
#else
#include "spin2014/SourceCageBinaryMaterial.h"
#include "spin2014/SlqpIterationStore.h"
#include "source_surface_test_utils.h"
#include "manufacturing/SlicerValidation.h"
#include <Eigen/Eigenvalues>
#include <array>
#include <chrono>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <system_error>
#include <type_traits>
#include <unistd.h>
using namespace spin2014;
using namespace surface_test;
namespace {
bool identical(const Eigen::VectorXd& a,const Eigen::VectorXd& b) {
    return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;
}
bool variableBlockChangedWithinPass(const Eigen::VectorXd& before,
    const Eigen::VectorXd& after,std::size_t offset,std::size_t count) {
    if(before.size()!=after.size()||offset>std::size_t(before.size())||
       count>std::size_t(before.size())-offset)return false;
    return (after.segment(static_cast<Eigen::Index>(offset),
        static_cast<Eigen::Index>(count))-
        before.segment(static_cast<Eigen::Index>(offset),
        static_cast<Eigen::Index>(count))).norm()>1e-9;
}
bool acceptedTraceBlockChanged(const Eigen::VectorXd& initial,
    const std::vector<SlqpIteration>& trace,std::size_t offset,std::size_t count) {
    // A SLQP row stores the state before its line search. Its acceptedStep
    // applies to the transition to the next row, so compare adjacent states
    // only when that transition was accepted. This remains valid when an
    // adaptive partition changes variable dimensions between solve passes.
    if(trace.empty()||!initial.isApprox(trace.front().state,0.))return false;
    for(std::size_t i=0;i+1<trace.size();++i) {
        const auto& from=trace[i];
        const auto& to=trace[i+1];
        if(from.acceptedStep<=0.)continue;
        if(from.state.size()!=initial.size()||to.state.size()!=initial.size())continue;
        if(variableBlockChangedWithinPass(from.state,to.state,offset,count))return true;
    }
    return false;
}
int jointIndexGuardTest() {
    Eigen::VectorXd before=Eigen::VectorXd::Zero(8+8);
    Eigen::VectorXd afterSplit=Eigen::VectorXd::Zero(12+8);
    afterSplit.head(12).setConstant(.25);
    afterSplit.segment(12,8).setConstant(.5);
    if(variableBlockChangedWithinPass(before,afterSplit,8,8))return 1;
    // The endpoint handle offset is 12, but comparing its block to the
    // pre-split full state is invalid: the two states have different layouts.
    // Joint-update success must come from accepted same-partition traces.
    if(variableBlockChangedWithinPass(before,afterSplit,12,8))return 1;
    Eigen::VectorXd passBefore=Eigen::VectorXd::Zero(12+8);
    Eigen::VectorXd passAfter=passBefore;
    passAfter.head(12).setConstant(.5);
    if(!variableBlockChangedWithinPass(passBefore,passAfter,0,12)||
       variableBlockChangedWithinPass(passBefore,passAfter,12,8))return 1;
    passAfter[12]=1.25;
    if(!variableBlockChangedWithinPass(passBefore,passAfter,12,8))return 1;
    SlqpIteration initialRow;initialRow.state=passBefore;initialRow.acceptedStep=.25;
    SlqpIteration acceptedIntermediate;acceptedIntermediate.state=passBefore;
    acceptedIntermediate.state[12]=1.25;acceptedIntermediate.acceptedStep=.1;
    SlqpIteration returnedEndpoint;returnedEndpoint.state=passBefore;
    const std::vector<SlqpIteration> acceptedThenReturned={initialRow,acceptedIntermediate,returnedEndpoint};
    if(!acceptedTraceBlockChanged(passBefore,acceptedThenReturned,12,8))return 1;
    SlqpIteration rejectedRow;rejectedRow.state=passBefore;rejectedRow.acceptedStep=0.;
    SlqpIteration rejectedEndpoint;rejectedEndpoint.state=passBefore;
    rejectedEndpoint.state[12]=1.25;
    if(acceptedTraceBlockChanged(passBefore,std::vector<SlqpIteration>{rejectedRow,rejectedEndpoint},12,8))return 1;
    Eigen::VectorXd unrelatedInitial=passBefore;unrelatedInitial[12]=-.5;
    if(acceptedTraceBlockChanged(unrelatedInitial,acceptedThenReturned,12,8))return 1;
    return 0;
}
class DurableProgressLog {
public:
    explicit DurableProgressLog(const std::string& path) {
        fd_=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0600);
        if(fd_<0)throw std::system_error(errno,std::generic_category(),"create joint progress log");
        append("{\"schema\":\"spin-it-joint-progress-v1\",\"record_type\":\"header\"}\n");
    }
    ~DurableProgressLog(){if(fd_>=0)::close(fd_);}
    DurableProgressLog(const DurableProgressLog&)=delete;
    DurableProgressLog& operator=(const DurableProgressLog&)=delete;
    void append(const std::string& line) {
        std::size_t offset=0;
        while(offset<line.size()) {
            const auto written=::write(fd_,line.data()+offset,line.size()-offset);
            if(written<0&&errno==EINTR)continue;
            if(written<=0)throw std::system_error(errno,std::generic_category(),"write joint progress log");
            offset+=static_cast<std::size_t>(written);
        }
        if(::fsync(fd_)!=0)throw std::system_error(errno,std::generic_category(),"sync joint progress log");
    }
    void iteration(const SlqpIteration& row,std::chrono::steady_clock::time_point start) {
        const auto encoded=encodeSlqpIteration(row);
        const auto recordHash=manufacturing::computeContentSha256(encoded);
        const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now()-start).count();
        std::ostringstream event;event.imbue(std::locale::classic());event<<std::setprecision(17)
            <<"{\"record_type\":\"iteration\",\"iteration\":"<<row.iteration
            <<",\"elapsed_ns\":"<<elapsed<<",\"objective\":"<<row.objective
            <<",\"equality_inf\":"<<row.equalityInfinityNorm
            <<",\"projected_gradient_inf\":"<<row.projectedGradientInfinityNorm
            <<",\"kkt_inf\":"<<row.kktInfinityNorm<<",\"accepted_step\":"<<row.acceptedStep
            <<",\"line_search_evaluations\":"<<row.lineSearchEvaluations
            <<",\"active_variables\":"<<row.activeVariables
            <<",\"membership_sha256\":\""<<row.membershipHash
            <<"\",\"full_iteration_sha256\":\""<<recordHash<<"\"}\n";
        append(event.str());
    }
private:
    int fd_=-1;
};
std::string jsonEscape(const std::string& value) {
    std::ostringstream out;
    for(const unsigned char c:value) {
        switch(c) {
        case '"': out<<"\\\"";break;
        case '\\': out<<"\\\\";break;
        case '\b': out<<"\\b";break;
        case '\f': out<<"\\f";break;
        case '\n': out<<"\\n";break;
        case '\r': out<<"\\r";break;
        case '\t': out<<"\\t";break;
        default:
            if(c<0x20)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec;
            else out<<char(c);
        }
    }
    return out.str();
}
std::string hashState(const Eigen::VectorXd& state) {
    std::ostringstream bytes;bytes.imbue(std::locale::classic());bytes<<std::hexfloat;
    for(Eigen::Index i=0;i<state.size();++i)bytes<<state[i]<<'\n';
    return manufacturing::computeContentSha256(bytes.str());
}
SourceCagePaperFunctional fixture(int level=4) {
    const auto outer=cube(.3125,.6875),design=cube(.425,.575);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,level,4,1000);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;
    pc.inertiaWeight=0.;pc.centerWeight=0.;
    return SourceCagePaperFunctional(SourceCageMomentFunctional(
        SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.17,mc),pc);
}
SourceCagePaperFunctional nontrivialJointFixture() {
    const auto outer=cube(.3,.7),design=cube(.34,.66);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,5,1000);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"frozen 26-handle fixture obtains actual BBW");
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;
    pc.inertiaWeight=.8;pc.centerWeight=.6;pc.fillWeight=.03;pc.shapeWeight=.04;pc.displacementWeight=.05;
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
    frame.rotation=Eigen::AngleAxisd(.17,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    return SourceCagePaperFunctional(SourceCageMomentFunctional(
        SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.3,mc),pc,frame);
}
void candidateChecks(const SourceCageBinaryMaterialAudit&,bool);
int nontrivialJointMaterialTest(const std::string& progressPath={}) {
    const auto f=nontrivialJointFixture();
    require(f.handleCount()==26&&f.betaCount()==8,"fixed joint fixture retains all handles and fill variables");
    Eigen::VectorXd beta=Eigen::VectorXd::LinSpaced(8,.17,.66);
    auto initial=f.initialState(beta,.23);
    for(std::size_t h=0;h<f.handleCount();++h) {
        initial[f.handleOffset()+4*h]=1.04+.0003*h;
        initial[f.handleOffset()+4*h+1]=.007+.0001*h;
        initial[f.handleOffset()+4*h+2]=-.011+.00007*h;
        initial[f.handleOffset()+4*h+3]=.004-.00006*h;
    }
    SourceCageAdaptiveSolverSettings settings;settings.maximumPasses=4;
    settings.solve.source.solver.maximumIterations=120;
    settings.solve.source.solver.maximumLineSearchSteps=24;
    settings.solve.source.solver.lbfgsMemory=8;
    // Geometry-evidence preflight for 496 primitives admits 400 full checker
    // calls within the 24 GiB retained-evidence envelope. This is a storage
    // ceiling, not a runtime guarantee: a measured baseline guard call takes
    // about 71 seconds at this primitive count. Optimizer iterations/line
    // searches and scientific tolerances remain unchanged. Exhausting either
    // independent evidence budget is RESOURCE_LIMIT, never convergence.
    settings.solve.budget.maximumGuardCalls=400;
    settings.solve.budget.maximumRetainedStageRecords=50000000;
    settings.solve.budget.maximumGuardedBytes=std::size_t(24)*1024*1024*1024;
    settings.budget.maximumKnownBufferBytes=std::size_t(24)*1024*1024*1024;
    const auto resources=planSourceCageGuardedSlqpResources(f,settings.solve);
    std::cout<<"joint-material resource-plan primitives="<<f.source().boundary().primitiveBoundaries().size()
        <<" guard-calls="<<settings.solve.budget.maximumGuardCalls
        <<" stages-per-call="<<resources.stagesPerGuard<<" stage-capacity="<<resources.stageCapacity
        <<" required-bytes="<<resources.requiredBytes<<" cap-bytes="<<settings.solve.budget.maximumGuardedBytes<<std::endl;
    std::unique_ptr<DurableProgressLog> progressLog;
    const auto progressStart=std::chrono::steady_clock::now();
    if(!progressPath.empty()) {
        progressLog=std::make_unique<DurableProgressLog>(progressPath);
        std::ostringstream start;start.imbue(std::locale::classic());start
            <<"{\"record_type\":\"joint_run_start\",\"functional_sha256\":\""<<f.configSha256()
            <<"\",\"upstream_sha256\":\""<<f.upstreamSha256()
            <<"\",\"initial_state_sha256\":\""<<hashState(initial)
            <<"\",\"handles\":"<<f.handleCount()<<",\"fill_variables\":"<<f.betaCount()
            <<",\"maximum_passes\":"<<settings.maximumPasses
            <<",\"maximum_iterations_per_pass\":"<<settings.solve.source.solver.maximumIterations
            <<",\"maximum_line_search_steps\":"<<settings.solve.source.solver.maximumLineSearchSteps
            <<",\"maximum_guard_calls\":"<<settings.solve.budget.maximumGuardCalls<<"}\n";
        progressLog->append(start.str());
        settings.solve.source.solver.iterationObserver=[&](const SlqpIteration& row) {
            progressLog->iteration(row,progressStart);
        };
    }
    std::cerr<<"joint-material solve started; durable SLQP iteration progress="
        <<(!progressPath.empty()?progressPath:"disabled")<<std::endl;
    const auto run=solveSourceCageAdaptive(f,initial,settings);
    const auto& facts=run.facts();
    if(progressLog) {
        std::ostringstream terminal;terminal.imbue(std::locale::classic());terminal
            <<"{\"record_type\":\"joint_run_terminal\",\"outcome\":"<<int(facts.outcome)
            <<",\"passes\":"<<facts.passes.size()<<",\"scientific_sha256\":\""
            <<facts.scientificSha256<<"\",\"terminal_reason\":\""
            <<jsonEscape(facts.terminalReason.data())<<"\"}\n";
        progressLog->append(terminal.str());
    }
    std::cout<<"joint-material outcome="<<int(facts.outcome)<<" passes="<<facts.passes.size()
        <<" reason="<<facts.terminalReason.data()<<" scientific="<<facts.scientificSha256<<'\n';
    if(facts.outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit &&
       !run.lastSolved() && facts.passes.size()==1 && !facts.passes.front().solve &&
       std::string(facts.terminalReason.data()).find("aggregate known-buffer reservation cap before work")!=std::string::npos) {
        std::cout<<"JOINT_METHOD_BLOCKED: complete four-pass provenance reservation exceeds 24 GiB before optimization; no scientific solve was attempted\n";
        return 0;
    }
    require(facts.outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete&&run.lastSolved(),
        "frozen nonzero TopMassAndLowering adaptive solve reaches a relaxed endpoint");
    bool acceptedFillUpdate=false,acceptedHandleUpdate=false;
    for(const auto& pass:facts.passes) {
        if(!pass.source||!pass.solve||
           !pass.solve->facts().source.numericalResultAvailable)continue;
        const auto& passFunctional=*pass.source;
        const auto& trace=pass.solve->facts().source.numerical.trace;
        acceptedFillUpdate=acceptedFillUpdate||acceptedTraceBlockChanged(
            pass.initialState,trace,0,passFunctional.betaCount());
        acceptedHandleUpdate=acceptedHandleUpdate||acceptedTraceBlockChanged(
            pass.initialState,trace,passFunctional.handleOffset(),
            4*passFunctional.handleCount());
    }
    require(acceptedFillUpdate,"accepted same-partition SLQP steps change fill variables");
    require(acceptedHandleUpdate,"accepted same-partition SLQP steps change all-handle state");
    const auto binary=validateSourceCageBinaryMaterial(run);
    require(binary.facts().outcome==SourceCageBinaryMaterialOutcome::Validated&&binary.validatedMaterial(),
        "frozen joint endpoint survives original beta rounding and full material checks");
    candidateChecks(binary,true);
    std::cout<<"METHOD_INTEGRATION nonzero TopMassAndLowering -> adaptive -> binary material PASS\n";
    return 0;
}
int jointGeometryProfileTest(const std::string& profilePath) {
    const auto f=nontrivialJointFixture();
    Eigen::VectorXd beta=Eigen::VectorXd::LinSpaced(8,.17,.66);
    auto initial=f.initialState(beta,.23);
    for(std::size_t h=0;h<f.handleCount();++h) {
        initial[f.handleOffset()+4*h]=1.04+.0003*h;
        initial[f.handleOffset()+4*h+1]=.007+.0001*h;
        initial[f.handleOffset()+4*h+2]=-.011+.00007*h;
        initial[f.handleOffset()+4*h+3]=.004-.00006*h;
    }
    DurableProgressLog log(profilePath);
    std::ostringstream start;start.imbue(std::locale::classic());start
        <<"{\"record_type\":\"geometry_profile_start\",\"functional_sha256\":\""
        <<f.configSha256()<<"\",\"upstream_sha256\":\""<<f.upstreamSha256()
        <<"\",\"initial_state_sha256\":\""<<hashState(initial)<<"\",\"handles\":"
        <<f.handleCount()<<",\"fill_variables\":"<<f.betaCount()<<"}\n";
    log.append(start.str());
    const auto begin=std::chrono::steady_clock::now();
    SourceCageGeometryChecker checker(f,SourceCageGeometrySettings{});
    const auto audit=checker.check(initial);
    const auto& facts=audit.facts();
    std::array<std::uint64_t,7> counts{},totalNs{},maxNs{};
    std::array<std::uint64_t,6> pairOutcomeCounts{},pairOutcomeNs{};
    for(const auto& row:facts.records) {
        const auto stage=static_cast<std::size_t>(row.stage);
        if(stage>=counts.size())throw std::runtime_error("unknown source geometry stage in profiler");
        ++counts[stage];totalNs[stage]+=row.elapsedNanoseconds;
        maxNs[stage]=std::max(maxNs[stage],row.elapsedNanoseconds);
        if(row.stage==SourceCageGeometryStage::PrimitivePair) {
            const auto outcome=static_cast<std::size_t>(row.outcome);
            if(outcome>=pairOutcomeCounts.size())throw std::runtime_error("unknown pair-stage outcome in profiler");
            ++pairOutcomeCounts[outcome];pairOutcomeNs[outcome]+=row.elapsedNanoseconds;
        }
    }
    for(std::size_t stage=0;stage<counts.size();++stage) {
        std::ostringstream row;row.imbue(std::locale::classic());row
            <<"{\"record_type\":\"geometry_stage_summary\",\"stage\":"<<stage
            <<",\"count\":"<<counts[stage]<<",\"total_ns\":"<<totalNs[stage]
            <<",\"max_ns\":"<<maxNs[stage]<<"}\n";
        log.append(row.str());
    }
    for(std::size_t outcome=0;outcome<pairOutcomeCounts.size();++outcome) {
        std::ostringstream row;row.imbue(std::locale::classic());row
            <<"{\"record_type\":\"pair_outcome_summary\",\"outcome\":"<<outcome
            <<",\"count\":"<<pairOutcomeCounts[outcome]<<",\"total_ns\":"<<pairOutcomeNs[outcome]<<"}\n";
        log.append(row.str());
    }
    const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-begin).count();
    std::ostringstream terminal;terminal.imbue(std::locale::classic());terminal
        <<"{\"record_type\":\"geometry_profile_terminal\",\"outcome\":\""
        <<sourceCageGeometryOutcomeToken(facts.outcome)<<"\",\"input_sha256\":\""
        <<facts.inputSha256<<"\",\"source_output_sha256\":\""<<facts.sourceOutputSha256
        <<"\",\"scientific_sha256\":\""<<facts.scientificSha256
        <<"\",\"records\":"<<facts.records.size()<<",\"pairs_visited\":"<<facts.pairsVisited
        <<",\"exact_pair_relations\":"<<facts.exactPairRelations
        <<",\"elapsed_ns\":"<<elapsed<<"}\n";
    log.append(terminal.str());
    std::cout<<"joint geometry profile outcome="<<sourceCageGeometryOutcomeToken(facts.outcome)
        <<" records="<<facts.records.size()<<" pairs="<<facts.pairsVisited
        <<" exact_pairs="<<facts.exactPairRelations<<" elapsed_ns="<<elapsed
        <<" sha256="<<facts.scientificSha256<<" log="<<profilePath<<'\n';
    return facts.outcome==SourceCageGeometryOutcome::ChecksPassed?0:1;
}
Eigen::VectorXd state(const SourceCagePaperFunctional& f,double beta) {
    auto x=f.initialState(Eigen::VectorXd::Constant(f.betaCount(),beta),.23);
    for(std::size_t h=0;h<f.handleCount();++h) {
        x[f.handleOffset()+4*h]=1.02;
        x.segment<3>(f.handleOffset()+4*h+1)<<.011,-.007,.019;
    }
    return x;
}
// KNOWN_SOLUTION_VALIDATION: fixed principal frame and zero objective weights.
// A nonzero Givens angle is retained, with the transverse basis counter-rotated.
SourceCagePaperFunctional stationary(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto m=f.evaluateValues(x).world.finalMoments.values();
    SourceCagePaperFrame frame;frame.origin=m.segment<3>(1)/m[0];
    Eigen::Matrix3d q;q<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(q-m[0]*frame.origin*frame.origin.transpose());
    Eigen::Matrix3d base=eig.eigenvectors().transpose();if(base.determinant()<0)base.row(0)*=-1;
    const double c=std::cos(x[f.angleOffset()]),s=std::sin(x[f.angleOffset()]);
    Eigen::Matrix3d g;g<<c,s,0,-s,c,0,0,0,1;frame.rotation=g*base;
    return SourceCagePaperFunctional(f.source(),f.config(),frame);
}
SourceCageAdaptiveSolverSettings solverSettings() {
    SourceCageAdaptiveSolverSettings s;s.maximumPasses=3;
    s.solve.source.solver.maximumIterations=1;s.solve.source.solver.maximumLineSearchSteps=2;
    s.solve.budget.maximumGuardCalls=8;return s;
}
SourceCageAdaptiveSolverAudit endpoint(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    auto audit=solveSourceCageAdaptive(stationary(f,x),x,solverSettings());
    require(audit.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete,"fixture must actually converge and complete");
    return audit;
}
void candidateChecks(const SourceCageBinaryMaterialAudit& audit,bool requireNonidentity=true) {
    const auto m=audit.candidate();require(bool(m),"candidate geometry retained");
    const auto& f=m->functional();const auto values=f.evaluateValues(m->binaryState());
    require(identical(m->relaxedState().tail(m->relaxedState().size()-f.betaCount()),
        m->binaryState().tail(m->binaryState().size()-f.betaCount())),"rounding preserves phi and all handles bit-for-bit");
    require((m->outerBoundary().triangles.array()==f.source().boundary().binding().source().outerFaces().array()).all(),"actual original outer connectivity retained");
    require((m->outerBoundary().vertices-values.world.deformedPoints.topRows(m->outerBoundary().vertices.rows())).norm()==0.,"emitted outer coordinates are actual binary deformation");
    if(requireNonidentity)
        require((m->outerBoundary().vertices-f.source().boundary().binding().source().outerVertices()).norm()>1e-3,"nonidentity deformation exercised");
    require((audit.facts().worldIntegratedMoments-values.world.finalMoments.values()).cwiseAbs().maxCoeff()<1e-12,"all ten independently emitted moments agree with actual functional");
    require((audit.facts().spinIntegratedMoments-values.spinMoments.values()).cwiseAbs().maxCoeff()<1e-12,"fixed frame independent transform agrees");
    require((audit.facts().rawConstraints-values.constraints).norm()<1e-12,"Givens constraints recomputed on actual binary output");
    require((audit.facts().centralInertia-values.spinMoments.inertiaAboutCenter()).norm()<1e-12,"central inertia uses independent binary moments in fixed nonidentity paper frame");
    require(audit.facts().centralInertiaPsd&&audit.facts().centralInertiaEigenvalues.minCoeff()>0.,"positive material has checked central inertia spectrum");
    require(identical(audit.facts().stateDelta,m->binaryState()-m->relaxedState()),"full raw rounding delta retained");
    Eigen::VectorXd normalized=values.constraints;
    for(Eigen::Index i=0;i<normalized.size();++i)normalized[i]/=f.source().density()*std::pow(f.source().boundary().binding().source().config().rootSide,i<2?4:5);
    require((audit.facts().normalizedConstraints-normalized).norm()<1e-12,"binary Givens residuals use fixed INPUT scales");
    require(audit.facts().geometry&&audit.facts().exactBoundary,"both actual checker and final boundary audit run");
}
}
int main(int argc,char** argv){try {
    if(argc==2&&std::string(argv[1])=="--joint-index-guard")return jointIndexGuardTest();
    if(argc==2&&std::string(argv[1])=="--joint-guard-once") {
        const auto f=nontrivialJointFixture();
        Eigen::VectorXd beta=Eigen::VectorXd::LinSpaced(8,.17,.66);
        auto initial=f.initialState(beta,.23);
        for(std::size_t h=0;h<f.handleCount();++h) {
            initial[f.handleOffset()+4*h]=1.04+.0003*h;
            initial[f.handleOffset()+4*h+1]=.007+.0001*h;
            initial[f.handleOffset()+4*h+2]=-.011+.00007*h;
            initial[f.handleOffset()+4*h+3]=.004-.00006*h;
        }
        std::cerr<<"joint-guard-once start primitives="<<f.source().boundary().primitiveBoundaries().size()<<std::endl;
        const auto start=std::chrono::steady_clock::now();
        const auto audit=SourceCageGeometryChecker(f,SourceCageGeometrySettings{}).check(initial);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto& facts=audit.facts();
        std::cout<<"joint-guard-once outcome="<<sourceCageGeometryOutcomeToken(facts.outcome)
            <<" seconds="<<seconds<<" primitives="<<facts.primitiveCount<<" pairs="<<facts.pairCount
            <<" visited="<<facts.pairsVisited<<" exact="<<facts.exactPairRelations
            <<" records="<<facts.records.size()<<" planned_bytes="<<facts.plannedCheckerBytes
            <<" reason="<<(facts.records.empty()?"no-record":facts.records.back().reason.data())
            <<" digest="<<facts.scientificSha256<<std::endl;
        require(facts.outcome==SourceCageGeometryOutcome::ChecksPassed,
            "frozen 26-handle geometry guard remains admitted under its registered checker-byte cap");
        std::array<std::uint64_t,7> stageNanoseconds{};
        std::array<std::size_t,7> stageRecords{};
        std::array<std::uint64_t,6> pairOutcomeNanoseconds{};
        std::array<std::size_t,6> pairOutcomeRecords{};
        for(const auto& record:facts.records) {
            const auto stage=static_cast<std::size_t>(record.stage);
            if(stage<stageNanoseconds.size()) {
                stageNanoseconds[stage]+=record.elapsedNanoseconds;
                ++stageRecords[stage];
            }
            if(record.stage==SourceCageGeometryStage::PrimitivePair) {
                const auto outcome=static_cast<std::size_t>(record.outcome);
                require(outcome<pairOutcomeNanoseconds.size(),"primitive-pair outcome is in range");
                pairOutcomeNanoseconds[outcome]+=record.elapsedNanoseconds;
                ++pairOutcomeRecords[outcome];
            }
        }
        std::uint64_t pairOutcomeTime=0;
        std::size_t pairOutcomeCount=0;
        for(std::size_t outcome=0;outcome<pairOutcomeNanoseconds.size();++outcome) {
            pairOutcomeTime+=pairOutcomeNanoseconds[outcome];
            pairOutcomeCount+=pairOutcomeRecords[outcome];
            if(pairOutcomeRecords[outcome])std::cout<<"joint-guard-pair-outcome outcome="<<outcome
                <<" records="<<pairOutcomeRecords[outcome]<<" elapsed_ns="<<pairOutcomeNanoseconds[outcome]<<std::endl;
        }
        require(pairOutcomeCount==facts.pairsVisited,"all visited pairs retain one diagnostic outcome");
        require(pairOutcomeTime==stageNanoseconds[static_cast<std::size_t>(SourceCageGeometryStage::PrimitivePair)],
                "pair outcome timing partitions the complete pair stage");
        for(std::size_t stage=0;stage<stageNanoseconds.size();++stage)
            if(stageRecords[stage])std::cout<<"joint-guard-stage stage="<<stage
                <<" records="<<stageRecords[stage]<<" elapsed_ns="<<stageNanoseconds[stage]<<std::endl;
        return facts.outcome==SourceCageGeometryOutcome::ChecksPassed?0:1;
    }
    if(argc==2&&std::string(argv[1])=="--joint-material")return nontrivialJointMaterialTest();
    if(argc==3&&std::string(argv[1])=="--joint-material-progress")return nontrivialJointMaterialTest(argv[2]);
    if(argc==3&&std::string(argv[1])=="--joint-material-geometry-profile")return jointGeometryProfileTest(argv[2]);
    static_assert(!std::is_default_constructible_v<SourceCageBinaryMaterialAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageBinaryMaterial>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<SourceCageBinaryMaterial>().binaryState())>>);
    const auto f=fixture();require(f.betaCount()==8,"eight actual source owners");
    for(double beta:{0.,.25,.5,.75,1.}) {
        const auto x=state(f,beta);const auto ep=endpoint(f,x);
        const auto audit=validateSourceCageBinaryMaterial(ep);candidateChecks(audit);
        require(audit.facts().outcome==SourceCageBinaryMaterialOutcome::Validated&&audit.validatedMaterial(),"stationary symmetric binary validates");
        require((audit.candidate()->binaryState().head(audit.candidate()->functional().betaCount()).array()==(beta>=.5?1.:0.)).all(),"zero/one and exact .5 tie round only beta");
        require(audit.facts().emptyCavity==(beta<.5),"empty cavity explicit");
        require(audit.candidate()->carvedBoundary().triangles.rows()==(beta<.5?0:96),"internal faces cancelled in carved union");
        require(audit.facts().exactBoundary->facts().boundaryValid,"actual emitted material exact audit valid");
        require(audit.facts().normalizedMomentResiduals.cwiseAbs().maxCoeff()<1e-10,"fixed input normalization gate");
        require(audit.facts().inputSha256.size()==64&&audit.facts().functionalSha256.size()==64&&audit.facts().partitionSha256.size()==64&&
            audit.facts().configSha256.size()==64&&audit.facts().outputSha256.size()==64,"complete immutable identity digests");
        if(beta==.5) {
            const auto repeated=validateSourceCageBinaryMaterial(ep);
            require(repeated.facts().scientificSha256==audit.facts().scientificSha256,"deterministic complete evidence digest");
        }
        std::cout<<"KNOWN_SOLUTION_VALIDATION beta="<<beta<<" PASS\n";
    }
    // METHOD_INTEGRATION: a real shape-regularized cage update must survive
    // adaptive completion, paper beta rounding, material export, and the
    // independent final-mesh moment/constraint checks.
    auto deformState=state(f,.25);
    for(std::size_t h=0;h<f.handleCount();++h) {
        deformState[f.handleOffset()+4*h]=1.02;
        deformState.segment<3>(f.handleOffset()+4*h+1).setZero();
    }
    const auto centered=stationary(f,deformState);
    auto shapeConfig=centered.config();shapeConfig.shapeWeight=.25;shapeConfig.displacementWeight=.25;
    const SourceCagePaperFunctional shaped(centered.source(),shapeConfig,centered.frame());
    auto shapeSettings=solverSettings();shapeSettings.solve.source.solver.maximumIterations=12;
    shapeSettings.solve.source.solver.maximumLineSearchSteps=12;
    shapeSettings.solve.budget.maximumGuardCalls=256;
    const auto shapeRun=solveSourceCageAdaptive(shaped,deformState,shapeSettings);
    require(shapeRun.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete&&shapeRun.lastSolved(),
        "nontrivial shape-regularized cage solve reaches a fresh adaptive endpoint");
    const auto& shapeTrace=shapeRun.lastSolved()->solve().facts().source.numerical.trace;
    bool acceptedShapeStep=false;
    for(const auto& row:shapeTrace)for(const auto& trial:row.objectiveTrials)acceptedShapeStep|=trial.accepted;
    require(acceptedShapeStep,"joint cage shape update has an actually accepted objective step");
    require((shapeRun.lastSolved()->state().tail(4*f.handleCount())-
        deformState.tail(4*f.handleCount())).norm()>1e-8,"optimized cage handle coordinates change");
    const auto shapeMaterial=validateSourceCageBinaryMaterial(shapeRun);
    require(shapeMaterial.facts().outcome==SourceCageBinaryMaterialOutcome::Validated&&
        shapeMaterial.validatedMaterial(),"accepted adaptive cage result reaches audited binary material output");
    candidateChecks(shapeMaterial,false);
    require(shapeMaterial.facts().exactBoundary->facts().boundaryValid&&
        shapeMaterial.facts().normalizedConstraints.cwiseAbs().maxCoeff()<=1e-6,
        "independent final material boundary and binary constraints pass");
    std::cout<<"METHOD_INTEGRATION accepted cage update -> adaptive -> binary material PASS\n";
    auto mixed=state(f,.4);
    const auto& leaves=f.source().boundary().binding().domain().leaves();
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==8)mixed[i]=.6;
    mixed[0]=.5;
    const auto mixedEndpoint=endpoint(f,mixed);auto mixedAudit=validateSourceCageBinaryMaterial(mixedEndpoint);candidateChecks(mixedAudit);
    require(mixedAudit.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!mixedAudit.validatedMaterial(),"genuine rounding-induced fixed-frame failure stays BINARY_VALIDATION_FAILED");
    require(mixedAudit.facts().rawConstraints.cwiseAbs().maxCoeff()>1e-6,"raw rounding constraint failure visible");
    require(mixedAudit.candidate()->binaryState()[0]==1.,"mixed exact .5 tie hollowed");
    // The same geometry at low density has small RAW first moments, while
    // dimensionless fixed-input residuals retain the genuine axis violation.
    const SourceCagePaperFunctional dilute(SourceCageMomentFunctional(f.source().boundary(),.001,f.source().config()),f.config());
    const auto diluteEndpoint=endpoint(dilute,mixed);
    const auto diluteAudit=validateSourceCageBinaryMaterial(diluteEndpoint);candidateChecks(diluteAudit);
    require(diluteAudit.facts().rawConstraints.cwiseAbs().maxCoeff()<1e-6,"low-density raw residual alone would pass");
    require(diluteAudit.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!diluteAudit.validatedMaterial(),
        "fixed input density normalization must reject despite small raw residual");
    require(diluteAudit.facts().normalizedConstraints.cwiseAbs().maxCoeff()>1e-6,"dimensionless rounding violation retained");
    SourceCageBinaryMaterialSettings smallGeometry;smallGeometry.geometry.maximumPoints=1;
    const auto bothFail=validateSourceCageBinaryMaterial(diluteEndpoint,smallGeometry);
    require(bothFail.facts().geometry&&bothFail.facts().geometry->facts().outcome==SourceCageGeometryOutcome::ResourceLimit&&
        bothFail.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!bothFail.validatedMaterial(),
        "established binary constraint failure stays failed when checker cannot complete");
    auto relaxedCopy=mixedAudit.facts();relaxedCopy.outcome=SourceCageBinaryMaterialOutcome::Validated;
    require(!mixedAudit.validatedMaterial(),"copied facts cannot authorize material");
    const auto& retained=mixedAudit.facts();const auto candidate=mixedAudit.candidate();const auto moved=std::move(mixedAudit);
    require(retained.scientificSha256==moved.facts().scientificSha256,"immutable facts survive move");
    rejects([&]{(void)mixedAudit.facts();},"moved-from access rejected");
    require(candidate->functional().evaluateValues(candidate->binaryState()).outputSha256==moved.facts().binaryFunctionalOutputSha256,"candidate retains owned source lifetime");
    auto diagonal=state(f,.2);
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==leaves[i].coordinates[1])diagonal[i]=.8;
    const auto nonmanifold=validateSourceCageBinaryMaterial(endpoint(f,diagonal));candidateChecks(nonmanifold);
    require(nonmanifold.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!nonmanifold.validatedMaterial()&&
        nonmanifold.facts().geometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed&&
        !nonmanifold.facts().exactBoundary->facts().boundaryValid,"edge contact cavity union honestly rejected by final boundary audit without repair");
    const auto complete=endpoint(f,state(f,1.));
    SourceCageBinaryMaterialSettings badGeometry;badGeometry.geometry.minimumOuterDomainClearance=.2;
    const auto geometryRejected=validateSourceCageBinaryMaterial(complete,badGeometry);
    require(geometryRejected.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!geometryRejected.validatedMaterial()&&
        geometryRejected.facts().geometry->facts().outcome==SourceCageGeometryOutcome::RejectedGeometry,"fresh explicit binary geometry settings can reject");
    for(double t:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        SourceCageBinaryMaterialSettings bad;bad.constraintTolerance=t;
        rejects([&]{(void)validateSourceCageBinaryMaterial(complete,bad);},"finite positive tolerance required");
    }
    for(std::size_t cap:{std::size_t(0),std::size_t(1),std::numeric_limits<std::size_t>::max()}) {
        SourceCageBinaryMaterialSettings bad;bad.maximumTriangles=cap;
        const auto limited=validateSourceCageBinaryMaterial(complete,bad);
        require(limited.facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit&&!limited.candidate()&&!limited.validatedMaterial(),"invalid or overflowing mesh budgets rejected before emission");
    }
    SourceCageBinaryMaterialSettings tiny;tiny.maximumOwnedBytes=1;
    require(validateSourceCageBinaryMaterial(complete,tiny).facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit,"owned buffer cap enforced");
    tiny=SourceCageBinaryMaterialSettings{};tiny.maximumExactDiagnosticBytes=1;
    require(validateSourceCageBinaryMaterial(complete,tiny).facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit,"final exact diagnostic cap enforced before allocation");
    const auto coarse=fixture(3);const auto cx=state(coarse,.25);const auto cf=stationary(coarse,cx);
    auto limit=solverSettings();limit.maximumPasses=1;
    const auto pending=solveSourceCageAdaptive(cf,cx,limit);
    require(bool(pending.pendingCheckpoint()),"actual changed checkpoint exists");
    const auto notAdmitted=validateSourceCageBinaryMaterial(pending);
    require(notAdmitted.facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted&&!notAdmitted.candidate(),"pending checkpoint never admitted as terminal");
    auto limitedSettings=solverSettings();limitedSettings.budget.maximumKnownBufferBytes=1;
    const auto resource=solveSourceCageAdaptive(cf,cx,limitedSettings);
    require(validateSourceCageBinaryMaterial(resource).facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted,"resource endpoint never admitted");
    SourceCagePaperFunctional nonstationary(f.source(),SourceCagePaperFunctionalConfig{},{});
    auto nx=nonstationary.initialState(Eigen::VectorXd::Constant(8,.25),0.);
    const auto stalled=solveSourceCageAdaptive(nonstationary,nx,solverSettings());
    require(stalled.facts().outcome==SourceCageAdaptiveSolverOutcome::SolveNotConverged,"actual nonconverged fixture");
    require(validateSourceCageBinaryMaterial(stalled).facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted,"nonconverged endpoint never admitted");
    std::cout<<"binary material regression PASS; bounded pre-seal object only; no sealing, feasibility, physics or paper grant\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
