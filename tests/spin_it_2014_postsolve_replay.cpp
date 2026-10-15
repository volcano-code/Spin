#if __has_include("spin2014/SourceCheckpointCodec.h")
#include "spin2014/SourceCheckpointCodec.h"
#include "spin2014/SourceMaterialBoundary.h"
#include "spin2014/SlqpIterationStore.h"
#include "manufacturing/SlicerValidation.h"
#include "source_surface_test_utils.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

using namespace spin2014;
using namespace surface_test;

namespace {
SourceAdaptiveSpinRequest fixture(const std::filesystem::path& store) {
    const auto outer=box({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    const auto inner=box({.125,-.125,.625},{.375,.125,.875});
    SourceSurfaceOctreeConfig c;c.maximumDepth=5;c.rootOrigin=Eigen::Vector3d(-1,-1,0);c.rootSide=2;
    SourceAdaptiveSpinRequest r;
    r.source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,inner.v,inner.f,c);
    r.initialLevel=2;r.maximumOptimizationLevel=5;r.maximumAdaptivePasses=8;r.laplacianWeight=1;
    r.termination=SourceAdaptiveTermination::PaperPostTransitionV1;
    r.slqp.maximumIterations=300;r.slqp.equalityTolerance=1e-10;r.slqp.projectedGradientTolerance=1e-8;
    r.slqp.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;
    r.iterationStoreDirectory=store.string();
    return r;
}

std::string materialCheckpoint(const SourceAdaptiveCheckpoint& checkpoint,
    const SourceMaterialBoundarySettings& boundary) {
    auto bytes=std::string("SCMAT001")+hashSourceMaterialBoundarySettings(boundary)+
        encodeSourceAdaptiveCheckpoint(checkpoint);
    return bytes+manufacturing::computeContentSha256(bytes);
}

void expectReject(const SourceAdaptiveSpinRequest& request,const std::string& checkpoint,
    const std::string& trace,const std::string& boundaryHash,const char* message) {
    bool rejected=false;
    try {(void)replaySourceAdaptivePostSolve(request,checkpoint,boundaryHash,trace,1);}
    catch(const std::exception&) { rejected=true; }
    require(rejected,message);
}

void writeTrace(const std::filesystem::path& path,const std::vector<SlqpIteration>& records,
    std::uint64_t cap) {
    SlqpIterationStore store(path.string(),cap);
    for(const auto& record:records)store.append(record);
    store.seal();
}
}

int main(int argc,char** argv) { try {
    require(spin2014G7PostSolvePrePassCount==20,
        "G7 stage-19 mapping must require twenty completed passes");
    require(argc==2,"output directory required");
    const auto root=std::filesystem::absolute(argv[1]) / ("run-"+std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count()));std::filesystem::create_directories(root);
    auto request=fixture(root/"stage-19");
    SourceAdaptiveStageControl first;first.maximumNewPasses=1;
    auto stage19=solveSourceAdaptiveSpinStages(request,first);
    require(stage19.checkpoint && stage19.result.passes.size()==1,
        "fixture did not create a committed stage boundary");
    request.iterationStoreDirectory=(root/"stage-20").string();
    SourceAdaptiveStageControl second;second.maximumNewPasses=1;second.resume=stage19.checkpoint;
    const auto normal=solveSourceAdaptiveSpinStages(request,second);
    require(normal.result.passes.size()==2,"fixture did not execute one post-boundary pass");
    const auto& expected=normal.result.passes.back();
    require(expected.solve.iterationStore && expected.solve.iterationStore->sealed(),
        "fixture did not retain sealed pass-20 trace");
    SourceMaterialBoundarySettings boundary;
    boundary.unionSettings.maximumInputCuboids=1000000;
    boundary.unionSettings.maximumOutputVertices=1000000;
    boundary.unionSettings.maximumOutputTriangles=1000000;
    boundary.unionSettings.maximumIntermediateComplexity=1000000;
    const auto boundaryHash=hashSourceMaterialBoundarySettings(boundary);
    const auto checkpoint=materialCheckpoint(*stage19.checkpoint,boundary);
    const auto tracePath=expected.solve.iterationStore->path();

    // Break caught: a diagnostic accidentally calls the optimizer or rebuilds
    // the transition from anything other than the sealed trace.
    const auto replay=replaySourceAdaptivePostSolve(request,checkpoint,boundaryHash,tracePath,1);
    const auto& actual=replay.result.passes.back();
    require(!replay.optimizerInvoked && replay.binaryValidation=="NOT_RUN" &&
        replay.scientificAcceptance=="NOT_GRANTED", "diagnostic-only scope was weakened");
    require(actual.beforeHash==expected.beforeHash && actual.afterHash==expected.afterHash &&
        actual.levelCounts==expected.levelCounts && actual.outputLevelCounts==expected.outputLevelCounts &&
        actual.split.beforeHash==expected.split.beforeHash && actual.split.afterHash==expected.split.afterHash &&
        actual.merge.beforeHash==expected.merge.beforeHash && actual.merge.afterHash==expected.merge.afterHash,
        "post-solve transition no longer matches the normal pass");
    require(actual.split.changes.size()==expected.split.changes.size() &&
        actual.merge.changes.size()==expected.merge.changes.size(),
        "post-solve split/merge change counts no longer match the normal pass");
    require(actual.split.removedBefore.values()==expected.split.removedBefore.values() &&
        actual.split.removedAfter.values()==expected.split.removedAfter.values() &&
        actual.merge.removedBefore.values()==expected.merge.removedBefore.values() &&
        actual.merge.removedAfter.values()==expected.merge.removedAfter.values(),
        "post-solve raw moments no longer match the normal pass");

    std::vector<SlqpIteration> records;forEachSlqpIteration(expected.solve,[&](const auto& t){records.push_back(t);});
    require(!records.empty(),"fixture trace is empty");
    auto firstState=records;firstState.front().state[0]=1.-firstState.front().state[0];
    writeTrace(root/"wrong-first.slqptrace",firstState,request.iterationStoreMaximumBytes);
    expectReject(request,checkpoint,(root/"wrong-first.slqptrace").string(),boundaryHash,
        "trace with a changed inherited state was accepted");
    auto noFinal=records;noFinal.back().finalState=false;
    writeTrace(root/"no-final.slqptrace",noFinal,request.iterationStoreMaximumBytes);
    expectReject(request,checkpoint,(root/"no-final.slqptrace").string(),boundaryHash,
        "trace without final-state marker was accepted");
    auto nonfinite=records;nonfinite.back().state[0]=std::numeric_limits<double>::quiet_NaN();
    writeTrace(root/"nonfinite.slqptrace",nonfinite,request.iterationStoreMaximumBytes);
    expectReject(request,checkpoint,(root/"nonfinite.slqptrace").string(),boundaryHash,
        "trace with nonfinite terminal state was accepted");
    std::ifstream input(tracePath,std::ios::binary);std::string sealed((std::istreambuf_iterator<char>(input)),{});
    sealed.resize(sealed.size()-1);std::ofstream(root/"truncated.slqptrace",std::ios::binary).write(sealed.data(),sealed.size());
    expectReject(request,checkpoint,(root/"truncated.slqptrace").string(),boundaryHash,
        "truncated or unsealed trace was accepted");
    auto changedRequest=request;changedRequest.initialFill=.25;
    expectReject(changedRequest,checkpoint,tracePath,boundaryHash,
        "checkpoint for changed request controls was accepted");
    require(bool(normal.checkpoint),"fixture did not retain a stage-20 boundary token for rejection coverage");
    const auto committedStage20=materialCheckpoint(*normal.checkpoint,boundary);
    expectReject(request,committedStage20,tracePath,boundaryHash,
        "already committed pass trace was accepted as a post-solve replay");
    std::cout<<"post-solve replay diagnostic contract passed\n";
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }}
#else
#include <iostream>
int main(){std::cerr<<"RED: post-solve diagnostic unavailable\n";return 1;}
#endif
