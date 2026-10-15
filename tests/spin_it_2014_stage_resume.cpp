#if __has_include("spin2014/SourceAdaptiveCheckpoint.h")
#include "spin2014/SourceAdaptiveCheckpoint.h"
#include "spin2014/SlqpIterationStore.h"
#include "source_surface_test_utils.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;

SourceAdaptiveSpinRequest request(bool angle=false) {
    auto outer=angle?cube(-.5,.5):box({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    auto inner=angle?cube(-.25,.25):box({.125,-.125,.625},{.375,.125,.875});
    SourceSurfaceOctreeConfig c;c.maximumDepth=5;c.rootOrigin=angle?Eigen::Vector3d(-1,-1,-1):Eigen::Vector3d(-1,-1,0);c.rootSide=2;
    SourceAdaptiveSpinRequest r;r.source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,inner.v,inner.f,c);
    r.initialLevel=2;r.maximumOptimizationLevel=angle?3:5;r.maximumAdaptivePasses=8;
    r.laplacianWeight=1;r.slqp.maximumIterations=300;r.slqp.equalityTolerance=1e-10;r.slqp.projectedGradientTolerance=1e-8;
    r.slqp.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;
    if(angle){r.objective=Spin2014Objective::YoYo;r.initialPhi=.37;}
    return r;
}
std::string traceBytes(const SlqpLbfgsResult& s) {
    std::string result;forEachSlqpIteration(s,[&](const auto& t){result+=encodeSlqpIteration(t);});return result;
}
void same(const SourceAdaptiveSpinResult& a,const SourceAdaptiveSpinResult& b) {
    require(a.status==b.status && a.traceHash==b.traceHash && a.outputHash==b.outputHash,"resumed scientific hashes/status differ");
    require(a.passes.size()==b.passes.size() && a.relaxedFill==b.relaxedFill && a.binaryFill==b.binaryFill && a.finalPhi==b.finalPhi,"resumed history or final state differs");
    for(std::size_t i=0;i<a.passes.size();++i){
        require(a.passes[i].index==int(i) && b.passes[i].index==int(i),"stage counter reset");
        require(traceBytes(a.passes[i].solve)==traceBytes(b.passes[i].solve),"full iteration bytes differ");
    }
}
int main(int argc,char** argv) {try {
    static_assert(!std::is_default_constructible_v<SourceAdaptiveCheckpoint>,"caller cannot forge checkpoint");
    const auto temp=std::filesystem::path(argc>1?argv[1]:"stage-resume-test-output");
    require(std::filesystem::create_directory(temp),"exclusive test directory required");
    std::size_t boundaries=0;
    for(bool angle:{false,true})for(auto policy:{SourceAdaptiveTermination::UnchangedPartitionKktV1,SourceAdaptiveTermination::PaperPostTransitionV1}) {
        auto r=request(angle);r.termination=policy;auto baseline=solveSourceAdaptiveSpin(r);
        const bool immediate=angle && policy==SourceAdaptiveTermination::PaperPostTransitionV1;
        if(immediate)require(baseline.passes.size()==1,"paper resolution should finish after first angle split");
        else require(baseline.passes.size()>=2,"fixture must exercise multiple actual stages");
        SourceAdaptiveStageControl control;control.maximumNewPasses=1;
        auto chunk=solveSourceAdaptiveSpinStages(r,control);
        if(immediate)require(!chunk.checkpoint,"complete first pass incorrectly paused");
        while(chunk.checkpoint) {
            ++boundaries;
            require(chunk.result.status==SourceAdaptiveSpinStatus::StageBoundaryPaused && !chunk.result.binaryCheck && !chunk.result.binaryCuboids,"paused state cannot grant binary acceptance");
            require(chunk.checkpoint->completedPasses()==chunk.result.passes.size(),"checkpoint history omitted");
            auto prior=chunk.checkpoint;
            // A returned diagnostic result is mutable, but is not the resume token.
            chunk.result.relaxedFill.setZero();chunk.result.passes.clear();
            auto reconstructed=request(angle);reconstructed.termination=policy;
            control.resume=prior;chunk=solveSourceAdaptiveSpinStages(reconstructed,control);
            require(chunk.result.passes.size()>prior->completedPasses(),"resume replayed zero stages");
        }
        same(baseline,chunk.result);
        std::cout<<"angle="<<angle<<" policy="<<sourceAdaptiveTerminationToken(policy)<<" passes="<<baseline.passes.size()<<" output="<<baseline.outputHash<<'\n';
    }
    auto r=request();SourceAdaptiveStageControl one;one.maximumNewPasses=1;
    auto paused=solveSourceAdaptiveSpinStages(r,one);require(bool(paused.checkpoint),"pause required");
    one.resume=paused.checkpoint;
    auto wrong=r;wrong.density=2;rejects([&]{solveSourceAdaptiveSpinStages(wrong,one);},"changed request accepted");
    wrong=r;wrong.maximumAdaptivePasses+=1;rejects([&]{solveSourceAdaptiveSpinStages(wrong,one);},"global pass budget reset");
    wrong=r;wrong.iterationStoreMaximumBytes-=1;rejects([&]{solveSourceAdaptiveSpinStages(wrong,one);},"trace budget changed");
    wrong=r;wrong.iterationStoreDirectory=(temp/"wrong-disk-mode").string();
    rejects([&]{solveSourceAdaptiveSpinStages(wrong,one);},"memory checkpoint switched to disk mode");
    require(!std::filesystem::exists(wrong.iterationStoreDirectory),"wrong mode created output");
    auto zero=one;zero.maximumNewPasses=0;rejects([&]{solveSourceAdaptiveSpinStages(r,zero);},"zero work budget accepted");
    wrong=r;wrong.slqp.maximumIterations=1;auto failed=solveSourceAdaptiveSpinStages(wrong,{});
    require(!failed.checkpoint && failed.result.status==SourceAdaptiveSpinStatus::RelaxationUnknown,"failed solve checkpointed");
    wrong=r;wrong.maximumAdaptivePasses=1;auto capped=solveSourceAdaptiveSpinStages(wrong,{});
    require(!capped.checkpoint && capped.result.status==SourceAdaptiveSpinStatus::AdaptivePassLimit,"global cap made resumable");
    wrong=r;wrong.maximumLeafCount=2;auto resource=solveSourceAdaptiveSpinStages(wrong,{});
    require(!resource.checkpoint && resource.result.status==SourceAdaptiveSpinStatus::ResourceLimitReached,"failed adaptation checkpointed");
    // Disk history must stay streamed, preserve cumulative bytes, and validate
    // the old immutable files before creating a new attempt directory.
    auto disk=request();disk.iterationStoreDirectory=(temp/"first").string();
    auto first=solveSourceAdaptiveSpinStages(disk,SourceAdaptiveStageControl{});
    require(bool(first.checkpoint) && first.checkpoint->storedTraceBytes()>0,"streamed pause missing");
    auto bytes=first.checkpoint->storedTraceBytes();SourceAdaptiveStageControl next;next.resume=first.checkpoint;
    auto memory=request();rejects([&]{solveSourceAdaptiveSpinStages(memory,next);},"disk checkpoint switched to memory mode");
    auto budget=disk;budget.iterationStoreMaximumBytes=bytes+48;
    budget.iterationStoreDirectory=(temp/"budget-first").string();
    auto budgetFirst=solveSourceAdaptiveSpinStages(budget,{});
    require(bool(budgetFirst.checkpoint),"budget fixture must save first complete stage");
    SourceAdaptiveStageControl budgetNext;budgetNext.resume=budgetFirst.checkpoint;
    budget.iterationStoreDirectory=(temp/"budget-next").string();
    rejects([&]{solveSourceAdaptiveSpinStages(budget,budgetNext);},"resume reset remaining trace-byte budget");
    require(std::filesystem::file_size(temp/"budget-next/pass-000001.slqptrace")<=48,"resumed trace exceeded remaining declared bytes");
    disk.iterationStoreDirectory=(temp/"next").string();auto second=solveSourceAdaptiveSpinStages(disk,next);
    require(bool(second.checkpoint) && second.checkpoint->storedTraceBytes()>bytes,"trace usage reset on resume");
    for(const auto& p:second.result.passes)require(p.solve.iterationStore && p.solve.trace.empty(),"history expanded to resident trace");
    auto baseline=solveSourceAdaptiveSpin(request());
    next.resume=second.checkpoint;next.maximumNewPasses=100;disk.iterationStoreDirectory=(temp/"last").string();
    auto done=solveSourceAdaptiveSpinStages(disk,next);same(baseline,done.result);
    // Retained tokens are not consumed; the same saved branch yields identical results.
    disk.iterationStoreDirectory=(temp/"repeat").string();same(done.result,solveSourceAdaptiveSpinStages(disk,next).result);
    auto path=first.result.passes[0].solve.iterationStore->path();
    std::vector<SlqpIteration> records;
    {SlqpIterationReader reader(path);SlqpIteration record;
        while(reader.next(record))records.push_back(record);}
    require(!records.empty(),"trace tamper fixture empty");
    // Produce a structurally valid file with recomputed record SHA-256/footer.
    // activeVariables is deliberately outside the old scientific trace hash.
    records.front().activeVariables+=1;
    std::filesystem::rename(path,path+".original");
    {SlqpIterationStore replacement(path,disk.iterationStoreMaximumBytes);
        for(const auto& record:records)replacement.append(record);replacement.seal();}
    next.resume=first.checkpoint;disk.iterationStoreDirectory=(temp/"rehashed-corrupt-attempt").string();
    rejects([&]{solveSourceAdaptiveSpinStages(disk,next);},"rehashed full-record tamper accepted");
    require(!std::filesystem::exists(disk.iterationStoreDirectory),"rehashed corrupt resume created output");
    std::filesystem::rename(path,path+".tampered");
    std::filesystem::rename(path+".original",path);
    std::filesystem::resize_file(path,std::filesystem::file_size(path)-1);
    next.resume=first.checkpoint;disk.iterationStoreDirectory=(temp/"corrupt-attempt").string();
    rejects([&]{solveSourceAdaptiveSpinStages(disk,next);},"corrupt prior trace accepted");
    require(!std::filesystem::exists(disk.iterationStoreDirectory),"corrupt resume created new output first");
    require(boundaries>=5,"insufficient distinct stage boundaries");
    std::cout<<"native complete-history stage pause/resume passed; NOT cross-process resume\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#else
#include <iostream>
int main(){std::cerr<<"RED: native complete-history stage continuation is unavailable\n";return 1;}
#endif
