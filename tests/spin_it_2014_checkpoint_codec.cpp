#if __has_include("spin2014/SourceCheckpointCodec.h")
#include "spin2014/SourceCheckpointCodec.h"
#include "spin2014/SlqpIterationStore.h"
#include "spin2014/SlqpAuditTrace.h"
#include "manufacturing/SlicerValidation.h"
#include "verification/Sha256Stream.h"
#include "../src/spin2014/SourceAdaptiveSpinEndpoint.h"
#include "source_surface_test_utils.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <vector>
using namespace spin2014;
using namespace surface_test;
SourceAdaptiveSpinRequest fixture(int id) {
    const bool angle=id>=2;
    auto outer=angle?cube(-.5,.5):box({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    auto inner=angle?cube(-.25,.25):box({.125,-.125,.625},{.375,.125,.875});
    SourceSurfaceOctreeConfig c;c.maximumDepth=5;c.rootOrigin=angle?Eigen::Vector3d(-1,-1,-1):Eigen::Vector3d(-1,-1,0);c.rootSide=2;
    SourceAdaptiveSpinRequest r;r.source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,inner.v,inner.f,c);
    r.initialLevel=2;r.maximumOptimizationLevel=angle?4:5;r.maximumAdaptivePasses=8;r.laplacianWeight=1;
    r.slqp.maximumIterations=300;r.slqp.equalityTolerance=1e-10;r.slqp.projectedGradientTolerance=1e-8;
    r.slqp.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;
    if(angle){r.objective=Spin2014Objective::YoYo;r.initialPhi=.37;}
    if(id%2)r.termination=SourceAdaptiveTermination::PaperPostTransitionV1;
    return r;
}
void science(const SourceAdaptiveSpinResult& r) {
    std::cout<<'F'<<r.outputHash<<'\n'<<r.traceHash<<'\n'<<r.passes.size()<<'\n';
    for(const auto& p:r.passes)forEachSlqpIteration(p.solve,[&](const auto& t){
        const auto b=encodeSlqpIteration(t);std::cout<<b.size()<<':';std::cout.write(b.data(),b.size());});
}
// Malformed producer for negative tests only. Recompute every enclosing hash,
// so checksum or scientific-hash mismatch cannot mask a missing semantic gate.
void forgedKkt(SourceAdaptiveSpinRequest r,int kind) {
    auto run=solveSourceAdaptiveSpinStages(r,{});require(bool(run.checkpoint),"forgery needs a paused real solve");
    auto& s=run.result.passes[0].solve;
    if(kind) {
        std::vector<SlqpIteration> records;forEachSlqpIteration(s,[&](const auto& t){records.push_back(t);});
        if(kind==1){records.back().kkt.valid=false;records.back().kkt.status="FORGED_TERMINAL_KKT";}
        else {
            const auto last=records.back();while(records.size()<=std::size_t(r.slqp.maximumIterations))records.push_back(last);
            for(std::size_t i=0;i<records.size();++i){records[i].iteration=int(i);records[i].finalState=i+1==records.size();}
        }
        auto replacement=std::make_shared<SlqpIterationStore>(r.iterationStoreDirectory+"/forged.slqptrace",r.iterationStoreMaximumBytes);
        for(const auto& record:records)replacement->append(record);replacement->seal();s.iterationStore=replacement;
    }else s.finalKkt.status="FORGED_SUMMARY_KKT";
    auto digest=hashAuditedSlqp(s,r.slqp);
    auto addU64=[](std::string& b,std::uint64_t n){for(int i=0;i<8;++i)b.push_back(char(n>>(8*i)));};
    auto number=[&](double x){std::uint64_t n;std::memcpy(&n,&x,8);addU64(digest,n);};
    number(s.finalPhi);for(Eigen::Index i=0;i<s.beta.size();++i)number(s.beta[i]);
    s.traceHash=manufacturing::computeContentSha256(digest);detail::finishSourceAdaptiveResult(run.result);
    std::string wire="SCSPV001";
    auto text=[&](const std::string& b){addU64(wire,b.size());wire+=b;};
    text(run.result.inputHash);addU64(wire,r.iterationStoreMaximumBytes);addU64(wire,run.result.passes.size());
    std::uint64_t size=0;
    for(const auto& p:run.result.passes){const auto& q=p.solve;
        verification::Sha256Stream hash;forEachSlqpIteration(q,[&](const auto& t){auto b=encodeSlqpIteration(t);hash.stream().write(b.data(),b.size());});
        SlqpIteration summary;summary.finalState=true;summary.state=q.finalState;summary.phi=q.finalPhi;summary.kkt=q.finalKkt;
        summary.objective=q.evaluation.objective;summary.gradient=q.evaluation.gradient;summary.constraints=q.evaluation.constraints;summary.constraintJacobian=q.evaluation.constraintJacobian;
        text(q.iterationStore->path());text(hash.digest());text(q.traceHash);text(q.diagnostic);text(encodeSlqpIteration(summary));text(p.beforeHash);text(p.afterHash);
        size+=q.iterationStore->bytesWritten();
    }
    text(run.result.traceHash);text(run.result.outputHash);addU64(wire,size);wire+=manufacturing::computeContentSha256(wire);
    bool rejected=false;
    try{(void)decodeSourceAdaptiveCheckpoint(r,wire);}
    catch(const std::exception& e){require(std::string(e.what()).find(kind==2?"historical iteration":"KKT")!=std::string::npos,"forgery test failed before semantic guard");rejected=true;}
    require(rejected,"forged semantic checkpoint accepted");std::cout<<"KKT forgery rejected\n";
}
int main(int argc,char** argv){try {
    require(argc>=4,"mode, fixture, new directory [checkpoint] required");
    std::string mode=argv[1];int id=std::stoi(argv[2]);require(id>=0&&id<4,"fixture range");
    auto r=fixture(id);r.iterationStoreDirectory=std::filesystem::absolute(argv[3]).string();
    if(mode=="forge-summary" || mode=="forge-terminal" || mode=="forge-count"){
        forgedKkt(r,mode=="forge-summary"?0:mode=="forge-terminal"?1:2);return 0;}
    if(mode=="reopen") {
        const auto cap=argc==5?std::stoull(argv[4]):r.iterationStoreMaximumBytes;
        auto store=SlqpIterationStore::openSealed(argv[3],cap);
        require(store->sealed() && store->count()>0 && store->retainedRecords()==0 && store->peakRetainedRecords()==0,"read-only reload not streamed");
        rejects([&]{store->append(SlqpIteration{});},"read-only append allowed");
        rejects([&]{store->back();},"read-only mutable record exposed");
        rejects([&]{store->seal();},"read-only seal allowed");
        std::cout<<"read-only verified\n";return 0;
    }
    if(mode=="visit") {
        std::filesystem::create_directories(argv[3]);
        const auto path=(std::filesystem::path(argv[3])/"visitor.slqptrace").string();
        std::vector<SlqpIteration> sampleRecords;
        {
            SlqpIterationStore writer(path,r.iterationStoreMaximumBytes);
            for(int i=0;i<3;++i) {SlqpIteration record;record.iteration=i;record.objective=double(i+1);sampleRecords.push_back(record);writer.append(record);}
            writer.seal();
        }
        std::vector<int> visited;
        auto store=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes,
            [&](const SlqpIteration& record){visited.push_back(record.iteration);});
        require(store->sealed() && store->count()==3 && visited==std::vector<int>({0,1,2}),
            "single-pass sealed visitor did not validate and visit every record in order");
        store->validateUnchangedSnapshot();
        SlqpLbfgsResult inMemory;inMemory.trace=sampleRecords;inMemory.status=SlqpStatus::Converged;
        const auto referenceTraceHash=hashAuditedSlqp(inMemory,r.slqp);
        const auto streamedTraceHash=hashAuditedSlqp(inMemory,r.slqp,
            [&](const std::function<void(const SlqpIteration&)>& consume) {
                const auto audited=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes,consume);
                require(audited->count()==sampleRecords.size(),"streamed hash visitor record count mismatch");
            });
        require(streamedTraceHash==referenceTraceHash,"streamed audited trace hash changed science identity");
        const auto damaged=(std::filesystem::path(argv[3])/"damaged.slqptrace").string();
        std::filesystem::copy_file(path,damaged);
        std::string bytes;
        {std::ifstream input(damaged,std::ios::binary);bytes.assign(std::istreambuf_iterator<char>(input),{});}
        require(bytes.size()>96,"trace fixture too short");bytes[96]^=1;
        {std::ofstream output(damaged,std::ios::binary|std::ios::trunc);output.write(bytes.data(),bytes.size());}
        std::size_t damagedVisits=0;bool rejected=false;
        try { (void)SlqpIterationStore::openSealedAndVisit(damaged,r.iterationStoreMaximumBytes,
            [&](const SlqpIteration&){++damagedVisits;}); }
        catch(const SlqpStorageError&) { rejected=true; }
        require(rejected && damagedVisits==0,"single-pass visitor exposed a checksum-invalid record");
        auto readBytes=[](const std::filesystem::path& source) {
            std::ifstream input(source,std::ios::binary);
            return std::string(std::istreambuf_iterator<char>(input),{});
        };
        auto writeBytes=[](const std::filesystem::path& destination,const std::string& data) {
            std::ofstream output(destination,std::ios::binary|std::ios::trunc);
            output.write(data.data(),static_cast<std::streamsize>(data.size()));
        };
        auto u64At=[](const std::string& data,std::size_t offset) {
            std::uint64_t value=0;for(int i=0;i<8;++i)value|=std::uint64_t(static_cast<unsigned char>(data.at(offset+i)))<<(8*i);return value;
        };
        auto putU64=[](std::string& data,std::size_t offset,std::uint64_t value) {
            for(int i=0;i<8;++i)data.at(offset+i)=char(value>>(8*i));
        };
        const auto validBytes=readBytes(path);
        require(validBytes.size()>168,"multi-record sealed trace fixture too short");
        const auto secondFrame=std::size_t(24+72+u64At(validBytes,24));
        const auto secondPayload=secondFrame+72;
        auto expectBadStore=[&](const std::string& name,const std::string& contents,std::size_t expectedVisits) {
            const auto badPath=(std::filesystem::path(argv[3])/name).string();writeBytes(badPath,contents);
            std::size_t visits=0;bool failedClosed=false;
            try{(void)SlqpIterationStore::openSealedAndVisit(badPath,r.iterationStoreMaximumBytes,
                [&](const SlqpIteration&){++visits;});}catch(const SlqpStorageError&){failedClosed=true;}
            require(failedClosed && visits==expectedVisits,"malformed sealed trace was accepted or exposed an invalid record");
        };
        auto laterCorrupt=validBytes;laterCorrupt.at(secondPayload)^=1;
        expectBadStore("damaged-later.slqptrace",laterCorrupt,1);
        auto badFooterMagic=validBytes;badFooterMagic.at(badFooterMagic.size()-24)^=1;
        expectBadStore("damaged-footer-magic.slqptrace",badFooterMagic,0);
        auto badFooterEnd=validBytes;putU64(badFooterEnd,badFooterEnd.size()-8,u64At(validBytes,validBytes.size()-8)-1);
        expectBadStore("damaged-footer-end.slqptrace",badFooterEnd,0);
        auto badFooterCount=validBytes;putU64(badFooterCount,badFooterCount.size()-16,2);
        expectBadStore("damaged-footer-count.slqptrace",badFooterCount,1);
        auto trailing=validBytes;trailing.push_back('x');
        expectBadStore("trailing-byte.slqptrace",trailing,0);
        auto writeVariant=[&](const std::filesystem::path& variant,double offset) {
            SlqpIterationStore writer(variant.string(),r.iterationStoreMaximumBytes);
            for(int i=0;i<3;++i){SlqpIteration record;record.iteration=i;record.objective=double(i+1)+offset;writer.append(record);}
            writer.seal();
        };
        const auto inPlacePath=std::filesystem::path(argv[3])/"in-place-rewrite.slqptrace";
        writeVariant(inPlacePath,0.0);
        auto inPlaceStore=SlqpIterationStore::openSealed(inPlacePath.string(),r.iterationStoreMaximumBytes);
        auto inPlaceBytes=readBytes(inPlacePath);auto alteredRecord=sampleRecords.front();alteredRecord.objective+=100.0;
        const auto alteredPayload=encodeSlqpIteration(alteredRecord);
        require(alteredPayload.size()==u64At(inPlaceBytes,24),"in-place rewrite fixture changed record length");
        inPlaceBytes.replace(96,alteredPayload.size(),alteredPayload);
        inPlaceBytes.replace(32,64,manufacturing::computeContentSha256(alteredPayload));
        const auto originalMtime=std::filesystem::last_write_time(inPlacePath);
        writeBytes(inPlacePath,inPlaceBytes);std::filesystem::last_write_time(inPlacePath,originalMtime);
        require(std::filesystem::last_write_time(inPlacePath)==originalMtime,"test failed to restore in-place file mtime");
        {SlqpIterationReader validRewrite(inPlacePath.string());SlqpIteration record;std::size_t count=0;while(validRewrite.next(record))++count;
            require(count==sampleRecords.size(),"in-place rewrite is not a valid self-consistent sealed trace");}
        bool inPlaceRejected=false;
        try{(void)inPlaceStore->openReader();}catch(const SlqpStorageError&){inPlaceRejected=true;}
        require(inPlaceRejected,"same-inode rewrite with recomputed digest and restored mtime bypassed snapshot checks");
        bool inPlaceSnapshotRejected=false;
        try{inPlaceStore->validateUnchangedSnapshot();}catch(const SlqpStorageError&){inPlaceSnapshotRejected=true;}
        require(inPlaceSnapshotRejected,
            "snapshot-only checkpoint validation accepted a rewritten sealed trace");
        const auto variant=std::filesystem::path(argv[3])/"same-size-variant.slqptrace";
        const auto backup=std::filesystem::path(argv[3])/"original-backup.slqptrace";
        writeVariant(variant,10.0);
        require(std::filesystem::file_size(path)==std::filesystem::file_size(variant),
            "path-replacement fixture must preserve sealed-file length");
        bool changedDuringVisit=false,returnedAfterSwap=false;
        try {
            const auto swapped=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes,
                [&](const SlqpIteration&) {
                    if(changedDuringVisit)return;
                    std::filesystem::rename(path,backup);std::filesystem::rename(variant,path);changedDuringVisit=true;
                });
            returnedAfterSwap=bool(swapped);
        } catch(const SlqpStorageError&) {}
        if(changedDuringVisit){std::filesystem::remove(path);std::filesystem::rename(backup,path);}
        require(changedDuringVisit && !returnedAfterSwap,
            "sealed visitor returned a store after its pathname was replaced by same-size valid content");
        writeVariant(variant,20.0);
        auto pinned=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes,[](const SlqpIteration&){});
        auto pinnedReader=pinned->openReader();
        const auto afterOpenBackup=std::filesystem::path(argv[3])/"after-open-backup.slqptrace";
        std::filesystem::rename(path,afterOpenBackup);std::filesystem::rename(variant,path);
        std::vector<double> pinnedObjectives;
        SlqpIteration pinnedRecord;bool pinnedReaderRejected=false;
        try{while(pinnedReader->next(pinnedRecord))pinnedObjectives.push_back(pinnedRecord.objective);}
        catch(const SlqpStorageError&){pinnedReaderRejected=true;}
        SlqpLbfgsResult pinnedResult;pinnedResult.iterationStore=pinned;bool replacementTraversalRejected=false;
        try{forEachSlqpIteration(pinnedResult,[](const SlqpIteration&){});}catch(const SlqpStorageError&){replacementTraversalRejected=true;}
        std::filesystem::remove(path);std::filesystem::rename(afterOpenBackup,path);
        require(pinnedReaderRejected && replacementTraversalRejected && pinnedObjectives.size()<=3 &&
            std::equal(pinnedObjectives.begin(),pinnedObjectives.end(),std::vector<double>({1.,2.,3.}).begin()),
            "read-only store followed a replaced path or accepted a changed descriptor snapshot");
        std::cout<<"single-pass sealed trace audit passed\n";return 0;
    }
    if(mode=="scanbench") {
        require(argc==6 || argc==7,"scanbench requires output directory, records, trials per record, optional repeats");
        const auto records=std::stoi(argv[4]),trialsPerRecord=std::stoi(argv[5]);
        require(records>=2 && records<=64 && trialsPerRecord>=1000 && trialsPerRecord<=100000,
            "scanbench fixture dimensions out of bounds");
        const int repeats=argc>6?std::stoi(argv[6]):7;
        require(repeats>=5 && repeats<=30,"scanbench repeat count out of bounds");
        const auto root=std::filesystem::path(argv[3]);std::filesystem::create_directories(root);
        const auto path=(root/"large-synthetic.slqptrace").string();
        verification::Sha256Stream expectedRaw;
        {
            SlqpIterationStore writer(path,r.iterationStoreMaximumBytes);
            for(int recordIndex=0;recordIndex<records;++recordIndex) {
                SlqpIteration record;record.iteration=recordIndex;record.objective=1.0/(recordIndex+1);
                record.objectiveTrials.reserve(trialsPerRecord);
                for(int i=0;i<trialsPerRecord;++i) {
                    RestorationTrial trial;trial.step=1.0/(i+1);trial.constraintL1=double((i+recordIndex)%997)/997.0;
                    trial.valid=(i%11)!=0;trial.accepted=(i%17)==0;trial.reason=trial.accepted?"ACCEPTED":"INSUFFICIENT_MERIT_DECREASE";
                    record.objectiveTrials.push_back(std::move(trial));
                }
                const auto encoded=encodeSlqpIteration(record);expectedRaw.stream().write(encoded.data(),encoded.size());
                writer.append(std::move(record));
            }
            writer.seal();
        }
        const auto expectedRawHash=expectedRaw.digest();
        SlqpLbfgsResult result;result.status=SlqpStatus::Converged;result.diagnostic="SYNTHETIC_SCAN_ONLY";
        result.formulation=Spin2014Formulation::LegacyAngleEliminated;
        result.finalKkt.valid=true;result.finalKkt.status="SYNTHETIC";
        const auto fileBytes=std::filesystem::file_size(path);
        const auto scanBaseline=[&]() {
            auto checked=SlqpIterationStore::openSealed(path,r.iterationStoreMaximumBytes);
            require(checked->count()==std::size_t(records),"baseline sealed count mismatch");
            verification::Sha256Stream raw;std::size_t rawCount=0;
            SlqpIterationReader rawReader(path);SlqpIteration item;
            while(rawReader.next(item)){const auto bytes=encodeSlqpIteration(item);raw.stream().write(bytes.data(),bytes.size());++rawCount;}
            require(rawCount==std::size_t(records) && raw.digest()==expectedRawHash,"baseline raw trace digest mismatch");
            const SlqpIterationTraversal traverse=[&](const auto& consume) {
                SlqpIterationReader reader(path);SlqpIteration record;std::size_t count=0;
                while(reader.next(record)){consume(record);++count;}
                require(count==std::size_t(records),"baseline paper-hash count mismatch");
            };
            return hashAuditedSlqp(result,r.slqp,traverse);
        };
        const auto scanCandidate=[&]() {
            verification::Sha256Stream raw;std::size_t count=0;
            const SlqpIterationTraversal traverse=[&](const auto& consume) {
                const auto store=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes,
                    [&](const SlqpIteration& record) {
                        const auto bytes=encodeSlqpIteration(record);raw.stream().write(bytes.data(),bytes.size());
                        consume(record);++count;
                    });
                require(store->count()==std::size_t(records),"candidate single-pass count mismatch");
            };
            const auto paperHash=hashAuditedSlqp(result,r.slqp,traverse);
            require(count==std::size_t(records) && raw.digest()==expectedRawHash,"candidate raw trace digest mismatch");
            return paperHash;
        };
        std::vector<double> baselineMs,candidateMs;std::string baselineHash,candidateHash;
        for(int repetition=0;repetition<repeats;++repetition) {
            const bool baselineFirst=(repetition%2)==0;
            for(int slot=0;slot<2;++slot) {
                const bool runBaseline=(slot==0)==baselineFirst;
                const auto started=std::chrono::steady_clock::now();
                const auto digest=runBaseline?scanBaseline():scanCandidate();
                const double elapsed=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
                auto& samples=runBaseline?baselineMs:candidateMs;auto& expected=runBaseline?baselineHash:candidateHash;
                if(expected.empty())expected=digest;else require(expected==digest,"repeated scan changed scientific hash");
                samples.push_back(elapsed);
                std::cout<<(runBaseline?"baseline":"candidate")<<','<<repetition<<','<<std::fixed<<std::setprecision(3)<<elapsed<<" ms\n";
            }
            require(baselineHash.empty() || candidateHash.empty() || baselineHash==candidateHash,
                "baseline and candidate paper hashes differ");
        }
        auto median=[](std::vector<double> values){std::sort(values.begin(),values.end());return values[values.size()/2];};
        const double baselineMedian=median(baselineMs),candidateMedian=median(candidateMs);
        require(baselineHash==candidateHash,"baseline/candidate science hash mismatch");
        std::cout<<"{"<<"\"scope\":\"SYNTHETIC_TRACE_SCAN_MICROBENCH_NOT_G7\","
            <<"\"trace_bytes\":"<<fileBytes<<",\"records\":"<<records<<",\"trials_per_record\":"<<trialsPerRecord<<","
            <<"\"repeats\":"<<repeats<<",\"baseline_median_ms\":"<<baselineMedian<<",\"candidate_median_ms\":"<<candidateMedian<<","
            <<"\"candidate_change_percent\":"<<100.0*(candidateMedian/baselineMedian-1.0)<<","
            <<"\"raw_trace_sha256\":\""<<expectedRawHash<<"\",\"paper_trace_sha256\":\""<<baselineHash<<"\"}\n";
        return 0;
    }
    if(mode=="baseline"){science(solveSourceAdaptiveSpin(r));return 0;}
    SourceAdaptiveStageControl control;
    if(mode=="resume" || mode=="roundtrip") {
        require(argc==5,"snapshot required");
        require(std::filesystem::file_size(argv[4])<=64ULL*1024*1024,"input metadata size");
        std::ifstream file(argv[4],std::ios::binary);std::string bytes((std::istreambuf_iterator<char>(file)),{});
        control.resume=decodeSourceAdaptiveCheckpoint(r,bytes);
        require(control.resume->completedPasses()>0 && control.resume->storedTraceBytes()>0,"history or disk budget lost");
        if(mode=="roundtrip") {
            require(encodeSourceAdaptiveCheckpoint(*control.resume)==bytes,"noncanonical state roundtrip");
            std::cout<<"roundtrip verified\n";return 0;
        }
    }else require(mode=="start","mode invalid");
    auto run=solveSourceAdaptiveSpinStages(r,control);
    if(run.checkpoint){const auto b=encodeSourceAdaptiveCheckpoint(*run.checkpoint);std::cout<<'P';std::cout.write(b.data(),b.size());}
    else science(run.result);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#else
#include <iostream>
int main(){std::cerr<<"RED: production source-stage checkpoint codec unavailable\n";return 1;}
#endif
