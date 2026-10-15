#include "spin2014/SourceCheckpointCodec.h"
#include "spin2014/DiagnosticStageTimer.h"
#include "spin2014/SlqpIterationStore.h"
#include "spin2014/SlqpAuditTrace.h"
#include "SourceAdaptiveSpinEndpoint.h"
#include "manufacturing/SlicerValidation.h"
#include "verification/Sha256Stream.h"
#include <cstring>
#include <chrono>
#include <filesystem>
#include <limits>

namespace spin2014 {
namespace {
constexpr std::uint64_t limit=64ULL*1024*1024;
const std::string magic="SCSPV001";
[[noreturn]] void fail(const char* why){throw SlqpStorageError(std::string("source checkpoint: ")+why);}
void require(bool ok,const char* why){if(!ok)fail(why);}
struct Writer {
    std::string bytes=magic;
    void u64(std::uint64_t n){require(bytes.size()<=limit-8,"metadata limit");for(int i=0;i<8;++i)bytes.push_back(char(n>>(8*i)));}
    void text(const std::string& s){u64(s.size());require(s.size()<=limit-bytes.size(),"metadata limit");bytes+=s;}
};
struct Reader {
    const std::string& bytes;std::size_t pos=8,end;
    explicit Reader(const std::string& b):bytes(b),end(b.size()>=64?b.size()-64:0) {
        require(b.size()<=limit && b.size()>=8+64 && b.compare(0,8,magic)==0,"metadata size/version");
        require(manufacturing::computeContentSha256(b.substr(0,end))==b.substr(end),"metadata checksum");
    }
    std::uint64_t u64(){require(pos<=end && end-pos>=8,"truncated metadata");std::uint64_t n=0;
        for(int i=0;i<8;++i)n|=std::uint64_t(static_cast<unsigned char>(bytes[pos++]))<<(8*i);return n;}
    std::string text(){auto n=u64();require(n<=end-pos,"invalid metadata string length");auto s=bytes.substr(pos,n);pos+=n;return s;}
};
std::string fullDigest(const SlqpLbfgsResult& solve) {
    verification::Sha256Stream hash;
    forEachSlqpIteration(solve,[&](const auto& t){auto b=encodeSlqpIteration(t);hash.stream().write(b.data(),b.size());});
    return hash.digest();
}
// A version-1 solution record uses only these semantically corresponding
// SlqpIteration fields. Canonical roundtrip rejects data in unused fields.
SlqpIteration solutionRecord(const SlqpLbfgsResult& s) {
    SlqpIteration t;t.finalState=true;t.state=s.finalState;t.phi=s.finalPhi;t.kkt=s.finalKkt;
    t.objective=s.evaluation.objective;t.gradient=s.evaluation.gradient;
    t.constraints=s.evaluation.constraints;t.constraintJacobian=s.evaluation.constraintJacobian;return t;
}
template<class A,class B> bool equal(const A& a,const B& b){return a.rows()==b.rows()&&a.cols()==b.cols()&&(a.array()==b.array()).all();}
void checkEvaluation(const SlqpIteration& t,const AdaptiveFillEvaluation& e) {
    require(t.objective==e.objective && equal(t.gradient,e.gradient) && equal(t.constraints,e.constraints) &&
        equal(t.constraintJacobian,e.constraintJacobian),"saved evaluation disagrees with source geometry");
}
std::string paperHashFromBase(std::string b,const SlqpLbfgsResult& s) {
    // Preserve the existing Givens wrapper's little-endian final-state suffix.
    auto number=[&](double x){std::uint64_t bits;std::memcpy(&bits,&x,8);for(int k=0;k<8;++k)b.push_back(char(bits>>(8*k)));};
    number(s.finalPhi);for(Eigen::Index i=0;i<s.beta.size();++i)number(s.beta[i]);
    return manufacturing::computeContentSha256(b);
}
std::string paperHash(const SlqpLbfgsResult& s,const SlqpLbfgsSettings& controls) {
    return paperHashFromBase(hashAuditedSlqp(s,controls),s);
}
std::string paperHash(const SlqpLbfgsResult& s,const SlqpLbfgsSettings& controls,
    const SlqpIterationTraversal& traverse) {
    return paperHashFromBase(hashAuditedSlqp(s,controls,traverse),s);
}
void counts(const std::vector<SourceAdaptiveLeaf>& leaves,const Eigen::VectorXd& beta,
    const SourceAdaptiveSpinRequest& r,std::vector<std::size_t>& levels,std::vector<std::size_t>& fractions,std::size_t& total) {
    levels.assign(r.maximumOptimizationLevel+1,0);fractions.assign(levels.size(),0);total=0;
    for(std::size_t i=0;i<leaves.size();++i){auto level=leaves[i].level;++levels[level];
        if(beta[i]>=r.fractionalEpsilon&&beta[i]<=1-r.fractionalEpsilon){++fractions[level];++total;}}
}
}
std::string encodeSourceAdaptiveCheckpoint(const SourceAdaptiveCheckpoint& saved) {
    require(saved.disk_,"durable codec requires streamed sealed histories");
    Writer w;w.text(saved.history_.inputHash);w.u64(saved.maximumBytes_);w.u64(saved.history_.passes.size());
    std::uint64_t bytes=0;
    for(std::size_t i=0;i<saved.history_.passes.size();++i) {
        const auto& p=saved.history_.passes[i];const auto& s=p.solve;
        require(s.iterationStore && s.iterationStore->sealed(),"missing sealed history");
        require(saved.traceDigests_.at(i).size()==64,"invalid cached history digest");
        s.iterationStore->validateUnchangedSnapshot();
        auto n=s.iterationStore->bytesWritten();
        require(n==s.iterationStore->bytesWritten() && n<=saved.maximumBytes_-bytes,"history bytes changed");bytes+=n;
        w.text(s.iterationStore->path());w.text(saved.traceDigests_[i]);w.text(s.traceHash);
        w.text(s.diagnostic);w.text(encodeSlqpIteration(solutionRecord(s)));
        w.text(p.beforeHash);w.text(p.afterHash);
    }
    require(bytes==saved.storedBytes_,"cumulative history size mismatch");
    w.text(saved.history_.traceHash);w.text(saved.history_.outputHash);w.u64(bytes);
    require(w.bytes.size()<=limit-64,"metadata limit");w.bytes+=manufacturing::computeContentSha256(w.bytes);return w.bytes;
}
struct DecodedCheckpoint {
    SourceAdaptiveSpinResult result;
    std::uint64_t storedBytes=0;
    std::vector<std::string> traceDigests;
};
DecodedCheckpoint decodeCheckpoint(const SourceAdaptiveSpinRequest& r,const std::string& bytes) {
    SPIN2014_DIAG_TIMER("checkpoint_decode_total",-1);
    validateSourceAdaptiveSpinRequest(r);Reader in(bytes);
    require(!r.iterationStoreDirectory.empty() && r.slqp.iterationStorePath.empty(),"decode requires source streaming mode");
    SourceAdaptiveSpinResult out;out.inputHash=computeSourceAdaptiveSpinRequestHash(r);
    require(in.text()==out.inputHash,"request/source identity changed");
    require(in.u64()==r.iterationStoreMaximumBytes,"trace budget changed");
    const auto count=in.u64();
    require(count>0 && count<std::uint64_t(r.maximumAdaptivePasses) && count<=(in.end-in.pos)/56,"invalid completed stage count");
    out.termination=r.termination;out.geometryDepth=r.source->config().maximumDepth;
    out.maximumOptimizationLevel=r.maximumOptimizationLevel;out.finalPhi=r.initialPhi;
    SourceAdaptiveDomain domain=[&]{SPIN2014_DIAG_TIMER("checkpoint_decode_initial_domain",-1);
        return SourceAdaptiveDomain(r.source,r.initialLevel,r.maximumOptimizationLevel,r.maximumLeafCount);}();
    Eigen::VectorXd beta=Eigen::VectorXd::Constant(domain.leaves().size(),r.initialFill);
    const bool angle=r.objective!=Spin2014Objective::StaticLowering;
    const auto problem=SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
    std::uint64_t used=0;std::vector<std::string> digests;
    for(std::uint64_t index=0;index<count;++index) {
        SPIN2014_DIAG_TIMER("checkpoint_decode_stage_total",static_cast<int>(index));
        SourceAdaptivePass p;p.index=int(index);p.before=domain.leaves();p.beforeHash=domain.partitionHash();p.inputBeta=beta;p.inputPhi=out.finalPhi;
        auto path=in.text(),digest=in.text(),solveHash=in.text(),diagnostic=in.text(),encodedSolution=in.text();
        const auto beforeHash=in.text(),afterHash=in.text();
        require(!path.empty() && path.find('\0')==std::string::npos && std::filesystem::path(path).is_absolute(),"invalid history path");
        require(beforeHash==p.beforeHash,"history input partition mismatch");
        auto& s=p.solve;
        s.status=SlqpStatus::Converged;s.diagnostic=diagnostic;
        s.algorithm=r.slqp.algorithm;s.formulation=Spin2014Formulation::PaperGivens;
        s.boundedVariableCount=beta.size();s.initialPhi=out.finalPhi;
        auto summary=decodeSlqpIteration(encodedSolution);s.finalState=summary.state;s.finalPhi=summary.phi;
        require(s.finalState.size()==beta.size()+(angle?1:0) && s.finalState.allFinite(),"invalid solution dimensions/nonfinite state");
        s.beta=s.finalState.head(beta.size());(void)domain.removedMoments(s.beta);
        require(s.finalPhi==(angle?s.finalState[beta.size()]:0.),"solution phi mismatch");
        AdaptiveFillFunctional f(domain,problem,r.laplacianWeight,r.density);
        s.evaluation=f.evaluatePaper(s.beta,s.finalPhi);checkEvaluation(summary,s.evaluation);
        // Re-evaluate the saved dual witness, never accept its boolean flag.
        const auto& k=summary.kkt;
        const auto actual=evaluateKkt(s.evaluation,s.finalState,beta.size(),k.lambda,k.alpha,k.eta);
        require(k.status=="WITNESS_LP_OPTIMAL_FIRST_ORDER" && k.valid && actual.valid && actual.equalityPrimal==k.equalityPrimal && actual.boundPrimal==k.boundPrimal &&
            actual.stationarity==k.stationarity && actual.dualSign==k.dualSign && actual.complementarity==k.complementarity &&
            actual.full==k.full,"invalid final KKT witness");
        require(actual.equalityPrimal<=r.slqp.equalityTolerance && actual.boundPrimal<=r.slqp.projectedGradientTolerance &&
            actual.stationarity<=r.slqp.projectedGradientTolerance && actual.dualSign<=r.slqp.projectedGradientTolerance &&
            actual.complementarity<=r.slqp.projectedGradientTolerance,"saved relaxation not converged");
        s.finalKkt=k;
        require(encodeSlqpIteration(solutionRecord(s))==encodedSolution,"noncanonical solution record");
        Eigen::VectorXd expected(beta.size()+(angle?1:0));expected.head(beta.size())=beta;if(angle)expected[beta.size()]=out.finalPhi;
        verification::Sha256Stream hash;std::size_t records=0;SlqpIteration last;
        const SlqpIterationTraversal traverse=[&](const auto& consume) {
            {SPIN2014_DIAG_TIMER("checkpoint_decode_open_and_audit_trace",static_cast<int>(index));
            s.iterationStore=SlqpIterationStore::openSealedAndVisit(path,r.iterationStoreMaximumBytes-used,
                [&](const SlqpIteration& t) {
                require(t.iteration==int(records) && records<std::size_t(r.slqp.maximumIterations) &&
                    t.state.size()==expected.size() && t.state.allFinite(),"invalid historical iteration/state");
                require(!records || !last.finalState,"records after terminal iteration");
                if(!records)require(equal(t.state,expected),"history does not begin at inherited state");
                (void)domain.removedMoments(t.state.head(beta.size()));
                require(t.phi==(angle?t.state[beta.size()]:0.),"historical phi mismatch");
                checkEvaluation(t,f.evaluatePaper(t.state.head(beta.size()),t.phi));
                const auto b=encodeSlqpIteration(t);hash.stream().write(b.data(),b.size());
                consume(t);last=t;++records;
            });
            used+=s.iterationStore->bytesWritten();
            }
            require(records>0 && records==s.iterationStore->count() && last.finalState && equal(last.state,s.finalState),"history terminal state mismatch");
            SlqpIteration terminalKkt,summaryKkt;terminalKkt.kkt=last.kkt;summaryKkt.kkt=s.finalKkt;
            require(encodeSlqpIteration(terminalKkt)==encodeSlqpIteration(summaryKkt),"history terminal KKT mismatch");
            require(hash.digest()==digest,"complete history digest mismatch");
        };
        s.traceHash=paperHash(s,r.slqp,traverse);require(s.traceHash==solveHash,"historical solve hash mismatch");
        beta=s.beta;out.finalPhi=s.finalPhi;
        counts(p.before,beta,r,p.levelCounts,p.fractionalLevelCounts,p.fractionalCount);
        {SPIN2014_DIAG_TIMER("checkpoint_decode_split",static_cast<int>(index));
            p.split=domain.splitFractionalLeaves(beta,r.fractionalEpsilon);}
        {SPIN2014_DIAG_TIMER("checkpoint_decode_merge",static_cast<int>(index));
            p.merge=domain.mergeBinarySiblings(beta,r.fractionalEpsilon);}
        p.after=domain.leaves();p.afterHash=domain.partitionHash();p.outputBeta=beta;
        require(p.afterHash==afterHash,"replayed transition partition mismatch");
        counts(p.after,beta,r,p.outputLevelCounts,p.outputFractionalLevelCounts,p.outputFractionalCount);
        p.paperResolutionComplete=true;
        for(int level=0;level<r.maximumOptimizationLevel;++level)
            if(p.outputFractionalLevelCounts[level])p.paperResolutionComplete=false;
        out.terminalRelaxationKktApplies=p.beforeHash==p.afterHash && s.beta==beta;
        const bool complete=p.paperResolutionComplete && (r.termination==SourceAdaptiveTermination::PaperPostTransitionV1 || out.terminalRelaxationKktApplies);
        require(!complete,"terminal stage cannot be a resumable checkpoint");
        out.passes.push_back(std::move(p));digests.push_back(digest);
    }
    const auto traceHash=in.text(),outputHash=in.text();require(in.u64()==used,"cumulative trace counter mismatch");
    require(in.pos==in.end,"trailing metadata");
    out.relaxedFill=beta;out.finalPartition=domain.leaves();out.finalPartitionHash=domain.partitionHash();
    out.status=SourceAdaptiveSpinStatus::StageBoundaryPaused;
    out.diagnostic="execution paused after complete solve and split/merge; not binary validated";
    detail::finishSourceAdaptiveResult(out);
    require(out.traceHash==traceHash && out.outputHash==outputHash,"reconstructed complete history hash mismatch");
    return {std::move(out),used,std::move(digests)};
}
std::shared_ptr<const SourceAdaptiveCheckpoint> decodeSourceAdaptiveCheckpoint(
    const SourceAdaptiveSpinRequest& r,const std::string& bytes) {
    auto decoded=decodeCheckpoint(r,bytes);
    return std::shared_ptr<const SourceAdaptiveCheckpoint>(new SourceAdaptiveCheckpoint(
        decoded.result,decoded.storedBytes,r.iterationStoreMaximumBytes,true,std::move(decoded.traceDigests)));
}

SourceAdaptivePostSolveReplay replaySourceAdaptivePostSolve(
    const SourceAdaptiveSpinRequest& r,const std::string& materialCheckpointBytes,
    const std::string& expectedMaterialSettingsHash,const std::string& tracePath,int expectedPrePassCount) {
    validateSourceAdaptiveSpinRequest(r);
    require(expectedMaterialSettingsHash.size()==64,"material settings hash size");
    require(materialCheckpointBytes.size()>136 && materialCheckpointBytes.compare(0,8,"SCMAT001")==0 &&
        materialCheckpointBytes.substr(8,64)==expectedMaterialSettingsHash &&
        manufacturing::computeContentSha256(materialCheckpointBytes.substr(0,materialCheckpointBytes.size()-64))==
            materialCheckpointBytes.substr(materialCheckpointBytes.size()-64),"material checkpoint settings/version mismatch");
    SourceAdaptivePostSolveReplay replay;
    replay.checkpointSha256=manufacturing::computeContentSha256(materialCheckpointBytes);
    const auto began=std::chrono::steady_clock::now();
    auto decoded=decodeCheckpoint(r,materialCheckpointBytes.substr(72,materialCheckpointBytes.size()-136));
    replay.timingsNanoseconds["checkpoint_decode_and_history_validation"]=
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-began).count();
    auto& out=replay.result=std::move(decoded.result);
    const int passIndex=static_cast<int>(out.passes.size());
    require(expectedPrePassCount>0 && passIndex==expectedPrePassCount && passIndex<r.maximumAdaptivePasses,
        "checkpoint is not the required pre-pass stage boundary");
    const auto requested=std::filesystem::weakly_canonical(std::filesystem::path(tracePath));
    for(const auto& prior:out.passes) {
        require(bool(prior.solve.iterationStore),"checkpoint history lacks immutable trace store");
        require(std::filesystem::weakly_canonical(std::filesystem::path(prior.solve.iterationStore->path()))!=requested,
            "trace is already a committed stage boundary");
    }
    std::vector<SourcePartitionKey> keys;keys.reserve(out.finalPartition.size());
    for(const auto& leaf:out.finalPartition)keys.push_back(leaf);
    SourceAdaptiveDomain domain=SourceAdaptiveDomain::restorePartition(r.source,r.initialLevel,
        r.maximumOptimizationLevel,r.maximumLeafCount,keys,out.finalPartitionHash);
    Eigen::VectorXd beta=out.relaxedFill;const bool angle=r.objective!=Spin2014Objective::StaticLowering;
    require(beta.size()==Eigen::Index(domain.leaves().size()) && beta.allFinite(),"checkpoint inherited beta invalid");
    SourceAdaptivePass pass;pass.index=passIndex;pass.before=domain.leaves();pass.beforeHash=domain.partitionHash();
    pass.inputBeta=beta;pass.inputPhi=out.finalPhi;
    const auto problem=SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
    AdaptiveFillFunctional functional(domain,problem,r.laplacianWeight,r.density);
    Eigen::VectorXd expected(beta.size()+(angle?1:0));expected.head(beta.size())=beta;if(angle)expected[beta.size()]=out.finalPhi;
    SlqpIteration last;std::size_t records=0;verification::Sha256Stream raw;
    const auto traceBegan=std::chrono::steady_clock::now();
    const auto store=SlqpIterationStore::openSealedAndVisit(tracePath,r.iterationStoreMaximumBytes-decoded.storedBytes,
        [&](const SlqpIteration& record) {
            require(record.iteration==int(records) && records<std::size_t(r.slqp.maximumIterations),
                "invalid post-solve iteration index/count");
            require(!records || !last.finalState,"records after terminal post-solve state");
            require(record.state.size()==expected.size() && record.state.allFinite(),"invalid post-solve state");
            if(!records)require(equal(record.state,expected),"post-solve trace does not begin at stage boundary state");
            require(record.phi==(angle?record.state[beta.size()]:0.),"post-solve phi mismatch");
            (void)domain.removedMoments(record.state.head(beta.size()));
            checkEvaluation(record,functional.evaluatePaper(record.state.head(beta.size()),record.phi));
            const auto bytes=encodeSlqpIteration(record);raw.stream().write(bytes.data(),bytes.size());
            last=record;++records;
        });
    replay.timingsNanoseconds["trace_digest_capture_and_validation"]=
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-traceBegan).count();
    require(records>0 && records==store->count() && last.finalState,"post-solve trace has no terminal final-state marker");
    const auto terminalEvaluation=functional.evaluatePaper(last.state.head(beta.size()),last.phi);
    checkEvaluation(last,terminalEvaluation);
    const auto terminalKkt=evaluateKkt(terminalEvaluation,last.state,beta.size(),last.kkt.lambda,last.kkt.alpha,last.kkt.eta);
    require(last.kkt.valid && terminalKkt.valid && last.kkt.status=="WITNESS_LP_OPTIMAL_FIRST_ORDER" &&
        terminalKkt.equalityPrimal==last.kkt.equalityPrimal && terminalKkt.boundPrimal==last.kkt.boundPrimal &&
        terminalKkt.stationarity==last.kkt.stationarity && terminalKkt.dualSign==last.kkt.dualSign &&
        terminalKkt.complementarity==last.kkt.complementarity && terminalKkt.full==last.kkt.full,
        "post-solve terminal KKT witness mismatch");
    require(terminalKkt.equalityPrimal<=r.slqp.equalityTolerance && terminalKkt.boundPrimal<=r.slqp.projectedGradientTolerance &&
        terminalKkt.stationarity<=r.slqp.projectedGradientTolerance && terminalKkt.dualSign<=r.slqp.projectedGradientTolerance &&
        terminalKkt.complementarity<=r.slqp.projectedGradientTolerance,"post-solve terminal KKT gate failed");
    auto& solve=pass.solve;solve.status=SlqpStatus::Converged;
    solve.diagnostic="FIRST_ORDER_KKT_TOLERANCES_MET; no minimum or physical certificate";
    solve.algorithm=r.slqp.algorithm;solve.formulation=Spin2014Formulation::PaperGivens;solve.boundedVariableCount=beta.size();
    solve.initialPhi=out.finalPhi;solve.finalState=last.state;solve.finalPhi=last.phi;solve.beta=last.state.head(beta.size());
    solve.evaluation=terminalEvaluation;solve.finalKkt=last.kkt;solve.iterationStore=store;
    const SlqpIterationTraversal traversal=[&](const auto& consume){forEachSlqpIteration(solve,consume);};
    const auto traceHashBegan=std::chrono::steady_clock::now();
    solve.traceHash=paperHash(solve,r.slqp,traversal);
    replay.timingsNanoseconds["trace_paper_hash_second_traversal"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-traceHashBegan).count();
    replay.traceRecordStreamSha256=raw.digest();replay.terminalState=last.state;replay.terminalKkt=last.kkt;
    beta=solve.beta;out.finalPhi=solve.finalPhi;
    const auto inputBegan=std::chrono::steady_clock::now();
    counts(pass.before,beta,r,pass.levelCounts,pass.fractionalLevelCounts,pass.fractionalCount);
    replay.timingsNanoseconds["input_level_accounting"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-inputBegan).count();
    const auto splitBegan=std::chrono::steady_clock::now();pass.split=domain.splitFractionalLeaves(beta,r.fractionalEpsilon);
    replay.timingsNanoseconds["split"]=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-splitBegan).count();
    const auto mergeBegan=std::chrono::steady_clock::now();pass.merge=domain.mergeBinarySiblings(beta,r.fractionalEpsilon);
    replay.timingsNanoseconds["merge"]=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-mergeBegan).count();
    pass.after=domain.leaves();pass.afterHash=domain.partitionHash();pass.outputBeta=beta;
    const auto outputBegan=std::chrono::steady_clock::now();
    counts(pass.after,beta,r,pass.outputLevelCounts,pass.outputFractionalLevelCounts,pass.outputFractionalCount);
    replay.timingsNanoseconds["output_level_accounting"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-outputBegan).count();
    pass.paperResolutionComplete=true;for(int level=0;level<r.maximumOptimizationLevel;++level)
        if(pass.outputFractionalLevelCounts[level])pass.paperResolutionComplete=false;
    out.terminalRelaxationKktApplies=pass.beforeHash==pass.afterHash && solve.beta==beta;
    out.passes.push_back(std::move(pass));out.relaxedFill=beta;out.finalPartition=domain.leaves();out.finalPartitionHash=domain.partitionHash();
    out.status=SourceAdaptiveSpinStatus::StageBoundaryPaused;
    out.diagnostic="SUPPLEMENTAL_DIAGNOSTIC_ONLY; binary validation not run";
    const auto resultHashBegan=std::chrono::steady_clock::now();detail::finishSourceAdaptiveResult(out);
    replay.timingsNanoseconds["scientific_result_hash"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-resultHashBegan).count();
    const auto diagnosticCaptureBegan=std::chrono::steady_clock::now();
    std::vector<std::string> traceDigests;traceDigests.reserve(out.passes.size());
    for(const auto& captured:out.passes)traceDigests.push_back(fullDigest(captured.solve));
    const SourceAdaptiveCheckpoint checkpoint(out,decoded.storedBytes+store->bytesWritten(),r.iterationStoreMaximumBytes,true,std::move(traceDigests));
    replay.timingsNanoseconds["checkpoint_capture_in_memory"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-diagnosticCaptureBegan).count();
    const auto encodeBegan=std::chrono::steady_clock::now();
    const auto encoded=encodeSourceAdaptiveCheckpoint(checkpoint);
    replay.encodedCheckpointSha256=manufacturing::computeContentSha256(encoded);
    replay.timingsNanoseconds["checkpoint_encode_in_memory"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-encodeBegan).count();
    const auto captureBegan=std::chrono::steady_clock::now();
    replay.diagnosticBytes="POSTSOLVE_DIAGNOSTIC_ONLY\n"+replay.checkpointSha256+"\n"+replay.traceRecordStreamSha256+"\n"+out.outputHash+"\n";
    replay.timingsNanoseconds["diagnostic_byte_capture"]=std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-captureBegan).count();
    return replay;
}
}
