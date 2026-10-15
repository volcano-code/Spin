#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "spin2014/DiagnosticStageTimer.h"
#include "SourceAdaptiveSpinEndpoint.h"
#include "spin2014/SlqpAuditTrace.h"
#include "spin2014/SlqpIterationStore.h"
#include "verification/Sha256Stream.h"
#include <filesystem>
#include "manufacturing/SlicerValidation.h"
#include <cmath>
#include <Eigen/LU>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace spin2014 {
const char* sourceAdaptiveTerminationToken(SourceAdaptiveTermination t) noexcept {
    switch(t) {
    case SourceAdaptiveTermination::UnchangedPartitionKktV1:return "UNCHANGED_PARTITION_KKT_V1";
    case SourceAdaptiveTermination::PaperPostTransitionV1:return "PAPER_POST_TRANSITION_V1";
    }
    return "UNKNOWN";
}
const char* sourceAdaptiveSpinStatusToken(SourceAdaptiveSpinStatus s) noexcept {
    if(s==SourceAdaptiveSpinStatus::StageBoundaryPaused)return "STAGE_BOUNDARY_PAUSED_NOT_VALIDATED";
    switch(s){case SourceAdaptiveSpinStatus::BinaryMomentValidated:return "BINARY_MOMENT_VALIDATED";case SourceAdaptiveSpinStatus::BinaryValidationFailed:return "BINARY_VALIDATION_FAILED";case SourceAdaptiveSpinStatus::InfeasibleCertified:return "INFEASIBLE_CERTIFIED";case SourceAdaptiveSpinStatus::RelaxationUnknown:return "RELAXATION_UNKNOWN";case SourceAdaptiveSpinStatus::AdaptivePassLimit:return "ADAPTIVE_PASS_LIMIT_UNKNOWN";case SourceAdaptiveSpinStatus::ResourceLimitReached:return "RESOURCE_LIMIT_REACHED_UNKNOWN";case SourceAdaptiveSpinStatus::StageBoundaryPaused:return "STAGE_BOUNDARY_PAUSED_NOT_VALIDATED";}return "UNKNOWN";
}
void validateSourceAdaptiveSpinRequest(const SourceAdaptiveSpinRequest& r) {
    if(r.termination!=SourceAdaptiveTermination::UnchangedPartitionKktV1 &&
       r.termination!=SourceAdaptiveTermination::PaperPostTransitionV1)
        throw std::invalid_argument("invalid source adaptive termination policy");
    if(!r.source||r.initialLevel<0||r.maximumOptimizationLevel<r.initialLevel||r.maximumOptimizationLevel>r.source->config().maximumDepth||r.maximumAdaptivePasses<=0||r.maximumLeafCount==0||r.maximumLeafCount>static_cast<std::size_t>(std::numeric_limits<int>::max()/6)||r.maximumOutputCuboids==0)throw std::invalid_argument("invalid source levels/resource guards");
    const double controls[]={r.initialFill,r.initialPhi,r.fractionalEpsilon,r.laplacianWeight,r.inertiaWeight,r.centerWeight,r.density,r.binaryConstraintTolerance};
    for(double c:controls)if(!std::isfinite(c))throw std::invalid_argument("nonfinite source request control");
    if(r.initialFill<0||r.initialFill>1||r.fractionalEpsilon<=0||r.fractionalEpsilon>=.5||r.laplacianWeight<0||r.inertiaWeight<0||r.centerWeight<0||r.density<=0||r.binaryConstraintTolerance<=0||(r.objective==Spin2014Objective::StaticLowering&&r.initialPhi!=0))throw std::invalid_argument("invalid source fill/weight/phi control");
    if(r.boundaryDensity.has_value()!=r.interiorDensity.has_value())
        throw std::invalid_argument("source dual-material request must specify both boundaryDensity and interiorDensity");
    for(double d:{r.resolvedBoundaryDensity(),r.resolvedInteriorDensity()}) {
        if(!std::isfinite(d)||d<=0)throw std::invalid_argument("source material densities must be finite and positive");
        for(int degree=3;degree<=5;++degree){double s=d*std::pow(r.source->config().rootSide,degree);if(!std::isfinite(s)||s<=0)throw std::invalid_argument("nonfinite/underflow source material scaling");}
    }
    if(!r.targetFrame.origin.allFinite()||!r.targetFrame.axes.allFinite() ||
       !(r.targetFrame.axes.transpose()*r.targetFrame.axes).isApprox(Eigen::Matrix3d::Identity(),1e-12) ||
       std::abs(r.targetFrame.axes.determinant()-1.)>1e-12)
        throw std::invalid_argument("source target frame must be finite, right-handed, and orthonormal");
    (void)SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
    const auto& s=r.slqp;if(s.algorithm!=SlqpAlgorithm::AuditedKktRestorationV1)throw std::invalid_argument("source adaptive requires AUDITED_KKT_RESTORATION_V1");
    validateSlqpAuditSettings(s.audit);validateSlqpDirectionSettings(s);
    const double sc[]={s.equalityTolerance,s.projectedGradientTolerance,s.activeSetTolerance,s.meritPenalty,s.lineSearchReduction,s.armijo};for(double c:sc)if(!std::isfinite(c))throw std::invalid_argument("nonfinite source SLQP control");
    if(s.maximumIterations<=0||s.maximumLineSearchSteps<=0||s.lbfgsMemory<=0||s.equalityTolerance<=0||s.projectedGradientTolerance<=0||s.activeSetTolerance<0||s.activeSetTolerance>=.5||s.meritPenalty<=0||s.lineSearchReduction<=0||s.lineSearchReduction>=1||s.armijo<=0||s.armijo>=1)throw std::invalid_argument("invalid source SLQP control");
}
std::string computeSourceAdaptiveSpinRequestHash(const SourceAdaptiveSpinRequest& r) {
    validateSourceAdaptiveSpinRequest(r);std::ostringstream o;o.imbue(std::locale::classic());o<<std::setprecision(17)<<"source-adaptive-request-v2 PAPER_GIVENS AUDITED_KKT_RESTORATION_V1 IMPLEMENTATION_DEFAULT "<<sourceAdaptiveTerminationToken(r.termination)<<' '
        <<r.source->partitionSha256()<<' '<<spin2014ObjectiveToken(r.objective)<<' '<<r.initialLevel<<' '<<r.maximumOptimizationLevel<<' '<<r.maximumAdaptivePasses<<' '<<r.maximumLeafCount<<' '<<r.maximumOutputCuboids<<' ';
    for(double c:{r.initialFill,r.initialPhi,r.fractionalEpsilon,r.laplacianWeight,r.inertiaWeight,r.centerWeight,r.density,r.binaryConstraintTolerance})o<<c<<',';
    // Preserve every historical single-density identity byte-for-byte. The
    // explicit dual-density contract receives a distinct request identity.
    if(r.boundaryDensity) o<<"source-dual-material-density-v1 "<<*r.boundaryDensity<<','<<*r.interiorDensity<<';';
    const auto identity=spin::SpinTargetFrame{};
    if(r.targetFrame.origin!=identity.origin || r.targetFrame.axes!=identity.axes) {
        o<<"target-frame-v1 ";
        for(int i=0;i<3;++i)o<<r.targetFrame.origin(i)<<',';
        for(int row=0;row<3;++row)for(int col=0;col<3;++col)o<<r.targetFrame.axes(row,col)<<',';
    }
    writeSlqpAuditControlsJson(o,r.slqp);return manufacturing::computeContentSha256(o.str());
}
namespace {
void number(std::ostream& o,double v){if(std::isfinite(v))o<<v;else if(std::isnan(v))o<<"ABSENT_NONFINITE_NAN";else o<<(v<0?"ABSENT_NEG_INF":"ABSENT_POS_INF");o<<',';}
void vector(std::ostream& o,const Eigen::VectorXd& v){o<<v.size()<<'[';for(Eigen::Index i=0;i<v.size();++i)number(o,v(i));o<<']';}
void transition(std::ostream& o,const SourcePartitionTransition& t){o<<t.beforeHash<<' '<<t.afterHash<<' ';vector(o,t.removedBefore.values());vector(o,t.removedAfter.values());
    auto key=[&](const SourcePartitionKey& k){o<<k.level<<':'<<k.coordinates[0]<<','<<k.coordinates[1]<<','<<k.coordinates[2]<<';';};
    for(const auto& c:t.changes){key(c.parent);o<<'[';for(const auto& k:c.children)key(k);o<<']';}}
std::optional<SourceAffineConstraintCertificate> affinePartitionCertificate(
    const AdaptiveFillFunctional& f, double phi, double tolerance) {
    // The first rows are center-of-mass and (for the paper Givens formulation)
    // transverse-moment equalities.  They are independent of phi.  The final
    // diagonalization row is not affine in phi and is therefore intentionally
    // excluded from this exact box certificate.
    const auto e=f.evaluatePaper(Eigen::VectorXd::Zero(f.betaCount()),phi);
    const Eigen::Index n=f.betaCount();
    const double certTolerance=std::max(1e-12,10.*tolerance);
    for(Eigen::Index row=0;row<e.constraints.size();++row) {
        if(f.paperHasAngle() && row==e.constraints.size()-1) continue;
        const auto a=e.constraintJacobian.row(row).head(n);
        if(!a.allFinite() || !std::isfinite(e.constraints(row))) continue;
        double lower=e.constraints(row),upper=e.constraints(row);
        for(Eigen::Index i=0;i<n;++i) {
            if(a(i)>=0.) upper+=a(i); else lower+=a(i);
        }
        if(!std::isfinite(lower)||!std::isfinite(upper)) continue;
        if(lower>certTolerance || upper<-certTolerance) {
            SourceAffineConstraintCertificate c;
            c.constraintIndex=static_cast<int>(row);c.lower=lower;c.upper=upper;c.tolerance=certTolerance;
            return c;
        }
    }
    return std::nullopt;
}
void finish(SourceAdaptiveSpinResult& r) {
    std::ostringstream trace;trace.imbue(std::locale::classic());trace<<std::hexfloat<<"source-adaptive-trace-v2 "<<r.inputHash;
    for(const auto& p:r.passes){trace<<p.index<<' '<<p.beforeHash<<' '<<p.afterHash<<' '<<p.solve.traceHash<<' ';vector(trace,p.inputBeta);vector(trace,p.outputBeta);number(trace,p.inputPhi);number(trace,p.solve.finalPhi);number(trace,p.solve.evaluation.objective);number(trace,p.solve.evaluation.paperObjective);number(trace,p.solve.evaluation.regularizationEnergy);vector(trace,p.solve.evaluation.finalMoments.values());transition(trace,p.split);transition(trace,p.merge);trace<<p.fractionalCount<<';';for(auto n:p.levelCounts)trace<<n<<',';trace<<';';for(auto n:p.fractionalLevelCounts)trace<<n<<',';}
    for(const auto& p:r.passes){trace<<"post="<<p.outputFractionalCount<<','<<p.paperResolutionComplete<<';';for(auto n:p.outputLevelCounts)trace<<n<<',';trace<<';';for(auto n:p.outputFractionalLevelCounts)trace<<n<<',';}
    r.traceHash=manufacturing::computeContentSha256(trace.str());std::ostringstream o;o.imbue(std::locale::classic());o<<std::hexfloat<<"source-adaptive-output-v2 NUMERICAL_CUBOID_LIST_ONLY "<<sourceAdaptiveTerminationToken(r.termination)<<' '<<r.terminalRelaxationKktApplies<<' '<<r.inputHash<<' '<<r.traceHash<<' '<<sourceAdaptiveSpinStatusToken(r.status)<<' '<<std::quoted(r.diagnostic)<<' '<<r.finalPartitionHash<<' '<<r.geometryDepth<<','<<r.maximumOptimizationLevel<<';';number(o,r.finalPhi);vector(o,r.relaxedFill);vector(o,r.binaryFill);vector(o,r.roundingDelta);
    o<<"moments="<<r.binaryMoments.has_value()<<';';if(r.binaryMoments)vector(o,r.binaryMoments->values());
    o<<"cuboids="<<r.binaryCuboids.has_value()<<';';if(r.binaryCuboids)o<<r.binaryCuboids->outputHash;
    o<<"check="<<r.binaryCheck.has_value()<<';';if(r.binaryCheck){const auto& c=*r.binaryCheck;o<<c.passed<<' '<<std::quoted(c.diagnostic)<<' ';o<<c.independentMoments.has_value();if(c.independentMoments)vector(o,c.independentMoments->values());vector(o,c.constraints);number(o,c.normalizedMomentError);number(o,c.rawConstraintResidual);number(o,c.normalizedConstraintResidual);number(o,c.minimumCentralInertiaEigenvalue);}
    if(r.infeasibilityCertificate) {
        const auto& c=*r.infeasibilityCertificate;
        o<<"affine_certificate="<<c.constraintIndex<<' ';number(o,c.lower);number(o,c.upper);number(o,c.tolerance);o<<std::quoted(c.scope)<<';';
    }
    r.outputHash=manufacturing::computeContentSha256(o.str());
}
}
void detail::finishSourceAdaptiveResult(SourceAdaptiveSpinResult& out) {finish(out);}
struct detail::SourceCheckpointAccess {
    static std::string fullTraceDigest(const SlqpLbfgsResult& solve) {
        verification::Sha256Stream hash;
        // The complete binary record includes fields not present in historical
        // scientific hashes. Do not confuse the base audited hash with the
        // Givens wrapper's additional final-state binding.
        forEachSlqpIteration(solve,[&](const SlqpIteration& record) {
            const auto bytes=encodeSlqpIteration(record);
            hash.stream().write(bytes.data(),bytes.size());
        });
        return hash.digest();
    }
    static SourceAdaptiveSpinResult restore(const SourceAdaptiveCheckpoint& saved,
        const SourceAdaptiveSpinRequest& r,std::uint64_t& bytes,
        std::vector<std::string>& validatedTraceDigests) {
        const auto& h=saved.history_;
        if(h.inputHash!=computeSourceAdaptiveSpinRequestHash(r) ||
           saved.maximumBytes_!=r.iterationStoreMaximumBytes ||
           saved.disk_!=!r.iterationStoreDirectory.empty())
            throw std::invalid_argument("checkpoint request/storage controls changed");
        if(h.passes.empty() || h.passes.size()>=std::size_t(r.maximumAdaptivePasses) ||
           saved.traceDigests_.size()!=h.passes.size())
            throw std::invalid_argument("checkpoint outside resumable stage range");
        std::uint64_t verified=0;
        validatedTraceDigests.clear();validatedTraceDigests.reserve(h.passes.size());
        for(std::size_t i=0;i<h.passes.size();++i) {
            const auto& p=h.passes[i];
            if(p.index!=int(i) || p.solve.status!=SlqpStatus::Converged ||
               saved.traceDigests_[i].size()!=64)
                throw SlqpStorageError("checkpoint prior stage/trace mismatch");
            if(saved.disk_) {
                if(!p.solve.iterationStore || !p.solve.iterationStore->sealed())
                    throw SlqpStorageError("checkpoint prior stage/trace mismatch");
                p.solve.iterationStore->validateUnchangedSnapshot();
                const auto n=p.solve.iterationStore->bytesWritten();
                if(n==0 || n>saved.maximumBytes_-verified)
                    throw SlqpStorageError("checkpoint prior trace size/budget mismatch");
                verified+=n;
            } else if(p.solve.iterationStore || slqpTraceCount(p.solve)==0) {
                throw SlqpStorageError("in-memory checkpoint prior trace mismatch");
            }
            validatedTraceDigests.push_back(saved.traceDigests_[i]);
        }
        if(verified!=saved.storedBytes_)throw SlqpStorageError("checkpoint cumulative trace mismatch");
        bytes=verified;return h;
    }
    static std::shared_ptr<const SourceAdaptiveCheckpoint> capture(
        const SourceAdaptiveSpinResult& out,std::uint64_t bytes,const SourceAdaptiveSpinRequest& r,
        const std::vector<std::string>& validatedPrefixDigests) {
        if(validatedPrefixDigests.size()>out.passes.size())
            throw SlqpStorageError("checkpoint validated trace prefix exceeds history");
        std::vector<std::string> digests=validatedPrefixDigests;digests.reserve(out.passes.size());
        for(std::size_t i=0;i<validatedPrefixDigests.size();++i) {
            const auto& store=out.passes[i].solve.iterationStore;
            if(!r.iterationStoreDirectory.empty()) {
                if(!store || !store->sealed())throw SlqpStorageError("checkpoint prefix lacks sealed history");
                store->validateUnchangedSnapshot();
            } else if(store || slqpTraceCount(out.passes[i].solve)==0) {
                throw SlqpStorageError("in-memory checkpoint prefix history mismatch");
            }
        }
        for(std::size_t i=validatedPrefixDigests.size();i<out.passes.size();++i)
            digests.push_back(fullTraceDigest(out.passes[i].solve));
        return std::shared_ptr<const SourceAdaptiveCheckpoint>(new SourceAdaptiveCheckpoint(
            out,bytes,r.iterationStoreMaximumBytes,!r.iterationStoreDirectory.empty(),std::move(digests)));
    }
};
SourceAdaptiveSpinResult detail::solveSourceAdaptiveSpinEndpoint(const SourceAdaptiveSpinRequest& r,
    SourceAdaptiveSpinEndpoint* endpoint,const SourceAdaptiveStageControl* control,
    std::shared_ptr<const SourceAdaptiveCheckpoint>* checkpoint) {
    validateSourceAdaptiveSpinRequest(r);SourceAdaptiveSpinResult out;out.termination=r.termination;out.inputHash=computeSourceAdaptiveSpinRequestHash(r);out.finalPhi=r.initialPhi;out.geometryDepth=r.source->config().maximumDepth;out.maximumOptimizationLevel=r.maximumOptimizationLevel;
    if(control && (control->maximumNewPasses<=0 || !checkpoint))
        throw std::invalid_argument("invalid native stage execution control");
    if(checkpoint)checkpoint->reset();
    if(!r.slqp.iterationStorePath.empty())throw std::invalid_argument("source solve requires iterationStoreDirectory, not a shared per-pass path");
    std::uint64_t storedBytes=0;std::vector<std::string> validatedTraceDigests;
    if(control && control->resume) {
        // Validate immutable earlier trace files before any new output creation.
        SPIN2014_DIAG_TIMER("endpoint_checkpoint_restore",-1);
        out=detail::SourceCheckpointAccess::restore(*control->resume,r,storedBytes,validatedTraceDigests);
        out.status=SourceAdaptiveSpinStatus::RelaxationUnknown;out.diagnostic.clear();
    }
    if(!r.iterationStoreDirectory.empty()) {
        std::error_code ec;
        if(!std::filesystem::create_directory(r.iterationStoreDirectory,ec))
            throw SlqpStorageError("SLQP iteration storage: exclusive directory creation failed: "+r.iterationStoreDirectory+": "+ec.message());
    }
    try {
        const int startPass=static_cast<int>(out.passes.size());
        auto makeDomain=[&] {
            if(!startPass)return SourceAdaptiveDomain(r.source,r.initialLevel,r.maximumOptimizationLevel,r.maximumLeafCount);
            std::vector<SourcePartitionKey> keys;keys.reserve(out.finalPartition.size());
            for(const auto& leaf:out.finalPartition)keys.push_back(leaf);
            return SourceAdaptiveDomain::restorePartition(r.source,r.initialLevel,
                r.maximumOptimizationLevel,r.maximumLeafCount,keys,out.finalPartitionHash);
        };
        SourceAdaptiveDomain domain=[&]{
            SPIN2014_DIAG_TIMER("endpoint_domain_restore",-1);
            return makeDomain();
        }();
        Eigen::VectorXd beta=startPass?out.relaxedFill:Eigen::VectorXd::Constant(domain.leaves().size(),r.initialFill);
        (void)domain.removedMoments(beta);
        auto snapshot=[&](int pass){SPIN2014_DIAG_TIMER("adaptive_snapshot",pass);
            out.relaxedFill=beta;out.finalPartition=domain.leaves();out.finalPartitionHash=domain.partitionHash();};snapshot(startPass-1);
        bool complete=false;
        const auto problem=SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
        for(int i=startPass;i<r.maximumAdaptivePasses;++i){SourceAdaptivePass p;
            {SPIN2014_DIAG_TIMER("adaptive_pass_input_snapshot",i);
                p.index=i;p.before=domain.leaves();p.beforeHash=domain.partitionHash();p.inputBeta=beta;p.inputPhi=out.finalPhi;}
            try {
                // Empty beta is not a source optimization, even when an angular
                // objective would give the shared solver a phi-only state.
                if(domain.leaves().empty())
                    throw std::invalid_argument("empty conservative source partition: no fill variables");
                AdaptiveFillFunctional f=[&]{SPIN2014_DIAG_TIMER("adaptive_functional_construction",i);
                    return AdaptiveFillFunctional(domain,problem,r.laplacianWeight,
                        r.resolvedBoundaryDensity(),r.resolvedInteriorDensity(),r.targetFrame);}();
                // Before invoking SLQP, certify the finite source partition's
                // affine equality box whenever a whole row excludes zero.
                // This is intentionally scoped to the current partition: a
                // later octree refinement may change the attainable interval.
                if(const auto certificate=affinePartitionCertificate(
                        f,out.finalPhi,r.slqp.equalityTolerance)) {
                    out.infeasibilityCertificate=certificate;
                    out.status=SourceAdaptiveSpinStatus::InfeasibleCertified;
                    std::ostringstream message;message.imbue(std::locale::classic());
                    message<<"discrete source partition affine constraint "
                           <<certificate->constraintIndex<<" interval ["
                           <<std::setprecision(17)<<certificate->lower<<", "
                           <<certificate->upper<<"] excludes zero; continuous/source-refinement problem not classified";
                    out.diagnostic=message.str();
                    finish(out);return out;
                }
                auto controls=r.slqp;
                if(!r.iterationStoreDirectory.empty()) {
                    std::ostringstream name;name<<"pass-"<<std::setw(6)<<std::setfill('0')<<i<<".slqptrace";
                    controls.iterationStorePath=(std::filesystem::path(r.iterationStoreDirectory)/name.str()).string();
                    controls.iterationStoreMaximumBytes=r.iterationStoreMaximumBytes-storedBytes;
                }
                {
                    SPIN2014_DIAG_TIMER("adaptive_slqp_solve",i);
                    p.solve=solveSlqpLbfgsPaper(f,beta,out.finalPhi,controls);
                }
                if(p.solve.iterationStore)storedBytes+=p.solve.iterationStore->bytesWritten();
            }catch(const std::bad_alloc&){throw;}catch(const std::length_error&){throw;}
            catch(const SlqpStorageError&){throw;}
            catch(const std::exception& e){
                // The source guard or an initial nonfinite solver evaluation
                // can reject the state before producing a trace. Keep the valid
                // partition/state and explicit invalid KKT; never infer success.
                p.solve.status=SlqpStatus::NumericalUnknown;p.solve.diagnostic=std::string("initial/source solve evaluation unavailable: ")+e.what();
                p.solve.beta=beta;p.solve.initialPhi=p.solve.finalPhi=out.finalPhi;
                p.solve.algorithm=r.slqp.algorithm;p.solve.formulation=Spin2014Formulation::PaperGivens;
                p.solve.boundedVariableCount=beta.size();p.solve.finalState.resize(beta.size()+(r.objective==Spin2014Objective::StaticLowering?0:1));
                p.solve.finalState.head(beta.size())=beta;if(p.solve.finalState.size()>beta.size())p.solve.finalState(beta.size())=out.finalPhi;
                p.solve.traceHash=hashAuditedSlqp(p.solve,r.slqp);
            }
            beta=p.solve.beta;out.finalPhi=p.solve.finalPhi;
            {
                SPIN2014_DIAG_TIMER("adaptive_input_level_accounting",i);
                p.levelCounts.resize(r.maximumOptimizationLevel+1);p.fractionalLevelCounts.resize(r.maximumOptimizationLevel+1);
                for(std::size_t k=0;k<domain.leaves().size();++k){int level=domain.leaves()[k].level;++p.levelCounts[level];if(beta(k)>=r.fractionalEpsilon&&beta(k)<=1-r.fractionalEpsilon){++p.fractionalCount;++p.fractionalLevelCounts[level];}}
            }
            bool resource=false;
            if(p.solve.status==SlqpStatus::Converged){try {
                // Copy before the pair of transitions makes split+merge atomic as
                // one adaptive update when either exceeds an allocation guard.
                auto nextDomain=domain;auto nextBeta=beta;
                SourcePartitionTransition split,merge;
                {SPIN2014_DIAG_TIMER("adaptive_split",i);split=nextDomain.splitFractionalLeaves(nextBeta,r.fractionalEpsilon);}
                {SPIN2014_DIAG_TIMER("adaptive_merge",i);merge=nextDomain.mergeBinarySiblings(nextBeta,r.fractionalEpsilon);}
                domain=std::move(nextDomain);beta.swap(nextBeta);p.split=std::move(split);p.merge=std::move(merge);
            }catch(const std::length_error& e){resource=true;out.diagnostic=e.what();}catch(const std::bad_alloc&){resource=true;out.diagnostic="RESOURCE_LIMIT: source adaptation allocation";}}
            p.after=domain.leaves();p.afterHash=domain.partitionHash();p.outputBeta=beta;
            {
                SPIN2014_DIAG_TIMER("adaptive_output_level_accounting",i);
                p.outputLevelCounts.resize(r.maximumOptimizationLevel+1);
                p.outputFractionalLevelCounts.resize(r.maximumOptimizationLevel+1);
                for(std::size_t k=0;k<p.after.size();++k) {
                    const int level=p.after[k].level;++p.outputLevelCounts[level];
                    if(beta(k)>=r.fractionalEpsilon && beta(k)<=1-r.fractionalEpsilon) {
                        ++p.outputFractionalCount;++p.outputFractionalLevelCounts[level];
                    }
                }
            }
            p.paperResolutionComplete=p.solve.status==SlqpStatus::Converged&&!resource;
            for(int level=0;level<r.maximumOptimizationLevel;++level)
                if(p.outputFractionalLevelCounts[level])p.paperResolutionComplete=false;
            out.terminalRelaxationKktApplies=p.solve.status==SlqpStatus::Converged &&
                p.solve.finalKkt.valid && p.beforeHash==p.afterHash && p.solve.beta==beta;
            complete=p.paperResolutionComplete;
            if(r.termination==SourceAdaptiveTermination::UnchangedPartitionKktV1)
                complete=complete&&out.terminalRelaxationKktApplies;
            const auto status=p.solve.status;out.passes.push_back(std::move(p));snapshot(i);
            if(resource){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;finish(out);return out;}
            if(status!=SlqpStatus::Converged){out.status=SourceAdaptiveSpinStatus::RelaxationUnknown;out.diagnostic="audited relaxation incomplete: "+out.passes.back().solve.diagnostic;finish(out);return out;}
            if(complete)break;
            if(control && i-startPass+1>=control->maximumNewPasses && i+1<r.maximumAdaptivePasses) {
                out.status=SourceAdaptiveSpinStatus::StageBoundaryPaused;
                out.diagnostic="execution paused after complete solve and split/merge; not binary validated";
                {SPIN2014_DIAG_TIMER("adaptive_finish_trace_hash",i);finish(out);}
                {SPIN2014_DIAG_TIMER("endpoint_checkpoint_capture",i);
                    *checkpoint=detail::SourceCheckpointAccess::capture(out,storedBytes,r,validatedTraceDigests);}
                return out;
            }
        }
        if(!complete){out.status=SourceAdaptiveSpinStatus::AdaptivePassLimit;out.diagnostic=std::string("adaptive pass limit under ")+sourceAdaptiveTerminationToken(r.termination);finish(out);return out;}
        Eigen::VectorXd binary=(beta.array()>=.5).cast<double>();
        // Build into local output so a cap failure cannot expose a partial list.
        auto cuboids=sourceBinaryCuboids(domain,binary,r.maximumOutputCuboids);
        out.binaryFill=binary;out.roundingDelta=binary-beta;out.binaryCuboids=std::move(cuboids);
        // The binary raw moments do not need a positive-mass objective evaluation.
        out.binaryMoments=(r.source->outerSurfaceMoments()-r.source->designSurfaceMoments()).scaledDensity(r.resolvedBoundaryDensity())
            +r.source->designSurfaceMoments().scaledDensity(r.resolvedInteriorDensity())
            -domain.removedMoments(binary).scaledDensity(r.resolvedInteriorDensity());
        out.binaryCheck=checkSourceBinaryMoments(*r.source,*out.binaryCuboids,*out.binaryMoments,r.objective,out.finalPhi,
            r.resolvedBoundaryDensity(),r.resolvedInteriorDensity(),r.binaryConstraintTolerance,r.targetFrame);
        out.status=out.binaryCheck->passed?SourceAdaptiveSpinStatus::BinaryMomentValidated:SourceAdaptiveSpinStatus::BinaryValidationFailed;out.diagnostic=out.binaryCheck->diagnostic;
        if(endpoint && out.status==SourceAdaptiveSpinStatus::BinaryMomentValidated) {
            endpoint->domain.emplace(std::move(domain));
            endpoint->binary=std::move(binary);
        }
    }catch(const std::length_error& e){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;out.diagnostic=e.what();}
    catch(const std::bad_alloc&){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;out.diagnostic="RESOURCE_LIMIT: source allocation";}
    finish(out);return out;
}
SourceAdaptiveSpinResult solveSourceAdaptiveSpin(const SourceAdaptiveSpinRequest& r) {
    return detail::solveSourceAdaptiveSpinEndpoint(r,nullptr);
}
SourceAdaptiveStageRun solveSourceAdaptiveSpinStages(const SourceAdaptiveSpinRequest& r,
    const SourceAdaptiveStageControl& control) {
    SourceAdaptiveStageRun run;
    run.result=detail::solveSourceAdaptiveSpinEndpoint(r,nullptr,&control,&run.checkpoint);
    return run;
}
}
