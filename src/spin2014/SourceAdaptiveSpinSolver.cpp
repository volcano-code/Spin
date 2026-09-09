#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "spin2014/SlqpAuditTrace.h"
#include "manufacturing/SlicerValidation.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace spin2014 {
const char* sourceAdaptiveSpinStatusToken(SourceAdaptiveSpinStatus s) noexcept {
    switch(s){case SourceAdaptiveSpinStatus::BinaryMomentValidated:return "BINARY_MOMENT_VALIDATED";case SourceAdaptiveSpinStatus::BinaryValidationFailed:return "BINARY_VALIDATION_FAILED";case SourceAdaptiveSpinStatus::RelaxationUnknown:return "RELAXATION_UNKNOWN";case SourceAdaptiveSpinStatus::AdaptivePassLimit:return "ADAPTIVE_PASS_LIMIT_UNKNOWN";case SourceAdaptiveSpinStatus::ResourceLimitReached:return "RESOURCE_LIMIT_REACHED_UNKNOWN";}return "UNKNOWN";
}
void validateSourceAdaptiveSpinRequest(const SourceAdaptiveSpinRequest& r) {
    if(!r.source||r.initialLevel<0||r.maximumOptimizationLevel<r.initialLevel||r.maximumOptimizationLevel>r.source->config().maximumDepth||r.maximumAdaptivePasses<=0||r.maximumLeafCount==0||r.maximumLeafCount>static_cast<std::size_t>(std::numeric_limits<int>::max()/6)||r.maximumOutputCuboids==0)throw std::invalid_argument("invalid source levels/resource guards");
    const double controls[]={r.initialFill,r.initialPhi,r.fractionalEpsilon,r.laplacianWeight,r.inertiaWeight,r.centerWeight,r.density,r.binaryConstraintTolerance};
    for(double c:controls)if(!std::isfinite(c))throw std::invalid_argument("nonfinite source request control");
    if(r.initialFill<0||r.initialFill>1||r.fractionalEpsilon<=0||r.fractionalEpsilon>=.5||r.laplacianWeight<0||r.inertiaWeight<0||r.centerWeight<0||r.density<=0||r.binaryConstraintTolerance<=0||(r.objective==Spin2014Objective::StaticLowering&&r.initialPhi!=0))throw std::invalid_argument("invalid source fill/weight/phi control");
    (void)SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
    for(int degree=3;degree<=5;++degree){double s=r.density*std::pow(r.source->config().rootSide,degree);if(!std::isfinite(s)||s<=0)throw std::invalid_argument("nonfinite/underflow fixed source scaling");}
    const auto& s=r.slqp;if(s.algorithm!=SlqpAlgorithm::AuditedKktRestorationV1)throw std::invalid_argument("source adaptive requires AUDITED_KKT_RESTORATION_V1");
    validateSlqpAuditSettings(s.audit);
    const double sc[]={s.equalityTolerance,s.projectedGradientTolerance,s.activeSetTolerance,s.meritPenalty,s.lineSearchReduction,s.armijo};for(double c:sc)if(!std::isfinite(c))throw std::invalid_argument("nonfinite source SLQP control");
    if(s.maximumIterations<=0||s.maximumLineSearchSteps<=0||s.lbfgsMemory<=0||s.equalityTolerance<=0||s.projectedGradientTolerance<=0||s.activeSetTolerance<0||s.activeSetTolerance>=.5||s.meritPenalty<=0||s.lineSearchReduction<=0||s.lineSearchReduction>=1||s.armijo<=0||s.armijo>=1)throw std::invalid_argument("invalid source SLQP control");
}
std::string computeSourceAdaptiveSpinRequestHash(const SourceAdaptiveSpinRequest& r) {
    validateSourceAdaptiveSpinRequest(r);std::ostringstream o;o.imbue(std::locale::classic());o<<std::setprecision(17)<<"source-adaptive-request-v1 PAPER_GIVENS AUDITED_KKT_RESTORATION_V1 IMPLEMENTATION_DEFAULT "
        <<r.source->partitionSha256()<<' '<<spin2014ObjectiveToken(r.objective)<<' '<<r.initialLevel<<' '<<r.maximumOptimizationLevel<<' '<<r.maximumAdaptivePasses<<' '<<r.maximumLeafCount<<' '<<r.maximumOutputCuboids<<' ';
    for(double c:{r.initialFill,r.initialPhi,r.fractionalEpsilon,r.laplacianWeight,r.inertiaWeight,r.centerWeight,r.density,r.binaryConstraintTolerance})o<<c<<',';
    writeSlqpAuditControlsJson(o,r.slqp);return manufacturing::computeContentSha256(o.str());
}
namespace {
void number(std::ostream& o,double v){if(std::isfinite(v))o<<v;else if(std::isnan(v))o<<"ABSENT_NONFINITE_NAN";else o<<(v<0?"ABSENT_NEG_INF":"ABSENT_POS_INF");o<<',';}
void vector(std::ostream& o,const Eigen::VectorXd& v){o<<v.size()<<'[';for(Eigen::Index i=0;i<v.size();++i)number(o,v(i));o<<']';}
void transition(std::ostream& o,const SourcePartitionTransition& t){o<<t.beforeHash<<' '<<t.afterHash<<' ';vector(o,t.removedBefore.values());vector(o,t.removedAfter.values());
    auto key=[&](const SourcePartitionKey& k){o<<k.level<<':'<<k.coordinates[0]<<','<<k.coordinates[1]<<','<<k.coordinates[2]<<';';};
    for(const auto& c:t.changes){key(c.parent);o<<'[';for(const auto& k:c.children)key(k);o<<']';}}
void finish(SourceAdaptiveSpinResult& r) {
    std::ostringstream trace;trace.imbue(std::locale::classic());trace<<std::hexfloat<<"source-adaptive-trace-v1 "<<r.inputHash;
    for(const auto& p:r.passes){trace<<p.index<<' '<<p.beforeHash<<' '<<p.afterHash<<' '<<p.solve.traceHash<<' ';vector(trace,p.inputBeta);vector(trace,p.outputBeta);number(trace,p.inputPhi);number(trace,p.solve.finalPhi);number(trace,p.solve.evaluation.objective);number(trace,p.solve.evaluation.paperObjective);number(trace,p.solve.evaluation.regularizationEnergy);vector(trace,p.solve.evaluation.finalMoments.values());transition(trace,p.split);transition(trace,p.merge);trace<<p.fractionalCount<<';';for(auto n:p.levelCounts)trace<<n<<',';trace<<';';for(auto n:p.fractionalLevelCounts)trace<<n<<',';}
    r.traceHash=manufacturing::computeContentSha256(trace.str());std::ostringstream o;o.imbue(std::locale::classic());o<<std::hexfloat<<"source-adaptive-output-v1 NUMERICAL_CUBOID_LIST_ONLY "<<r.inputHash<<' '<<r.traceHash<<' '<<sourceAdaptiveSpinStatusToken(r.status)<<' '<<std::quoted(r.diagnostic)<<' '<<r.finalPartitionHash<<' '<<r.geometryDepth<<','<<r.maximumOptimizationLevel<<';';number(o,r.finalPhi);vector(o,r.relaxedFill);vector(o,r.binaryFill);vector(o,r.roundingDelta);
    o<<"moments="<<r.binaryMoments.has_value()<<';';if(r.binaryMoments)vector(o,r.binaryMoments->values());
    o<<"cuboids="<<r.binaryCuboids.has_value()<<';';if(r.binaryCuboids)o<<r.binaryCuboids->outputHash;
    o<<"check="<<r.binaryCheck.has_value()<<';';if(r.binaryCheck){const auto& c=*r.binaryCheck;o<<c.passed<<' '<<std::quoted(c.diagnostic)<<' ';o<<c.independentMoments.has_value();if(c.independentMoments)vector(o,c.independentMoments->values());vector(o,c.constraints);number(o,c.normalizedMomentError);number(o,c.rawConstraintResidual);number(o,c.normalizedConstraintResidual);number(o,c.minimumCentralInertiaEigenvalue);}
    r.outputHash=manufacturing::computeContentSha256(o.str());
}
}
SourceAdaptiveSpinResult solveSourceAdaptiveSpin(const SourceAdaptiveSpinRequest& r) {
    validateSourceAdaptiveSpinRequest(r);SourceAdaptiveSpinResult out;out.inputHash=computeSourceAdaptiveSpinRequestHash(r);out.finalPhi=r.initialPhi;out.geometryDepth=r.source->config().maximumDepth;out.maximumOptimizationLevel=r.maximumOptimizationLevel;
    try {
        SourceAdaptiveDomain domain(r.source,r.initialLevel,r.maximumOptimizationLevel,r.maximumLeafCount);
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(domain.leaves().size(),r.initialFill);
        auto snapshot=[&]{out.relaxedFill=beta;out.finalPartition=domain.leaves();out.finalPartitionHash=domain.partitionHash();};snapshot();
        bool complete=false;
        const auto problem=SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight);
        for(int i=0;i<r.maximumAdaptivePasses;++i){SourceAdaptivePass p;p.index=i;p.before=domain.leaves();p.beforeHash=domain.partitionHash();p.inputBeta=beta;p.inputPhi=out.finalPhi;
            try {
                AdaptiveFillFunctional f(domain,problem,r.laplacianWeight,r.density);
                p.solve=solveSlqpLbfgsPaper(f,beta,out.finalPhi,r.slqp);
            }catch(const std::bad_alloc&){throw;}catch(const std::length_error&){throw;}
            catch(const std::exception& e){
                // The existing solver can reject an empty prefix or an initial
                // nonfinite evaluation before producing a trace. Keep the valid
                // partition/state and explicit invalid KKT; never infer success.
                p.solve.status=SlqpStatus::NumericalUnknown;p.solve.diagnostic=std::string("initial/source solve evaluation unavailable: ")+e.what();
                p.solve.beta=beta;p.solve.initialPhi=p.solve.finalPhi=out.finalPhi;
                p.solve.algorithm=r.slqp.algorithm;p.solve.formulation=Spin2014Formulation::PaperGivens;
                p.solve.boundedVariableCount=beta.size();p.solve.finalState.resize(beta.size()+(r.objective==Spin2014Objective::StaticLowering?0:1));
                p.solve.finalState.head(beta.size())=beta;if(p.solve.finalState.size()>beta.size())p.solve.finalState(beta.size())=out.finalPhi;
                p.solve.traceHash=hashAuditedSlqp(p.solve,r.slqp);
            }
            beta=p.solve.beta;out.finalPhi=p.solve.finalPhi;
            p.levelCounts.resize(r.maximumOptimizationLevel+1);p.fractionalLevelCounts.resize(r.maximumOptimizationLevel+1);
            for(std::size_t k=0;k<domain.leaves().size();++k){int level=domain.leaves()[k].level;++p.levelCounts[level];if(beta(k)>=r.fractionalEpsilon&&beta(k)<=1-r.fractionalEpsilon){++p.fractionalCount;++p.fractionalLevelCounts[level];}}
            bool resource=false;
            if(p.solve.status==SlqpStatus::Converged){try {
                // Copy before the pair of transitions makes split+merge atomic as
                // one adaptive update when either exceeds an allocation guard.
                auto nextDomain=domain;auto nextBeta=beta;
                auto split=nextDomain.splitFractionalLeaves(nextBeta,r.fractionalEpsilon);
                auto merge=nextDomain.mergeBinarySiblings(nextBeta,r.fractionalEpsilon);
                domain=std::move(nextDomain);beta.swap(nextBeta);p.split=std::move(split);p.merge=std::move(merge);
            }catch(const std::length_error& e){resource=true;out.diagnostic=e.what();}catch(const std::bad_alloc&){resource=true;out.diagnostic="RESOURCE_LIMIT: source adaptation allocation";}}
            p.after=domain.leaves();p.afterHash=domain.partitionHash();p.outputBeta=beta;
            complete=p.solve.status==SlqpStatus::Converged&&!resource&&p.beforeHash==p.afterHash;
            if(complete)for(int level=0;level<r.maximumOptimizationLevel;++level)if(p.fractionalLevelCounts[level])complete=false;
            const auto status=p.solve.status;out.passes.push_back(std::move(p));snapshot();
            if(resource){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;finish(out);return out;}
            if(status!=SlqpStatus::Converged){out.status=SourceAdaptiveSpinStatus::RelaxationUnknown;out.diagnostic="audited relaxation incomplete: "+out.passes.back().solve.diagnostic;finish(out);return out;}
            if(complete)break;
        }
        if(!complete){out.status=SourceAdaptiveSpinStatus::AdaptivePassLimit;out.diagnostic="adaptive pass limit; no unchanged-partition final KKT";finish(out);return out;}
        Eigen::VectorXd binary=(beta.array()>=.5).cast<double>();
        // Build into local output so a cap failure cannot expose a partial list.
        auto cuboids=sourceBinaryCuboids(domain,binary,r.maximumOutputCuboids);
        out.binaryFill=binary;out.roundingDelta=binary-beta;out.binaryCuboids=std::move(cuboids);
        // The binary raw moments do not need a positive-mass objective evaluation.
        out.binaryMoments=r.source->outerSurfaceMoments().scaledDensity(r.density)-domain.removedMoments(binary).scaledDensity(r.density);
        out.binaryCheck=checkSourceBinaryMoments(*r.source,*out.binaryCuboids,*out.binaryMoments,r.objective,out.finalPhi,r.density,r.binaryConstraintTolerance);
        out.status=out.binaryCheck->passed?SourceAdaptiveSpinStatus::BinaryMomentValidated:SourceAdaptiveSpinStatus::BinaryValidationFailed;out.diagnostic=out.binaryCheck->diagnostic;
    }catch(const std::length_error& e){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;out.diagnostic=e.what();}
    catch(const std::bad_alloc&){out.status=SourceAdaptiveSpinStatus::ResourceLimitReached;out.diagnostic="RESOURCE_LIMIT: source allocation";}
    finish(out);return out;
}
}
