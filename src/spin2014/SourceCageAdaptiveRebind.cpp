#include "spin2014/SourceCageAdaptiveRebind.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace spin2014 {
namespace {
using Stage=SourceCageAdaptiveRebindStage;
using Status=SourceCageAdaptiveRebindStageOutcome;
using Outcome=SourceCageAdaptiveRebindOutcome;
[[noreturn]] void resource(const char* why) {
    throw std::length_error(std::string("RESOURCE_LIMIT: source adaptive rebind ")+why);
}
std::size_t add(std::size_t a,std::size_t b) {if(a>SIZE_MAX-b)resource("addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b) {if(b&&a>SIZE_MAX/b)resource("multiplication overflow");return a*b;}
class Hash {
    verification::Sha256 h_;
public:
    void integer(std::uint64_t n) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);h_.update(b,8);}
    void string(const std::string& s) {integer(s.size());h_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double d) {std::uint64_t n;std::memcpy(&n,&d,8);integer(n);}
    template<class T> void values(const T& x) {integer(x.size());for(Eigen::Index i=0;i<x.size();++i)real(x[i]);}
    void key(const SourcePartitionKey& k) {integer(k.level);for(auto c:k.coordinates)integer(c);}
    std::string finish() {return h_.finish();}
};
std::string stateHash(const Eigen::VectorXd& x) {Hash h;h.string("source-adaptive-full-state-v1");h.values(x);return h.finish();}
std::string functionalHash(const SourceCagePaperFunctional& f) {
    Hash h;h.string("source-adaptive-functional-v1");
    for(const auto* s:{&f.configSha256(),&f.upstreamSha256(),&f.source().configSha256(),
        &f.source().upstreamSha256(),&f.source().boundary().configSha256(),&f.source().boundary().outputSha256(),
        &f.source().boundary().binding().configSha256(),&f.source().boundary().binding().outputSha256(),
        &f.source().boundary().binding().domain().partitionHash()})h.string(*s);
    return h.finish();
}
std::string configHash(const SourceCageAdaptiveRebindSettings& s) {
    Hash h;h.string("source-adaptive-rebind-v1 accounting-v1 fractional-epsilon IMPLEMENTATION_DEFAULT=.05");
    h.real(s.fractionalEpsilon);h.integer(s.budget.maximumOwnedBytes);
    const auto& g=s.geometry;h.real(g.minimumOuterDomainClearance);
    for(auto n:{g.maximumPoints,g.maximumTriangles,g.maximumPrimitives,g.maximumPairs,g.maximumRecords,
        g.maximumCheckerBytes,g.maximumNefComplexity,g.maximumConnectedShells,g.maximumExactStringBytes})h.integer(n);
    return h.finish();
}
void transitionHash(Hash& h,const SourcePartitionTransition& t) {
    h.string(t.beforeHash);h.string(t.afterHash);h.values(t.removedBefore.values());h.values(t.removedAfter.values());
    h.integer(t.changes.size());for(const auto& c:t.changes) {
        h.key(c.parent);h.integer(c.children.size());for(const auto& child:c.children)h.key(child);
    }
}
void materialHash(Hash& h,const SourceCageAdaptiveMaterialFacts& m) {
    h.integer(m.available);h.values(m.worldMoments.values());h.values(m.spinMoments.values());
    h.string(m.outputSha256);h.string(m.worldOutputSha256);
}
void finish(SourceCageAdaptiveRebindFacts& f) {
    Hash h;h.string("actual-source-adaptive-rebind-scientific-v1 REST-transitions; fresh-deformed-material; no-solver-conclusion");
    h.integer(unsigned(f.outcome));h.integer(unsigned(f.lastStage));h.integer(f.plannedOwnedBytes);
    for(const auto* s:{&f.configSha256,&f.oldStateSha256,&f.newStateSha256,&f.oldFunctionalSha256,&f.newFunctionalSha256,
        &f.oldPartitionSha256,&f.newPartitionSha256,&f.oldBindingSha256,&f.newBindingSha256})h.string(*s);
    for(const auto& r:f.stages) {h.integer(unsigned(r.outcome));h.string(r.reason.data());h.string(r.reasonSha256.data());}
    for(const auto* leaves:{&f.oldLeaves,&f.newLeaves}) {h.integer(leaves->size());for(const auto& k:*leaves)h.key(k);}
    h.values(f.oldBeta);h.values(f.newBeta);transitionHash(h,f.split);transitionHash(h,f.merge);
    for(const auto* audit:{&f.inputGeometry,&f.outputGeometry}) {
        h.integer(bool(*audit));if(*audit)h.string((*audit)->facts().scientificSha256);
    }
    materialHash(h,f.beforeMaterial);materialHash(h,f.afterMaterial);
    h.integer(f.deformedMomentDifferenceAvailable);h.values(f.deformedMomentDifference);f.scientificSha256=h.finish();
}
template<std::size_t N> void copy(std::array<char,N>& out,const char* text) {
    const auto n=std::min(N-1,std::strlen(text));std::memcpy(out.data(),text,n);out[n]=0;
}
void fail(SourceCageAdaptiveRebindFacts& f,Outcome o,const char* why) {
    f.outcome=o;auto& r=f.stages[unsigned(f.lastStage)];
    r.outcome=o==Outcome::ResourceLimit?Status::ResourceLimit:
        o==Outcome::RejectedGeometry?Status::RejectedGeometry:Status::Unknown;
    copy(r.reason,why);Hash h;h.string(why);copy(r.reasonSha256,h.finish().c_str());
}
// Each adaptive split is one level (at most 8 descendants), and each primitive
// can split into at most 8 Full cuboids. Merge chains are bounded by source depth.
// Two vector payloads cover geometric capacity growth of transition evidence.
// Mutable domain copies/rebuilt immutable layers retain their upstream caps.
std::size_t evidencePlan(const SourceCagePaperFunctional& f,const SourceCageAdaptiveRebindSettings& s) {
    const auto b=f.betaCount(),next=mul(b,8),p=f.source().boundary().binding().primitives().size();
    const auto changes=mul(next,add(std::size_t(f.source().boundary().binding().source().config().maximumDepth),1));
    std::size_t bytes=add(sizeof(SourceCageAdaptiveRebindFacts),4096); // shared controls and bounded digest strings
    bytes=add(bytes,mul(add(b,next),sizeof(SourcePartitionKey)+sizeof(double)));
    bytes=add(bytes,mul(mul(changes,2),add(sizeof(SourcePartitionChange),mul(8,sizeof(SourcePartitionKey)))));
    // local beta and checkpoint state; moves avoid duplicated full state storage.
    bytes=add(bytes,mul(add(next,add(next,add(f.hasAngle()?1:0,mul(f.handleCount(),4)))),sizeof(double)));
    for(auto primitives:{p,mul(p,8)}) {
        const auto pairs=primitives?mul(primitives/2,primitives-1+(primitives%2)):0;
        const auto records=add(3,add(mul(primitives,3),pairs));
        bytes=add(bytes,mul(records,mul(2,sizeof(SourceCageGeometryRecord))));
        bytes=add(bytes,mul(primitives,add(s.geometry.maximumExactStringBytes,1)));
        bytes=add(bytes,4096); // audit object, fixed metadata and digest strings
    }
    return bytes;
}
// Full evaluation temporaries, including upstream moment and joint arrays,
// decomposition, deformed points, per-handle rows and shape residual/gradient.
std::size_t evaluationPlan(const SourceCagePaperFunctional& f) {
    const auto n=f.stateCount(),b=f.betaCount(),h=f.handleCount();
    auto values=mul(n,add(f.constraintCount(),11));
    values=add(values,add(mul(b,10),mul(h,40)));
    values=add(values,add(mul(n,4),add(mul(h,10),f.constraintCount())));
    values=add(values,mul(std::size_t(f.source().boundary().points().rows()),3));
    if(f.shapeFunctional())values=add(values,add(mul(f.source().boundary().outerPointCount(),6),mul(h,4)));
    return add(mul(values,sizeof(double)),4096);
}
void account(SourceCageAdaptiveRebindFacts& out,std::size_t bytes) {
    out.plannedOwnedBytes=std::max(out.plannedOwnedBytes,bytes);
    if(bytes>out.settings.budget.maximumOwnedBytes)resource("owned evidence/state/evaluation budget exceeded");
}
SourceCageAdaptiveMaterialFacts material(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto e=f.evaluate(x);SourceCageAdaptiveMaterialFacts m;
    if(!e.world.finalMoments.isFinite()||!e.spinMoments.isFinite()||!e.jointGradient.allFinite()||
        !e.jointConstraintJacobian.allFinite()||!e.spinMomentJacobian.allFinite())
        throw std::overflow_error("nonfinite rebuilt material or derivatives");
    m.worldMoments=e.world.finalMoments;m.spinMoments=e.spinMoments;
    m.outputSha256=e.outputSha256;m.worldOutputSha256=e.world.outputSha256;m.available=true;return m;
}
void snapshot(const SourceAdaptiveDomain& d,std::vector<SourcePartitionKey>& keys) {
    keys.reserve(d.leaves().size());for(const auto& l:d.leaves())keys.push_back(l);
}
bool geometryPass(SourceCageAdaptiveRebindFacts& f,const SourceCageGeometryAudit& a) {
    switch(a.facts().outcome) {
    case SourceCageGeometryOutcome::ChecksPassed:return true;
    case SourceCageGeometryOutcome::RejectedGeometry:fail(f,Outcome::RejectedGeometry,"actual represented geometry rejected");return false;
    case SourceCageGeometryOutcome::ResourceLimit:fail(f,Outcome::ResourceLimit,"actual geometry checker resource limit");return false;
    default:fail(f,Outcome::Unknown,"actual geometry checker unknown");return false;
    }
}
}
struct SourceCageAdaptiveRebindAudit::Data {
    SourceCageAdaptiveRebindFacts facts;
    std::shared_ptr<const SourceCageAdaptiveCheckpoint> checkpoint;
};
const SourceCageAdaptiveRebindFacts& SourceCageAdaptiveRebindAudit::facts() const {
    if(!data_)throw std::logic_error("moved-from adaptive rebind audit");return data_->facts;
}
std::shared_ptr<const SourceCageAdaptiveCheckpoint> SourceCageAdaptiveRebindAudit::checkpoint() const {
    if(!data_)throw std::logic_error("moved-from adaptive rebind audit");return data_->checkpoint;
}
SourceCageAdaptiveRebindAudit rebindSourceCageAdaptive(const SourceCagePaperFunctional& original,
    const Eigen::VectorXd& x,const SourceCageAdaptiveRebindSettings& settings) {
    if(std::size_t(x.size())!=original.stateCount()||!x.allFinite())
        throw std::invalid_argument("adaptive rebind requires finite actual full state dimensions");
    if((x.head(original.betaCount()).array()<0).any()||(x.head(original.betaCount()).array()>1).any())
        throw std::invalid_argument("adaptive rebind beta must be in [0,1]");
    for(std::size_t h=0;h<original.handleCount();++h)if(x[original.handleOffset()+4*h]<=0)
        throw std::invalid_argument("adaptive rebind handle scales must be positive");
    if(!std::isfinite(settings.fractionalEpsilon)||settings.fractionalEpsilon<=0||settings.fractionalEpsilon>=.5)
        throw std::invalid_argument("adaptive rebind epsilon must be in (0,.5)");
    if(!settings.budget.maximumOwnedBytes||settings.budget.maximumOwnedBytes>std::size_t(PTRDIFF_MAX))
        throw std::invalid_argument("invalid adaptive rebind owned byte cap");
    // Constructor validates checker API settings before any transaction work.
    const SourceCageGeometryChecker inputChecker(original,settings.geometry);
    auto data=std::make_shared<SourceCageAdaptiveRebindAudit::Data>();auto& out=data->facts;
    out.settings=settings;out.configSha256=configHash(settings);out.oldStateSha256=stateHash(x);
    out.oldFunctionalSha256=functionalHash(original);
    const auto& binding=original.source().boundary().binding();
    out.oldPartitionSha256=binding.domain().partitionHash();out.oldBindingSha256=binding.outputSha256();
    auto begin=[&](Stage stage) {out.lastStage=stage;out.stages[unsigned(stage)].outcome=Status::Attempted;};
    auto passed=[&] {out.stages[unsigned(out.lastStage)].outcome=Status::Passed;};
    auto run=[&] {
        begin(Stage::Preflight);const auto evidence=evidencePlan(original,settings);
        account(out,add(evidence,evaluationPlan(original)));
        snapshot(binding.domain(),out.oldLeaves);out.oldBeta=x.head(original.betaCount());passed();
        begin(Stage::InputGeometry);
        out.inputGeometry=std::make_shared<const SourceCageGeometryAudit>(inputChecker.check(x));
        if(!geometryPass(out,*out.inputGeometry))return;passed();
        begin(Stage::InputMaterial);out.beforeMaterial=material(original,x);passed();
        begin(Stage::Split);auto domain=binding.domain();Eigen::VectorXd beta=out.oldBeta;
        out.split=domain.splitFractionalLeaves(beta,settings.fractionalEpsilon);passed();
        begin(Stage::Merge);out.merge=domain.mergeBinarySiblings(beta,settings.fractionalEpsilon);passed();
        snapshot(domain,out.newLeaves);out.newBeta=beta;out.newPartitionSha256=domain.partitionHash();
        begin(Stage::Rebuild);
        const SourceCageBinding rebound(domain,binding.bbw(),binding.config());
        const SourceConformingBoundary boundary(rebound,original.source().boundary().config());
        const SourceCageMomentFunctional moments(boundary,original.source().density(),original.source().config());
        const SourceCagePaperFunctional next(moments,original.config(),original.frame());
        out.newBindingSha256=rebound.outputSha256();out.newFunctionalSha256=functionalHash(next);
        account(out,add(evidence,evaluationPlan(next)));
        Eigen::VectorXd state(next.stateCount());state.head(next.betaCount())=beta;
        if(next.hasAngle())state[next.angleOffset()]=x[original.angleOffset()];
        state.segment(next.handleOffset(),4*next.handleCount())=x.segment(original.handleOffset(),4*original.handleCount());
        out.newStateSha256=stateHash(state);passed();
        begin(Stage::OutputGeometry);
        out.outputGeometry=std::make_shared<const SourceCageGeometryAudit>(SourceCageGeometryChecker(next,settings.geometry).check(state));
        if(!geometryPass(out,*out.outputGeometry))return;passed();
        begin(Stage::OutputMaterial);out.afterMaterial=material(next,state);
        out.deformedMomentDifference=out.afterMaterial.worldMoments.values()-out.beforeMaterial.worldMoments.values();
        if(!out.deformedMomentDifference.allFinite())throw std::overflow_error("nonfinite deformed moment difference");
        out.deformedMomentDifferenceAvailable=true;passed();
        // The only checkpoint constructor call; nothing escapes before finish.
        data->checkpoint=std::shared_ptr<const SourceCageAdaptiveCheckpoint>(new SourceCageAdaptiveCheckpoint(next,std::move(state)));
        out.outcome=Outcome::Published;
    };
    try {run();}
    catch(const std::bad_alloc&) {fail(out,Outcome::ResourceLimit,"RESOURCE_LIMIT: allocation failure");}
    catch(const std::length_error& e) {fail(out,Outcome::ResourceLimit,e.what());}
    catch(const std::exception& e) {fail(out,Outcome::Unknown,e.what());}
    catch(...) {fail(out,Outcome::Unknown,"unknown actual adaptive rebind exception");}
    if(out.outcome!=Outcome::Published)data->checkpoint.reset();
    finish(out); // Diagnostic allocation failure may propagate; no audit escapes.
    return SourceCageAdaptiveRebindAudit(std::move(data));
}
}
