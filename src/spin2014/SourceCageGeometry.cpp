#include "spin2014/SourceCageGeometry.h"
#include "verification/ExactMaterialRegion.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <numeric>
#include <new>
#include <stdexcept>

namespace spin2014 {
namespace {
using namespace verification;
using Stage=SourceCageGeometryStage;
using Status=SourceCageGeometryStageOutcome;
using Outcome=SourceCageGeometryOutcome;
[[noreturn]] void resource(const char* s) {throw std::length_error(std::string("RESOURCE_LIMIT: source geometry ")+s);}
std::size_t add(std::size_t a,std::size_t b) {if(a>SIZE_MAX-b)resource("addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b) {if(b&&a>SIZE_MAX/b)resource("multiplication overflow");return a*b;}
std::size_t pairs(std::size_t n) {return n?mul(n/2,n-1+(n%2)):0;}
class Hash {
    Sha256 h;
public:
    void integer(std::uint64_t n) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);h.update(b,8);}
    void string(const std::string& s) {integer(s.size());h.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double v) {std::uint64_t n;std::memcpy(&n,&v,8);integer(n);}
    template<class T> void vector(const std::vector<T>& v) {integer(v.size());for(const auto& n:v)integer(n);}
    template<std::size_t N> void pairVector(const std::vector<std::array<int,N>>& v) {integer(v.size());for(const auto& a:v)for(auto n:a)integer(n);}
    void optional(const std::optional<bool>& v) {integer(v?(*v?2:1):0);}
    std::string finish() {return h.finish();}
};
template<std::size_t N> void copy(std::array<char,N>& a,const std::string& s) {
    const auto n=std::min(N-1,s.size());std::memcpy(a.data(),s.data(),n);a[n]=0;
}
void reason(SourceCageGeometryRecord& r,const std::string& s) {copy(r.reason,s);Hash h;h.string(s);copy(r.reasonSha256,h.finish());}
std::string boundaryHash(const ExactBoundaryFacts& f) {
    Hash h;h.string(f.algorithm);h.string(f.backend);h.string(f.inputDigest);h.string(f.orientation);h.string(f.intersectionError);
    for(bool b:{f.boundaryValid,f.closed,f.edgeManifold,f.vertexManifold,f.consistentlyOriented,f.componentOrientationsValid,f.intersectionCompleted,f.intersectionFree})h.integer(b);
    h.integer(f.exactPairTests);h.vector(f.unreferencedVertices);h.vector(f.degenerateFaces);h.vector(f.nonManifoldVertices);
    h.pairVector(f.boundaryEdges);h.pairVector(f.nonManifoldEdges);h.pairVector(f.inconsistentEdges);h.pairVector(f.intersectingPairs);
    h.vector(f.faceComponents);h.vector(f.componentVolumeSigns);h.integer(f.componentSignedVolumesExact.size());
    for(const auto& s:f.componentSignedVolumesExact)h.string(s);return h.finish();
}
std::string regionHash(const ExactMaterialRegion& r) {
    const auto& f=r.facts();Hash h;h.string(boundaryHash(r.boundaryFacts()));
    for(const auto& s:{f.algorithm,f.backend,f.semantics,f.inputDigest,f.regionDigest,f.error})h.string(s);
    h.integer(f.completed);h.integer(f.valid);h.integer(f.explicitEmpty);h.vector(f.parents);h.vector(f.depths);h.vector(f.originalSigns);h.vector(f.checkerReversedComponents);
    h.integer(f.materialComponents.has_value());if(f.materialComponents)h.integer(*f.materialComponents);
    h.integer(f.boundedVoidComponents.has_value());if(f.boundedVoidComponents)h.integer(*f.boundedVoidComponents);return h.finish();
}
std::string relationHash(const ExactMaterialRelationFacts& f) {
    Hash h;for(const auto& s:{f.algorithm,f.semantics,f.aDigest,f.bDigest,f.relationDigest,f.error})h.string(s);
    h.integer(f.completed);h.integer(f.aExplicitEmpty);h.integer(f.bExplicitEmpty);
    h.optional(f.intersectionNonempty);h.optional(f.positiveVolumeOverlap);h.optional(f.aSubsetB);h.optional(f.bSubsetA);h.optional(f.equal);return h.finish();
}
std::string clearanceHash(const ExactBoundaryClearanceFacts& f) {
    Hash h;for(const auto& s:{f.algorithm,f.backend,f.semantics,f.aDigest,f.bDigest,f.aInputDigest,f.bInputDigest,f.minimumClearanceExact,f.clearanceDigest,f.error,f.state})h.string(s);
    h.integer(f.completed);h.integer(f.aExplicitEmpty);h.integer(f.bExplicitEmpty);h.optional(f.applicable);h.optional(f.atLeastMinimum);h.optional(f.contact);
    h.integer(f.squaredDistanceExact.has_value());if(f.squaredDistanceExact)h.string(*f.squaredDistanceExact);
    h.integer(f.witnessFaces.has_value());if(f.witnessFaces)for(auto n:*f.witnessFaces)h.integer(n);return h.finish();
}
bool resourceError(const std::string& s) {return s.find("RESOURCE_LIMIT")!=s.npos||s.find("bad_alloc")!=s.npos||s.find("bad allocation")!=s.npos;}
void validate(const SourceCageGeometrySettings& s) {
    if(!std::isfinite(s.minimumOuterDomainClearance)||s.minimumOuterDomainClearance<0)throw std::invalid_argument("invalid minimum outer/domain clearance");
    for(auto n:{s.maximumPoints,s.maximumTriangles,s.maximumPrimitives,s.maximumPairs,s.maximumRecords,s.maximumConnectedShells,s.maximumExactStringBytes})
        if(!n||n>std::size_t(INT_MAX/3))resource("zero or unrepresentable cap");
    if(!s.maximumCheckerBytes||s.maximumCheckerBytes>std::size_t(PTRDIFF_MAX)||!s.maximumNefComplexity||s.maximumNefComplexity>std::size_t(INT_MAX))resource("zero or unrepresentable byte/Nef cap");
    // Validate configured products too, before any size-dependent allocation.
    mul(s.maximumRecords,sizeof(SourceCageGeometryRecord));mul(s.maximumTriangles,s.maximumTriangles);mul(s.maximumPrimitives,s.maximumExactStringBytes);
}
std::string inputHash(const SourceCagePaperFunctional& f,const SourceCageGeometrySettings& s,const Eigen::VectorXd& x) {
    Hash h;h.string("source-cage-geometry-v1 accounting-v1");h.string(f.configSha256());h.string(f.upstreamSha256());
    h.string(f.source().configSha256());h.string(f.source().upstreamSha256());h.string(f.source().boundary().outputSha256());
    h.integer(x.size());for(Eigen::Index i=0;i<x.size();++i)h.real(x[i]);h.real(s.minimumOuterDomainClearance);
    for(auto n:{s.maximumPoints,s.maximumTriangles,s.maximumPrimitives,s.maximumPairs,s.maximumRecords,s.maximumCheckerBytes,s.maximumNefComplexity,s.maximumConnectedShells,s.maximumExactStringBytes})h.integer(n);
    return h.finish();
}
// Known maximum diagnostic vector payload for one indexed mesh, including every
// possible intersecting face pair. Exact rational/native predicate scratch is
// separately opaque; per-component rational volume strings are bounded below.
std::size_t diagnostics(std::size_t p,std::size_t f,const SourceCageGeometrySettings& s) {
    auto bytes=add(mul(p,16),add(mul(f,128),mul(pairs(f),sizeof(std::array<int,2>))));
    // Binary64 signed tetrahedral volumes have bounded-degree dyadic arithmetic;
    // reserve 8192 bytes per component even when the retained-string cap is tiny.
    bytes=add(bytes,mul(std::min(f,s.maximumConnectedShells),add(sizeof(std::string),std::max(std::size_t(8192),s.maximumExactStringBytes))));
    // Geometric-growth vector capacity and the audit->region facts copy may
    // coexist. Four payloads conservatively cover both, excluding allocator nodes.
    return add(mul(bytes,4),mul(f,f)); // signed hierarchy: conservative one byte per slot
}
void plan(const SourceCagePaperFunctional& f,const SourceCageGeometrySettings& s,SourceCageGeometryFacts& out) {
    const auto& b=f.source().boundary();const auto& outer=b.binding().source().outerFaces();
    out.pointCount=b.points().rows();out.primitiveCount=b.primitiveBoundaries().size();out.pairCount=pairs(out.primitiveCount);
    out.triangleCount=outer.rows();std::size_t bytes=mul(out.pointCount,3*sizeof(double));
    // evaluateValues output, raw full-state decomposition, constraints and full
    // one-row BBW weight scratch; upstream immutable source is shared, not copied.
    bytes=add(bytes,mul(add(f.stateCount(),add(f.betaCount(),add(f.constraintCount(),mul(f.handleCount(),2)))),sizeof(double)));
    if(f.shapeFunctional())bytes=add(bytes,mul(b.outerPointCount(),3*sizeof(double))); // value-only shape residual
    auto account=[&](std::size_t p,std::size_t t) {
        // Local mesh, region's retained mesh copy, local index mapping and DSU.
        bytes=add(bytes,add(mul(p,2*3*sizeof(double)+3*sizeof(std::size_t)),mul(t,2*3*sizeof(int))));
        bytes=add(bytes,diagnostics(p,t,s));
    };
    account(b.outerPointCount(),outer.rows());
    for(const auto& refs:b.primitiveBoundaries()) {out.triangleCount=add(out.triangleCount,refs.size());account(mul(refs.size(),3),refs.size());}
    out.plannedRecords=add(3,add(mul(out.primitiveCount,3),out.pairCount));
    bytes=add(bytes,mul(out.plannedRecords,sizeof(SourceCageGeometryRecord)));
    bytes=add(bytes,mul(out.primitiveCount,add(s.maximumExactStringBytes,1024))); // retained distance and region/vector metadata
    out.plannedCheckerBytes=bytes;
    if(out.pointCount>s.maximumPoints||out.pointCount>std::size_t(INT_MAX/3)||out.triangleCount>s.maximumTriangles||out.triangleCount>std::size_t(INT_MAX/3)||
       !out.primitiveCount||out.primitiveCount>s.maximumPrimitives||out.pairCount>s.maximumPairs||out.plannedRecords>s.maximumRecords||bytes>s.maximumCheckerBytes)resource("preflight cap exceeded");
}
struct Mesh {Eigen::MatrixXd v;Eigen::MatrixXi f;Eigen::Vector3d lo,hi;};
Mesh primitiveMesh(const SourceConformingBoundary& b,const Eigen::MatrixXd& points,std::size_t id) {
    const auto& refs=b.primitiveBoundaries()[id];std::vector<std::size_t> ids;ids.reserve(mul(refs.size(),3));
    for(const auto& ref:refs)for(auto i:b.triangles().at(ref.triangle))ids.push_back(i);
    std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    Mesh m;m.v.resize(ids.size(),3);m.f.resize(refs.size(),3);
    for(std::size_t i=0;i<ids.size();++i)m.v.row(i)=points.row(ids[i]);
    for(std::size_t i=0;i<refs.size();++i) {
        auto t=b.triangles().at(refs[i].triangle);if(refs[i].orientation==-1)std::swap(t[1],t[2]);
        else if(refs[i].orientation!=1)throw std::invalid_argument("invalid immutable primitive orientation");
        for(int j=0;j<3;++j)m.f(i,j)=int(std::lower_bound(ids.begin(),ids.end(),t[j])-ids.begin());
    }
    m.lo=m.v.colwise().minCoeff();m.hi=m.v.colwise().maxCoeff();return m;
}
void shellPreflight(const Mesh& m,const SourceCageGeometrySettings& s) {
    std::vector<std::size_t> roots(m.v.rows());std::iota(roots.begin(),roots.end(),0);
    auto root=[&](std::size_t i) {while(roots[i]!=i){roots[i]=roots[roots[i]];i=roots[i];}return i;};
    for(Eigen::Index i=0;i<m.f.rows();++i)for(int j=1;j<3;++j)roots[root(m.f(i,j))]=root(m.f(i,0));
    std::size_t n=0;for(std::size_t i=0;i<roots.size();++i)n+=root(i)==i;
    if(n>s.maximumConnectedShells)resource("connected shell cap exceeded");
}
void boundStrings(const ExactMaterialRegion& r,const SourceCageGeometrySettings& s) {
    for(const auto& v:r.boundaryFacts().componentSignedVolumesExact)if(v.size()>s.maximumExactStringBytes)resource("exact component volume string cap exceeded");
}
void finish(SourceCageGeometryFacts& f) {
    Hash h;h.string(f.semantics);h.string(f.inputSha256);h.string(f.sourceOutputSha256);h.integer(unsigned(f.outcome));
    for(auto n:{f.pointCount,f.triangleCount,f.primitiveCount,f.pairCount,f.plannedRecords,f.plannedCheckerBytes,f.pairsVisited,f.exactPairRelations})h.integer(n);
    h.string(hashExactNefBudgetFacts(f.nefBudget));h.integer(f.records.size());
    for(const auto& r:f.records) {
        h.integer(unsigned(r.stage));h.integer(unsigned(r.outcome));h.integer(r.first);h.integer(r.second);h.integer(r.betaOwner);h.integer(r.secondBetaOwner);
        h.string(r.evidenceSha256.data());h.string(r.reasonSha256.data());h.string(r.reason.data());h.string(r.squaredDistanceExact);
        h.integer(r.witnessFaces.has_value());if(r.witnessFaces)for(auto n:*r.witnessFaces)h.integer(n);
    }f.scientificSha256=h.finish();
}
}
const char* sourceCageGeometryOutcomeToken(SourceCageGeometryOutcome o) {
    switch(o){case Outcome::ChecksPassed:return "CHECKS_PASSED";case Outcome::RejectedGeometry:return "REJECTED_GEOMETRY";case Outcome::ResourceLimit:return "RESOURCE_LIMIT";default:return "UNKNOWN";}
}
SourceCageGeometryChecker::SourceCageGeometryChecker(const SourceCagePaperFunctional& f,const SourceCageGeometrySettings& s):source_(f),settings_(s) {validate(s);}
SourceCageGeometryAudit SourceCageGeometryChecker::check(const Eigen::VectorXd& x) const {
    if(std::size_t(x.size())!=source_.stateCount())throw std::invalid_argument("source geometry full-state dimension mismatch");
    SourceCageGeometryAudit audit;auto& out=audit.facts_;out.inputSha256=inputHash(source_,settings_,x);
    auto budget=std::make_shared<ExactNefBudget>(settings_.maximumNefComplexity);
    // A single preflight event exists even when the requested tape cannot fit.
    out.records.emplace_back();
    auto event=[&](Stage stage,std::size_t a=SourceCageGeometryRecord::absent,std::size_t b=SourceCageGeometryRecord::absent)->SourceCageGeometryRecord& {
        out.records.emplace_back();auto& r=out.records.back();r.stage=stage;r.first=a;r.second=b;
        const auto& p=source_.source().boundary().binding().primitives();
        if(a!=SourceCageGeometryRecord::absent)r.betaOwner=p.at(a).betaIndex;
        if(b!=SourceCageGeometryRecord::absent)r.secondBetaOwner=p.at(b).betaIndex;return r;
    };
    auto failure=[&](const std::string& error,bool completed) {
        auto& r=out.records.back();reason(r,error);
        if(budget->facts().exceeded||resourceError(error)) {r.outcome=Status::ResourceLimit;out.outcome=Outcome::ResourceLimit;}
        else if(completed) {r.outcome=Status::Rejected;out.outcome=Outcome::RejectedGeometry;}
        else {r.outcome=Status::Unknown;out.outcome=Outcome::Unknown;}
    };
    auto run=[&] {
        plan(source_,settings_,out);out.records.reserve(out.plannedRecords);out.records.back().outcome=Status::Passed;
        event(Stage::SourceEvaluation);const auto evaluated=source_.evaluateValues(x);
        out.sourceOutputSha256=evaluated.outputSha256;copy(out.records.back().evidenceSha256,evaluated.outputSha256);out.records.back().outcome=Status::Passed;
        const auto& points=evaluated.world.deformedPoints;const auto& b=source_.source().boundary();
        event(Stage::OuterRegion);Mesh outer;outer.v=points.topRows(b.outerPointCount());outer.f=b.binding().source().outerFaces();shellPreflight(outer,settings_);
        const auto region=makeExactMaterialRegionWithBudget(outer.v,outer.f,budget);copy(out.records.back().evidenceSha256,regionHash(region));boundStrings(region,settings_);
        if(!region.facts().valid||!region.facts().completed) {failure(region.facts().error,region.facts().completed);return;}
        out.records.back().outcome=Status::Passed;
        std::vector<ExactMaterialRegion> cells;cells.reserve(out.primitiveCount);std::vector<std::array<Eigen::Vector3d,2>> boxes;boxes.reserve(out.primitiveCount);
        for(std::size_t i=0;i<out.primitiveCount;++i) {
            event(Stage::PrimitiveRegion,i);const auto mesh=primitiveMesh(b,points,i);shellPreflight(mesh,settings_);
            auto cell=makeExactMaterialRegionWithBudget(mesh.v,mesh.f,budget);copy(out.records.back().evidenceSha256,regionHash(cell));boundStrings(cell,settings_);
            if(!cell.facts().valid||!cell.facts().completed) {failure(cell.facts().error,cell.facts().completed);return;}
            if(cell.boundaryFacts().componentVolumeSigns!=std::vector<int>{1}) {failure("primitive must have one positive nonzero closed component",true);return;}
            out.records.back().outcome=Status::Passed;
            event(Stage::PrimitiveSubset,i);const auto rel=relateExactMaterialRegionsWithBudget(cell,region,budget);copy(out.records.back().evidenceSha256,relationHash(rel));
            if(!rel.completed||!rel.aSubsetB.has_value()) {failure(rel.error,false);return;}
            if(!*rel.aSubsetB) {failure("primitive not subset of actual outer material",true);return;}out.records.back().outcome=Status::Passed;
            auto& cr=event(Stage::OuterDomainClearance,i);const auto c=checkExactBoundaryClearance(region,cell,settings_.minimumOuterDomainClearance);copy(cr.evidenceSha256,clearanceHash(c));
            if(!c.completed||!c.applicable.value_or(false)||!c.contact.has_value()||!c.atLeastMinimum.has_value()||!c.squaredDistanceExact.has_value()) {failure(c.error,false);return;}
            if(c.squaredDistanceExact->size()>settings_.maximumExactStringBytes)resource("exact clearance string cap exceeded");
            cr.squaredDistanceExact=*c.squaredDistanceExact;cr.witnessFaces=c.witnessFaces;
            if(*c.contact||!*c.atLeastMinimum) {failure("outer/domain contact or insufficient conservative clearance",true);return;}cr.outcome=Status::Passed;
            cells.push_back(std::move(cell));boxes.push_back({mesh.lo,mesh.hi});
        }
        for(std::size_t i=0;i<cells.size();++i)for(std::size_t j=i+1;j<cells.size();++j) {
            auto& r=event(Stage::PrimitivePair,i,j);++out.pairsVisited;
            bool disjoint=false;for(int axis=0;axis<3;++axis)disjoint|=boxes[i][1][axis]<boxes[j][0][axis]||boxes[j][1][axis]<boxes[i][0][axis];
            if(disjoint) {Hash h;h.string("exact-safe binary64 strict AABB separation");for(auto k:{i,j})for(const auto& v:boxes[k])for(int a=0;a<3;++a)h.real(v[a]);copy(r.evidenceSha256,h.finish());r.outcome=Status::BroadPhaseDisjoint;continue;}
            ++out.exactPairRelations;const auto rel=relateExactMaterialRegionsWithBudget(cells[i],cells[j],budget);copy(r.evidenceSha256,relationHash(rel));
            if(!rel.completed||!rel.positiveVolumeOverlap.has_value()) {failure(rel.error,false);return;}
            if(*rel.positiveVolumeOverlap) {failure("positive-volume primitive pair overlap",true);return;}r.outcome=Status::Passed;
        }
        out.outcome=Outcome::ChecksPassed;
    };
    try {run();}
    catch(const std::bad_alloc&) {failure("RESOURCE_LIMIT: allocation failure",false);}
    catch(const std::length_error& e) {failure(std::string("RESOURCE_LIMIT: ")+e.what(),false);}
    catch(const std::exception& e) {failure(e.what(),false);}
    catch(...) {failure("unknown source/geometry exception",false);}
    out.nefBudget=budget->facts();finish(out);return audit;
}
}
