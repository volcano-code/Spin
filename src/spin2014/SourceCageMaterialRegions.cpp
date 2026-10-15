#include "spin2014/SourceCageMaterialRegions.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>

namespace spin2014 {
namespace {
using Outcome=SourceCageMaterialRegionsOutcome;
using namespace verification;
[[noreturn]] void resource(const char* s){throw std::length_error(s);}
std::size_t add(std::size_t a,std::size_t b){if(a>SIZE_MAX-b)resource("region count addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b){if(b&&a>SIZE_MAX/b)resource("region count multiplication overflow");return a*b;}
class Hash {
    Sha256 h_;
public:
    void put(std::uint64_t n){unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);h_.update(b,8);}
    void put(int n){put(std::uint64_t(n));}
    void put(bool n){put(std::uint64_t(n));}
    void put(double d){std::uint64_t n;std::memcpy(&n,&d,8);put(n);}
    void put(const std::string& s){put(std::uint64_t(s.size()));h_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    template<class T> void put(const std::vector<T>& v){put(std::uint64_t(v.size()));for(const auto& x:v)put(T(x));}
    template<class T,std::size_t N> void put(const std::array<T,N>& v){put(std::uint64_t(N));for(const auto& x:v)put(x);}
    template<class T> void put(const std::optional<T>& v){put(bool(v));if(v)put(*v);}
    template<class T> void matrix(const Eigen::MatrixBase<T>& m){put(std::uint64_t(m.rows()));put(std::uint64_t(m.cols()));for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)put(double(m(i,j)));}
    void values(const Eigen::VectorXd& x){put(std::uint64_t(x.size()));for(Eigen::Index i=0;i<x.size();++i)put(x[i]);}
    template<class... T> void fields(const T&... x){(put(x),...);}
    std::string finish(){return h_.finish();}
};
void hash(Hash& h,const ExactBoundaryFacts& f){
    h.fields(f.algorithm,f.backend,f.inputDigest,f.orientation,f.intersectionError,f.boundaryValid,f.closed,f.edgeManifold,f.vertexManifold,
        f.consistentlyOriented,f.componentOrientationsValid,f.intersectionCompleted,f.intersectionFree,f.exactPairTests,
        f.unreferencedVertices,f.degenerateFaces,f.nonManifoldVertices,f.boundaryEdges,f.nonManifoldEdges,f.inconsistentEdges,
        f.intersectingPairs,f.faceComponents,f.componentVolumeSigns,f.componentSignedVolumesExact);
}
void hash(Hash& h,const SourceCageRegionEvidence& e){
    const auto& f=e.region;
    h.fields(f.algorithm,f.backend,f.semantics,f.inputDigest,f.regionDigest,f.error,f.completed,f.valid,f.explicitEmpty,
        f.parents,f.depths,f.originalSigns,f.checkerReversedComponents,f.materialComponents,f.boundedVoidComponents,bool(e.boundary));
    if(e.boundary)hash(h,*e.boundary);
}
void hash(Hash& h,const ExactMaterialPartitionFacts& f){
    h.fields(f.algorithm,f.backend,f.semantics,f.emptyProtocol,f.wholeDigest,f.partitionDigest,f.error,f.partDigests,
        f.completed,f.wholeExplicitEmpty,f.partExplicitEmpty,f.unionEqualsWhole,f.noPositiveVolumeOverlap,
        f.missingWholeContent,f.excessUnionContent,f.completePartition,f.partSubsetWhole,f.positiveVolumeOverlapPairs);
}
void hash(Hash& h,const ExactMaterialRelationFacts& f){
    h.fields(f.algorithm,f.semantics,f.aDigest,f.bDigest,f.relationDigest,f.error,f.completed,f.aExplicitEmpty,f.bExplicitEmpty,
        f.intersectionNonempty,f.positiveVolumeOverlap,f.aSubsetB,f.bSubsetA,f.equal,f.symmetricDifferenceVolumeExact);
}
void hash(Hash& h,const ExactBoundaryClearanceFacts& f){
    h.fields(f.algorithm,f.backend,f.semantics,f.aDigest,f.bDigest,f.aInputDigest,f.bInputDigest,f.minimumClearanceExact,
        f.clearanceDigest,f.error,f.state,f.completed,f.aExplicitEmpty,f.bExplicitEmpty,f.applicable,f.atLeastMinimum,
        f.contact,f.squaredDistanceExact,f.witnessFaces);
}
std::string configHash(const SourceCageMaterialRegionsSettings& s){
    Hash h;h.fields(std::string("source-cage-material-regions-config-v1"),s.minimumOuterWall,s.maximumPrimitives,s.maximumPoints,
        s.maximumTriangles,s.maximumOwnedBytes,s.maximumExactDiagnosticBytes,s.maximumNefComplexity);return h.finish();
}
std::string stateHash(const Eigen::VectorXd& x){Hash h;h.put(std::string("source-adaptive-full-state-v1"));h.values(x);return h.finish();}
std::string functionalHash(const SourceCagePaperFunctional& f){
    Hash h;h.put(std::string("source-adaptive-functional-v1"));
    h.fields(f.configSha256(),f.upstreamSha256(),f.source().configSha256(),f.source().upstreamSha256(),
        f.source().boundary().configSha256(),f.source().boundary().outputSha256(),f.source().boundary().binding().configSha256(),
        f.source().boundary().binding().outputSha256(),f.source().boundary().binding().domain().partitionHash());return h.finish();
}
void meshHash(Hash& h,const SourceCageBinaryMesh& m){h.matrix(m.vertices);h.matrix(m.triangles);}
std::string binaryOutputHash(const SourceCageBinaryMaterial& m){
    Hash h;h.put(std::string("actual-deformed-binary-emitted-triangles-v1"));
    meshHash(h,m.outerBoundary());meshHash(h,m.carvedBoundary());meshHash(h,m.materialBoundary());return h.finish();
}
void finish(SourceCageMaterialRegionsFacts& f,bool product){
    Hash h;h.fields(f.semantics,f.designMembershipProof,f.fullManufacturingThickness,int(f.outcome),std::string(f.reason.data()),
        f.configSha256,f.inputSha256,f.sourceSha256,f.functionalSha256,f.partitionSha256,f.bindingSha256,f.binaryStateSha256,
        f.binaryOutputSha256,f.evaluationSha256,f.admitted,f.identitiesChecked,f.representedDesignMembershipEstablished,
        f.selectedCount,f.plannedPoints,f.plannedTriangles,f.plannedOwnedBytes,f.plannedExactDiagnosticBytes,product);
    hash(h,f.outer);hash(h,f.carved);hash(h,f.material);hash(h,f.materialPartition);hash(h,f.selectedPartition);
    h.put(std::uint64_t(f.selected.size()));
    for(const auto& p:f.selected){
        h.fields(p.primitiveId,p.betaOwner,p.sourceKey.level,p.binaryBeta,p.meshSha256);
        for(auto n:p.sourceKey.coordinates)h.put(std::uint64_t(n));
        hash(h,p.region);hash(h,p.subsetOuter);
    }
    hash(h,f.outerClearance);h.put(hashExactNefBudgetFacts(f.nefBudget));f.scientificSha256=h.finish();
}
void terminal(SourceCageMaterialRegionsFacts& f,Outcome o,const char* s){
    f.outcome=o;const auto n=std::min(f.reason.size()-1,std::strlen(s));std::memcpy(f.reason.data(),s,n);f.reason[n]=0;
}
SourceCageRegionEvidence evidence(const ExactMaterialRegion& r){
    SourceCageRegionEvidence e;e.region=r.facts();if(!e.region.explicitEmpty)e.boundary=r.boundaryFacts();return e;
}
void plan(const SourceCageBinaryMaterial& m,SourceCageMaterialRegionsFacts& out){
    const auto& s=out.settings;
    for(auto n:{s.maximumPrimitives,s.maximumPoints,s.maximumTriangles,s.maximumNefComplexity})
        if(!n||n>std::size_t(INT_MAX/3))resource("invalid region native count cap");
    for(auto n:{s.maximumOwnedBytes,s.maximumExactDiagnosticBytes})
        if(!n||n>std::size_t(PTRDIFF_MAX))resource("invalid region byte cap");
    mul(s.maximumTriangles,s.maximumTriangles);
    const auto& f=m.functional();const auto& b=f.source().boundary();
    const auto& primitives=b.binding().primitives();
    if(primitives.size()!=b.primitiveBoundaries().size())throw std::logic_error("immutable primitive boundary identity mismatch");
    std::size_t p=0,t=0,selected=0;
    for(const auto* mesh:{&m.outerBoundary(),&m.carvedBoundary(),&m.materialBoundary()}){
        p=add(p,std::size_t(mesh->vertices.rows()));t=add(t,std::size_t(mesh->triangles.rows()));
    }
    for(std::size_t i=0;i<primitives.size();++i){
        const auto owner=primitives[i].betaIndex;
        if(owner>=f.betaCount())throw std::logic_error("immutable primitive has invalid beta owner");
        const auto beta=m.binaryState()[owner];
        if(beta!=0.&&beta!=1.)throw std::logic_error("immutable upstream state is not binary");
        if(beta==1.){
            selected=add(selected,1);
            const auto faces=b.primitiveBoundaries()[i].size();
            t=add(t,faces);p=add(p,mul(faces,3)); // Conservative unique point upper bound.
        }
    }
    out.selectedCount=selected;out.plannedPoints=p;out.plannedTriangles=t;
    if(selected>s.maximumPrimitives||p>s.maximumPoints||t>s.maximumTriangles||
        std::size_t(b.points().rows())>s.maximumPoints)resource("selected witness/aggregate region mesh cap");
    // Retained evidence, native-region mesh copies, extraction and map scratch,
    // full value-only evaluation and state/BBW/shape scratch. Native exact
    // rational allocator internals are explicitly outside this accounting.
    auto bytes=add(32768,mul(p,6*sizeof(double)+4*sizeof(std::size_t)));
    bytes=add(bytes,mul(t,6*sizeof(int)+4*sizeof(std::size_t)));
    bytes=add(bytes,mul(add(selected,4),4096));
    bytes=add(bytes,mul(std::size_t(b.points().rows()),3*sizeof(double)+sizeof(std::size_t)));
    bytes=add(bytes,mul(add(mul(f.stateCount(),6),add(mul(f.handleCount(),2),mul(f.constraintCount(),4))),sizeof(double)));
    if(f.shapeFunctional())bytes=add(bytes,mul(b.outerPointCount(),3*sizeof(double)));
    out.plannedOwnedBytes=bytes;
    const auto pairs=t?mul(t/2,t-1+t%2):0;
    // Conservative retained pair/edge/component diagnostics (four capacities),
    // nesting matrices and partition witness arrays, before any exact call.
    out.plannedExactDiagnosticBytes=mul(4,add(add(mul(p,16),mul(t,128+8192+sizeof(std::string))),
        add(mul(pairs,sizeof(std::array<int,2>)+2),mul(mul(add(selected,4),add(selected,4)),128))));
    if(bytes>s.maximumOwnedBytes||out.plannedExactDiagnosticBytes>s.maximumExactDiagnosticBytes)
        resource("region owned buffer/exact diagnostic cap before allocation");
}
bool regionGate(const SourceCageRegionEvidence& e,SourceCageMaterialRegionsFacts& out){
    if(e.region.completed&&e.region.valid)return true;
    terminal(out,e.region.completed?Outcome::Rejected:Outcome::Unknown,e.region.error.c_str());return false;
}
bool partitionGate(const ExactMaterialPartitionFacts& p,SourceCageMaterialRegionsFacts& out){
    if(!p.completed){terminal(out,Outcome::Unknown,p.error.c_str());return false;}
    const bool subsets=p.partSubsetWhole&&std::all_of(p.partSubsetWhole->begin(),p.partSubsetWhole->end(),[](bool b){return b;});
    if(p.completePartition==true&&p.unionEqualsWhole==true&&p.noPositiveVolumeOverlap==true&&
        p.missingWholeContent==false&&p.excessUnionContent==false&&subsets)return true;
    terminal(out,Outcome::Rejected,"exact partition rejected: union, overlap, missing/excess or subset gate");return false;
}
SourceCageBinaryMesh primitiveMesh(const SourceConformingBoundary& b,const Eigen::MatrixXd& points,std::size_t id){
    const auto& refs=b.primitiveBoundaries().at(id);
    std::vector<std::size_t> map(std::size_t(points.rows()),SIZE_MAX);
    for(const auto& ref:refs)for(auto index:b.triangles().at(ref.triangle))map.at(index)=0;
    std::size_t count=0;for(auto& index:map)if(index!=SIZE_MAX)index=count++;
    SourceCageBinaryMesh m;m.vertices.resize(count,3);m.triangles.resize(refs.size(),3);
    for(std::size_t i=0;i<map.size();++i)if(map[i]!=SIZE_MAX)m.vertices.row(map[i])=points.row(i);
    for(std::size_t i=0;i<refs.size();++i){
        auto face=b.triangles().at(refs[i].triangle);
        if(refs[i].orientation==-1)std::swap(face[1],face[2]);
        else if(refs[i].orientation!=1)throw std::logic_error("immutable primitive orientation invalid");
        for(int j=0;j<3;++j)m.triangles(i,j)=int(map.at(face[j]));
    }
    return m;
}
}
struct SourceCageMaterialRegionsAudit::Data{
    SourceCageMaterialRegionsFacts facts;
    std::shared_ptr<const SourceCageMaterialRegions> validated;
};
const SourceCageMaterialRegionsFacts& SourceCageMaterialRegionsAudit::facts()const{
    if(!data_)throw std::logic_error("moved-from region audit has no evidence");return data_->facts;
}
std::shared_ptr<const SourceCageMaterialRegions> SourceCageMaterialRegionsAudit::validatedRegions()const{(void)facts();return data_->validated;}
std::shared_ptr<const SourceCageBinaryMaterial> SourceCageMaterialRegions::binaryMaterial()const{
    if(!material_)throw std::logic_error("moved-from region product has no material");return material_;
}
const SourceCageBinaryMaterialAudit& SourceCageMaterialRegions::binaryAudit()const{
    (void)binaryMaterial();return audit_;
}
const char* sourceCageMaterialRegionsOutcomeToken(Outcome o){
    switch(o){case Outcome::Validated:return "EXACT_MATERIAL_REGIONS_VALIDATED_PRE_SEAL";case Outcome::NotAdmitted:return "NOT_ADMITTED";
    case Outcome::Rejected:return "EXACT_MATERIAL_REGIONS_REJECTED";case Outcome::ResourceLimit:return "RESOURCE_LIMIT";default:return "UNKNOWN";}
}
SourceCageMaterialRegionsAudit validateSourceCageMaterialRegions(const SourceCageBinaryMaterialAudit& audit,
    const SourceCageMaterialRegionsSettings& settings){
    if(!std::isfinite(settings.minimumOuterWall)||settings.minimumOuterWall<0.)
        throw std::invalid_argument("minimum outer wall must be finite nonnegative");
    auto data=std::make_shared<SourceCageMaterialRegionsAudit::Data>();auto& out=data->facts;
    out.settings=settings;out.configSha256=configHash(settings);out.inputSha256=audit.facts().scientificSha256;
    ExactNefBudgetPtr budget;
    auto run=[&]{
        const auto m=audit.validatedMaterial();
        if(audit.facts().outcome!=SourceCageBinaryMaterialOutcome::Validated||!m){
            terminal(out,Outcome::NotAdmitted,"requires genuine immutable validated binary material");return;
        }
        const auto& f=m->functional();const auto& b=f.source().boundary();const auto& upstream=audit.facts();
        out.sourceSha256=b.binding().source().sourceSha256();out.functionalSha256=functionalHash(f);
        out.partitionSha256=b.binding().domain().partitionHash();out.bindingSha256=b.binding().outputSha256();
        out.binaryStateSha256=stateHash(m->binaryState());out.binaryOutputSha256=binaryOutputHash(*m);
        if(out.functionalSha256!=upstream.functionalSha256||out.partitionSha256!=upstream.partitionSha256||
            out.binaryStateSha256!=upstream.binaryStateSha256||out.binaryOutputSha256!=upstream.outputSha256){
            terminal(out,Outcome::NotAdmitted,"immutable binary source/functional/partition/state/output identity mismatch");return;
        }
        out.admitted=true;plan(*m,out);
        const auto values=f.evaluateValues(m->binaryState());out.evaluationSha256=values.outputSha256;
        if(values.hasJacobians||values.world.hasJacobians||out.evaluationSha256!=upstream.binaryFunctionalOutputSha256||
            out.selectedCount!=upstream.selectedPrimitives){
            terminal(out,Outcome::NotAdmitted,"binary value-only evaluation or selected-count identity mismatch");return;
        }
        out.identitiesChecked=true;
        budget=std::make_shared<ExactNefBudget>(settings.maximumNefComplexity);
        auto outer=makeExactMaterialRegionWithBudget(m->outerBoundary().vertices,m->outerBoundary().triangles,budget);
        out.outer=evidence(outer);if(!regionGate(out.outer,out))return;
        auto carved=upstream.emptyCavity?emptyExactMaterialRegionWithBudget(budget):
            makeExactMaterialRegionWithBudget(m->carvedBoundary().vertices,m->carvedBoundary().triangles,budget);
        out.carved=evidence(carved);if(!regionGate(out.carved,out))return;
        auto material=makeExactMaterialRegionWithBudget(m->materialBoundary().vertices,m->materialBoundary().triangles,budget);
        out.material=evidence(material);if(!regionGate(out.material,out))return;
        if(outer.facts().explicitEmpty||material.facts().explicitEmpty||outer.facts().materialComponents.value_or(0)<=0||
            material.facts().materialComponents.value_or(0)<=0){
            terminal(out,Outcome::Rejected,"outer and material must be positive nonempty regions");return;
        }
        out.materialPartition=checkExactMaterialPartitionWithBudget(outer,{material,carved},budget);
        if(!partitionGate(out.materialPartition,out))return;
        std::vector<ExactMaterialRegion> selected;selected.reserve(std::max(std::size_t(1),out.selectedCount));
        out.selected.reserve(out.selectedCount);
        const auto& primitives=b.binding().primitives();
        for(std::size_t id=0;id<primitives.size();++id){
            const auto& p=primitives[id];if(m->binaryState()[p.betaIndex]!=1.)continue;
            SourceCageSelectedRegionEvidence e;e.primitiveId=id;e.betaOwner=p.betaIndex;e.sourceKey=p.sourceKey;e.binaryBeta=1.;
            auto mesh=primitiveMesh(b,values.world.deformedPoints,id);
            Hash h;h.put(std::string("actual-selected-oriented-primitive-mesh-v1"));meshHash(h,mesh);e.meshSha256=h.finish();
            selected.push_back(makeExactMaterialRegionWithBudget(mesh.vertices,mesh.triangles,budget));
            e.region=evidence(selected.back());out.selected.push_back(std::move(e));auto& saved=out.selected.back();
            if(!regionGate(saved.region,out))return;
            saved.subsetOuter=relateExactMaterialRegionsWithBudget(selected.back(),outer,budget);
            if(!saved.subsetOuter.completed){terminal(out,Outcome::Unknown,saved.subsetOuter.error.c_str());return;}
            if(saved.subsetOuter.aSubsetB!=true){terminal(out,Outcome::Rejected,"selected represented design primitive is not subset outer");return;}
        }
        if(selected.empty())selected.push_back(emptyExactMaterialRegionWithBudget(budget));
        if(!selected.back().facts().valid){terminal(out,Outcome::Unknown,"explicit empty witness construction incomplete");return;}
        out.selectedPartition=checkExactMaterialPartitionWithBudget(carved,selected,budget);
        if(!partitionGate(out.selectedPartition,out))return;
        out.representedDesignMembershipEstablished=true;
        out.outerClearance=checkExactBoundaryClearance(outer,carved,settings.minimumOuterWall);
        const auto& c=out.outerClearance;
        if(!c.completed){terminal(out,Outcome::Unknown,c.error.c_str());return;}
        const bool clear=upstream.emptyCavity?
            (c.state=="NOT_APPLICABLE"&&c.applicable==false&&!c.contact&&!c.squaredDistanceExact&&!c.witnessFaces&&!c.atLeastMinimum):
            (c.applicable==true&&c.contact==false&&c.atLeastMinimum==true&&c.squaredDistanceExact&&c.witnessFaces);
        if(!clear){terminal(out,Outcome::Rejected,"actual exact outer-to-carved clearance gate rejected");return;}
        terminal(out,Outcome::Validated,"exact represented regions, partitions, constructive design membership and outer clearance passed; no full thickness or seal");
        data->validated=std::shared_ptr<const SourceCageMaterialRegions>(new SourceCageMaterialRegions(audit,m));
    };
    try{run();}
    catch(const std::length_error& e){terminal(out,Outcome::ResourceLimit,e.what());}
    catch(const std::bad_alloc&){terminal(out,Outcome::ResourceLimit,"region allocation failure");}
    catch(const std::exception& e){terminal(out,Outcome::Unknown,e.what());}
    catch(...){terminal(out,Outcome::Unknown,"unknown exact region validation failure");}
    if(budget){out.nefBudget=budget->facts();if(out.nefBudget.exceeded){
        data->validated.reset();terminal(out,Outcome::ResourceLimit,"shared exact Nef complexity budget exceeded; partial evidence retained");
    }}
    finish(out,bool(data->validated));return SourceCageMaterialRegionsAudit(data);
}
}
