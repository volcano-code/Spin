#include "spin2014/SourceCageComponentClearance.h"
#include "verification/Sha256.h"
#include <CGAL/Gmpq.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>

namespace spin2014 {
namespace {
using namespace verification;
using Outcome=SourceCageComponentClearanceOutcome;
[[noreturn]] void resource(const char* s){throw std::length_error(s);}
std::size_t add(std::size_t a,std::size_t b){if(a>SIZE_MAX-b)resource("component count addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b){if(b&&a>SIZE_MAX/b)resource("component count multiplication overflow");return a*b;}
std::size_t pairCount(std::size_t n){return n?mul(n/2,n-1+n%2):0;}
class Hash {
    Sha256 h_;
public:
    void raw(const void* p,std::size_t n){h_.update(static_cast<const unsigned char*>(p),n);}
    void put(std::uint64_t n){unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);raw(b,8);}
    void put(int n){put(std::uint64_t(n));}
    void put(bool n){put(std::uint64_t(n));}
    void put(double d){std::uint64_t n;std::memcpy(&n,&d,8);put(n);}
    void put(const std::string& s){put(std::uint64_t(s.size()));raw(s.data(),s.size());}
    template<class T> void put(const std::vector<T>& v){put(std::uint64_t(v.size()));for(const auto& x:v)put(T(x));}
    template<class T,std::size_t N> void put(const std::array<T,N>& v){put(std::uint64_t(N));for(const auto& x:v)put(x);}
    template<class T> void put(const std::optional<T>& v){put(bool(v));if(v)put(*v);}
    template<class... T> void fields(const T&... x){(put(x),...);}
    std::string finish(){return h_.finish();}
};
// The exact backend's public input identity protocol, recomputed from actual
// arrays. No copied digest is accepted as proof of the material bytes.
std::string meshHash(const SourceCageBinaryMesh& m){
    Hash h;constexpr char prefix[]="exact-boundary-input-v1";h.raw(prefix,sizeof(prefix));
    const unsigned char mode=static_cast<unsigned char>(ExactBoundaryOrientation::AnyNonzero);h.raw(&mode,1);
    h.fields(std::uint64_t(m.vertices.rows()),std::uint64_t(m.triangles.rows()));
    for(Eigen::Index i=0;i<m.vertices.rows();++i)for(int j=0;j<3;++j)h.put(m.vertices(i,j));
    for(Eigen::Index i=0;i<m.triangles.rows();++i)for(int j=0;j<3;++j)h.put(std::uint64_t(m.triangles(i,j)));
    return h.finish();
}
void hash(Hash& h,const ExactBoundaryFacts& f){
    h.fields(f.algorithm,f.backend,f.inputDigest,f.orientation,f.intersectionError,f.boundaryValid,f.closed,f.edgeManifold,f.vertexManifold,
        f.consistentlyOriented,f.componentOrientationsValid,f.intersectionCompleted,f.intersectionFree,f.exactPairTests,
        f.unreferencedVertices,f.degenerateFaces,f.nonManifoldVertices,f.boundaryEdges,f.nonManifoldEdges,f.inconsistentEdges,
        f.intersectingPairs,f.faceComponents,f.componentVolumeSigns,f.componentSignedVolumesExact);
}
void hash(Hash& h,const SourceCageRegionEvidence& e){const auto& f=e.region;
    h.fields(f.algorithm,f.backend,f.semantics,f.inputDigest,f.regionDigest,f.error,f.completed,f.valid,f.explicitEmpty,
        f.parents,f.depths,f.originalSigns,f.checkerReversedComponents,f.materialComponents,f.boundedVoidComponents,bool(e.boundary));
    if(e.boundary)hash(h,*e.boundary);
}
void hash(Hash& h,const ExactBoundaryClearanceFacts& f){
    h.fields(f.algorithm,f.backend,f.semantics,f.aDigest,f.bDigest,f.aInputDigest,f.bInputDigest,f.minimumClearanceExact,
        f.clearanceDigest,f.error,f.state,f.completed,f.aExplicitEmpty,f.bExplicitEmpty,f.applicable,f.atLeastMinimum,
        f.contact,f.squaredDistanceExact,f.witnessFaces);
}
std::string configHash(const SourceCageComponentClearanceSettings& s){Hash h;
    h.fields(std::string("source-cage-component-clearance-config-v1"),s.minimumComponentSeparation,s.maximumComponents,s.maximumPoints,
        s.maximumTriangles,s.maximumPairs,s.maximumOwnedBytes,s.maximumExactDiagnosticBytes,s.maximumNefComplexity);return h.finish();
}
void finish(SourceCageComponentClearanceFacts& f,bool product){Hash h;
    h.fields(std::string("source-cage-component-clearance-v1"),int(f.outcome),std::string(f.reason.data()),f.scope,f.semantics,
        f.fullManufacturingThickness,f.state,f.admitted,f.identitiesChecked,f.coverageChecked,f.componentCount,f.plannedPoints,
        f.plannedTriangles,f.plannedPairCount,f.plannedOwnedBytes,f.plannedExactDiagnosticBytes,f.sourceFaceComponents,
        f.applicable,f.atLeastMinimum,f.contact,f.minimumSquaredDistanceExact,f.minimumPairIndex,f.minimumComponentIds,
        f.minimumLocalWitnessFaces,f.minimumParentWitnessFaces,f.inputSha256,f.materialInputSha256,f.configSha256,product);
    h.put(std::uint64_t(f.components.size()));for(const auto& c:f.components){
        h.fields(c.componentId,c.sourceVertexIds,c.sourceFaceIds,c.originalSign,c.parent,c.depth,c.checkerReversed,c.originalInputSha256,c.checkerInputSha256);hash(h,c.region);
    }
    h.put(std::uint64_t(f.pairs.size()));for(const auto& p:f.pairs){h.fields(p.componentIds,p.parentWitnessFaces);hash(h,p.clearance);}
    h.put(hashExactNefBudgetFacts(f.nefBudget));f.scientificSha256=h.finish();
}
void terminal(SourceCageComponentClearanceFacts& f,Outcome o,const char* s){f.outcome=o;
    const auto n=std::min(f.reason.size()-1,std::strlen(s));std::memcpy(f.reason.data(),s,n);f.reason[n]=0;
}
void plan(const SourceCageBinaryMesh& m,const SourceCageMaterialRegionsFacts& upstream,SourceCageComponentClearanceFacts& out){
    const auto& s=out.settings;
    for(auto n:{s.maximumComponents,s.maximumPoints,s.maximumTriangles,s.maximumPairs,s.maximumNefComplexity})
        if(!n||n>std::size_t(INT_MAX/3))resource("invalid component native count cap");
    for(auto n:{s.maximumOwnedBytes,s.maximumExactDiagnosticBytes})if(!n||n>std::size_t(PTRDIFF_MAX))resource("invalid component byte cap");
    if(!upstream.material.boundary||m.vertices.rows()<=0||m.triangles.rows()<=0||m.vertices.cols()!=3||m.triangles.cols()!=3)
        throw std::logic_error("immutable material boundary evidence missing");
    const auto& b=*upstream.material.boundary;const auto& r=upstream.material.region;
    const auto p=std::size_t(m.vertices.rows()),t=std::size_t(m.triangles.rows()),n=b.componentVolumeSigns.size();
    out.plannedPoints=p;out.plannedTriangles=t;out.componentCount=n;
    if(!n||n>s.maximumComponents||p>s.maximumPoints||t>s.maximumTriangles)resource("component/point/triangle cap");
    if(b.faceComponents.size()!=t||b.componentSignedVolumesExact.size()!=n||r.parents.size()!=n||r.depths.size()!=n||r.originalSigns.size()!=n)
        throw std::logic_error("immutable component diagnostics size mismatch");
    out.plannedPairCount=pairCount(n);if(out.plannedPairCount>s.maximumPairs)resource("component pair cap");
    // Mesh copies, global-to-local/coverage scratch, component maps, input
    // faceComponents, native wrapper storage, and bounded pair containers.
    out.plannedOwnedBytes=add(32768,add(mul(p,12*sizeof(double)+8*sizeof(std::size_t)),
        add(mul(t,12*sizeof(int)+8*sizeof(std::size_t)),add(mul(n,8192),mul(out.plannedPairCount,8192)))));
    // Conservatively budget edge/pair diagnostics and rational serialization;
    // exact allocator internals and GMP intermediates use the external guard.
    out.plannedExactDiagnosticBytes=mul(4,add(mul(p,128),add(mul(t,8192+256),
        add(mul(pairCount(t),sizeof(std::array<int,2>)+2),mul(out.plannedPairCount,32768)))));
    if(out.plannedOwnedBytes>s.maximumOwnedBytes||out.plannedExactDiagnosticBytes>s.maximumExactDiagnosticBytes)
        resource("component owned/exact diagnostic budget before allocation");
}
// Bound all retained variable diagnostic strings, including actual rational
// outputs, before copying them out of the backend's externally guarded result.
struct DiagnosticBytes {
    std::size_t used=0,cap;
    void put(const std::string& s){used=add(used,add(s.size(),1));if(used>cap)resource("actual exact diagnostic byte cap");}
    template<class T> void put(const std::vector<T>& v){used=add(used,mul(v.size(),sizeof(T)));if(used>cap)resource("actual exact diagnostic array cap");}
    template<class... T> void fields(const T&... x){(put(x),...);}
    void region(const ExactMaterialRegion& value){const auto& r=value.facts();const auto& b=value.boundaryFacts();
        fields(r.algorithm,r.backend,r.semantics,r.inputDigest,r.regionDigest,r.error,r.parents,r.depths,r.originalSigns,r.checkerReversedComponents);
        fields(b.algorithm,b.backend,b.inputDigest,b.orientation,b.intersectionError,b.unreferencedVertices,b.degenerateFaces,b.nonManifoldVertices,
            b.boundaryEdges,b.nonManifoldEdges,b.inconsistentEdges,b.intersectingPairs,b.faceComponents,b.componentVolumeSigns,b.componentSignedVolumesExact);
        for(const auto& s:b.componentSignedVolumesExact)put(s);
    }
    void pair(const ExactBoundaryClearanceFacts& c){fields(c.algorithm,c.backend,c.semantics,c.aDigest,c.bDigest,c.aInputDigest,c.bInputDigest,c.minimumClearanceExact,c.clearanceDigest,c.error,c.state);if(c.squaredDistanceExact)put(*c.squaredDistanceExact);}
};
}
struct SourceCageComponentClearanceAudit::Data {SourceCageComponentClearanceFacts facts;std::shared_ptr<const SourceCageComponentClearance> validated;};
const SourceCageComponentClearanceFacts& SourceCageComponentClearanceAudit::facts()const{if(!data_)throw std::logic_error("moved-from component audit has no evidence");return data_->facts;}
std::shared_ptr<const SourceCageComponentClearance> SourceCageComponentClearanceAudit::validatedClearance()const{(void)facts();return data_->validated;}
const SourceCageMaterialRegionsAudit& SourceCageComponentClearance::regionsAudit()const{(void)audit_.facts();return audit_;}
std::shared_ptr<const SourceCageMaterialRegions> SourceCageComponentClearance::regions()const{return regionsAudit().validatedRegions();}
const char* sourceCageComponentClearanceOutcomeToken(Outcome o){switch(o){case Outcome::Validated:return "BOUNDARY_COMPONENT_SEPARATION_VALIDATED_PRE_SEAL";case Outcome::NotAdmitted:return "NOT_ADMITTED";case Outcome::Rejected:return "BOUNDARY_COMPONENT_SEPARATION_REJECTED";case Outcome::ResourceLimit:return "RESOURCE_LIMIT";default:return "UNKNOWN";}}
SourceCageComponentClearanceAudit validateSourceCageComponentClearance(const SourceCageMaterialRegionsAudit& audit,const SourceCageComponentClearanceSettings& settings){
    if(!std::isfinite(settings.minimumComponentSeparation)||settings.minimumComponentSeparation<0.)throw std::invalid_argument("minimum component separation must be finite nonnegative");
    auto data=std::make_shared<SourceCageComponentClearanceAudit::Data>();auto& out=data->facts;out.settings=settings;
    out.configSha256=configHash(settings);out.inputSha256=audit.facts().scientificSha256;ExactNefBudgetPtr budget;
    auto run=[&]{
        const auto regions=audit.validatedRegions();if(audit.facts().outcome!=SourceCageMaterialRegionsOutcome::Validated||!regions){terminal(out,Outcome::NotAdmitted,"requires genuine immutable validated regions");return;}
        out.admitted=true;const auto material=regions->binaryMaterial();const auto& mesh=material->materialBoundary();const auto& upstream=audit.facts();
        plan(mesh,upstream,out);const auto& boundary=*upstream.material.boundary;const auto& original=upstream.material.region;
        out.materialInputSha256=meshHash(mesh);
        if(out.materialInputSha256!=original.inputDigest||out.materialInputSha256!=boundary.inputDigest||!boundary.boundaryValid||!original.valid){terminal(out,Outcome::NotAdmitted,"actual material hash or immutable boundary identity mismatch");return;}
        out.identitiesChecked=true;out.sourceFaceComponents=boundary.faceComponents;
        const auto n=out.componentCount,p=out.plannedPoints,t=out.plannedTriangles;
        std::vector<std::size_t> owner(p,SIZE_MAX),local(p,SIZE_MAX),faceCounts(n,0),vertexCounts(n,0);
        for(std::size_t f=0;f<t;++f){const int id=boundary.faceComponents[f];if(id<0||std::size_t(id)>=n)throw std::logic_error("invalid exact face component ID");++faceCounts[id];
            for(int j=0;j<3;++j){const int v=mesh.triangles(f,j);if(v<0||std::size_t(v)>=p)throw std::logic_error("invalid actual face vertex index");
                if(owner[v]==SIZE_MAX){owner[v]=std::size_t(id);++vertexCounts[id];}else if(owner[v]!=std::size_t(id))throw std::logic_error("vertex shared across closed boundary components");}
        }
        if(std::find(owner.begin(),owner.end(),SIZE_MAX)!=owner.end())throw std::logic_error("actual material vertex coverage incomplete");
        out.components.reserve(n);
        for(std::size_t i=0;i<n;++i){if(!faceCounts[i]||!vertexCounts[i])throw std::logic_error("empty indexed boundary component");
            SourceCageBoundaryComponentEvidence c;c.componentId=i;c.originalSign=boundary.componentVolumeSigns[i];c.parent=original.parents[i];c.depth=original.depths[i];
            if((c.originalSign!=1&&c.originalSign!=-1)||c.originalSign!=original.originalSigns[i]||c.depth<0||c.parent< -1||c.parent>=int(n)||c.parent==int(i))throw std::logic_error("invalid original component sign/nesting");
            c.checkerReversed=c.originalSign==-1;c.sourceVertexIds.reserve(vertexCounts[i]);c.sourceFaceIds.reserve(faceCounts[i]);out.components.push_back(std::move(c));
        }
        for(std::size_t v=0;v<p;++v){auto& c=out.components[owner[v]];local[v]=c.sourceVertexIds.size();c.sourceVertexIds.push_back(v);}
        for(std::size_t f=0;f<t;++f)out.components[boundary.faceComponents[f]].sourceFaceIds.push_back(f);
        std::size_t coveredP=0,coveredT=0;for(const auto& c:out.components){coveredP=add(coveredP,c.sourceVertexIds.size());coveredT=add(coveredT,c.sourceFaceIds.size());}
        if(coveredP!=p||coveredT!=t)throw std::logic_error("component source coverage totals mismatch");out.coverageChecked=true;
        budget=std::make_shared<ExactNefBudget>(settings.maximumNefComplexity);std::vector<ExactMaterialRegion> standalone;standalone.reserve(n);
        // Four-capacity allowance covers backend evidence, retained copies and
        // serialization scratch as in preflight; never spend the whole cap on
        // a single logical copy of diagnostics.
        DiagnosticBytes diagnostics{0,settings.maximumExactDiagnosticBytes/4};
        for(auto& c:out.components){SourceCageBinaryMesh checker;checker.vertices.resize(c.sourceVertexIds.size(),3);checker.triangles.resize(c.sourceFaceIds.size(),3);
            for(std::size_t v=0;v<c.sourceVertexIds.size();++v)checker.vertices.row(v)=mesh.vertices.row(c.sourceVertexIds[v]);
            for(std::size_t f=0;f<c.sourceFaceIds.size();++f)for(int j=0;j<3;++j)checker.triangles(f,j)=int(local[mesh.triangles(c.sourceFaceIds[f],j)]);
            c.originalInputSha256=meshHash(checker);if(c.checkerReversed)checker.triangles.col(1).swap(checker.triangles.col(2));c.checkerInputSha256=meshHash(checker);
            standalone.push_back(makeExactMaterialRegionWithBudget(checker.vertices,checker.triangles,budget));const auto& value=standalone.back();diagnostics.region(value);
            c.region.region=value.facts();c.region.boundary=value.boundaryFacts();
            if(!value.facts().completed){terminal(out,Outcome::Unknown,value.facts().error.c_str());return;}
            if(!value.facts().valid||value.facts().inputDigest!=c.checkerInputSha256||value.facts().materialComponents!=1||value.facts().boundedVoidComponents!=0||value.facts().originalSigns!=std::vector<int>{1}||!value.boundaryFacts().closed){terminal(out,Outcome::Rejected,"standalone exact component is not one bounded positive region");return;}
        }
        out.pairs.reserve(out.plannedPairCount);bool incomplete=false,allClear=true,contact=false;std::optional<CGAL::Gmpq> minimum;std::optional<std::size_t> minimumIndex;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j){SourceCageComponentPairEvidence pair;pair.componentIds={i,j};
            pair.clearance=checkExactBoundaryClearance(standalone[i],standalone[j],settings.minimumComponentSeparation);diagnostics.pair(pair.clearance);const auto& c=pair.clearance;
            const bool measured=c.completed&&c.state=="COMPLETE"&&c.applicable==true&&c.squaredDistanceExact&&c.witnessFaces&&c.contact&&c.atLeastMinimum;
            if(!measured){incomplete=true;}else{
                pair.parentWitnessFaces=std::array<std::size_t,2>{out.components[i].sourceFaceIds.at((*c.witnessFaces)[0]),out.components[j].sourceFaceIds.at((*c.witnessFaces)[1])};
                const CGAL::Gmpq distance(c.squaredDistanceExact->c_str());if(distance<CGAL::Gmpq(0))throw std::logic_error("negative exact squared distance");
                if(!minimum||distance<*minimum){minimum=distance;minimumIndex=out.pairs.size();}
                contact=contact||*c.contact;allClear=allClear&&!*c.contact&&*c.atLeastMinimum;
            }
            out.pairs.push_back(std::move(pair));
        }
        if(out.pairs.size()!=out.plannedPairCount||incomplete){terminal(out,Outcome::Unknown,"not every component pair has complete exact measurement");return;}
        if(n==1){out.state="NOT_APPLICABLE";out.applicable=false;}
        else{
            if(!minimumIndex)throw std::logic_error("missing complete global minimum");const auto& pair=out.pairs[*minimumIndex];
            diagnostics.put(*pair.clearance.squaredDistanceExact);out.state="MEASURED";out.applicable=true;out.contact=contact;out.atLeastMinimum=allClear;
            out.minimumPairIndex=minimumIndex;out.minimumSquaredDistanceExact=pair.clearance.squaredDistanceExact;out.minimumComponentIds=pair.componentIds;
            out.minimumLocalWitnessFaces=pair.clearance.witnessFaces;out.minimumParentWitnessFaces=pair.parentWitnessFaces;
        }
        if(!allClear){terminal(out,Outcome::Rejected,"complete distinct boundary component separation rejected");return;}
        terminal(out,Outcome::Validated,"complete distinct boundary component separation passed; full manufacturing thickness NOT_ESTABLISHED");
        data->validated=std::shared_ptr<const SourceCageComponentClearance>(new SourceCageComponentClearance(audit));
    };
    try{run();}catch(const std::length_error& e){terminal(out,Outcome::ResourceLimit,e.what());}catch(const std::bad_alloc&){terminal(out,Outcome::ResourceLimit,"component allocation failure");}
    catch(const std::exception& e){terminal(out,Outcome::Unknown,e.what());}catch(...){terminal(out,Outcome::Unknown,"unknown component clearance failure");}
    if(budget){out.nefBudget=budget->facts();if(out.nefBudget.exceeded){data->validated.reset();terminal(out,Outcome::ResourceLimit,"shared component Nef budget exceeded");}}
    finish(out,bool(data->validated));return SourceCageComponentClearanceAudit(data);
}
}
