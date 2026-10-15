#include "detail/RankTwoArcRegistrySeams.h"
#include "detail/RankTwoRegistryContinuation.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <set>
#include <sstream>
#include <tuple>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using C=detail::CertifiedAlgebraicReal;
[[noreturn]] void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]] void topology(const std::string& why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
std::string hash(const std::string& bytes){
    verification::Sha256 value;
    value.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return value.finish();
}
void coordinateRecord(std::ostream& stream,const C& point){
    const auto& c=point.certificate();
    stream<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& coefficient:c.coefficients)stream<<coefficient<<',';
    stream<<'\n';
}
bool finite(const RankTwoArcRecord& arc){
    return arc.ends[0].location==RankTwoArcEndLocation::Finite
        && arc.ends[1].location==RankTwoArcEndLocation::Finite;
}
CORE::Expr evaluateExact(const std::array<Q,25>& p,const std::array<C,2>& xy){
    CORE::Expr value(0);
    for(int i=4;i>=0;--i){CORE::Expr y(0);
        for(int j=4;j>=0;--j){std::ostringstream q;q.imbue(std::locale::classic());q<<p[5*i+j];
            y=y*xy[1].expression()+CORE::Expr(CORE::BigRat(q.str().c_str()));}
        value=value*xy[0].expression()+y;
    }
    return value;
}
bool boundEndpointProof(const RankTwoRegistryPoint& p,const std::string& n,
    const std::string& proof,const std::string& gate){
    if(p.endpointNeighborhoodIdentity==n&&p.endpointArcIdentity==proof&&p.gateIdentity==gate)return true;
    return std::any_of(p.additionalEndpointBindings.begin(),p.additionalEndpointBindings.end(),[&](const auto& b){
        return b.neighborhoodIdentity==n&&b.proofIdentity==proof&&b.gateIdentity==gate;});
}
void bindFanProof(RankTwoRegistryPoint& p,const std::string& proof,const std::string& fan){
    if(p.endpointArcIdentity==proof){p.endpointFanIdentity=fan;return;}
    for(auto& b:p.additionalEndpointBindings)if(b.proofIdentity==proof){b.fanIdentity=fan;return;}
    invalid("fan proof absent from actual point bindings");
}
}
struct RankTwoArcRegistry::Impl {
    std::string authority,record,hash;
    std::vector<RankTwoRegistryPoint> points;
    std::vector<RankTwoRegistryEdge> edges;
    std::vector<RankTwoRegistryArc> arcs;
    std::vector<RankTwoRegistryRegion> regions;
    std::vector<RankTwoEndpointFan> fans;
    std::vector<RankTwoEndpointNeighborhood> neighborhoods;
    std::vector<RankTwoArcSubdivision> subdivisions;
    std::vector<RankTwoOutsideIncidence> outsideIncidences;
    std::shared_ptr<const RankTwoSharedFiniteSeams> seamAuthority;
    std::string seamIdentity;
    std::size_t seamFace=0;
    RankTwoArcSamplingPlan samplingPlan;
    std::shared_ptr<const RankTwoRegistrySeamImport> seamImport;
    int isolationBudget=512;
};
RankTwoArcRegistry::RankTwoArcRegistry(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoArcRegistry::arrangementIdentity()const{return impl_->authority;}
const std::vector<RankTwoRegistryPoint>& RankTwoArcRegistry::points()const{return impl_->points;}
const std::vector<RankTwoRegistryEdge>& RankTwoArcRegistry::edges()const{return impl_->edges;}
const std::vector<RankTwoRegistryArc>& RankTwoArcRegistry::arcs()const{return impl_->arcs;}
const std::vector<RankTwoRegistryRegion>& RankTwoArcRegistry::regions()const{return impl_->regions;}
const std::vector<RankTwoEndpointFan>& RankTwoArcRegistry::endpointFans()const{return impl_->fans;}
const std::vector<RankTwoEndpointNeighborhood>& RankTwoArcRegistry::endpointNeighborhoods()const{return impl_->neighborhoods;}
const std::vector<RankTwoArcSubdivision>& RankTwoArcRegistry::subdivisions()const{return impl_->subdivisions;}
const std::vector<RankTwoOutsideIncidence>& RankTwoArcRegistry::outsideIncidences()const{return impl_->outsideIncidences;}
const std::string& RankTwoArcRegistry::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoArcRegistry::scientificHash()const{return impl_->hash;}

std::vector<std::size_t> RankTwoArcRegistry::loopVertices(const std::string& region,
    bool inner,std::size_t index)const{
    const RankTwoRegistryRegion* selected=nullptr;
    for(const auto& r:impl_->regions)if(r.identity==region)selected=&r;
    if(!selected)invalid("region absent from shared registry");
    const auto& loops=inner?selected->innerLoops:selected->outerLoops;
    if(index>=loops.size())invalid("loop absent from shared region");
    const auto& loop=loops[index];
    if(!loop.finite)throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
        "infinite native loop requires certified finite clipping, not truncation");
    if(loop.edges.empty())topology("finite native loop has no shared edges");
    std::vector<std::size_t> points;points.reserve(loop.edges.size());
    for(std::size_t i=0;i<loop.edges.size();++i){
        const auto& ref=loop.edges[i];const auto& next=loop.edges[(i+1)%loop.edges.size()];
        const auto& edge=impl_->edges.at(ref.edgeIndex);
        const auto& successor=impl_->edges.at(next.edgeIndex);
        const auto begin=edge.points[ref.reversed?1:0],end=edge.points[ref.reversed?0:1];
        if(begin==end || end!=successor.points[next.reversed?1:0])
            topology("shared native ring is disconnected or has a zero dependency edge");
        points.push_back(begin);
    }
    return points;
}

RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const std::vector<RankTwoEndpointNeighborhood>&,const RankTwoSourceLimits&);
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement&,
    const RankTwoArcSamplingPlan&,const std::vector<RankTwoEndpointFan>&,
    const std::vector<RankTwoEndpointNeighborhood>&,
    const std::vector<RankTwoArcSubdivision>&,const RankTwoSourceLimits&);
RankTwoArcRegistry buildRankTwoArcRegistryWithSubdivisions(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const std::vector<RankTwoArcSubdivision>& subdivisions,const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistry(arrangement,requested,{},neighborhoods,subdivisions,limits);
}
RankTwoArcRegistry rebuildRankTwoArcRegistryWithSubdivisions(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& authority,const RankTwoArcSamplingPlan& requested,
    const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const std::vector<RankTwoArcSubdivision>& subdivisions,const RankTwoSourceLimits& limits){
    if(arrangement.scientificHash()!=authority.arrangementIdentity())invalid("continuation has foreign native arrangement");
    auto normalized=requested;for(auto& entry:normalized)std::sort(entry.second.begin(),entry.second.end());
    if(normalized!=authority.impl_->samplingPlan)invalid("continuation changed frozen primitive sampling plan");
    std::set<std::string> original,current,neighborhoodFans;
    for(const auto& n:authority.endpointNeighborhoods())original.insert(n.scientificHash());
    for(const auto& n:neighborhoods){if(!current.insert(n.scientificHash()).second)invalid("duplicate continuation neighborhood");
        for(const auto& f:n.fans())neighborhoodFans.insert(f.scientificHash());}
    if(original!=current)invalid("continuation changed complete endpoint authorities");
    for(const auto& old:authority.subdivisions()){
        const auto found=std::find_if(subdivisions.begin(),subdivisions.end(),[&](const auto& s){
            return old.arcIdentity==s.arcIdentity&&old.pointIdentities==s.pointIdentities&&old.fraction==s.fraction;});
        if(found==subdivisions.end())invalid("continuation dropped an existing moving subdivision dependency");
    }
    std::vector<RankTwoEndpointFan> standalone;
    for(const auto& f:authority.endpointFans())if(!neighborhoodFans.count(f.scientificHash()))standalone.push_back(f);
    return buildRankTwoArcRegistryImported(arrangement,requested,standalone,neighborhoods,subdivisions,
        authority.impl_->seamImport.get(),limits,&authority);
}
RankTwoArcRegistry buildRankTwoArcRegistryWithNeighborhoods(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistry(arrangement,requested,{},neighborhoods,limits);
}
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistry(arrangement,requested,{},limits);
}
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointFan>& requestedFans,
    const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistry(arrangement,requested,requestedFans,{},limits);
}
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointFan>& requestedFans,
    const std::vector<RankTwoEndpointNeighborhood>& requestedNeighborhoods,const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistry(arrangement,requested,requestedFans,requestedNeighborhoods,{},limits);
}
RankTwoArcRegistry buildRankTwoArcRegistry(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointFan>& requestedFans,
    const std::vector<RankTwoEndpointNeighborhood>& requestedNeighborhoods,
    const std::vector<RankTwoArcSubdivision>& requestedSubdivisions,const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistryImported(arrangement,requested,requestedFans,requestedNeighborhoods,requestedSubdivisions,nullptr,limits);
}
RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointFan>& requestedFans,
    const std::vector<RankTwoEndpointNeighborhood>& requestedNeighborhoods,
    const std::vector<RankTwoArcSubdivision>& requestedSubdivisions,
    const RankTwoRegistrySeamImport* imported,const RankTwoSourceLimits& limits){
    return buildRankTwoArcRegistryImported(arrangement,requested,requestedFans,requestedNeighborhoods,
        requestedSubdivisions,imported,limits,nullptr);
}
RankTwoArcRegistry buildRankTwoArcRegistryImported(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::vector<RankTwoEndpointFan>& requestedFans,
    const std::vector<RankTwoEndpointNeighborhood>& requestedNeighborhoods,
    const std::vector<RankTwoArcSubdivision>& requestedSubdivisions,
    const RankTwoRegistrySeamImport* imported,const RankTwoSourceLimits& limits,
    const RankTwoArcRegistry* reuse){
    detail::requireExactAlgebraicEnvironment();
    detail::rankTwoDiagnostic("registry:validate-plan");
    if(limits.maximumIsolationBits<1 || limits.maximumIsolationBits>512 || !limits.maximumVertices)
        invalid("positive registry vertex and valid isolation budgets required");
    const auto& records=arrangementRecords(arrangement);
    std::map<std::string,RankTwoEndpointFan> fans;
    for(const auto& fan:requestedFans){
        if(fan.gate().arrangementIdentity()!=arrangement.scientificHash()
            || !fans.emplace(fan.scientificHash(),fan).second)
            invalid("endpoint fan is foreign or duplicated in registry");
    }
    std::map<std::string,const RankTwoArcRecord*> arcs;
    std::map<std::string,const RankTwoVertexRecord*> events;
    for(const auto& arc:records.arcs)if(!arcs.emplace(arc.identity,&arc).second)
        topology("duplicate authoritative arc identity");
    for(const auto& v:records.vertices)if(!events.emplace(v.identity,&v).second)
        topology("duplicate authoritative event identity");
    std::map<std::string,RankTwoEndpointNeighborhood> neighborhoods;
    std::map<std::string,std::string> fanNeighborhood;
    std::map<std::string,const RankTwoRegionRecord*> originalRegions;
    for(const auto& r:records.regions)if(!originalRegions.emplace(r.identity,&r).second)
        topology("duplicate authoritative region identity");
    std::vector<RankTwoOutsideIncidence> outside;
    std::set<std::pair<std::string,std::string>> outsideKeys;
    for(const auto& n:requestedNeighborhoods){
        if(n.arrangementIdentity()!=arrangement.scientificHash() || !events.count(n.endpointIdentity())
            || !neighborhoods.emplace(n.endpointIdentity(),n).second)
            invalid("complete endpoint neighborhood is foreign or duplicated");
        for(const auto& id:n.incidentArcs()){
            const auto arc=arcs.find(id);
            if(arc==arcs.end())invalid("complete neighborhood names foreign arc");
            if(!finite(*arc->second)){
                const auto& native=*arc->second;
                const auto ray=std::find_if(n.boundaryRays().begin(),n.boundaryRays().end(),
                    [&](const auto& r){return r.proof.arcIdentity==id;});
                const auto left=originalRegions.find(native.leftRegion),right=originalRegions.find(native.rightRegion);
                const bool ownsD=std::any_of(native.owners.begin(),native.owners.end(),
                    [](const auto& owner){return owner.polynomial=="D";});
                // The immutable DCEL is split at ALL D intersections. A ray
                // not owned by D has constant nonzero D sign on its open arc.
                // Both adjacent original cell signs and an actual certified
                // on-arc gate independently select the negative side. Never
                // infer this from Outside (which also includes cone filtering).
                if(ray==n.boundaryRays().end()||left==originalRegions.end()||right==originalRegions.end()
                    ||ownsD||left->second->originalSigns[2]!=-1||right->second->originalSigns[2]!=-1
                    ||ray->proof.arrangementIdentity!=arrangement.scientificHash()
                    ||ray->proof.endpointIdentity!=n.endpointIdentity()
                    ||ray->proof.gateIdentity!=ray->gate.identity()
                    ||ray->proof.strip.chordDistanceUpperBound!=0
                    ||evaluateExact(records.originalFhD[2],ray->proof.gateCoordinates).sign()!=-1)
                    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
                        "complete neighborhood unbounded incidence lacks exact original D-negative ray exclusion");
                outsideKeys.emplace(n.scientificHash(),id);
                outside.push_back({id,n.endpointIdentity(),n.scientificHash(),ray->proof.identity,
                    ray->gate.identity(),native.leftRegion,native.rightRegion});
            }
        }
        for(const auto& fan:n.fans()){
            if(!fans.emplace(fan.scientificHash(),fan).second)
                invalid("complete neighborhood fan duplicates another supplied proof");
            fanNeighborhood.emplace(fan.scientificHash(),n.scientificHash());
        }
    }
    auto plan=requested;
    for(auto& entry:plan){
        const auto found=arcs.find(entry.first);
        if(found==arcs.end() || !finite(*found->second))invalid("sample plan names foreign or unbounded arc");
        std::sort(entry.second.begin(),entry.second.end());
        for(std::size_t i=0;i<entry.second.size();++i)
            if(entry.second[i]<=0 || entry.second[i]>=1 || (i && entry.second[i]==entry.second[i-1]))
                invalid("sample plan has endpoint or duplicated interior fraction");
    }
    // A proof is not a second geometric vertex. Alias ONLY a literal native
    // affine dependency (same arc/branch, both actual ends, axis and fraction),
    // never two roots which merely coincide for the current coefficients.
    std::map<std::string,std::string> primitiveAliasProofs;
    const auto aliasPrimitive=[&](const RankTwoSharedGateArcCertificate& proof,const RankTwoMovingGate& gate){
        const auto found=arcs.find(proof.arcIdentity);
        if(found==arcs.end()||!finite(*found->second)
            ||(imported&&imported->chains_.count(proof.arcIdentity)))return;
        const auto& arc=*found->second;
        if(gate.projectionAxis()!=(arc.vertical?1:0))return;
        Q f;
        if(gate.endpointIdentity()==*arc.ends[0].vertexIdentity&&gate.anchorIdentity()==*arc.ends[1].vertexIdentity)
            f=gate.fraction();
        else if(gate.endpointIdentity()==*arc.ends[1].vertexIdentity&&gate.anchorIdentity()==*arc.ends[0].vertexIdentity)
            f=Q(1)-gate.fraction();
        else return;
        if(!std::binary_search(plan.at(arc.identity).begin(),plan.at(arc.identity).end(),f))return;
        std::ostringstream row;row.imbue(std::locale::classic());
        row<<"rank-two-moving-arc-sample-v1\n"<<arrangement.scientificHash()<<'\n'<<arc.identity<<'\n'
            <<*arc.ends[0].vertexIdentity<<'\n'<<*arc.ends[1].vertexIdentity<<'\n'
            <<f.numerator()<<'/'<<f.denominator()<<'\n'<<arc.vertical<<' '<<arc.branch<<'\n';
        if(!primitiveAliasProofs.emplace(proof.identity,hash(row.str())).second)
            invalid("duplicate literal primitive endpoint proof");
    };
    for(const auto& f:fans)for(const auto& p:f.second.orderedArcs())aliasPrimitive(p,f.second.gate());
    for(const auto& n:neighborhoods)for(const auto& r:n.second.boundaryRays())aliasPrimitive(r.proof,r.gate);
    std::size_t count=events.size();
    if(count>limits.maximumVertices)throw RankTwoSourceError("RESOURCE_LIMIT","registry event count exceeds vertex budget");
    for(const auto& entry:arcs)if(finite(*entry.second)){
        const auto p=plan.find(entry.first);
        if(p==plan.end())invalid("sample plan omitted a finite shared arc");
        if(p->second.size()>limits.maximumVertices-count)
            throw RankTwoSourceError("RESOURCE_LIMIT","registry moving points exceed vertex budget");
        count+=p->second.size();
    }
    for(const auto& entry:fans){
        std::size_t additions=0;
        for(const auto& proof:entry.second.orderedArcs())
            if(!primitiveAliasProofs.count(proof.identity)
                &&(!imported||!imported->endpointAliasProofs_.count(proof.identity)))++additions;
        if(additions>limits.maximumVertices-count)
            throw RankTwoSourceError("RESOURCE_LIMIT","registry common-gate points exceed vertex budget");
        count+=additions;
    }
    for(const auto& entry:neighborhoods){
        std::size_t additions=0;
        for(const auto& ray:entry.second.boundaryRays())
            if(!outsideKeys.count({entry.second.scientificHash(),ray.proof.arcIdentity})
                &&!primitiveAliasProofs.count(ray.proof.identity)
                &&(!imported||!imported->endpointAliasProofs_.count(ray.proof.identity)))++additions;
        if(additions>limits.maximumVertices-count)
            throw RankTwoSourceError("RESOURCE_LIMIT","registry boundary rays exceed vertex budget");
        count+=additions;
    }
    std::map<std::string,const RankTwoRegistryPoint*> inherited;
    if(reuse){
        if(reuse->arrangementIdentity()!=arrangement.scientificHash()
            ||reuse->impl_->samplingPlan!=plan||reuse->impl_->seamImport.get()!=imported)
            invalid("coordinate reuse changed native, primitive or shared seam authority");
        std::set<std::string> oldFans,newFans,oldNeighborhoods,newNeighborhoods;
        for(const auto& f:reuse->endpointFans())oldFans.insert(f.scientificHash());
        for(const auto& f:fans)newFans.insert(f.first);
        for(const auto& n:reuse->endpointNeighborhoods())oldNeighborhoods.insert(n.scientificHash());
        for(const auto& n:neighborhoods)newNeighborhoods.insert(n.second.scientificHash());
        if(oldFans!=newFans||oldNeighborhoods!=newNeighborhoods)
            invalid("coordinate reuse changed complete endpoint authority");
        if(limits.maximumIsolationBits<reuse->impl_->isolationBudget)
            throw RankTwoSourceError("RESOURCE_LIMIT","coordinate reuse cannot bypass a tighter isolation budget");
        for(const auto& old:reuse->subdivisions())if(std::none_of(requestedSubdivisions.begin(),requestedSubdivisions.end(),
            [&](const auto& s){return old.arcIdentity==s.arcIdentity&&old.pointIdentities==s.pointIdentities&&old.fraction==s.fraction;}))
            invalid("coordinate reuse dropped a sealed moving dependency");
        for(const auto& p:reuse->points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision)
            if(!inherited.emplace(p.identity,&p).second)invalid("duplicate inherited subdivision identity");
    }
    if(requestedSubdivisions.size()>limits.maximumVertices-count)
        throw RankTwoSourceError("RESOURCE_LIMIT","registry subdivision DAG exceeds vertex budget");
    count+=requestedSubdivisions.size();
    auto out=std::make_shared<RankTwoArcRegistry::Impl>();out->authority=arrangement.scientificHash();
    out->samplingPlan=plan;
    out->isolationBudget=limits.maximumIsolationBits;
    std::sort(outside.begin(),outside.end(),[](const auto& a,const auto& b){
        return std::tie(a.neighborhoodIdentity,a.arcIdentity)<std::tie(b.neighborhoodIdentity,b.arcIdentity);});
    out->outsideIncidences=std::move(outside);
    if(imported){
        if(imported->faceAuthority_!=out->authority||!imported->authority_)invalid("foreign immutable seam import");
        out->seamAuthority=imported->authority_;out->seamIdentity=imported->masterIdentity_;out->seamFace=imported->face_;
        out->seamImport=std::make_shared<const RankTwoRegistrySeamImport>(*imported);
        const auto seam=[&](const std::string& arc){return imported->chains_.count(arc)!=0;};
        for(const auto& s:requestedSubdivisions)if(seam(s.arcIdentity))
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","local subdivision would change published common seam; new master required");
        for(const auto& f:fans)for(const auto& p:f.second.orderedArcs())if(seam(p.arcIdentity)&&!imported->endpointBindings_.count(p.identity))
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","local fan gate would change published common seam; certified composition required");
        for(const auto& n:neighborhoods)for(const auto& ray:n.second.boundaryRays())if(seam(ray.proof.arcIdentity)&&!imported->endpointBindings_.count(ray.proof.identity))
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","local endpoint ray would change published common seam; certified composition required");
    }
    std::map<std::string,RankTwoRegistryPoint> pending;
    for(const auto& entry:events){
        const auto& v=*entry.second;
        pending.emplace(v.identity,RankTwoRegistryPoint{v.identity,v.identity,{},RankTwoRegistryPointKind::Event,
            {C::reload(v.coordinates[0],limits.maximumIsolationBits),C::reload(v.coordinates[1],limits.maximumIsolationBits)},
            {},Q(0)});
    }
    std::map<std::string,std::vector<std::string>> chains;
    for(const auto& entry:arcs){
        const auto& arc=*entry.second;
        if(!finite(arc))continue;
        for(const auto& end:arc.ends)if(!end.vertexIdentity || !events.count(*end.vertexIdentity))
            topology("finite authoritative arc lacks a known endpoint");
        auto& chain=chains[arc.identity];chain.push_back(*arc.ends[0].vertexIdentity);
        if(imported&&imported->chains_.count(arc.identity)){
            chain.clear();const auto& supplied=imported->chains_.at(arc.identity);
            const auto extra=imported->endpointCounts_.find(arc.identity);
            if(supplied.size()!=plan.at(arc.identity).size()+2+(extra==imported->endpointCounts_.end()?0:extra->second))
                invalid("master chain count differs from native plan and published endpoint gates");
            for(const auto& p:supplied){
                const auto existing=pending.find(p.identity);
                if(p.kind==RankTwoRegistryPointKind::Event){
                    if(existing==pending.end()||existing->second.kind!=RankTwoRegistryPointKind::Event)
                        invalid("master endpoint absent from native events");
                    for(int d=0;d<2;++d)if(existing->second.coordinates[d].expression()!=p.coordinates[d].expression())
                        invalid("master and native endpoint coordinates differ");
                    existing->second.coordinates=p.coordinates; // actual master's immutable expressions
                }else if(!pending.emplace(p.identity,p).second)topology("master seam point identity collision");
                chain.push_back(p.identity);
            }
            continue;
        }
        for(const auto& fraction:plan.at(arc.identity)){
            const auto phase="registry:sample:"+arc.identity;
            detail::rankTwoDiagnostic(phase.c_str());
            auto sample=sampleRankTwoArc(arrangement,arc.identity,fraction,limits.maximumIsolationBits);
            const auto id=sample.identity;
            if(!pending.emplace(id,RankTwoRegistryPoint{id,{},arc.identity,RankTwoRegistryPointKind::ArcSample,
                std::move(sample.coordinates),sample.endpointIdentities,fraction}).second)
                topology("registry sample dependency identity collides with another point");
            chain.push_back(id);
        }
        chain.push_back(*arc.ends[1].vertexIdentity);
    }
    struct FanEdge {std::string arc,start,gate,fan,neighborhood;};
    const auto bindPrimitive=[&](const std::string& id,const RankTwoSharedGateArcCertificate& proof,
        const RankTwoMovingGate& gate,const std::string& neighborhood,const std::string& fan){
        const auto point=pending.find(id);
        if(point==pending.end()||point->second.kind!=RankTwoRegistryPointKind::ArcSample
            ||point->second.arcIdentity!=proof.arcIdentity)
            invalid("literal endpoint alias lost its original primitive sample");
        for(int d=0;d<2;++d)if(point->second.coordinates[d].expression()!=proof.gateCoordinates[d].expression())
            invalid("literal primitive dependency differs from actual certified gate");
        point->second.additionalEndpointBindings.push_back({neighborhood,proof.identity,gate.identity(),fan});
    };
    detail::rankTwoDiagnostic("registry:insert-fans");
    std::vector<FanEdge> fanEdges;
    for(const auto& entry:fans){
        const auto& fan=entry.second;const auto& gate=fan.gate();
        out->fans.push_back(fan);
        const auto found=fanNeighborhood.find(fan.scientificHash());
        const std::string neighborhood=found==fanNeighborhood.end()?std::string{}:found->second;
        for(const auto& proof:fan.orderedArcs()){
            const auto arc=arcs.find(proof.arcIdentity);
            if(arc==arcs.end() || !events.count(proof.endpointIdentity)
                || proof.arrangementIdentity!=out->authority || proof.gateIdentity!=gate.identity())
                invalid("fan point has a different actual arc, endpoint or gate authority");
            const auto published=imported?imported->endpointBindings_.find(proof.identity):std::map<std::string,std::string>::const_iterator{};
            const bool shared=imported&&published!=imported->endpointBindings_.end();
            const auto primitive=primitiveAliasProofs.find(proof.identity);
            const bool localAlias=primitive!=primitiveAliasProofs.end();
            const auto id=shared?published->second:localAlias?primitive->second:hash("rank-two-registry-fan-point-v1\n"+proof.identity+"\n");
            if(shared){const auto point=pending.find(id);
                if(point==pending.end()||point->second.kind!=RankTwoRegistryPointKind::SharedSeam
                    ||point->second.arcIdentity!=proof.arcIdentity
                    ||!boundEndpointProof(point->second,neighborhood,proof.identity,gate.identity()))
                    invalid("published fan point lost actual complete native proof authority");
                for(int d=0;d<2;++d)if(point->second.coordinates[d].expression()!=proof.gateCoordinates[d].expression())
                    invalid("published common endpoint differs from actual fan proof");
                bindFanProof(point->second,proof.identity,fan.scientificHash());
            }else if(localAlias)bindPrimitive(id,proof,gate,neighborhood,fan.scientificHash());
            else if(!pending.emplace(id,RankTwoRegistryPoint{id,{},proof.arcIdentity,RankTwoRegistryPointKind::SharedGate,
                proof.gateCoordinates,{gate.endpointIdentity(),gate.anchorIdentity()},gate.fraction(),
                gate.identity(),fan.scientificHash(),gate.projectionAxis(),neighborhood,
                neighborhood.empty()?std::string{}:proof.identity}).second)
                topology("registry fan point identity collides with another dependency");
            if(finite(*arc->second)){
                if(!shared&&!localAlias)chains.at(proof.arcIdentity).push_back(id);
                fanEdges.push_back({proof.arcIdentity,proof.endpointIdentity,id,fan.scientificHash(),neighborhood});
            }
        }
    }
    for(const auto& entry:neighborhoods){
        const auto& n=entry.second;out->neighborhoods.push_back(n);
        for(const auto& ray:n.boundaryRays()){
            const auto& proof=ray.proof;const auto& gate=ray.gate;
            if(outsideKeys.count({n.scientificHash(),proof.arcIdentity}))continue;
            if(proof.arrangementIdentity!=out->authority || proof.endpointIdentity!=n.endpointIdentity()
                || proof.gateIdentity!=gate.identity() || proof.strip.chordDistanceUpperBound!=0)
                invalid("complete boundary ray has inconsistent authority or is not unchanged");
            const auto published=imported?imported->endpointBindings_.find(proof.identity):std::map<std::string,std::string>::const_iterator{};
            const bool shared=imported&&published!=imported->endpointBindings_.end();
            const auto primitive=primitiveAliasProofs.find(proof.identity);
            const bool localAlias=primitive!=primitiveAliasProofs.end();
            const auto id=shared?published->second:localAlias?primitive->second:hash("rank-two-registry-boundary-ray-point-v1\n"+proof.identity+"\n");
            if(shared){const auto point=pending.find(id);
                if(point==pending.end()||point->second.kind!=RankTwoRegistryPointKind::SharedSeam
                    ||point->second.arcIdentity!=proof.arcIdentity
                    ||!boundEndpointProof(point->second,n.scientificHash(),proof.identity,gate.identity()))
                    invalid("published endpoint point lost actual native proof authority");
                for(int d=0;d<2;++d)if(point->second.coordinates[d].expression()!=proof.gateCoordinates[d].expression())
                    invalid("published common endpoint differs from actual ray proof");
            }else if(localAlias)bindPrimitive(id,proof,gate,n.scientificHash(),{});
            else if(!pending.emplace(id,RankTwoRegistryPoint{id,{},proof.arcIdentity,RankTwoRegistryPointKind::BoundaryRay,
                proof.gateCoordinates,{gate.endpointIdentity(),gate.anchorIdentity()},gate.fraction(),
                gate.identity(),{},gate.projectionAxis(),n.scientificHash(),proof.identity}).second)
                topology("boundary-ray dependency identity collides with another point");
            if(!shared&&!localAlias)chains.at(proof.arcIdentity).push_back(id);
            fanEdges.push_back({proof.arcIdentity,proof.endpointIdentity,id,{},n.scientificHash()});
        }
    }
    std::map<std::string,RankTwoArcSubdivision> splits;
    for(const auto& s:requestedSubdivisions){
        const auto arc=arcs.find(s.arcIdentity);
        if(arc==arcs.end() || !finite(*arc->second) || s.fraction<=0 || s.fraction>=1
            || s.pointIdentities[0].empty() || s.pointIdentities[1].empty()
            || s.pointIdentities[0]==s.pointIdentities[1])
            invalid("subdivision needs an actual finite arc, distinct parents and interior fraction");
        std::ostringstream row;row.imbue(std::locale::classic());
        row<<"rank-two-moving-arc-subdivision-v1\n"<<out->authority<<'\n'<<s.arcIdentity<<'\n'
            <<s.pointIdentities[0]<<'\n'<<s.pointIdentities[1]<<'\n'
            <<s.fraction.numerator()<<'/'<<s.fraction.denominator()<<'\n';
        if(!splits.emplace(hash(row.str()),s).second)invalid("duplicate subdivision dependency");
    }
    // Request order cannot select scientific vertex order or refinement state.
    // At each wave, stable IDs order all currently evaluable DAG dependencies.
    // Preserve rational affine dependence on the ORIGINAL native endpoints
    // (and opaque gates). Primitive ArcSample projections already interpolate
    // those endpoints; treating their correlated roots as independent leaves
    // invents conjugate combinations and unnecessary resultant degree.
    // Sibling samples are not independent conjugate variables: repeated
    // elimination between their coordinate polynomials would invent growing
    // degree even though their projections remain in the same affine span.
    using ProjectionWeights=std::map<std::string,Q>;
    std::map<std::string,ProjectionWeights> projectionWeights;
    for(const auto& entry:pending){
        const auto& p=entry.second;
        if(p.kind!=RankTwoRegistryPointKind::ArcSample)continue;
        const auto found=arcs.find(p.arcIdentity);
        if(found==arcs.end() || !finite(*found->second) || p.fraction<=0 || p.fraction>=1)
            invalid("primitive sample lacks finite native affine authority");
        const auto& arc=*found->second;
        for(int end=0;end<2;++end){
            if(p.endpointIdentities[end]!=*arc.ends[end].vertexIdentity
                || !pending.count(p.endpointIdentities[end])
                || pending.at(p.endpointIdentities[end]).kind!=RankTwoRegistryPointKind::Event)
                invalid("primitive sample affine endpoints differ from actual native arc");
        }
        projectionWeights.emplace(entry.first,ProjectionWeights{
            {p.endpointIdentities[0],Q(1)-p.fraction},
            {p.endpointIdentities[1],p.fraction}});
    }
    const auto weights=[&](const std::string& id)->ProjectionWeights{
        const auto found=projectionWeights.find(id);
        if(found!=projectionWeights.end())return found->second;
        if(pending.at(id).kind==RankTwoRegistryPointKind::ArcSubdivision)
            invalid("subdivision projection dependency was not evaluated");
        return {{id,Q(1)}};
    };
    while(!splits.empty()){
        bool progress=false;
        for(auto i=splits.begin();i!=splits.end();){
            const auto& s=i->second;
            const auto l=pending.find(s.pointIdentities[0]),r=pending.find(s.pointIdentities[1]);
            if(l==pending.end() || r==pending.end()){++i;continue;}
            const auto& arc=*arcs.at(s.arcIdentity);
            const auto onArc=[&](const RankTwoRegistryPoint& p){
                if(p.kind!=RankTwoRegistryPointKind::Event)return p.arcIdentity==arc.identity;
                return p.eventIdentity==arc.ends[0].vertexIdentity || p.eventIdentity==arc.ends[1].vertexIdentity;
            };
            const auto subdivisionPhase="registry:subdivision:"+i->first;
            detail::rankTwoDiagnostic((subdivisionPhase+":begin").c_str());
            if(!onArc(l->second) || !onArc(r->second))invalid("subdivision parents do not belong to same actual native arc");
            const int axis=arc.vertical?1:0;
            const auto known=inherited.find(i->first);
            if(known!=inherited.end()&&(known->second->arcIdentity!=s.arcIdentity
                ||known->second->endpointIdentities!=s.pointIdentities||known->second->fraction!=s.fraction))
                invalid("inherited coordinate does not bind the same actual subdivision");
            if(known==inherited.end()&&!(l->second.coordinates[axis].expression()<r->second.coordinates[axis].expression()))
                invalid("subdivision parents not native-monotone increasing");
            ProjectionWeights combined;
            for(const auto& p:weights(l->first))combined[p.first]+=(1-s.fraction)*p.second;
            for(const auto& p:weights(r->first))combined[p.first]+=s.fraction*p.second;
            Q total(0);for(const auto& p:combined){
                if(p.second<=0)invalid("moving projection has nonpositive original-leaf weight");
                total+=p.second;
            }
            if(total!=1 || combined.empty())invalid("moving projection affine weights do not sum to one");
            if(known!=inherited.end()){
                if(!pending.emplace(i->first,*known->second).second)
                    topology("inherited subdivision collides with existing source feature");
                detail::rankTwoDiagnostic((subdivisionPhase+":reused-coordinate").c_str());
            }else{
                const auto first=combined.begin();Q sum=first->second;
                auto u=pending.at(first->first).coordinates[axis];
                for(auto p=std::next(first);p!=combined.end();++p){
                    u=detail::sampleRankTwoAlgebraicCoordinate(u,pending.at(p->first).coordinates[axis],
                        p->second/(sum+p->second),limits.maximumIsolationBits);
                    sum+=p->second;
                }
                auto xy=detail::evaluateRankTwoArcProjection(arrangement,arc.identity,u,limits.maximumIsolationBits);
                if(!pending.emplace(i->first,RankTwoRegistryPoint{i->first,{},arc.identity,
                    RankTwoRegistryPointKind::ArcSubdivision,std::move(xy),s.pointIdentities,s.fraction}).second)
                    topology("subdivision identity collides with existing source feature");
                detail::rankTwoDiagnostic((subdivisionPhase+":coordinates-ready").c_str());
            }
            projectionWeights.emplace(i->first,std::move(combined));
            chains.at(arc.identity).push_back(i->first);out->subdivisions.push_back(s);
            i=splits.erase(i);progress=true;
        }
        if(!progress)invalid("subdivision DAG has missing parents or cyclic dependencies");
    }
    for(auto& entry:chains){
        const auto phase="registry:sort:"+entry.first;detail::rankTwoDiagnostic(phase.c_str());
        const int axis=arcs.at(entry.first)->vertical?1:0;
        std::sort(entry.second.begin(),entry.second.end(),[&](const auto& a,const auto& b){
            return pending.at(a).coordinates[axis].expression()<pending.at(b).coordinates[axis].expression();
        });
    }
    std::map<std::string,std::size_t> pointIndex,arcIndex;
    detail::rankTwoDiagnostic("registry:publish-edges");
    for(auto& entry:pending){pointIndex.emplace(entry.first,out->points.size());out->points.push_back(std::move(entry.second));}
    if(out->points.size()!=count)topology("registry point count differs from verified plan");
    for(const auto& entry:arcs){
        const auto& native=*entry.second;
        RankTwoRegistryArc saved{native.identity,native.supportPolynomial,native.leftRegion,native.rightRegion,
            finite(native),native.vertical,native.branch,native.ends,native.owners,{}};
        if(saved.finite){
            const auto& chain=chains.at(saved.identity);
            const int axis=native.vertical?1:0;
            for(std::size_t i=1;i<chain.size();++i){
                const std::array<std::size_t,2> indices{pointIndex.at(chain[i-1]),pointIndex.at(chain[i])};
                if(!(out->points[indices[0]].coordinates[axis].expression()
                    <out->points[indices[1]].coordinates[axis].expression()))
                    topology("registry chain is not strictly native-monotone: arc="+saved.identity
                        +" first="+chain[i-1]+" second="+chain[i]);
                const auto id=hash("rank-two-registry-edge-v1\n"+out->authority+"\n"+saved.identity
                    +"\n"+chain[i-1]+"\n"+chain[i]+"\n");
                saved.edges.push_back(out->edges.size());
                out->edges.push_back({id,saved.identity,indices,saved.owners,{}});
            }
        }
        arcIndex.emplace(saved.identity,out->arcs.size());out->arcs.push_back(std::move(saved));
    }
    for(const auto& proof:fanEdges){
        detail::rankTwoDiagnostic("registry:bind-fan-edges");
        const auto& chain=chains.at(proof.arc);
        const auto& arc=out->arcs.at(arcIndex.at(proof.arc));
        bool atStart=chain.front()==proof.start;
        if(!atStart && chain.back()!=proof.start)topology("fan endpoint is not a shared chain endpoint");
        if((atStart?chain[1]:chain[chain.size()-2])!=proof.gate)
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
                "independent sample splits common-gate incident chord; simultaneous fan proof invalid");
        auto& edge=out->edges.at(atStart?arc.edges.front():arc.edges.back());
        if(!edge.endpointFanIdentity.empty() && edge.endpointFanIdentity!=proof.fan)
            invalid("shared endpoint edge has conflicting fan authorities");
        edge.endpointFanIdentity=proof.fan;
        if(!edge.endpointNeighborhoodIdentity.empty() && edge.endpointNeighborhoodIdentity!=proof.neighborhood)
            invalid("shared endpoint edge has conflicting complete neighborhood authorities");
        edge.endpointNeighborhoodIdentity=proof.neighborhood;
    }
    for(const auto& n:out->neighborhoods){
        std::set<std::string> bound;
        for(const auto& edge:out->edges)if(edge.endpointNeighborhoodIdentity==n.scientificHash())
            if(!bound.insert(edge.arcIdentity).second)topology("complete star binds same incidence twice");
        for(const auto& excluded:out->outsideIncidences)if(excluded.neighborhoodIdentity==n.scientificHash())
            if(!bound.insert(excluded.arcIdentity).second)topology("outside incidence also bound to finite carrier edge");
        if(bound!=std::set<std::string>(n.incidentArcs().begin(),n.incidentArcs().end()))
            topology("published shared edges omit actual complete endpoint incidence");
    }
    std::map<std::string,const RankTwoRegionRecord*> regions;
    detail::rankTwoDiagnostic("registry:publish-loops");
    for(const auto& r:records.regions)if(!regions.emplace(r.identity,&r).second)
        topology("duplicate authoritative region identity");
    for(const auto& entry:regions){
        const auto& original=*entry.second;
        RankTwoRegistryRegion saved{original.identity,original.kind,original.unbounded,
            original.originalSigns,original.normalSigns,{},{},{}};
        const auto loops=[&](const auto& originals,auto& published){
            for(const auto& sourceLoop:originals){
                RankTwoRegistryLoop loop;loop.finite=true;loop.sourceArcs=sourceLoop;
                for(const auto& ref:sourceLoop){
                    if(ref.atInfinity){loop.finite=false;continue;}
                    const auto found=arcIndex.find(ref.identity);
                    if(found==arcIndex.end())topology("native region names an unknown shared arc");
                    const auto& arc=out->arcs[found->second];
                    if((ref.reversed?arc.rightRegion:arc.leftRegion)!=original.identity)
                        topology("directed native arc disagrees with region side");
                    if(!arc.finite){loop.finite=false;continue;}
                    if(ref.reversed)for(auto e=arc.edges.rbegin();e!=arc.edges.rend();++e)
                        loop.edges.push_back({*e,true});
                    else for(const auto e:arc.edges)loop.edges.push_back({e,false});
                }
                if(!loop.finite)loop.edges.clear(); // retained source ring is unresolved, never a partial polygon
                published.push_back(std::move(loop));
            }
        };
        loops(original.outerLoops,saved.outerLoops);loops(original.innerLoops,saved.innerLoops);
        for(const auto& vertex:original.isolatedVertices){
            if(!events.count(vertex))topology("native isolated feature is not an event");
            saved.isolatedPoints.push_back(pointIndex.at(vertex));
        }
        out->regions.push_back(std::move(saved));
    }
    RankTwoArcRegistry result(out);
    for(const auto& r:out->regions)for(bool inner:{false,true}){
        const auto& loops=inner?r.innerLoops:r.outerLoops;
        for(std::size_t i=0;i<loops.size();++i)if(loops[i].finite)result.loopVertices(r.identity,inner,i);
    }
    std::ostringstream record,scientific;record.imbue(std::locale::classic());scientific.imbue(std::locale::classic());
    const bool multipleProofs=std::any_of(out->points.begin(),out->points.end(),[](const auto& p){return !p.additionalEndpointBindings.empty();});
    const char* version=multipleProofs?"rank-two-shared-arc-registry-v6":
        !out->outsideIncidences.empty()?"rank-two-shared-arc-registry-v5":
        out->seamAuthority?"rank-two-shared-arc-registry-v4":!out->subdivisions.empty()?"rank-two-shared-arc-registry-v3":
        out->neighborhoods.empty()?"rank-two-shared-arc-registry-v1":"rank-two-shared-arc-registry-v2";
    record<<version<<"\nAUTHORITY "<<out->authority<<'\n';
    scientific<<version<<'\n'<<out->authority<<'\n';
    if(out->seamAuthority){
        const auto row="MASTER_SEAMS "+out->seamIdentity+" FACE "+std::to_string(out->seamFace)+"\n";
        record<<row;scientific<<row;
    }
    for(const auto& x:out->outsideIncidences){
        const auto row="OUTSIDE_INCIDENCE "+x.arcIdentity+" "+x.endpointIdentity+" "+x.neighborhoodIdentity+" "
            +x.endpointArcIdentity+" "+x.gateIdentity+" "+x.leftRegion+" "+x.rightRegion+" D -1 -1 OPEN_NATIVE_ARC\n";
        record<<row;scientific<<row;
    }
    for(std::size_t i=0;i<out->points.size();++i){
        const auto& p=out->points[i];
        std::ostringstream row;row.imbue(std::locale::classic());
        const char* kind=p.kind==RankTwoRegistryPointKind::Event?"EVENT":
            p.kind==RankTwoRegistryPointKind::ArcSample?"ARC_SAMPLE":
            p.kind==RankTwoRegistryPointKind::SharedGate?"SHARED_GATE":
            p.kind==RankTwoRegistryPointKind::BoundaryRay?"BOUNDARY_RAY":
            p.kind==RankTwoRegistryPointKind::SharedSeam?"SHARED_SEAM":"ARC_SUBDIVISION";
        row<<"POINT "<<i<<' '<<p.identity<<' '<<kind
            <<' '<<p.eventIdentity<<' '<<p.arcIdentity<<' '<<p.endpointIdentities[0]<<' '<<p.endpointIdentities[1]
            <<' '<<p.fraction<<'\n';
        record<<row.str();scientific<<row.str();
        std::ostringstream gate;gate.imbue(std::locale::classic());
        gate<<"GATE "<<p.gateIdentity<<' '<<p.endpointFanIdentity<<' '<<p.projectionAxis<<'\n';
        record<<gate.str();scientific<<gate.str();
        if(!out->neighborhoods.empty() || !out->subdivisions.empty()){
            const auto row="NEIGHBORHOOD_GATE "+p.endpointNeighborhoodIdentity+" "+p.endpointArcIdentity+"\n";
            record<<row;scientific<<row;
        }
        for(const auto& c:p.coordinates){coordinateRecord(record,c);scientific<<c.certificate().rootIdentity<<'\n';}
        for(const auto& binding:p.additionalEndpointBindings){
            const auto row="POINT_ENDPOINT_BINDING "+p.identity+" "+binding.neighborhoodIdentity+" "+binding.proofIdentity+" "
                +binding.gateIdentity+" "+(binding.fanIdentity.empty()?"NONE":binding.fanIdentity)+"\n";
            record<<row;scientific<<row;
        }
    }
    for(std::size_t i=0;i<out->edges.size();++i){
        const auto& e=out->edges[i];std::ostringstream row;
        row.imbue(std::locale::classic());row<<"EDGE "<<i<<' '<<e.identity<<' '<<e.arcIdentity<<' '<<e.points[0]<<' '<<e.points[1];
        for(const auto& owner:e.owners)row<<' '<<owner.polynomial<<':'<<owner.multiplicity;
        row<<" FAN "<<(e.endpointFanIdentity.empty()?"NONE":e.endpointFanIdentity);
        if(!out->neighborhoods.empty() || !out->subdivisions.empty())row<<" NEIGHBORHOOD "
            <<(e.endpointNeighborhoodIdentity.empty()?"NONE":e.endpointNeighborhoodIdentity);
        row<<'\n';record<<row.str();scientific<<row.str();
    }
    // The immutable arrangement hash binds support polynomials, branches and
    // all native ends/owners, including unsampled infinity. Published indices
    // and all directed consumers are additionally bound below.
    for(const auto& arc:out->arcs){
        std::ostringstream row;row.imbue(std::locale::classic());
        row<<"ARC "<<arc.identity<<' '<<arc.finite<<' '<<arc.vertical<<' '<<arc.branch;
        for(const auto e:arc.edges)row<<' '<<e;row<<'\n';record<<row.str();scientific<<row.str();
    }
    for(const auto& region:out->regions){
        std::ostringstream row;row.imbue(std::locale::classic());
        row<<"REGION "<<region.identity<<' '<<int(region.kind)<<' '<<region.unbounded;
        for(const auto s:region.originalSigns)row<<' '<<s;row<<" NORMAL";
        for(const auto s:region.normalSigns)row<<' '<<s;row<<'\n';
        const auto append=[&](const auto& loops,const char* kind){
            for(const auto& loop:loops){
                row<<"LOOP "<<kind<<' '<<loop.finite;
                for(const auto& ref:loop.sourceArcs)row<<' '<<(ref.atInfinity?"INFINITY":(ref.reversed?"-":"+")+ref.identity);
                row<<" EDGES";for(const auto& ref:loop.edges)row<<' '<<(ref.reversed?"-":"+")<<ref.edgeIndex;
                row<<'\n';
            }
        };
        append(region.outerLoops,"OUTER");append(region.innerLoops,"INNER");
        row<<"ISOLATED";for(const auto p:region.isolatedPoints)row<<' '<<p;row<<'\n';
        record<<row.str();scientific<<row.str();
    }
    for(const auto& fan:out->fans){
        const auto row="FAN_BEGIN "+fan.scientificHash()+"\n";
        record<<row<<fan.canonicalRecord()<<"FAN_END\n";
        scientific<<row;
    }
    for(const auto& n:out->neighborhoods){
        const auto row="NEIGHBORHOOD_BEGIN "+n.scientificHash()+"\n";
        record<<row<<n.canonicalRecord()<<"NEIGHBORHOOD_END\n";scientific<<row;
    }
    record<<"chord_topology NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();out->hash=hash(scientific.str());
    detail::requireExactAlgebraicEnvironment();return result;
}
} // namespace spin
