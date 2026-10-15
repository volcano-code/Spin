#include "detail/RankTwoFiniteCellTopology.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Edge=std::pair<std::size_t,std::size_t>;
[[noreturn]]void invalid(const std::string& why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
struct Budget {
    std::uint64_t used=0,limit;
    void tick(){if(used==limit)throw RankTwoSourceError("RESOURCE_LIMIT","finite-cell topology step budget exhausted");++used;}
};
struct EdgeUse {std::size_t triangle,first,second;};
using Uses=std::map<Edge,std::vector<EdgeUse>>;
Uses uses(const std::vector<std::array<std::size_t,3>>& triangles,Budget& budget){
    Uses out;std::set<std::array<std::size_t,3>> seen;
    for(std::size_t i=0;i<triangles.size();++i){
        auto key=triangles[i];std::sort(key.begin(),key.end());
        if(key[0]==key[1]||key[1]==key[2]||!seen.insert(key).second)invalid("duplicate/degenerate combinatorial triangle");
        for(std::size_t j=0;j<3;++j){budget.tick();auto x=triangles[i][j],y=triangles[i][(j+1)%3];
            auto& edge=out[std::minmax(x,y)];edge.push_back({i,x,y});
            if(edge.size()>2)invalid("edge has more than two incident triangles");
            if(edge.size()==2&&(edge[0].first!=y||edge[0].second!=x))invalid("internal triangle edge has equal rather than opposite orientations");
        }
    }
    return out;
}
// Trace the actual union boundary; no native label or expected count drives
// the result. Convex-domain checks occur only after this extraction.
std::vector<detail::RankTwoComplexBoundary> traceBoundary(const Uses& edges,Budget& budget){
    std::map<std::size_t,std::size_t> next,previous;
    for(const auto& item:edges)if(item.second.size()==1){budget.tick();const auto& e=item.second[0];
        if(!next.emplace(e.first,e.second).second||!previous.emplace(e.second,e.first).second)
            invalid("actual union boundary branches or touches itself at a vertex");
    }
    std::set<std::size_t> remaining;for(const auto& p:next){remaining.insert(p.first);if(!previous.count(p.first))invalid("actual boundary is open");}
    std::vector<detail::RankTwoComplexBoundary> out;
    while(!remaining.empty()){
        const auto start=*remaining.begin();auto v=start;std::vector<std::size_t> ring;
        do{budget.tick();if(!remaining.erase(v))invalid("actual boundary revisits another cycle");ring.push_back(v);
            const auto f=next.find(v);if(f==next.end())invalid("actual boundary has missing successor");v=f->second;
        }while(v!=start);
        out.push_back({false,std::move(ring)});
    }
    return out;
}
}
namespace detail {
RankTwoComplexInspection inspectRankTwoTriangleComplex(
    const std::vector<std::array<std::size_t,3>>& triangles,
    const std::vector<RankTwoComplexBoundary>& boundaries,std::uint64_t maximumSteps){
    if(!maximumSteps)throw RankTwoSourceError("RESOURCE_LIMIT","positive topology step budget required");
    if(triangles.empty())invalid("empty actual triangle complex");
    Budget budget{0,maximumSteps};const auto edges=uses(triangles,budget);
    std::map<Edge,Edge> wanted;std::set<std::size_t> boundaryVertices;
    for(const auto& ring:boundaries){
        if(ring.points.size()<3)invalid("boundary cycle has fewer than three vertices");
        for(std::size_t i=0;i<ring.points.size();++i){budget.tick();const auto x=ring.points[i],y=ring.points[(i+1)%ring.points.size()];
            if(!boundaryVertices.insert(x).second)invalid("boundary cycles are not simple and disjoint");
            if(!wanted.emplace(std::minmax(x,y),Edge{x,y}).second)invalid("boundary edge has duplicate consumer");
        }
    }
    std::vector<std::vector<std::size_t>> adjacent(triangles.size());std::map<std::size_t,std::vector<Edge>> links;
    for(const auto& item:edges){budget.tick();const auto& list=item.second;const auto w=wanted.find(item.first);
        if(list.size()==1){if(w==wanted.end()||w->second!=Edge{list[0].first,list[0].second})invalid("actual boundary differs from oriented native cycles");}
        else{if(w!=wanted.end())invalid("native boundary is actually an internal edge");
            adjacent[list[0].triangle].push_back(list[1].triangle);adjacent[list[1].triangle].push_back(list[0].triangle);}
    }
    for(const auto& item:wanted)if(!edges.count(item.first))invalid("native boundary missing from actual triangles");
    for(const auto& t:triangles)for(std::size_t j=0;j<3;++j){budget.tick();links[t[j]].emplace_back(t[(j+1)%3],t[(j+2)%3]);}
    RankTwoComplexInspection out;
    for(const auto& vertex:links){
        const bool boundary=boundaryVertices.count(vertex.first);std::map<std::size_t,std::size_t> next,previous;
        std::set<std::size_t> neighbors;
        for(const auto& e:vertex.second){budget.tick();neighbors.insert(e.first);neighbors.insert(e.second);
            if(!next.emplace(e.first,e.second).second||!previous.emplace(e.second,e.first).second)
                invalid("vertex link branches or has duplicate oriented edges");}
        std::vector<std::size_t> starts,ends;
        for(const auto v:neighbors){budget.tick();if(!previous.count(v))starts.push_back(v);if(!next.count(v))ends.push_back(v);}
        if(boundary?(starts.size()!=1||ends.size()!=1):(!starts.empty()||!ends.empty()))
            invalid("vertex link is not a boundary path or interior cycle");
        const auto start=boundary?starts[0]:*neighbors.begin();auto v=start;std::set<std::size_t> seen;std::vector<std::size_t> order;
        while(seen.insert(v).second){budget.tick();order.push_back(v);const auto f=next.find(v);
            if(f==next.end())break;v=f->second;}
        if(seen!=neighbors||(boundary?(order.back()!=ends[0]):(v!=start)))invalid("vertex link is disconnected or pinched");
        out.links.push_back({vertex.first,boundary,std::move(order)});
    }
    std::vector<std::size_t> component(triangles.size(),triangles.size());
    for(std::size_t i=0;i<triangles.size();++i)if(component[i]==triangles.size()){
        const auto index=out.components.size();std::vector<std::size_t> pending{i},faces;component[i]=index;
        while(!pending.empty()){budget.tick();const auto f=pending.back();pending.pop_back();faces.push_back(f);
            for(const auto n:adjacent[f]){budget.tick();if(component[n]==triangles.size()){component[n]=index;pending.push_back(n);}}}
        std::sort(faces.begin(),faces.end());std::set<std::size_t> vertices;std::set<Edge> used;
        for(const auto f:faces){budget.tick();const auto& t=triangles[f];vertices.insert(t.begin(),t.end());
            for(std::size_t j=0;j<3;++j)used.insert(std::minmax(t[j],t[(j+1)%3]));}
        out.components.push_back({vertices.size(),used.size(),0,
            std::int64_t(vertices.size())-std::int64_t(used.size())+std::int64_t(faces.size()),std::move(faces),{}});
    }
    for(std::size_t i=0;i<boundaries.size();++i){budget.tick();const auto& ring=boundaries[i];
        const auto index=component[edges.at(std::minmax(ring.points[0],ring.points[1]))[0].triangle];
        for(std::size_t j=0;j<ring.points.size();++j){budget.tick();
            if(component[edges.at(std::minmax(ring.points[j],ring.points[(j+1)%ring.points.size()]))[0].triangle]!=index)
                invalid("one boundary cycle assigned to multiple components");}
        out.components[index].boundaries.push_back(i);out.components[index].holes+=ring.inner;
    }
    for(const auto& c:out.components){
        if(c.boundaries.size()!=c.holes+1||c.eulerCharacteristic!=1-std::int64_t(c.holes))
            invalid("component lacks exactly one outer boundary or has incorrect genus/hole assignment");
    }
    std::ostringstream record;record.imbue(std::locale::classic());record<<"rank-two-orientable-triangle-complex-v1\n";
    for(std::size_t i=0;i<boundaries.size();++i){record<<"BOUNDARY "<<i<<' '<<(boundaries[i].inner?"INNER":"OUTER");for(const auto v:boundaries[i].points)record<<' '<<v;record<<'\n';}
    for(const auto& l:out.links){record<<"LINK "<<l.vertex<<' '<<(l.boundary?"BOUNDARY":"INTERIOR");for(const auto v:l.neighbors)record<<' '<<v;record<<'\n';}
    for(std::size_t i=0;i<out.components.size();++i){const auto& c=out.components[i];
        record<<"COMPONENT "<<i<<' '<<c.vertexCount<<' '<<c.edgeCount<<' '<<c.triangles.size()<<' '<<c.holes<<' '<<c.eulerCharacteristic<<'\n';
        record<<"COMPONENT_TRIANGLES "<<i;for(const auto f:c.triangles)record<<' '<<f;record<<'\n';
        record<<"COMPONENT_BOUNDARIES "<<i;for(const auto b:c.boundaries)record<<' '<<b;record<<'\n';}
    out.workSteps=budget.used;out.canonicalRecord=record.str();return out;
}
}
struct RankTwoFiniteCellTopology::Impl {
    RankTwoCoveredFiniteMesh mesh;std::size_t regions=0,holes=0;std::string record,hash;
    explicit Impl(RankTwoCoveredFiniteMesh m):mesh(std::move(m)){}
};
RankTwoFiniteCellTopology::RankTwoFiniteCellTopology(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoCoveredFiniteMesh& RankTwoFiniteCellTopology::coveredMesh()const{return impl_->mesh;}
std::size_t RankTwoFiniteCellTopology::regionCount()const{return impl_->regions;}
std::size_t RankTwoFiniteCellTopology::nativeHoleCount()const{return impl_->holes;}
const std::string& RankTwoFiniteCellTopology::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoFiniteCellTopology::scientificHash()const{return impl_->hash;}
RankTwoFiniteCellTopology certifyRankTwoFiniteCellTopology(const RankTwoArrangement& a,
    const RankTwoCoveredFiniteMesh& mesh,const RankTwoFiniteCellTopologySettings& settings){
    const auto& p=mesh.coverage().partition();const auto& r=p.registry();const auto& t=mesh.triangulation();
    if(a.scientificHash()!=p.arrangementIdentity()||a.scientificHash()!=mesh.coverage().arrangementIdentity()
        ||p.scientificHash()!=t.partitionIdentity()||r.scientificHash()!=t.registryIdentity())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","cell topology requires identical native and finite authorities");
    if(!settings.maximumTopologySteps)throw RankTwoSourceError("RESOURCE_LIMIT","positive cell topology budget required");
    auto out=std::make_shared<RankTwoFiniteCellTopology::Impl>(mesh);std::uint64_t remaining=settings.maximumTopologySteps;
    std::ostringstream record;record.imbue(std::locale::classic());record<<"rank-two-finite-cell-topology-v1\nAUTHORITY "
        <<a.scientificHash()<<' '<<mesh.scientificHash()<<' '<<p.scientificHash()<<' '<<r.scientificHash()<<' '<<t.scientificHash()<<'\n';
    std::vector<std::array<std::size_t,3>> whole;std::map<std::size_t,std::vector<std::array<std::size_t,3>>> cells;
    for(const auto& triangle:t.triangles()){whole.push_back(triangle.points);cells[triangle.registryRegionIndex].push_back(triangle.points);}
    if(cells.size()!=p.regionIndices().size())invalid("actual finite cell set differs from native domain");
    for(const auto region:p.regionIndices()){
        const auto& native=r.regions().at(region);std::vector<detail::RankTwoComplexBoundary> boundaries;
        std::set<std::size_t> loopVertices;
        for(bool inner:{false,true}){const auto& loops=inner?native.innerLoops:native.outerLoops;
            for(std::size_t i=0;i<loops.size();++i){const auto vertices=r.loopVertices(native.identity,inner,i);
                for(const auto v:vertices)if(!loopVertices.insert(v).second)
                    throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","native cell closure has touching/non-Jordan boundaries; no finite manifold grant");
                boundaries.push_back({inner,vertices});}}
        const auto inspected=detail::inspectRankTwoTriangleComplex(cells.at(region),boundaries,remaining);
        remaining-=inspected.workSteps;
        // Isolated native points remain actual named interior vertices, not
        // artificial punctures or points silently dropped from the carrier.
        for(const auto point:native.isolatedPoints){const auto f=std::find_if(inspected.links.begin(),inspected.links.end(),
            [point](const auto& link){return link.vertex==point;});
            if(f==inspected.links.end()||f->boundary)invalid("isolated native event not retained as actual interior vertex");}
        record<<"CELL_BEGIN "<<region<<' '<<native.identity<<' '<<int(native.kind)<<'\n'<<inspected.canonicalRecord<<"CELL_END\n";
        ++out->regions;out->holes+=native.innerLoops.size();
    }
    Budget tracing{0,remaining};const auto unionEdges=uses(whole,tracing);const auto unionBoundary=traceBoundary(unionEdges,tracing);
    remaining-=tracing.used;
    // D>0 intersected with the actual linear halfplanes is convex. All cells
    // including Empty must assemble into one unpunctured finite disk. This
    // does not assemble the union of distinct normal-owner cones.
    if(unionBoundary.size()!=1)invalid("complete convex finite domain is disconnected or has a spurious hole");
    const auto domain=detail::inspectRankTwoTriangleComplex(whole,unionBoundary,remaining);
    if(domain.components.size()!=1||domain.components[0].holes||domain.components[0].eulerCharacteristic!=1)
        invalid("complete finite domain is not an orientable disk");
    record<<"DOMAIN_BEGIN\n"<<domain.canonicalRecord<<"DOMAIN_END\nCOVERED_MESH_BEGIN\n"<<mesh.canonicalRecord()
        <<"COVERED_MESH_END\nfinite_manifold_cell_correspondence PASSED\n"
        <<"global_isotopy NOT_GRANTED\nall_normal_owners NOT_GRANTED\nprofile_precision NOT_GRANTED\n"
        <<"strict_thickness NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-finite-cell-topology-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    return RankTwoFiniteCellTopology(std::move(out));
}
}
