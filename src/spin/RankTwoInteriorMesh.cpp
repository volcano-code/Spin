#include "detail/RankTwoInteriorMesh.h"
#include "detail/RankTwoMovingMeshEdge.h"
#include "detail/RankTwoObservedCartesian.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <deque>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using K=detail::RankTwoObservedCartesian;using Q=CGAL::Gmpq;
using Pair=std::pair<std::size_t,std::size_t>;
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
Pair edge(std::size_t a,std::size_t b){return std::minmax(a,b);}
void boundary(std::map<Pair,int>& b,const RankTwoFiniteTriangle& t){
    for(int j=0;j<3;++j){auto a=t.points[j],v=t.points[(j+1)%3];b[edge(a,v)]+=a<v?1:-1;}}
void trim(std::map<Pair,int>& b){for(auto i=b.begin();i!=b.end();)if(!i->second)i=b.erase(i);else ++i;}
void coord(std::ostream& s,const detail::CertifiedAlgebraicReal& value){const auto& c=value.certificate();
    s<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& k:c.coefficients)s<<k<<',';s<<'\n';}
}
struct RankTwoInteriorMesh::Impl {
    RankTwoCoveredFiniteMesh base;
    std::vector<RankTwoInteriorMeshPoint> points;
    std::vector<RankTwoFiniteTriangle> triangles;
    std::vector<RankTwoRegionTriangulation> summaries;
    std::string record,hash;
    explicit Impl(RankTwoCoveredFiniteMesh b):base(std::move(b)){}
};
RankTwoInteriorMesh::RankTwoInteriorMesh(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoCoveredFiniteMesh& RankTwoInteriorMesh::baseMesh()const{return impl_->base;}
const std::vector<RankTwoInteriorMeshPoint>& RankTwoInteriorMesh::points()const{return impl_->points;}
const std::vector<RankTwoFiniteTriangle>& RankTwoInteriorMesh::triangles()const{return impl_->triangles;}
const std::vector<RankTwoRegionTriangulation>& RankTwoInteriorMesh::regionSummaries()const{return impl_->summaries;}
const std::string& RankTwoInteriorMesh::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoInteriorMesh::scientificHash()const{return impl_->hash;}

RankTwoInteriorMesh insertRankTwoMovingInteriorNodes(const RankTwoArrangement& a,
    const RankTwoCoveredFiniteMesh& base,const std::vector<RankTwoMovingInteriorPoint>& nodes,
    const RankTwoInteriorMeshSettings& settings){
    detail::requireExactAlgebraicEnvironment();const auto& p=base.coverage().partition();const auto& r=p.registry();
    if(a.scientificHash()!=base.coverage().arrangementIdentity() || a.scientificHash()!=p.arrangementIdentity())
        invalid("moving insertion requires the SAME actual covered domain");
    if(!settings.maximumPredicates)throw RankTwoSourceError("RESOURCE_LIMIT","positive incidence proof budget required");
    std::uint64_t predicates=0;
    const auto consume=[&]{if(predicates==settings.maximumPredicates)
        throw RankTwoSourceError("RESOURCE_LIMIT","moving mesh incidence budget exhausted");++predicates;};
    std::set<std::string> identities;
    std::map<std::string,std::size_t> regions;
    for(const auto i:p.regionIndices())regions.emplace(r.regions().at(i).identity,i);
    std::vector<const RankTwoMovingInteriorPoint*> ordered;
    for(const auto& n:nodes){
        if(n.arrangementIdentity()!=a.scientificHash() || !regions.count(n.definition().regionIdentity)
            || !identities.insert(n.scientificHash()).second)invalid("foreign/duplicate moving node or excluded native owner");
        const auto kind=r.regions().at(regions.at(n.definition().regionIdentity)).kind;
        if(kind!=n.regionKind() || (kind!=RankTwoRegionKind::Full && kind!=RankTwoRegionKind::Cut))
            invalid("moving node active native label mismatch");
        ordered.push_back(&n);
    }
    std::sort(ordered.begin(),ordered.end(),[](const auto* x,const auto* y){return x->scientificHash()<y->scientificHash();});
    auto out=std::make_shared<RankTwoInteriorMesh::Impl>(base);
    for(const auto& v:r.points())out->points.push_back({v.identity,v.coordinates});
    out->triangles=base.triangulation().triangles();
    std::set<Pair> constraints;
    for(const auto i:p.edgeIndices()){const auto& e=r.edges().at(i);constraints.insert(edge(e.points[0],e.points[1]));}
    const auto point=[&](std::size_t i){const auto& xy=out->points.at(i).coordinates;return K::Point_2(xy[0].expression(),xy[1].expression());};
    const auto orientation=[&](std::size_t u,std::size_t v,std::size_t w){consume();
        return K().orientation_2_object()(point(u),point(v),point(w));};
    std::ostringstream steps;steps.imbue(std::locale::classic());
    for(const auto* node:ordered){
        const auto index=out->points.size(),owner=regions.at(node->definition().regionIdentity);
        out->points.push_back({node->scientificHash(),node->coordinates()});
        detail::rankTwoDiagnostic(("interior-mesh:locate-begin "+std::to_string(index)).c_str());
        std::vector<std::size_t> incident;
        std::map<std::size_t,std::array<CGAL::Orientation,3>> signs;
        // Root-isolator rectangles are conservative exclusion ONLY. All
        // admitted candidate faces are tested by exact orientation, never by
        // interval midpoints or numerical distance tolerances.
        for(std::size_t t=0;t<out->triangles.size();++t){consume();const auto& face=out->triangles[t];bool excluded=false;
            for(int d=0;d<2;++d){Q lo(out->points[face.points[0]].coordinates[d].certificate().lower);
                Q hi(out->points[face.points[0]].coordinates[d].certificate().upper);
                for(int j=1;j<3;++j){const auto& c=out->points[face.points[j]].coordinates[d].certificate();
                    lo=std::min(lo,Q(c.lower));hi=std::max(hi,Q(c.upper));}
                const auto& c=node->coordinates()[d].certificate();
                if(Q(c.upper)<lo || Q(c.lower)>hi)excluded=true;}
            if(excluded)continue;
            std::array<CGAL::Orientation,3> s;bool contains=true;
            for(int j=0;j<3;++j){s[j]=orientation(face.points[j],face.points[(j+1)%3],index);
                if(s[j]==CGAL::RIGHT_TURN){contains=false;break;}}
            if(contains){incident.push_back(t);signs.emplace(t,s);}
        }
        if(incident.empty() || incident.size()>2)topology("moving native point has no unique face/internal-edge incidence");
        for(const auto t:incident)if(out->triangles[t].registryRegionIndex!=owner)
            topology("analytic moving node is not inside its actual finite native polygon");
        std::vector<RankTwoFiniteTriangle> children;
        if(incident.size()==1){const auto t=incident[0];const auto& f=out->triangles[t];
            for(int j=0;j<3;++j){if(signs.at(t)[j]!=CGAL::LEFT_TURN)
                    topology("single-face moving insertion touches a native constraint or existing vertex");
                children.push_back({{f.points[j],f.points[(j+1)%3],index},owner});}
        }else{
            std::optional<Pair> shared;
            for(const auto t:incident){const auto& f=out->triangles[t];int zeros=0;
                for(int j=0;j<3;++j)if(signs.at(t)[j]==CGAL::COLLINEAR){++zeros;const auto e=edge(f.points[j],f.points[(j+1)%3]);
                    if(shared && *shared!=e)topology("inconsistent shared-edge incidence");shared=e;}
                if(zeros!=1)topology("moving insertion coincides with a vertex or ambiguous edge");}
            if(!shared || constraints.count(*shared))topology("moving internal node cannot replace a native boundary/interface constraint");
            for(const auto t:incident){const auto& f=out->triangles[t];
                for(int j=0;j<3;++j)if(edge(f.points[j],f.points[(j+1)%3])==*shared){
                    children.push_back({{f.points[j],index,f.points[(j+2)%3]},owner});
                    children.push_back({{index,f.points[(j+1)%3],f.points[(j+2)%3]},owner});}}
        }
        std::map<Pair,int> before,after;
        for(const auto t:incident)boundary(before,out->triangles[t]);
        for(const auto& child:children){
            if(orientation(child.points[0],child.points[1],child.points[2])!=CGAL::LEFT_TURN)
                topology("zero/reversed actual stellar child");boundary(after,child);}
        trim(before);trim(after);if(before!=after)topology("stellar replacement does not cover the same directed incident boundary");
        steps<<"INSERT "<<index<<' '<<owner<<' '<<node->scientificHash()<<' '<<incident.size()<<' '<<children.size()<<'\n';
        for(const auto t:incident){const auto& f=out->triangles[t];steps<<"REPLACED "<<owner<<' '<<f.points[0]<<' '<<f.points[1]<<' '<<f.points[2]<<'\n';}
        for(const auto& f:children)steps<<"CHILD "<<owner<<' '<<f.points[0]<<' '<<f.points[1]<<' '<<f.points[2]<<'\n';
        std::vector<RankTwoFiniteTriangle> next;
        for(std::size_t t=0;t<out->triangles.size();++t)if(std::find(incident.begin(),incident.end(),t)==incident.end())next.push_back(out->triangles[t]);
        next.insert(next.end(),children.begin(),children.end());out->triangles=std::move(next);
        detail::rankTwoDiagnostic(("interior-mesh:locate-end "+std::to_string(index)).c_str());
    }
    for(auto& t:out->triangles)std::rotate(t.points.begin(),std::min_element(t.points.begin(),t.points.end()),t.points.end());
    std::sort(out->triangles.begin(),out->triangles.end(),[](const auto& x,const auto& y){return x.points!=y.points?x.points<y.points:x.registryRegionIndex<y.registryRegionIndex;});
    std::set<std::size_t> used;
    for(const auto region:p.regionIndices()){
        std::vector<std::size_t> own;std::set<std::size_t> vertices;std::map<Pair,std::vector<std::size_t>> edges;
        std::map<Pair,int> actual,wanted;
        for(std::size_t t=0;t<out->triangles.size();++t)if(out->triangles[t].registryRegionIndex==region){consume();
            own.push_back(t);const auto& f=out->triangles[t];vertices.insert(f.points.begin(),f.points.end());used.insert(f.points.begin(),f.points.end());boundary(actual,f);
            for(int j=0;j<3;++j)edges[edge(f.points[j],f.points[(j+1)%3])].push_back(t);}
        const auto& cell=r.regions().at(region);
        for(bool inner:{false,true})for(const auto& loop:inner?cell.innerLoops:cell.outerLoops)for(const auto& ref:loop.edges){
            const auto& e=r.edges().at(ref.edgeIndex);auto u=e.points[0],v=e.points[1];if(ref.reversed)std::swap(u,v);wanted[edge(u,v)]+=u<v?1:-1;}
        trim(actual);trim(wanted);if(own.empty() || actual!=wanted)topology("moving subdivision changed native covered rings");
        std::map<std::size_t,std::vector<std::size_t>> adjacent;
        for(const auto& e:edges){consume();if(e.second.size()>2)topology("nonmanifold moving triangle edge");
            if(e.second.size()==2){adjacent[e.second[0]].push_back(e.second[1]);adjacent[e.second[1]].push_back(e.second[0]);}}
        std::set<std::size_t> visited;std::size_t components=0;
        for(const auto t:own)if(visited.insert(t).second){++components;std::deque<std::size_t> queue{t};
            while(!queue.empty()){consume();const auto first=queue.front();queue.pop_front();for(const auto next:adjacent[first])
                if(visited.insert(next).second)queue.push_back(next);}}
        const auto chi=static_cast<std::int64_t>(vertices.size())-static_cast<std::int64_t>(edges.size())+static_cast<std::int64_t>(own.size());
        const auto holes=static_cast<std::int64_t>(components)-chi;
        if(holes<0 || components!=cell.outerLoops.size() || std::size_t(holes)!=cell.innerLoops.size())topology("moving mesh connectivity/Euler holes changed");
        out->summaries.push_back({region,own.size(),vertices.size(),edges.size(),components,std::size_t(holes),chi});
    }
    std::set<std::size_t> required(p.pointIndices().begin(),p.pointIndices().end());
    for(std::size_t i=r.points().size();i<out->points.size();++i)required.insert(i);
    if(used!=required)topology("actual moving mesh omitted/shared an undeclared point");
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-interior-finite-mesh-v1\nAUTHORITY "<<a.scientificHash()<<' '<<base.scientificHash()<<'\n'
        <<"POLICY EXACT_STELLAR_SUBDIVISION SAME_NATIVE_RINGS\nBUDGET "<<settings.maximumPredicates<<" PREDICATES "<<predicates<<'\n';
    for(std::size_t i=0;i<out->points.size();++i){const auto& v=out->points[i];
        record<<"POINT "<<i<<' '<<(i<r.points().size()?"NATIVE":"INTERIOR")<<' '<<v.identity<<'\n';coord(record,v.coordinates[0]);coord(record,v.coordinates[1]);}
    for(const auto* n:ordered)record<<"NODE_BEGIN "<<n->scientificHash()<<'\n'<<n->canonicalRecord()<<"NODE_END\n";
    record<<steps.str();
    for(const auto& t:out->triangles)record<<"TRIANGLE "<<t.registryRegionIndex<<' '<<t.points[0]<<' '<<t.points[1]<<' '<<t.points[2]<<'\n';
    for(const auto& s:out->summaries)record<<"REGION_SUMMARY "<<s.registryRegionIndex<<' '<<s.triangleCount<<' '<<s.vertexCount<<' '<<s.edgeCount<<' '<<s.connectedComponents<<' '<<s.holes<<' '<<s.eulerCharacteristic<<'\n';
    record<<"BASE_MESH_BEGIN\n"<<base.canonicalRecord()<<"BASE_MESH_END\ninterior_node_insertion PASSED\n"
        <<"strict_thickness NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-interior-finite-mesh-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoInteriorMesh(std::move(out));
}
RankTwoInteriorMesh insertRankTwoMovingEdgeNodes(const RankTwoArrangement& a,
    const RankTwoInteriorMesh& parent,const std::vector<RankTwoMovingMeshEdgePoint>& nodes,
    const RankTwoInteriorMeshSettings& settings){
    detail::requireExactAlgebraicEnvironment();const auto& p=parent.baseMesh().coverage().partition();const auto& r=p.registry();
    if(a.scientificHash()!=p.arrangementIdentity())invalid("moving edge insertion needs same parent authority");
    if(!settings.maximumPredicates)throw RankTwoSourceError("RESOURCE_LIMIT","positive edge subdivision budget required");
    std::uint64_t predicates=0;const auto consume=[&]{if(predicates==settings.maximumPredicates)
        throw RankTwoSourceError("RESOURCE_LIMIT","edge subdivision budget exhausted");++predicates;};
    std::set<std::string> identities;std::set<Pair> requested;
    std::vector<const RankTwoMovingMeshEdgePoint*> ordered;
    for(const auto& n:nodes){const auto& d=n.definition();
        if(n.arrangementIdentity()!=a.scientificHash() || n.parentMeshIdentity()!=parent.scientificHash()
            || !identities.insert(n.scientificHash()).second || !requested.insert(edge(d.pointIndices[0],d.pointIndices[1])).second)
            invalid("duplicate/foreign moving edge or repeated subdivision of one edge");
        ordered.push_back(&n);}
    std::sort(ordered.begin(),ordered.end(),[](const auto* x,const auto* y){return x->scientificHash()<y->scientificHash();});
    auto out=std::make_shared<RankTwoInteriorMesh::Impl>(*parent.impl_);out->summaries.clear();
    const auto point=[&](std::size_t i){const auto& xy=out->points.at(i).coordinates;return K::Point_2(xy[0].expression(),xy[1].expression());};
    const auto orientation=[&](const RankTwoFiniteTriangle& t){consume();return K().orientation_2_object()(point(t.points[0]),point(t.points[1]),point(t.points[2]));};
    // Monotonically allocated IDs reproduce the old survivor/append order,
    // without reindexing every unrelated face after each stellar split.
    // This is CURRENT exact combinatorial incidence, not stale parent face
    // indices, rounded spatial lookup or a caller-supplied successful subset.
    std::map<std::size_t,RankTwoFiniteTriangle> current;
    std::map<Pair,std::set<std::size_t>> incidence;
    std::size_t nextFace=0;
    const auto addFace=[&](const RankTwoFiniteTriangle& face){
        const auto id=nextFace++;current.emplace(id,face);
        for(int j=0;j<3;++j){consume();
            if(!incidence[edge(face.points[j],face.points[(j+1)%3])].insert(id).second)
                topology("repeated edge in indexed moving face");}
    };
    const auto removeFace=[&](std::size_t id){
        const auto& face=current.at(id);
        for(int j=0;j<3;++j){consume();const auto key=edge(face.points[j],face.points[(j+1)%3]);
            auto found=incidence.find(key);
            if(found==incidence.end() || found->second.erase(id)!=1)
                topology("stale indexed moving edge incidence");
            if(found->second.empty())incidence.erase(found);}
        current.erase(id);
    };
    for(const auto& face:out->triangles)addFace(face);
    std::ostringstream steps;steps.imbue(std::locale::classic());
    for(const auto* n:ordered){const auto& d=n->definition();const auto target=edge(d.pointIndices[0],d.pointIndices[1]);
        const auto index=out->points.size(),owner=n->registryRegionIndex();
        std::vector<std::size_t> incident;std::vector<RankTwoFiniteTriangle> children;
        consume();const auto found=incidence.find(target);
        if(found==incidence.end() || found->second.size()!=2)
            topology("moving internal edge lost its two actual incident faces");
        // Copy before any mutation; ordered IDs give exactly the former
        // ascending current-face traversal for adjacent-edge repairs too.
        incident.assign(found->second.begin(),found->second.end());
        for(const auto i:incident){consume();const auto& f=current.at(i);bool matched=false;
            for(int j=0;j<3;++j)if(edge(f.points[j],f.points[(j+1)%3])==target){
                if(matched)topology("ambiguous indexed moving edge");matched=true;
                if(f.registryRegionIndex!=owner)topology("moving edge subdivision crossed native owner");
                children.push_back({{f.points[j],index,f.points[(j+2)%3]},owner});
                children.push_back({{index,f.points[(j+1)%3],f.points[(j+2)%3]},owner});}
            if(!matched)topology("indexed moving edge is absent from actual face");}
        out->points.push_back({n->scientificHash(),n->coordinates()});
        std::map<Pair,int> before,after;
        for(const auto t:incident)boundary(before,current.at(t));
        for(const auto& child:children){if(orientation(child)!=CGAL::LEFT_TURN)topology("moving edge has zero/reversed child");boundary(after,child);}
        trim(before);trim(after);if(before!=after)topology("moving edge changed directed incident boundary");
        steps<<"INSERT_EDGE "<<index<<' '<<owner<<' '<<n->scientificHash()<<' '<<target.first<<' '<<target.second<<" 2 4\n";
        for(const auto t:incident){const auto& f=current.at(t);steps<<"REPLACED "<<owner<<' '<<f.points[0]<<' '<<f.points[1]<<' '<<f.points[2]<<'\n';}
        for(const auto& f:children)steps<<"CHILD "<<owner<<' '<<f.points[0]<<' '<<f.points[1]<<' '<<f.points[2]<<'\n';
        for(const auto t:incident)removeFace(t);
        for(const auto& child:children)addFace(child);
    }
    out->triangles.clear();out->triangles.reserve(current.size());
    for(const auto& item:current)out->triangles.push_back(item.second);
    for(auto& t:out->triangles)std::rotate(t.points.begin(),std::min_element(t.points.begin(),t.points.end()),t.points.end());
    std::sort(out->triangles.begin(),out->triangles.end(),[](const auto& x,const auto& y){return x.points!=y.points?x.points<y.points:x.registryRegionIndex<y.registryRegionIndex;});
    std::set<std::size_t> used;
    for(const auto region:p.regionIndices()){
        std::vector<std::size_t> own;std::set<std::size_t> vertices;std::map<Pair,std::vector<std::size_t>> edges;
        std::map<Pair,int> actual,wanted;
        for(std::size_t t=0;t<out->triangles.size();++t)if(out->triangles[t].registryRegionIndex==region){consume();
            const auto& f=out->triangles[t];own.push_back(t);vertices.insert(f.points.begin(),f.points.end());used.insert(f.points.begin(),f.points.end());boundary(actual,f);
            for(int j=0;j<3;++j)edges[edge(f.points[j],f.points[(j+1)%3])].push_back(t);}
        const auto& cell=r.regions().at(region);
        for(bool inner:{false,true})for(const auto& loop:inner?cell.innerLoops:cell.outerLoops)for(const auto& ref:loop.edges){
            const auto& e=r.edges().at(ref.edgeIndex);auto u=e.points[0],v=e.points[1];if(ref.reversed)std::swap(u,v);wanted[edge(u,v)]+=u<v?1:-1;}
        trim(actual);trim(wanted);if(own.empty() || actual!=wanted)topology("edge subdivision changed original native rings");
        std::map<std::size_t,std::vector<std::size_t>> adjacent;
        for(const auto& e:edges){consume();if(e.second.size()>2)topology("nonmanifold edge subdivision");
            if(e.second.size()==2){adjacent[e.second[0]].push_back(e.second[1]);adjacent[e.second[1]].push_back(e.second[0]);}}
        std::set<std::size_t> visited;std::size_t components=0;
        for(const auto t:own)if(visited.insert(t).second){++components;std::deque<std::size_t> queue{t};
            while(!queue.empty()){consume();const auto first=queue.front();queue.pop_front();for(const auto next:adjacent[first])
                if(visited.insert(next).second)queue.push_back(next);}}
        const auto chi=static_cast<std::int64_t>(vertices.size())-static_cast<std::int64_t>(edges.size())+static_cast<std::int64_t>(own.size());
        const auto holes=static_cast<std::int64_t>(components)-chi;
        if(holes<0 || components!=cell.outerLoops.size() || std::size_t(holes)!=cell.innerLoops.size())topology("edge subdivision changed native connectivity/holes");
        out->summaries.push_back({region,own.size(),vertices.size(),edges.size(),components,std::size_t(holes),chi});
    }
    std::set<std::size_t> required(p.pointIndices().begin(),p.pointIndices().end());
    for(std::size_t i=r.points().size();i<out->points.size();++i)required.insert(i);
    if(used!=required)topology("moving edge mesh omitted an original or newly inserted point");
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-edge-refined-finite-mesh-v1\nAUTHORITY "<<a.scientificHash()<<' '<<parent.scientificHash()<<'\n'
        <<"POLICY EXACT_TWO_SIDED_EDGE_SUBDIVISION SAME_NATIVE_RINGS\nBUDGET "<<settings.maximumPredicates<<" PREDICATES "<<predicates<<'\n';
    for(std::size_t i=0;i<out->points.size();++i){const auto& v=out->points[i];record<<"POINT "<<i<<' '
        <<(i<parent.points().size()?"PARENT":"EDGE")<<' '<<v.identity<<'\n';coord(record,v.coordinates[0]);coord(record,v.coordinates[1]);}
    for(const auto* n:ordered)record<<"NODE_BEGIN "<<n->scientificHash()<<'\n'<<n->canonicalRecord()<<"NODE_END\n";
    record<<steps.str();for(const auto& t:out->triangles)record<<"TRIANGLE "<<t.registryRegionIndex<<' '<<t.points[0]<<' '<<t.points[1]<<' '<<t.points[2]<<'\n';
    for(const auto& s:out->summaries)record<<"REGION_SUMMARY "<<s.registryRegionIndex<<' '<<s.triangleCount<<' '<<s.vertexCount<<' '<<s.edgeCount<<' '<<s.connectedComponents<<' '<<s.holes<<' '<<s.eulerCharacteristic<<'\n';
    record<<"PARENT_MESH_BEGIN\n"<<parent.canonicalRecord()<<"PARENT_MESH_END\nedge_node_insertion PASSED\n"
        <<"strict_thickness NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\n"
        <<"source_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-edge-refined-finite-mesh-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoInteriorMesh(std::move(out));
}
} // namespace spin
