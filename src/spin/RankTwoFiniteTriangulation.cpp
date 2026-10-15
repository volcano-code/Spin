#include "detail/RankTwoFiniteTriangulation.h"
#include "detail/RankTwoDiagnostics.h"
#include "detail/RankTwoObservedCartesian.h"
#include "verification/Sha256.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <algorithm>
#include <deque>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using K=detail::RankTwoObservedCartesian;
constexpr auto outside=std::numeric_limits<std::size_t>::max(),unassigned=outside-1;
struct FaceInfo {std::size_t label=unassigned;};
using Vb=CGAL::Triangulation_vertex_base_with_info_2<std::size_t,K>;
using Cb=CGAL::Constrained_triangulation_face_base_2<K>;
using Fb=CGAL::Triangulation_face_base_with_info_2<FaceInfo,K,Cb>;
using Tds=CGAL::Triangulation_data_structure_2<Vb,Fb>;
using CT=CGAL::Constrained_triangulation_2<K,Tds,CGAL::No_intersection_tag>;
using Pair=std::pair<std::size_t,std::size_t>;
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
[[noreturn]]void resource(const char* why){throw RankTwoSourceError("RESOURCE_LIMIT",why);}
std::string sha(const std::string& bytes){verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());return h.finish();}
Pair ordered(std::size_t a,std::size_t b){return a<b?Pair{a,b}:Pair{b,a};}
void boundary(std::map<Pair,std::int64_t>& m,std::size_t a,std::size_t b){++m[{a,b}];--m[{b,a}];}
void trim(std::map<Pair,std::int64_t>& m){for(auto p=m.begin();p!=m.end();)if(!p->second)p=m.erase(p);else ++p;}
}
struct RankTwoFiniteTriangulation::Impl {
    std::string partition,registry,record,hash;
    std::vector<RankTwoFiniteTriangle> triangles;
    std::vector<RankTwoRegionTriangulation> summaries;
};
RankTwoFiniteTriangulation::RankTwoFiniteTriangulation(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoFiniteTriangulation::partitionIdentity()const{return impl_->partition;}
const std::string& RankTwoFiniteTriangulation::registryIdentity()const{return impl_->registry;}
const std::vector<RankTwoFiniteTriangle>& RankTwoFiniteTriangulation::triangles()const{return impl_->triangles;}
const std::vector<RankTwoRegionTriangulation>& RankTwoFiniteTriangulation::regionSummaries()const{return impl_->summaries;}
const std::string& RankTwoFiniteTriangulation::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoFiniteTriangulation::scientificHash()const{return impl_->hash;}

RankTwoFiniteTriangulation triangulateRankTwoFinitePartition(const RankTwoArrangement& a,
    const RankTwoFinitePartition& p,const RankTwoChordIncidence& incidence,const RankTwoRegistryRotation& rotation,
    const RankTwoFiniteTriangulationSettings& settings){
    detail::requireExactAlgebraicEnvironment();const auto& r=p.registry();
    if(a.scientificHash()!=p.arrangementIdentity() || incidence.arrangementIdentity()!=a.scientificHash()
        || incidence.registryIdentity()!=r.scientificHash() || rotation.arrangementIdentity()!=a.scientificHash()
        || rotation.registryIdentity()!=r.scientificHash() || rotation.incidenceIdentity()!=incidence.scientificHash()
        ||(!incidence.partitionIdentity().empty() && incidence.partitionIdentity()!=p.scientificHash())
        ||rotation.partitionIdentity()!=incidence.partitionIdentity())
        invalid("finite triangulation requires the SAME partition, incidence and native rotation authorities");
    if(!settings.maximumTriangles || !settings.maximumAdjacencySteps)resource("positive finite mesh budgets required");
    std::uint64_t steps=0;const auto consume=[&]{if(steps==settings.maximumAdjacencySteps)
        resource("finite triangulation adjacency budget exhausted");++steps;};
    auto out=std::make_shared<RankTwoFiniteTriangulation::Impl>();out->partition=p.scientificHash();out->registry=r.scientificHash();
    detail::RankTwoPredicateObservation observation;
    const bool observing=detail::rankTwoDiagnosticTarget.sink!=nullptr;
    CT ct(K(observing?&observation:nullptr));std::map<std::size_t,CT::Vertex_handle> handles;std::set<CT::Vertex_handle> inserted;
    auto points=p.pointIndices();std::sort(points.begin(),points.end(),[&](auto x,auto y){return r.points()[x].identity<r.points()[y].identity;});
    for(const auto index:points){consume();const auto& point=r.points().at(index);
        if(detail::rankTwoDiagnosticTarget.sink){const auto message="finite-cdt:insert-begin "+std::to_string(index)+" "+point.identity;
            detail::rankTwoDiagnostic(message.c_str());}
        const K::Point_2 actualPoint(point.coordinates[0].expression(),point.coordinates[1].expression());
        if(observing)observation.bind(actualPoint,index);
        const auto v=ct.insert(actualPoint);
        if(!inserted.insert(v).second)topology("distinct source point identities coincide at finite CT insertion");
        v->info()=index;handles.emplace(index,v);
        if(detail::rankTwoDiagnosticTarget.sink){const auto message="finite-cdt:insert-end "+std::to_string(index);
            detail::rankTwoDiagnostic(message.c_str());}
    }
    std::set<std::size_t> included(p.regionIndices().begin(),p.regionIndices().end());
    if(included.empty()){
        if(!points.empty() || !p.edgeIndices().empty())topology("empty finite domain has anonymous boundary features");
    }else if(ct.dimension()!=2)topology("included finite domain has no two-dimensional shared point span");
    auto constraints=p.edgeIndices();std::sort(constraints.begin(),constraints.end(),[&](auto x,auto y){return r.edges()[x].identity<r.edges()[y].identity;});
    for(const auto e:constraints){consume();const auto& edge=r.edges().at(e);
        if(detail::rankTwoDiagnosticTarget.sink){const auto message="finite-cdt:constraint-begin "+std::to_string(e);
            detail::rankTwoDiagnostic(message.c_str());}
        observation.constraint=e;
        ct.insert_constraint(handles.at(edge.points[0]),handles.at(edge.points[1]));
        observation.constraint.reset();
        if(ct.number_of_vertices()!=points.size())topology("No_intersection_tag CT created an undeclared source point");
        if(detail::rankTwoDiagnosticTarget.sink){const auto message="finite-cdt:constraint-end "+std::to_string(e);
            detail::rankTwoDiagnostic(message.c_str());}
    }
    std::map<std::string,std::size_t> regionIndex;
    for(std::size_t i=0;i<r.regions().size();++i)regionIndex.emplace(r.regions()[i].identity,i);
    std::map<std::string,const RankTwoRegistryArc*> arcs;for(const auto& arc:r.arcs())arcs.emplace(arc.identity,&arc);
    std::deque<CT::Face_handle> queue;
    const auto assign=[&](CT::Face_handle f,std::size_t label){consume();
        if(f->info().label==unassigned){f->info().label=label;queue.push_back(f);}
        else if(f->info().label!=label)topology("finite constraint sides/nesting assign conflicting original native cells");
    };
    if(!included.empty()){
        for(auto f=ct.all_faces_begin();f!=ct.all_faces_end();++f)if(ct.is_infinite(f))assign(f,outside);
        for(const auto e:constraints){consume();const auto& edge=r.edges().at(e);CT::Face_handle f;int i=0;
            if(!ct.is_edge(handles.at(edge.points[0]),handles.at(edge.points[1]),f,i) || !f->is_constrained(i))
                topology("original shared constraint was split/lost in finite CT");
            if(ct.is_infinite(f)){const auto other=f->neighbor(i);i=other->index(f);f=other;}
            if(ct.is_infinite(f))topology("actual finite constraint has no incident finite face");
            const auto& va=handles.at(edge.points[0])->point();const auto& vb=handles.at(edge.points[1])->point();
            const auto o=CGAL::orientation(va,vb,f->vertex(i)->point());
            if(o==CGAL::COLLINEAR)topology("finite constraint side has zero orientation");
            const auto neighbor=f->neighbor(i),left=o==CGAL::LEFT_TURN?f:neighbor,right=o==CGAL::LEFT_TURN?neighbor:f;
            const auto& arc=*arcs.at(edge.arcIdentity);const auto l=regionIndex.at(arc.leftRegion),rr=regionIndex.at(arc.rightRegion);
            assign(left,included.count(l)?l:outside);assign(right,included.count(rr)?rr:outside);
        }
        while(!queue.empty()){const auto f=queue.front();queue.pop_front();
            for(int i=0;i<3;++i){consume();if(!f->is_constrained(i))assign(f->neighbor(i),f->info().label);}
        }
    }
    for(auto f=ct.finite_faces_begin();f!=ct.finite_faces_end();++f){consume();
        if(f->info().label==unassigned)topology("finite CT face has no native ownership seed");
        if(f->info().label==outside)continue;
        if(!included.count(f->info().label))topology("excluded native cell promoted to finite material");
        if(out->triangles.size()==settings.maximumTriangles)resource("finite triangle budget exhausted");
        std::array<std::size_t,3> triangle{f->vertex(0)->info(),f->vertex(1)->info(),f->vertex(2)->info()};
        if(CGAL::orientation(f->vertex(0)->point(),f->vertex(1)->point(),f->vertex(2)->point())!=CGAL::LEFT_TURN)
            topology("zero or reversed actual finite CT face");
        const auto first=std::min_element(triangle.begin(),triangle.end());std::rotate(triangle.begin(),first,triangle.end());
        out->triangles.push_back({triangle,f->info().label});
    }
    std::sort(out->triangles.begin(),out->triangles.end(),[](const auto& x,const auto& y){
        return x.points!=y.points?x.points<y.points:x.registryRegionIndex<y.registryRegionIndex;});
    std::set<std::size_t> usedPoints;
    for(const auto region:p.regionIndices()){
        std::vector<std::size_t> own;std::set<std::size_t> vertices;std::map<Pair,std::vector<std::size_t>> edges;
        std::map<Pair,std::int64_t> actualBoundary,expectedBoundary;
        for(std::size_t t=0;t<out->triangles.size();++t)if(out->triangles[t].registryRegionIndex==region){
            consume();own.push_back(t);const auto& v=out->triangles[t].points;
            vertices.insert(v.begin(),v.end());usedPoints.insert(v.begin(),v.end());
            for(int j=0;j<3;++j){edges[ordered(v[j],v[(j+1)%3])].push_back(t);boundary(actualBoundary,v[j],v[(j+1)%3]);}
        }
        if(own.empty())topology("included native finite cell has no actual triangles");
        const auto& cell=r.regions().at(region);
        for(bool inner:{false,true})for(const auto& loop:inner?cell.innerLoops:cell.outerLoops)
            for(const auto& ref:loop.edges){consume();const auto& e=r.edges().at(ref.edgeIndex);
                boundary(expectedBoundary,e.points[ref.reversed?1:0],e.points[ref.reversed?0:1]);}
        trim(actualBoundary);trim(expectedBoundary);
        if(actualBoundary!=expectedBoundary)topology("triangles do not cover every exact original directed outer/inner ring");
        std::map<std::size_t,std::vector<std::size_t>> neighbors;
        for(const auto& e:edges){consume();if(e.second.size()>2)topology("nonmanifold finite region triangle edge");
            if(e.second.size()==2){neighbors[e.second[0]].push_back(e.second[1]);neighbors[e.second[1]].push_back(e.second[0]);}}
        std::set<std::size_t> visited;std::size_t components=0;
        for(const auto t:own)if(!visited.count(t)){++components;std::deque<std::size_t> q{t};visited.insert(t);
            while(!q.empty()){consume();const auto first=q.front();q.pop_front();
                for(const auto next:neighbors[first])if(visited.insert(next).second)q.push_back(next);}}
        const auto chi=static_cast<std::int64_t>(vertices.size())-static_cast<std::int64_t>(edges.size())+static_cast<std::int64_t>(own.size());
        const auto holes=static_cast<std::int64_t>(components)-chi;
        if(holes<0 || components!=cell.outerLoops.size() || std::size_t(holes)!=cell.innerLoops.size())
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","actual finite mesh connectivity/Euler holes disagree with original native cell");
        out->summaries.push_back({region,own.size(),vertices.size(),edges.size(),components,std::size_t(holes),chi});
    }
    if(usedPoints!=std::set<std::size_t>(points.begin(),points.end()))topology("included isolated/shared finite feature was omitted from actual triangles");
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-finite-native-triangulation-v1\nAUTHORITY "<<a.scientificHash()<<' '<<r.scientificHash()<<' '
        <<p.scientificHash()<<' '<<incidence.scientificHash()<<' '<<rotation.scientificHash()<<'\n';
    record<<"POLICY CARTESIAN_CORE_EXPR NO_INTERSECTION_TAG NATIVE_SIDES EXACT_EXPRESSION_IDENTITY_ZERO\nBUDGET "<<settings.maximumTriangles<<' '
        <<settings.maximumAdjacencySteps<<" STEPS "<<steps<<'\n';
    for(const auto& t:out->triangles)record<<"TRIANGLE "<<t.registryRegionIndex<<' '<<t.points[0]<<' '<<t.points[1]<<' '<<t.points[2]<<'\n';
    for(const auto& s:out->summaries)record<<"REGION_SUMMARY "<<s.registryRegionIndex<<' '<<s.triangleCount<<' '<<s.vertexCount<<' '
        <<s.edgeCount<<' '<<s.connectedComponents<<' '<<s.holes<<' '<<s.eulerCharacteristic<<'\n';
    record<<"PARTITION_BEGIN\n"<<p.canonicalRecord()<<"PARTITION_END\nINCIDENCE_BEGIN\n"<<incidence.canonicalRecord()
        <<"INCIDENCE_END\nROTATION_BEGIN\n"<<rotation.canonicalRecord()<<"ROTATION_END\n";
    record<<"actual_finite_triangulation PASSED\ncurve_isotopy NOT_GRANTED\nstrict_thickness NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();out->hash=sha("rank-two-finite-native-triangulation-certificate-v1\n"+out->record);
    detail::requireExactAlgebraicEnvironment();return RankTwoFiniteTriangulation(out);
}
RankTwoFiniteTriangulation reloadRankTwoFiniteTriangulation(const RankTwoArrangement& a,const RankTwoFinitePartition& p,
    const RankTwoChordIncidence& incidence,const RankTwoRegistryRotation& rotation,const std::string& saved,
    const RankTwoFiniteTriangulationSettings& settings){
    auto result=triangulateRankTwoFinitePartition(a,p,incidence,rotation,settings);
    if(result.canonicalRecord()!=saved)invalid("saved finite triangles are not the actual complete native-owned mesh");return result;
}
} // namespace spin
