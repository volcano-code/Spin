#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoMovingMeshEdge.h")
#include "spin/detail/RankTwoMovingMeshEdge.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong moving-edge failure class");return;}throw std::runtime_error("invalid moving edge accepted");}
spin::RankTwoInteriorPointRequest centerRequest(const spin::RankTwoArrangement& a){
    const auto& n=spin::arrangementRecords(a);spin::RankTwoInteriorPointRequest r;
    const spin::RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:n.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0 && (!left || x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0 && (!right || x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left && right,"actual disk native events missing");r.eventIdentities={left->identity,right->identity};
    for(const auto& cell:n.regions)if(cell.kind==spin::RankTwoRegionKind::Full){require(r.regionIdentity.empty(),"ambiguous frozen full disk cell");r.regionIdentity=cell.identity;}
    const auto eventX=[&](const std::string& id){for(const auto& v:n.vertices)if(v.identity==id)return C::reload(v.coordinates[0]).expression();throw std::runtime_error("native event absent");};
    for(const auto& arc:n.arcs){if(arc.vertical || arc.ends[0].location!=spin::RankTwoArcEndLocation::Finite || arc.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        if(!(eventX(*arc.ends[0].vertexIdentity)<0 && eventX(*arc.ends[1].vertexIdentity)>0))continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous lower disk fiber");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous upper disk fiber");r.boundaryArcIdentities[1]=arc.identity;}}
    return r;
}
void tests(const char* path){using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(1));RankTwoArcSamplingPlan sampling;
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
        sampling.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto coverage=refineRankTwoFiniteArcCoverage(a,sampling,{},Q(1,32));
    const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    const auto center=buildRankTwoMovingInteriorPoint(a,centerRequest(a));
    const auto parent=insertRankTwoMovingInteriorNodes(a,base,{center});
    using Pair=std::pair<std::size_t,std::size_t>;
    std::map<Pair,std::size_t> oldEdges;
    for(const auto& t:parent.triangles())for(int j=0;j<3;++j)++oldEdges[std::minmax(t.points[j],t.points[(j+1)%3])];
    const auto d=[&](const auto& point){const auto& xy=point.coordinates;return E(1)-xy[0].expression()*xy[0].expression()-xy[1].expression()*xy[1].expression();};
    std::vector<RankTwoMovingMeshEdgePoint> nodes;
    for(const auto& e:oldEdges)if(e.second==2 && d(parent.points()[e.first.first]).sign()==0 && d(parent.points()[e.first.second]).sign()==0){
        RankTwoMovingMeshEdgeRequest request;request.pointIndices={e.first.first,e.first.second};
        const auto node=buildRankTwoMovingMeshEdgePoint(a,parent,request);
        require(node.parentMeshIdentity()==parent.scientificHash(),"moving edge lost actual parent identity");nodes.push_back(node);}
    require(nodes.size()==21,"frozen disk zero-edge inventory changed");
    rejects([&]{buildRankTwoMovingMeshEdgePoint(a,parent,RankTwoMovingMeshEdgeRequest{});},"INPUT_IDENTITY_MISMATCH");
    auto third=nodes.front().definition();third.fraction=Q(1,3);
    const auto forward=buildRankTwoMovingMeshEdgePoint(a,parent,third);
    std::swap(third.pointIndices[0],third.pointIndices[1]);third.fraction=Q(2,3);
    const auto reverse=buildRankTwoMovingMeshEdgePoint(a,parent,third);
    require(forward.scientificHash()==reverse.scientificHash(),"reversed endpoint DAG changed actual affine point identity");
    for(int j=0;j<2;++j)require(forward.coordinates()[j].expression()==reverse.coordinates()[j].expression(),"reversed endpoint DAG changed actual affine coordinates");
    const auto repaired=insertRankTwoMovingEdgeNodes(a,parent,nodes);
    require(repaired.points().size()==46 && repaired.triangles().size()==66,"moving edge samples are metadata or incident faces were deleted");
    std::set<std::size_t> used;std::map<Pair,std::size_t> edges;std::map<Pair,int> boundary;
    for(const auto& t:repaired.triangles()){
        require(t.registryRegionIndex==parent.regionSummaries()[0].registryRegionIndex,"moving edge changed native label");
        used.insert(t.points.begin(),t.points.end());
        for(int j=0;j<3;++j){const auto u=t.points[j],v=t.points[(j+1)%3];++edges[std::minmax(u,v)];boundary[std::minmax(u,v)]+=u<v?1:-1;}}
    require(used.size()==46,"actual repaired mesh omits certified point");
    std::vector<E> height;for(const auto& p:repaired.points()){const auto D=d(p);require(D.sign()>=0,"edge point outside disk");height.push_back(CORE::sqrt(D));}
    for(std::size_t j=parent.points().size();j<repaired.points().size();++j)require(height[j].sign()>0,"new edge sample has no positive thickness");
    for(const auto& t:repaired.triangles()){
        std::array<E,3> lower,upper;std::array<bool,3> flags;
        for(int j=0;j<3;++j){lower[j]=-height[t.points[j]];upper[j]=height[t.points[j]];
            flags[j]=edges.at(std::minmax(t.points[j],t.points[(j+1)%3]))==1;
            if(flags[j])require(height[t.points[j]].sign()==0 && height[t.points[(j+1)%3]].sign()==0,"invented material boundary");}
        detail::certifyRankTwoTriangleThickness(lower,upper,flags);}
    for(const auto& e:oldEdges)if(e.second==1)require(edges[e.first]==1,"moving edge changed original covered disk boundary");
    for(const auto& e:edges)require(e.second==1 || e.second==2,"nonmanifold repaired edge");
    require(repaired.regionSummaries().size()==1 && repaired.regionSummaries()[0].holes==0 && repaired.regionSummaries()[0].eulerCharacteristic==1,"disk topology changed");
    RankTwoMovingMeshEdgeSettings zero;zero.maximumNativeChecks=0;
    RankTwoMovingMeshEdgeRequest r;r.pointIndices={oldEdges.begin()->first.first,oldEdges.begin()->first.second};
    rejects([&]{buildRankTwoMovingMeshEdgePoint(a,parent,r,zero);},"RESOURCE_LIMIT");
    const auto exterior=std::find_if(oldEdges.begin(),oldEdges.end(),[](const auto& e){return e.second==1;});
    require(exterior!=oldEdges.end(),"disk boundary edge missing");r.pointIndices={exterior->first.first,exterior->first.second};
    rejects([&]{buildRankTwoMovingMeshEdgePoint(a,parent,r);},"SOURCE_TOPOLOGY_INVALID");
    r.fraction=0;rejects([&]{buildRankTwoMovingMeshEdgePoint(a,parent,r);},"INPUT_IDENTITY_MISMATCH");
    auto duplicate=nodes;duplicate.push_back(nodes.front());rejects([&]{insertRankTwoMovingEdgeNodes(a,parent,duplicate);},"INPUT_IDENTITY_MISMATCH");
    const auto foreignParent=insertRankTwoMovingInteriorNodes(a,base,{});
    rejects([&]{insertRankTwoMovingEdgeNodes(a,foreignParent,nodes);},"INPUT_IDENTITY_MISMATCH");
    RankTwoInteriorMeshSettings exhausted;exhausted.maximumPredicates=0;
    rejects([&]{insertRankTwoMovingEdgeNodes(a,parent,nodes,exhausted);},"RESOURCE_LIMIT");
    std::ofstream out(path);require(bool(out),"cannot save actual moving edge mesh");
    out<<"rank-two-moving-mesh-edge-test-v1\nCASE full-disk EDGE_MESH "<<repaired.scientificHash()<<'\n'<<repaired.canonicalRecord()<<"CASE_END\n"
       <<"ACTUAL_DISK_THICKNESS PASS 66 TRIANGLES 21 INSERTED_EDGE_POINTS\nclosed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"incomplete moving edge mesh evidence");
    // Break caught: repeatedly scanning every unrelated face for each edge
    // consumes the fixed budget before the final topology checks. This is
    // the actual21-edge disk, not a fake incidence list. 24 parent faces,
    // 21 two-sided splits and66 final faces are independently fixed above.
    RankTwoInteriorMeshSettings bounded;bounded.maximumPredicates=1000;
    const auto withinBudget=insertRankTwoMovingEdgeNodes(a,parent,nodes,bounded);
    require(withinBudget.points().size()==46 && withinBudget.triangles().size()==66,
        "fixed-budget insertion omitted a repair or changed actual mesh");
    for(std::size_t i=0;i<repaired.triangles().size();++i)
        require(withinBudget.triangles()[i].points==repaired.triangles()[i].points
            &&withinBudget.triangles()[i].registryRegionIndex==repaired.triangles()[i].registryRegionIndex,
            "budget-index path changed deterministic triangle/owner sequence");
    auto reversed=nodes;std::reverse(reversed.begin(),reversed.end());
    const auto reverseBudget=insertRankTwoMovingEdgeNodes(a,parent,reversed,bounded);
    require(reverseBudget.canonicalRecord()==withinBudget.canonicalRecord(),
        "input order changed indexed edge subdivision or predicate accounting");
    RankTwoInteriorMeshSettings one;one.maximumPredicates=1;
    rejects([&]{insertRankTwoMovingEdgeNodes(a,parent,nodes,one);},"RESOURCE_LIMIT");
    out<<"CASE fixed-budget-disk EDGE_MESH "<<withinBudget.scientificHash()<<'\n'
       <<withinBudget.canonicalRecord()<<"CASE_END\n"
       <<"ACTUAL_EDGE_BUDGET PASS 1000; ALL21 REPAIRS;46POINTS;66TRIANGLES;REVERSED_ORDER_IDENTICAL;ONE_REJECTED\n";
    out.flush();require(bool(out),"incomplete fixed-budget evidence");
    std::cout<<"PASS actual disk46points/66triangles;21moving edge points; ALL actual triangles/internal edges positive; closed3D/source motion NOT_GRANTED\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoMovingMeshEdge.h")
    if(argc!=2)throw std::invalid_argument("moving-edge evidence path required");
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: certified moving internal-edge nodes in the actual common mesh");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
