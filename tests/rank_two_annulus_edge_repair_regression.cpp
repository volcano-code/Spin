#include "spin/detail/RankTwoMaterialLifts.h"
#include "spin/detail/RankTwoMovingMeshEdge.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>

namespace {
using Q=CGAL::Gmpq;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
spin::RankTwoInteriorPointRequest request(const spin::RankTwoArrangement& a,spin::RankTwoRegionKind kind,int sign){
    using namespace spin;const auto& n=arrangementRecords(a);RankTwoInteriorPointRequest r;
    const RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:n.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0 && (!left || x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0 && (!right || x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left && right,"native adjacent events missing");r.eventIdentities={left->identity,right->identity};
    for(const auto& c:n.regions)if(c.kind==kind && c.originalSigns[1]==sign){require(r.regionIdentity.empty(),"ambiguous region");r.regionIdentity=c.identity;}
    const auto x=[&](const std::string& id){for(const auto& v:n.vertices)if(v.identity==id)return C::reload(v.coordinates[0]).expression();throw std::runtime_error("native event missing");};
    for(const auto& arc:n.arcs){if(arc.vertical || arc.ends[0].location!=RankTwoArcEndLocation::Finite || arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        if(!(x(*arc.ends[0].vertexIdentity)<0 && x(*arc.ends[1].vertexIdentity)>0))continue;
        if(kind==RankTwoRegionKind::Cut && C::reload(arc.interiorWitness[1]).expression().sign()<=0)continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous lower boundary");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous upper boundary");r.boundaryArcIdentities[1]=arc.identity;}}
    return r;
}
using Pair=std::pair<std::size_t,std::size_t>;
std::set<Pair> zeroEdges(const spin::RankTwoInteriorMesh& mesh,const spin::RankTwoMaterialLifts& lifts){
    std::set<Pair> result;
    for(const auto& w:lifts.triangles())for(int j=0;j<3;++j)if(w.zeroInternalEdges[j]){
        const auto& t=mesh.triangles().at(w.triangleIndex);result.insert(std::minmax(t.points[j],t.points[(j+1)%3]));}
    return result;
}
std::size_t failed(const spin::RankTwoMaterialLifts& lifts){std::size_t count=0;for(const auto& w:lifts.triangles())count+=w.active && !w.strictPassed;return count;}
void tests(const char* path,const std::string& mode){using namespace spin;
    const bool asymmetric=mode=="asymmetric",all=mode=="all";
    detail::ScopedRankTwoDiagnostics observations([](const char* phase,void*)noexcept{
        const std::string s(phase);if(s.find("annulus-edge:")==0)std::cerr<<s<<'\n';},nullptr);
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(3,2));RankTwoArcSamplingPlan samples;
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
        samples.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto coverage=refineRankTwoFiniteArcCoverage(a,samples,{},Q(1,32));
    const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);std::vector<RankTwoMovingInteriorPoint> inside;
    for(const auto [kind,s]:{std::pair{RankTwoRegionKind::Full,-1},std::pair{RankTwoRegionKind::Cut,-1},std::pair{RankTwoRegionKind::Cut,1}})
        inside.push_back(buildRankTwoMovingInteriorPoint(a,request(a,kind,s)));
    const auto parent=insertRankTwoMovingInteriorNodes(a,base,inside);
    const auto old=buildRankTwoMaterialLifts(source,a,parent);const auto inventory=zeroEdges(parent,old);
    require(inventory.size()==232 && failed(old)==258,"frozen annulus initial negative scope changed");
    std::ofstream out(path);require(bool(out),"cannot save annulus repair evidence");
    out<<"rank-two-annulus-edge-repair-test-v2\nPOLICY "<<(all?"ALL_ACTUAL_ZERO_EDGES_INDEX_ORDER":asymmetric?"FROZEN_ACTUAL_ASYMMETRIC_EDGE_10_215":"FIRST_ACTUAL_ZERO_EDGE_INDEX_ORDER")<<" FIXED_WEIGHT 1/2\n"
       <<"CASE annulus-parent MATERIAL_LIFTS "<<old.scientificHash()<<'\n'<<old.canonicalRecord()<<"CASE_END\n";
    out.flush();require(bool(out),"cannot preserve actual parent before edge attempts");
    // Freeze selection BEFORE any attempted point construction: lexicographic
    // actual parent indices, exact weight1/2. No successful-candidate filtering.
    const auto first=asymmetric?Pair(10,215):*inventory.begin();
    require(inventory.count(first) && (asymmetric || first==Pair(10,187)),"actual frozen zero edge changed");
    std::vector<RankTwoMovingMeshEdgePoint> nodes;
    const auto selected=all?inventory:std::set<Pair>{first};
    for(const auto e:selected){
        RankTwoMovingMeshEdgeRequest r;r.pointIndices={e.first,e.second};
        std::cerr<<"annulus-edge:frozen-point:"<<e.first<<':'<<e.second<<'\n';
        nodes.push_back(buildRankTwoMovingMeshEdgePoint(a,parent,r));
        out<<"COMPLETED_NODE "<<e.first<<' '<<e.second<<' '<<nodes.back().scientificHash()<<'\n';
        out.flush();require(bool(out),"incomplete actual edge checkpoint");
    }
    std::cerr<<"annulus-edge:first-insertion\n";
    const auto refined=insertRankTwoMovingEdgeNodes(a,parent,nodes);
    std::cerr<<"annulus-edge:actual-refined-lifts\n";
    const auto result=buildRankTwoMaterialLifts(source,a,refined);
    require(refined.points().size()==2183+nodes.size() && refined.triangles().size()==2856+2*nodes.size(),"edge insertion omitted a side or changed unrelated faces");
    const auto remaining=zeroEdges(refined,result);
    // Parent face225 has TWO zero internal sides; after this one repair its
    // other side10--1653 is still zero. Parent face233 is fully repaired.
    // Therefore258->257, not258->256. Do not erase the remaining failure.
    // The asymmetric edge's parent faces227,236 EACH have three zero sides.
    // Its four children still contain the four other unrepaired zero sides;
    // failed-face count258->260 is honest, not a reason to erase those faces.
    const auto expectedFailed=all?0u:asymmetric?260u:257u;
    require(!remaining.count(first) && remaining.size()==(all?0u:231u) && failed(result)==expectedFailed,"actual shared edge repair did not restore strict local thickness");
    require(result.allActiveStrictPassed()==all,"partial repair promoted incomplete annulus or complete repair retains pinch");
    require(parent.regionSummaries().size()==refined.regionSummaries().size(),"region omitted");
    for(std::size_t i=0;i<parent.regionSummaries().size();++i){const auto& p=parent.regionSummaries()[i];const auto& q=refined.regionSummaries()[i];
        require(p.registryRegionIndex==q.registryRegionIndex && p.connectedComponents==q.connectedComponents && p.holes==q.holes && p.eulerCharacteristic==q.eulerCharacteristic,"actual regional holes changed");}
    out<<"CASE "<<(all?"annulus-all-edges":asymmetric?"annulus-asymmetric-edge":"annulus-single-edge")<<" MATERIAL_LIFTS "<<result.scientificHash()<<'\n'<<result.canonicalRecord()<<"CASE_END\n"
       <<"SCOPE INITIAL_ZERO_EDGES 232 FINAL_ZERO_EDGES "<<remaining.size()<<" INITIAL_FAILED 258 FINAL_FAILED "<<expectedFailed<<"\n"
       <<"all_annulus_strict "<<(all?"METHOD_FIXTURE_PASSED":"NOT_GRANTED")<<"\nclosed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"incomplete annulus evidence");
    std::cout<<"PASS actual annulus frozen shared edge repair; "<<remaining.size()<<"zero edges/"<<expectedFailed<<"failed faces remain; no full carrier grant\n";
}
}
int main(int argc,char** argv){try{if(argc!=2 && !(argc==3 && (std::string(argv[2])=="asymmetric" || std::string(argv[2])=="all")))throw std::invalid_argument("annulus edge evidence path and optional asymmetric/all mode required");tests(argv[1],argc==3?argv[2]:"single");return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
