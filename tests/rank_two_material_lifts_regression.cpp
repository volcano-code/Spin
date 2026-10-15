#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoMaterialLifts.h")
#include "spin/detail/RankTwoMaterialLifts.h"
#include "spin/detail/RankTwoMovingMeshEdge.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong lift failure class");return;}throw std::runtime_error("invalid lift authority/budget accepted");}
spin::RankTwoArcSamplingPlan sampling(const spin::RankTwoArrangement& a){spin::RankTwoArcSamplingPlan result;
    for(const auto& arc:spin::arrangementRecords(a).arcs)if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
        result.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});return result;}
spin::RankTwoInteriorPointRequest request(const spin::RankTwoArrangement& a,spin::RankTwoRegionKind kind,int hsign){
    using namespace spin;const auto& n=arrangementRecords(a);RankTwoInteriorPointRequest r;
    const RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:n.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0 && (!left || x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0 && (!right || x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left && right,"actual native adjacent events absent");r.eventIdentities={left->identity,right->identity};
    for(const auto& c:n.regions)if(c.kind==kind && c.originalSigns[1]==hsign){require(r.regionIdentity.empty(),"ambiguous native active region");r.regionIdentity=c.identity;}
    const auto eventX=[&](const std::string& id){for(const auto& v:n.vertices)if(v.identity==id)return C::reload(v.coordinates[0]).expression();throw std::runtime_error("native event absent");};
    for(const auto& arc:n.arcs){if(arc.vertical || arc.ends[0].location!=RankTwoArcEndLocation::Finite || arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        if(!(eventX(*arc.ends[0].vertexIdentity)<0 && eventX(*arc.ends[1].vertexIdentity)>0))continue;
        if(kind==RankTwoRegionKind::Cut && C::reload(arc.interiorWitness[1]).expression().sign()<=0)continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous lower native gap");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous upper native gap");r.boundaryArcIdentities[1]=arc.identity;}}
    return r;
}
void check(const spin::PreparedRankTwoSource& source,const spin::RankTwoArrangement& a,
    const spin::RankTwoInteriorMesh& mesh,const char* name,std::ostream& out,bool requirePassed){
    using namespace spin;std::cerr<<"material-test:"<<name<<":build\n";const auto lifts=buildRankTwoMaterialLifts(source,a,mesh);
    std::cerr<<"material-test:"<<name<<":independent-native-kernel\n";
    const auto& registry=mesh.baseMesh().coverage().partition().registry();
    require(lifts.sourceIdentity()==source.inputHash() && lifts.meshIdentity()==mesh.scientificHash(),"lift detached from actual source/mesh");
    require(lifts.triangles().size()==mesh.triangles().size(),"Empty region or failed active triangle omitted");
    std::size_t active=0,failed=0,empty=0,fullCut=0,cutEmpty=0;
    using Pair=std::pair<std::size_t,std::size_t>;std::map<Pair,std::vector<std::size_t>> uses;
    for(const auto& t:mesh.triangles())for(int j=0;j<3;++j)uses[std::minmax(t.points[j],t.points[(j+1)%3])].push_back(t.registryRegionIndex);
    for(std::size_t i=0;i<mesh.triangles().size();++i){const auto& t=mesh.triangles()[i];const auto& w=lifts.triangles()[i];
        require(w.triangleIndex==i && w.registryRegionIndex==t.registryRegionIndex,"lift witness uses stale triangle index");
        const auto kind=registry.regions()[t.registryRegionIndex].kind;
        if(kind==RankTwoRegionKind::Empty){++empty;require(!w.active,"Empty region promoted to material");continue;}
        ++active;require(w.active,"active material region omitted");
        std::array<E,3> lo,hi;
        for(int j=0;j<3;++j){lo[j]=lifts.lower(t.registryRegionIndex,t.points[j]);hi[j]=lifts.upper(t.points[j]);
            require(hi[j].cmp(lo[j])>=0,"negative actual corner thickness");
            const auto& adjacent=uses.at(std::minmax(t.points[j],t.points[(j+1)%3]));
            if(adjacent.size()==2){const auto first=registry.regions()[adjacent[0]].kind,second=registry.regions()[adjacent[1]].kind;
                if((first==RankTwoRegionKind::Full && second==RankTwoRegionKind::Cut) || (second==RankTwoRegionKind::Full && first==RankTwoRegionKind::Cut)){
                    ++fullCut;require(!w.materialBoundary[j],"Full/Cut shared interface incorrectly exempted from thickness");
                    for(const auto p:{t.points[j],t.points[(j+1)%3]})require(lifts.lower(adjacent[0],p)==lifts.lower(adjacent[1],p),"Full/Cut lower lifts disagree");}
                if(first==RankTwoRegionKind::Empty || second==RankTwoRegionKind::Empty){++cutEmpty;require(w.materialBoundary[j],"Cut/Empty true material boundary mislabeled internal");
                    require(hi[j]==lo[j] && lifts.upper(t.points[(j+1)%3])==lifts.lower(t.registryRegionIndex,t.points[(j+1)%3]),"Cut/Empty boundary does not collapse");}}}
        bool passed=true;try{detail::certifyRankTwoTriangleThickness(lo,hi,w.materialBoundary);}catch(const RankTwoSourceError& e){require(e.code()=="SOURCE_TOPOLOGY_INVALID","unexpected kernel failure class");passed=false;}
        require(w.strictPassed==passed,"lift wrapper falsified actual kernel acceptance");failed+=!passed;
    }
    require(lifts.allActiveStrictPassed()==(active>0 && failed==0),"whole thickness grant ignores failed actual triangle");
    if(requirePassed)require(active==66 && failed==0,"repaired real disk lost strict thickness");
    if(std::string(name)=="annulus")require(empty>0 && fullCut>0 && cutEmpty>0,"native Empty/Cut/Full boundary coverage missing");
    if(std::string(name)=="unrepaired-disk")require(failed>0,"zero disk ears/edges silently granted thickness");
    RankTwoMaterialLiftSettings zero;zero.maximumPredicates=0;
    rejects([&]{buildRankTwoMaterialLifts(source,a,mesh,zero);},"RESOURCE_LIMIT");
    const auto foreign=PreparedRankTwoSource::fromExact(source.exactCoefficients(),"foreign-lift-input-frame");
    rejects([&]{buildRankTwoMaterialLifts(foreign,a,mesh);},"INPUT_IDENTITY_MISMATCH");
    out<<"CASE "<<name<<" MATERIAL_LIFTS "<<lifts.scientificHash()<<'\n'<<lifts.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" active="<<active<<" failed="<<failed<<" empty="<<empty<<" FullCut="<<fullCut<<" CutEmpty="<<cutEmpty<<'\n';
}
void tests(const char* path){using namespace spin;std::ofstream out(path);require(bool(out),"cannot save actual lifts");out<<"rank-two-material-lifts-test-v1\n";
    detail::ScopedRankTwoDiagnostics observations([](const char* phase,void*)noexcept{
        if(std::string(phase).find("material-lift:")==0)std::cerr<<phase<<'\n';},nullptr);
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));const auto coverage=refineRankTwoFiniteArcCoverage(a,sampling(a),{},Q(1,32));
    const auto base=triangulateRankTwoCoveredFiniteDomain(a,coverage);
    const auto center=buildRankTwoMovingInteriorPoint(a,request(a,RankTwoRegionKind::Full,-1));
    const auto parent=insertRankTwoMovingInteriorNodes(a,base,{center});check(disk,a,parent,"unrepaired-disk",out,false);
    using Pair=std::pair<std::size_t,std::size_t>;std::map<Pair,int> edges;
    for(const auto& t:parent.triangles())for(int j=0;j<3;++j)++edges[std::minmax(t.points[j],t.points[(j+1)%3])];
    std::vector<RankTwoMovingMeshEdgePoint> nodes;
    for(const auto& e:edges)if(e.second==2 && e.first.first!=24 && e.first.second!=24){RankTwoMovingMeshEdgeRequest r;r.pointIndices={e.first.first,e.first.second};nodes.push_back(buildRankTwoMovingMeshEdgePoint(a,parent,r));}
    const auto repaired=insertRankTwoMovingEdgeNodes(a,parent,nodes);check(disk,a,repaired,"repaired-disk",out,true);
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));const auto covered=refineRankTwoFiniteArcCoverage(b,sampling(b),{},Q(1,32));
    const auto ring=triangulateRankTwoCoveredFiniteDomain(b,covered);std::vector<RankTwoMovingInteriorPoint> interior;
    for(const auto [kind,sign]:{std::pair{RankTwoRegionKind::Full,-1},std::pair{RankTwoRegionKind::Cut,-1},std::pair{RankTwoRegionKind::Cut,1}})
        interior.push_back(buildRankTwoMovingInteriorPoint(b,request(b,kind,sign)));
    const auto actual=insertRankTwoMovingInteriorNodes(b,ring,interior);check(ann,b,actual,"annulus",out,false);
    out<<"closed_3d_source NOT_GRANTED\nsource_motion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";out.flush();require(bool(out),"incomplete lift evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoMaterialLifts.h")
    if(argc!=2)throw std::invalid_argument("material lift evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: source-bound actual mesh lifts and native material-boundary classification");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
