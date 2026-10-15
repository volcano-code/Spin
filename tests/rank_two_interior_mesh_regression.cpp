#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoInteriorMesh.h")
#include "spin/detail/RankTwoInteriorMesh.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){bool caught=false;
    try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);}
spin::RankTwoArcSamplingPlan plan(const spin::RankTwoArrangement& a){spin::RankTwoArcSamplingPlan p;
    for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)
            p.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});return p;}
spin::RankTwoInteriorPointRequest request(const spin::RankTwoArrangement& a,spin::RankTwoRegionKind kind,int hsign){
    const auto& n=spin::arrangementRecords(a);spin::RankTwoInteriorPointRequest r;
    const spin::RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:n.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0 && (!left || x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0 && (!right || x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left && right,"frozen mesh fixture lacks adjacent native events");r.eventIdentities={left->identity,right->identity};
    for(const auto& region:n.regions)if(region.kind==kind && region.originalSigns[1]==hsign){
        require(r.regionIdentity.empty(),"fixture must identify one active native region");r.regionIdentity=region.identity;}
    require(!r.regionIdentity.empty(),"fixture active native region missing");
    const auto x=spin::detail::sampleRankTwoAlgebraicCoordinate(C::reload(left->coordinates[0]),C::reload(right->coordinates[0]),Q(1,2)).expression();
    const auto eventX=[&](const std::string& id){for(const auto& v:n.vertices)if(v.identity==id)return C::reload(v.coordinates[0]).expression();throw std::runtime_error("missing native endpoint");};
    for(const auto& arc:n.arcs){if(arc.vertical || arc.ends[0].location!=spin::RankTwoArcEndLocation::Finite || arc.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        if(!(eventX(*arc.ends[0].vertexIdentity)<x && x<eventX(*arc.ends[1].vertexIdentity)))continue;
        if(kind==spin::RankTwoRegionKind::Cut && C::reload(arc.interiorWitness[1]).expression().sign()<=0)continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous frozen lower fiber");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous frozen upper fiber");r.boundaryArcIdentities[1]=arc.identity;}}
    return r;
}
void check(const spin::RankTwoArrangement& a,const spin::RankTwoCoveredFiniteMesh& base,
    const std::vector<spin::RankTwoMovingInteriorPoint>& nodes,std::ostream& out,const char* name,bool annulus){
    using namespace spin;const auto mesh=insertRankTwoMovingInteriorNodes(a,base,nodes);
    const auto& p=base.coverage().partition();const auto& r=p.registry();
    require(mesh.baseMesh().scientificHash()==base.scientificHash(),"interior insertion dropped the owned refined domain");
    require(mesh.points().size()==r.points().size()+nodes.size() && mesh.triangles().size()==base.triangulation().triangles().size()+2*nodes.size(),
        "moving nodes were attached as metadata instead of actual incident retriangulation");
    std::set<std::size_t> used;std::map<std::size_t,std::map<std::pair<std::size_t,std::size_t>,int>> boundary;
    for(const auto& t:mesh.triangles()){
        used.insert(t.points.begin(),t.points.end());
        for(int j=0;j<3;++j){auto u=t.points[j],v=t.points[(j+1)%3];boundary[t.registryRegionIndex][std::minmax(u,v)]+=u<v?1:-1;}}
    auto wantedUsed=std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end());
    for(std::size_t i=r.points().size();i<mesh.points().size();++i){wantedUsed.insert(i);require(used.count(i),"new moving point absent from actual triangles");}
    require(used==wantedUsed,"unpublished or omitted shared source point in interior mesh");
    std::size_t holes=0,empty=0,cut=0,full=0;
    for(const auto& s:mesh.regionSummaries()){
        const auto& region=r.regions().at(s.registryRegionIndex);holes+=s.holes;
        empty+=region.kind==RankTwoRegionKind::Empty;cut+=region.kind==RankTwoRegionKind::Cut;full+=region.kind==RankTwoRegionKind::Full;
        require(s.holes==region.innerLoops.size(),"moving insertion filled a native hole");
        std::map<std::pair<std::size_t,std::size_t>,int> wanted;
        for(bool inner:{false,true})for(const auto& loop:inner?region.innerLoops:region.outerLoops)
            for(const auto& ref:loop.edges){const auto& e=r.edges().at(ref.edgeIndex);auto u=e.points[0],v=e.points[1];if(ref.reversed)std::swap(u,v);wanted[std::minmax(u,v)]+=u<v?1:-1;}
        auto actual=boundary.at(s.registryRegionIndex);
        for(auto i=actual.begin();i!=actual.end();)if(!i->second)i=actual.erase(i);else ++i;
        for(auto i=wanted.begin();i!=wanted.end();)if(!i->second)i=wanted.erase(i);else ++i;
        require(actual==wanted,"moving insertion changed an original directed covered ring");}
    if(annulus)require(holes==3 && empty==1 && cut==2 && full==1,"interior insertion used parity or lost annulus labels");
    else {
        require(holes==0 && full==1 && used.size()==25 && mesh.triangles().size()==24,"actual disk subdivision differs from Euler-derived count");
        // A center point is positive, but unchanged boundary-only ears and
        // diagonals still fail. Do not mark the whole carrier as passed.
        std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
        for(const auto& t:mesh.triangles())for(int j=0;j<3;++j)++edges[std::minmax(t.points[j],t.points[(j+1)%3])];
        std::size_t zeroTriangles=0,zeroInternalEdges=0,positiveTriangles=0,strictPassed=0,rejectedPositive=0;
        for(const auto& t:mesh.triangles()){
            std::array<E,3> lower,upper;std::array<bool,3> materialBoundary;
            bool allZero=true;
            for(int j=0;j<3;++j){const auto& xy=mesh.points().at(t.points[j]).coordinates;
                const E D=E(1)-xy[0].expression()*xy[0].expression()-xy[1].expression()*xy[1].expression();
                require(D.sign()>=0,"actual inserted disk corner leaves the projection domain");upper[j]=CORE::sqrt(D);lower[j]=-upper[j];allZero=allZero && D.sign()==0;
                materialBoundary[j]=edges.at(std::minmax(t.points[j],t.points[(j+1)%3]))==1;}
            bool zeroInternal=false;
            for(int j=0;j<3;++j)zeroInternal=zeroInternal || (!materialBoundary[j] &&
                (upper[j]-lower[j]+upper[(j+1)%3]-lower[(j+1)%3]).sign()==0);
            if(allZero){++zeroTriangles;rejects([&]{detail::certifyRankTwoTriangleThickness(lower,upper,materialBoundary);},"SOURCE_TOPOLOGY_INVALID","unchanged zero ear promoted to carrier");}
            else {++positiveTriangles;
                if(zeroInternal){++rejectedPositive;rejects([&]{detail::certifyRankTwoTriangleThickness(lower,upper,materialBoundary);},"SOURCE_TOPOLOGY_INVALID","positive face with zero internal seam promoted to carrier");}
                else {++strictPassed;detail::certifyRankTwoTriangleThickness(lower,upper,materialBoundary);}}}
        for(const auto& entry:edges)if(entry.second==2){
            const auto& u=mesh.points()[entry.first.first].coordinates;const auto& v=mesh.points()[entry.first.second].coordinates;
            const E du=E(1)-u[0].expression()*u[0].expression()-u[1].expression()*u[1].expression();
            const E dv=E(1)-v[0].expression()*v[0].expression()-v[1].expression()*v[1].expression();
            if(du.sign()==0 && dv.sign()==0)++zeroInternalEdges;}
        require(positiveTriangles>0 && zeroTriangles>0 && zeroInternalEdges>0 && rejectedPositive>0,"single-node insertion incorrectly claimed complete thickness repair");
        out<<"THICKNESS_SCOPE full-disk ACTIVE "<<mesh.triangles().size()<<" POSITIVE "<<positiveTriangles<<" ZERO_TRIANGLES "<<zeroTriangles<<" ZERO_INTERNAL_EDGES "<<zeroInternalEdges
           <<" STRICT_PASSED "<<strictPassed<<" REJECTED_POSITIVE "<<rejectedPositive<<" NOT_GRANTED\n";
    }
    RankTwoInteriorMeshSettings zero;zero.maximumPredicates=0;
    rejects([&]{insertRankTwoMovingInteriorNodes(a,base,nodes,zero);},"RESOURCE_LIMIT","zero incidence proof budget granted interior mesh");
    auto duplicate=nodes;duplicate.push_back(nodes.front());
    rejects([&]{insertRankTwoMovingInteriorNodes(a,base,duplicate);},"INPUT_IDENTITY_MISMATCH","duplicate moving point identity accepted");
    out<<"CASE "<<name<<" INTERIOR_MESH "<<mesh.scientificHash()<<'\n'<<mesh.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" points="<<used.size()<<" triangles="<<mesh.triangles().size()<<" holes="<<holes<<" labels="<<empty<<'/'<<cut<<'/'<<full<<'\n';
}
void tests(const char* path){using namespace spin;std::ostringstream out;out<<"rank-two-interior-mesh-test-v1\n";
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-coverage-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));const auto first=refineRankTwoFiniteArcCoverage(a,plan(a),{},Q(1,32));
    const auto base=triangulateRankTwoCoveredFiniteDomain(a,first);
    const auto center=buildRankTwoMovingInteriorPoint(a,request(a,RankTwoRegionKind::Full,-1));
    check(a,base,{center},out,"full-disk",false);
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-coverage-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));const auto second=refineRankTwoFiniteArcCoverage(b,plan(b),{},Q(1,32));
    const auto ring=triangulateRankTwoCoveredFiniteDomain(b,second);
    std::vector<RankTwoMovingInteriorPoint> nodes;
    for(const auto [kind,sign]:{std::pair{RankTwoRegionKind::Full,-1},std::pair{RankTwoRegionKind::Cut,-1},std::pair{RankTwoRegionKind::Cut,1}})
        nodes.push_back(buildRankTwoMovingInteriorPoint(b,request(b,kind,sign)));
    check(b,ring,nodes,out,"annulus",true);
    rejects([&]{insertRankTwoMovingInteriorNodes(b,base,{center});},"INPUT_IDENTITY_MISMATCH","foreign covered domain granted insertion");
    rejects([&]{insertRankTwoMovingInteriorNodes(a,base,nodes);},"INPUT_IDENTITY_MISMATCH","foreign node authority granted insertion");
    std::ofstream saved(path);require(bool(saved),"cannot save actual interior mesh");saved<<out.str();saved.flush();require(bool(saved),"incomplete actual interior mesh output");}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoInteriorMesh.h")
    if(argc!=2)throw std::invalid_argument("interior mesh evidence path required");
    const auto start=std::chrono::steady_clock::now();const auto observe=+[](const char* phase,void* context)noexcept{try{
        const auto& t=*static_cast<const std::chrono::steady_clock::time_point*>(context);
        std::cerr<<"DIAGNOSTIC_NS "<<std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-t).count()<<' '<<phase<<'\n';}catch(...){}};
    spin::detail::ScopedRankTwoDiagnostics observation(observe,const_cast<std::chrono::steady_clock::time_point*>(&start));tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: moving native nodes in the actual covered common triangle mesh");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
