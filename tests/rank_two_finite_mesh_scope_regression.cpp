#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include "spin/detail/RankTwoChordSweep.h"
#if defined(RANK_TWO_FINITE_PARTITION_INCIDENCE_V1)
#include "spin/detail/RankTwoFiniteTriangulation.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(F action,const char* code,const char* why){
    bool caught=false;try{action();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
spin::RankTwoArcRegistry registry(const spin::RankTwoArrangement& a,int denominator){
    spin::RankTwoArcSamplingPlan plan;
    for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite){
            std::vector<Q> fractions;for(int i=1;i<denominator;++i)fractions.emplace_back(i,denominator);
            plan.emplace(arc.identity,std::move(fractions));
        }
    return spin::buildRankTwoArcRegistry(a,plan);
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-mesh-scope-disk-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(source,global,cone);
    const auto r=registry(a,4);const auto p=buildRankTwoFinitePartition(a,r);
    const auto n=p.pointIndices().size(),m=p.edgeIndices().size();
    require(n>=3 && m>=3 && n<r.points().size() && m<r.edges().size(),"fixture does not actually restrict the native graph");
    const auto incidence=certifyRankTwoChordIncidenceSweep(a,p);
    require(incidence.partitionIdentity()==p.scientificHash(),"incidence scope lost the immutable partition");
    require(incidence.pointPairCount()==n*(n-1)/2 && incidence.pointEdgeCount()==m*(n-2)
        && incidence.edgePairCount()==m*(m-1)/2,"partition sweep silently omitted in-domain pairs");
    const auto rotation=certifyRankTwoRegistryRotation(a,p,incidence);
    require(rotation.partitionIdentity()==p.scientificHash() && rotation.pointCount()==n
        && rotation.finiteDirectionCount()==2*m,"rotation consumed wrong finite domain");
    require(rotation.omittedOutsidePartitionDirectionCount()>0,"excluded finite directions incorrectly reported as infinite");
    const auto mesh=triangulateRankTwoFinitePartition(a,p,incidence,rotation);
    std::set<std::size_t> used;for(const auto& t:mesh.triangles())used.insert(t.points.begin(),t.points.end());
    require(used==std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end()),"scoped triangulation lost actual points");
    require(mesh.regionSummaries().size()==p.regionIndices().size(),"scoped triangulation lost actual cells");
    rejects([&]{certifyRankTwoRegistryRotation(a,r,incidence);},"INPUT_IDENTITY_MISMATCH","scoped certificate promoted to the full registry");
    rejects([&]{certifyRankTwoChordIncidenceSweep(global,p);},"INPUT_IDENTITY_MISMATCH","foreign arrangement accepted");
    const auto otherRegistry=registry(a,2);const auto other=buildRankTwoFinitePartition(a,otherRegistry);
    rejects([&]{certifyRankTwoRegistryRotation(a,other,incidence);},"INPUT_IDENTITY_MISMATCH","stale partition incidence accepted");
    rejects([&]{triangulateRankTwoFinitePartition(a,other,incidence,rotation);},"INPUT_IDENTITY_MISMATCH","foreign mesh scope accepted");
    rejects([&]{certifyRankTwoChordIncidenceSweep(a,p,{0});},"RESOURCE_LIMIT","empty scoped predicate budget granted incidence");
    rejects([&]{certifyRankTwoRegistryRotation(a,p,incidence,{0});},"RESOURCE_LIMIT","empty scoped rotation budget granted rotation");
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-mesh-scope-crossing-v1");
    const auto b=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    const auto coarse=registry(b,2);const auto allFinite=buildRankTwoFinitePartition(b,coarse);
    rejects([&]{certifyRankTwoChordIncidenceSweep(b,allFinite);},"ARRANGEMENT_UNRESOLVED","actual in-domain crossing hidden by scope");
    rejects([&]{certifyRankTwoChordIncidenceSweep(b,coarse);},"ARRANGEMENT_UNRESOLVED","legacy full-registry crossing rejection changed");
    std::ofstream out(path);require(bool(out),"cannot save partition scope evidence");
    out<<"rank-two-finite-mesh-scope-test-v1\nPARTITION "<<p.scientificHash()<<'\n'<<p.canonicalRecord()
       <<"INCIDENCE "<<incidence.scientificHash()<<'\n'<<incidence.canonicalRecord()
       <<"ROTATION "<<rotation.scientificHash()<<'\n'<<rotation.canonicalRecord()
       <<"TRIANGULATION "<<mesh.scientificHash()<<'\n'<<mesh.canonicalRecord();
    out.flush();require(bool(out),"incomplete scope evidence");
    std::cout<<"PASS immutable partition incidence/rotation/mesh; retained crossings and stale/full promotion rejected; points="<<n<<" edges="<<m<<'\n';
}
}
#endif
int main(int argc,char** argv){try{
#if defined(RANK_TWO_FINITE_PARTITION_INCIDENCE_V1)
    if(argc!=2)throw std::invalid_argument("scope evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: incidence and native rotation bound to immutable finite partition");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
