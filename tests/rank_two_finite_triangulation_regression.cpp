#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoFiniteTriangulation.h")
#include "spin/detail/RankTwoFiniteTriangulation.h"
#include "spin/detail/RankTwoChordSweep.h"
#include "spin/detail/RankTwoDiagnostics.h"
#if __has_include("spin/detail/RankTwoObservedCartesian.h")
#include "spin/detail/RankTwoObservedCartesian.h"
#define RANK_TWO_PREDICATE_OBSERVER_AVAILABLE 1
#endif
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
template<class F>auto timed(const std::string& label,F action){
    const auto begin=std::chrono::steady_clock::now();std::cerr<<"STAGE_BEGIN "<<label<<'\n';
    auto value=action();const auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-begin).count();
    std::cerr<<"STAGE_END "<<label<<" elapsed_ns="<<ns<<'\n';return value;
}
spin::RankTwoArcRegistry registry(const spin::RankTwoArrangement& a,int denominator=2){
    spin::RankTwoArcSamplingPlan plan;for(const auto& arc:spin::arrangementRecords(a).arcs)
        if(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite){
            std::vector<Q> fractions;for(int i=1;i<denominator;++i)fractions.emplace_back(i,denominator);
            plan.emplace(arc.identity,std::move(fractions));
        }
    return spin::buildRankTwoArcRegistry(a,plan);
}
E det(const spin::RankTwoArcRegistry& r,const std::array<std::size_t,3>& p){
    const auto& a=r.points().at(p[0]).coordinates;const auto& b=r.points().at(p[1]).coordinates;const auto& c=r.points().at(p[2]).coordinates;
    return (b[0].expression()-a[0].expression())*(c[1].expression()-a[1].expression())
        -(b[1].expression()-a[1].expression())*(c[0].expression()-a[0].expression());
}
void check(const spin::RankTwoArrangement& a,const char* name,std::ostream& saved,bool annulus,const char* path){
    // Frozen *polygonal* method resolution, not a curve error/isotopy grant.
    // For a circle of squared radius s and an x-affine n-way arc subdivision,
    // every chord has squared radial distance >= s*(1-1/n). The annulus has
    // successive squared radii (3-sqrt(2))/2, 2, (3+sqrt(2))/2, 9/4.
    // n=128 satisfies all strict nested-circle separation inequalities before
    // CT runs. In particular 9/4*(127/128) > (3+sqrt(2))/2. It is not chosen
    // by retrying CDT or changing the source coefficients.
    using namespace spin;const std::string label=name;
    const auto r=timed(label+"/registry",[&]{return registry(a,annulus?128:2);});
    const auto p=timed(label+"/partition",[&]{return buildRankTwoFinitePartition(a,r);});
    const auto incidence=timed(label+"/incidence",[&]{return certifyRankTwoChordIncidenceSweep(a,r);});
    const auto rotation=timed(label+"/rotation",[&]{return certifyRankTwoRegistryRotation(a,r,incidence);});
    // Incomplete diagnostic antecedents cannot be mistaken for the final mesh
    // artifact, which is only written after the WHOLE test succeeds.
    if(annulus){std::ofstream checkpoint(std::string(path)+".annulus-pretriangulation.txt");require(bool(checkpoint),"cannot save diagnostic antecedents");
        checkpoint<<"rank-two-finite-pretriangulation-diagnostic-v1\nINCOMPLETE NOT_ACCEPTED\nPARTITION "<<p.scientificHash()<<'\n'<<p.canonicalRecord()
            <<"INCIDENCE "<<incidence.scientificHash()<<'\n'<<incidence.canonicalRecord()<<"ROTATION "<<rotation.scientificHash()<<'\n'<<rotation.canonicalRecord();
        checkpoint.flush();require(bool(checkpoint),"incomplete diagnostic antecedent write");}
    const auto result=timed(label+"/triangulation",[&]{return triangulateRankTwoFinitePartition(a,p,incidence,rotation);});
    require(!result.triangles().empty() && result.regionSummaries().size()==p.regionIndices().size(),"finite region triangulation missing");
    require(result.partitionIdentity()==p.scientificHash() && result.registryIdentity()==r.scientificHash(),"triangulation authority changed");
    std::set<std::size_t> used;std::map<std::size_t,std::size_t> counts;std::set<std::array<std::size_t,3>> unique;
    for(const auto& t:result.triangles()){
        require(det(r,t.points).sign()>0,"zero/reversed actual finite triangle");
        require(std::find(p.regionIndices().begin(),p.regionIndices().end(),t.registryRegionIndex)!=p.regionIndices().end(),"triangle labelled by excluded/parity cell");
        require(unique.insert(t.points).second,"duplicate finite triangle");
        for(auto point:t.points)require(std::find(p.pointIndices().begin(),p.pointIndices().end(),point)!=p.pointIndices().end(),"anonymous library point accepted");
        used.insert(t.points.begin(),t.points.end());++counts[t.registryRegionIndex];
    }
    require(used==std::set<std::size_t>(p.pointIndices().begin(),p.pointIndices().end()),"actual finite feature omitted");
    std::size_t computedHoles=0,empty=0,full=0,cut=0;
    for(const auto& s:result.regionSummaries()){
        const auto& cell=r.regions().at(s.registryRegionIndex);
        require(s.triangleCount==counts.at(s.registryRegionIndex) && s.connectedComponents==cell.outerLoops.size(),"triangle ownership/component coverage wrong");
        require(s.holes==cell.innerLoops.size() && s.eulerCharacteristic==static_cast<std::int64_t>(s.connectedComponents)-static_cast<std::int64_t>(s.holes),"actual triangle complex lost a native hole");
        computedHoles+=s.holes;empty+=cell.kind==RankTwoRegionKind::Empty;full+=cell.kind==RankTwoRegionKind::Full;cut+=cell.kind==RankTwoRegionKind::Cut;
    }
    if(annulus)require(computedHoles==3 && empty==1 && full==1 && cut==2,"annulus region labels replaced by parity filling");
    else require(computedHoles==0 && full==1 && result.triangles().size()==2,"full disk finite triangulation incorrect");
    // Hashes/labels are never caller-entered completion. Reload must recreate
    // the entire same sorted mesh and reject missing/reversed/foreign data.
    require(timed(label+"/reload",[&]{return reloadRankTwoFiniteTriangulation(a,p,incidence,rotation,result.canonicalRecord());}).scientificHash()==result.scientificHash(),"reload changed finite mesh");
    auto bad=result.canonicalRecord();const auto at=bad.find("TRIANGLE ");require(at!=std::string::npos,"actual triangle record absent");
    bad.erase(at,bad.find('\n',at)-at+1);
    rejects([&]{reloadRankTwoFiniteTriangulation(a,p,incidence,rotation,bad);},"INPUT_IDENTITY_MISMATCH","omitted triangle accepted");
    RankTwoFiniteTriangulationSettings budget;budget.maximumTriangles=0;
    rejects([&]{triangulateRankTwoFinitePartition(a,p,incidence,rotation,budget);},"RESOURCE_LIMIT","zero triangle budget granted mesh");
    saved<<"CASE "<<name<<" TRIANGULATION "<<result.scientificHash()<<'\n'<<result.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" triangles="<<result.triangles().size()<<" computed_holes="<<computedHoles
        <<" labels Empty/Cut/Full="<<empty<<'/'<<cut<<'/'<<full<<'\n';
}
void tests(const char* path){
    using namespace spin;std::ostringstream saved;saved<<"rank-two-finite-triangulation-test-v1\n";
#if defined(RANK_TWO_PREDICATE_OBSERVER_AVAILABLE)
    {
        using K=CGAL::Cartesian<CORE::Expr>;using P=K::Point_2;
        detail::RankTwoPredicateObservation state;std::vector<std::string> observations;
        const auto sink=[](const char* message,void* context)noexcept{
            static_cast<std::vector<std::string>*>(context)->emplace_back(message);
        };
        detail::ScopedRankTwoDiagnostics nested(sink,&observations);
        const P p(E(0),E(0)),q(E(2),E(0)),left(E(0),E(2)),right(E(0),E(-2)),middle(E(1),E(0));
        state.bind(p,10);state.bind(q,11);state.bind(left,12);state.bind(right,13);state.bind(middle,14);
        state.constraint=7;detail::RankTwoObservedCartesian traits(&state);
        const auto observed=traits.orientation_2_object();const auto original=K().orientation_2_object();
        for(const P& third:{left,right,middle})require(observed(P(p),P(q),P(third))==original(p,q,third),"predicate observer changed exact orientation");
        require(observations.size()==6,"predicate observer lost begin/end pairs");
        require(observations.front()=="finite-cdt:predicate-begin ORIENTATION2 7 10 11 12","predicate observer lost actual coordinate identities");
        require(observations.back()=="finite-cdt:predicate-end ORIENTATION2 7 10 11 14 0","zero orientation observation changed");
        state.constraint.reset();observed(p,q,left);require(observations.size()==6,"observer logged outside selected recovery scope");
    }
    std::cout<<"PASS observation forwards original Cartesian CORE predicates\n";
#else
    throw std::runtime_error("MISSING_FEATURE: read-only exact constraint predicate observation");
#endif
    const auto observe=[](const char* message,void*)noexcept{if(std::string(message).rfind("finite-cdt:",0)==0)std::cerr<<"CDT_OBSERVATION "<<message<<'\n';};
    detail::ScopedRankTwoDiagnostics observation(observe,nullptr);
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"finite-triangulation-full-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(disk,Q(1));check(a,"full-disk",saved,false,path);
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"finite-triangulation-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));
    const auto coarse=registry(b);
    rejects([&]{certifyRankTwoChordIncidenceSweep(b,coarse);},"ARRANGEMENT_UNRESOLVED","crossing coarse annulus entered No_intersection_tag CT");
    const E root2=CORE::sqrt(E(2));const std::array<E,4> radii{(E(3)-root2)/2,E(2),(E(3)+root2)/2,E(9)/4};
    for(std::size_t i=1;i<radii.size();++i)require((radii[i]*E(127)/128-radii[i-1]).sign()>0,"frozen nested-circle chord separation not proven");
    check(b,"annulus",saved,true,path);
    const auto r=registry(a),foreign=timed("foreign/registry",[&]{return registry(b,128);});const auto p=buildRankTwoFinitePartition(a,r);
    const auto good=certifyRankTwoChordIncidenceSweep(a,r),other=certifyRankTwoChordIncidenceSweep(b,foreign);
    const auto correct=certifyRankTwoRegistryRotation(a,r,good),wrong=certifyRankTwoRegistryRotation(b,foreign,other);
    rejects([&]{triangulateRankTwoFinitePartition(a,p,other,correct);},"INPUT_IDENTITY_MISMATCH","foreign incidence granted triangulation");
    rejects([&]{triangulateRankTwoFinitePartition(a,p,good,wrong);},"INPUT_IDENTITY_MISMATCH","foreign rotation granted triangulation");
    std::ofstream out(path);require(bool(out),"cannot write triangulation evidence");out<<saved.str();out.flush();require(bool(out),"incomplete triangulation evidence");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoFiniteTriangulation.h")
    if(argc!=2)throw std::invalid_argument("finite triangulation evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: hole-preserving finite CDT with native region ownership and no anonymous points");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
