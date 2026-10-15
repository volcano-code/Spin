#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoMixedArcIntervals.h")
#include "spin/detail/RankTwoMixedArcIntervals.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "spin/detail/RankTwoDiagnostics.h"
#if __has_include("spin/detail/RankTwoArcMembership.h")
#include "spin/detail/RankTwoArcMembership.h"
#endif
#include <algorithm>
#include <chrono>
#include <set>
namespace {
using Q=CGAL::Gmpq;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"endpoint-neighborhood-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(source,global,cone);
    std::string event;RankTwoArcSamplingPlan plan;
    for(const auto& v:arrangementRecords(a).vertices)
        if(C::reload(v.coordinates[0]).expression()==1 && C::reload(v.coordinates[1]).expression()==0)event=v.identity;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    RankTwoEndpointNeighborhoodSettings ns;ns.maximumChordDistance=Q(1,32);ns.heightMargin=Q(1,8);
    const auto n=certifyRankTwoEndpointNeighborhood(a,event,ns);
    const auto initial=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n});
    // Freeze ALL actual eight finite incident arcs, from the shared event to
    // each original native midpoint. No selection based on solver outcomes.
    std::vector<RankTwoMixedArcInterval> intervals;
    for(const auto& arc:initial.arcs()){
        const bool first=arc.ends[0].vertexIdentity==event,last=arc.ends[1].vertexIdentity==event;
        if(!first && !last)continue;
        require(arc.finite && first!=last,"frozen star needs exact finite arc bounds");
        const RankTwoRegistryPoint* middle=nullptr;
        for(const auto& p:initial.points())if(p.kind==RankTwoRegistryPointKind::ArcSample && p.arcIdentity==arc.identity){
            require(middle==nullptr,"unexpected multiple original midpoint samples");middle=&p;
        }
        require(middle!=nullptr,"incident native midpoint absent");
        intervals.push_back({arc.identity,first?std::array<std::string,2>{event,middle->identity}
            :std::array<std::string,2>{middle->identity,event}});
    }
    require(intervals.size()==8,"incomplete frozen native star");
#if __has_include("spin/detail/RankTwoArcMembership.h")
    for(const auto& i:intervals)for(const auto& p:initial.points())
        if(p.kind==RankTwoRegistryPointKind::ArcSample && p.arcIdentity==i.arcIdentity)
            require(detail::isOnRankTwoNativeArc(a,i.arcIdentity,p.coordinates),"actual native midpoint rejected");
    detail::AlgebraicKernel kernel;
    const auto rational=[&](int value){
        std::vector<detail::AlgebraicReal> roots;
        const std::vector<CORE::BigInt> coefficients{CORE::BigInt(-value),CORE::BigInt(1)};
        kernel.solve_1_object()(detail::AlgebraicPolynomial(coefficients.begin(),coefficients.end()),
            false,std::back_inserter(roots));
        require(roots.size()==1,"rational fixture ambiguity");return C::fromRoot(roots[0]);
    };
    const std::string upperH="7b6b470063cd4201ee24ed859be307244295edfd3da11bf4ab3057bb5dc41df1";
    require(detail::isOnRankTwoNativeArc(a,upperH,{rational(1),rational(0)}),"actual common event rejected");
    require(!detail::isOnRankTwoNativeArc(a,upperH,{rational(1),rational(-1)}),
        "another valid h zero at same event abscissa accepted as native endpoint");
    require(!detail::isOnRankTwoNativeArc(a,upperH,{rational(2),rational(0)}),"point outside native arc accepted");
    std::vector<detail::AlgebraicReal> ys;
    kernel.solve_1_object()(detail::AlgebraicPolynomial(-1,1,1),false,std::back_inserter(ys));
    require(ys.size()==2,"frozen h fiber root count");
    require(detail::isOnRankTwoNativeArc(a,upperH,{rational(0),C::fromRoot(ys[1])}),"actual h start event rejected");
    require(!detail::isOnRankTwoNativeArc(a,upperH,{rational(0),C::fromRoot(ys[0])}),"opposite actual h branch accepted");
#else
    throw std::runtime_error("MISSING_FEATURE: exact native fiber membership without Cartesian product elimination");
#endif
    detail::rankTwoDiagnostic("mixed-test:initial-eight-intervals-ready");
    const auto result=refineRankTwoMixedArcIntervals(a,plan,{n},intervals,Q(1,32));
    detail::rankTwoDiagnostic("mixed-test:first-result-ready");
    const auto& r=result.registry();
    require(result.intervals().size()==8 && result.rounds().size()>1 && !r.subdivisions().empty(),
        "eight intervals were merely initialized without actual common refinement");
    require(result.rounds().front().selectedEdges==16 && result.rounds().front().splitEdges>0,
        "frozen coarse star not checked in full");
    std::set<std::string> expected,covered,changedArcs;
    for(const auto& interval:result.intervals()){
        const auto arc=std::find_if(r.arcs().begin(),r.arcs().end(),[&](const auto& x){return x.identity==interval.arcIdentity;});
        require(arc!=r.arcs().end(),"selected arc lost");
        bool inside=false,finished=false;
        for(const auto index:arc->edges){const auto& e=r.edges().at(index);
            const auto& left=r.points().at(e.points[0]);const auto& right=r.points().at(e.points[1]);
            if(left.identity==interval.pointIdentities[0])inside=true;
            if(inside){require(expected.emplace(e.identity).second,"overlapping requested interval silently coalesced");}
            if(inside && right.identity==interval.pointIdentities[1]){finished=true;break;}
        }
        require(finished,"actual final interval has an unproved gap");
    }
    for(const auto& id:result.endpointEdges())require(expected.count(id) && covered.emplace(id).second,"invalid endpoint coverage");
    require(result.endpointEdges().size()==8,"common endpoint star partially lost");
    for(const auto& p:result.proofs()){
        require(p.registryIdentity()==r.scientificHash() && p.arrangementIdentity()==a.scientificHash()
            && expected.count(p.edgeIdentity()) && covered.emplace(p.edgeIdentity()).second,
            "one interval retained a proof from an earlier independent registry");
        require(p.strip().chordDistanceUpperBound<=Q(1,32),"mixed error relaxed");
    }
    require(expected==covered && result.rounds().back().selectedEdges==expected.size()
        && result.rounds().back().splitEdges==0,"common final partition not completely covered");
    for(const auto& p:r.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision){
        const auto arc=std::find_if(r.arcs().begin(),r.arcs().end(),[&](const auto& z){return z.identity==p.arcIdentity;});
        require(arc!=r.arcs().end(),"subdivision native owner lost");
        bool checked=false;
        for(const auto& owner:arc->owners){
            if(owner.polynomial=="h" || owner.polynomial=="D" || owner.polynomial=="F"){
                checked=true;changedArcs.insert(p.arcIdentity);
            }
            // The two frozen straight normal-boundary intervals also require
            // refinement when their boxes must exclude the tangent curves.
            if(owner.polynomial=="NORMAL_EDGE:0:1")checked=true;
        }
        require(checked,"undeclared subdivision support");
        require(detail::isOnRankTwoNativeArc(a,p.arcIdentity,p.coordinates),
            "subdivision left the original exact native fiber/branch");
        // The independent saved-evidence audit checks every literal h/D/F
        // polynomial and its correlated roots via rational resultants/Sturm.
        // Do not rebuild Cartesian CORE expression products here: doing so
        // repeats elimination for hundreds of already certified native points.
    }
    detail::rankTwoDiagnostic("mixed-test:all-native-points-ready");
    require(changedArcs.size()>1,"multiple actual curved arcs were not refined together");
    for(const auto& h:result.rounds())require(h.selectedEdges==h.endpointEdges+h.regularEdges+h.splitEdges
        && h.endpointEdges==8,"round dropped another interval's failed edge");
    auto reverse=intervals;std::reverse(reverse.begin(),reverse.end());
    const auto reordered=refineRankTwoMixedArcIntervals(a,plan,{n},reverse,Q(1,32));
    require(reordered.canonicalRecord()==result.canonicalRecord(),"caller interval order changed scientific mapping");
    auto duplicate=intervals;duplicate.push_back(intervals.front());
    rejects([&]{refineRankTwoMixedArcIntervals(a,plan,{n},duplicate,Q(1,32));},"INPUT_IDENTITY_MISMATCH","duplicate coverage accepted");
    auto overlap=intervals;const auto& i=intervals.front();
    // A strict subinterval shares a real selected edge, not just a vertex.
    const auto gate=std::find_if(initial.points().begin(),initial.points().end(),[&](const auto& p){return p.arcIdentity==i.arcIdentity
        && !p.endpointNeighborhoodIdentity.empty();});
    require(gate!=initial.points().end(),"common gate absent");
    overlap.push_back({i.arcIdentity,i.pointIdentities[0]==event?std::array<std::string,2>{event,gate->identity}
        :std::array<std::string,2>{gate->identity,event}});
    rejects([&]{refineRankTwoMixedArcIntervals(a,plan,{n},overlap,Q(1,32));},"INPUT_IDENTITY_MISMATCH","overlapping interval accepted");
    rejects([&]{refineRankTwoMixedArcIntervals(a,plan,{n},{},Q(1,32));},"INPUT_IDENTITY_MISMATCH","empty coverage accepted");
    auto backwards=intervals;std::swap(backwards.front().pointIdentities[0],backwards.front().pointIdentities[1]);
    rejects([&]{refineRankTwoMixedArcIntervals(a,plan,{n},backwards,Q(1,32));},"INPUT_IDENTITY_MISMATCH","reversed interval accepted");
    RankTwoSourceLimits zero;zero.maximumRefinementDepth=0;
    rejects([&]{refineRankTwoMixedArcIntervals(a,plan,{n},intervals,Q(1,32),{},zero);},"RESOURCE_LIMIT","unproved shared star granted success");
    std::ofstream out(path);require(bool(out),"cannot write mixed intervals evidence");
    out<<"rank-two-mixed-intervals-test-v1\nMIXED_INTERVALS "<<result.scientificHash()<<'\n'<<result.canonicalRecord();
    out.flush();require(bool(out),"incomplete mixed intervals evidence");
    std::cout<<"PASS eight actual native intervals, one common final registry, all proof rebuilds, multiple curved arcs and canonical order; global isotopy NOT_GRANTED\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoMixedArcIntervals.h")
    if(argc!=2)throw std::invalid_argument("mixed intervals evidence path required");
    const auto start=std::chrono::steady_clock::now();
    const auto observe=+[](const char* phase,void* context)noexcept{
        try{const std::string p(phase);
            if(p.rfind("mixed-",0)!=0 && p.rfind("regular-edge:",0)!=0)return;
            const auto& origin=*static_cast<const std::chrono::steady_clock::time_point*>(context);
            std::cerr<<"DIAGNOSTIC_NS "<<std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-origin).count()<<' '<<p<<'\n';
        }catch(...){}
    };
    spin::detail::ScopedRankTwoDiagnostics observer(observe,const_cast<std::chrono::steady_clock::time_point*>(&start));
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: simultaneous mixed native intervals in one common final registry");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
