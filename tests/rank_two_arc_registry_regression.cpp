#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcRegistry.h")
#include "spin/detail/RankTwoArcRegistry.h"
#include "spin/detail/RankTwoDiagnostics.h"
#if __has_include("spin/detail/RankTwoArcRegistryFans.h")
#include "spin/detail/RankTwoArcRegistryFans.h"
#endif
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <map>
#include <sstream>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=spin::detail::CertifiedAlgebraicReal;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
bool onOriginal(const std::array<Q,25>& coefficients,const std::array<C,2>& coordinates){
    using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
    using P=Traits::Polynomial_2;
    CGAL::Gmpz denominator(1);
    for(const auto& c:coefficients)mpz_lcm(denominator.mpz(),denominator.mpz(),c.denominator().mpz());
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);P p(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(coefficients[5*i+j]!=0){
        const auto& c=coefficients[5*i+j];std::ostringstream integer;
        integer<<c.numerator()*(denominator/c.denominator());
        p+=P(CORE::BigInt(integer.str()))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    auto& k=Traits::CKvA_2::instance().kernel();
    const auto point=k.construct_algebraic_real_2_object()(
        spin::detail::rankTwoAlgebraicHandle(coordinates[0]),spin::detail::rankTwoAlgebraicHandle(coordinates[1]));
    return k.is_zero_at_2_object()(CGAL::make_square_free(p),point);
}
void check(const spin::RankTwoArrangement& authority,std::ostream& out,bool hasHoles,bool hasInfinity,bool testOrder=false){
    using namespace spin;
    const auto& native=arrangementRecords(authority);
    RankTwoArcSamplingPlan plan;
    std::size_t finite=0;
    for(const auto& arc:native.arcs){
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite
            && arc.ends[1].location==RankTwoArcEndLocation::Finite){
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});++finite;
        }
    }
    const auto registry=buildRankTwoArcRegistry(authority,plan);
    require(registry.arrangementIdentity()==authority.scientificHash(),"registry authority lost");
    require(registry.points().size()==native.vertices.size()+finite
        && registry.edges().size()==2*finite && registry.arcs().size()==native.arcs.size()
        && registry.regions().size()==native.regions.size(),"shared feature omitted or duplicated");
    std::map<std::string,std::size_t> events;
    for(std::size_t i=0;i<registry.points().size();++i){
        const auto& p=registry.points()[i];
        if(p.kind==RankTwoRegistryPointKind::Event){
            require(events.emplace(p.eventIdentity,i).second,"duplicate event seam");
            const RankTwoVertexRecord* original=nullptr;
            for(const auto& v:native.vertices)if(v.identity==p.eventIdentity)original=&v;
            require(original && p.coordinates[0].expression()==C::reload(original->coordinates[0]).expression()
                && p.coordinates[1].expression()==C::reload(original->coordinates[1]).expression(),
                "shared event coordinate differs from original authority");
        }else{
            require(p.fraction==Q(1,2),"moving fraction lost");
            const RankTwoArcRecord* parent=nullptr;
            for(const auto& a:native.arcs)if(a.identity==p.arcIdentity)parent=&a;
            require(parent!=nullptr && p.endpointIdentities[0]==parent->ends[0].vertexIdentity
                && p.endpointIdentities[1]==parent->ends[1].vertexIdentity,"moving point DAG changed");
            const int axis=parent->vertical?1:0;
            const auto endpoint=[&](int k){
                for(const auto& v:native.vertices)if(v.identity==*parent->ends[k].vertexIdentity)
                    return C::reload(v.coordinates[axis]).expression();
                throw std::runtime_error("native endpoint missing");
            };
            require(p.coordinates[axis].expression()==(endpoint(0)+endpoint(1))/2,
                "registry substituted a root-interval midpoint");
            require(p.coordinates[1].expression().sign()==C::reload(parent->interiorWitness[1]).expression().sign(),
                "registered point switched annulus upper/lower branch");
            for(const auto& owner:parent->owners){
                const int material=owner.polynomial=="F"?0:owner.polynomial=="h"?1:owner.polynomial=="D"?2:-1;
                if(material>=0)require(onOriginal(native.originalFhD[material],p.coordinates),
                    "registered coordinate left its original material curve");
                else{
                    const int normal=owner.polynomial=="NORMAL_EDGE:0:1"?2:
                        owner.polynomial=="NORMAL_EDGE:0:2"?1:owner.polynomial=="NORMAL_EDGE:1:2"?0:-1;
                    require(normal>=0 && std::size_t(normal)<native.normalHalfplanes.size()
                        && onOriginal(native.normalHalfplanes[normal],p.coordinates),
                        "registered normal coordinate left its original plane");
                }
            }
        }
    }
    std::vector<unsigned> left(registry.edges().size()),right(registry.edges().size());
    std::size_t holes=0,infinity=0;
    for(const auto& region:registry.regions()){
        const RankTwoRegionRecord* original=nullptr;
        for(const auto& r:native.regions)if(r.identity==region.identity)original=&r;
        require(original!=nullptr && original->kind==region.kind
            && original->normalSigns==region.normalSigns && original->originalSigns==region.originalSigns
            && original->outerLoops.size()==region.outerLoops.size()
            && original->innerLoops.size()==region.innerLoops.size(),"DCEL region/hole labels lost");
        holes+=region.innerLoops.size();
        for(bool inner:{false,true}){
            const auto& loops=inner?region.innerLoops:region.outerLoops;
            const auto& sourceLoops=inner?original->innerLoops:original->outerLoops;
            for(std::size_t l=0;l<loops.size();++l){
                const auto& loop=loops[l];
                require(loop.sourceArcs.size()==sourceLoops[l].size(),"hole arc silently omitted");
                for(std::size_t j=0;j<loop.sourceArcs.size();++j)
                    require(loop.sourceArcs[j].identity==sourceLoops[l][j].identity
                        && loop.sourceArcs[j].reversed==sourceLoops[l][j].reversed
                        && loop.sourceArcs[j].atInfinity==sourceLoops[l][j].atInfinity,
                        "original directed ring reordered or reversed");
                if(!loop.finite){
                    ++infinity;
                    rejects([&]{registry.loopVertices(region.identity,inner,l);},
                        "ARRANGEMENT_UNRESOLVED","infinite ring truncated into a closed polygon");
                    continue;
                }
                const auto points=registry.loopVertices(region.identity,inner,l);
                require(points.size()==loop.edges.size() && points.size()>=4,
                    "shared midpoint not consumed by adjacent ring");
                for(std::size_t j=0;j<loop.edges.size();++j){
                    const auto ref=loop.edges[j];
                    const auto& edge=registry.edges().at(ref.edgeIndex);
                    const auto begin=edge.points[ref.reversed?1:0],end=edge.points[ref.reversed?0:1];
                    require(points[j]==begin && points[(j+1)%points.size()]==end,
                        "ring has a gap, duplicate seam or independently sampled edge");
                    const RankTwoArcRecord* arc=nullptr;
                    for(const auto& a:native.arcs)if(a.identity==edge.arcIdentity)arc=&a;
                    require(arc && (ref.reversed?arc->rightRegion:arc->leftRegion)==region.identity,
                        "native side of region loop changed");
                    ++(ref.reversed?right:left)[ref.edgeIndex];
                }
            }
        }
    }
    if((holes>0)!=hasHoles || (infinity>0)!=hasInfinity)
        std::cerr<<"COVERAGE authority="<<authority.scientificHash()<<" vertices="<<native.vertices.size()
            <<" finite_arcs="<<finite<<" holes="<<holes<<" infinite_loops="<<infinity
            <<" expected_holes="<<hasHoles<<" expected_infinity="<<hasInfinity<<'\n'
            <<registry.canonicalRecord();
    require((holes>0)==hasHoles && (infinity>0)==hasInfinity,"missing hole/infinity fixture coverage");
    if(finite==native.arcs.size())for(std::size_t i=0;i<left.size();++i)
        require(left[i]==1 && right[i]==1,"two adjacent faces did not consume one shared reversed edge");
    auto reversedPlan=plan;
    for(auto& entry:reversedPlan)entry.second={Q(3,4),Q(1,4),Q(1,2)};
    auto sortedPlan=plan;
    for(auto& entry:sortedPlan)entry.second={Q(1,4),Q(1,2),Q(3,4)};
    // A second order-independent plan is checked on the simple disk only;
    // annulus/local registry still has a real moving midpoint on every arc.
    if(testOrder){
        const auto a=buildRankTwoArcRegistry(authority,reversedPlan);
        const auto b=buildRankTwoArcRegistry(authority,sortedPlan);
        require(a.scientificHash()==b.scientificHash() && a.canonicalRecord()==b.canonicalRecord(),
            "plan input order changed scientific registry or published indices");
    }
    auto missing=plan;missing.erase(missing.begin());
    rejects([&]{buildRankTwoArcRegistry(authority,missing);},"INPUT_IDENTITY_MISMATCH",
        "caller omitted a finite shared arc");
    auto duplicate=plan;duplicate.begin()->second.push_back(Q(1,2));
    rejects([&]{buildRankTwoArcRegistry(authority,duplicate);},"INPUT_IDENTITY_MISMATCH",
        "duplicate moving point manufactured a zero edge");
    auto invalid=plan;invalid.begin()->second.push_back(Q(0));
    rejects([&]{buildRankTwoArcRegistry(authority,invalid);},"INPUT_IDENTITY_MISMATCH",
        "endpoint alias accepted as an interior sample");
    RankTwoSourceLimits tiny;tiny.maximumVertices=native.vertices.size()+finite-1;
    rejects([&]{buildRankTwoArcRegistry(authority,plan,tiny);},"RESOURCE_LIMIT",
        "registry silently exceeded vertex budget");
    rejects([&]{registry.loopVertices("not-a-region",false,0);},"INPUT_IDENTITY_MISMATCH",
        "foreign region consumed local shared edges");
    out<<"REGISTRY "<<registry.scientificHash()<<' '<<native.vertices.size()<<' '
        <<finite<<' '<<holes<<' '<<infinity<<'\n'<<registry.canonicalRecord();
}
void tests(const char* path){
    using namespace spin;
    std::ofstream out(path);require(bool(out),"cannot open registry evidence");
    out<<"rank-two-shared-arc-registry-test-v1\n";
    const auto disk=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"shared-registry-disk-v1");
    // The unbounded exterior has a fictitious outer ring even though every
    // real circle arc is finite. Its finite INNER ring is still shared with
    // the disk; never delete the fictitious ring to fit a polygon-only API.
    check(analyzeRankTwoGlobalArrangement(disk,Q(1)),out,true,true,true);
    const auto annulus=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-registry-annulus-v1");
    const auto global=analyzeRankTwoGlobalArrangement(annulus,Q(3,2));
    check(global,out,true,true);
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    check(analyzeRankTwoNormalFace(annulus,global,cone),out,false,true);
#if __has_include("spin/detail/RankTwoArcRegistryFans.h")
    const auto equator=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"shared-registry-fans-v1");
    detail::rankTwoDiagnostic("registry-fixture:equator-arrangement");
    const auto full=analyzeRankTwoGlobalArrangement(equator,Q(1));
    RankTwoNormalFace normal;normal.vertexIds={0,1,2};normal.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto overlay=analyzeRankTwoNormalFace(equator,full,normal);
    const auto& records=arrangementRecords(overlay);
    std::string start,negativeAnchor,positiveAnchor;
    RankTwoArcSamplingPlan plan;
    for(const auto& v:records.vertices){
        const E x=C::reload(v.coordinates[0]).expression(),y=C::reload(v.coordinates[1]).expression();
        if(x==1 && y==0)start=v.identity;
        if(x>1 && y==-E(1)/2)negativeAnchor=v.identity;
        if(x==0 && y==1)positiveAnchor=v.identity;
    }
    for(const auto& arc:records.arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite
        && arc.ends[1].location==RankTwoArcEndLocation::Finite)plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    require(!start.empty() && !negativeAnchor.empty() && !positiveAnchor.empty(),"fan registry events absent");
    const auto negative=RankTwoMovingGate::fromVertices(overlay,start,negativeAnchor,1,Q(1,512));
    const auto positive=RankTwoMovingGate::fromVertices(overlay,start,positiveAnchor,1,Q(1,1024));
    const std::vector<RankTwoEndpointFan> fans{
        certifyRankTwoEndpointFan(overlay,negative,{Q(-1,1024),Q(0),Q(99,100),Q(101,100)},Q(1,1048576)),
        certifyRankTwoEndpointFan(overlay,positive,{Q(0),Q(1,1024),Q(99,100),Q(101,100)},Q(1,1048576))};
    detail::rankTwoDiagnostic("registry-fixture:fan-publication");
    const auto shared=buildRankTwoArcRegistry(overlay,plan,fans);
    require(shared.endpointFans().size()==2,"whole incident proof group omitted");
    std::size_t gates=0,covered=0;
    for(const auto& point:shared.points())if(point.kind==RankTwoRegistryPointKind::SharedGate){
        ++gates;require((point.gateIdentity==negative.identity() || point.gateIdentity==positive.identity())
            && point.endpointIdentities[0]==start && point.projectionAxis==1
            && point.endpointIdentities[1]==(point.gateIdentity==negative.identity()?negativeAnchor:positiveAnchor)
            && point.fraction==(point.gateIdentity==negative.identity()?Q(1,512):Q(1,1024))
            && point.coordinates[1].expression()==(point.gateIdentity==negative.identity()?-E(1)/1024:E(1)/1024),
            "shared gate point lost its moving DAG or exact coordinates");
    }
    for(const auto& edge:shared.edges())if(!edge.endpointFanIdentity.empty()){
        ++covered;const RankTwoEndpointFan* fan=nullptr;
        for(const auto& f:shared.endpointFans())if(f.scientificHash()==edge.endpointFanIdentity)fan=&f;
        require(fan!=nullptr,"edge attached a nonexistent or stale fan proof");
        const RankTwoSharedGateArcCertificate* proof=nullptr;
        for(const auto& a:fan->orderedArcs())if(a.arcIdentity==edge.arcIdentity)proof=&a;
        require(proof!=nullptr,"covered shared edge not in complete fan");
        const auto& a=shared.points()[edge.points[0]],&b=shared.points()[edge.points[1]];
        const auto& event=a.kind==RankTwoRegistryPointKind::Event?a:b;
        const auto& gate=a.kind==RankTwoRegistryPointKind::SharedGate?a:b;
        require(event.eventIdentity==proof->endpointIdentity && gate.gateIdentity==proof->gateIdentity
            && gate.coordinates[0].expression()==proof->gateCoordinates[0].expression()
            && gate.coordinates[1].expression()==proof->gateCoordinates[1].expression(),
            "proof attached to a different published chord");
    }
    require(gates==6 && covered==6,"one fan arc was independently regenerated or omitted");
    auto split=plan;std::string d;
    for(const auto& arc:records.arcs)if(arc.ends[1].vertexIdentity==start)
        for(const auto& owner:arc.owners)if(owner.polynomial=="D")d=arc.identity;
    require(!d.empty(),"missing equator D arc for fan-splitting counterexample");
    split[d].push_back(Q(1)-Q(1,1073741824));
    detail::rankTwoDiagnostic("registry-fixture:reject-split-fan");
    rejects([&]{buildRankTwoArcRegistry(overlay,split,fans);},"ARRANGEMENT_UNRESOLVED",
        "independent interior sample split a simultaneous common-gate chord");
    rejects([&]{buildRankTwoArcRegistry(overlay,plan,{fans[0],fans[0]});},"INPUT_IDENTITY_MISMATCH",
        "duplicate fan manufactured repeated shared gate points");
    RankTwoArcSamplingPlan fullPlan;
    for(const auto& arc:arrangementRecords(full).arcs)if(arc.ends[0].location==RankTwoArcEndLocation::Finite
        && arc.ends[1].location==RankTwoArcEndLocation::Finite)fullPlan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    rejects([&]{buildRankTwoArcRegistry(full,fullPlan,fans);},"INPUT_IDENTITY_MISMATCH",
        "proof from another normal authority attached to shared edges");
    out<<"FAN_REGISTRY "<<shared.scientificHash()<<' '<<gates<<' '<<covered<<'\n'<<shared.canonicalRecord();
#else
    throw std::runtime_error("MISSING_FEATURE: shared registry edges bound to complete common-gate fans");
#endif
    out<<"CHECKED 4 REGISTRIES FAN_EDGE_BINDINGS 6\nchord_topology NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"registry evidence incomplete");
    std::cout<<"PASS shared actual-arc points/edges, reversed native region consumers, retained holes/infinity and point DAGs\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcRegistry.h")
    const bool trace=argc==3 && std::string(argv[2])=="--trace";
    if(argc!=2 && !trace)throw std::invalid_argument("registry evidence path and optional --trace required");
    if(trace){
        const auto sink=[](const char* phase,void*)noexcept{std::cerr<<"PHASE "<<phase<<'\n';};
        spin::detail::ScopedRankTwoDiagnostics observer(sink,nullptr);tests(argv[1]);
    }else tests(argv[1]);
    return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: one shared registry for actual directed boundary rings");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
