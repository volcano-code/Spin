#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include "spin/detail/RankTwoDiagnostics.h"
#if __has_include("spin/detail/RankTwoCertifiedVertices.h")
#include "spin/detail/RankTwoCertifiedVertices.h"
#endif
#if __has_include("spin/detail/RankTwoRegistryEdge.h")
#include "spin/detail/RankTwoRegistryEdge.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <variant>
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=spin::detail::CertifiedAlgebraicReal;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class A,class R>
auto attemptedProof(const A& a,const R& r,const std::string& id,std::size_t index,int axis,
    const spin::detail::RankTwoStripBox& b,const Q& e,spin::detail::RankTwoRegularStripMemo* m,int)
    ->decltype(tryCertifyRankTwoRegistryEdgeAtIndex(a,r,id,index,axis,b,e,m)){
    return tryCertifyRankTwoRegistryEdgeAtIndex(a,r,id,index,axis,b,e,m);
}
template<class A,class R>
std::variant<spin::RankTwoRegularEdge,std::string> attemptedProof(const A&,const R&,
    const std::string&,std::size_t,int,const spin::detail::RankTwoStripBox&,const Q&,
    spin::detail::RankTwoRegularStripMemo*,long){
    throw std::runtime_error("MISSING_FEATURE: typed guarded native edge attempt");
}
template<class F> void rejects(F f,const char* code,const char* why){
    bool rejected=false;try{f();}catch(const spin::RankTwoSourceError& e){rejected=e.code()==code;}
    require(rejected,why);
}
// Break caught: a stale/foreign index must not select a different current
// edge, and the indexed path must still perform the complete native proof.
// ADL feature detection keeps RED an actual runtime assertion, not a build error.
template<class A,class R>
auto indexedProof(const A& a,const R& r,const std::string& id,std::size_t index,
    int axis,const spin::detail::RankTwoStripBox& box,const Q& error,int)
    ->decltype(certifyRankTwoRegistryEdgeAtIndex(a,r,id,index,axis,box,error,nullptr)){
    return certifyRankTwoRegistryEdgeAtIndex(a,r,id,index,axis,box,error,nullptr);
}
template<class A,class R>
spin::RankTwoRegularEdge indexedProof(const A&,const R&,const std::string&,
    std::size_t,int,const spin::detail::RankTwoStripBox&,const Q&,long){
    throw std::runtime_error("MISSING_FEATURE: guarded current-registry indexed regular-edge proof");
}
void retainedNativeVertices(const spin::RankTwoArrangement& a){
#if __has_include("spin/detail/RankTwoCertifiedVertices.h")
    std::size_t reconstructed=0;
    const auto observe=[](const char* phase,void* context) noexcept {
        if(std::strcmp(phase,"algebraic-root:reload-sturm")==0)
            ++*static_cast<std::size_t*>(context);
    };
    const auto& vertices=spin::arrangementRecords(a).vertices;
    require(!vertices.empty(),"owner regression has no actual native vertices");
    for(const auto& vertex:vertices){
        // The independent serialized reload validates the original complete
        // polynomial/ordinal/isolator, before observing owner access.
        const auto x=C::reload(vertex.coordinates[0]);
        const auto y=C::reload(vertex.coordinates[1]);
        spin::detail::ScopedRankTwoDiagnostics measured(observe,&reconstructed);
        const auto& p=spin::certifiedRankTwoVertex(a,vertex.identity);
        require(p[0].expression()==x.expression() && p[1].expression()==y.expression(),
            "retained native coordinates differ from authenticated records");
        for(int axis=0;axis<2;++axis){
            const auto& actual=p[axis].certificate();
            const auto& expected=vertex.coordinates[axis];
            require(actual.coefficients==expected.coefficients
                && actual.rootIndex==expected.rootIndex && actual.lower==expected.lower
                && actual.upper==expected.upper && actual.rootIdentity==expected.rootIdentity,
                "owner lost the original complete root certificate");
        }
    }
    require(reconstructed==0,"current immutable event values were reloaded per access");
    rejects([&]{spin::certifiedRankTwoVertex(a,"foreign-native-vertex");},
        "INPUT_IDENTITY_MISMATCH","caller-selected unknown event acquired native authority");
#else
    throw std::runtime_error("MISSING_FEATURE: arrangement-owned certified native event values");
#endif
}
struct Selected {std::string arc;spin::RankTwoArcSamplingPlan plan;};
Selected selected(const spin::RankTwoArrangement& a,bool vertical,const Q& center){
    Selected out;
    for(const auto& arc:spin::arrangementRecords(a).arcs){
        if(arc.ends[0].location!=spin::RankTwoArcEndLocation::Finite
            || arc.ends[1].location!=spin::RankTwoArcEndLocation::Finite)continue;
        out.plan[arc.identity]={Q(1,2)};
        bool owner=false;for(const auto& o:arc.owners)
            owner|=o.polynomial==(vertical?"NORMAL_EDGE:0:2":"F");
        const E y=C::reload(arc.interiorWitness[1]).expression();
        // The native exact dyadic witness of the intended positive vertical
        // interval is exactly 1, not a strictly larger rounded midpoint.
        // Include that witness without changing the arc, sample fractions,
        // strip box or requested chord accuracy.
        if(arc.vertical==vertical && owner && (vertical?(y>=E(1) && y<E(13)/10):(y>0 && y<1))){
            require(out.arc.empty(),"frozen registry-edge fixture is ambiguous");
            out.arc=arc.identity;
        }
    }
    if(out.arc.empty())for(const auto& arc:spin::arrangementRecords(a).arcs){
        if(!arc.vertical)continue;
        std::cerr<<"unselected vertical arc "<<arc.identity<<" owners";
        for(const auto& owner:arc.owners)std::cerr<<' '<<owner.polynomial;
        std::cerr<<" y="<<C::reload(arc.interiorWitness[1]).expression().doubleValue()
            <<" finite="<<(arc.ends[0].location==spin::RankTwoArcEndLocation::Finite)
            <<','<<(arc.ends[1].location==spin::RankTwoArcEndLocation::Finite)<<'\n';
    }
    require(!out.arc.empty(),"frozen actual arc unavailable");
    out.plan.at(out.arc)={center-Q(1,65536),center+Q(1,65536)};
    return out;
}
std::string middleEdge(const spin::RankTwoArcRegistry& registry,const std::string& arc){
    std::string id;
    for(const auto& e:registry.edges())if(e.arcIdentity==arc){
        const auto& a=registry.points()[e.points[0]],&b=registry.points()[e.points[1]];
        if(a.kind==spin::RankTwoRegistryPointKind::ArcSample
            && b.kind==spin::RankTwoRegistryPointKind::ArcSample){
            require(id.empty(),"frozen adjacent samples selected multiple edges");id=e.identity;
        }
    }
    require(!id.empty(),"adjacent exact samples not published");return id;
}
void check(const spin::RankTwoArrangement& authority,const Selected& chosen,int axis,
    const spin::detail::RankTwoStripBox& box,std::ostream& out){
    using namespace spin;
    const auto registry=buildRankTwoArcRegistry(authority,chosen.plan);
    const auto edgeId=middleEdge(registry,chosen.arc);
    const auto proof=certifyRankTwoRegistryEdge(authority,registry,edgeId,axis,box,Q(1,1048576));
    require(proof.arrangementIdentity()==authority.scientificHash()
        && proof.registryIdentity()==registry.scientificHash() && proof.edgeIdentity()==edgeId,
        "regular proof bound a different actual edge");
    require(proof.strip().chordDistanceUpperBound<=Q(1,1048576),"regular actual chord exceeds locked error");
    require(proof.obstacles().size()+1==authority.factors().size()
        && proof.strip().obstacleSigns.size()==proof.obstacles().size(),"native factor/normal obstacle omitted");
    require(proof.checkedVertexCount()==arrangementRecords(authority).vertices.size(),
        "regular proof did not consider all native events");
    require(proof.pointIdentities()[0]!=proof.pointIdentities()[1]
        && proof.belowRegion()!=proof.aboveRegion(),"regular proof collapsed chord or adjacent regions");
    for(const auto& e:registry.edges())if(e.identity==edgeId){
        require(proof.pointIdentities()[0]==registry.points()[e.points[0]].identity
            && proof.pointIdentities()[1]==registry.points()[e.points[1]].identity,
            "proof endpoints were regenerated rather than read from shared registry");
    }
    const auto again=certifyRankTwoRegistryEdge(authority,registry,edgeId,axis,box,Q(1,1048576));
    require(proof.scientificHash()==again.scientificHash() && proof.canonicalRecord()==again.canonicalRecord(),
        "same regular proof changed identity or scientific output");
    std::size_t index=registry.edges().size();
    for(std::size_t i=0;i<registry.edges().size();++i)if(registry.edges()[i].identity==edgeId)index=i;
    require(index<registry.edges().size() && registry.edges().size()>1,"actual indexed fixture incomplete");
    const auto fast=indexedProof(authority,registry,edgeId,index,axis,box,Q(1,1048576),0);
    require(fast.canonicalRecord()==proof.canonicalRecord() && fast.scientificHash()==proof.scientificHash(),
        "indexed entry changed actual native proof, chart, events or authority");
    spin::detail::RankTwoRegularStripMemo memo;
    const auto typed=attemptedProof(authority,registry,edgeId,index,axis,box,Q(1,1048576),&memo,0);
    require(typed.index()==0 && std::get<0>(typed).canonicalRecord()==proof.canonicalRecord(),
        "typed attempt skipped or changed actual native guards");
    rejects([&]{attemptedProof(authority,registry,edgeId,registry.edges().size(),axis,box,
        Q(1,1048576),&memo,0);},"INPUT_IDENTITY_MISMATCH","typed attempt swallowed invalid authority");
    if(axis==0){
        for(int repeat=0;repeat<2;++repeat){
            const auto rejected=attemptedProof(authority,registry,edgeId,index,axis,box,
                Q(1,1073741824),&memo,0);
            require(rejected.index()==1 && std::get<1>(rejected)=="regular graph chord bound exceeds requested accuracy",
                "typed native rejection granted proof or lost witness");
        }
    }
    rejects([&]{indexedProof(authority,registry,edgeId,registry.edges().size(),axis,box,Q(1,1048576),0);},
        "INPUT_IDENTITY_MISMATCH","out-of-range index acquired a native proof");
    rejects([&]{indexedProof(authority,registry,edgeId,(index+1)%registry.edges().size(),axis,box,Q(1,1048576),0);},
        "INPUT_IDENTITY_MISMATCH","stale index selected another current edge");
    rejects([&]{indexedProof(authority,registry,"foreign-edge",index,axis,box,Q(1,1048576),0);},
        "INPUT_IDENTITY_MISMATCH","foreign identity acquired indexed authority");
    rejects([&]{indexedProof(authority,registry,edgeId,index,2,box,Q(1,1048576),0);},
        "INPUT_IDENTITY_MISMATCH","invalid chart bypassed indexed proof guard");
    rejects([&]{certifyRankTwoRegistryEdge(authority,registry,"not-an-edge",axis,box,Q(1,1048576));},
        "INPUT_IDENTITY_MISMATCH","foreign shared edge acquired a regular proof");
    if(axis==0)rejects([&]{certifyRankTwoRegistryEdge(authority,registry,edgeId,axis,box,Q(1,1073741824));},
        "ARRANGEMENT_UNRESOLVED","accuracy failure silently relaxed its chord bound");
    if(axis==0)rejects([&]{indexedProof(authority,registry,edgeId,index,axis,box,Q(1,1073741824),0);},
        "ARRANGEMENT_UNRESOLVED","indexed entry relaxed its chord bound");
    if(axis==1 && !chosen.plan.empty()){
        bool vertical=false;for(const auto& arc:registry.arcs())if(arc.identity==chosen.arc)vertical=arc.vertical;
        if(vertical){
            require(proof.strip().chordDistanceUpperBound==0,"vertical native line not represented exactly");
            rejects([&]{certifyRankTwoRegistryEdge(authority,registry,edgeId,0,
                {box[2],box[3],box[0],box[1]},Q(1,1048576));},
                "ARRANGEMENT_UNRESOLVED","zero-projection chart accepted for vertical edge");
        }
    }
    out<<"REGULAR_EDGE "<<proof.scientificHash()<<'\n'<<proof.canonicalRecord();
}
void tests(const char* path){
    using namespace spin;
    std::ofstream out(path);require(bool(out),"cannot open regular-edge evidence");
    out<<"rank-two-registry-regular-edge-test-v1\n";
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"registry-edge-annulus-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(3,2));
    retainedNativeVertices(global);
    const auto ordinary=selected(global,false,Q(1,2));
    check(global,ordinary,0,{Q(-1,4096),Q(1,4096),Q(89,100),Q(9,10)},out);
    // An isometric chart swap on a noncritical, non-axis-aligned quartic
    // branch. Bounds were frozen from the analytic annulus before any solve.
    check(global,selected(global,false,Q(9,16)),1,
        {Q(88344,100000),Q(88349,100000),Q(1110,10000),Q(1116,10000)},out);
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto overlay=analyzeRankTwoNormalFace(source,global,cone);
    retainedNativeVertices(overlay);
#if __has_include("spin/detail/RankTwoCertifiedVertices.h")
    bool checkedForeign=false;
    for(const auto& v:arrangementRecords(overlay).vertices){
        bool present=false;
        for(const auto& old:arrangementRecords(global).vertices)present|=old.identity==v.identity;
        if(!present){
            rejects([&]{certifiedRankTwoVertex(global,v.identity);},
                "INPUT_IDENTITY_MISMATCH","overlay event was reused through foreign authority");
            checkedForeign=true;break;
        }
    }
    require(checkedForeign,"foreign authority fixture lacks a new overlay event");
#endif
    check(overlay,selected(overlay,true,Q(1,2)),1,
        {Q(114,100),Q(116,100),Q(-1,100000),Q(1,100000)},out);
    const auto registry=buildRankTwoArcRegistry(global,ordinary.plan);
    rejects([&]{certifyRankTwoRegistryEdge(overlay,registry,middleEdge(registry,ordinary.arc),0,
        {Q(-1,4096),Q(1,4096),Q(89,100),Q(9,10)},Q(1,1048576));},
        "INPUT_IDENTITY_MISMATCH","foreign arrangement supplied incomplete normal exclusion");
    const auto edgeId=middleEdge(registry,ordinary.arc);std::size_t index=0;
    while(registry.edges().at(index).identity!=edgeId)++index;
    rejects([&]{indexedProof(overlay,registry,edgeId,index,0,
        {Q(-1,4096),Q(1,4096),Q(89,100),Q(9,10)},Q(1,1048576),0);},
        "INPUT_IDENTITY_MISMATCH","indexed entry accepted foreign arrangement authority");
    rejects([&]{certifyRankTwoRegistryEdge(global,registry,middleEdge(registry,ordinary.arc),0,
        {Q(-1),Q(1),Q(1,2),Q(1)},Q(1,1048576));},
        "ARRANGEMENT_UNRESOLVED","coarse/intersecting graph strip accepted as whole-chord proof");
    out<<"CHECKED 3 ACTUAL_REGULAR_EDGES\nfull_chord_coverage NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"regular-edge evidence incomplete");
    std::cout<<"PASS actual shared-edge strips, complete factor exclusion, native events and swapped/vertical charts\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoRegistryEdge.h")
    if(argc!=2)throw std::invalid_argument("regular-edge evidence path required");
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: regular strips bound to actual registry chords and all native factors");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
