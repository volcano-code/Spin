#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoEndpointArc.h")
#include "spin/detail/RankTwoEndpointArc.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "fixtures/rank_two_source_cases.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=spin::detail::CertifiedAlgebraicReal;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* code,const char* why){
    bool rejected=false;try{f();}catch(const spin::RankTwoSourceError& e){rejected=e.code()==code;}
    require(rejected,why);
}
void tests(const char* output){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"endpoint-arc-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(3,2));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto local=analyzeRankTwoNormalFace(source,global,cone);
    const auto& records=arrangementRecords(local);
    const auto vertex=[&](const std::string& id)->const RankTwoVertexRecord&{
        for(const auto& v:records.vertices)if(v.identity==id)return v;
        throw std::runtime_error("missing actual vertex");
    };
    std::ofstream evidence(output);require(bool(evidence),"cannot write endpoint evidence");
    evidence<<"rank-two-endpoint-arc-chart-v1\n";
    unsigned checked=0;
    for(const auto& arc:records.arcs){
        bool isD=false;for(const auto& owner:arc.owners)isD|=owner.polynomial=="D";
        if(!isD || !arc.ends[0].vertexIdentity || !arc.ends[1].vertexIdentity)continue;
        const auto& l=vertex(*arc.ends[0].vertexIdentity);
        const auto& r=vertex(*arc.ends[1].vertexIdentity);
        // Both upper/lower quarter circles: gate orientation is not assumed.
        const E leftY=C::reload(l.coordinates[1]).expression();
        const int side=leftY.sign();
        if(C::reload(l.coordinates[0]).expression()!=0
            || leftY*E(side)!=E(3)/2
            || C::reload(r.coordinates[0]).expression()!=E(3)/2
            || C::reload(r.coordinates[1]).expression()!=0)continue;
        const detail::RankTwoStripBox box{side>0?Q(0):Q(-1,512),
            side>0?Q(1,512):Q(0),Q(149,100),Q(151,100)};
        const auto proof=certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(1,1024),box,Q(1,1048576));
        const auto& xy=proof.gateCoordinates;
        require(xy[1].expression()==E(3)*E(side)/2048,"moving gate fraction/axis was lost");
        require((xy[0].expression()*xy[0].expression()+xy[1].expression()*xy[1].expression()-E(9)/4).sign()==0,
            "gate not on actual circle");
        require(xy[0].expression()>0 && xy[1].expression()*E(side)>0,"wrong quarter/branch chosen after chart swap");
        require(proof.endpointIdentity==r.identity && proof.otherEndpointIdentity==l.identity
            && proof.arcIdentity==arc.identity && proof.arrangementIdentity==local.scientificHash()
            && proof.projectionAxis==1 && proof.fraction==Q(1,1024),
            "moving point provenance not bound to actual arc endpoints");
        require(proof.strip.chordDistanceUpperBound<=Q(1,1048576)
            && proof.checkedVertexCount==records.vertices.size(),
            "error bound or full event exclusion missing");
        const auto again=certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(1,1024),box,Q(1,1048576));
        require(proof.identity==again.identity && xy[0].expression()==again.gateCoordinates[0].expression(),
            "endpoint chart replay differs");
        // Native x chart fails at this smooth x-critical endpoint; it must
        // not be silently accepted with a zero denominator.
        rejects([&]{certifyRankTwoEndpointArc(local,arc.identity,1,0,Q(1,1024),
            {Q(149,100),Q(151,100),Q(-1,512),Q(1,512)},Q(1,1048576));},
            "ARRANGEMENT_UNRESOLVED","critical native projection accepted as regular");
        // Same supporting circle, wrong half: exact polynomial equality
        // alone is insufficient to establish membership in this arc.
        rejects([&]{certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(1,1024),
            {box[0],box[1],Q(-151,100),Q(-149,100)},Q(1,1048576));},
            "ARRANGEMENT_UNRESOLVED","opposite branch accepted");
        rejects([&]{certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(1,2),box,Q(1,1048576));},
            "ARRANGEMENT_UNRESOLVED","out-of-strip gate accepted");
        rejects([&]{certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(0),box,Q(1,1048576));},
            "INPUT_IDENTITY_MISMATCH","zero fraction accepted as new sample");
        rejects([&]{certifyRankTwoEndpointArc(local,arc.identity,1,2,Q(1,1024),box,Q(1,1048576));},
            "INPUT_IDENTITY_MISMATCH","invalid chart accepted");
        evidence<<"ARC "<<proof.arrangementIdentity<<' '<<proof.arcIdentity<<' '
            <<proof.endpointIdentity<<' '<<proof.otherEndpointIdentity<<' '
            <<proof.identity<<' '<<proof.projectionAxis<<' '<<proof.fraction<<' '
            <<proof.strip.inputIdentity<<' '<<proof.strip.chordDistanceUpperBound<<' '
            <<proof.checkedVertexCount<<'\n';
        for(const auto& v:xy){
            const auto& c=v.certificate();
            evidence<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
            for(const auto& a:c.coefficients)evidence<<a<<',';
            evidence<<'\n';
        }
        ++checked;
    }
    require(checked==2,"one gate orientation missing or duplicated");
    unsigned quartics=0,equalProjections=0;
    for(const auto& arc:records.arcs){
        bool isF=false;for(const auto& owner:arc.owners)isF|=owner.polynomial=="F";
        if(!isF || !arc.ends[0].vertexIdentity || !arc.ends[1].vertexIdentity)continue;
        const auto& l=vertex(*arc.ends[0].vertexIdentity);const auto& r=vertex(*arc.ends[1].vertexIdentity);
        const E rx=C::reload(r.coordinates[0]).expression();
        if(C::reload(l.coordinates[0]).expression()!=0 || C::reload(l.coordinates[1]).expression()<=0
            || rx<=E(89)/100 || rx>=E(9)/10 || C::reload(r.coordinates[1]).expression()!=0)continue;
        const auto proof=certifyRankTwoEndpointArc(local,arc.identity,1,1,Q(1,1024),
            {Q(0),Q(1,512),Q(89,100),Q(9,10)},Q(1,1048576));
        require(proof.gateCoordinates[1].expression()==C::reload(l.coordinates[1]).expression()/1024,
            "quartic critical point dependency lost");
        using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
        using P=Traits::Polynomial_2;
        const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1),s=x*x+y*y;
        auto& k=Traits::CKvA_2::instance().kernel();
        const auto point=k.construct_algebraic_real_2_object()(
            detail::rankTwoAlgebraicHandle(proof.gateCoordinates[0]),
            detail::rankTwoAlgebraicHandle(proof.gateCoordinates[1]));
        require(k.is_zero_at_2_object()(4*s*s-12*s+7,point),"gate is not on actual annulus quartic");
        require(proof.gateCoordinates[0].expression()>E(89)/100
            && proof.gateCoordinates[0].expression()<E(9)/10,"quartic gate switched inner/outer branch");
        evidence<<"QUARTIC "<<proof.identity<<' '<<proof.arcIdentity<<' '
            <<proof.strip.chordDistanceUpperBound<<' '<<proof.checkedVertexCount<<'\n';
        for(const auto& v:proof.gateCoordinates){
            const auto& c=v.certificate();
            evidence<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
            for(const auto& a:c.coefficients)evidence<<a<<',';
            evidence<<'\n';
        }
        ++quartics;
    }
    require(quartics==1,"actual quartic endpoint chart missing or duplicated");
    for(const auto& arc:arrangementRecords(global).arcs){
        bool isD=false;for(const auto& owner:arc.owners)isD|=owner.polynomial=="D";
        if(!isD)continue;
        rejects([&]{certifyRankTwoEndpointArc(global,arc.identity,1,1,Q(1,1024),
            {Q(0),Q(1,512),Q(149,100),Q(151,100)},Q(1,1048576));},
            "ARRANGEMENT_UNRESOLVED","equal endpoint ordinates produced a fabricated moving gate");
        ++equalProjections;
    }
    require(equalProjections==2,"global semicircle rejection coverage missing");
    rejects([&]{certifyRankTwoEndpointArc(local,"unknown",0,0,Q(1,2),
        {Q(0),Q(1),Q(0),Q(1)},Q(1));},
        "INPUT_IDENTITY_MISMATCH","unknown arc accepted");
    evidence<<"CHECKED "<<checked+quartics<<" EQUAL_PROJECTION_REJECTIONS "<<equalProjections
        <<"\nfull_incident_fan NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    evidence.flush();require(bool(evidence),"endpoint evidence write failed");
    std::cout<<"PASS actual smooth critical endpoint: alternate chart, moving gate, native branch, event exclusion; complete fan/mesh NOT_EXECUTED\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoEndpointArc.h")
    if(argc!=2)throw std::invalid_argument("evidence path required");
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: certified endpoint arc chart with moving gate");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
