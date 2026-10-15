#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoEndpointFan.h")
#include "spin/detail/RankTwoEndpointFan.h"
#include "spin/detail/RankTwoNormalOverlay.h"
namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using C=spin::detail::CertifiedAlgebraicReal;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,why);
}
void tests(const char* path){
    using namespace spin;
    // Frozen structural-equator method input. h=y-D, F=h^2-D;
    // every intersection between F/D/h has y=0. This includes actual
    // tangency of F and D at (1,0), plus transverse h and a normal edge.
    const RankTwoExactCoefficients q{1,0,-1,1,-1,-1,0,0,0,0};
    const auto source=PreparedRankTwoSource::fromExact(q,"endpoint-fan-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};
    cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto local=analyzeRankTwoNormalFace(source,global,cone);
    const auto& records=arrangementRecords(local);
    std::string start,negativeAnchor,positiveAnchor;
    for(const auto& v:records.vertices){
        const E x=C::reload(v.coordinates[0]).expression(),y=C::reload(v.coordinates[1]).expression();
        if(x==1 && y==0)start=v.identity;
        if(x>1 && y==-E(1)/2)negativeAnchor=v.identity;
        if(x==0 && y==1)positiveAnchor=v.identity;
    }
    require(!start.empty() && !negativeAnchor.empty() && !positiveAnchor.empty(),
        "frozen actual events missing");
    // Same registered error bound. The coarse gate is expected to fail;
    // one dyadic split is frozen before execution, rather than relaxing it.
    const auto coarse=RankTwoMovingGate::fromVertices(local,start,negativeAnchor,1,Q(1,256));
    rejects([&]{certifyRankTwoEndpointFan(local,coarse,
        {Q(-1,512),Q(0),Q(99,100),Q(101,100)},Q(1,1048576));},
        "ARRANGEMENT_UNRESOLVED","coarse F chord received an inaccurate fan certificate");
    const auto negative=RankTwoMovingGate::fromVertices(local,start,negativeAnchor,1,Q(1,512));
    const auto positive=RankTwoMovingGate::fromVertices(local,start,positiveAnchor,1,Q(1,1024));
    std::ofstream out(path);require(bool(out),"cannot write fan evidence");
    out<<"rank-two-common-endpoint-fan-v1\n";
    for(int side:{-1,1}){
        const auto& gate=side<0?negative:positive;
        const detail::RankTwoStripBox box{side<0?Q(-1,1024):Q(0),
            side<0?Q(0):Q(1,1024),Q(99,100),Q(101,100)};
        const auto fan=certifyRankTwoEndpointFan(local,gate,box,Q(1,1048576));
        const auto& arcs=fan.orderedArcs();
        require(arcs.size()==3 && fan.supportRecords().size()==local.factors().size(),
            "incident curve or normal owner omitted");
        std::vector<std::string> order;
        unsigned boundary=0,excluded=0;
        for(const auto& support:fan.supportRecords()){
            if(support.kind==RankTwoFanSupportKind::AtStartBoundary){
                ++boundary;require(support.owners.size()==1
                    && support.owners[0].polynomial=="NORMAL_EDGE:0:1",
                    "start-plane normal owner not retained");
            }
            if(support.kind==RankTwoFanSupportKind::Excluded)++excluded;
        }
        require(boundary==1 && excluded==2,"normal-boundary coverage not certified");
        for(const auto& arc:arcs){
            require(arc.gateIdentity==gate.identity()
                && arc.endpointIdentity==start && arc.arrangementIdentity==local.scientificHash(),
                "arcs used distinct gates or different input identities");
            require(arc.gateCoordinates[1].expression()==E(side)/1024
                && arc.strip.chordDistanceUpperBound<=Q(1,1048576),
                "common moving ordinate or whole-arc error lost");
            std::string owner;
            for(const auto& candidate:records.arcs)if(candidate.identity==arc.arcIdentity)
                for(const auto& tag:candidate.owners)
                    if(tag.polynomial=="F" || tag.polynomial=="D" || tag.polynomial=="h")owner=tag.polynomial;
            require(!owner.empty(),"gate not tied to actual material arc");
            order.push_back(owner);
        }
        require(order==(side<0?std::vector<std::string>{"F","D","h"}:
            std::vector<std::string>{"h","F","D"}),"tangent incident order changed");
        for(std::size_t i=1;i<arcs.size();++i){
            require(arcs[i-1].gateCoordinates[0].expression()<arcs[i].gateCoordinates[0].expression(),
                "overlapping incident chords received fan validity");
            require(arcs[i-1].aboveRegion==arcs[i].belowRegion,
                "actual DCEL region adjacency not preserved");
        }
        const auto again=certifyRankTwoEndpointFan(local,gate,box,Q(1,1048576));
        require(again.scientificHash()==fan.scientificHash()
            && again.canonicalRecord()==fan.canonicalRecord(),"fan repeat changed");
        out<<"FAN "<<side<<' '<<fan.scientificHash()<<' '<<gate.identity()<<' '
            <<arcs.size()<<' '<<fan.supportRecords().size()<<' '<<boundary<<' '<<excluded<<'\n';
        out<<fan.canonicalRecord();
    }
    rejects([&]{RankTwoMovingGate::fromVertices(local,start,"missing",1,Q(1,2));},
        "INPUT_IDENTITY_MISMATCH","unknown moving event accepted");
    rejects([&]{certifyRankTwoEndpointFan(global,negative,
        {Q(-1,1024),Q(0),Q(99,100),Q(101,100)},Q(1,1048576));},
        "INPUT_IDENTITY_MISMATCH","gate from another arrangement accepted");
    rejects([&]{RankTwoMovingGate::fromVertices(local,start,start,1,Q(1,2));},
        "ARRANGEMENT_UNRESOLVED","constant ordinate manufactured a common gate");
    // Frozen normal triangles, specified in projected (x,y) coordinates.
    // Their actual gradient rays are (-2*x,-2*y,1). The common gate crosses
    // one finite and one unbounded native normal arc in each configuration;
    // the second fan also requires a vertical native arc in the swapped chart.
    const auto circles=PreparedRankTwoSource::fromExact(
        {2,0,0,1,-1,-1,0,0,0,0},"normal-endpoint-fan-method-v1");
    const auto circleGlobal=analyzeRankTwoGlobalArrangement(circles,Q(3,2));
    for(int chart:{0,1}){
        RankTwoNormalFace normal;normal.vertexIds={0,1,2};
        if(chart==0)normal.rays={{{6,6,1},{-6,6,1},{0,-12,1}}};
        else normal.rays={{{0,-6,1},{0,0,1},{-6,-18,1}}};
        const auto overlay=analyzeRankTwoNormalFace(circles,circleGlobal,normal);
        const auto& native=arrangementRecords(overlay);
        std::string event,anchor;
        for(const auto& v:native.vertices){
            const E x=C::reload(v.coordinates[0]).expression(),y=C::reload(v.coordinates[1]).expression();
            if(x==0 && y==(chart==0?6:3))event=v.identity;
            if(x==3 && y==(chart==0?-3:9))anchor=v.identity;
        }
        require(!event.empty() && !anchor.empty(),"normal fan actual vertices absent");
        const auto gate=RankTwoMovingGate::fromVertices(overlay,event,anchor,chart,Q(1,1024));
        const detail::RankTwoStripBox box=chart==0?
            detail::RankTwoStripBox{Q(0),Q(3,1024),Q(59,10),Q(61,10)}:
            detail::RankTwoStripBox{Q(3),Q(3)+Q(6,1024),Q(-1,100),Q(1,100)};
        const auto fan=certifyRankTwoEndpointFan(overlay,gate,box,Q(1,1048576));
        require(fan.orderedArcs().size()==2
            && fan.supportRecords().size()==overlay.factors().size(),"normal fan coverage lost");
        unsigned infinite=0,vertical=0,excluded=0;
        for(const auto& support:fan.supportRecords()){
            if(support.kind==RankTwoFanSupportKind::Excluded)++excluded;
            else require(support.kind==RankTwoFanSupportKind::Graph,
                "normal graph incorrectly treated as start-plane boundary");
        }
        for(const auto& proof:fan.orderedArcs()){
            require(proof.gateIdentity==gate.identity() && proof.strip.chordDistanceUpperBound==0,
                "straight normal arc lost shared gate or zero chord error");
            const RankTwoArcRecord* arc=nullptr;
            for(const auto& a:native.arcs)if(a.identity==proof.arcIdentity)arc=&a;
            require(arc!=nullptr && arc->owners.size()==1
                && arc->owners[0].polynomial.find("NORMAL_EDGE:")==0,"normal source owner lost");
            for(const auto& end:arc->ends)infinite+=end.location!=RankTwoArcEndLocation::Finite;
            vertical+=arc->vertical;
            const E u=proof.gateCoordinates[chart].expression(),v=proof.gateCoordinates[1-chart].expression();
            require(u==(chart==0?E(3)/1024:E(3)+E(6)/1024),"normal gate not moving from actual events");
            require(chart==0?(v==6-3*u || v==6+3*u):(v==0 || v==(u-3)/2),
                "normal gate differs from independent analytic line");
        }
        require(infinite==1 && vertical==unsigned(chart==1) && excluded==4,
            "frozen fan failed to exercise native infinity/vertical membership");
        require(fan.orderedArcs()[0].aboveRegion==fan.orderedArcs()[1].belowRegion,
            "normal sector adjacency lost");
        out<<"NORMAL_FAN "<<chart<<' '<<fan.scientificHash()<<' '<<infinite<<' '<<vertical<<'\n';
        out<<fan.canonicalRecord();
    }
    out<<"CHECKED 4 FANS 10 ARC_SEGMENTS\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"fan write incomplete");
    std::cout<<"PASS actual common-gate fans: equator tangency, both sides, all material/normal owners, ordered chords and DCEL regions\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoEndpointFan.h")
    if(argc!=2)throw std::invalid_argument("fan evidence path required");
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: complete incident fan with one moving gate");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
