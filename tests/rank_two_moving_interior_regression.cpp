#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoMovingInterior.h")
#include "spin/detail/RankTwoMovingInterior.h"
#include "spin/detail/RankTwoCarrierKernels.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){bool caught=false;
    try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);}
spin::RankTwoInteriorPointRequest request(const spin::RankTwoArrangement& a,spin::RankTwoRegionKind kind,
    int hsign,Q xf=Q(1,2),Q yf=Q(1,2)){
    const auto& native=spin::arrangementRecords(a);spin::RankTwoInteriorPointRequest r;
    r.abscissaFraction=xf;r.ordinateFraction=yf;
    const spin::RankTwoVertexRecord *left=nullptr,*right=nullptr;
    for(const auto& v:native.vertices){const auto x=C::reload(v.coordinates[0]).expression();
        if(x.sign()<0 && (!left || x>C::reload(left->coordinates[0]).expression()))left=&v;
        if(x.sign()>0 && (!right || x<C::reload(right->coordinates[0]).expression()))right=&v;}
    require(left && right,"fixture missing two native bounding events");r.eventIdentities={left->identity,right->identity};
    for(const auto& cell:native.regions)if(cell.kind==kind && (hsign==0 || cell.originalSigns[1]==hsign)){
        require(r.regionIdentity.empty(),"fixture expected one selected active native cell");r.regionIdentity=cell.identity;}
    require(!r.regionIdentity.empty(),"fixture active region absent");
    const auto selectedX=spin::detail::sampleRankTwoAlgebraicCoordinate(
        C::reload(left->coordinates[0]),C::reload(right->coordinates[0]),xf).expression();
    const auto eventX=[&](const std::string& id){for(const auto& v:native.vertices)
        if(v.identity==id)return C::reload(v.coordinates[0]).expression();
        throw std::runtime_error("frozen native arc end absent");};
    for(const auto& arc:native.arcs)if(!arc.vertical && arc.ends[0].location==spin::RankTwoArcEndLocation::Finite
        && arc.ends[1].location==spin::RankTwoArcEndLocation::Finite){
        if(!(eventX(*arc.ends[0].vertexIdentity)<selectedX && selectedX<eventX(*arc.ends[1].vertexIdentity)))continue;
        // A concentric annulus has two vertical gaps in the same native face.
        // This frozen test specifies the positive-y gap, not a new region or
        // a boundary selected after observing an optimization result.
        if(kind==spin::RankTwoRegionKind::Cut && C::reload(arc.interiorWitness[1]).expression().sign()<=0)continue;
        if(arc.leftRegion==r.regionIdentity){require(r.boundaryArcIdentities[0].empty(),"ambiguous lower arc in frozen fixture");r.boundaryArcIdentities[0]=arc.identity;}
        if(arc.rightRegion==r.regionIdentity){require(r.boundaryArcIdentities[1].empty(),"ambiguous upper arc in frozen fixture");r.boundaryArcIdentities[1]=arc.identity;}}
    require(!r.boundaryArcIdentities[0].empty() && !r.boundaryArcIdentities[1].empty(),"selected cell lacks two actual native fibers");
    return r;
}
void save(std::ostream& out,const char* name,const spin::RankTwoMovingInteriorPoint& p){
    out<<"CASE "<<name<<" POINT "<<p.scientificHash()<<'\n'<<p.canonicalRecord()<<"CASE_END\n";
    std::cout<<"PASS "<<name<<" exact moving interior definition\n";
}
void tests(const char* path){
    using namespace spin;std::ostringstream out;out<<"rank-two-moving-interior-test-v1\n";
    const auto source=PreparedRankTwoSource::fromExact({2,0,0,1,1,1,0,0,0,0},"moving-interior-full-disk-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(1));const auto centered=request(a,RankTwoRegionKind::Full,-1);
    const auto center=buildRankTwoMovingInteriorPoint(a,centered);
    require(center.coordinates()[0].expression().sign()==0 && center.coordinates()[1].expression().sign()==0,
        "native midpoint definitions failed their analytic disk value");save(out,"full-disk-center",center);
    for(const Q R:{Q(1),Q(3,2)}){
        const auto authority=analyzeRankTwoGlobalArrangement(source,R);
        const auto req=request(authority,RankTwoRegionKind::Full,-1,Q(1,4),Q(1,3));
        const auto point=buildRankTwoMovingInteriorPoint(authority,req);
        const auto x=point.coordinates()[0].expression(),y=point.coordinates()[1].expression();
        // Hand-derived: x=-R/2; y=-R*sqrt(3)/6. No builder supplies expected values.
        std::ostringstream numerator,denominator;numerator<<R.numerator();denominator<<R.denominator();
        const E rr(CORE::BigRat(CORE::BigInt(numerator.str()),CORE::BigInt(denominator.str())));
        require((2*x+rr).sign()==0 && y.sign()<0 && (12*y*y-rr*rr).sign()==0,
            "point replaced by a fixed origin, rational witness or wrong fiber branch");
        const E d=rr*rr-x*x-y*y;require(d.sign()>0,"native disk internal point lacks positive material");
        require(reloadRankTwoMovingInteriorPoint(authority,req,point.canonicalRecord()).scientificHash()==point.scientificHash(),
            "moving definition reload changed scientific identity");
        save(out,R==Q(1)?"full-disk-offcenter":"scaled-disk-offcenter",point);
    }
    const auto ann=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"moving-interior-annulus-v1");
    const auto b=analyzeRankTwoGlobalArrangement(ann,Q(3,2));
    for(const auto [kind,sign]:{std::pair{RankTwoRegionKind::Full,-1},std::pair{RankTwoRegionKind::Cut,-1},std::pair{RankTwoRegionKind::Cut,1}}){
        const auto req=request(b,kind,sign);const auto point=buildRankTwoMovingInteriorPoint(b,req);
        const auto x=point.coordinates()[0].expression(),y=point.coordinates()[1].expression();
        const E h=x*x+y*y-2,D=E(9)/4-x*x-y*y;require(D.sign()>0,"annulus interior outside ball projection");
        const E upper=CORE::sqrt(D),lower=kind==RankTwoRegionKind::Full?-upper:h;
        require((upper-lower).sign()>0,"actual active annulus point has zero/negative thickness");
        save(out,kind==RankTwoRegionKind::Full?"annulus-full":sign<0?"annulus-cut-lower":"annulus-cut-upper",point);
    }
    auto bad=centered;bad.abscissaFraction=0;
    rejects([&]{buildRankTwoMovingInteriorPoint(a,bad);},"INPUT_IDENTITY_MISMATCH","event endpoint accepted as interior");
    bad=centered;std::swap(bad.boundaryArcIdentities[0],bad.boundaryArcIdentities[1]);
    rejects([&]{buildRankTwoMovingInteriorPoint(a,bad);},"SOURCE_TOPOLOGY_INVALID","reversed native fiber gap accepted");
    bad=centered;bad.regionIdentity="caller-invented-cell";
    rejects([&]{buildRankTwoMovingInteriorPoint(a,bad);},"INPUT_IDENTITY_MISMATCH","invented native region accepted");
    auto wide=request(b,RankTwoRegionKind::Full,-1);
    const auto& n=arrangementRecords(b);
    for(const auto& arc:n.arcs)if(arc.owners.size()==1 && arc.owners[0].polynomial=="D"){
        wide.eventIdentities={*arc.ends[0].vertexIdentity,*arc.ends[1].vertexIdentity};break;}
    rejects([&]{buildRankTwoMovingInteriorPoint(b,wide);},"ARRANGEMENT_UNRESOLVED","nonadjacent event slab accepted");
    auto corrupt=center.canonicalRecord();corrupt+="caller_passed true\n";
    rejects([&]{reloadRankTwoMovingInteriorPoint(a,centered,corrupt);},"INPUT_IDENTITY_MISMATCH","altered moving evidence accepted");
    RankTwoInteriorPointSettings zero;zero.maximumFiberRoots=0;
    rejects([&]{buildRankTwoMovingInteriorPoint(a,centered,zero);},"RESOURCE_LIMIT","zero proof budget granted internal point");
    rejects([&]{reloadRankTwoMovingInteriorPoint(b,request(b,RankTwoRegionKind::Full,-1),center.canonicalRecord());},
        "INPUT_IDENTITY_MISMATCH","foreign arrangement reused interior definition");
    std::ofstream saved(path);require(bool(saved),"cannot save interior evidence");saved<<out.str();saved.flush();require(bool(saved),"incomplete interior evidence write");
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoMovingInterior.h")
    if(argc!=2)throw std::invalid_argument("moving interior evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: native adjacent-event/fiber-defined moving internal nodes");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
