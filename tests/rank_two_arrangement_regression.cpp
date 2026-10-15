#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>
#include <set>
#if __has_include("spin/detail/RankTwoConstruction.h")
#include "spin/detail/RankTwoConstruction.h"
#include "fixtures/rank_two_source_cases.h"
#if __has_include("spin/detail/RankTwoNormalOverlay.h")
#include "spin/detail/RankTwoNormalOverlay.h"
#endif
#if __has_include("spin/detail/RankTwoNormalEvents.h")
#include "spin/detail/RankTwoNormalEvents.h"
#endif
#if __has_include("spin/detail/RankTwoArrangementRecords.h")
#include "spin/detail/RankTwoArrangementRecords.h"
#endif
#if __has_include("spin/detail/RankTwoNormalVertices.h")
#include "spin/detail/RankTwoNormalVertices.h"
#endif

namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void normalVertexTests(std::ostream& evidence) {
#if __has_include("spin/detail/RankTwoNormalVertices.h")
    using namespace spin;
    // h=x^2+y^2, G=(-2x,-2y,1): at x=3/4,y=0 the surface
    // radius sqrt(x^2+x^4)=15/16 is rational, giving an exact F=0
    // normal vertex rather than a nearby open-edge event.
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::disk(),"normal-endpoint");
    const auto global=analyzeRankTwoGlobalArrangement(source,CGAL::Gmpq(15,16));
    RankTwoNormalFace first,second;
    first.vertexIds={7,8,9};first.rays={{{CGAL::Gmpq(-3,2),0,1},{0,1,0},{0,0,-1}}};
    second.vertexIds={10,7,8};second.rays={{{1,0,0},{CGAL::Gmpq(-3,2),0,1},{0,1,0}}};
    const auto a=analyzeRankTwoNormalFace(source,global,first);
    const auto b=analyzeRankTwoNormalFace(source,global,second);
    const auto find=[](const auto& overlay,std::size_t id)->const RankTwoNormalVertexRecord& {
        for(const auto& v:normalVertexRecords(overlay))if(v.normalVertex==id)return v;
        throw std::runtime_error("missing normal vertex record");
    };
    const auto& va=find(a,7);const auto& vb=find(b,7);
    require(va.location==RankTwoNormalVertexLocation::FinitePositive && va.point
        && (*va.point)[0]==CGAL::Gmpq(3,4) && (*va.point)[1]==0,"wrong finite normal preimage");
    require(va.originalSigns==std::array<int,3>{0,1,1},"radius endpoint signs lost");
    require(va.identity==vb.identity && va.point==vb.point && va.originalSigns==vb.originalSigns,
        "normal endpoint split across incident faces");
    require(find(a,8).location==RankTwoNormalVertexLocation::Ideal && !find(a,8).point,
        "ideal normal received finite coordinates");
    require(find(a,9).location==RankTwoNormalVertexLocation::Negative && !find(a,9).point,
        "negative normal ray was admitted as finite positive");
    for(const auto* overlay:{&a,&b}){
        std::size_t matches=0;
        for(const auto& vertex:arrangementRecords(*overlay).vertices)if(vertex.identity==va.identity){
            ++matches;
            require(detail::CertifiedAlgebraicReal::reload(vertex.coordinates[0]).expression()==CORE::Expr(3)/4
                && detail::CertifiedAlgebraicReal::reload(vertex.coordinates[1]).expression()==CORE::Expr(0),
                "normal endpoint ID not attached to actual DCEL point");
        }
        require(matches==1,"normal endpoint absent from typed DCEL");
        for(const auto& event:normalBoundaryEvents(*overlay))
            require(event.identity!=va.identity && event.parameter>CORE::Expr(0)
                && event.parameter<CORE::Expr(1),"normal endpoint polluted open-edge root ordinals");
    }
    require(normalVertexRecords(global).empty(),"global has fabricated normal vertex records");
    evidence<<"ENDPOINT_FIRST\n"<<a.canonicalRecord()<<"ENDPOINT_SECOND\n"<<b.canonicalRecord();
#else
    (void)evidence;
    throw std::runtime_error("MISSING_FEATURE: cross-face normal vertex identity and ideal classification");
#endif
}
void structuredRecordTests(std::ostream& evidence) {
#if __has_include("spin/detail/RankTwoArrangementRecords.h")
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test");
    const auto global=analyzeRankTwoGlobalArrangement(source,rank_two_cases::radius());
    const auto& records=arrangementRecords(global);
    // Four disjoint nested circles give five regions and V-E=0. Additional
    // regular vertices at another branch's critical x are valid subdivisions;
    // do not require each circle to be represented by exactly two arcs.
    require(records.vertices.size()>=8 && records.vertices.size()==records.arcs.size()
            && records.regions.size()==5,
        "structured arrangement omitted a circle or region");
    std::map<std::string,const RankTwoArcRecord*> arcs;
    std::set<std::string> vertices,regions;
    for(const auto& v:records.vertices){
        require(vertices.insert(v.identity).second,"duplicate structured vertex");
        for(const auto& c:v.coordinates)detail::CertifiedAlgebraicReal::reload(c);
    }
    for(const auto& f:records.regions)require(regions.insert(f.identity).second,"duplicate structured region");
    for(const auto& a:records.arcs){
        require(arcs.emplace(a.identity,&a).second,"duplicate structured arc");
        for(const auto& end:a.ends)require(end.vertexIdentity && vertices.count(*end.vertexIdentity),"dangling finite arc endpoint");
        require(regions.count(a.leftRegion) && regions.count(a.rightRegion),"dangling arc region");
        require(!a.owners.empty(),"arc dropped original-polynomial owners");
    }
    int cuts=0,fulls=0;
    std::map<std::string,int> uses;
    for(const auto& f:records.regions){
        require(f.originalSigns[2]!=0 && f.originalSigns[0]!=0 && f.originalSigns[1]!=0,"missing exact region signs");
        if(f.kind==RankTwoRegionKind::Cut){
            ++cuts;
            require(f.outerLoops.size()==1 && f.innerLoops.size()==1,"CUT annulus hole missing from structured data");
        }
        if(f.kind==RankTwoRegionKind::Full)++fulls;
        auto check=[&](const auto& loops){for(const auto& loop:loops){
            require(!loop.empty(),"empty boundary loop");
            for(std::size_t i=0;i<loop.size();++i){
                const auto& ref=loop[i];const auto& next=loop[(i+1)%loop.size()];
                if(ref.atInfinity)continue;
                require(arcs.count(ref.identity) && arcs.count(next.identity),"dangling loop arc");
                const auto& a=*arcs.at(ref.identity);const auto& b=*arcs.at(next.identity);
                require((ref.reversed?a.rightRegion:a.leftRegion)==f.identity,"halfedge orientation/region mismatch");
                require(a.ends[ref.reversed?0:1].vertexIdentity==b.ends[next.reversed?1:0].vertexIdentity,
                    "boundary loop endpoints do not join");
                ++uses[ref.identity];
            }
        }};
        check(f.outerLoops);check(f.innerLoops);
    }
    require(cuts==2 && fulls==1,"original sign classification changed in structured view");
    for(const auto& a:arcs)require(uses[a.first]==2,"lost DCEL incidence");
    RankTwoNormalFace face;
    face.vertexIds={1,2,3};face.rays={{{8,0,1},{-8,0,1},{0,-12,1}}};
    const auto overlay=analyzeRankTwoNormalFace(source,global,face);
    const auto restored=reloadRankTwoNormalFace(source,global,face,overlay.canonicalRecord());
    require(restored.scientificHash()==overlay.scientificHash(),"normal overlay reload mismatch");
    require(normalBoundaryEvents(restored).size()==normalBoundaryEvents(overlay).size(),"reload dropped shared events");
    for(std::size_t i=0;i<normalBoundaryEvents(overlay).size();++i){
        const auto& a=normalBoundaryEvents(overlay)[i];const auto& b=normalBoundaryEvents(restored)[i];
        require(a.identity==b.identity && a.parameter==b.parameter,"reload changed exact event");
    }
    evidence<<"RELOADABLE_OVERLAY_BEGIN\n"<<overlay.canonicalRecord()<<"RELOADABLE_OVERLAY_END\n";
    auto reject=[&](const RankTwoNormalFace& f,const std::string& record){
        bool failed=false;
        try {reloadRankTwoNormalFace(source,global,f,record);}
        catch(const RankTwoSourceError& e){failed=e.code()=="INPUT_IDENTITY_MISMATCH";}
        require(failed,"forged normal arrangement accepted");
    };
    reject(face,overlay.canonicalRecord()+"passed=true\n");
    reject(face,overlay.canonicalRecord().substr(0,overlay.canonicalRecord().size()/2));
    auto wrong=face;wrong.rays[2][1]=-10;reject(wrong,overlay.canonicalRecord());
#else
    (void)evidence;
    throw std::runtime_error("MISSING_FEATURE: structured arrangement and checked normal reload");
#endif
}
void sharedEventTests(std::ostream& evidence) {
#if __has_include("spin/detail/RankTwoNormalEvents.h")
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test");
    const auto global=analyzeRankTwoGlobalArrangement(source,rank_two_cases::radius());
    RankTwoNormalFace above,below;
    above.vertexIds={1,2,3};above.rays={{{8,0,1},{-8,0,1},{0,-12,1}}};
    // Local edge 2->1 is reversed; the third ray is on the other side.
    below.vertexIds={2,1,4};below.rays={{{-8,0,1},{8,0,1},{0,12,1}}};
    const auto first=analyzeRankTwoNormalFace(source,global,above);
    const auto second=analyzeRankTwoNormalFace(source,global,below);
    std::vector<RankTwoNormalBoundaryEvent> roots[2];
    int side=0;
    for(const auto* overlay:{&first,&second}){
        for(const auto& event:normalBoundaryEvents(*overlay))
            if(event.edge==std::array<std::size_t,2>{1,2} && event.radiusOrdinal>=0)
                roots[side].push_back(event);
        ++side;
    }
    require(roots[0].size()==4 && roots[1].size()==4,"shared annulus edge must have four radius roots");
    // On y=0: F=x^4-3x^2+7/4. Normal interpolation t=(x+4)/8.
    // These independent closed-form roots also detect reversed ordinals and
    // accidental inclusion of negative normal combinations.
    const CORE::Expr inner=CORE::sqrt((CORE::Expr(3)-CORE::sqrt(CORE::Expr(2)))/2);
    const CORE::Expr outer=CORE::sqrt((CORE::Expr(3)+CORE::sqrt(CORE::Expr(2)))/2);
    const std::array<CORE::Expr,4> xs={-outer,-inner,inner,outer};
    for(int i=0;i<4;++i){
        require(roots[0][i].radiusOrdinal==i && roots[1][i].radiusOrdinal==i,"wrong canonical radius ordinal");
        require(roots[0][i].identity==roots[1][i].identity,"cracked cross-face root identity");
        require(roots[0][i].parameter==(xs[i]+4)/8,"wrong positive-cone normal parameter");
        require(roots[1][i].parameter==roots[0][i].parameter,"local edge orientation leaked into canonical parameter");
        evidence<<"SHARED_ROOT "<<i<<' '<<roots[0][i].identity<<'\n';
    }
#else
    (void)evidence;
    throw std::runtime_error("MISSING_FEATURE: shared normal edge event identities");
#endif
}
void normalFaceTests(std::ostream& evidence) {
#if __has_include("spin/detail/RankTwoNormalOverlay.h")
    using namespace spin;
    // G=(-2x,-2y,1): these rays invert to the PRE-REGISTERED triangle
    // (-4,-4),(4,-4),(0,6), strictly containing the radius-1.5 disk.
    // Dropping a face because its edges have no F roots must fail this test.
    RankTwoNormalFace triangle;
    triangle.vertexIds={10,11,12};
    triangle.rays={{{8,8,1},{-8,8,1},{0,-12,1}}};
    for(int kind=0;kind<3;++kind){
        auto coefficients=kind==0?rank_two_cases::disk():
            (kind==1?rank_two_cases::annulus():rank_two_cases::twoLobes());
        auto face=triangle;
        if(kind==2)for(auto& ray:face.rays)ray[1]=-ray[1];
        const auto source=PreparedRankTwoSource::fromExact(coefficients,"test");
        const auto global=analyzeRankTwoGlobalArrangement(source,rank_two_cases::radius());
        const auto result=analyzeRankTwoNormalFace(source,global,face);
        require(result.cutComponents()==(kind==2?2:1),"normal face lost hidden component");
        require(result.cutBoundaryLoops()==(kind==0?1:2),"normal face lost hidden loop/hole");
        // h=2+x²-y² in the saddle case is negative near each y pole of
        // the radius-1.5 disk. F>0 there gives TWO separate FULL slivers.
        const int fullCount[3]={0,1,2};
        require(result.fullComponents()==fullCount[kind],"normal face lost FULL interior");
        require(result.faceCount()==result.certifiedFaceCount(),"overlay face missing witness");
        std::swap(face.vertexIds[0],face.vertexIds[1]);
        std::swap(face.rays[0],face.rays[1]);
        const auto reversed=analyzeRankTwoNormalFace(source,global,face);
        evidence<<"NORMAL_BEGIN "<<kind<<" "<<result.scientificHash()<<'\n'<<result.canonicalRecord()<<"NORMAL_END\n";
        evidence<<"REVERSED_BEGIN "<<kind<<" "<<reversed.scientificHash()<<'\n'<<reversed.canonicalRecord()<<"REVERSED_END\n";
        require(result.canonicalRecord()==reversed.canonicalRecord(),"face permutation changed exact overlay");
    }
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::disk(),"test");
    const auto global=analyzeRankTwoGlobalArrangement(source,rank_two_cases::radius());
    // All these rays point into z<0, while G_z=1 everywhere. Its inverse
    // cone is empty even though its three supporting lines bound a triangle.
    for(auto& ray:triangle.rays)for(auto& q:ray)q=-q;
    const auto empty=analyzeRankTwoNormalFace(source,global,triangle);
    require(empty.cutComponents()==0 && empty.fullComponents()==0,"inequality signs discarded");
    bool rejected=false;
    triangle.rays[2]=triangle.rays[1];
    try {analyzeRankTwoNormalFace(source,global,triangle);}
    catch(const RankTwoSourceError& e){rejected=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(rejected,"singular normal triangle accepted");
    rejected=false;
    try {analyzeRankTwoNormalFace(PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test"),global,triangle);}
    catch(const RankTwoSourceError& e){rejected=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(rejected,"global cache from wrong source accepted");
#else
    (void)evidence;
    throw std::runtime_error("MISSING_FEATURE: certified fixed-normal-face overlay");
#endif
}
void tests(std::ostream& evidence) {
    using namespace spin;
    normalVertexTests(evidence);
    const auto disk=analyzeRankTwoGlobalArrangement(
        PreparedRankTwoSource::fromExact(rank_two_cases::disk(),"test"),rank_two_cases::radius());
    require(disk.cutComponents()==1 && disk.cutBoundaryLoops()==1,"hidden disk incomplete");
    const auto annulus=analyzeRankTwoGlobalArrangement(
        PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test"),rank_two_cases::radius());
    require(annulus.cutComponents()==1 && annulus.cutBoundaryLoops()==2,"annulus hole lost");
    require(annulus.fullComponents()==1,"full-ray annulus interior lost");
    const auto lobes=analyzeRankTwoGlobalArrangement(
        PreparedRankTwoSource::fromExact(rank_two_cases::twoLobes(),"test"),rank_two_cases::radius());
    require(lobes.cutComponents()==2 && lobes.cutBoundaryLoops()==2,"two-lobe topology lost");
    auto tangent=rank_two_cases::disk();tangent[0]=CGAL::Gmpq(5,4);
    const auto even=analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact(tangent,"test"),1);
    require(even.cutComponents()==0,"square-free zero set changed original F sign");
    bool doubleFactor=false;
    for(const auto& factor:even.factors())for(const auto& owner:factor.owners)
        if(owner.polynomial=="F" && owner.multiplicity==2)doubleFactor=true;
    require(doubleFactor,"even multiplicity not preserved");
    auto isolated=rank_two_cases::disk();isolated[0]=-1;
    const auto point=analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact(isolated,"test"),1);
    require(point.cutComponents()==0,"isolated root invented positive-area region");
    require(point.isolatedPoints()>0,"isolated point dropped");
    auto common=rank_two_cases::disk();common[0]=1;
    const auto shared=analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact(common,"test"),1);
    bool threeOwners=false;
    for(const auto& factor:shared.factors())if(factor.owners.size()==3)threeOwners=true;
    require(threeOwners,"shared F/h/D factor lost owners");
    for(const auto* result:{&disk,&annulus,&lobes,&even,&point,&shared}){
        require(result->faceCount()==result->certifiedFaceCount(),"open cell missing exact sample");
        require(!result->canonicalRecord().empty(),"missing DCEL evidence");
        evidence<<"BEGIN "<<result->scientificHash()<<'\n'<<result->canonicalRecord()<<"END\n";
    }
    const auto repeat=analyzeRankTwoGlobalArrangement(
        PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test"),rank_two_cases::radius());
    require(annulus.canonicalRecord()==repeat.canonicalRecord(),"arrangement identity nondeterministic");
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"test");
    const auto reloaded=RankTwoArrangement::reload(source,rank_two_cases::radius(),annulus.canonicalRecord());
    require(reloaded.scientificHash()==annulus.scientificHash(),"saved partition failed reload");
    bool mismatch=false;
    try {RankTwoArrangement::reload(source,rank_two_cases::radius(),annulus.canonicalRecord()+"passed=true\n");}
    catch(const RankTwoSourceError& e){mismatch=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(mismatch,"hand-written appended acceptance field accepted");
    mismatch=false;
    try {RankTwoArrangement::reload(source,rank_two_cases::radius(),
                                  annulus.canonicalRecord().substr(0,annulus.canonicalRecord().size()/2));}
    catch(const RankTwoSourceError& e){mismatch=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(mismatch,"truncated DCEL evidence accepted");
    mismatch=false;
    try {RankTwoArrangement::reload(PreparedRankTwoSource::fromExact(rank_two_cases::disk(),"test"),
                                  rank_two_cases::radius(),annulus.canonicalRecord());}
    catch(const RankTwoSourceError& e){mismatch=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(mismatch,"partition from another authority accepted");
    RankTwoSourceLimits tiny;tiny.maximumArrangementCells=1;
    bool refused=false;
    try {analyzeRankTwoGlobalArrangement(PreparedRankTwoSource::fromExact(rank_two_cases::disk(),"test"),
                                       rank_two_cases::radius(),tiny);}
    catch(const RankTwoSourceError& e){refused=e.code()=="RESOURCE_LIMIT";}
    require(refused,"cell budget ignored");
    normalFaceTests(evidence);
    sharedEventTests(evidence);
    structuredRecordTests(evidence);
    std::cout<<"PASS global disk/annulus/two-lobe regions, original signs, all open cells, isolated roots\n";
    std::cout<<annulus.scientificHash()<<'\n';
}
}
#endif
int main(int argc,char** argv) {
    try {
#if __has_include("spin/detail/RankTwoConstruction.h")
        if(argc!=2)throw std::invalid_argument("explicit evidence file required");
        std::ifstream previous(argv[1]);
        if(previous.good())throw std::invalid_argument("refusing to overwrite arrangement evidence");
        std::ofstream evidence(argv[1]);
        if(!evidence)throw std::runtime_error("cannot create arrangement evidence");
        tests(evidence);evidence.close();
        if(!evidence)throw std::runtime_error("failed to save arrangement evidence");
        return 0;
#else
        throw std::runtime_error("MISSING_FEATURE: complete global algebraic partition");
#endif
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
