#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cstring>
#include <chrono>
#if __has_include("spin/detail/RankTwoDiagnostics.h")
#include "spin/detail/RankTwoDiagnostics.h"
#endif
#if __has_include("spin/detail/RankTwoWholeArrangement.h")
#include "spin/detail/RankTwoWholeArrangement.h"
#include "spin/detail/RankTwoArrangementRecords.h"
#include "fixtures/rank_two_source_cases.h"
#if __has_include("spin/detail/RankTwoSettings.h")
#include "spin/detail/RankTwoSettings.h"
#endif
#if __has_include("spin/detail/RankTwoSharedIntervals.h")
#include "spin/detail/RankTwoSharedIntervals.h"
#endif
namespace {
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
void tests(const char* output){
    using namespace spin;
    using Q=CGAL::Gmpq;
    PreparedRankTwoNormalMesh::Rays rays={{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}};
    const PreparedRankTwoNormalMesh::Faces faces={{{0,2,4}},{{2,1,4}},{{1,3,4}},{{3,0,4}},
        {{2,0,5}},{{1,2,5}},{{3,1,5}},{{0,3,5}}};
    const auto normal=PreparedRankTwoNormalMesh::fromExact(rays,faces);
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"whole-test");
#if __has_include("spin/detail/RankTwoSettings.h")
    QuadricTriangulationSettings settings;
    settings.rankTwoArrangementSource=true;settings.icosphereSubdivisions=0;settings.truncationRadius=1.5;
    const auto configured=analyzeRankTwoArrangement(source,settings);
    const auto& configuredRecords=wholeArrangementRecords(configured);
    require(configuredRecords.normals.rays().size()==12 && configuredRecords.normalFaces.size()==20
        && configuredRecords.sharedEdges.size()==30,"settings path changed the requested normal grid");
    const auto rejectSettings=[&](const QuadricTriangulationSettings& s){
        bool refused=false;try{analyzeRankTwoArrangement(source,s);}
        catch(const RankTwoSourceError& e){refused=e.code()=="INPUT_IDENTITY_MISMATCH" || e.code()=="RESOURCE_LIMIT";}
        require(refused,"invalid/conflicting rank-two settings silently accepted");
    };
    auto bad=settings;bad.rankTwoArrangementSource=false;rejectSettings(bad);
    bad=settings;bad.fixedNormalParaboloid=true;rejectSettings(bad);
    bad=settings;bad.algebraicFiniteMap=true;rejectSettings(bad);
    bad=settings;bad.paperProjectiveHyperboloid=true;rejectSettings(bad);
    bad=settings;bad.truncationRadius=0;rejectSettings(bad);
    bad=settings;bad.maximumIcosphereSubdivisions=2;rejectSettings(bad);
    bad=settings;bad.icosphereSubdivisions=20;rejectSettings(bad);
    bad=settings;bad.icosphereSubdivisions=-1;rejectSettings(bad);
    bad=settings;bad.truncationRadius=std::numeric_limits<double>::quiet_NaN();rejectSettings(bad);
    bad=settings;bad.rankTwoLimits.maximumVertices=11;rejectSettings(bad);
    bad=settings;bad.rankTwoLimits.maximumArrangementCells=49;rejectSettings(bad);
#else
    throw std::runtime_error("MISSING_FEATURE: approved settings-only rank-two arrangement entry");
#endif
    const auto result=analyzeRankTwoArrangement(source,Q(3,2),normal);
    const auto& all=wholeArrangementRecords(result);
#if __has_include("spin/detail/RankTwoDiagnostics.h")
    std::array<unsigned,5> observed{};
    const auto count=[](const char* event,void* data) noexcept {
        auto& n=*static_cast<std::array<unsigned,5>*>(data);
        ++n[0];if(std::strcmp(event,"overlay:begin")==0)++n[1];
        if(std::strcmp(event,"overlay:end")==0)++n[2];
        if(std::strcmp(event,"events-certificates:begin")==0)n[3]=1;
        if(std::strcmp(event,"events-certificates:end")==0)n[3]=0;
        if(n[3] && std::strcmp(event,"algebraic-root:certify")==0)++n[4];
    };
    {
        detail::ScopedRankTwoDiagnostics observer(count,&observed);
        const auto local=analyzeRankTwoNormalFace(source,all.global,normal.face(0));
        require(local.canonicalRecord()==all.normalFaces[0].canonicalRecord(),"diagnostics changed scientific record");
    }
    require(observed[0]>2 && observed[1]==1 && observed[2]==1,"missing or unbalanced overlay diagnostics");
    require(observed[3]==0 && observed[4]==0,
        "event collection recertifies roots already certified in the same finalize");
    const auto countAfter=observed[0];detail::rankTwoDiagnostic("outside-scope");
    require(observed[0]==countAfter,"diagnostic sink leaked after scope");
#else
    throw std::runtime_error("MISSING_FEATURE: scoped read-only rank-two phase diagnostics");
#endif
    require(all.normalFaces.size()==8 && all.sharedEdges.size()==12 && all.normalVertices.size()==6,
        "whole normal atlas incomplete");
    require(all.global.cutComponents()==1 && all.global.cutBoundaryLoops()==2,"global hole lost in atlas");
    std::size_t exactCells=0;
    const auto cells=[](const auto& a){const auto& r=arrangementRecords(a);return r.regions.size()+r.vertices.size()+r.arcs.size();};
    exactCells+=cells(all.global);
    for(const auto& face:all.normalFaces){
        exactCells+=cells(face);
        require(face.faceCount()==face.certifiedFaceCount(),"whole atlas has unclassified face");
    }
    require(all.totalArrangementCells==exactCells,"whole budget counts do not match actual objects");
    std::size_t openRadiusEvents=0;
    for(const auto& edge:all.sharedEdges){
        require(edge.incidentFaces[0]!=edge.incidentFaces[1],"shared edge lacks two owners");
        for(const auto& event:edge.events){
            require(event.edge==edge.edge,"shared event assigned to wrong edge");
            if(event.radiusOrdinal>=0)++openRadiusEvents;
        }
    }
    require(openRadiusEvents==8,"annulus quadrant edges should have eight open F roots");
    const auto reloaded=reloadRankTwoArrangement(source,Q(3,2),normal,result.canonicalRecord());
    require(result.canonicalRecord()==reloaded.canonicalRecord(),"whole-grid replay changed evidence");
    bool failed=false;
    try{reloadRankTwoArrangement(source,Q(3,2),normal,result.canonicalRecord()+"passed=true\n");}
    catch(const RankTwoSourceError& e){failed=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(failed,"whole-grid forged pass accepted");
    RankTwoSourceLimits small;small.maximumArrangementCells=exactCells-1;failed=false;
    try{analyzeRankTwoArrangement(source,Q(3,2),normal,small);}
    catch(const RankTwoSourceError& e){failed=e.code()=="RESOURCE_LIMIT";}
    require(failed,"global total budget was enforced only per face");
    // Canonical shared registries do not depend on which incident face ran first.
    auto reversed=faces;std::reverse(reversed.begin(),reversed.end());
    const auto second=analyzeRankTwoArrangement(source,Q(3,2),PreparedRankTwoNormalMesh::fromExact(rays,reversed));
    const auto& other=wholeArrangementRecords(second);
    for(std::size_t i=0;i<all.sharedEdges.size();++i){
        require(all.sharedEdges[i].edge==other.sharedEdges[i].edge,"shared edge ordering changed");
        const auto& a=all.sharedEdges[i].events;const auto& b=other.sharedEdges[i].events;
        require(a.size()==b.size(),"face visitation changed event count");
        for(std::size_t j=0;j<a.size();++j)require(a[j].identity==b[j].identity && a[j].parameter==b[j].parameter,
            "face visitation changed exact shared event");
    }
    for(std::size_t i=0;i<all.normalVertices.size();++i)
        require(all.normalVertices[i].identity==other.normalVertices[i].identity,"normal vertex registry depends on face traversal");
    failed=false;
    try{arrangementRecords(result);}
    catch(const RankTwoSourceError&){failed=true;}
    require(failed,"whole object exposed an empty local DCEL as valid");
    std::string coincidentRecord;
#if __has_include("spin/detail/RankTwoSharedIntervals.h")
    // q=z-xy is regular rank two, but h=xy vanishes along ALL four
    // positive equator-to-north normal edges. Isolated roots alone cannot
    // certify these shared zero strata or their material-side labels.
    const auto saddle=PreparedRankTwoSource::fromExact({0,0,0,1,0,0,0,0,0,Q(-1,2)},"coincident-test");
    const auto coincident=analyzeRankTwoArrangement(saddle,Q(1),normal);
    const auto& c=wholeArrangementRecords(coincident);
    std::size_t positiveEdges=0;
    for(const auto& edge:c.sharedEdges){
        const auto& intervals=sharedIntervals(edge);
        if(edge.edge[1]!=4){require(intervals.empty(),"ideal/negative edge has positive finite interval");continue;}
        ++positiveEdges;
        require(intervals.size()>=2,"shared h=0 edge lost its inside/outside split");
        require(intervals.front().lowerParameter==CORE::Expr(0) && intervals.back().upperParameter==CORE::Expr(1),
            "shared interval misses ideal or finite normal endpoint");
        CORE::Expr last(0);bool inside=false,outside=false;
        for(const auto& part:intervals){
            require(part.lowerParameter==last && part.lowerParameter<part.upperParameter,"gap/overlap in shared interval cover");
            last=part.upperParameter;
            require(part.originalSigns[1]==0,"coincident h sign lost");
            bool owner=false;for(const auto& o:part.originalOwners)if(o.polynomial=="h" && o.multiplicity==1)owner=true;
            require(owner,"coincident h multiplicity lost");
            require(!part.sourceArcs[0].empty() && !part.sourceArcs[1].empty(),"shared interval not bound to both actual DCEL arcs");
            inside|=part.originalSigns[2]>0;outside|=part.originalSigns[2]<0;
        }
        require(inside && outside,"shared interval signs not independently classified");
    }
    require(positiveEdges==4,"wrong number of h=0 normal edges");
    coincidentRecord=coincident.canonicalRecord();
    // An invertible rational image of an octahedron. Several different
    // positive normal cones now share x=0 as their supporting graph line;
    // their intervals and events must NOT be copied by line identity.
    const PreparedRankTwoNormalMesh::Rays tilted={
        {{Q(4,5),0,Q(3,5)}},{{Q(-4,5),0,Q(3,5)}},{{0,Q(-4,5),Q(3,5)}},
        {{Q(-4,5),0,Q(-3,5)}},{{Q(4,5),0,Q(-3,5)}},{{0,Q(4,5),Q(-3,5)}}};
    const PreparedRankTwoNormalMesh::Faces tiltedFaces={
        {{0,1,2}},{{0,5,1}},{{0,2,4}},{{0,4,5}},{{3,2,1}},{{3,1,5}},{{3,4,2}},{{3,5,4}}};
    const auto difficult=analyzeRankTwoArrangement(saddle,Q(1),PreparedRankTwoNormalMesh::fromExact(tilted,tiltedFaces));
    const auto findEdge=[&](std::array<std::size_t,2> key)->const RankTwoSharedEdgeRecord&{
        for(const auto& e:wholeArrangementRecords(difficult).sharedEdges)if(e.edge==key)return e;
        throw std::runtime_error("missing tilted shared edge");
    };
    const auto& center=findEdge({0,1});
    require(center.events.size()==3,"coincident height branch lost its center event");
    const std::array<CORE::Expr,3> exactParameters={CORE::Expr(1)/8,CORE::Expr(1)/2,CORE::Expr(7)/8};
    for(int i=0;i<3;++i)require(center.events[i].parameter==exactParameters[i],"wrong exact finite-edge event parameter");
    require(center.events[1].zeroPolynomials==std::vector<std::string>{"h"}
        && center.events[1].radiusOrdinal==-1,"height crossing falsely counted as radius root");
    const auto& unbounded=findEdge({0,4});
    require(unbounded.events.empty(),"events copied from another cone on the same line");
    const auto& intervals=sharedIntervals(unbounded);
    require(!intervals.empty() && intervals.front().lowerParameter==CORE::Expr(0)
        && intervals.back().upperParameter==CORE::Expr(1)/2,"ideal edge limit not recovered exactly");
    for(const auto& part:intervals)require(part.sourceReversed[0] && part.sourceReversed[1],
        "unbounded shared parameter orientation lost");
    require(sharedIntervals(findEdge({3,4})).empty() && findEdge({3,4}).events.empty(),
        "negative cone reused positive supporting-line data");
    coincidentRecord+=difficult.canonicalRecord();
#else
    throw std::runtime_error("MISSING_FEATURE: shared normal-edge arc intervals and coincident owners");
#endif
    std::ifstream old(output);require(!old.good(),"refuse overwrite whole evidence");
    std::ofstream saved(output);saved<<result.canonicalRecord()<<second.canonicalRecord()<<coincidentRecord;
#if __has_include("spin/detail/RankTwoSettings.h")
    saved<<configured.canonicalRecord();
#endif
    saved.close();
    require(bool(saved),"cannot save whole evidence");
    std::cout<<"PASS whole normal atlas: 8 faces,12 shared edges,8 open F roots,6 normal vertices\n"
        <<result.scientificHash()<<'\n';
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoWholeArrangement.h")
    if(argc!=2 && (argc!=3 || std::string(argv[2])!="--trace"))
        throw std::invalid_argument("new whole arrangement evidence path and optional --trace required");
#if __has_include("spin/detail/RankTwoDiagnostics.h")
    auto origin=std::chrono::steady_clock::now();
    const auto trace=[](const char* phase,void* context) noexcept {
        const auto start=*static_cast<std::chrono::steady_clock::time_point*>(context);
        const auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"TRACE "<<ns<<' '<<phase<<'\n'<<std::flush;
    };
    spin::detail::ScopedRankTwoDiagnostics diagnostic(argc==3?+trace:nullptr,&origin);
#endif
    tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: whole-grid arrangement and shared registry");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
