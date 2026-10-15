#include "detail/RankTwoRegistryRotation.h"
#include "detail/RankTwoCarrierKernels.h"
#include "verification/Sha256.h"
#include <CGAL/Cartesian.h>
#include <algorithm>
#include <locale>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Direction=std::pair<std::string,int>; // native arc and outgoing end
using K=CGAL::Cartesian<CORE::Expr>;
[[noreturn]]void invalid(const std::string& why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void topology(const std::string& why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
std::string text(const Direction& d){return d.first+":"+std::to_string(d.second);}
std::string hash(const std::string& bytes){verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());return h.finish();}
}
struct RankTwoRegistryRotation::Impl {
    std::string record,hash;
    std::string arrangement,registry,incidence,partition;
    std::uint64_t points=0,directions=0,omitted=0,excluded=0;
    std::size_t maximumDegree=0;
};
RankTwoRegistryRotation::RankTwoRegistryRotation(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoRegistryRotation::arrangementIdentity()const{return impl_->arrangement;}
const std::string& RankTwoRegistryRotation::registryIdentity()const{return impl_->registry;}
const std::string& RankTwoRegistryRotation::incidenceIdentity()const{return impl_->incidence;}
const std::string& RankTwoRegistryRotation::partitionIdentity()const{return impl_->partition;}
std::uint64_t RankTwoRegistryRotation::pointCount()const{return impl_->points;}
std::uint64_t RankTwoRegistryRotation::finiteDirectionCount()const{return impl_->directions;}
std::uint64_t RankTwoRegistryRotation::omittedInfiniteDirectionCount()const{return impl_->omitted;}
std::uint64_t RankTwoRegistryRotation::omittedOutsidePartitionDirectionCount()const{return impl_->excluded;}
std::size_t RankTwoRegistryRotation::maximumFiniteDegree()const{return impl_->maximumDegree;}
const std::string& RankTwoRegistryRotation::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoRegistryRotation::scientificHash()const{return impl_->hash;}

RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const RankTwoChordIncidence& incidence,
    const RankTwoRegistryRotationSettings& settings){
    std::vector<std::size_t> points(registry.points().size()),edges(registry.edges().size());
    std::iota(points.begin(),points.end(),0);std::iota(edges.begin(),edges.end(),0);
    return RankTwoRegistryRotation::certify(arrangement,registry,incidence,"",points,edges,settings);
}
RankTwoRegistryRotation RankTwoRegistryRotation::certify(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const RankTwoChordIncidence& incidence,const std::string& partition,
    const std::vector<std::size_t>& selectedPoints,const std::vector<std::size_t>& selectedEdges,
    const RankTwoRegistryRotationSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(registry.arrangementIdentity()!=arrangement.scientificHash()
        ||incidence.arrangementIdentity()!=arrangement.scientificHash()
        ||incidence.registryIdentity()!=registry.scientificHash()
        ||incidence.partitionIdentity()!=partition)
        invalid("finite rotation requires the same private arrangement, registry and incidence");
    const auto& native=arrangementRecords(arrangement);
    const auto& points=registry.points();const auto& edges=registry.edges();
    const std::set<std::size_t> pointSet(selectedPoints.begin(),selectedPoints.end()),edgeSet(selectedEdges.begin(),selectedEdges.end());
    if(pointSet.size()!=selectedPoints.size() || edgeSet.size()!=selectedEdges.size()
        ||(!pointSet.empty() && *pointSet.rbegin()>=points.size())
        ||(!edgeSet.empty() && *edgeSet.rbegin()>=edges.size()))invalid("invalid immutable native rotation scope");
    std::uint64_t used=0;
    const auto consume=[&]{if(used==settings.maximumRotationSteps)
        throw RankTwoSourceError("RESOURCE_LIMIT","complete native rotation step budget exhausted");++used;};
    std::map<std::string,const RankTwoArcRecord*> arcs;
    std::map<std::string,std::set<Direction>> nativeStars;
    std::map<std::string,std::map<Direction,Direction>> nativeNext;
    for(const auto& v:native.vertices)nativeStars.emplace(v.identity,std::set<Direction>{});
    std::ostringstream rows;rows.imbue(std::locale::classic());
    for(const auto& arc:native.arcs){
        if(!arcs.emplace(arc.identity,&arc).second)topology("duplicate native arc in rotation");
        rows<<"NATIVE_ARC "<<arc.identity;
        for(int k=0;k<2;++k){
            const auto& end=arc.ends[k];rows<<' '<<int(end.location)<<' '
                <<(end.vertexIdentity?*end.vertexIdentity:"NONE");
            if(end.location==RankTwoArcEndLocation::Finite){
                if(!end.vertexIdentity || !nativeStars.count(*end.vertexIdentity))
                    topology("finite native direction has no original event");
                nativeStars.at(*end.vertexIdentity).insert({arc.identity,k});
            }
        }
        rows<<'\n';
    }
    // Face is to the left of each directed boundary walk. At a finite corner,
    // the twin of the incoming halfedge rotates CLOCKWISE to the outgoing one.
    // Use all walks, including those whose other ends lie at infinity.
    for(const auto& region:native.regions)for(bool inner:{false,true}){
        const auto& loops=inner?region.innerLoops:region.outerLoops;
        for(std::size_t l=0;l<loops.size();++l){
            const auto& loop=loops[l];rows<<"NATIVE_LOOP "<<region.identity<<' '
                <<(inner?"INNER":"OUTER")<<' '<<l;
            for(const auto& ref:loop)rows<<' '<<(ref.atInfinity?"INFINITY":(ref.reversed?"-":"+")+ref.identity);
            rows<<'\n';
            for(std::size_t j=0;j<loop.size();++j){
                consume();const auto& current=loop[j];
                if(current.atInfinity)continue;
                const auto& arc=*arcs.at(current.identity);
                const int end=current.reversed?0:1;
                if(arc.ends[end].location!=RankTwoArcEndLocation::Finite)continue;
                const auto& next=loop[(j+1)%loop.size()];
                if(next.atInfinity)topology("finite native event connected to fictitious infinity");
                const auto& nextArc=*arcs.at(next.identity);const int start=next.reversed?1:0;
                if(nextArc.ends[start].location!=RankTwoArcEndLocation::Finite
                    ||nextArc.ends[start].vertexIdentity!=arc.ends[end].vertexIdentity)
                    topology("native face walk changes event at its corner");
                auto& successor=nativeNext[*arc.ends[end].vertexIdentity];
                if(!successor.emplace(Direction{arc.identity,end},Direction{next.identity,start}).second)
                    topology("native outgoing direction has two face successors");
            }
        }
    }
    for(const auto& event:nativeStars){
        const auto& star=event.second;const auto& successor=nativeNext[event.first];
        if(star.size()!=successor.size())topology("native event rotation omitted a direction");
        if(star.empty())continue;
        std::set<Direction> visited;auto direction=*star.begin();
        do {consume();if(!star.count(direction) || !visited.insert(direction).second)
                topology("native finite event rotation is not one cycle");
            const auto found=successor.find(direction);
            if(found==successor.end())topology("native direction lacks a successor");
            direction=found->second;
        }while(direction!=*star.begin());
        if(visited.size()!=star.size())topology("native event split into multiple rotation cycles");
    }
    std::vector<std::vector<std::size_t>> actualStars(points.size());
    for(const auto e:selectedEdges)for(const auto p:edges[e].points){
        if(!pointSet.count(p))topology("finite edge references point outside its private scope");
        actualStars[p].push_back(e);
    }
    std::vector<K::Point_2> positions;
    for(const auto& p:points)positions.emplace_back(p.coordinates[0].expression(),p.coordinates[1].expression());
    const auto cmp=[&](std::size_t a,std::size_t b,int k){
        return int(detail::rankTwoAlgebraicHandle(points[a].coordinates[k]).compare(
            detail::rankTwoAlgebraicHandle(points[b].coordinates[k])));
    };
    const auto other=[&](std::size_t e,std::size_t p){
        const auto& pair=edges[e].points;
        if(pair[0]!=p && pair[1]!=p)topology("rotation edge not incident to point");
        return pair[pair[0]==p?1:0];
    };
    std::map<std::string,const RankTwoRegistryArc*> sharedArcs;
    for(const auto& arc:registry.arcs())sharedArcs.emplace(arc.identity,&arc);
    auto out=std::make_shared<RankTwoRegistryRotation::Impl>();out->points=selectedPoints.size();out->partition=partition;
    out->arrangement=arrangement.scientificHash();out->registry=registry.scientificHash();out->incidence=incidence.scientificHash();
    for(const auto p:selectedPoints){
        consume();const auto& point=points[p];auto order=actualStars[p];
        std::map<std::size_t,std::size_t> expected;
        std::size_t nativeDegree=0;
        if(point.kind==RankTwoRegistryPointKind::Event){
            const auto found=nativeStars.find(point.eventIdentity);
            if(found==nativeStars.end())topology("registry event absent from native rotation");
            const auto& star=found->second;const auto& next=nativeNext.at(point.eventIdentity);
            nativeDegree=star.size();std::map<Direction,std::size_t> retained;
            for(const auto& d:star){
                const auto& arc=*sharedArcs.at(d.first);
                if(!arc.finite){++out->omitted;continue;}
                if(arc.edges.empty())topology("finite native direction lacks a shared chord");
                const auto edge=d.second?arc.edges.back():arc.edges.front();
                if(!edgeSet.count(edge)){++out->excluded;continue;}
                if(edges.at(edge).points[d.second?1:0]!=p)
                    topology("shared event edge disagrees with native outgoing end");
                retained.emplace(d,edge);
            }
            if(retained.size()!=order.size())topology("actual event has a nonnative finite direction");
            for(const auto& d:star)rows<<"NATIVE_STEP "<<p<<' '<<text(d)<<" NEXT "<<text(next.at(d))<<'\n';
            for(const auto& entry:retained){
                auto d=next.at(entry.first);std::size_t traversed=0;
                while(!retained.count(d)){consume();d=next.at(d);
                    if(++traversed>star.size())topology("native projected successor never returns");}
                expected.emplace(entry.second,retained.at(d));
            }
        }else{
            if(order.size()!=2 || edges[order[0]].arcIdentity!=point.arcIdentity
                ||edges[order[1]].arcIdentity!=point.arcIdentity)
                topology("moving source point lacks its degree-two same-arc star");
            expected.emplace(order[0],order[1]);expected.emplace(order[1],order[0]);
        }
        std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){
            consume();if(a==b)return false;
            const auto x=other(a,p),y=other(b,p);
            const auto half=[&](std::size_t q){const int sy=cmp(q,p,1);
                return sy>0 ||(sy==0 && cmp(q,p,0)>0)?0:1;};
            const int hx=half(x),hy=half(y);if(hx!=hy)return hx<hy;
            const int sign=int(CGAL::orientation(positions[p],positions[x],positions[y]));
            if(!sign)topology("indistinguishable same-half finite chord directions");
            return sign>0;
        });
        rows<<"STAR "<<p<<' '<<(point.kind==RankTwoRegistryPointKind::Event?"EVENT":"SAMPLE")
            <<" NATIVE "<<nativeDegree<<" FINITE "<<order.size()<<"\nGEOMETRIC_CCW "<<p;
        for(const auto e:order)rows<<' '<<e;rows<<"\nEXPECTED_CW "<<p;
        for(const auto& entry:expected)rows<<' '<<entry.first<<':'<<entry.second;rows<<'\n';
        for(std::size_t i=0;i<order.size();++i){
            consume();const auto predecessor=order[(i+order.size()-1)%order.size()];
            if(expected.at(order[i])!=predecessor)
                topology("actual finite chord cyclic order differs from native DCEL at "+point.identity);
        }
        out->directions+=order.size();out->maximumDegree=std::max(out->maximumDegree,order.size());
    }
    if(out->directions!=2*selectedEdges.size())topology("finite rotation did not consume both edge ends");
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<(partition.empty()?"rank-two-registry-rotation-v1":"rank-two-partition-registry-rotation-v1")
        <<"\nAUTHORITY "<<arrangement.scientificHash()<<' '<<registry.scientificHash()<<' '<<incidence.scientificHash();
    if(!partition.empty())record<<"\nPARTITION_SCOPE "<<partition<<"\nEXCLUDED_FINITE_DIRECTIONS "<<out->excluded;
    record<<"\nMAX_STEPS "<<settings.maximumRotationSteps
        <<"\nCOUNTS "<<out->points<<' '<<out->directions<<' '<<out->maximumDegree<<' '
        <<out->omitted<<' '<<used<<"\nNATIVE_ARRANGEMENT_BEGIN\n"<<arrangement.canonicalRecord()
        <<"NATIVE_ARRANGEMENT_END\nCHORD_INCIDENCE_BEGIN\n"<<incidence.canonicalRecord()
        <<"CHORD_INCIDENCE_END\n"<<rows.str()
        <<"finite_star_native_rotation PASSED\nnative_arc_isotopy NOT_GRANTED\nunbounded_arc_geometry NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();out->hash=hash("rank-two-registry-rotation-certificate-v1\n"+out->record);
    detail::requireExactAlgebraicEnvironment();return RankTwoRegistryRotation(std::move(out));
}
} // namespace spin
