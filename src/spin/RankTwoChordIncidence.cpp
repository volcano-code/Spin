#include "detail/RankTwoChordIncidence.h"
#include "detail/RankTwoChordSweep.h"
#include "detail/RankTwoCarrierKernels.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <CGAL/Cartesian.h>
#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <locale>
#include <sstream>

namespace spin {
namespace {
using K=CGAL::Cartesian<CORE::Expr>;
[[noreturn]]void unresolved(const std::string& why){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",why);}
std::string hash(const std::string& bytes){verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());return h.finish();}
}
struct RankTwoChordIncidence::Impl {
    std::string arrangement,registry,partition,record,hash;
    std::uint64_t pointPairs=0,pointEdges=0,edgePairs=0;
    std::uint64_t groups=0;
};
RankTwoChordIncidence::RankTwoChordIncidence(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoChordIncidence::arrangementIdentity()const{return impl_->arrangement;}
const std::string& RankTwoChordIncidence::registryIdentity()const{return impl_->registry;}
const std::string& RankTwoChordIncidence::partitionIdentity()const{return impl_->partition;}
std::uint64_t RankTwoChordIncidence::pointPairCount()const{return impl_->pointPairs;}
std::uint64_t RankTwoChordIncidence::pointEdgeCount()const{return impl_->pointEdges;}
std::uint64_t RankTwoChordIncidence::edgePairCount()const{return impl_->edgePairs;}
std::uint64_t RankTwoChordIncidence::predicateGroupCount()const{return impl_->groups;}
const std::string& RankTwoChordIncidence::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoChordIncidence::scientificHash()const{return impl_->hash;}

RankTwoChordIncidence certifyRankTwoChordIncidence(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const RankTwoChordIncidenceSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(registry.arrangementIdentity()!=arrangement.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","chord incidence requires the actual registry arrangement");
    auto out=std::make_shared<RankTwoChordIncidence::Impl>();
    out->arrangement=arrangement.scientificHash();out->registry=registry.scientificHash();
    const auto& points=registry.points();const auto& edges=registry.edges();
    std::vector<K::Point_2> positions;positions.reserve(points.size());
    for(const auto& p:points)positions.emplace_back(p.coordinates[0].expression(),p.coordinates[1].expression());
    const auto compare=[&](std::size_t a,std::size_t b,int axis){
        return int(detail::rankTwoAlgebraicHandle(points[a].coordinates[axis]).compare(
            detail::rankTwoAlgebraicHandle(points[b].coordinates[axis])));
    };
    const auto ordered=[&](std::size_t p,std::size_t a,std::size_t b){
        return compare(p,a,0)*compare(p,b,0)<=0 && compare(p,a,1)*compare(p,b,1)<=0;
    };
    const auto orientation=[&](std::size_t a,std::size_t b,std::size_t c){
        return int(CGAL::orientation(positions[a],positions[b],positions[c]));
    };
    std::uint64_t used=0;
    const auto consume=[&]{
        if(used==settings.maximumIncidenceTests)
            throw RankTwoSourceError("RESOURCE_LIMIT","complete chord incidence test budget exhausted at registry "+out->registry);
        ++used;
    };
    std::ostringstream rows;rows.imbue(std::locale::classic());
    detail::rankTwoDiagnostic("chord-incidence:point-pairs");
    for(std::size_t a=0;a<points.size();++a)for(std::size_t b=a+1;b<points.size();++b){
        consume();++out->pointPairs;const int x=compare(a,b,0),y=compare(a,b,1);
        if(x==0 && y==0)unresolved("distinct shared point identities coincide: "+points[a].identity+" "+points[b].identity);
        rows<<"DISTINCT_POINT "<<a<<' '<<b<<' '<<x<<' '<<y<<'\n';
    }
    detail::rankTwoDiagnostic("chord-incidence:point-edges");
    for(std::size_t e=0;e<edges.size();++e){
        const auto a=edges[e].points[0],b=edges[e].points[1];
        if(a>=points.size() || b>=points.size() || a==b)unresolved("invalid actual shared chord: "+edges[e].identity);
        for(std::size_t p=0;p<points.size();++p){
            if(p==a || p==b)continue;
            consume();++out->pointEdges;const int o=orientation(a,b,p);
            const bool between=o==0 && ordered(p,a,b);
            if(between)unresolved("nonincident shared point lies on actual chord: "+points[p].identity+" "+edges[e].identity);
            rows<<"POINT_EDGE "<<p<<' '<<e<<" ORIENTATION "<<o<<" ORDERED "<<between<<'\n';
        }
    }
    detail::rankTwoDiagnostic("chord-incidence:edge-pairs");
    for(std::size_t i=0;i<edges.size();++i)for(std::size_t j=i+1;j<edges.size();++j){
        consume();++out->edgePairs;
        const auto a=edges[i].points[0],b=edges[i].points[1],c=edges[j].points[0],d=edges[j].points[1];
        int common=-1,shared=0;
        for(const auto p:{a,b})if(p==c || p==d){common=static_cast<int>(p);++shared;}
        if(shared==2)unresolved("duplicate actual shared chords: "+edges[i].identity+" "+edges[j].identity);
        const std::array<int,4> signs{orientation(a,b,c),orientation(a,b,d),orientation(c,d,a),orientation(c,d,b)};
        const bool intersects=(signs[0]*signs[1]<0 && signs[2]*signs[3]<0)
            || (signs[0]==0 && ordered(c,a,b)) || (signs[1]==0 && ordered(d,a,b))
            || (signs[2]==0 && ordered(a,c,d)) || (signs[3]==0 && ordered(b,c,d));
        if(!shared && intersects)unresolved("nonincident actual coarse chords intersect: "+edges[i].identity+" "+edges[j].identity);
        if(shared){
            const auto p=static_cast<std::size_t>(common),x=a==p?b:a,y=c==p?d:c;
            if(orientation(p,x,y)==0 && (ordered(x,p,y) || ordered(y,p,x)))
                unresolved("shared-endpoint actual chords overlap: "+edges[i].identity+" "+edges[j].identity);
        }
        rows<<"EDGE_PAIR "<<i<<' '<<j<<" SHARED "<<common<<" SIGNS";
        for(const auto s:signs)rows<<' '<<s;
        rows<<" RESULT "<<(shared?"SHARED_ENDPOINT":"DISJOINT")<<'\n';
    }
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-chord-incidence-v1\nAUTHORITY "<<out->arrangement<<' '<<out->registry
        <<"\nMAX_TESTS "<<settings.maximumIncidenceTests<<"\nCOUNTS "<<out->pointPairs<<' '
        <<out->pointEdges<<' '<<out->edgePairs<<"\nSHARED_REGISTRY_BEGIN\n"
        <<registry.canonicalRecord()<<"SHARED_REGISTRY_END\n"<<rows.str()
        <<"actual_finite_chord_incidence PASSED\nnative_arc_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->groups=used;
    out->record=record.str();out->hash=hash("rank-two-chord-incidence-certificate-v1\n"+out->record);
    detail::requireExactAlgebraicEnvironment();return RankTwoChordIncidence(std::move(out));
}

RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const RankTwoChordSweepSettings& settings){
    std::vector<std::size_t> points(registry.points().size()),edges(registry.edges().size());
    std::iota(points.begin(),points.end(),0);std::iota(edges.begin(),edges.end(),0);
    return RankTwoChordIncidence::sweep(arrangement,registry,"",points,edges,settings);
}
RankTwoChordIncidence RankTwoChordIncidence::sweep(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,const std::string& partition,
    const std::vector<std::size_t>& selectedPoints,const std::vector<std::size_t>& selectedEdges,
    const RankTwoChordSweepSettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(registry.arrangementIdentity()!=arrangement.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","exact chord sweep requires the actual registry arrangement");
    const auto& points=registry.points();const auto& edges=registry.edges();
    const std::size_t n=selectedPoints.size(),m=selectedEdges.size();
    const std::set<std::size_t> pointSet(selectedPoints.begin(),selectedPoints.end()),edgeSet(selectedEdges.begin(),selectedEdges.end());
    if(pointSet.size()!=n || edgeSet.size()!=m
        ||(!pointSet.empty() && *pointSet.rbegin()>=points.size())
        ||(!edgeSet.empty() && *edgeSet.rbegin()>=edges.size()))
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","invalid immutable incidence scope");
    for(const auto e:selectedEdges)for(const auto p:edges[e].points)if(!pointSet.count(p))
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","partition edge endpoint missing from its point scope");
    auto out=std::make_shared<RankTwoChordIncidence::Impl>();
    out->arrangement=arrangement.scientificHash();out->registry=registry.scientificHash();
    out->partition=partition;
    const auto product=[](std::uint64_t a,std::uint64_t b){
        if(a && b>std::numeric_limits<std::uint64_t>::max()/a)
            throw RankTwoSourceError("RESOURCE_LIMIT","exact sweep coverage count overflow");
        return a*b;
    };
    out->pointPairs=n?product(n,n-1)/2:0;
    out->edgePairs=m?product(m,m-1)/2:0;
    if(m && n<2)unresolved("actual chord has no two distinct endpoints");
    out->pointEdges=m?product(m,n-2):0;
    std::vector<K::Point_2> positions;positions.reserve(points.size());
    for(const auto& p:points)positions.emplace_back(p.coordinates[0].expression(),p.coordinates[1].expression());
    std::ostringstream rows,trace;rows.imbue(std::locale::classic());trace.imbue(std::locale::classic());
    const auto consume=[&]{
        if(out->groups==settings.maximumPredicateGroups)
            throw RankTwoSourceError("RESOURCE_LIMIT","exact sweep predicate budget exhausted at registry "+out->registry);
        ++out->groups;
    };
    const auto compare=[&](std::size_t a,std::size_t b,int axis){
        consume();const int value=int(detail::rankTwoAlgebraicHandle(points[a].coordinates[axis]).compare(
            detail::rankTwoAlgebraicHandle(points[b].coordinates[axis])));
        trace<<"C "<<a<<' '<<b<<' '<<axis<<' '<<value<<'\n';return value;
    };
    const auto orientation=[&](std::size_t a,std::size_t b,std::size_t c){
        consume();
        // Same immutable registry index means EXACTLY the same point. Avoid
        // asking CORE to cancel two copies of a high-degree root expression
        // to establish this structural zero. No proximity/tolerance test, and
        // no shortcut for different indices, even if their isolators overlap.
        const int value=(a==b || a==c || b==c)?0:
            int(CGAL::orientation(positions[a],positions[b],positions[c]));
        trace<<"O "<<a<<' '<<b<<' '<<c<<' '<<value<<'\n';return value;
    };
    const auto ordered=[&](std::size_t p,std::size_t a,std::size_t b){
        return compare(p,a,0)*compare(p,b,0)<=0 && compare(p,a,1)*compare(p,b,1)<=0;
    };
    auto pointOrder=selectedPoints;
    detail::rankTwoDiagnostic("chord-sweep:point-order");
    std::sort(pointOrder.begin(),pointOrder.end(),[&](std::size_t a,std::size_t b){
        const int x=compare(a,b,0);if(x)return x<0;
        const int y=compare(a,b,1);return y?y<0:a<b;
    });
    rows<<"POINT_ORDER "<<n;for(const auto p:pointOrder)rows<<' '<<p;rows<<'\n';
    for(std::size_t i=1;i<n;++i){
        const auto a=pointOrder[i-1],b=pointOrder[i];const int x=compare(a,b,0),y=compare(a,b,1);
        if(x==0 && y==0)unresolved("distinct shared point identities coincide: "+points[a].identity+" "+points[b].identity);
        rows<<"DISTINCT_NEIGHBOR "<<a<<' '<<b<<' '<<x<<' '<<y<<'\n';
    }
    // Bounds are selected ACTUAL algebraic endpoints, never floating AABBs.
    std::vector<std::array<std::size_t,4>> bounds(edges.size());
    detail::rankTwoDiagnostic("chord-sweep:edge-bounds");
    for(const auto e:selectedEdges){
        const auto a=edges[e].points[0],b=edges[e].points[1];
        if(a>=points.size() || b>=points.size() || a==b)unresolved("invalid actual shared chord: "+edges[e].identity);
        const bool x=compare(a,b,0)<=0,y=compare(a,b,1)<=0;
        bounds[e]={x?a:b,x?b:a,y?a:b,y?b:a};
        rows<<"EDGE_BOUNDS "<<e;for(const auto p:bounds[e])rows<<' '<<p;rows<<'\n';
    }
    std::uint64_t pe=0,ee=0;
    detail::rankTwoDiagnostic("chord-sweep:point-windows");
    for(const auto e:selectedEdges){
        const auto a=edges[e].points[0],b=edges[e].points[1];
        const auto& box=bounds[e];
        const auto lower=[&](std::size_t bound,bool inclusive){
            std::size_t lo=0,hi=n;
            while(lo<hi){const auto mid=lo+(hi-lo)/2;
                const int sign=compare(pointOrder[mid],bound,0);
                if(sign<0 || (inclusive && sign==0))lo=mid+1;else hi=mid;
            }return lo;
        };
        const auto begin=lower(box[0],false),end=lower(box[1],true);
        rows<<"POINT_WINDOW "<<e<<' '<<begin<<' '<<end<<'\n';
        for(std::size_t i=begin;i<end;++i){
            const auto p=pointOrder[i];if(p==a || p==b)continue;++pe;
            if(compare(p,box[2],1)<0){rows<<"POINT_EDGE "<<p<<' '<<e<<" Y BELOW\n";continue;}
            if(compare(p,box[3],1)>0){rows<<"POINT_EDGE "<<p<<' '<<e<<" Y ABOVE\n";continue;}
            const int o=orientation(a,b,p);const bool on=o==0 && ordered(p,a,b);
            if(on)unresolved("nonincident shared point lies on actual chord: "+points[p].identity+" "+edges[e].identity);
            rows<<"POINT_EDGE "<<p<<' '<<e<<" ORIENTATION "<<o<<" ORDERED "<<on<<'\n';
        }
    }
    auto edgeOrder=selectedEdges;
    detail::rankTwoDiagnostic("chord-sweep:edge-order");
    std::sort(edgeOrder.begin(),edgeOrder.end(),[&](std::size_t a,std::size_t b){
        const int sign=compare(bounds[a][0],bounds[b][0],0);return sign?sign<0:a<b;
    });
    rows<<"EDGE_ORDER "<<m;for(const auto e:edgeOrder)rows<<' '<<e;rows<<'\n';
    detail::rankTwoDiagnostic("chord-sweep:edge-windows");
    for(std::size_t i=0;i<m;++i){
        if(i%32==0)detail::rankTwoDiagnostic(("chord-sweep:edge-window "+std::to_string(i)).c_str());
        const auto first=edgeOrder[i];std::size_t lo=i+1,hi=m;
        while(lo<hi){const auto mid=lo+(hi-lo)/2;
            if(compare(bounds[edgeOrder[mid]][0],bounds[first][1],0)<=0)lo=mid+1;else hi=mid;
        }
        rows<<"EDGE_WINDOW "<<i<<' '<<i+1<<' '<<lo<<'\n';
        for(std::size_t j=i+1;j<lo;++j){
            ++ee;const auto second=edgeOrder[j];
            const auto a=edges[first].points[0],b=edges[first].points[1];
            const auto c=edges[second].points[0],d=edges[second].points[1];
            int common=-1,shared=0;for(const auto p:{a,b})if(p==c || p==d){common=int(p);++shared;}
            if(shared==2)unresolved("duplicate actual shared chords: "+edges[first].identity+" "+edges[second].identity);
            rows<<"EDGE_PAIR "<<first<<' '<<second<<" SHARED "<<common;
            if(compare(bounds[first][3],bounds[second][2],1)<0){rows<<" Y FIRST_BELOW\n";continue;}
            if(compare(bounds[second][3],bounds[first][2],1)<0){rows<<" Y SECOND_BELOW\n";continue;}
            const std::array<int,4> signs{orientation(a,b,c),orientation(a,b,d),orientation(c,d,a),orientation(c,d,b)};
            const bool intersects=(signs[0]*signs[1]<0 && signs[2]*signs[3]<0)
                || (signs[0]==0 && ordered(c,a,b)) || (signs[1]==0 && ordered(d,a,b))
                || (signs[2]==0 && ordered(a,c,d)) || (signs[3]==0 && ordered(b,c,d));
            if(!shared && intersects)unresolved("nonincident actual coarse chords intersect: "+edges[first].identity+" "+edges[second].identity);
            if(shared){const auto p=std::size_t(common),x=a==p?b:a,y=c==p?d:c;
                if(orientation(p,x,y)==0 && (ordered(x,p,y) || ordered(y,p,x)))
                    unresolved("shared-endpoint actual chords overlap: "+edges[first].identity+" "+edges[second].identity);
            }
            rows<<" SIGNS";for(const auto s:signs)rows<<' '<<s;
            rows<<" RESULT "<<(shared?"SHARED_ENDPOINT":"DISJOINT")<<'\n';
        }
    }
    if(pe>out->pointEdges || ee>out->edgePairs)unresolved("sweep duplicated candidate coverage");
    std::ostringstream record;record.imbue(std::locale::classic());
    detail::rankTwoDiagnostic("chord-sweep:serialize");
    record<<(partition.empty()?"rank-two-chord-sweep-v1":"rank-two-partition-chord-sweep-v1")
        <<"\nAUTHORITY "<<out->arrangement<<' '<<out->registry;
    if(!partition.empty()){
        record<<"\nPARTITION_SCOPE "<<partition<<"\nSCOPE_POINTS";
        for(const auto p:selectedPoints)record<<' '<<p;record<<"\nSCOPE_EDGES";
        for(const auto e:selectedEdges)record<<' '<<e;
    }
    record<<"\nMAX_PREDICATE_GROUPS "<<settings.maximumPredicateGroups<<"\nCOUNTS "
        <<out->pointPairs<<' '<<out->pointEdges<<' '<<out->edgePairs<<"\nWORK "<<out->groups<<' '
        <<(n?n-1:0)<<' '<<pe<<' '<<ee<<"\nSHARED_REGISTRY_BEGIN\n"<<registry.canonicalRecord()
        <<"SHARED_REGISTRY_END\n"<<rows.str()<<"PREDICATE_TRACE_BEGIN\n"<<trace.str()
        <<"PREDICATE_TRACE_END\nactual_finite_chord_incidence PASSED\nnative_arc_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();out->hash=hash("rank-two-chord-sweep-certificate-v1\n"+out->record);
    detail::requireExactAlgebraicEnvironment();return RankTwoChordIncidence(std::move(out));
}
} // namespace spin
