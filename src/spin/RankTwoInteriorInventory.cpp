#include "detail/RankTwoInteriorInventory.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using C=detail::CertifiedAlgebraicReal;using Q=CGAL::Gmpq;
[[noreturn]]void reject(const char* code,const char* why){throw RankTwoSourceError(code,why);}
bool active(RankTwoRegionKind k){return k==RankTwoRegionKind::Cut||k==RankTwoRegionKind::Full;}
struct Event {const RankTwoVertexRecord* native;C x;};
struct Fiber {const RankTwoArcRecord* native;C y;};
}
struct RankTwoInteriorInventory::Impl {
    std::string authority,record,hash;
    std::size_t slabs=0;
    std::vector<RankTwoMovingInteriorPoint> points;
};
RankTwoInteriorInventory::RankTwoInteriorInventory(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoInteriorInventory::arrangementIdentity()const{return impl_->authority;}
std::size_t RankTwoInteriorInventory::slabCount()const{return impl_->slabs;}
const std::vector<RankTwoMovingInteriorPoint>& RankTwoInteriorInventory::points()const{return impl_->points;}
const std::string& RankTwoInteriorInventory::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoInteriorInventory::scientificHash()const{return impl_->hash;}

RankTwoInteriorInventory enumerateRankTwoMovingInteriors(const RankTwoArrangement& a,
    const RankTwoInteriorInventorySettings& settings){
    detail::requireExactAlgebraicEnvironment();
    if(!settings.maximumCandidateChecks)reject("RESOURCE_LIMIT","positive complete native inventory budget required");
    const int bits=settings.pointProof.maximumIsolationBits;
    if(bits<1||bits>512)reject("INPUT_IDENTITY_MISMATCH","valid native interior isolation budget required");
    std::uint64_t checks=0;
    const auto consume=[&]{if(checks==settings.maximumCandidateChecks)reject("RESOURCE_LIMIT","complete native inventory check budget exhausted");++checks;};
    const auto& n=arrangementRecords(a);
    std::map<std::string,RankTwoRegionKind> kinds;std::set<std::string> wanted,covered;
    for(const auto& r:n.regions){consume();kinds.emplace(r.identity,r.kind);if(active(r.kind)){
        if(r.unbounded)reject("ARRANGEMENT_UNRESOLVED","active native region unexpectedly unbounded");wanted.insert(r.identity);}}
    const auto isActive=[&](const std::string& id){consume();const auto at=kinds.find(id);
        if(at==kinds.end())reject("INPUT_IDENTITY_MISMATCH","native arc region absent from authority");return active(at->second);};
    const auto compare=[&](const C& x,const C& y){consume();return detail::rankTwoAlgebraicHandle(x).compare(detail::rankTwoAlgebraicHandle(y));};
    std::vector<Event> events;
    for(const auto& v:n.vertices){consume();events.push_back({&v,C::reload(v.coordinates[0],bits)});}
    std::sort(events.begin(),events.end(),[&](const Event& x,const Event& y){const auto c=compare(x.x,y.x);
        return c==CGAL::SMALLER||(c==CGAL::EQUAL&&x.native->identity<y.native->identity);});
    std::map<std::string,std::size_t> indices;
    for(std::size_t i=0;i<events.size();++i)indices.emplace(events[i].native->identity,i);
    std::vector<std::size_t> distinct;
    for(std::size_t i=0;i<events.size();++i)if(distinct.empty()||compare(events[distinct.back()].x,events[i].x)!=CGAL::EQUAL)distinct.push_back(i);
    std::vector<RankTwoInteriorPointRequest> requests;
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-native-interior-inventory-v1\nAUTHORITY "<<a.scientificHash()<<' '<<a.sourceIdentity()<<' '<<a.truncationRadius()
        <<"\nPOLICY ALL_ACTIVE_NATIVE_VERTICAL_GAPS FIXED_X_Y_WEIGHTS 1/2 1/2\n";
    for(std::size_t i=0;i<distinct.size();++i)record<<"EVENT "<<i<<' '<<events[distinct[i]].native->identity<<'\n';
    for(std::size_t slab=1;slab<distinct.size();++slab){
        const auto& l=events[distinct[slab-1]];const auto& r=events[distinct[slab]];
        const C x=detail::sampleRankTwoAlgebraicCoordinate(l.x,r.x,Q(1,2),bits);
        std::vector<Fiber> fibers;
        for(const auto& arc:n.arcs){consume();if(arc.vertical)continue;
            if(!isActive(arc.leftRegion)&&!isActive(arc.rightRegion))continue;
            for(const auto& end:arc.ends)if(end.location!=RankTwoArcEndLocation::Finite||!end.vertexIdentity)
                reject("ARRANGEMENT_UNRESOLVED","active native boundary requires certified finite clipping");
            const auto endpoint=[&](int k)->const C&{consume();const auto at=indices.find(*arc.ends[k].vertexIdentity);
                if(at==indices.end())reject("INPUT_IDENTITY_MISMATCH","native boundary event missing");return events[at->second].x;};
            if(compare(endpoint(0),x)!=CGAL::SMALLER||compare(x,endpoint(1))!=CGAL::SMALLER)continue;
            fibers.push_back({&arc,detail::evaluateRankTwoArcProjection(a,arc.identity,x,bits)[1]});
        }
        std::sort(fibers.begin(),fibers.end(),[&](const Fiber& u,const Fiber& v){return compare(u.y,v.y)==CGAL::SMALLER;});
        const auto first=requests.size();
        if(!fibers.empty()&&(isActive(fibers.front().native->rightRegion)||isActive(fibers.back().native->leftRegion)))
            reject("ARRANGEMENT_UNRESOLVED","active native region lacks a bounded outer fiber");
        for(std::size_t i=1;i<fibers.size();++i){const auto& lower=*fibers[i-1].native;const auto& upper=*fibers[i].native;
            if(compare(fibers[i-1].y,fibers[i].y)!=CGAL::SMALLER)reject("ARRANGEMENT_UNRESOLVED","coincident native fibers inside a regular slab");
            const bool lo=isActive(lower.leftRegion),hi=isActive(upper.rightRegion);
            if(!lo&&!hi)continue;
            if(lower.leftRegion!=upper.rightRegion||!lo||!hi)reject("SOURCE_TOPOLOGY_INVALID","consecutive native fibers disagree on active gap ownership");
            RankTwoInteriorPointRequest req;req.regionIdentity=lower.leftRegion;
            req.eventIdentities={l.native->identity,r.native->identity};req.boundaryArcIdentities={lower.identity,upper.identity};
            req.abscissaFraction=req.ordinateFraction=Q(1,2);requests.push_back(req);covered.insert(req.regionIdentity);
            record<<"GAP "<<requests.size()-1<<' '<<slab-1<<' '<<req.regionIdentity<<' '<<lower.identity<<' '<<upper.identity<<'\n';
        }
        record<<"SLAB "<<slab-1<<' '<<l.native->identity<<' '<<r.native->identity<<' '<<requests.size()-first<<'\n';
    }
    if(wanted!=covered)reject("ARRANGEMENT_UNRESOLVED","not every active native region has a complete interior inventory");
    if(requests.size()>settings.maximumPoints)reject("RESOURCE_LIMIT","complete native gap inventory exceeds point budget");
    auto out=std::make_shared<RankTwoInteriorInventory::Impl>();out->authority=a.scientificHash();
    out->slabs=distinct.empty()?0:distinct.size()-1;out->points.reserve(requests.size());
    // Freeze ALL definitions before native point construction. Existing point
    // factory independently verifies adjacent events, all roots and signs.
    for(std::size_t i=0;i<requests.size();++i){consume();detail::rankTwoDiagnostic("native-interior-inventory:point-proof");
        out->points.push_back(buildRankTwoMovingInteriorPoint(a,requests[i],settings.pointProof));
        record<<"NODE_BEGIN "<<i<<' '<<out->points.back().scientificHash()<<'\n'<<out->points.back().canonicalRecord()<<"NODE_END\n";}
    record<<"COUNTS EVENTS "<<distinct.size()<<" SLABS "<<out->slabs<<" GAPS "<<requests.size()<<" ACTIVE_REGIONS "<<covered.size()
        <<"\nBUDGET "<<settings.maximumPoints<<' '<<settings.maximumCandidateChecks<<" CHECKS "<<checks
        <<"\ncomplete_native_gap_inventory PASSED\ncarrier_insertion NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-native-interior-inventory-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoInteriorInventory(std::move(out));
}
RankTwoInteriorInventory reloadRankTwoMovingInteriors(const RankTwoArrangement& a,const std::string& saved,
    const RankTwoInteriorInventorySettings& settings){
    const auto fresh=enumerateRankTwoMovingInteriors(a,settings);
    if(fresh.canonicalRecord()!=saved)reject("INPUT_IDENTITY_MISMATCH","saved native inventory differs from actual complete reconstruction");return fresh;
}
} // namespace spin
