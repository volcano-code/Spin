#include "detail/RankTwoFinitePartition.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
[[noreturn]] void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]] void topology(const char* why){throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",why);}
Q evaluate(const std::array<Q,25>& p,const std::array<Q,2>& w){
    Q value(0);for(int i=4;i>=0;--i){Q y(0);for(int j=4;j>=0;--j)y=y*w[1]+p[5*i+j];value=value*w[0]+y;}return value;
}
int sign(const Q& x){return x<0?-1:x>0?1:0;}
std::string sha(const std::string& text){verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(text.data()),text.size());return h.finish();}
void polynomial(std::ostream& out,const std::array<Q,25>& p){for(const auto& c:p)out<<' '<<c;out<<'\n';}
bool sameArcWalk(const std::vector<RankTwoDirectedArcRecord>& a,const std::vector<RankTwoDirectedArcRecord>& b){
    if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();++i)
        if(a[i].identity!=b[i].identity || a[i].reversed!=b[i].reversed || a[i].atInfinity!=b[i].atInfinity)return false;
    return true;
}
}
struct RankTwoFinitePartition::Impl {
    explicit Impl(RankTwoArcRegistry r):registry(std::move(r)){}
    RankTwoArcRegistry registry;
    std::string authority,record,hash;
    std::vector<std::size_t> regions,edges,points;
    std::vector<RankTwoExcludedRegion> excludedRegions;
    std::vector<std::string> excludedArcs;
};
RankTwoFinitePartition::RankTwoFinitePartition(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoFinitePartition::arrangementIdentity()const{return impl_->authority;}
const std::string& RankTwoFinitePartition::registryIdentity()const{return impl_->registry.scientificHash();}
const RankTwoArcRegistry& RankTwoFinitePartition::registry()const{return impl_->registry;}
const std::vector<std::size_t>& RankTwoFinitePartition::regionIndices()const{return impl_->regions;}
const std::vector<std::size_t>& RankTwoFinitePartition::edgeIndices()const{return impl_->edges;}
const std::vector<std::size_t>& RankTwoFinitePartition::pointIndices()const{return impl_->points;}
const std::vector<RankTwoExcludedRegion>& RankTwoFinitePartition::excludedRegions()const{return impl_->excludedRegions;}
const std::vector<std::string>& RankTwoFinitePartition::excludedArcs()const{return impl_->excludedArcs;}
const std::string& RankTwoFinitePartition::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoFinitePartition::scientificHash()const{return impl_->hash;}

RankTwoFinitePartition buildRankTwoFinitePartition(const RankTwoArrangement& a,
    const RankTwoArcRegistry& r,const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();
    if(r.arrangementIdentity()!=a.scientificHash())invalid("finite domain and registry have different native authorities");
    const auto& n=arrangementRecords(a);
    if(n.regions.size()>limits.maximumArrangementCells || r.points().size()>limits.maximumVertices)
        throw RankTwoSourceError("RESOURCE_LIMIT","complete finite partition exceeds registered cell/point budget");
    const auto& d=n.originalFhD[2];
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(i+j>2 && d[5*i+j]!=0)
        invalid("projection domain is not the authority's bounded quadratic");
    const Q m00=-d[10],m01=-d[6]/2,m11=-d[2],det=m00*m11-m01*m01;
    if(m00<=0 || det<=0)invalid("projection quadratic is not negative definite");
    const std::array<Q,2> center{(m11*d[5]-m01*d[1])/(2*det),(m00*d[1]-m01*d[5])/(2*det)};
    const Q peak=evaluate(d,center);if(peak<=0)invalid("projection ellipse has no nonempty bounded interior");
    if(n.regions.size()!=r.regions().size() || n.arcs.size()!=r.arcs().size())
        topology("registry omitted a native cell or unsampled arc");
    std::map<std::string,const RankTwoRegionRecord*> originals;
    std::map<std::string,std::size_t> indices;
    for(const auto& cell:n.regions)if(!originals.emplace(cell.identity,&cell).second)topology("duplicate native cell");
    auto out=std::make_shared<RankTwoFinitePartition::Impl>(r);out->authority=a.scientificHash();
    std::set<std::size_t> kept;
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-finite-native-partition-v1\nAUTHORITY "<<a.scientificHash()<<' '<<r.scientificHash()<<'\n';
    record<<"PROJECTION_ELLIPSE "<<m00<<' '<<m01<<' '<<m11<<' '<<det<<' '<<center[0]<<' '<<center[1]<<' '<<peak<<'\n';
    for(std::size_t i=0;i<3;++i){record<<"ORIGINAL_POLYNOMIAL "<<i;polynomial(record,n.originalFhD[i]);}
    for(std::size_t i=0;i<n.normalHalfplanes.size();++i){record<<"NORMAL_POLYNOMIAL "<<i;polynomial(record,n.normalHalfplanes[i]);}
    for(std::size_t i=0;i<r.regions().size();++i){
        const auto& cell=r.regions()[i];const auto found=originals.find(cell.identity);
        if(found==originals.end() || !indices.emplace(cell.identity,i).second)topology("unknown or duplicate registry cell");
        const auto& source=*found->second;
        if(cell.kind!=source.kind || cell.unbounded!=source.unbounded || cell.originalSigns!=source.originalSigns
            || cell.normalSigns!=source.normalSigns || cell.normalSigns.size()!=n.normalHalfplanes.size())
            topology("finite registry changed original region labels");
        for(int k=0;k<3;++k)if(!source.originalSigns[k] || sign(evaluate(n.originalFhD[k],source.witness))!=source.originalSigns[k])
            topology("original sign witness cannot be independently recomputed");
        RankTwoExcludedRegion excluded{i,source.originalSigns[2]<0,{}};
        for(std::size_t k=0;k<n.normalHalfplanes.size();++k){
            if(!source.normalSigns[k] || sign(evaluate(n.normalHalfplanes[k],source.witness))!=source.normalSigns[k])
                topology("normal cone witness cannot be independently recomputed");
            if(source.normalSigns[k]<0)excluded.negativeNormalHalfplanes.push_back(k);
        }
        const auto expected=excluded.outsideProjection || !excluded.negativeNormalHalfplanes.empty()?RankTwoRegionKind::Outside:
            source.originalSigns[0]<0?RankTwoRegionKind::Cut:source.originalSigns[1]<0?RankTwoRegionKind::Full:RankTwoRegionKind::Empty;
        if(expected!=source.kind)topology("native label disagrees with original polynomial signs");
        record<<"REGION_WITNESS "<<i<<' '<<cell.identity<<' '<<source.witness[0]<<' '<<source.witness[1]<<'\n';
        if(expected==RankTwoRegionKind::Outside){out->excludedRegions.push_back(excluded);
            record<<"EXCLUDE_REGION "<<i<<' '<<excluded.outsideProjection<<" NORMAL";
            for(auto j:excluded.negativeNormalHalfplanes)record<<' '<<j;record<<'\n';continue;}
        if(cell.unbounded)topology("positive ellipse cell incorrectly remains unbounded");
        if(cell.outerLoops.empty())topology("bounded included cell has no original outer loop");
        kept.insert(i);out->regions.push_back(i);record<<"INCLUDE_REGION "<<i<<' '<<int(cell.kind)<<'\n';
    }
    std::set<std::size_t> edgeSet,pointSet;
    std::map<std::pair<std::size_t,bool>,unsigned> consumers;
    for(const auto i:out->regions){
        const auto& cell=r.regions()[i];const auto& source=*originals.at(cell.identity);
        for(bool inner:{false,true}){
            const auto& loops=inner?cell.innerLoops:cell.outerLoops;
            const auto& nativeLoops=inner?source.innerLoops:source.outerLoops;
            if(loops.size()!=nativeLoops.size())topology("native finite domain loop omitted");
            for(std::size_t l=0;l<loops.size();++l){const auto& loop=loops[l];
                if(!loop.finite || !sameArcWalk(loop.sourceArcs,nativeLoops[l]))
                    topology("included ellipse cell has an altered/infinite native ring");
                const auto vertices=r.loopVertices(cell.identity,inner,l);
                record<<"DOMAIN_LOOP "<<i<<' '<<(inner?"INNER":"OUTER")<<' '<<l;
                for(const auto& ref:loop.edges){
                    const auto& e=r.edges().at(ref.edgeIndex);
                    edgeSet.insert(ref.edgeIndex);pointSet.insert(e.points[0]);pointSet.insert(e.points[1]);
                    if(++consumers[{ref.edgeIndex,ref.reversed}]!=1)topology("duplicate finite seam consumer");
                    record<<' '<<(ref.reversed?"-":"+")<<ref.edgeIndex;
                }
                record<<'\n';if(vertices.size()<2)topology("included native ring has no actual shared boundary");
            }
        }
        for(const auto p:cell.isolatedPoints){pointSet.insert(p);record<<"DOMAIN_ISOLATED "<<i<<' '<<p<<'\n';}
    }
    std::map<std::string,const RankTwoArcRecord*> arcs;
    for(const auto& arc:n.arcs)arcs.emplace(arc.identity,&arc);
    for(const auto& arc:r.arcs()){
        const auto found=arcs.find(arc.identity);if(found==arcs.end())topology("unknown finite registry arc");
        const auto& native=*found->second;
        if(arc.leftRegion!=native.leftRegion || arc.rightRegion!=native.rightRegion)topology("native arc side changed");
        const auto left=indices.at(native.leftRegion),right=indices.at(native.rightRegion);
        const bool l=kept.count(left),rr=kept.count(right);
        if(!l && !rr){out->excludedArcs.push_back(arc.identity);
            record<<"EXCLUDE_ARC "<<arc.identity<<' '<<arc.finite<<' '<<left<<' '<<right<<'\n';continue;}
        if(!arc.finite || native.ends[0].location!=RankTwoArcEndLocation::Finite
            || native.ends[1].location!=RankTwoArcEndLocation::Finite)
            topology("unbounded native branch survived exact ellipse/cone restriction");
        if(arc.edges.empty())topology("included finite native arc has no shared edges");
        if(l!=rr){bool boundary=false;for(const auto& owner:native.owners)
            boundary=boundary || owner.polynomial=="D" || owner.polynomial.rfind("NORMAL_EDGE:",0)==0;
            if(!boundary)topology("domain boundary lacks projection/normal authority");}
        record<<"INCLUDE_ARC "<<arc.identity<<' '<<left<<' '<<right<<'\n';
        for(const auto e:arc.edges){
            if(!edgeSet.count(e) || consumers[{e,false}]!=unsigned(l) || consumers[{e,true}]!=unsigned(rr))
                topology("included finite arc has missing, reversed or duplicated region consumption");
            record<<"DOMAIN_EDGE "<<e<<' '<<l<<' '<<rr<<'\n';
        }
    }
    out->edges.assign(edgeSet.begin(),edgeSet.end());out->points.assign(pointSet.begin(),pointSet.end());
    record<<"DOMAIN_POINTS";for(auto p:out->points)record<<' '<<p;record<<'\n';
    record<<"NATIVE_BEGIN\n"<<a.canonicalRecord()<<"NATIVE_END\nREGISTRY_BEGIN\n"<<r.canonicalRecord()<<"REGISTRY_END\n";
    record<<"finite_native_domain PASSED\nchord_isotopy NOT_GRANTED\ntriangulation NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();out->hash=sha("rank-two-finite-native-partition-certificate-v1\n"+out->record);
    detail::requireExactAlgebraicEnvironment();return RankTwoFinitePartition(out);
}
RankTwoFinitePartition reloadRankTwoFinitePartition(const RankTwoArrangement& a,const RankTwoArcRegistry& r,
    const std::string& saved,const RankTwoSourceLimits& limits){
    auto rebuilt=buildRankTwoFinitePartition(a,r,limits);
    if(saved!=rebuilt.canonicalRecord())invalid("saved finite partition is not the original complete domain extraction");
    return rebuilt;
}
} // namespace spin
