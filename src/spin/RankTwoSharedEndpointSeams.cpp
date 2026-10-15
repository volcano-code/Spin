#include "detail/RankTwoSharedEndpointSeams.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;using C=detail::CertifiedAlgebraicReal;
std::string hash(const std::string& s){verification::Sha256 h;h.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());return h.finish();}
[[noreturn]]void invalid(const char* s){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",s);}
[[noreturn]]void unresolved(const char* s){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",s);}
bool same(const std::array<C,2>& a,const std::array<C,2>& b){return a[0].expression()==b[0].expression()&&a[1].expression()==b[1].expression();}
std::string fraction(const Q& q){std::ostringstream o;o.imbue(std::locale::classic());o<<q.numerator()<<'/'<<q.denominator();return o.str();}
std::string binding(const RankTwoSharedEndpointBinding& b){std::ostringstream o;o.imbue(std::locale::classic());
    o<<b.face<<' '<<b.intervalIdentity<<' '<<b.arcIdentity<<' '<<b.endpointIdentity<<' '
     <<b.neighborhoodIdentity<<' '<<b.proofIdentity<<' '<<b.gateIdentity<<' '<<b.anchorIdentity
     <<' '<<b.projectionAxis<<' '<<fraction(b.fraction);return o.str();}
void root(std::ostream& o,const char* axis,const C& c){const auto& r=c.certificate();o<<axis<<" ROOT "<<r.rootIdentity<<" INDEX "<<r.rootIndex<<" POLY ";
    for(std::size_t k=0;k<r.coefficients.size();++k)o<<(k?",":"")<<r.coefficients[k];o<<" LOWER "<<r.lower<<" UPPER "<<r.upper<<'\n';}
}
struct RankTwoSharedEndpointSeams::Impl {
    RankTwoSharedFiniteSeams master;
    std::vector<RankTwoSharedEndpointPoint> points;
    std::vector<RankTwoSharedEndpointBinding> bindings;
    std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>> neighborhoods;
    std::string record,science;
    explicit Impl(const RankTwoSharedFiniteSeams& m):master(m){}
};
RankTwoSharedEndpointSeams::RankTwoSharedEndpointSeams(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoSharedFiniteSeams& RankTwoSharedEndpointSeams::master()const{return impl_->master;}
const std::string& RankTwoSharedEndpointSeams::masterIdentity()const{return impl_->master.scientificHash();}
const std::vector<RankTwoSharedEndpointPoint>& RankTwoSharedEndpointSeams::points()const{return impl_->points;}
const std::vector<RankTwoSharedEndpointBinding>& RankTwoSharedEndpointSeams::bindings()const{return impl_->bindings;}
const std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>>& RankTwoSharedEndpointSeams::neighborhoods()const{return impl_->neighborhoods;}
const std::string& RankTwoSharedEndpointSeams::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoSharedEndpointSeams::scientificHash()const{return impl_->science;}
RankTwoSharedEndpointSeams publishRankTwoSharedEndpointSeams(const PreparedRankTwoSource& s,
    const RankTwoSharedFiniteSeams& m,const std::vector<std::pair<std::size_t,RankTwoEndpointNeighborhood>>& requested,
    const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();const auto& w=wholeArrangementRecords(m.wholeArrangement());
    if(w.global.sourceIdentity()!=s.inputHash())invalid("endpoint publication source differs from master");
    if(limits.maximumIsolationBits<1||limits.maximumIsolationBits>512)invalid("invalid endpoint publication isolation budget");
    if(!limits.maximumVertices||!limits.maximumArrangementCells)throw RankTwoSourceError("RESOURCE_LIMIT","positive endpoint publication budgets required");
    if(requested.empty())invalid("endpoint publication requires actual neighborhoods");
    auto out=std::make_shared<RankTwoSharedEndpointSeams::Impl>(m);out->neighborhoods=requested;
    std::sort(out->neighborhoods.begin(),out->neighborhoods.end(),[](const auto& x,const auto& y){
        return x.first!=y.first?x.first<y.first:x.second.endpointIdentity()<y.second.endpointIdentity();});
    std::set<std::pair<std::size_t,std::string>> seen;
    struct Entry {RankTwoSharedEndpointBinding binding;std::array<C,2> xy;};
    std::map<std::string,std::vector<Entry>> groups;std::size_t checks=0;bool hasSharedFan=false;
    for(const auto& item:out->neighborhoods){const auto face=item.first;const auto& n=item.second;
        if(face>=w.normalFaces.size()||n.arrangementIdentity()!=w.normalFaces[face].scientificHash()
            ||!seen.emplace(face,n.endpointIdentity()).second)invalid("foreign or duplicated endpoint neighborhood");
        const auto publish=[&](const RankTwoSharedGateArcCertificate& p,const RankTwoMovingGate& g){
            const RankTwoSharedFiniteInterval* interval=nullptr;int side=-1;
            for(const auto& i:m.intervals())if(i.finite)for(int k=0;k<2;++k)
                if(i.incidentFaces[k]==face&&i.sourceArcs[k]==p.arcIdentity){
                    if(interval)unresolved("partial native seam interval composition not implemented");interval=&i;side=k;}
            if(!interval)return false;
            if(++checks>limits.maximumArrangementCells)throw RankTwoSourceError("RESOURCE_LIMIT","endpoint publication check budget exhausted");
            const auto chain=m.consumerPoints(interval->identity,side);
            if(chain.size()<2||(m.points()[chain.front()].identity!=p.endpointIdentity&&m.points()[chain.back()].identity!=p.endpointIdentity))
                unresolved("shared gate endpoint is not an actual master interval endpoint");
            const auto& first=m.points()[chain.front()].coordinates;const auto& last=m.points()[chain.back()].coordinates;
            const int axis=first[0].expression()==last[0].expression()?1:0;
            const auto& value=p.gateCoordinates[axis].expression();
            if(!(value>std::min(first[axis].expression(),last[axis].expression())&&value<std::max(first[axis].expression(),last[axis].expression())))
                unresolved("endpoint gate lies outside original finite shared interval");
            // Check the whole straight component, not only an axis coordinate.
            const auto determinant=(last[0].expression()-first[0].expression())*(p.gateCoordinates[1].expression()-first[1].expression())
                -(last[1].expression()-first[1].expression())*(p.gateCoordinates[0].expression()-first[0].expression());
            if(determinant!=0)unresolved("actual certified gate left original shared line");
            if(p.arrangementIdentity!=n.arrangementIdentity()||p.endpointIdentity!=n.endpointIdentity()
                ||p.gateIdentity!=g.identity()||p.strip.chordDistanceUpperBound!=0)
                invalid("actual straight shared gate lost its complete endpoint authority");
            // Definition-pair grouping: interval + endpoint only. Equality is
            // proved below, while both independent anchors remain in identity.
            const auto key=interval->identity+";"+p.endpointIdentity;
            groups[key].push_back({{face,interval->identity,p.arcIdentity,p.endpointIdentity,n.scientificHash(),
                p.identity,g.identity(),g.anchorIdentity(),g.projectionAxis(),g.fraction(),0},p.gateCoordinates});
            return true;
        };
        for(const auto& ray:n.boundaryRays())publish(ray.proof,ray.gate);
        for(const auto& fan:n.fans())for(const auto& p:fan.orderedArcs())
            if(publish(p,fan.gate()))hasSharedFan=true;
        // Actual local non-seam endpoints still belong to the complete finite
        // inventory. Retain their sealed neighborhoods without inventing a
        // common gate; only actual master incidences publish shared points.
    }
    if(groups.empty())unresolved("publication has no actual finite shared endpoint definitions");
    for(auto& row:groups){auto& e=row.second;std::sort(e.begin(),e.end(),[](const auto& a,const auto& b){return a.binding.face<b.binding.face;});
        const auto interval=std::find_if(m.intervals().begin(),m.intervals().end(),[&](const auto& i){return i.identity==e.front().binding.intervalIdentity;});
        if(e.size()!=2||e[0].binding.face==e[1].binding.face||interval==m.intervals().end()
            ||std::set<std::size_t>{e[0].binding.face,e[1].binding.face}!=std::set<std::size_t>{interval->incidentFaces[0],interval->incidentFaces[1]})
            unresolved("common endpoint gate lacks exactly both incident owner definitions");
        if(!same(e[0].xy,e[1].xy))unresolved("independent endpoint gates differ; common certified gate search required");
        if(out->points.size()==limits.maximumVertices)throw RankTwoSourceError("RESOURCE_LIMIT","endpoint gate point budget exceeded");
        const auto id=hash("rank-two-common-endpoint-definition-pair-v1\n"+m.scientificHash()+"\n"+binding(e[0].binding)+"\n"+binding(e[1].binding)+"\n");
        RankTwoSharedEndpointPoint point{id,e[0].binding.intervalIdentity,e[0].binding.endpointIdentity,e[0].xy,{}};
        // Reuse a master sample ONLY when both actual gate definitions are
        // literally its affine dependency. Equal CURRENT roots alone are
        // insufficient: they may have different velocities after source motion.
        for(const auto index:interval->points){const auto& existing=m.points().at(index);
            if(existing.endpoint||!same(existing.coordinates,point.coordinates))continue;
            const int axis=m.points().at(interval->points.front()).coordinates[0].expression()
                ==m.points().at(interval->points.back()).coordinates[0].expression()?1:0;
            for(const auto& owner:e){const auto& b=owner.binding;
                const bool forward=b.endpointIdentity==existing.parentIdentities[0]
                    &&b.anchorIdentity==existing.parentIdentities[1]&&b.fraction==existing.fraction;
                const bool reversed=b.endpointIdentity==existing.parentIdentities[1]
                    &&b.anchorIdentity==existing.parentIdentities[0]&&b.fraction==Q(1)-existing.fraction;
                if(b.projectionAxis!=axis||(!forward&&!reversed))
                    unresolved("coincident master/gate points have different moving definitions");
            }
            if(point.masterAlias)invalid("duplicate existing master dependency at common gate");
            point.identity=existing.identity;point.coordinates=existing.coordinates;
            point.masterAlias=true;point.definitionPairIdentity=id;
        }
        for(auto& x:e){x.binding.pointIndex=out->points.size();point.consumers.push_back(out->bindings.size());out->bindings.push_back(std::move(x.binding));}
        out->points.push_back(std::move(point));
    }
    const bool aliases=std::any_of(out->points.begin(),out->points.end(),[](const auto& p){return p.masterAlias;});
    std::ostringstream o;o.imbue(std::locale::classic());o<<(hasSharedFan?"rank-two-shared-endpoint-seams-v3\n":
        aliases?"rank-two-shared-endpoint-seams-v2\n":"rank-two-shared-endpoint-seams-v1\n")
        <<"MASTER "<<m.scientificHash()<<"\nPOINTS "<<out->points.size()<<'\n';
    for(const auto& p:out->points){o<<"POINT "<<p.identity<<" INTERVAL "<<p.intervalIdentity<<" ENDPOINT "<<p.endpointIdentity<<" CONSUMERS";
        for(const auto k:p.consumers)o<<' '<<k;o<<'\n';root(o,"X",p.coordinates[0]);root(o,"Y",p.coordinates[1]);}
    o<<"BINDINGS "<<out->bindings.size()<<'\n';for(const auto& b:out->bindings)o<<"BINDING "<<binding(b)<<" POINT "<<b.pointIndex<<'\n';
    if(aliases)for(std::size_t index=0;index<out->points.size();++index){const auto& p=out->points[index];if(p.masterAlias){
        o<<"MASTER_ALIAS "<<p.identity<<" DEFINITION_PAIR "<<p.definitionPairIdentity;
        if(hasSharedFan)o<<" GROUP "<<index;o<<'\n';}}
    for(const auto& n:out->neighborhoods)o<<"NEIGHBORHOOD_BEGIN "<<n.first<<' '<<n.second.scientificHash()<<'\n'<<n.second.canonicalRecord()<<"NEIGHBORHOOD_END\n";
    o<<"common_current_point_equality PASSED\nsource_motion NOT_GRANTED\nfull_coverage NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\nstrict_thickness NOT_GRANTED\n";
    out->record=o.str();out->science=hash(out->record);return RankTwoSharedEndpointSeams(std::move(out));
}
RankTwoArcRegistry buildRankTwoArcRegistryWithSharedEndpointSeams(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const RankTwoSharedEndpointSeams& publication,
    std::size_t face,const RankTwoSourceLimits& limits){
    auto imported=makeRankTwoRegistrySeamImport(a,plan,publication.master(),face);
    imported.endpointAuthority_=std::make_shared<const RankTwoSharedEndpointSeams>(publication);
    imported.masterIdentity_=publication.scientificHash();
    std::vector<RankTwoEndpointNeighborhood> neighborhoods;
    for(const auto& n:publication.neighborhoods())if(n.first==face)neighborhoods.push_back(n.second);
    if(neighborhoods.empty())invalid("requested face has no published endpoint authority");
    for(const auto& b:publication.bindings())if(b.face==face){
        const auto& p=publication.points().at(b.pointIndex);
        if(!imported.endpointBindings_.emplace(b.proofIdentity,p.identity).second)
            invalid("duplicate actual endpoint proof in import");
        auto found=imported.chains_.find(b.arcIdentity);
        if(found==imported.chains_.end())invalid("endpoint proof has no imported finite seam");
        if(p.masterAlias){
            const auto existing=std::find_if(found->second.begin(),found->second.end(),[&](const auto& v){return v.identity==p.identity;});
            if(existing==found->second.end()||existing->kind!=RankTwoRegistryPointKind::SharedSeam
                ||existing->arcIdentity!=b.arcIdentity||!same(existing->coordinates,p.coordinates))
                invalid("structural gate alias lost original actual master sample");
            if(existing->endpointNeighborhoodIdentity.empty()){
                existing->gateIdentity=b.gateIdentity;existing->projectionAxis=b.projectionAxis;
                existing->endpointNeighborhoodIdentity=b.neighborhoodIdentity;existing->endpointArcIdentity=b.proofIdentity;
            }else{
                // Publication already proved BOTH definitions are literally
                // this immutable master's parents/axis/fraction. Preserve
                // both endpoint authorities without overwriting or duplicating
                // the actual geometric point. Mere root equality still rejects.
                existing->additionalEndpointBindings.push_back({b.neighborhoodIdentity,b.proofIdentity,b.gateIdentity,{}});
            }
            imported.endpointAliasProofs_.insert(b.proofIdentity);
        }else{
            found->second.push_back({p.identity,{},b.arcIdentity,RankTwoRegistryPointKind::SharedSeam,p.coordinates,
                {b.endpointIdentity,b.anchorIdentity},b.fraction,b.gateIdentity,{},b.projectionAxis,b.neighborhoodIdentity,b.proofIdentity});
            ++imported.endpointCounts_[b.arcIdentity];
        }
    }
    return buildRankTwoArcRegistryImported(a,plan,{},neighborhoods,{},&imported,limits);
}
} // namespace spin
