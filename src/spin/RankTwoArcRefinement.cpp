#include "detail/RankTwoArcRefinement.h"
#include "detail/RankTwoRegularRegistry.h"
#include "detail/RankTwoMixedArcRefinement.h"
#include "detail/RankTwoMixedArcIntervals.h"
#include "detail/RankTwoRegistryContinuation.h"
#include "detail/RankTwoCarrierKernels.h"
#include "detail/RankTwoDiagnostics.h"
#include "detail/RankTwoStripMemo.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <iomanip>
#include <locale>
#include <optional>
#include <set>
#include <sstream>

namespace spin {
namespace {
using Q=CGAL::Gmpq;
using Z=CGAL::Gmpz;
using AR=detail::AlgebraicReal;
template<class T>std::string text(const T& x){std::ostringstream s;s.imbue(std::locale::classic());s<<x;return s.str();}
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
std::string hash(const std::string& bytes){verification::Sha256 s;
    s.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());return s.finish();}
CORE::BigRat bound(const Q& q){return CORE::BigRat(CORE::BigInt(text(q.numerator())),CORE::BigInt(text(q.denominator())));}
Z floorGrid(AR root,const Q& step,int bits){
    // Floor(root/step) is an exact value property. A prior query may refine a
    // shared native isolator, but cannot choose a different grid cell here.
    for(int i=0;i<=4*bits;++i){
        const Q low=Q(text(root.low()))/step,high=Q(text(root.high()))/step;
        Z lo,hi;mpz_fdiv_q(lo.mpz(),low.numerator().mpz(),low.denominator().mpz());
        mpz_fdiv_q(hi.mpz(),high.numerator().mpz(),high.denominator().mpz());
        if(lo==hi)return lo;
        if(hi==lo+1)return root.compare(bound(Q(hi)*step))==CGAL::SMALLER?lo:hi;
        root.refine();
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","regular box exact floor exceeded isolation operation budget");
}
std::optional<RankTwoRegularEdge> search(const RankTwoArrangement& arrangement,
    const RankTwoArcRegistry& registry,std::size_t edgeIndex,int preferred,
    const Q& error,const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits,
    std::ostringstream& attempts,std::size_t& count,detail::RankTwoRegularStripMemo& memo){
    const auto& edge=registry.edges().at(edgeIndex);
    const auto& a=registry.points()[edge.points[0]],&b=registry.points()[edge.points[1]];
    std::array<AR,4> roots{detail::rankTwoAlgebraicHandle(a.coordinates[0]),
        detail::rankTwoAlgebraicHandle(a.coordinates[1]),detail::rankTwoAlgebraicHandle(b.coordinates[0]),
        detail::rankTwoAlgebraicHandle(b.coordinates[1])};
    const std::array<bool,2> nonzero{roots[0].compare(roots[2])!=CGAL::EQUAL,
        roots[1].compare(roots[3])!=CGAL::EQUAL};
    Q step=settings.initialBoxGridStep;
    for(int level=0;level<=settings.maximumBoxHalvings;++level,step/=2){
        std::array<Z,4> tile;
        for(int i=0;i<4;++i)tile[i]=floorGrid(roots[i],step,limits.maximumIsolationBits);
        for(const int axis:{preferred,1-preferred}){
            if(!nonzero[axis])continue;
            const Q lo=Q(std::min(tile[axis],tile[2+axis])-1)*step;
            const Q hi=Q(std::max(tile[axis],tile[2+axis])+2)*step;
            for(const int padding:{1,4,16}){
                const Q bottom=Q(std::min(tile[1-axis],tile[3-axis])-padding)*step;
                const Q top=Q(std::max(tile[1-axis],tile[3-axis])+1+padding)*step;
                ++count;
                const auto trialPhase="regular-edge:attempt:"+edge.identity+":"+std::to_string(axis)+":"
                    +text(lo)+":"+text(hi)+":"+text(bottom)+":"+text(top);
                detail::rankTwoDiagnostic(trialPhase.c_str());
                try{
                    auto result=tryCertifyRankTwoRegistryEdgeAtIndex(arrangement,registry,edge.identity,edgeIndex,axis,
                        {lo,hi,bottom,top},error,&memo);
                    if(auto* witness=std::get_if<std::string>(&result)){
                        attempts<<"BOX_TRIAL "<<edge.identity<<' '<<axis<<' '<<step<<' '<<padding
                            <<" ARRANGEMENT_UNRESOLVED "<<std::quoted(*witness)<<'\n';
                        continue;
                    }
                    attempts<<"BOX_TRIAL "<<edge.identity<<' '<<axis<<' '<<step<<' '<<padding<<" ACCEPTED\n";
                    return std::move(std::get<RankTwoRegularEdge>(result));
                }catch(const RankTwoSourceError& e){
                    attempts<<"BOX_TRIAL "<<edge.identity<<' '<<axis<<' '<<step<<' '<<padding
                        <<' '<<e.code()<<' '<<std::quoted(e.witness())<<'\n';
                    if(e.code()!="ARRANGEMENT_UNRESOLVED")throw;
                }
            }
        }
    }
    return {};
}
}
struct RankTwoRegularArcRefinement::Impl {
    std::string arc,record,hash;
    RankTwoArcRegistry registry;
    std::vector<RankTwoRegularEdge> proofs;
    std::vector<RankTwoArcRefinementRound> rounds;
    Impl(std::string a,RankTwoArcRegistry r,std::vector<RankTwoRegularEdge> p,
        std::vector<RankTwoArcRefinementRound> history):arc(std::move(a)),registry(std::move(r)),
        proofs(std::move(p)),rounds(std::move(history)){}
};
RankTwoRegularArcRefinement::RankTwoRegularArcRefinement(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoRegularArcRefinement::arcIdentity()const{return impl_->arc;}
const RankTwoArcRegistry& RankTwoRegularArcRefinement::registry()const{return impl_->registry;}
const std::vector<RankTwoRegularEdge>& RankTwoRegularArcRefinement::proofs()const{return impl_->proofs;}
const std::vector<RankTwoArcRefinementRound>& RankTwoRegularArcRefinement::rounds()const{return impl_->rounds;}
const std::string& RankTwoRegularArcRefinement::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoRegularArcRefinement::scientificHash()const{return impl_->hash;}

RankTwoRegularArcRefinement refineRankTwoArcInterval(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const std::string& arcId,const Q& begin,const Q& end,
    const Q& error,const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();
    if(begin<=0 || begin>=end || end>=1 || error<=0 || settings.initialBoxGridStep<=0
        || settings.maximumBoxHalvings<0 || settings.maximumBoxHalvings>64
        || limits.maximumRefinementDepth<0 || limits.maximumRefinementDepth>64
        || limits.maximumIsolationBits<1 || limits.maximumIsolationBits>512)
        invalid("ordered interior fractions, positive accuracy/grid and supported budgets required");
    const auto wanted=requested.find(arcId);
    if(wanted==requested.end() || std::count(wanted->second.begin(),wanted->second.end(),begin)!=1
        || std::count(wanted->second.begin(),wanted->second.end(),end)!=1)
        invalid("complete shared plan must already contain the named interval endpoints exactly once");
    const RankTwoArcRecord* native=nullptr;
    for(const auto& a:arrangementRecords(arrangement).arcs)if(a.identity==arcId)native=&a;
    if(!native || native->ends[0].location!=RankTwoArcEndLocation::Finite
        || native->ends[1].location!=RankTwoArcEndLocation::Finite)
        invalid("regular subdivision requires an actual finite native arc");
    auto plan=requested;std::vector<RankTwoArcRefinementRound> rounds;
    std::ostringstream trials;trials.imbue(std::locale::classic());
    detail::RankTwoRegularStripMemo memo;
    for(int depth=0;depth<=limits.maximumRefinementDepth;++depth){
        const auto phase="arc-refinement:round:"+std::to_string(depth);detail::rankTwoDiagnostic(phase.c_str());
        auto registry=buildRankTwoArcRegistry(arrangement,plan,limits);
        const RankTwoRegistryArc* selected=nullptr;
        for(const auto& a:registry.arcs())if(a.identity==arcId)selected=&a;
        if(!selected)invalid("selected native arc absent after shared publication");
        RankTwoArcRefinementRound round{registry.scientificHash()};
        std::vector<RankTwoRegularEdge> proofs;std::vector<Q> splits;Q next=begin;
        for(const auto index:selected->edges){
            const auto& edge=registry.edges()[index];
            const auto& a=registry.points()[edge.points[0]],&b=registry.points()[edge.points[1]];
            if(a.kind!=RankTwoRegistryPointKind::ArcSample || b.kind!=RankTwoRegistryPointKind::ArcSample)continue;
            if(a.fraction<begin || b.fraction>end)continue;
            if(a.fraction!=next || a.fraction>=b.fraction)invalid("regular interval coverage is not contiguous native order");
            next=b.fraction;++round.selectedEdges;
            auto proof=search(arrangement,registry,index,native->vertical?1:0,error,settings,limits,trials,round.boxTrials,memo);
            if(proof){++round.certifiedEdges;proofs.push_back(std::move(*proof));}
            else{++round.splitEdges;splits.push_back((a.fraction+b.fraction)/2);}
        }
        if(next!=end || !round.selectedEdges)invalid("regular interval has incomplete actual shared-edge coverage");
        rounds.push_back(round);
        if(splits.empty()){
            auto out=std::make_shared<RankTwoRegularArcRefinement::Impl>(arcId,std::move(registry),std::move(proofs),std::move(rounds));
            std::ostringstream record,scientific;record.imbue(std::locale::classic());scientific.imbue(std::locale::classic());
            scientific<<"rank-two-regular-interval-refinement-v1\n"<<arrangement.scientificHash()<<'\n'<<arcId<<'\n'
                <<begin<<' '<<end<<' '<<error<<'\n'<<settings.initialBoxGridStep<<' '<<settings.maximumBoxHalvings
                <<" PAD 1 4 16 NATIVE_AXIS_FIRST\n"<<limits.maximumRefinementDepth<<' '
                <<limits.maximumIsolationBits<<' '<<limits.maximumVertices<<'\n'<<out->registry.scientificHash()<<'\n';
            record<<"REFINEMENT_AUTHORITY "<<arrangement.scientificHash()<<' '<<arcId<<' '
                <<out->registry.scientificHash()<<"\nINTERVAL "<<begin<<' '<<end<<" ERROR "<<error
                <<"\nBOX_POLICY "<<settings.initialBoxGridStep<<' '<<settings.maximumBoxHalvings
                <<" PAD 1 4 16 NATIVE_AXIS_FIRST\nLIMITS "<<limits.maximumRefinementDepth<<' '
                <<limits.maximumIsolationBits<<' '<<limits.maximumVertices<<'\n';
            for(std::size_t i=0;i<out->rounds.size();++i){const auto& r=out->rounds[i];
                std::ostringstream row;row.imbue(std::locale::classic());
                row<<"ROUND "<<i<<' '<<r.registryIdentity<<' '<<r.selectedEdges<<' '
                    <<r.certifiedEdges<<' '<<r.splitEdges<<' '<<r.boxTrials<<'\n';
                record<<row.str();scientific<<row.str();}
            record<<trials.str();scientific<<trials.str();
            record<<"SHARED_REGISTRY_BEGIN\n"<<out->registry.canonicalRecord()<<"SHARED_REGISTRY_END\n";
            for(const auto& p:out->proofs){
                scientific<<p.scientificHash()<<'\n';
                record<<"CERTIFIED_EDGE "<<p.scientificHash()<<'\n'<<p.canonicalRecord();}
            record<<"full_arc_completion NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
            out->record=record.str();out->hash=hash(scientific.str());detail::requireExactAlgebraicEnvironment();
            return RankTwoRegularArcRefinement(std::move(out));
        }
        if(depth==limits.maximumRefinementDepth)throw RankTwoSourceError("RESOURCE_LIMIT",
            "regular interval subdivision exhausted at registry "+registry.scientificHash()
            +" with "+std::to_string(round.splitEdges)+" unproved chords");
        // No accepted proof survives a changed registry. Re-publish exact
        // fractions and reconstruct every final proof in the next round.
        auto& fractions=plan.at(arcId);fractions.insert(fractions.end(),splits.begin(),splits.end());
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","regular interval refinement did not finish");
}
struct RankTwoRegularRegistryRefinement::Impl {
    RankTwoArcRegistry registry;
    std::vector<RankTwoRegularEdge> proofs;
    std::vector<RankTwoArcRefinementRound> rounds;
    std::string record,hash;
    Impl(RankTwoArcRegistry r,std::vector<RankTwoRegularEdge> p,
        std::vector<RankTwoArcRefinementRound> h):registry(std::move(r)),proofs(std::move(p)),rounds(std::move(h)){}
};
RankTwoRegularRegistryRefinement::RankTwoRegularRegistryRefinement(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoArcRegistry& RankTwoRegularRegistryRefinement::registry()const{return impl_->registry;}
const std::vector<RankTwoRegularEdge>& RankTwoRegularRegistryRefinement::proofs()const{return impl_->proofs;}
const std::vector<RankTwoArcRefinementRound>& RankTwoRegularRegistryRefinement::rounds()const{return impl_->rounds;}
const std::string& RankTwoRegularRegistryRefinement::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoRegularRegistryRefinement::scientificHash()const{return impl_->hash;}

RankTwoRegularRegistryRefinement refineRankTwoRegularRegistry(const RankTwoArrangement& arrangement,
    const RankTwoArcSamplingPlan& requested,const Q& error,const RankTwoArcRefinementSettings& settings,
    const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();
    if(error<=0 || settings.initialBoxGridStep<=0 || settings.maximumBoxHalvings<0
        || settings.maximumBoxHalvings>64 || limits.maximumRefinementDepth<0
        || limits.maximumRefinementDepth>64 || limits.maximumIsolationBits<1 || limits.maximumIsolationBits>512)
        invalid("positive regular registry accuracy/grid and supported budgets required");
    const auto& native=arrangementRecords(arrangement);
    std::map<std::string,const RankTwoArcRecord*> arcs;
    for(const auto& a:native.arcs){
        if(a.ends[0].location!=RankTwoArcEndLocation::Finite || a.ends[1].location!=RankTwoArcEndLocation::Finite)
            throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","whole regular registry requires unresolved real infinite arc "+a.identity);
        if(!arcs.emplace(a.identity,&a).second)invalid("duplicate native arc in regular registry authority");
    }
    if(arcs.empty())throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED","no finite regular arcs to certify");
    auto plan=requested;std::vector<RankTwoArcRefinementRound> history;
    std::ostringstream trials;trials.imbue(std::locale::classic());
    detail::RankTwoRegularStripMemo memo;
    for(int depth=0;depth<=limits.maximumRefinementDepth;++depth){
        const auto phase="regular-registry:round:"+std::to_string(depth);detail::rankTwoDiagnostic(phase.c_str());
        auto registry=buildRankTwoArcRegistry(arrangement,plan,limits);
        RankTwoArcRefinementRound round{registry.scientificHash()};
        std::vector<RankTwoRegularEdge> proofs;
        std::vector<std::pair<std::string,Q>> split;
        for(std::size_t index=0;index<registry.edges().size();++index){
            const auto& edge=registry.edges()[index];
            ++round.selectedEdges;
            const auto* arc=arcs.at(edge.arcIdentity);
            const auto& a=registry.points()[edge.points[0]],&b=registry.points()[edge.points[1]];
            const auto fraction=[&](const RankTwoRegistryPoint& point,int end){
                if(point.kind==RankTwoRegistryPointKind::ArcSample && point.arcIdentity==arc->identity)
                    return point.fraction;
                if(point.kind==RankTwoRegistryPointKind::Event && point.eventIdentity==arc->ends[end].vertexIdentity)
                    return Q(end);
                invalid("whole regular refinement cannot replace a simultaneous fan or unknown endpoint");
            };
            const Q begin=fraction(a,0),end=fraction(b,1);
            if(begin>=end)invalid("whole regular edge fraction order is not native increasing");
            auto proof=search(arrangement,registry,index,arc->vertical?1:0,error,settings,limits,trials,round.boxTrials,memo);
            if(proof){++round.certifiedEdges;proofs.push_back(std::move(*proof));}
            else{++round.splitEdges;split.emplace_back(arc->identity,(begin+end)/2);}
        }
        if(round.selectedEdges==0 || round.selectedEdges!=round.certifiedEdges+round.splitEdges)
            invalid("whole regular registry edge accounting failed");
        history.push_back(round);
        if(split.empty()){
            auto out=std::make_shared<RankTwoRegularRegistryRefinement::Impl>(std::move(registry),std::move(proofs),std::move(history));
            std::ostringstream record;record.imbue(std::locale::classic());
            record<<"rank-two-whole-regular-registry-v1\nAUTHORITY "<<arrangement.scientificHash()<<' '
                <<out->registry.scientificHash()<<"\nERROR "<<error<<"\nBOX_POLICY "<<settings.initialBoxGridStep<<' '
                <<settings.maximumBoxHalvings<<" PAD 1 4 16 NATIVE_AXIS_FIRST\nLIMITS "<<limits.maximumRefinementDepth
                <<' '<<limits.maximumIsolationBits<<' '<<limits.maximumVertices<<'\n';
            for(std::size_t i=0;i<out->rounds.size();++i){const auto& r=out->rounds[i];
                record<<"ROUND "<<i<<' '<<r.registryIdentity<<' '<<r.selectedEdges<<' '<<r.certifiedEdges
                    <<' '<<r.splitEdges<<' '<<r.boxTrials<<'\n';}
            record<<trials.str()<<"NATIVE_ARRANGEMENT_BEGIN\n"<<arrangement.canonicalRecord()
                <<"NATIVE_ARRANGEMENT_END\nSHARED_REGISTRY_BEGIN\n"<<out->registry.canonicalRecord()<<"SHARED_REGISTRY_END\n";
            for(const auto& proof:out->proofs)record<<"CERTIFIED_EDGE "<<proof.scientificHash()<<'\n'<<proof.canonicalRecord();
            record<<"all_finite_regular_edges PASSED\nsimultaneous_endpoint_fans NOT_GRANTED\n"
                <<"native_arc_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
            out->record=record.str();out->hash=hash("rank-two-whole-regular-registry-certificate-v1\n"+out->record);
            detail::requireExactAlgebraicEnvironment();return RankTwoRegularRegistryRefinement(std::move(out));
        }
        if(depth==limits.maximumRefinementDepth)throw RankTwoSourceError("RESOURCE_LIMIT",
            "whole regular registry refinement exhausted at "+registry.scientificHash()+" with "
            +std::to_string(round.splitEdges)+" unproved edges; simultaneous endpoint fans may be required");
        // No proof from a changed registry can survive: reconstruct every
        // accepted and split edge against the next common publication.
        for(const auto& s:split)plan.at(s.first).push_back(s.second);
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","whole regular registry refinement did not finish");
}
struct RankTwoMixedArcRefinement::Impl {
    std::string arc,record,hash;
    std::array<std::string,2> ends;
    RankTwoArcRegistry registry;
    std::vector<RankTwoRegularEdge> proofs;
    std::vector<std::string> endpointEdges;
    std::vector<RankTwoMixedArcRefinementRound> rounds;
    Impl(std::string id,std::array<std::string,2> e,RankTwoArcRegistry r,
        std::vector<RankTwoRegularEdge> p,std::vector<std::string> b,
        std::vector<RankTwoMixedArcRefinementRound> h):arc(std::move(id)),ends(std::move(e)),
        registry(std::move(r)),proofs(std::move(p)),endpointEdges(std::move(b)),rounds(std::move(h)){}
};
RankTwoMixedArcRefinement::RankTwoMixedArcRefinement(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoMixedArcRefinement::arcIdentity()const{return impl_->arc;}
const std::array<std::string,2>& RankTwoMixedArcRefinement::endpointIdentities()const{return impl_->ends;}
const RankTwoArcRegistry& RankTwoMixedArcRefinement::registry()const{return impl_->registry;}
const std::vector<RankTwoRegularEdge>& RankTwoMixedArcRefinement::proofs()const{return impl_->proofs;}
const std::vector<std::string>& RankTwoMixedArcRefinement::endpointEdges()const{return impl_->endpointEdges;}
const std::vector<RankTwoMixedArcRefinementRound>& RankTwoMixedArcRefinement::rounds()const{return impl_->rounds;}
const std::string& RankTwoMixedArcRefinement::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoMixedArcRefinement::scientificHash()const{return impl_->hash;}

namespace {
std::vector<std::size_t> mixedIntervalEdges(const RankTwoArcRegistry& r,const std::string& arcId,
    const std::array<std::string,2>& ends){
    const RankTwoRegistryArc* arc=nullptr;
    for(const auto& a:r.arcs())if(a.identity==arcId)arc=&a;
    if(!arc)invalid("mixed interval names no actual native arc");
    if(!arc->finite)throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
        "mixed interval requires certified finite native clipping");
    std::vector<std::string> chain;
    for(const auto i:arc->edges){
        const auto& e=r.edges().at(i);
        if(chain.empty())chain.push_back(r.points().at(e.points[0]).identity);
        else if(chain.back()!=r.points().at(e.points[0]).identity)invalid("native shared chain has a gap");
        chain.push_back(r.points().at(e.points[1]).identity);
    }
    const auto begin=std::find(chain.begin(),chain.end(),ends[0]);
    const auto end=std::find(chain.begin(),chain.end(),ends[1]);
    if(begin==chain.end() || end==chain.end() || begin>=end)
        invalid("mixed interval endpoints must be actual same-arc points in increasing native order");
    return {arc->edges.begin()+(begin-chain.begin()),arc->edges.begin()+(end-chain.begin())};
}
// This binding is against the ACTUAL final registry. Its factory preserved
// each full neighborhood, but no strip may silently replace the bound chord.
std::string mixedEndpointBinding(const RankTwoArcRegistry& r,const RankTwoRegistryEdge& edge,const Q& error){
    const RankTwoEndpointNeighborhood* neighborhood=nullptr;
    for(const auto& n:r.endpointNeighborhoods())if(n.scientificHash()==edge.endpointNeighborhoodIdentity)neighborhood=&n;
    if(!neighborhood)invalid("mixed edge has no complete endpoint authority");
    const RankTwoSharedGateArcCertificate* proof=nullptr;
    if(!edge.endpointFanIdentity.empty()){
        for(const auto& fan:neighborhood->fans())if(fan.scientificHash()==edge.endpointFanIdentity)
            for(const auto& p:fan.orderedArcs())if(p.arcIdentity==edge.arcIdentity)proof=&p;
    }else for(const auto& ray:neighborhood->boundaryRays())if(ray.proof.arcIdentity==edge.arcIdentity)proof=&ray.proof;
    if(!proof || proof->arrangementIdentity!=r.arrangementIdentity())invalid("mixed endpoint arc proof absent or foreign");
    const auto& a=r.points().at(edge.points[0]),&b=r.points().at(edge.points[1]);
    const RankTwoRegistryPoint* gate=nullptr;
    if(a.kind==RankTwoRegistryPointKind::Event && a.eventIdentity==neighborhood->endpointIdentity())gate=&b;
    else if(b.kind==RankTwoRegistryPointKind::Event && b.eventIdentity==neighborhood->endpointIdentity())gate=&a;
    const auto bound=[&]{
        if(!gate)return false;
        if(gate->endpointNeighborhoodIdentity==neighborhood->scientificHash()
            &&gate->endpointArcIdentity==proof->identity&&gate->gateIdentity==proof->gateIdentity
            &&gate->endpointFanIdentity==edge.endpointFanIdentity)return true;
        return std::any_of(gate->additionalEndpointBindings.begin(),gate->additionalEndpointBindings.end(),[&](const auto& b){
            return b.neighborhoodIdentity==neighborhood->scientificHash()&&b.proofIdentity==proof->identity
                &&b.gateIdentity==proof->gateIdentity&&b.fanIdentity==edge.endpointFanIdentity;});
    };
    if(!gate || gate->arcIdentity!=edge.arcIdentity || !bound())
        invalid("mixed endpoint proof does not bind the actual event/gate chord");
    for(int k=0;k<2;++k)if(gate->coordinates[k].expression()!=proof->gateCoordinates[k].expression())
        invalid("mixed endpoint proof gate coordinates changed");
    if(proof->strip.chordDistanceUpperBound>error)throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",
        "complete endpoint chord does not satisfy requested mixed interval error");
    std::ostringstream row;row.imbue(std::locale::classic());
    row<<"ENDPOINT_EDGE "<<edge.identity<<' '<<neighborhood->scientificHash()<<' '
        <<proof->identity<<' '<<proof->strip.chordDistanceUpperBound<<' '<<a.identity<<' '<<b.identity<<'\n';
    return row.str();
}
}
RankTwoMixedArcRefinement refineRankTwoMixedArcInterval(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const std::string& arcId,const std::array<std::string,2>& ends,const Q& error,
    const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();
    if(error<=0 || settings.initialBoxGridStep<=0 || settings.maximumBoxHalvings<0
        || settings.maximumBoxHalvings>64 || limits.maximumRefinementDepth<0
        || limits.maximumRefinementDepth>64 || limits.maximumIsolationBits<1 || limits.maximumIsolationBits>512)
        invalid("positive mixed interval accuracy/grid and supported budgets required");
    const RankTwoArcRecord* native=nullptr;
    for(const auto& arc:arrangementRecords(a).arcs)if(arc.identity==arcId)native=&arc;
    if(!native)invalid("mixed interval arc absent from native authority");
    std::vector<RankTwoArcSubdivision> subdivisions;
    std::vector<RankTwoMixedArcRefinementRound> rounds;
    std::ostringstream trials;trials.imbue(std::locale::classic());
    detail::RankTwoRegularStripMemo memo;
    for(int depth=0;depth<=limits.maximumRefinementDepth;++depth){
        const auto phase="mixed-arc:round:"+std::to_string(depth);detail::rankTwoDiagnostic(phase.c_str());
        auto registry=buildRankTwoArcRegistryWithSubdivisions(a,plan,neighborhoods,subdivisions,limits);
        const auto selected=mixedIntervalEdges(registry,arcId,ends);
        RankTwoMixedArcRefinementRound round{registry.scientificHash()};
        std::vector<RankTwoRegularEdge> proofs;
        std::vector<std::string> endpointEdges;
        std::vector<RankTwoArcSubdivision> split;
        std::ostringstream bindings;bindings.imbue(std::locale::classic());
        for(const auto index:selected){
            const auto& edge=registry.edges().at(index);++round.selectedEdges;
            if(!edge.endpointNeighborhoodIdentity.empty()){
                bindings<<mixedEndpointBinding(registry,edge,error);
                ++round.endpointEdges;endpointEdges.push_back(edge.identity);continue;
            }
            auto proof=search(a,registry,index,native->vertical?1:0,error,settings,limits,trials,round.boxTrials,memo);
            if(proof){++round.regularEdges;proofs.push_back(std::move(*proof));}
            else{
                ++round.splitEdges;
                split.push_back({arcId,{registry.points()[edge.points[0]].identity,
                    registry.points()[edge.points[1]].identity},Q(1,2)});
            }
        }
        if(!round.selectedEdges || round.selectedEdges!=round.endpointEdges+round.regularEdges+round.splitEdges)
            invalid("mixed interval edge accounting incomplete");
        rounds.push_back(round);
        if(split.empty()){
            auto out=std::make_shared<RankTwoMixedArcRefinement::Impl>(arcId,ends,std::move(registry),
                std::move(proofs),std::move(endpointEdges),std::move(rounds));
            std::ostringstream record;record.imbue(std::locale::classic());
            record<<"rank-two-mixed-arc-refinement-v1\nAUTHORITY "<<a.scientificHash()<<' '
                <<out->registry.scientificHash()<<' '<<arcId<<"\nINTERVAL "<<ends[0]<<' '<<ends[1]
                <<"\nERROR "<<error<<"\nBOX_POLICY "<<settings.initialBoxGridStep<<' '
                <<settings.maximumBoxHalvings<<" PAD 1 4 16 NATIVE_AXIS_FIRST\nLIMITS "
                <<limits.maximumRefinementDepth<<' '<<limits.maximumIsolationBits<<' '<<limits.maximumVertices<<'\n';
            for(std::size_t i=0;i<out->rounds.size();++i){const auto& h=out->rounds[i];
                record<<"ROUND "<<i<<' '<<h.registryIdentity<<' '<<h.selectedEdges<<' '
                    <<h.endpointEdges<<' '<<h.regularEdges<<' '<<h.splitEdges<<' '<<h.boxTrials<<'\n';}
            record<<trials.str()<<"SHARED_REGISTRY_BEGIN\n"<<out->registry.canonicalRecord()<<"SHARED_REGISTRY_END\n";
            record<<bindings.str();
            for(const auto& p:out->proofs)record<<"CERTIFIED_EDGE "<<p.scientificHash()<<'\n'<<p.canonicalRecord();
            record<<"selected_interval_error_coverage PASSED\nall_other_edges NOT_GRANTED\n"
                <<"global_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
            out->record=record.str();out->hash=hash("rank-two-mixed-arc-certificate-v1\n"+out->record);
            detail::requireExactAlgebraicEnvironment();return RankTwoMixedArcRefinement(std::move(out));
        }
        if(depth==limits.maximumRefinementDepth)throw RankTwoSourceError("RESOURCE_LIMIT",
            "mixed interval refinement exhausted at "+registry.scientificHash()+" with "
            +std::to_string(round.splitEdges)+" unproved regular edges");
        // Every regular proof is discarded, including the accepted ones. The
        // complete neighborhood binds its unchanged chord to the new registry.
        subdivisions.insert(subdivisions.end(),split.begin(),split.end());
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","mixed interval refinement did not finish");
}
struct RankTwoMixedArcIntervals::Impl {
    std::vector<RankTwoMixedArcInterval> intervals;
    RankTwoArcRegistry registry;
    std::vector<RankTwoRegularEdge> proofs;
    std::vector<std::string> endpointEdges;
    std::vector<RankTwoMixedArcRefinementRound> rounds;
    std::string record,hash;
    Impl(std::vector<RankTwoMixedArcInterval> i,RankTwoArcRegistry r,
        std::vector<RankTwoRegularEdge> p,std::vector<std::string> e,
        std::vector<RankTwoMixedArcRefinementRound> h):intervals(std::move(i)),registry(std::move(r)),
        proofs(std::move(p)),endpointEdges(std::move(e)),rounds(std::move(h)){}
};
RankTwoMixedArcIntervals::RankTwoMixedArcIntervals(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::vector<RankTwoMixedArcInterval>& RankTwoMixedArcIntervals::intervals()const{return impl_->intervals;}
const RankTwoArcRegistry& RankTwoMixedArcIntervals::registry()const{return impl_->registry;}
const std::vector<RankTwoRegularEdge>& RankTwoMixedArcIntervals::proofs()const{return impl_->proofs;}
const std::vector<std::string>& RankTwoMixedArcIntervals::endpointEdges()const{return impl_->endpointEdges;}
const std::vector<RankTwoMixedArcRefinementRound>& RankTwoMixedArcIntervals::rounds()const{return impl_->rounds;}
const std::string& RankTwoMixedArcIntervals::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoMixedArcIntervals::scientificHash()const{return impl_->hash;}

RankTwoMixedArcIntervals refineRankTwoMixedArcIntervals(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const std::vector<RankTwoMixedArcInterval>& requested,const Q& error,
    const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    const auto initial=buildRankTwoArcRegistryWithSubdivisions(a,plan,neighborhoods,{},limits);
    return refineRankTwoMixedArcIntervalsFromRegistry(a,initial,plan,neighborhoods,requested,error,settings,limits);
}
RankTwoMixedArcIntervals refineRankTwoMixedArcIntervalsFromRegistry(const RankTwoArrangement& a,
    const RankTwoArcRegistry& initial,const RankTwoArcSamplingPlan& plan,
    const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const std::vector<RankTwoMixedArcInterval>& requested,const Q& error,
    const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    detail::requireExactAlgebraicEnvironment();
    if(requested.empty() || error<=0 || settings.initialBoxGridStep<=0 || settings.maximumBoxHalvings<0
        || settings.maximumBoxHalvings>64 || limits.maximumRefinementDepth<0
        || limits.maximumRefinementDepth>64 || limits.maximumIsolationBits<1 || limits.maximumIsolationBits>512)
        invalid("nonempty actual intervals, positive accuracy/grid and supported budgets required");
    auto intervals=requested;
    std::sort(intervals.begin(),intervals.end(),[](const auto& left,const auto& right){
        if(left.arcIdentity!=right.arcIdentity)return left.arcIdentity<right.arcIdentity;
        return left.pointIdentities<right.pointIdentities;
    });
    std::map<std::string,const RankTwoArcRecord*> native;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(!native.emplace(arc.identity,&arc).second)invalid("duplicate native arc authority");
    for(const auto& i:intervals)if(!native.count(i.arcIdentity))invalid("selected interval has no native arc authority");
    std::vector<RankTwoArcSubdivision> subdivisions=initial.subdivisions();
    auto authority=initial;
    std::vector<RankTwoMixedArcRefinementRound> rounds;
    std::ostringstream trials;trials.imbue(std::locale::classic());
    detail::RankTwoRegularStripMemo memo;
    for(int depth=0;depth<=limits.maximumRefinementDepth;++depth){
        const auto phase="mixed-intervals:round:"+std::to_string(depth);detail::rankTwoDiagnostic(phase.c_str());
        auto registry=rebuildRankTwoArcRegistryWithSubdivisions(a,authority,plan,neighborhoods,subdivisions,limits);
        std::vector<std::size_t> selected;std::set<std::string> identities;
        detail::rankTwoDiagnostic("mixed-intervals:registry-ready");
        // Resolve ALL intervals against the current common publication before
        // certifying any edge. Duplicate or overlapping coverage is never a
        // successful set union that hides a malformed request.
        for(const auto& i:intervals)for(const auto k:mixedIntervalEdges(registry,i.arcIdentity,i.pointIdentities)){
            if(!identities.emplace(registry.edges().at(k).identity).second)
                invalid("duplicate or edge-overlapping native intervals");
            selected.push_back(k);
        }
        RankTwoMixedArcRefinementRound round{registry.scientificHash()};
        std::vector<RankTwoRegularEdge> proofs;std::vector<std::string> endpointEdges;
        std::vector<RankTwoArcSubdivision> split;
        std::ostringstream bindings;bindings.imbue(std::locale::classic());
        for(const auto k:selected){
            const auto& edge=registry.edges().at(k);++round.selectedEdges;
            const auto edgePhase="mixed-intervals:edge:"+edge.identity;detail::rankTwoDiagnostic(edgePhase.c_str());
            if(!edge.endpointNeighborhoodIdentity.empty()){
                bindings<<mixedEndpointBinding(registry,edge,error);
                ++round.endpointEdges;endpointEdges.push_back(edge.identity);continue;
            }
            auto proof=search(a,registry,k,native.at(edge.arcIdentity)->vertical?1:0,error,
                settings,limits,trials,round.boxTrials,memo);
            if(proof){++round.regularEdges;proofs.push_back(std::move(*proof));}
            else{
                ++round.splitEdges;
                split.push_back({edge.arcIdentity,{registry.points().at(edge.points[0]).identity,
                    registry.points().at(edge.points[1]).identity},Q(1,2)});
            }
            detail::rankTwoDiagnostic("mixed-intervals:edge-ready");
        }
        if(!round.selectedEdges || round.selectedEdges!=round.endpointEdges+round.regularEdges+round.splitEdges)
            invalid("common mixed interval accounting incomplete");
        rounds.push_back(round);
        const auto roundPhase="mixed-intervals:round-counts:"+std::to_string(round.selectedEdges)+":"
            +std::to_string(round.regularEdges)+":"+std::to_string(round.splitEdges);
        detail::rankTwoDiagnostic(roundPhase.c_str());
        if(split.empty()){
            auto out=std::make_shared<RankTwoMixedArcIntervals::Impl>(std::move(intervals),std::move(registry),
                std::move(proofs),std::move(endpointEdges),std::move(rounds));
            std::ostringstream record;record.imbue(std::locale::classic());
            record<<"rank-two-mixed-intervals-refinement-v1\nAUTHORITY "<<a.scientificHash()<<' '
                <<out->registry.scientificHash()<<"\nINTERVAL_COUNT "<<out->intervals.size()<<'\n';
            for(const auto& i:out->intervals)record<<"INTERVAL "<<i.arcIdentity<<' '
                <<i.pointIdentities[0]<<' '<<i.pointIdentities[1]<<'\n';
            record<<"ERROR "<<error<<"\nBOX_POLICY "<<settings.initialBoxGridStep<<' '
                <<settings.maximumBoxHalvings<<" PAD 1 4 16 NATIVE_AXIS_FIRST\nLIMITS "
                <<limits.maximumRefinementDepth<<' '<<limits.maximumIsolationBits<<' '<<limits.maximumVertices<<'\n';
            for(std::size_t j=0;j<out->rounds.size();++j){const auto& h=out->rounds[j];
                record<<"ROUND "<<j<<' '<<h.registryIdentity<<' '<<h.selectedEdges<<' '
                    <<h.endpointEdges<<' '<<h.regularEdges<<' '<<h.splitEdges<<' '<<h.boxTrials<<'\n';}
            record<<trials.str()<<"SHARED_REGISTRY_BEGIN\n"<<out->registry.canonicalRecord()<<"SHARED_REGISTRY_END\n"
                <<bindings.str();
            for(const auto& p:out->proofs)record<<"CERTIFIED_EDGE "<<p.scientificHash()<<'\n'<<p.canonicalRecord();
            record<<"selected_intervals_error_coverage PASSED\nall_other_edges NOT_GRANTED\n"
                <<"global_neighborhood_compatibility NOT_GRANTED\nglobal_isotopy NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
            out->record=record.str();out->hash=hash("rank-two-mixed-intervals-certificate-v1\n"+out->record);
            detail::requireExactAlgebraicEnvironment();return RankTwoMixedArcIntervals(std::move(out));
        }
        if(depth==limits.maximumRefinementDepth)throw RankTwoSourceError("RESOURCE_LIMIT",
            "common mixed intervals exhausted at "+registry.scientificHash()+" with "
            +std::to_string(round.splitEdges)+" unproved ordinary chords");
        // All proofs (also unchanged arcs) are destroyed with this round.
        subdivisions.insert(subdivisions.end(),split.begin(),split.end());
        // Retain only the SAME factory-owned moving coordinates. All strip
        // witnesses are still rebuilt against the next complete registry.
        authority=std::move(registry);
    }
    throw RankTwoSourceError("RESOURCE_LIMIT","common mixed intervals did not finish");
}
} // namespace spin
