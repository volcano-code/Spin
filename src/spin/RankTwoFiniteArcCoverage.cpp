#include "detail/RankTwoFiniteCoverageContinuation.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {
namespace {
[[noreturn]]void invalid(const char* why){throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH",why);}
[[noreturn]]void unresolved(const char* why){throw RankTwoSourceError("ARRANGEMENT_UNRESOLVED",why);}
std::vector<RankTwoMixedArcInterval> wholeIntervals(const RankTwoFinitePartition& p){
    const auto& r=p.registry();std::set<std::string> ids;
    for(const auto i:p.edgeIndices())ids.insert(r.edges().at(i).arcIdentity);
    std::vector<RankTwoMixedArcInterval> out;
    for(const auto& a:r.arcs())if(ids.erase(a.identity)){
        if(!a.finite || !a.ends[0].vertexIdentity || !a.ends[1].vertexIdentity)
            unresolved("included finite domain arc has no actual two finite endpoints");
        out.push_back({a.identity,{*a.ends[0].vertexIdentity,*a.ends[1].vertexIdentity}});
    }
    if(!ids.empty())invalid("finite edge lost its native arc");
    if(out.empty())unresolved("finite domain has no nonempty arc coverage");
    std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){
        return a.arcIdentity<b.arcIdentity;
    });return out;
}
}
struct RankTwoFiniteArcCoverage::Impl {
    RankTwoFinitePartition partition;
    RankTwoMixedArcIntervals refinement;
    std::string record,hash;
    Impl(RankTwoFinitePartition p,RankTwoMixedArcIntervals r):partition(std::move(p)),refinement(std::move(r)){}
};
RankTwoFiniteArcCoverage::RankTwoFiniteArcCoverage(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const std::string& RankTwoFiniteArcCoverage::arrangementIdentity()const{return impl_->partition.arrangementIdentity();}
const RankTwoFinitePartition& RankTwoFiniteArcCoverage::partition()const{return impl_->partition;}
const RankTwoMixedArcIntervals& RankTwoFiniteArcCoverage::refinement()const{return impl_->refinement;}
std::size_t RankTwoFiniteArcCoverage::edgeCount()const{return impl_->partition.edgeIndices().size();}
const std::string& RankTwoFiniteArcCoverage::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoFiniteArcCoverage::scientificHash()const{return impl_->hash;}

RankTwoFiniteArcCoverage certifyRankTwoFiniteArcCoverage(const RankTwoArrangement& a,
    const RankTwoFinitePartition& p,const RankTwoMixedArcIntervals& m){
    detail::requireExactAlgebraicEnvironment();
    const auto& r=p.registry();
    if(p.arrangementIdentity()!=a.scientificHash() || m.registry().arrangementIdentity()!=a.scientificHash()
        || p.registryIdentity()!=m.registry().scientificHash())
        invalid("finite coverage partition/refinement have different actual native authorities");
    const auto intervals=wholeIntervals(p);
    if(intervals.size()!=m.intervals().size())unresolved("selected local intervals do not cover every finite native arc");
    for(std::size_t i=0;i<intervals.size();++i)
        if(intervals[i].arcIdentity!=m.intervals()[i].arcIdentity
            || intervals[i].pointIdentities!=m.intervals()[i].pointIdentities)
            unresolved("local proof stops before an original finite native endpoint");
    std::map<std::string,std::size_t> wanted;
    for(const auto i:p.edgeIndices())if(!wanted.emplace(r.edges().at(i).identity,i).second)
        invalid("duplicate finite native edge identity");
    std::set<std::string> covered;
    for(const auto& proof:m.proofs()){
        const auto at=wanted.find(proof.edgeIdentity());
        if(at==wanted.end() || proof.arrangementIdentity()!=a.scientificHash()
            || proof.registryIdentity()!=r.scientificHash() || !covered.insert(proof.edgeIdentity()).second)
            invalid("regular witness missing, duplicate, foreign or outside actual finite domain");
        const auto& edge=r.edges().at(at->second);
        for(int k=0;k<2;++k){
            const auto& point=r.points().at(edge.points[k]);
            if(proof.pointIdentities()[k]!=point.identity)invalid("strip witness binds a different final point");
            for(int axis=0;axis<2;++axis)
                if(proof.endpointCoordinates()[k][axis].certificate().rootIdentity
                    !=point.coordinates[axis].certificate().rootIdentity)
                    invalid("strip witness coordinates differ from actual final registry");
        }
    }
    for(const auto& id:m.endpointEdges()){
        const auto at=wanted.find(id);
        if(at==wanted.end() || !covered.insert(id).second)invalid("invalid or overlapping finite endpoint witness");
        const auto& edge=r.edges().at(at->second);
        if(edge.endpointNeighborhoodIdentity.empty())invalid("finite endpoint lacks complete neighborhood authority");
        // The private mixed factory already checked the actual event/gate,
        // full-star membership, coordinates and bound against THIS registry.
        const auto found=std::find_if(r.endpointNeighborhoods().begin(),r.endpointNeighborhoods().end(),
            [&](const auto& n){return n.scientificHash()==edge.endpointNeighborhoodIdentity;});
        if(found==r.endpointNeighborhoods().end() || found->arrangementIdentity()!=a.scientificHash())
            invalid("finite endpoint authority is absent or foreign");
    }
    if(covered.size()!=wanted.size())unresolved("some actual finite domain chords have no current local proof");
    if(m.rounds().empty() || m.rounds().back().registryIdentity!=r.scientificHash()
        || m.rounds().back().splitEdges || m.rounds().back().selectedEdges!=wanted.size())
        invalid("finite coverage is not the completed final common refinement round");
    auto out=std::make_shared<RankTwoFiniteArcCoverage::Impl>(p,m);
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-whole-finite-arc-coverage-v1\nAUTHORITY "<<a.scientificHash()<<' '
        <<p.scientificHash()<<' '<<r.scientificHash()<<' '<<m.scientificHash()<<'\n';
    for(const auto& i:intervals)record<<"WHOLE_ARC "<<i.arcIdentity<<' '
        <<i.pointIdentities[0]<<' '<<i.pointIdentities[1]<<'\n';
    for(const auto& e:wanted)record<<"COVERED_EDGE "<<e.second<<' '<<e.first<<' '
        <<(r.edges().at(e.second).endpointNeighborhoodIdentity.empty()?"REGULAR":"ENDPOINT")<<'\n';
    record<<"PARTITION_BEGIN\n"<<p.canonicalRecord()<<"PARTITION_END\nMIXED_BEGIN\n"
        <<m.canonicalRecord()<<"MIXED_END\nall_finite_domain_local_edges PASSED\n"
        <<"global_neighborhood_compatibility NOT_GRANTED\nglobal_isotopy NOT_GRANTED\n"
        <<"profile_precision NOT_GRANTED\nstrict_thickness NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    out->record=record.str();const auto bytes="rank-two-whole-finite-arc-coverage-certificate-v1\n"+out->record;
    verification::Sha256 h;h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=h.finish();
    detail::requireExactAlgebraicEnvironment();return RankTwoFiniteArcCoverage(std::move(out));
}
RankTwoFiniteArcCoverage refineRankTwoFiniteArcCoverage(const RankTwoArrangement& a,
    const RankTwoArcSamplingPlan& plan,const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,
    const CGAL::Gmpq& error,const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    const auto initial=buildRankTwoArcRegistryWithNeighborhoods(a,plan,neighborhoods,limits);
    return refineRankTwoFiniteArcCoverageFromRegistry(a,initial,plan,neighborhoods,error,settings,limits);
}
RankTwoFiniteArcCoverage refineRankTwoFiniteArcCoverageFromRegistry(const RankTwoArrangement& a,
    const RankTwoArcRegistry& initial,const RankTwoArcSamplingPlan& plan,
    const std::vector<RankTwoEndpointNeighborhood>& neighborhoods,const CGAL::Gmpq& error,
    const RankTwoArcRefinementSettings& settings,const RankTwoSourceLimits& limits){
    const auto partition=buildRankTwoFinitePartition(a,initial,limits);
    const auto intervals=wholeIntervals(partition);
    // All selected arcs are WHOLE immutable native arcs. Refinement changes
    // the shared sampling DAG, not D, native labels, holes or exclusions.
    const auto mixed=refineRankTwoMixedArcIntervalsFromRegistry(a,initial,plan,neighborhoods,intervals,error,settings,limits);
    const auto final=buildRankTwoFinitePartition(a,mixed.registry(),limits);
    return certifyRankTwoFiniteArcCoverage(a,final,mixed);
}
} // namespace spin
