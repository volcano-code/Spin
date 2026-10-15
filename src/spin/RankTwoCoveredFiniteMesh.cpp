#include "detail/RankTwoCoveredFiniteMesh.h"
#include "detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#include <locale>
#include <set>
#include <sstream>

namespace spin {
struct RankTwoCoveredFiniteMesh::Impl {
    RankTwoFiniteArcCoverage coverage;
    RankTwoChordIncidence incidence;
    RankTwoRegistryRotation rotation;
    RankTwoFiniteTriangulation triangulation;
    std::string record,hash;
    Impl(RankTwoFiniteArcCoverage c,RankTwoChordIncidence i,RankTwoRegistryRotation r,
        RankTwoFiniteTriangulation t):coverage(std::move(c)),incidence(std::move(i)),
        rotation(std::move(r)),triangulation(std::move(t)){}
};
RankTwoCoveredFiniteMesh::RankTwoCoveredFiniteMesh(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
const RankTwoFiniteArcCoverage& RankTwoCoveredFiniteMesh::coverage()const{return impl_->coverage;}
const RankTwoChordIncidence& RankTwoCoveredFiniteMesh::incidence()const{return impl_->incidence;}
const RankTwoRegistryRotation& RankTwoCoveredFiniteMesh::rotation()const{return impl_->rotation;}
const RankTwoFiniteTriangulation& RankTwoCoveredFiniteMesh::triangulation()const{return impl_->triangulation;}
const std::string& RankTwoCoveredFiniteMesh::canonicalRecord()const{return impl_->record;}
const std::string& RankTwoCoveredFiniteMesh::scientificHash()const{return impl_->hash;}
RankTwoCoveredFiniteMesh triangulateRankTwoCoveredFiniteDomain(const RankTwoArrangement& a,
    const RankTwoFiniteArcCoverage& c,const RankTwoChordSweepSettings& ss,
    const RankTwoRegistryRotationSettings& rs,const RankTwoFiniteTriangulationSettings& ts){
    detail::requireExactAlgebraicEnvironment();const auto& p=c.partition();const auto& r=p.registry();
    if(c.arrangementIdentity()!=a.scientificHash() || p.arrangementIdentity()!=a.scientificHash()
        || c.refinement().registry().scientificHash()!=r.scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","covered mesh requires the SAME final owned finite domain");
    if(!ss.maximumPredicateGroups || !rs.maximumRotationSteps || !ts.maximumTriangles || !ts.maximumAdjacencySteps)
        throw RankTwoSourceError("RESOURCE_LIMIT","positive covered-mesh budgets required");
    detail::rankTwoDiagnostic("covered-finite-mesh:incidence-begin");
    const auto incidence=certifyRankTwoChordIncidenceSweep(a,p,ss);
    detail::rankTwoDiagnostic("covered-finite-mesh:rotation-begin");
    const auto rotation=certifyRankTwoRegistryRotation(a,p,incidence,rs);
    detail::rankTwoDiagnostic("covered-finite-mesh:triangulation-begin");
    const auto mesh=triangulateRankTwoFinitePartition(a,p,incidence,rotation,ts);
    // The triangulator already proves the complete directed native boundary
    // chain, Euler holes, exact CCW faces and all owned point usage. Retain
    // the coverage-edge consumer map as an independently reloadable witness.
    std::set<std::pair<std::size_t,std::size_t>> edges;
    for(const auto& t:mesh.triangles())for(std::size_t j=0;j<3;++j)
        edges.insert(std::minmax(t.points[j],t.points[(j+1)%3]));
    std::ostringstream record;record.imbue(std::locale::classic());
    record<<"rank-two-covered-finite-mesh-v1\nAUTHORITY "<<a.scientificHash()<<' '
        <<c.scientificHash()<<' '<<p.scientificHash()<<' '<<r.scientificHash()<<' '
        <<incidence.scientificHash()<<' '<<rotation.scientificHash()<<' '<<mesh.scientificHash()<<'\n';
    for(const auto e:p.edgeIndices()){
        const auto& edge=r.edges().at(e);
        if(!edges.count(std::minmax(edge.points[0],edge.points[1])))
            throw RankTwoSourceError("SOURCE_TOPOLOGY_INVALID","covered native edge omitted by actual finite mesh");
        record<<"MESH_COVERED_EDGE "<<e<<' '<<edge.identity<<' '<<edge.points[0]<<' '<<edge.points[1]<<'\n';
    }
    record<<"COVERAGE_BEGIN\n"<<c.canonicalRecord()<<"COVERAGE_END\nTRIANGULATION_BEGIN\n"
        <<mesh.canonicalRecord()<<"TRIANGULATION_END\nactual_covered_polygon_mesh PASSED\n"
        <<"global_isotopy NOT_GRANTED\nprofile_precision NOT_GRANTED\nstrict_thickness NOT_GRANTED\nfull_carrier_completion NOT_GRANTED\n";
    auto out=std::make_shared<RankTwoCoveredFiniteMesh::Impl>(c,incidence,rotation,mesh);out->record=record.str();
    const auto bytes="rank-two-covered-finite-mesh-certificate-v1\n"+out->record;
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());out->hash=hash.finish();
    detail::rankTwoDiagnostic("covered-finite-mesh:complete");detail::requireExactAlgebraicEnvironment();
    return RankTwoCoveredFiniteMesh(std::move(out));
}
} // namespace spin
