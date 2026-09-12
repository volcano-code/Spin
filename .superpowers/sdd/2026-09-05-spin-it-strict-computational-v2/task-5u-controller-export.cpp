#include "spin2014/SourceCageGeometry.h"
#include "source_surface_test_utils.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
template<class T> void array(const T& a) {std::cout<<'[';for(Eigen::Index i=0;i<a.size();++i){if(i)std::cout<<',';std::cout<<a(i);}std::cout<<']';}
template<class T> void matrix(const T& a) {std::cout<<'[';for(Eigen::Index i=0;i<a.rows();++i){if(i)std::cout<<',';array(a.row(i));}std::cout<<']';}
int main() {
    using namespace spin2014;
    std::cout<<std::setprecision(17);
    const auto outer=surface_test::cube(5./16,11./16),design=surface_test::cube(.34,.66);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=3;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,3);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),SourceCageFemConfig{}),SourceCageBbwConfig{});
    if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
    SourceCageBinding binding(domain,bbw,SourceCageBindingConfig{});
    SourceConformingBoundary boundary(binding,SourceConformingBoundaryConfig{});
    SourceCageMomentFunctional material(boundary,1.25,SourceCageMomentFunctionalConfig{});
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopLoweringOnly;
    SourceCagePaperFunctional paper(material,pc);
    if(paper.betaCount()!=8||paper.handleCount()!=26)throw std::runtime_error("unexpected full source dimensions");
    SourceCageGeometrySettings settings;const SourceCageGeometryChecker checker(paper,settings);
    for(int k=0;k<3;++k) {
        const char* id=k==0?"identity":k==1?"affine":"beta_phi";
        Eigen::VectorXd beta(8);for(int i=0;i<8;++i)beta[i]=k==2?.31+.009*i:.13+.017*i;
        auto x=paper.initialState(beta,k==2?.71:.23);
        if(k==1) {
            const auto& cage=bbw.fem().grid().cage();
            for(std::size_t h=0;h<paper.handleCount();++h) {
                x[paper.handleOffset()+4*h]=2.;
                x.segment<3>(paper.handleOffset()+4*h+1)=Eigen::Vector3d(-.125,.25,.5)+cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
            }
        }
        const auto audit=checker.check(x);const auto& f=audit.facts();
        if(f.outcome!=SourceCageGeometryOutcome::ChecksPassed) {
            for(const auto& r:f.records)std::cerr<<int(r.stage)<<' '<<int(r.outcome)<<' '<<r.reason.data()<<'\n';
            throw std::runtime_error(std::string("actual geometry audit did not pass: ")+id);
        }
        const auto e=paper.evaluateValues(x);
        if(f.sourceOutputSha256!=e.outputSha256)throw std::runtime_error("source mismatch");
        std::cout<<"{\"id\":"<<std::quoted(id)<<",\"scope\":\"SYNTHETIC_REGRESSION\",\"feasibility\":\"UNKNOWN\",\"sealed\":false,\"converged\":false"
                 <<",\"status\":"<<std::quoted(sourceCageGeometryOutcomeToken(f.outcome))
                 <<",\"input_hash\":"<<std::quoted(f.inputSha256)<<",\"output_hash\":"<<std::quoted(f.scientificSha256)
                 <<",\"source_output_hash\":"<<std::quoted(f.sourceOutputSha256)
                 <<",\"density\":"<<material.density()<<",\"minimum_clearance\":"<<settings.minimumOuterDomainClearance<<",\"state\":";array(x);
        std::cout<<",\"beta\":";array(beta);std::cout<<",\"points\":";matrix(e.world.deformedPoints);
        std::cout<<",\"rest_points\":";matrix(boundary.points());std::cout<<",\"outer_faces\":";matrix(outer.f);
        std::cout<<",\"world_moments\":";array(e.world.finalMoments.values());std::cout<<",\"primitives\":[";
        for(std::size_t p=0;p<boundary.primitiveBoundaries().size();++p) {
            if(p)std::cout<<',';const SourceCageGeometryRecord* record=nullptr;
            for(const auto& r:f.records)if(r.stage==SourceCageGeometryStage::OuterDomainClearance&&r.first==p) {
                if(record)throw std::runtime_error("duplicate clearance event");record=&r;
            }
            if(!record||record->outcome!=SourceCageGeometryStageOutcome::Passed)throw std::runtime_error("missing passed clearance");
            std::cout<<"{\"id\":"<<p<<",\"owner\":"<<binding.primitives()[p].betaIndex<<",\"stage\":\"PASSED\",\"distance_squared_exact\":"<<std::quoted(record->squaredDistanceExact)<<",\"faces\":[";
            const auto& refs=boundary.primitiveBoundaries()[p];
            for(std::size_t i=0;i<refs.size();++i) {
                if(i)std::cout<<',';auto t=boundary.triangles()[refs[i].triangle];if(refs[i].orientation<0)std::swap(t[1],t[2]);
                std::cout<<'['<<t[0]<<','<<t[1]<<','<<t[2]<<']';
            }std::cout<<"]}";
        }
        std::cout<<"],\"pair_count\":"<<f.pairCount<<",\"pairs_visited\":"<<f.pairsVisited<<",\"events\":[";
        for(std::size_t i=0;i<f.records.size();++i) {
            if(i)std::cout<<',';const auto& r=f.records[i];
            std::cout<<"{\"stage\":"<<int(r.stage)<<",\"outcome\":"<<int(r.outcome)<<",\"first\":"<<r.first<<",\"second\":"<<r.second<<",\"evidence_hash\":"<<std::quoted(r.evidenceSha256.data())<<'}';
        }std::cout<<"]}\n";
    }
}
