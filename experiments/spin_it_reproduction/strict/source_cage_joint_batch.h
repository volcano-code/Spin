#pragma once
// Included by the existing asymmetric fixture after its geometry/oracle helpers.
// This is an experiment adapter; it cannot construct a validated material grant.
#include "spin2014/SourceCageBinaryMaterial.h"
#include "spin2014/SlqpAuditTrace.h"
#include <igl/point_mesh_squared_distance.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <vector>

namespace joint_batch {
namespace fs=std::filesystem;
inline std::string quote(const char* text) {
    std::ostringstream o;o<<'"';
    for(const unsigned char* p=reinterpret_cast<const unsigned char*>(text);*p;++p) {
        if(*p=='"'||*p=='\\')o<<'\\'<<char(*p);
        else if(*p<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(*p)<<std::dec;
        else o<<char(*p);
    }o<<'"';return o.str();
}
inline std::ofstream file(const fs::path& p) {
    std::ofstream o; o.exceptions(std::ios::failbit|std::ios::badbit);
    o.open(p); o.imbue(std::locale::classic());o<<std::setprecision(17);return o;
}
template<class T> void vector(std::ostream& o,const Eigen::MatrixBase<T>& v) {
    o<<'[';for(Eigen::Index i=0;i<v.size();++i) {
        if(i)o<<',';if(std::isfinite(v.derived().coeff(i)))o<<v.derived().coeff(i);else o<<"null";
    }o<<']';
}
inline void obj(const fs::path& p,const Eigen::MatrixXd& v,const Eigen::MatrixXi& f) {
    auto o=file(p);for(Eigen::Index i=0;i<v.rows();++i)o<<"v "<<v.row(i)<<'\n';
    for(Eigen::Index i=0;i<f.rows();++i)o<<"f "<<f(i,0)+1<<' '<<f(i,1)+1<<' '<<f(i,2)+1<<'\n';
}
inline void trace(const fs::path& p,const spin2014::SlqpLbfgsResult& r) {
    auto o=file(p);
    for(const auto& t:r.trace) {
        o<<"{\"iteration\":"<<t.iteration;spin2014::writeSlqpAuditJson(o,t);
        o<<",\"objective\":"<<t.objective<<",\"accepted_step\":"<<t.acceptedStep
         <<",\"line_search_evaluations\":"<<t.lineSearchEvaluations<<"}\n";
    }
}
inline void evaluation(std::ostream& o,const spin2014::SourceCagePaperEvaluation& e) {
    o<<"{\"objective\":"<<e.totalObjective<<",\"paper_objective\":"<<e.paperObjective
     <<",\"constraints\":";vector(o,e.constraints);
    o<<",\"world_moments\":";vector(o,e.world.finalMoments.values());
    o<<",\"spin_moments\":";vector(o,e.spinMoments.values());o<<'}';
}
inline void outerMetrics(std::ostream& o,const Eigen::MatrixXd& original,
    const Eigen::MatrixXd& final,const Eigen::MatrixXi& faces) {
    const Eigen::VectorXd displacement=(final-original).rowwise().norm();
    Eigen::VectorXd a,b;Eigen::VectorXi index;Eigen::MatrixXd closest;
    igl::point_mesh_squared_distance(original,final,faces,a,index,closest);
    igl::point_mesh_squared_distance(final,original,faces,b,index,closest);
    // Same triangle connectivity supplies a continuous barycentric matching:
    // maximum corresponding-vertex displacement bounds both directed distances.
    o<<",\"outer_rms_displacement\":"<<std::sqrt(displacement.squaredNorm()/displacement.size())
     <<",\"outer_max_displacement\":"<<displacement.maxCoeff()
     <<",\"hausdorff_vertex_sample_lower_bound\":"<<std::sqrt(std::max(a.maxCoeff(),b.maxCoeff()))
     <<",\"hausdorff_correspondence_upper_bound\":"<<displacement.maxCoeff();
}
struct PrimitiveWitnessMesh {
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    std::vector<std::size_t> sourceIndices;
};
inline PrimitiveWitnessMesh primitiveWitnessMesh(const spin2014::SourceConformingBoundary& boundary,
    const Eigen::MatrixXd& points,std::size_t primitive) {
    const auto& refs=boundary.primitiveBoundaries().at(primitive);
    std::vector<std::size_t> indices;
    indices.reserve(refs.size()*3);
    for(const auto& ref:refs)
        for(const auto index:boundary.triangles().at(ref.triangle))indices.push_back(index);
    std::sort(indices.begin(),indices.end());
    indices.erase(std::unique(indices.begin(),indices.end()),indices.end());
    PrimitiveWitnessMesh mesh;
    mesh.sourceIndices=indices;
    mesh.vertices.resize(static_cast<Eigen::Index>(indices.size()),3);
    mesh.triangles.resize(static_cast<Eigen::Index>(refs.size()),3);
    for(std::size_t row=0;row<indices.size();++row)mesh.vertices.row(static_cast<Eigen::Index>(row))=points.row(static_cast<Eigen::Index>(indices[row]));
    for(std::size_t row=0;row<refs.size();++row) {
        auto face=boundary.triangles().at(refs[row].triangle);
        if(refs[row].orientation==-1)std::swap(face[1],face[2]);
        else if(refs[row].orientation!=1)throw std::invalid_argument("invalid immutable primitive orientation");
        for(int column=0;column<3;++column)mesh.triangles(static_cast<Eigen::Index>(row),column)=
            static_cast<int>(std::lower_bound(indices.begin(),indices.end(),face[column])-indices.begin());
    }
    return mesh;
}
inline void matrix(std::ostream& o,const Eigen::MatrixXd& value) {
    o<<'[';
    for(Eigen::Index row=0;row<value.rows();++row) {
        if(row)o<<',';vector(o,value.row(row).transpose());
    }
    o<<']';
}
inline void writePrimitiveSubsetWitnesses(const fs::path& directory,
    const spin2014::SourceCageAdaptiveSolverPass& pass,const Eigen::VectorXd& initial) {
    using namespace spin2014;
    if(!pass.solve||!pass.source)return;
    const auto& guarded=pass.solve->facts();
    const auto& source=guarded.source;
    const auto& boundary=pass.source->source().boundary();
    auto witness=file(directory/"primitive_subset_witnesses.jsonl");
    bool emitted=false;
    for(const auto& attempt:guarded.attempts) {
        if(emitted||!attempt.checkCompleted||attempt.callbackRecord>=source.recordCount||
            attempt.callbackRecord>=static_cast<std::size_t>(source.callbackStates.cols()))continue;
        const auto last=std::min(guarded.stages.size(),attempt.firstStage+attempt.stageCount);
        for(std::size_t index=attempt.firstStage;index<last;++index) {
            const auto& stage=guarded.stages[index];
            if(stage.stage!=SourceCageGeometryStage::PrimitiveSubset||
                stage.outcome!=SourceCageGeometryStageOutcome::Rejected||stage.first!=0||
                std::string(stage.reason.data())!="primitive not subset of actual outer material")continue;
            const Eigen::VectorXd state=source.callbackStates.col(static_cast<Eigen::Index>(attempt.callbackRecord));
            const auto evaluated=pass.source->evaluateValues(state);
            const auto primitive=primitiveWitnessMesh(boundary,evaluated.world.deformedPoints,stage.first);
            const auto outer=evaluated.world.deformedPoints.topRows(boundary.outerPointCount());
            const auto outerName="primitive_subset_callback_"+std::to_string(attempt.callbackRecord)+"_outer.obj";
            const auto primitiveName="primitive_subset_callback_"+std::to_string(attempt.callbackRecord)+"_primitive_0.obj";
            obj(directory/outerName,outer,boundary.binding().source().outerFaces());
            obj(directory/primitiveName,primitive.vertices,primitive.triangles);
            const Eigen::MatrixXd allWeights=boundary.binding().pointWeights(0,boundary.binding().points().rows());
            Eigen::MatrixXd weights(static_cast<Eigen::Index>(primitive.sourceIndices.size()),allWeights.cols());
            for(std::size_t row=0;row<primitive.sourceIndices.size();++row)
                weights.row(static_cast<Eigen::Index>(row))=allWeights.row(static_cast<Eigen::Index>(primitive.sourceIndices[row]));
            witness<<"{\"schema\":\"source-cage-primitive-subset-witness-v1\",\"callback\":"<<attempt.callbackRecord
                <<",\"geometry_stage\":"<<index<<",\"primitive\":"<<stage.first
                <<",\"state_dimension\":"<<state.size()<<",\"state\":";vector(witness,state);
            witness<<",\"initial_state_distance\":"<<(state-initial).norm()
                <<",\"state_sha256\":\""<<source.records[attempt.callbackRecord].stateSha256.data()
                <<"\",\"geometry_input_sha256\":\""<<attempt.inputSha256.data()
                <<"\",\"geometry_output_sha256\":\""<<attempt.sourceOutputSha256.data()
                <<"\",\"replayed_evaluation_output_sha256\":\""<<evaluated.outputSha256.data()
                <<"\",\"outer_obj\":\""<<outerName<<"\",\"primitive_obj\":\""<<primitiveName
                <<"\",\"primitive_vertex_count\":"<<primitive.vertices.rows()
                <<",\"primitive_triangle_count\":"<<primitive.triangles.rows()
                <<",\"primitive_source_indices\":[";
            for(std::size_t row=0;row<primitive.sourceIndices.size();++row) {if(row)witness<<',';witness<<primitive.sourceIndices[row];}
            witness<<"],\"weight_rows\":"<<weights.rows()<<",\"weight_columns\":"<<weights.cols()<<",\"weights\":";
            matrix(witness,weights);witness<<",\"reason\":"<<quote(stage.reason.data())<<"}\n";
            emitted=true;
            break;
        }
    }
}
inline int run(const fs::path& directory,bool preflight) {
    using namespace spin2014;
    if(!fs::create_directory(directory))throw std::invalid_argument("run directory must be new");
    const auto boundary=asymmetricBoundary();
    const SourceCageMomentFunctional moments(boundary,1.3,{});
    SourceCagePaperFunctionalConfig config;
    config.objective=Spin2014Objective::TopMassAndLowering;
    config.inertiaWeight=.8;config.centerWeight=.6;config.fillWeight=.03;
    config.shapeWeight=.04;config.displacementWeight=.05;
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.5,.5,.5);
    const SourceCagePaperFunctional functional(moments,config,frame);
    const Eigen::VectorXd beta=Eigen::VectorXd::LinSpaced(functional.betaCount(),.2,.65);
    const Eigen::VectorXd initial=functional.initialState(beta,.17);
    const auto initialEvaluation=functional.evaluate(initial);
    SourceCageAdaptiveSolverSettings settings;
    // Existing defaults are frozen in the preregistration below. No continuation.
    {auto o=file(directory/"input.json");
     o<<"{\"schema\":\"source-cage-joint-batch-v1\",\"classification\":\"SYNTHETIC_REGRESSION\","
      <<"\"handle_count\":"<<functional.handleCount()<<",\"beta_count\":"<<functional.betaCount()
      <<",\"initial_constraint_inf\":"<<initialEvaluation.constraints.lpNorm<Eigen::Infinity>()
      <<",\"functional_sha256\":\""<<functional.configSha256()<<"\",\"upstream_sha256\":\""<<functional.upstreamSha256()
      <<"\",\"boundary_sha256\":\""<<boundary.outputSha256()<<"\",\"bbw_sha256\":\""<<boundary.binding().bbw().outputSha256()
      <<"\",\"initial_state\":";vector(o,initial);
     o<<",\"maximum_passes\":"<<settings.maximumPasses<<",\"controls\":";
     writeSlqpAuditControlsJson(o,settings.solve.source.solver);o<<",\"initial\":";evaluation(o,initialEvaluation);
     o<<",\"binary_constraint_tolerance\":1e-6,\"rounding\":\"beta>=0.5 hollowed\","
      <<"\"hollowing_scope\":\"fixed initial partition, fixed handles; no adaptive completion grant\"}\n";}
    obj(directory/"original_outer.obj",boundary.binding().source().outerVertices(),boundary.binding().source().outerFaces());
    if(preflight)return 0;

    // Same initial state, objective, angle, source and controls, restricted to
    // beta/phi. This baseline is explicitly a fixed-partition diagnostic.
    EqualityBoxProblem hollow;hollow.betaCount=functional.betaCount();hollow.variableCount=functional.handleOffset();
    hollow.evaluate=[&](const Eigen::VectorXd& y) {
        Eigen::VectorXd x=initial;x.head(y.size())=y;
        const auto e=functional.evaluate(x);AdaptiveFillEvaluation r;
        r.objective=e.totalObjective;r.paperObjective=e.paperObjective;
        r.regularizationEnergy=e.totalObjective-e.paperObjective;
        r.gradient=e.jointGradient.head(y.size());r.constraints=e.constraints;
        r.constraintJacobian=e.jointConstraintJacobian.leftCols(y.size());r.finalMoments=e.spinMoments;return r;
    };
    const auto hollowResult=solveSlqpLbfgs(hollow,initial.head(hollow.variableCount),settings.solve.source.solver);
    trace(directory/"hollowing_trace.jsonl",hollowResult);
    {auto o=file(directory/"hollowing.json");Eigen::VectorXd x=initial;x.head(hollow.variableCount)=hollowResult.finalState;
     o<<"{\"scope\":\"FIXED_PARTITION_DIAGNOSTIC\",\"status\":\""<<slqpStatusToken(hollowResult.status)
      <<"\",\"trace_sha256\":\""<<hollowResult.traceHash<<"\",\"final\":";evaluation(o,functional.evaluateValues(x));o<<"}\n";}
    std::cout<<"[BATCH] hollowing="<<slqpStatusToken(hollowResult.status)<<std::endl;
    const auto adaptive=solveSourceCageAdaptive(functional,initial,settings);
    const auto& facts=adaptive.facts();
    for(const auto& pass:facts.passes) {
        const auto prefix="pass_"+std::to_string(pass.passNumber);
        {auto meta=file(directory/(prefix+"_metadata.json"));
         meta<<"{\"pass\":"<<pass.passNumber<<",\"starting_beta_values\":"<<pass.startingBetaValues
             <<",\"cumulative_reserved_bytes\":"<<pass.cumulativeKnownBufferBytes
             <<",\"old_partition\":\""<<pass.oldPartitionSha256<<"\",\"new_partition\":\""<<pass.newPartitionSha256
             <<"\",\"old_state\":\""<<pass.oldStateSha256<<"\",\"solved_state\":\""<<pass.solvedStateSha256
             <<"\",\"new_state\":\""<<pass.newStateSha256<<"\",\"rebind\":";
         if(!pass.rebind)meta<<"null";
         else {
             const auto& r=pass.rebind->facts();
             meta<<"{\"outcome\":"<<int(r.outcome)<<",\"last_stage\":"<<int(r.lastStage)
                 <<",\"scientific_sha256\":\""<<r.scientificSha256<<"\",\"old_beta\":";vector(meta,r.oldBeta);
             meta<<",\"new_beta\":";vector(meta,r.newBeta);
             meta<<",\"moment_difference_available\":"<<r.deformedMomentDifferenceAvailable<<",\"moment_difference\":";
             vector(meta,r.deformedMomentDifference);meta<<",\"stages\":[";
             for(std::size_t i=0;i<r.stages.size();++i) {
                 if(i)meta<<',';const auto& s=r.stages[i];
                 meta<<"{\"stage\":"<<i<<",\"outcome\":"<<int(s.outcome)<<",\"reason\":"<<quote(s.reason.data())
                     <<",\"reason_sha256\":\""<<s.reasonSha256.data()<<"\"}";
             }meta<<"]}";
         }meta<<"}\n";}
        if(!pass.solve)continue;
        const auto& g=pass.solve->facts();
        if(g.source.numericalResultAvailable)trace(directory/(prefix+"_trace.jsonl"),g.source.numerical);
        auto o=file(directory/(prefix+"_callbacks.jsonl"));
        for(std::size_t i=0;i<g.source.recordCount;++i) {
            const auto& r=g.source.records[i];o<<"{\"sequence\":"<<r.sequence<<",\"purpose\":"<<int(r.purpose)
             <<",\"status\":"<<int(r.status)<<",\"state_sha256\":\""<<r.stateSha256.data()
             <<"\",\"output_sha256\":\""<<r.outputSha256.data()<<"\",\"reason_sha256\":\""<<r.reasonSha256.data()
             <<"\",\"reason\":"<<quote(r.reason.data())<<",\"state\":";
            vector(o,g.source.callbackStates.col(i));o<<"}\n";
        }
        auto geometry=file(directory/(prefix+"_geometry.jsonl"));
        for(const auto& r:g.stages)geometry<<"{\"stage\":"<<int(r.stage)<<",\"outcome\":"<<int(r.outcome)
            <<",\"first\":"<<r.first<<",\"second\":"<<r.second<<",\"evidence_sha256\":\""<<r.evidenceSha256.data()
            <<"\",\"reason_sha256\":\""<<r.reasonSha256.data()<<"\",\"reason\":"<<quote(r.reason.data())<<"}\n";
        auto attempts=file(directory/(prefix+"_guard_attempts.jsonl"));
        for(const auto& a:g.attempts)attempts<<"{\"callback\":"<<a.callbackRecord<<",\"first_stage\":"<<a.firstStage
            <<",\"stage_count\":"<<a.stageCount<<",\"completed\":"<<a.checkCompleted
            <<",\"outcome\":"<<int(a.outcome)<<",\"scientific_sha256\":\""<<a.scientificSha256.data()<<"\"}\n";
        writePrimitiveSubsetWitnesses(directory,pass,initial);
    }
    const auto binary=validateSourceCageBinaryMaterial(adaptive);
    const auto& bf=binary.facts();
    if(const auto candidate=binary.candidate()) {
        for(const auto& entry:std::vector<std::pair<std::string,const SourceCageBinaryMesh*>>{
            {"binary_outer.obj",&candidate->outerBoundary()},{"binary_cavity.obj",&candidate->carvedBoundary()},
            {"binary_material.obj",&candidate->materialBoundary()}})
            obj(directory/entry.first,entry.second->vertices,entry.second->triangles);
    }
    auto o=file(directory/"result.json");
    o<<"{\"schema\":\"source-cage-joint-result-v1\",\"classification\":\"SYNTHETIC_REGRESSION\","
     <<"\"adaptive_outcome\":"<<int(facts.outcome)<<",\"adaptive_reason\":"<<quote(facts.terminalReason.data())
     <<",\"adaptive_reason_sha256\":\""<<facts.terminalReasonSha256<<"\",\"cumulative_reserved_bytes\":"<<facts.resources.cumulativeKnownBufferBytes
     <<",\"adaptive_sha256\":\""<<facts.scientificSha256
     <<"\",\"passes\":"<<facts.passes.size()<<",\"binary_status\":\""<<sourceCageBinaryMaterialOutcomeToken(bf.outcome)
     <<"\",\"binary_reason\":"<<quote(bf.reason.data())
     <<",\"binary_sha256\":\""<<bf.scientificSha256<<"\",\"binary_candidate_available\":"<<bool(binary.candidate())
     <<",\"binary_moments\":";vector(o,bf.worldIntegratedMoments);
    o<<",\"binary_normalized_constraints\":";vector(o,bf.normalizedConstraints);
    o<<",\"binary_moments_available\":"<<bf.momentsAvailable<<",\"binary_moment_residuals\":";vector(o,bf.normalizedMomentResiduals);
    o<<",\"binary_geometry\":";
    if(!bf.geometry)o<<"null";
    else {const auto& g=bf.geometry->facts();o<<"{\"outcome\":"<<int(g.outcome)<<",\"scientific_sha256\":\""<<g.scientificSha256<<"\",\"records\":[";
        for(std::size_t i=0;i<g.records.size();++i){if(i)o<<',';const auto& r=g.records[i];
            o<<"{\"stage\":"<<int(r.stage)<<",\"outcome\":"<<int(r.outcome)<<",\"first\":"<<r.first<<",\"second\":"<<r.second
             <<",\"reason\":"<<quote(r.reason.data())<<",\"evidence_sha256\":\""<<r.evidenceSha256.data()<<"\"}";}
        o<<"]}";}
    o<<",\"binary_boundary\":";
    if(!bf.exactBoundary)o<<"null";
    else {const auto& e=bf.exactBoundary->facts();
        o<<"{\"valid\":"<<e.boundaryValid<<",\"closed\":"<<e.closed<<",\"edge_manifold\":"<<e.edgeManifold
         <<",\"vertex_manifold\":"<<e.vertexManifold<<",\"oriented\":"<<e.consistentlyOriented
         <<",\"intersection_completed\":"<<e.intersectionCompleted<<",\"intersection_free\":"<<e.intersectionFree
         <<",\"intersection_error\":"<<quote(e.intersectionError.c_str())<<",\"input_digest\":\""<<e.inputDigest<<"\"}";}
    if(!facts.passes.empty()&&facts.passes.back().solve&&facts.passes.back().solve->facts().source.numericalResultAvailable) {
        const auto& pass=facts.passes.back();const auto& raw=pass.solve->facts().source;
        const auto& f=*pass.source;const auto& x=raw.numerical.finalState;
        const auto e=f.evaluateValues(x);const auto& b=f.source().boundary();
        const auto integrated=independentlyIntegrate(b,e,x.head(f.betaCount()),f.source().density());
        o<<",\"solver_status\":\""<<slqpStatusToken(raw.numerical.status)<<"\",\"trace_sha256\":\""<<raw.numerical.traceHash
         <<"\",\"geometry_status\":"<<int(pass.solve->facts().finalStatus)<<",\"iterations\":"<<raw.numerical.trace.size()
         <<",\"fill_change_last_pass\":"<<(x.head(f.betaCount())-pass.initialState.head(f.betaCount())).norm()
         <<",\"handle_change_from_initial\":"<<(x.tail(4*f.handleCount())-initial.tail(4*f.handleCount())).norm()
         <<",\"final_beta_count\":"<<f.betaCount()
         <<",\"relaxed_independent_moment_error\":"<<(integrated-e.world.finalMoments.values()).cwiseAbs().maxCoeff()/f.source().density()
         <<",\"final_state\":";vector(o,x);o<<",\"final\":";evaluation(o,e);
        const Eigen::MatrixXd outer=e.world.deformedPoints.topRows(b.outerPointCount());
        obj(directory/"last_accepted_outer_DIAGNOSTIC.obj",outer,b.binding().source().outerFaces());
        outerMetrics(o,b.binding().source().outerVertices(),outer,b.binding().source().outerFaces());
    }
    o<<",\"sealed\":false,\"original_paper_case\":false}\n";
    std::cout<<"[BATCH] adaptive_outcome="<<int(facts.outcome)<<" binary="<<sourceCageBinaryMaterialOutcomeToken(bf.outcome)<<std::endl;
    return 0; // Successful capture is not scientific convergence.
}
}
