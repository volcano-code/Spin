#include "SourceSpinMaterialBindings.h"
#include "ExactMeshSequenceIngestion.h"
#include "spin2014/SourceMaterialBoundary.h"
#include "spin2014/SourceAdaptiveDomain.h"
#include "spin2014/DiagnosticStageTimer.h"
#include "spin2014/SourceCheckpointCodec.h"
#include "manufacturing/SlicerValidation.h"
#include "spin2014/SlqpAuditTrace.h"
#include "spin2014/SlqpIterationStore.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <cmath>
#include <climits>
#include <filesystem>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace py = pybind11;
using namespace spin2014;
namespace {
double scalar(const py::handle& v) {
    return exact_mesh_ingestion::exactVertexScalar(v,exact_mesh_ingestion::numpyScalarType());
}
void require(bool ok,const char* message) {if(!ok)throw py::value_error(message);}
struct SettingsReader {
    py::dict values;
    std::set<std::string> consumed;
    explicit SettingsReader(const py::dict& v):values(v){}
    py::handle get(const char* key) {
        if(!values.contains(key))throw py::value_error(std::string("missing setting: ")+key);
        consumed.insert(key);return values[key];
    }
    double real(const char* key) {return scalar(get(key));}
    std::size_t count(const char* key,std::size_t lo=1,std::size_t hi=9007199254740992ULL) {
        auto v=get(key); const double number=scalar(v);
        auto integer=exact_mesh_ingestion::integerObject(v,"counts must be integer scalars");
        require(number>=double(lo)&&number<=double(hi),"count outside supported range");
        return py::cast<std::size_t>(integer);
    }
    std::string token(const char* key) {
        auto v=get(key);if(!py::isinstance<py::str>(v))throw py::type_error("setting token must be string");
        return py::cast<std::string>(v);
    }
    void finish() {require(consumed.size()==values.size(),"unknown settings keys");}
};

SourceSurfaceOctreeConfig sourceConfig(const py::dict& value) {
    SettingsReader d(value);SourceSurfaceOctreeConfig s;
    auto origin=d.get("root_origin");
    require(py::isinstance<py::list>(origin)||py::isinstance<py::tuple>(origin),"root_origin requires explicit three scalars");
    auto seq=py::reinterpret_borrow<py::sequence>(origin);require(seq.size()==3,"root_origin requires three scalars");
    for(int a=0;a<3;++a)s.rootOrigin(a)=scalar(seq[a]);
    s.rootSide=d.real("root_side");require(s.rootSide>0,"root_side must be positive");
    s.maximumDepth=int(d.count("maximum_depth",0,20));s.maximumNodes=d.count("maximum_nodes");d.finish();return s;
}
SlqpLbfgsSettings slqpConfig(const py::dict& value) {
    SettingsReader d(value);SlqpLbfgsSettings s;
    require(d.token("algorithm")=="AUDITED_KKT_RESTORATION_V1","unsupported SLQP algorithm");
    s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;
    s.maximumIterations=int(d.count("maximum_iterations",1,INT_MAX));
    s.maximumLineSearchSteps=int(d.count("maximum_line_search_steps",1,INT_MAX));
    s.lbfgsMemory=int(d.count("lbfgs_memory",1,INT_MAX));
    s.equalityTolerance=d.real("equality_tolerance");s.projectedGradientTolerance=d.real("projected_gradient_tolerance");
    s.activeSetTolerance=d.real("active_set_tolerance");s.meritPenalty=d.real("merit_penalty");
    s.lineSearchReduction=d.real("line_search_reduction");s.armijo=d.real("armijo");
    s.audit.restorationTrustBound=d.real("restoration_trust_bound");
    s.audit.lpIterationLimit=int(d.count("lp_iteration_limit",1,INT_MAX));
    s.audit.lpTimeLimit=d.real("lp_time_limit");s.audit.lpFeasibilityTolerance=d.real("lp_feasibility_tolerance");
    if(value.contains("direction_policy")) {
        const auto policy=d.token("direction_policy");
        require(policy=="METRIC_TANGENT_RESCUE_V1" || policy=="FINITE_BOX_QP_V1" || policy=="FINITE_BOX_PRIMAL_V2","unsupported direction policy");
        const bool box=policy!="METRIC_TANGENT_RESCUE_V1";
        s.directionPolicy=policy=="FINITE_BOX_PRIMAL_V2" ? SlqpDirectionPolicy::FiniteBoxPrimalV2 :
            box ? SlqpDirectionPolicy::FiniteBoxQpV1 : SlqpDirectionPolicy::MetricTangentRescueV1;
        const char* faceKey=box?"finite_box_maximum_faces":"metric_tangent_maximum_faces";
        if(value.contains(faceKey) && py::isinstance<py::str>(value[faceKey])) {
            require(policy=="FINITE_BOX_PRIMAL_V2" && d.token(faceKey)=="BOUNDED_VARIABLES",
                    "unsupported face capacity policy");
            s.metricTangent.maximumFaces=-1;
        } else s.metricTangent.maximumFaces=int(d.count(faceKey,1,100000));
        s.metricTangent.maximumIterations=int(d.count(box?"finite_box_maximum_iterations":"metric_tangent_maximum_iterations",1,INT_MAX));
        s.metricTangent.residualTolerance=d.real(box?"finite_box_residual_tolerance":"metric_tangent_residual_tolerance");
    }
    d.finish();validateSlqpDirectionSettings(s);
    validateSlqpAuditSettings(s.audit);
    require(s.equalityTolerance>0&&s.projectedGradientTolerance>0&&s.activeSetTolerance>=0&&s.activeSetTolerance<.5&&s.meritPenalty>0&&s.lineSearchReduction>0&&s.lineSearchReduction<1&&s.armijo>0&&s.armijo<1,"invalid SLQP controls");
    return s;
}
SourceAdaptiveSpinRequest solverConfig(const py::dict& value,int geometryDepth) {
    SettingsReader d(value);SourceAdaptiveSpinRequest r;
    const auto token=d.token("objective");bool found=false;
    for(auto o:{Spin2014Objective::StaticLowering,Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering})
        if(token==spin2014ObjectiveToken(o)){r.objective=o;found=true;}
    require(found,"invalid source objective");require(d.token("formulation")=="PAPER_GIVENS","unsupported formulation");
    // Missing key is explicitly the legacy API, never an implicit paper upgrade.
    if(value.contains("termination_policy")) {
        const auto policy=d.token("termination_policy");
        if(policy=="PAPER_POST_TRANSITION_V1")r.termination=SourceAdaptiveTermination::PaperPostTransitionV1;
        else require(policy=="UNCHANGED_PARTITION_KKT_V1","unsupported termination policy");
    }
    r.initialLevel=int(d.count("initial_level",0,20));r.maximumOptimizationLevel=int(d.count("maximum_optimization_level",0,20));
    r.maximumAdaptivePasses=int(d.count("maximum_adaptive_passes",1,INT_MAX));
    r.maximumLeafCount=d.count("maximum_leaf_count",1,INT_MAX/6);r.maximumOutputCuboids=d.count("maximum_output_cuboids",1,INT_MAX);
    r.initialFill=d.real("initial_fill");r.initialPhi=d.real("initial_phi");r.fractionalEpsilon=d.real("fractional_epsilon");
    r.laplacianWeight=d.real("laplacian_weight");r.inertiaWeight=d.real("inertia_weight");r.centerWeight=d.real("center_weight");
    r.binaryConstraintTolerance=d.real("binary_constraint_tolerance");
    if(value.contains("boundary_density") || value.contains("interior_density")) {
        require(value.contains("boundary_density") && value.contains("interior_density"),
            "dual material solver settings require both boundary_density and interior_density");
        r.boundaryDensity=d.real("boundary_density");r.interiorDensity=d.real("interior_density");
    }
    const auto origin=d.get("target_origin");
    require(py::isinstance<py::list>(origin)||py::isinstance<py::tuple>(origin),"target_origin requires three scalars");
    const auto originSequence=py::reinterpret_borrow<py::sequence>(origin);
    require(originSequence.size()==3,"target_origin requires three scalars");
    for(int i=0;i<3;++i)r.targetFrame.origin(i)=scalar(originSequence[i]);
    const auto axes=d.get("target_axes");
    require(py::isinstance<py::list>(axes)||py::isinstance<py::tuple>(axes),"target_axes requires a 3x3 matrix");
    const auto axesRows=py::reinterpret_borrow<py::sequence>(axes);
    require(axesRows.size()==3,"target_axes requires a 3x3 matrix");
    for(int row=0;row<3;++row) {
        const auto values=py::reinterpret_borrow<py::sequence>(axesRows[row]);
        require(values.size()==3,"target_axes requires a 3x3 matrix");
        for(int col=0;col<3;++col)r.targetFrame.axes(row,col)=scalar(values[col]);
    }
    r.slqp=slqpConfig(py::cast<py::dict>(d.get("slqp")));d.finish();
    require(r.initialLevel<=r.maximumOptimizationLevel&&r.maximumOptimizationLevel<=geometryDepth,"invalid source optimization levels");
    require(r.initialFill>=0&&r.initialFill<=1&&r.fractionalEpsilon>0&&r.fractionalEpsilon<.5&&r.laplacianWeight>=0&&r.inertiaWeight>=0&&r.centerWeight>=0&&r.binaryConstraintTolerance>0&&(r.objective!=Spin2014Objective::StaticLowering||r.initialPhi==0),"invalid source solve controls");return r;
}
SourceMaterialBoundarySettings boundaryConfig(const py::dict& value,double wall) {
    SettingsReader d(value);SourceMaterialBoundarySettings s;s.minimumOuterWall=wall;
    s.unionSettings.maximumInputCuboids=d.count("maximum_input_cuboids");
    s.unionSettings.maximumOutputVertices=d.count("maximum_output_vertices");
    s.unionSettings.maximumOutputTriangles=d.count("maximum_output_triangles");
    s.unionSettings.maximumIntermediateComplexity=d.count("maximum_intermediate_complexity");d.finish();
    validateSourceMaterialBoundarySettings(s);return s;
}

py::object number(double x) {return std::isfinite(x)?py::object(py::float_(x)):py::object(py::none());}
template<class V> py::list numbers(const V& v) {
    py::list out;for(Eigen::Index i=0;i<v.size();++i)out.append(number(v(i)));return out;
}
py::dict boundarySettings(const verification::ExactCuboidUnionSettings& s) {
    py::dict d;
    d["maximum_input_cuboids"]=s.maximumInputCuboids;
    d["maximum_output_vertices"]=s.maximumOutputVertices;
    d["maximum_output_triangles"]=s.maximumOutputTriangles;
    d["maximum_intermediate_complexity"]=s.maximumIntermediateComplexity;
    return d;
}
std::string settingsJson(const SourceSurfaceOctreeConfig& source,
    const SourceAdaptiveSpinRequest& r,const SourceMaterialBoundarySettings& boundary) {
    py::dict d,s,c;
    s["root_origin"]=numbers(source.rootOrigin);s["root_side"]=source.rootSide;
    s["maximum_depth"]=source.maximumDepth;s["maximum_nodes"]=source.maximumNodes;
    d["source_settings"]=s;
    c["objective"]=spin2014ObjectiveToken(r.objective);c["formulation"]="PAPER_GIVENS";
    c["termination_policy"]=sourceAdaptiveTerminationToken(r.termination);
    c["initial_level"]=r.initialLevel;c["maximum_optimization_level"]=r.maximumOptimizationLevel;
    c["maximum_adaptive_passes"]=r.maximumAdaptivePasses;c["maximum_leaf_count"]=r.maximumLeafCount;
    c["maximum_output_cuboids"]=r.maximumOutputCuboids;c["initial_fill"]=r.initialFill;c["initial_phi"]=r.initialPhi;
    c["fractional_epsilon"]=r.fractionalEpsilon;c["laplacian_weight"]=r.laplacianWeight;
    c["inertia_weight"]=r.inertiaWeight;c["center_weight"]=r.centerWeight;c["binary_constraint_tolerance"]=r.binaryConstraintTolerance;
    if(r.boundaryDensity) {c["boundary_density"]=*r.boundaryDensity;c["interior_density"]=*r.interiorDensity;}
    c["target_origin"]=numbers(r.targetFrame.origin);
    py::list axes;for(int row=0;row<3;++row) {py::list values;for(int col=0;col<3;++col)values.append(number(r.targetFrame.axes(row,col)));axes.append(values);}c["target_axes"]=axes;
    std::ostringstream controls;controls.imbue(std::locale::classic());controls<<std::setprecision(17);
    writeSlqpAuditControlsJson(controls,r.slqp);
    const auto json=py::module_::import("json");
    auto slqp=py::cast<py::dict>(json.attr("loads")(controls.str()));
    py::dict engine;
    for(const char* key:{"provenance","threads","parallel","random_seed","solver","presolve"}) {
        engine[key]=slqp[key];slqp.attr("pop")(key);
    }
    slqp["algorithm"]=slqpAlgorithmToken(r.slqp.algorithm);c["slqp"]=slqp;
    d["solver_settings"]=c;d["boundary_settings"]=boundarySettings(boundary.unionSettings);
    d["density"]=r.density;d["minimum_outer_wall"]=boundary.minimumOuterWall;
    d["audit_engine"]=engine;d["authority"]="IMPLEMENTATION_DEFAULT";
    return py::cast<std::string>(json.attr("dumps")(d));
}
py::dict key(const SourcePartitionKey& k) {
    py::dict d;d["level"]=k.level;d["coordinates"]=k.coordinates;return d;
}
py::list leaves(const std::vector<SourceAdaptiveLeaf>& list) {
    py::list out;for(const auto& l:list){auto d=key(l);d["moments"]=numbers(l.moments.values());out.append(d);}return out;
}
py::dict sourceFacts(const SourceSurfaceOctree& source) {
    py::dict out;out["source_hash"]=source.sourceSha256();out["config_hash"]=source.configSha256();
    out["partition_hash"]=source.partitionSha256();out["outer_surface_moments"]=numbers(source.outerSurfaceMoments().values());
    out["design_surface_moments"]=numbers(source.designSurfaceMoments().values());out["carveable_moments"]=numbers(source.carveableMoments().values());
    out["unresolved_boundary_volume"]=number(source.unresolvedBoundaryVolume());const auto& f=source.stats();py::dict stats;
    stats["internal_nodes"]=f.internalNodes;stats["boundary_nodes"]=f.boundaryNodes;stats["full_nodes"]=f.fullNodes;
    stats["empty_nodes"]=f.emptyNodes;stats["maximum_visited_level"]=f.maximumVisitedLevel;stats["levels"]=f.levels;out["stats"]=stats;return out;
}
const char* sourceCellKindToken(SourceSurfaceCellKind kind) {
    switch(kind) {
    case SourceSurfaceCellKind::Internal:return "INTERNAL";
    case SourceSurfaceCellKind::Boundary:return "BOUNDARY";
    case SourceSurfaceCellKind::Full:return "FULL";
    case SourceSurfaceCellKind::Empty:return "EMPTY";
    }
    throw std::logic_error("unknown source surface cell kind");
}
py::list sourceCells(const SourceSurfaceOctree& source) {
    constexpr std::size_t maximumInspectableCells=100000;
    const auto& cells=source.cells();
    if(cells.size()>maximumInspectableCells)
        throw std::length_error("source-domain inspection cell-output limit exceeded");
    py::list out;for(const auto& cell:cells) {
        py::dict record;record["index"]=cell.index;record["level"]=cell.level;record["coordinates"]=cell.coordinates;
        record["kind"]=sourceCellKindToken(cell.kind);record["moments"]=numbers(cell.moments.values());
        record["unresolved_volume"]=number(cell.unresolvedVolume);record["parent"]=cell.parent;
        py::list children;if(cell.firstChild!=SourceSurfaceCell::noIndex)
            for(std::size_t child=0;child<8;++child)children.append(cell.firstChild+child);
        record["children"]=children;out.append(record);
    }
    return out;
}
py::list sparseTriplets(const Eigen::SparseMatrix<double>& matrix) {
    py::list out;for(int column=0;column<matrix.outerSize();++column)
        for(Eigen::SparseMatrix<double>::InnerIterator entry(matrix,column);entry;++entry)
            out.append(py::make_tuple(entry.row(),entry.col(),number(entry.value())));
    return out;
}
py::dict transition(const SourcePartitionTransition& t) {
    py::dict d;d["before_hash"]=t.beforeHash;d["after_hash"]=t.afterHash;
    d["removed_before"]=numbers(t.removedBefore.values());d["removed_after"]=numbers(t.removedAfter.values());
    py::list changes;for(const auto& c:t.changes){py::dict item;item["parent"]=key(c.parent);py::list children;for(const auto& k:c.children)children.append(key(k));item["children"]=children;changes.append(item);}d["changes"]=changes;return d;
}
py::dict kkt(const KktResidual& k) {
    py::dict d;d["valid"]=k.valid;d["status"]=k.status;d["lambda"]=numbers(k.lambda);d["alpha"]=numbers(k.alpha);d["eta"]=numbers(k.eta);
    d["equality_primal"]=number(k.equalityPrimal);d["bound_primal"]=number(k.boundPrimal);d["stationarity"]=number(k.stationarity);
    d["dual_sign"]=number(k.dualSign);d["complementarity"]=number(k.complementarity);d["full"]=number(k.full);return d;
}
std::string vectorDigest(const Eigen::VectorXd& v) {
    std::string bytes;bytes.reserve(sizeof(std::uint64_t)*(std::size_t(v.size())+1));
    const auto append=[&](std::uint64_t value) {for(int byte=0;byte<8;++byte)bytes.push_back(char(value>>(8*byte)));};
    append(static_cast<std::uint64_t>(v.size()));
    for(Eigen::Index i=0;i<v.size();++i) {std::uint64_t bits;std::memcpy(&bits,&v[i],sizeof(bits));append(bits);}
    return manufacturing::computeContentSha256(bytes);
}
py::dict compactKkt(const KktResidual& k) {
    py::dict d;d["valid"]=k.valid;d["status"]=k.status;
    d["lambda_dimension"]=k.lambda.size();d["lambda_sha256"]=vectorDigest(k.lambda);
    d["alpha_dimension"]=k.alpha.size();d["alpha_sha256"]=vectorDigest(k.alpha);
    d["eta_dimension"]=k.eta.size();d["eta_sha256"]=vectorDigest(k.eta);
    d["equality_primal"]=number(k.equalityPrimal);d["bound_primal"]=number(k.boundPrimal);
    d["stationarity"]=number(k.stationarity);d["dual_sign"]=number(k.dualSign);
    d["complementarity"]=number(k.complementarity);d["full"]=number(k.full);return d;
}
std::string traceJson(const SlqpIteration& t) {
    // Existing serializer is authoritative for the full per-iteration audit.
    std::ostringstream o;o.imbue(std::locale::classic());o<<std::setprecision(17)<<"{\"iteration\":"<<t.iteration;
    writeSlqpAuditJson(o,t);
    // These fields were Python float objects in the legacy dictionary API.
    // Keep decimal representation (including -0.0), not JSON integer -0.
    auto field=[&](const char* name,double v){o<<",\""<<name<<"\":";if(std::isfinite(v))o<<std::showpoint<<v<<std::noshowpoint;else o<<"null";};
    field("objective",t.objective);field("equality_infinity_norm",t.equalityInfinityNorm);
    field("projected_gradient_infinity_norm",t.projectedGradientInfinityNorm);field("kkt_infinity_norm",t.kktInfinityNorm);
    field("merit_value",t.meritValue);field("accepted_step",t.acceptedStep);field("phi",t.phi);
    o<<",\"line_search_evaluations\":"<<t.lineSearchEvaluations<<",\"active_variables\":"<<t.activeVariables
     <<",\"lbfgs_history_size\":"<<t.lbfgsHistorySize<<'}';return o.str();
}
py::dict trace(const SlqpIteration& t) {
    return py::cast<py::dict>(py::module_::import("json").attr("loads")(traceJson(t)));
}
py::dict passFacts(const SourceAdaptivePass& p,bool includeTrace) {
        py::dict a;a["index"]=p.index;a["before"]=leaves(p.before);a["after"]=leaves(p.after);
        a["before_hash"]=p.beforeHash;a["after_hash"]=p.afterHash;a["input_beta"]=numbers(p.inputBeta);a["output_beta"]=numbers(p.outputBeta);a["input_phi"]=number(p.inputPhi);
        a["split"]=transition(p.split);a["merge"]=transition(p.merge);a["level_counts"]=p.levelCounts;a["fractional_level_counts"]=p.fractionalLevelCounts;a["fractional_count"]=p.fractionalCount;
        a["output_level_counts"]=p.outputLevelCounts;a["output_fractional_level_counts"]=p.outputFractionalLevelCounts;
        a["output_fractional_count"]=p.outputFractionalCount;a["paper_resolution_complete"]=p.paperResolutionComplete;
        a["kkt_partition_hash"]=p.beforeHash;
        const auto& s=p.solve;a["status"]=slqpStatusToken(s.status);a["diagnostic"]=s.diagnostic;a["trace_hash"]=s.traceHash;
        a["initial_phi"]=number(s.initialPhi);a["final_phi"]=number(s.finalPhi);a["final_state"]=numbers(s.finalState);a["beta"]=numbers(s.beta);
        a["bounded_variable_count"]=s.boundedVariableCount;a["final_kkt"]=kkt(s.finalKkt);a["algorithm"]=slqpAlgorithmToken(s.algorithm);a["formulation"]="PAPER_GIVENS";
        // Initial rejection leaves a default evaluation; mark it unavailable,
        // rather than publishing those default zeros as measured quantities.
        if(slqpTraceCount(s)==0)a["evaluation"]=py::none();
        else {py::dict e;e["objective"]=number(s.evaluation.objective);e["paper_objective"]=number(s.evaluation.paperObjective);e["regularization_energy"]=number(s.evaluation.regularizationEnergy);e["final_moments"]=numbers(s.evaluation.finalMoments.values());e["gradient"]=numbers(s.evaluation.gradient);e["constraints"]=numbers(s.evaluation.constraints);py::list jac;for(int i=0;i<s.evaluation.constraintJacobian.rows();++i)jac.append(numbers(s.evaluation.constraintJacobian.row(i)));e["constraint_jacobian"]=jac;a["evaluation"]=e;}
        a["evaluation_state"]=slqpTraceCount(s)==0?"UNAVAILABLE":"AVAILABLE";
        if(includeTrace){py::list records;forEachSlqpIteration(s,[&](const auto& t){records.append(trace(t));});a["trace"]=records;}
        return a;
}
// Never expose replay.result through solveFacts(): historical G7 traces are
// multi-gigabyte and materializing them into Python would defeat the bounded
// diagnostic.  This is deliberately a terminal-pass summary only.
py::dict postSolveSummary(const SourceAdaptiveSpinResult& r) {
    require(!r.passes.empty(),"post-solve replay has no terminal pass");
    const auto& p=r.passes.back();py::dict d;
    d["completed_passes"]=r.passes.size();d["status"]=sourceAdaptiveSpinStatusToken(r.status);
    d["input_hash"]=r.inputHash;d["output_hash"]=r.outputHash;d["final_partition_hash"]=r.finalPartitionHash;
    d["before_hash"]=p.beforeHash;d["after_hash"]=p.afterHash;
    d["level_counts"]=p.levelCounts;d["output_level_counts"]=p.outputLevelCounts;
    d["fractional_level_counts"]=p.fractionalLevelCounts;d["output_fractional_level_counts"]=p.outputFractionalLevelCounts;
    d["fractional_count"]=p.fractionalCount;d["output_fractional_count"]=p.outputFractionalCount;
    py::dict split,merge;
    auto transitionSummary=[](const SourcePartitionTransition& t) {py::dict x;x["before_hash"]=t.beforeHash;x["after_hash"]=t.afterHash;
        x["change_count"]=t.changes.size();x["removed_before"]=numbers(t.removedBefore.values());x["removed_after"]=numbers(t.removedAfter.values());return x;};
    d["split"]=transitionSummary(p.split);d["merge"]=transitionSummary(p.merge);
    d["terminal_trace_hash"]=p.solve.traceHash;
    return d;
}
py::dict solveFacts(const SourceAdaptiveSpinResult& r,bool includePasses) {
    py::dict d;d["status"]=sourceAdaptiveSpinStatusToken(r.status);d["diagnostic"]=r.diagnostic;
    d["termination_policy"]=sourceAdaptiveTerminationToken(r.termination);
    d["terminal_relaxation_kkt_applies"]=r.terminalRelaxationKktApplies;
    d["input_hash"]=r.inputHash;d["trace_hash"]=r.traceHash;d["output_hash"]=r.outputHash;d["final_partition_hash"]=r.finalPartitionHash;
    d["geometry_depth"]=r.geometryDepth;d["maximum_optimization_level"]=r.maximumOptimizationLevel;d["final_phi"]=number(r.finalPhi);
    d["relaxed_fill"]=numbers(r.relaxedFill);d["binary_fill"]=numbers(r.binaryFill);d["rounding_delta"]=numbers(r.roundingDelta);d["final_partition"]=leaves(r.finalPartition);
    d["infeasibility_certificate"]=py::none();
    if(r.infeasibilityCertificate) {
        const auto& c=*r.infeasibilityCertificate;py::dict certificate;
        certificate["scope"]=c.scope;certificate["constraint_index"]=c.constraintIndex;
        certificate["lower"]=number(c.lower);certificate["upper"]=number(c.upper);
        certificate["tolerance"]=number(c.tolerance);d["infeasibility_certificate"]=certificate;
    }
    if(includePasses){py::list passes;for(const auto& p:r.passes)passes.append(passFacts(p,true));d["passes"]=passes;}
    d["binary_moments"]=r.binaryMoments?py::object(numbers(r.binaryMoments->values())):py::object(py::none());
    d["binary_check"]=py::none();if(r.binaryCheck){const auto& c=*r.binaryCheck;py::dict b;b["passed"]=c.passed;b["diagnostic"]=c.diagnostic;b["independent_moments"]=c.independentMoments?py::object(numbers(c.independentMoments->values())):py::object(py::none());b["constraints"]=numbers(c.constraints);b["normalized_moment_error"]=number(c.normalizedMomentError);b["raw_constraint_residual"]=number(c.rawConstraintResidual);b["normalized_constraint_residual"]=number(c.normalizedConstraintResidual);b["minimum_central_inertia_eigenvalue"]=number(c.minimumCentralInertiaEigenvalue);d["binary_check"]=b;}
    d["binary_cuboids"]=py::none();if(r.binaryCuboids){const auto& b=*r.binaryCuboids;py::dict c;c["source_hash"]=b.sourceHash;c["partition_hash"]=b.partitionHash;c["output_hash"]=b.outputHash;c["maximum_cuboids"]=b.maximumCuboids;py::list cubes;for(const auto& cube:b.cuboids){py::dict x;x["cell"]=key(cube.cell);x["lower"]=numbers(cube.lower);x["upper"]=numbers(cube.upper);cubes.append(x);}c["cuboids"]=cubes;d["binary_cuboids"]=c;}
    return d;
}

#define FIELD(name, member) d[name]=f.member
py::dict budget(const verification::ExactNefBudgetFacts& f) {
    py::dict d;FIELD("definition",definition);FIELD("maximum_complexity",maximumComplexity);FIELD("current_complexity",currentComplexity);FIELD("peak_complexity",peakComplexity);FIELD("preflight_checks",preflightChecks);FIELD("last_requested_complexity",lastRequestedComplexity);FIELD("exceeded",exceeded);FIELD("last_operation",lastOperation);FIELD("observations_by_operation",observationsByOperation);return d;
}
py::dict region(const verification::ExactMaterialRegionFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("backend",backend);FIELD("semantics",semantics);FIELD("input_digest",inputDigest);FIELD("region_digest",regionDigest);FIELD("error",error);FIELD("completed",completed);FIELD("valid",valid);FIELD("explicit_empty",explicitEmpty);FIELD("parents",parents);FIELD("depths",depths);FIELD("original_signs",originalSigns);FIELD("checker_reversed_components",checkerReversedComponents);FIELD("material_components",materialComponents);FIELD("bounded_void_components",boundedVoidComponents);return d;
}
py::dict clearance(const verification::ExactBoundaryClearanceFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("backend",backend);FIELD("semantics",semantics);FIELD("a_digest",aDigest);FIELD("b_digest",bDigest);FIELD("a_input_digest",aInputDigest);FIELD("b_input_digest",bInputDigest);FIELD("minimum_clearance_exact",minimumClearanceExact);FIELD("clearance_digest",clearanceDigest);FIELD("error",error);FIELD("state",state);FIELD("completed",completed);FIELD("a_explicit_empty",aExplicitEmpty);FIELD("b_explicit_empty",bExplicitEmpty);FIELD("applicable",applicable);FIELD("at_least_minimum",atLeastMinimum);FIELD("contact",contact);FIELD("squared_distance_exact",squaredDistanceExact);FIELD("witness_faces",witnessFaces);return d;
}
py::dict partition(const verification::ExactMaterialPartitionFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("backend",backend);FIELD("semantics",semantics);FIELD("empty_protocol",emptyProtocol);FIELD("whole_digest",wholeDigest);FIELD("partition_digest",partitionDigest);FIELD("error",error);FIELD("part_digests",partDigests);FIELD("completed",completed);FIELD("whole_explicit_empty",wholeExplicitEmpty);FIELD("part_explicit_empty",partExplicitEmpty);FIELD("union_equals_whole",unionEqualsWhole);FIELD("no_positive_volume_overlap",noPositiveVolumeOverlap);FIELD("missing_whole_content",missingWholeContent);FIELD("excess_union_content",excessUnionContent);FIELD("complete_partition",completePartition);FIELD("part_subset_whole",partSubsetWhole);FIELD("positive_volume_overlap_pairs",positiveVolumeOverlapPairs);return d;
}
py::dict relation(const verification::ExactMaterialRelationFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("semantics",semantics);FIELD("a_digest",aDigest);FIELD("b_digest",bDigest);FIELD("relation_digest",relationDigest);FIELD("error",error);FIELD("completed",completed);FIELD("a_explicit_empty",aExplicitEmpty);FIELD("b_explicit_empty",bExplicitEmpty);FIELD("intersection_nonempty",intersectionNonempty);FIELD("positive_volume_overlap",positiveVolumeOverlap);FIELD("a_subset_b",aSubsetB);FIELD("b_subset_a",bSubsetA);FIELD("equal",equal);return d;
}
py::dict boundary(const verification::ExactBoundaryFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("backend",backend);FIELD("input_digest",inputDigest);FIELD("orientation",orientation);FIELD("intersection_error",intersectionError);FIELD("boundary_valid",boundaryValid);FIELD("closed",closed);FIELD("edge_manifold",edgeManifold);FIELD("vertex_manifold",vertexManifold);FIELD("consistently_oriented",consistentlyOriented);FIELD("component_orientations_valid",componentOrientationsValid);FIELD("intersection_completed",intersectionCompleted);FIELD("intersection_free",intersectionFree);FIELD("exact_pair_tests",exactPairTests);FIELD("unreferenced_vertices",unreferencedVertices);FIELD("degenerate_faces",degenerateFaces);FIELD("non_manifold_vertices",nonManifoldVertices);FIELD("boundary_edges",boundaryEdges);FIELD("non_manifold_edges",nonManifoldEdges);FIELD("inconsistent_edges",inconsistentEdges);FIELD("intersecting_pairs",intersectingPairs);FIELD("face_components",faceComponents);FIELD("component_volume_signs",componentVolumeSigns);FIELD("component_signed_volumes_exact",componentSignedVolumesExact);return d;
}
py::dict unionFacts(const verification::ExactCuboidUnionFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("backend",backend);FIELD("complexity_definition",complexityDefinition);FIELD("status",status);FIELD("diagnostic",diagnostic);FIELD("input_digest",inputDigest);FIELD("output_digest",outputDigest);FIELD("mesh_digest",meshDigest);FIELD("input_cuboids",inputCuboids);FIELD("boolean_operations",booleanOperations);FIELD("peak_complexity",peakComplexity);FIELD("simple",simple);FIELD("roundtrip_equal",roundtripEqual);d["budget"]=budget(f.budget);
    py::dict c;c["vertices"]=f.finalComplexity.vertices;c["halfedges"]=f.finalComplexity.halfedges;c["halffacets"]=f.finalComplexity.halffacets;c["volumes"]=f.finalComplexity.volumes;d["final_complexity"]=c;
    d["settings"]=boundarySettings(f.settings);return d;
}
py::dict geometry(const SourceMaterialBoundaryFacts& f) {
    py::dict d;FIELD("algorithm",algorithm);FIELD("status",status);FIELD("diagnostic",diagnostic);FIELD("input_digest",inputDigest);FIELD("product_digest",productDigest);FIELD("source_hash",sourceHash);FIELD("partition_hash",partitionHash);FIELD("binary_hash",binaryHash);FIELD("outer_mesh_hash",outerMeshHash);FIELD("material_mesh_hash",materialMeshHash);
    const bool attempted=f.status!="NOT_ATTEMPTED_UNKNOWN";d["state"]=attempted?"ATTEMPTED":"UNAVAILABLE";
    d["density"]=attempted?number(f.density):py::object(py::none());
    d["boundary_density"]=attempted?number(f.boundaryDensity):py::object(py::none());
    d["interior_density"]=attempted?number(f.interiorDensity):py::object(py::none());
    d["material_components_exported"]=f.materialComponentsExported;
    d["minimum_outer_wall"]=attempted?number(f.minimumOuterWall):py::object(py::none());
    d["union"]=unionFacts(f.unionFacts);d["resource_budget"]=budget(f.resourceBudget);d["outer_region"]=region(f.outerRegion);d["design_region"]=region(f.designRegion);d["material_region"]=region(f.materialRegion);d["carved_region"]=region(f.carvedRegion);
    d["boundary_material_region"]=region(f.boundaryMaterialRegion);d["interior_material_region"]=region(f.interiorMaterialRegion);
    d["material_boundary"]=boundary(f.materialBoundary);d["partition"]=partition(f.partition);d["dual_material_partition"]=partition(f.dualMaterialPartition);
    d["carved_subset_design"]=relation(f.carvedSubsetDesign);d["clearance"]=clearance(f.clearance);
    const bool integrated=std::isfinite(f.normalizedMomentError);
    d["moments_state"]=integrated?"AVAILABLE":"UNAVAILABLE";
    d["source_subtraction_moments"]=integrated?py::object(numbers(f.sourceSubtractionMoments.values())):py::object(py::none());
    d["integrated_moments"]=integrated?py::object(numbers(f.integratedMoments.values())):py::object(py::none());
    d["boundary_material_moments"]=integrated?py::object(numbers(f.boundaryMaterialMoments.values())):py::object(py::none());
    d["interior_material_moments"]=integrated?py::object(numbers(f.interiorMaterialMoments.values())):py::object(py::none());
    d["normalized_moment_residuals"]=numbers(f.normalizedMomentResiduals);d["normalized_moment_error"]=number(f.normalizedMomentError);d["minimum_central_inertia_eigenvalue"]=number(f.minimumCentralInertiaEigenvalue);return d;
}
#undef FIELD

struct OwnedRun {
    SourceAdaptiveSpinMaterialResult result;
    std::shared_ptr<const SourceSurfaceOctree> source;
    std::string settings;
    bool diskIterations=false;
    std::uint64_t iterationMaximumBytes=1073741824;
    std::string checkpoint;
};
py::dict diagnostics(const OwnedRun& r,bool includePasses) {
    py::dict d;d["schema"]="source-spin2014-native-material-v1";d["output_digest"]=r.result.outputDigest();d["solve"]=solveFacts(r.result.solve(),includePasses);d["geometry"]=geometry(r.result.boundary().facts());
    d["settings"]=py::module_::import("json").attr("loads")(r.settings);
    d["resource_budget_scope"]="NATIVE_SOURCE_MATERIAL_CHAIN_ONLY_EXCLUDES_INDEPENDENT_PYTHON_SEAL";
    const auto& s=*r.source;py::dict source;source["source_hash"]=s.sourceSha256();source["config_hash"]=s.configSha256();source["partition_hash"]=s.partitionSha256();source["outer_surface_moments"]=numbers(s.outerSurfaceMoments().values());source["design_surface_moments"]=numbers(s.designSurfaceMoments().values());source["carveable_moments"]=numbers(s.carveableMoments().values());source["unresolved_boundary_volume"]=number(s.unresolvedBoundaryVolume());
    const auto& f=s.stats();py::dict stats;stats["internal_nodes"]=f.internalNodes;stats["boundary_nodes"]=f.boundaryNodes;stats["full_nodes"]=f.fullNodes;stats["empty_nodes"]=f.emptyNodes;stats["maximum_visited_level"]=f.maximumVisitedLevel;stats["levels"]=f.levels;source["stats"]=stats;d["source"]=source;
    d["nonfinite_encoding"]="NULL_WITH_ENCLOSING_AVAILABILITY_OR_STATUS";return d;
}
// Export memory is one metadata/iteration record, independent of trace length.
// Disk mode opens one sequential reader per pass, retaining the native owner.
struct DiagnosticsIterator {
    std::shared_ptr<const OwnedRun> owner;
    bool encoded=false;
    bool compactStageAudit=false;
    std::size_t passIndex=0,traceIndex=0;
    bool headerPending=true,passPending=true;
    std::unique_ptr<SlqpIterationReader> reader;
    py::object encodeMetadata(const py::dict& event) const {
        if(!encoded)return event;
        const auto data=py::cast<std::string>(py::module_::import("json").attr("dumps")(
            event,py::arg("sort_keys")=true,py::arg("separators")=py::make_tuple(",",":"),py::arg("allow_nan")=false));
        return py::bytes(data+"\n");
    }
    py::object next() {
        py::dict event;
        const auto& passes=owner->result.solve().passes;
        if(headerPending){
            event["kind"]="header";event["payload"]=diagnostics(*owner,false);
            py::list counts;for(const auto& p:passes)counts.append(slqpTraceCount(p.solve));
            event["trace_counts"]=counts;headerPending=false;return encodeMetadata(event);
        }
        while(passIndex<passes.size()){
            const auto& p=passes[passIndex];
            if(passPending){
                event["kind"]="pass";event["pass_index"]=passIndex;event["payload"]=passFacts(p,false);
                passPending=false;return encodeMetadata(event);
            }
            if(traceIndex<slqpTraceCount(p.solve)){
                event["kind"]="iteration";event["pass_index"]=passIndex;event["trace_index"]=traceIndex;
                SlqpIteration record;
                const SlqpIteration* current;
                if(p.solve.iterationStore){
                    if(!reader)reader=std::make_unique<SlqpIterationReader>(p.solve.iterationStore->path());
                    if(!reader->next(record))throw SlqpStorageError("SLQP iteration storage: early end");
                    current=&record;
                }else current=&p.solve.trace[traceIndex];
                if(compactStageAudit) {
                    if(current->iteration!=static_cast<int>(traceIndex))
                        throw SlqpStorageError("SLQP iteration storage: iteration sequence mismatch");
                    const bool terminal=traceIndex+1==slqpTraceCount(p.solve);
                    if(current->finalState!=terminal)
                        throw SlqpStorageError("SLQP iteration storage: terminal marker mismatch");
                    auto bytes="{\"kind\":\"iteration\",\"pass_index\":"+std::to_string(passIndex)+
                        ",\"trace_index\":"+std::to_string(traceIndex)+",\"payload\":"+
                        (terminal?traceJson(*current):std::string("{}"))+"}\n";
                    ++traceIndex;return py::bytes(bytes);
                }
                if(encoded) {
                    auto bytes="{\"kind\":\"iteration\",\"pass_index\":"+std::to_string(passIndex)+
                        ",\"trace_index\":"+std::to_string(traceIndex)+",\"payload\":"+traceJson(*current)+"}\n";
                    ++traceIndex;return py::bytes(bytes);
                }
                event["payload"]=trace(*current);
                ++traceIndex;return event;
            }
            if(reader){SlqpIteration extra;if(reader->next(extra))throw SlqpStorageError("SLQP iteration storage: excess records");reader.reset();}
            ++passIndex;traceIndex=0;passPending=true;
        }
        throw py::stop_iteration();
    }
};
py::dict traceStorage(const OwnedRun& r){
    py::dict d;std::size_t count=0,retained=0,peak=0,leafCount=0,betaScalars=0;
    std::uint64_t bytes=0;py::list files;
    for(const auto& p:r.result.solve().passes){
        count+=slqpTraceCount(p.solve);
        leafCount+=p.before.size()+p.after.size();betaScalars+=p.inputBeta.size()+p.outputBeta.size();
        if(p.solve.iterationStore){const auto& s=*p.solve.iterationStore;retained+=s.retainedRecords();peak=std::max(peak,s.peakRetainedRecords());bytes+=s.bytesWritten();files.append(s.path());}
        else retained+=p.solve.trace.size();
    }
    d["mode"]=r.diskIterations?"DISK_ITERATIONS_V1":"MEMORY";
    d["record_count"]=count;d["retained_iteration_records"]=retained;
    d["peak_retained_iteration_records"]=r.diskIterations?peak:count;
    d["bytes_written"]=bytes;d["files"]=files;
    d["maximum_bytes"]=r.iterationMaximumBytes;d["byte_budget_scope"]="CUMULATIVE_SOURCE_RUN_ITERATION_FILES";
    d["resident_pass_count"]=r.result.solve().passes.size();
    d["resident_pass_partition_leaves"]=leafCount;d["resident_pass_beta_scalars"]=betaScalars;
    d["resident_pass_partition_and_beta_payload_bytes"]=leafCount*sizeof(SourceAdaptiveLeaf)+betaScalars*sizeof(double);
    d["resident_pass_memory_scope"]="PARTITION_AND_BETA_PAYLOAD_ONLY_EXCLUDES_OTHER_PASS_FIELDS_CAPACITY_AND_ALLOCATOR_OVERHEAD";
    d["iteration_memory_scope"]="PENDING_HISTORY_ONLY_EXCLUDES_ACTIVE_SOLVER_RECORD_READER_WRITER_BUFFERS_AND_EXPLICIT_DIAGNOSTICS";
    d["scientific_hash_io"]=r.diskIterations?"ONE_FULL_SEQUENTIAL_ITERATION_READBACK_PER_SOLVE":"MEMORY_TRAVERSAL";
    d["maximum_record_bytes"]=67108864;
    return d;
}
py::object mesh(const OwnedRun& r,int which) {
    const auto& p=r.result.boundary().product();if(!p)return py::none();
    // Explicit evaluated copies: returned writable arrays cannot mutate owner.
    if(which==0)return py::make_tuple(Eigen::MatrixXd(p->outerVertices()),Eigen::MatrixXi(p->outerFaces()));
    if(which==1)return py::make_tuple(Eigen::MatrixXd(p->materialVertices()),Eigen::MatrixXi(p->materialFaces()));
    if(which==2)return py::make_tuple(Eigen::MatrixXd(p->cavityVertices()),Eigen::MatrixXi(p->cavityFaces()));
    if(!p->facts().materialComponentsExported)return py::none();
    if(which==3)return py::make_tuple(Eigen::MatrixXd(p->boundaryMaterialVertices()),Eigen::MatrixXi(p->boundaryMaterialFaces()));
    return py::make_tuple(Eigen::MatrixXd(p->interiorMaterialVertices()),Eigen::MatrixXi(p->interiorMaterialFaces()));
}
} // namespace

void registerSourceSpinMaterialBindings(py::module_& m,Eigen::MatrixXd (*vertices)(const py::handle&),Eigen::MatrixXi (*faces)(const py::handle&)) {
    py::class_<DiagnosticsIterator>(m,"NativeSourceSpinDiagnosticsIterator")
        .def("__iter__",[](DiagnosticsIterator& it)->DiagnosticsIterator&{return it;},py::return_value_policy::reference_internal)
        .def("__next__",&DiagnosticsIterator::next);
    py::class_<OwnedRun,std::shared_ptr<OwnedRun>>(m,"NativeSourceSpinMaterialResult")
        .def_property_readonly("checkpoint",[](const OwnedRun& r)->py::object {
            return r.checkpoint.empty()?py::object(py::none()):py::object(py::bytes(r.checkpoint));
        })
        .def_property_readonly("trace_storage",&traceStorage)
        .def_property_readonly("diagnostics",[](const OwnedRun& r){return diagnostics(r,true);})
        .def("iter_diagnostics",[](std::shared_ptr<OwnedRun> r){return DiagnosticsIterator{std::move(r)};},
            "Lazy complete header/pass/iteration payloads; retains native owner; does not grant a seal.")
        .def("iter_diagnostics_json",[](std::shared_ptr<OwnedRun> r){return DiagnosticsIterator{std::move(r),true};},
            "Lossless native JSON-line bytes; use write_encoded_archive for sequence and integrity checks.")
        .def("iter_stage_audit_json",[](std::shared_ptr<OwnedRun> r){return DiagnosticsIterator{std::move(r),true,true};},
            "Compact audit events: validates every native record, emits only terminal iteration payloads.")
        .def_property_readonly("has_geometry",[](const OwnedRun& r){return bool(r.result.boundary().product());})
        .def_property_readonly("outer",[](const OwnedRun& r){return mesh(r,0);})
        .def_property_readonly("material",[](const OwnedRun& r){return mesh(r,1);})
        .def_property_readonly("cavity",[](const OwnedRun& r){return mesh(r,2);})
        .def_property_readonly("boundary_material",[](const OwnedRun& r){return mesh(r,3);})
        .def_property_readonly("interior_material",[](const OwnedRun& r){return mesh(r,4);});
    m.def("inspect_source_spin2014_domain",[vertices,faces](const py::object& ov,const py::object& of,
            const py::object& dv,const py::object& df,const py::dict& sourceSettings,
            const py::object& initialLevel,const py::object& maximumOptimizationLevel,
            const py::object& maximumLeafCount,bool includeDetails) {
        const auto config=sourceConfig(sourceSettings);py::dict levels;
        levels["initial_level"]=initialLevel;levels["maximum_optimization_level"]=maximumOptimizationLevel;
        levels["maximum_leaf_count"]=maximumLeafCount;SettingsReader reader(levels);
        const int initial=int(reader.count("initial_level",0,config.maximumDepth));
        const int maximum=int(reader.count("maximum_optimization_level",initial,config.maximumDepth));
        const auto budget=reader.count("maximum_leaf_count",1,std::size_t(INT_MAX/6));reader.finish();
        const auto outerV=vertices(ov);const auto outerF=faces(of);const auto designV=vertices(dv);const auto designF=faces(df);
        std::shared_ptr<const SourceSurfaceOctree> source;
        {py::gil_scoped_release release;source=std::make_shared<const SourceSurfaceOctree>(outerV,outerF,designV,designF,config);}
        SourceAdaptiveDomain domain(source,initial,maximum,budget);domain.validateCoverage();py::dict out;
        out["schema"]="source-spin2014-domain-inspection-v1";out["source"]=sourceFacts(*source);
        out["initial_level"]=initial;out["maximum_optimization_level"]=maximum;out["maximum_leaf_count"]=budget;
        out["inspection_detail"]=includeDetails?"FULL":"SUMMARY";
        out["partition_hash"]=domain.partitionHash();out["source_cell_count"]=source->cells().size();
        out["leaf_count"]=domain.leaves().size();out["edge_count"]=domain.edges().size();
        py::dict leafLevels;std::map<int,std::size_t> levelCounts;
        for(const auto& leaf:domain.leaves())++levelCounts[leaf.level];
        for(const auto& item:levelCounts)leafLevels[py::str(std::to_string(item.first))]=py::int_(item.second);
        out["leaf_level_counts"]=leafLevels;
        if(includeDetails) {
            out["source_cells"]=sourceCells(*source);out["leaves"]=leaves(domain.leaves());
            py::list edges;for(const auto& edge:domain.edges())edges.append(py::make_tuple(edge.first,edge.second));out["edges"]=edges;
            out["laplacian_triplets"]=sparseTriplets(domain.uniformLaplacian());
        }
        return out;
    },py::arg("outer_vertices"),py::arg("outer_faces"),py::arg("design_vertices"),py::arg("design_faces"),py::kw_only(),
    py::arg("source_settings"),py::arg("initial_level"),py::arg("maximum_optimization_level"),py::arg("maximum_leaf_count"),
    py::arg("include_details")=true,
    "Read-only source-octree/adaptive-domain inspection; optionally summarizes without expanding cells, leaves, edges, or sparse Laplacian, never solves or mutates a partition.");
    m.def("solve_source_spin2014_material",[vertices,faces](const py::object& ov,const py::object& of,const py::object& dv,const py::object& df,const py::dict& sourceSettings,const py::dict& solverSettings,const py::dict& boundarySettings,const py::object& density,const py::object& wall,const py::object& storeDirectory,std::uint64_t storeMaximumBytes,const py::object& maximumNewPasses,const py::object& resumeCheckpoint){
        const auto config=sourceConfig(sourceSettings);auto request=solverConfig(solverSettings,config.maximumDepth);
        if(!storeDirectory.is_none()) {
            require(py::isinstance<py::str>(storeDirectory),"iteration_store_directory must be a string or None");
            request.iterationStoreDirectory=py::cast<std::string>(storeDirectory);
            require(!request.iterationStoreDirectory.empty(),"iteration_store_directory must not be empty");
            require(request.iterationStoreDirectory.find('\0')==std::string::npos,"iteration_store_directory contains NUL");
        }
        request.iterationStoreMaximumBytes=storeMaximumBytes;
        const bool staged=!maximumNewPasses.is_none();
        SourceAdaptiveStageControl control;std::string resume;
        if(staged) {
            py::dict value;value["maximum_new_passes"]=maximumNewPasses;
            SettingsReader reader(value);control.maximumNewPasses=int(reader.count("maximum_new_passes",1,INT_MAX));
            require(!request.iterationStoreDirectory.empty(),"staged material requires disk iteration storage");
        }
        if(!resumeCheckpoint.is_none()) {
            require(staged,"resume requires stage control");
            require(py::isinstance<py::bytes>(resumeCheckpoint),"checkpoint requires bytes");
            require(py::len(resumeCheckpoint)>0 && py::len(resumeCheckpoint)<=64ULL*1024*1024,"checkpoint metadata size");
            resume=py::cast<std::string>(resumeCheckpoint);
        }
        request.density=scalar(density);require(request.density>0,"density must be positive");
        const auto bounds=boundaryConfig(boundarySettings,scalar(wall));
        for(double density:{request.resolvedBoundaryDensity(),request.resolvedInteriorDensity()})
            for(int power=3;power<=5;++power){double scale=density*std::pow(config.rootSide,power);require(std::isfinite(scale)&&scale>0,"invalid source scaling");}
        const auto outerV=vertices(ov);const auto outerF=faces(of);const auto designV=vertices(dv);const auto designF=faces(df);
        auto copiedSettings=settingsJson(config,request,bounds);
        std::shared_ptr<const SourceSurfaceOctree> source;
        // Native work has no Python callbacks and consumes owned copied inputs.
        py::gil_scoped_release release;
        {SPIN2014_DIAG_TIMER("binding_source_octree_build",-1);
            source=std::make_shared<const SourceSurfaceOctree>(outerV,outerF,designV,designF,config);}
        request.source=source;
        // Material tokens bind geometry export policy too, unlike the lower-level
        // solver-only codec. Verify before creating any new iteration directory.
        const auto materialSettingsHash=hashSourceMaterialBoundarySettings(bounds);
        if(!resume.empty()) {
            if(resume.size()<=136 || resume.substr(0,8)!="SCMAT001" ||
               resume.substr(8,64)!=materialSettingsHash ||
               manufacturing::computeContentSha256(resume.substr(0,resume.size()-64))!=resume.substr(resume.size()-64))
                throw std::invalid_argument("material checkpoint settings/version mismatch");
            {SPIN2014_DIAG_TIMER("binding_checkpoint_decode",-1);
                control.resume=decodeSourceAdaptiveCheckpoint(request,resume.substr(72,resume.size()-136));}
        }
        std::shared_ptr<const SourceAdaptiveCheckpoint> checkpoint;
        auto result=[&]{SPIN2014_DIAG_TIMER("binding_material_stages",-1);
            return solveSourceAdaptiveSpinMaterialStages(request,bounds,staged?&control:nullptr,checkpoint);}();
        std::string encoded;
        if(checkpoint) {
            {SPIN2014_DIAG_TIMER("binding_checkpoint_encode",-1);
                encoded="SCMAT001"+materialSettingsHash+encodeSourceAdaptiveCheckpoint(*checkpoint);}
            // Integrity only, not authentication against a writer able to
            // recompute hashes. The supervisor must bind immutable identities.
            encoded+=manufacturing::computeContentSha256(encoded);
            if(encoded.size()>64ULL*1024*1024)throw std::length_error("material checkpoint metadata size");
        }
        return OwnedRun{std::move(result),std::move(source),std::move(copiedSettings),!request.iterationStoreDirectory.empty(),storeMaximumBytes,std::move(encoded)};
    },py::arg("outer_vertices"),py::arg("outer_faces"),py::arg("design_vertices"),py::arg("design_faces"),py::kw_only(),py::arg("source_settings"),py::arg("solver_settings"),py::arg("boundary_settings"),py::arg("density"),py::arg("minimum_outer_wall"),py::arg("iteration_store_directory")=py::none(),py::arg("iteration_store_max_bytes")=1073741824ULL,py::arg("maximum_new_passes")=py::none(),py::arg("resume_checkpoint")=py::none(),
    "Actual source solver/private-terminal material wrapper, called once; no caller status or moments accepted.");
    m.def("replay_source_spin2014_postsolve",[vertices,faces](const py::object& ov,const py::object& of,const py::object& dv,const py::object& df,const py::dict& sourceSettings,const py::dict& solverSettings,const py::dict& boundarySettings,const py::object& density,const py::object& wall,const py::bytes& checkpoint,const std::string& tracePath,std::uint64_t storeMaximumBytes){
        const auto config=sourceConfig(sourceSettings);auto request=solverConfig(solverSettings,config.maximumDepth);
        request.density=scalar(density);require(request.density>0,"density must be positive");
        const auto bounds=boundaryConfig(boundarySettings,scalar(wall));
        for(double density:{request.resolvedBoundaryDensity(),request.resolvedInteriorDensity()})
            for(int power=3;power<=5;++power){const double scale=density*std::pow(config.rootSide,power);require(std::isfinite(scale)&&scale>0,"invalid source scaling");}
        const auto outerV=vertices(ov);const auto outerF=faces(of);
        const auto designV=vertices(dv);const auto designF=faces(df);
        require(!tracePath.empty() && tracePath.find('\0')==std::string::npos,"immutable trace path required");
        request.iterationStoreDirectory=std::filesystem::absolute(std::filesystem::path(tracePath)).parent_path().string();
        request.iterationStoreMaximumBytes=storeMaximumBytes;
        const auto materialHash=hashSourceMaterialBoundarySettings(bounds);const auto wire=std::string(checkpoint);
        SourceAdaptivePostSolveReplay replay;
        {py::gil_scoped_release release;
            auto source=std::make_shared<const SourceSurfaceOctree>(outerV,outerF,designV,designF,config);
            request.source=std::move(source);
            replay=replaySourceAdaptivePostSolve(request,wire,materialHash,tracePath,spin2014G7PostSolvePrePassCount);
        }
        py::dict out;out["scope"]=replay.scope;out["optimizer_invoked"]=replay.optimizerInvoked;
        out["binary_validation"]=replay.binaryValidation;out["scientific_acceptance"]=replay.scientificAcceptance;
        out["checkpoint_sha256"]=replay.checkpointSha256;out["trace_record_stream_sha256"]=replay.traceRecordStreamSha256;
        out["encoded_checkpoint_sha256"]=replay.encodedCheckpointSha256;
        out["terminal_state_dimension"]=replay.terminalState.size();out["terminal_state_sha256"]=vectorDigest(replay.terminalState);
        out["terminal_kkt"]=compactKkt(replay.terminalKkt);
        out["timings_ns"]=replay.timingsNanoseconds;out["summary"]=postSolveSummary(replay.result);
        out["diagnostic_bytes"]=py::bytes(replay.diagnosticBytes);return out;
    },py::arg("outer_vertices"),py::arg("outer_faces"),py::arg("design_vertices"),py::arg("design_faces"),py::kw_only(),py::arg("source_settings"),py::arg("solver_settings"),py::arg("boundary_settings"),py::arg("density"),py::arg("minimum_outer_wall"),py::arg("checkpoint"),py::arg("trace_path"),py::arg("iteration_store_max_bytes")=1073741824ULL,
    "Read-only post-SLQP diagnostic replay. It validates immutable evidence and never resumes, solves, rounds, exports, or grants scientific acceptance.");
}
