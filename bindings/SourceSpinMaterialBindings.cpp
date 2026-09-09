#include "SourceSpinMaterialBindings.h"
#include "ExactMeshSequenceIngestion.h"
#include "spin2014/SourceMaterialBoundary.h"
#include "spin2014/SlqpAuditTrace.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <cmath>
#include <climits>
#include <iomanip>
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
    s.audit.lpTimeLimit=d.real("lp_time_limit");s.audit.lpFeasibilityTolerance=d.real("lp_feasibility_tolerance");d.finish();
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
    r.initialLevel=int(d.count("initial_level",0,20));r.maximumOptimizationLevel=int(d.count("maximum_optimization_level",0,20));
    r.maximumAdaptivePasses=int(d.count("maximum_adaptive_passes",1,INT_MAX));
    r.maximumLeafCount=d.count("maximum_leaf_count",1,INT_MAX/6);r.maximumOutputCuboids=d.count("maximum_output_cuboids",1,INT_MAX);
    r.initialFill=d.real("initial_fill");r.initialPhi=d.real("initial_phi");r.fractionalEpsilon=d.real("fractional_epsilon");
    r.laplacianWeight=d.real("laplacian_weight");r.inertiaWeight=d.real("inertia_weight");r.centerWeight=d.real("center_weight");
    r.binaryConstraintTolerance=d.real("binary_constraint_tolerance");
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
    c["initial_level"]=r.initialLevel;c["maximum_optimization_level"]=r.maximumOptimizationLevel;
    c["maximum_adaptive_passes"]=r.maximumAdaptivePasses;c["maximum_leaf_count"]=r.maximumLeafCount;
    c["maximum_output_cuboids"]=r.maximumOutputCuboids;c["initial_fill"]=r.initialFill;c["initial_phi"]=r.initialPhi;
    c["fractional_epsilon"]=r.fractionalEpsilon;c["laplacian_weight"]=r.laplacianWeight;
    c["inertia_weight"]=r.inertiaWeight;c["center_weight"]=r.centerWeight;c["binary_constraint_tolerance"]=r.binaryConstraintTolerance;
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
py::dict trace(const SlqpIteration& t) {
    // Existing serializer is authoritative for the full per-iteration audit.
    std::ostringstream o;o.imbue(std::locale::classic());o<<std::setprecision(17)<<"{\"iteration\":"<<t.iteration;
    writeSlqpAuditJson(o,t);o<<'}';
    auto d=py::cast<py::dict>(py::module_::import("json").attr("loads")(o.str()));
    d["objective"]=number(t.objective);d["equality_infinity_norm"]=number(t.equalityInfinityNorm);
    d["projected_gradient_infinity_norm"]=number(t.projectedGradientInfinityNorm);d["kkt_infinity_norm"]=number(t.kktInfinityNorm);
    d["merit_value"]=number(t.meritValue);d["accepted_step"]=number(t.acceptedStep);d["phi"]=number(t.phi);
    d["line_search_evaluations"]=t.lineSearchEvaluations;d["active_variables"]=t.activeVariables;d["lbfgs_history_size"]=t.lbfgsHistorySize;return d;
}
py::dict solveFacts(const SourceAdaptiveSpinResult& r) {
    py::dict d;d["status"]=sourceAdaptiveSpinStatusToken(r.status);d["diagnostic"]=r.diagnostic;
    d["input_hash"]=r.inputHash;d["trace_hash"]=r.traceHash;d["output_hash"]=r.outputHash;d["final_partition_hash"]=r.finalPartitionHash;
    d["geometry_depth"]=r.geometryDepth;d["maximum_optimization_level"]=r.maximumOptimizationLevel;d["final_phi"]=number(r.finalPhi);
    d["relaxed_fill"]=numbers(r.relaxedFill);d["binary_fill"]=numbers(r.binaryFill);d["rounding_delta"]=numbers(r.roundingDelta);d["final_partition"]=leaves(r.finalPartition);
    py::list passes;
    for(const auto& p:r.passes) {
        py::dict a;a["index"]=p.index;a["before"]=leaves(p.before);a["after"]=leaves(p.after);
        a["before_hash"]=p.beforeHash;a["after_hash"]=p.afterHash;a["input_beta"]=numbers(p.inputBeta);a["output_beta"]=numbers(p.outputBeta);a["input_phi"]=number(p.inputPhi);
        a["split"]=transition(p.split);a["merge"]=transition(p.merge);a["level_counts"]=p.levelCounts;a["fractional_level_counts"]=p.fractionalLevelCounts;a["fractional_count"]=p.fractionalCount;
        const auto& s=p.solve;a["status"]=slqpStatusToken(s.status);a["diagnostic"]=s.diagnostic;a["trace_hash"]=s.traceHash;
        a["initial_phi"]=number(s.initialPhi);a["final_phi"]=number(s.finalPhi);a["final_state"]=numbers(s.finalState);a["beta"]=numbers(s.beta);
        a["bounded_variable_count"]=s.boundedVariableCount;a["final_kkt"]=kkt(s.finalKkt);a["algorithm"]=slqpAlgorithmToken(s.algorithm);a["formulation"]="PAPER_GIVENS";
        // Initial rejection leaves a default evaluation; mark it unavailable,
        // rather than publishing those default zeros as measured quantities.
        if(s.trace.empty())a["evaluation"]=py::none();
        else {py::dict e;e["objective"]=number(s.evaluation.objective);e["paper_objective"]=number(s.evaluation.paperObjective);e["regularization_energy"]=number(s.evaluation.regularizationEnergy);e["final_moments"]=numbers(s.evaluation.finalMoments.values());e["gradient"]=numbers(s.evaluation.gradient);e["constraints"]=numbers(s.evaluation.constraints);py::list jac;for(int i=0;i<s.evaluation.constraintJacobian.rows();++i)jac.append(numbers(s.evaluation.constraintJacobian.row(i)));e["constraint_jacobian"]=jac;a["evaluation"]=e;}
        a["evaluation_state"]=s.trace.empty()?"UNAVAILABLE":"AVAILABLE";
        py::list records;for(const auto& t:s.trace)records.append(trace(t));a["trace"]=records;passes.append(a);
    }d["passes"]=passes;
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
    d["density"]=attempted?number(f.density):py::object(py::none());d["minimum_outer_wall"]=attempted?number(f.minimumOuterWall):py::object(py::none());
    d["union"]=unionFacts(f.unionFacts);d["resource_budget"]=budget(f.resourceBudget);d["outer_region"]=region(f.outerRegion);d["design_region"]=region(f.designRegion);d["material_region"]=region(f.materialRegion);d["carved_region"]=region(f.carvedRegion);d["material_boundary"]=boundary(f.materialBoundary);d["partition"]=partition(f.partition);d["carved_subset_design"]=relation(f.carvedSubsetDesign);d["clearance"]=clearance(f.clearance);
    const bool integrated=std::isfinite(f.normalizedMomentError);
    d["moments_state"]=integrated?"AVAILABLE":"UNAVAILABLE";
    d["source_subtraction_moments"]=integrated?py::object(numbers(f.sourceSubtractionMoments.values())):py::object(py::none());
    d["integrated_moments"]=integrated?py::object(numbers(f.integratedMoments.values())):py::object(py::none());
    d["normalized_moment_residuals"]=numbers(f.normalizedMomentResiduals);d["normalized_moment_error"]=number(f.normalizedMomentError);d["minimum_central_inertia_eigenvalue"]=number(f.minimumCentralInertiaEigenvalue);return d;
}
#undef FIELD

struct OwnedRun {
    SourceAdaptiveSpinMaterialResult result;
    std::shared_ptr<const SourceSurfaceOctree> source;
    std::string settings;
};
py::dict diagnostics(const OwnedRun& r) {
    py::dict d;d["schema"]="source-spin2014-native-material-v1";d["output_digest"]=r.result.outputDigest();d["solve"]=solveFacts(r.result.solve());d["geometry"]=geometry(r.result.boundary().facts());
    d["settings"]=py::module_::import("json").attr("loads")(r.settings);
    d["resource_budget_scope"]="NATIVE_SOURCE_MATERIAL_CHAIN_ONLY_EXCLUDES_INDEPENDENT_PYTHON_SEAL";
    const auto& s=*r.source;py::dict source;source["source_hash"]=s.sourceSha256();source["config_hash"]=s.configSha256();source["partition_hash"]=s.partitionSha256();source["outer_surface_moments"]=numbers(s.outerSurfaceMoments().values());source["design_surface_moments"]=numbers(s.designSurfaceMoments().values());source["carveable_moments"]=numbers(s.carveableMoments().values());source["unresolved_boundary_volume"]=number(s.unresolvedBoundaryVolume());
    const auto& f=s.stats();py::dict stats;stats["internal_nodes"]=f.internalNodes;stats["boundary_nodes"]=f.boundaryNodes;stats["full_nodes"]=f.fullNodes;stats["empty_nodes"]=f.emptyNodes;stats["maximum_visited_level"]=f.maximumVisitedLevel;stats["levels"]=f.levels;source["stats"]=stats;d["source"]=source;
    d["nonfinite_encoding"]="NULL_WITH_ENCLOSING_AVAILABILITY_OR_STATUS";return d;
}
py::object mesh(const OwnedRun& r,int which) {
    const auto& p=r.result.boundary().product();if(!p)return py::none();
    // Explicit evaluated copies: returned writable arrays cannot mutate owner.
    if(which==0)return py::make_tuple(Eigen::MatrixXd(p->outerVertices()),Eigen::MatrixXi(p->outerFaces()));
    if(which==1)return py::make_tuple(Eigen::MatrixXd(p->materialVertices()),Eigen::MatrixXi(p->materialFaces()));
    return py::make_tuple(Eigen::MatrixXd(p->cavityVertices()),Eigen::MatrixXi(p->cavityFaces()));
}
} // namespace

void registerSourceSpinMaterialBindings(py::module_& m,Eigen::MatrixXd (*vertices)(const py::handle&),Eigen::MatrixXi (*faces)(const py::handle&)) {
    py::class_<OwnedRun>(m,"NativeSourceSpinMaterialResult")
        .def_property_readonly("diagnostics",&diagnostics)
        .def_property_readonly("has_geometry",[](const OwnedRun& r){return bool(r.result.boundary().product());})
        .def_property_readonly("outer",[](const OwnedRun& r){return mesh(r,0);})
        .def_property_readonly("material",[](const OwnedRun& r){return mesh(r,1);})
        .def_property_readonly("cavity",[](const OwnedRun& r){return mesh(r,2);});
    m.def("solve_source_spin2014_material",[vertices,faces](const py::object& ov,const py::object& of,const py::object& dv,const py::object& df,const py::dict& sourceSettings,const py::dict& solverSettings,const py::dict& boundarySettings,const py::object& density,const py::object& wall){
        const auto config=sourceConfig(sourceSettings);auto request=solverConfig(solverSettings,config.maximumDepth);
        request.density=scalar(density);require(request.density>0,"density must be positive");
        const auto bounds=boundaryConfig(boundarySettings,scalar(wall));
        for(int power=3;power<=5;++power){double scale=request.density*std::pow(config.rootSide,power);require(std::isfinite(scale)&&scale>0,"invalid source scaling");}
        const auto outerV=vertices(ov);const auto outerF=faces(of);const auto designV=vertices(dv);const auto designF=faces(df);
        auto copiedSettings=settingsJson(config,request,bounds);
        std::shared_ptr<const SourceSurfaceOctree> source;
        // Native work has no Python callbacks and consumes owned copied inputs.
        py::gil_scoped_release release;
        source=std::make_shared<const SourceSurfaceOctree>(outerV,outerF,designV,designF,config);request.source=source;
        return OwnedRun{solveSourceAdaptiveSpinMaterial(request,bounds),std::move(source),std::move(copiedSettings)};
    },py::arg("outer_vertices"),py::arg("outer_faces"),py::arg("design_vertices"),py::arg("design_faces"),py::kw_only(),py::arg("source_settings"),py::arg("solver_settings"),py::arg("boundary_settings"),py::arg("density"),py::arg("minimum_outer_wall"),
    "Actual source solver/private-terminal material wrapper, called once; no caller status or moments accepted.");
}
