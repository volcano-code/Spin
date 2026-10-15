#include "Spin2014AuthorityBindings.h"
#include "spin2014/QuadricAuthorityNewtonState.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <gmpxx.h>
namespace py=pybind11;
namespace {
py::list modelVertexList(const Eigen::MatrixXd& vertices) {
    py::list result;
    for(Eigen::Index i=0;i<vertices.rows();++i) {
        py::list row;
        for(Eigen::Index j=0;j<vertices.cols();++j) {
            const mpq_class exact(vertices(i,j));
            row.append(exact.get_num().get_str()+"/"+exact.get_den().get_str());
        }
        result.append(std::move(row));
    }
    return result;
}
py::list faceList(const Eigen::MatrixXi& faces) {
    py::list result;
    for(Eigen::Index row=0;row<faces.rows();++row) {
        py::list triangle;
        for(Eigen::Index col=0;col<faces.cols();++col)triangle.append(faces(row,col));
        result.append(std::move(triangle));
    }
    return result;
}
py::dict materialToPython(const spin::QuadricAuthorityTotalMomentEvaluation& owned,
    const spin::QuadricAuthorityMaterialContext& context,
    const spin::QuadricMomentEvaluatorSettings& evaluation) {
    const auto& material=owned.material;
    if(owned.contextIdentity!=context.identity() || !owned.derivative.defined)
        throw std::invalid_argument("authority material serialization rejected unbound context or derivative");
    py::dict result;
    result["context_identity"]=context.identity();
    result["request_identity"]=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation);
    // These are the original represented INPUT coordinates, not coefficient
    // or material display arrays. Preserve their exact dyadic values.
    result["outer_vertices_model_authority"]=modelVertexList(context.request().outerSurface.vertices);
    result["outer_faces"]=faceList(context.request().outerSurface.faces);
    result["design_vertices_model_authority"]=modelVertexList(context.request().designDomain.vertices);
    result["design_faces"]=faceList(context.request().designDomain.faces);
    result["model_to_target_authority"]=context.modelToTarget();
    result["target_to_model_authority"]=context.targetToModel();
    result["fixed_material_density"]=context.fixedDensity();
    result["design_material_density"]=context.request().materialDensity;
    result["fixed_moments_target"]=context.fixedMomentsTarget();
    result["variable_moments_target"]=material.variableMomentsExactTarget;
    result["total_moments_target"]=owned.totalMomentsTarget;
    result["material_vertices_target_authority"]=material.materialExactVerticesTarget;
    result["material_faces"]=faceList(material.materialBoundaryTarget.faces);
    py::list births;
    for(Eigen::Index i=0;i<material.materialBirthFaces.size();++i)births.append(material.materialBirthFaces(i));
    result["material_birth_faces"]=std::move(births);
    result["design_vertices_target_authority"]=material.designExactVerticesTarget;
    result["source_vertices_target_authority"]=material.quadricTriangulation.exactVertices;
    result["source_coordinate_lower"]=material.quadricTriangulation.coordinateLower;
    result["source_coordinate_upper"]=material.quadricTriangulation.coordinateUpper;
    result["source_faces"]=faceList(material.quadricTriangulation.closedBoundary.faces);
    result["source_construction_identity"]=material.quadricTriangulation.rationalConstructionIdentity;
    result["provenance_topology_hash"]=material.provenanceTopologyHash;
    result["material_topology_hash"]=material.topologyHash;
    result["material_vertex_provenance"]=material.vertexProvenance;
    result["surface_type"]=static_cast<int>(material.quadricTriangulation.type);
    result["effective_normal_subdivisions"]=material.quadricTriangulation.effectiveIcosphereSubdivisions;
    result["source_truncation_radius"]=material.quadricTriangulation.truncationRadius;
    result["source_interface_face_count"]=material.quadricTriangulation.interfaceFaceCount;
    result["source_cap_face_count"]=material.quadricTriangulation.capFaceCount;
    result["native_moment_jacobian"]=owned.derivative.differential.coefficientJacobian;
    result["material_scientific_hash"]=owned.scientificHash;
    result["boolean_backend_identity"]=material.booleanBackendIdentity;
    return result;
}
py::dict pointToPython(const spin2014::QuadricAuthorityNewtonPoint& point,
    const spin::QuadricAuthorityMaterialContext& context,
    const spin::QuadricMomentEvaluatorSettings& evaluation) {
    if(point.evaluatorRequestIdentity!=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation))
        throw std::invalid_argument("authority point serialization rejected a different evaluator recipe");
    auto result=materialToPython(point.material,context,evaluation);
    result["schema"]="spin2014-authority-actual-point-v1";
    result["scope"]="NUMERICAL_POINT_NOT_FINAL_SEAL";
    result["coefficient_policy"]=point.state.coefficientPolicy;
    result["precision_bits"]=point.state.precisionBits;
    result["geometry_authority"]=point.state.geometryAuthority;
    result["dual_authority"]=point.state.dualAuthority;
    result["objective"]=point.state.objective.objective;
    result["objective_gradient"]=point.state.objective.objectiveGradient;
    result["objective_hessian"]=point.state.objective.objectiveHessian;
    result["constraint_jacobian"]=point.state.objective.constraintJacobian;
    result["constraint_hessians"]=point.state.objective.constraintHessians;
    result["objective_identity"]=point.state.objective.objectiveIdentity;
    result["constraints"]=point.state.objective.constraints;
    result["raw_equation_residual"]=point.state.rawEquationResidual;
    result["multiplier_mu"]=point.state.multiplierMu;
    result["multipliers"]=point.state.multipliers;
    result["matrix_authority"]=point.state.matrix;
    result["rhs_authority"]=point.state.rightHandSide;
    result["evaluation_map_identity"]=point.state.evaluationMapIdentity;
    result["material_scientific_hash"]=point.material.scientificHash;
    result["state_scientific_hash"]=point.state.scientificHash;
    return result;
}
}
namespace {
py::dict authorityEventToPython(const spin2014::QuadricAuthorityNewtonEvent& event) {
    py::dict result;
    result["schema"]="spin2014-authority-live-event-v1";
    result["scope"]="OBSERVATION_NOT_CHECKPOINT_OR_FINAL_SEAL";
    result["sequence"]=event.sequence;result["restart"]=event.restart;
    result["evaluation_index"]=event.evaluationIndex;result["type"]=event.type;
    result["status"]=event.status;result["diagnostic"]=event.diagnostic;
    result["geometry_authority"]=event.geometryAuthority;result["dual_authority"]=event.dualAuthority;
    result["state_identity"]=event.stateIdentity;result["material_identity"]=event.materialIdentity;
    result["trial_identity"]=event.trialIdentity;result["raw_residual"]=event.rawResidual;
    result["arc_length"]=event.arcLength;result["accepted"]=event.accepted;
    return result;
}
}
void registerSpin2014AuthorityBindings(py::module_& module) {
    module.def("solve_spin2014_authority_quadric_hollowing", [](
        spin::SpinDesignRequest request,spin2014::QuadricObjectiveSettings objective,
        double fixedDensity,spin::QuadricMomentEvaluatorSettings evaluation,
        unsigned budget,const std::string& linearLimit,py::object eventSink) {
        if(!eventSink.is_none() && !PyCallable_Check(eventSink.ptr()))
            throw py::type_error("event_sink must be callable or None");
        const auto context=[&]{py::gil_scoped_release release;
            return spin::makeQuadricAuthorityMaterialContext(request,fixedDensity);}();
        const auto contextualEvent=[&](const spin2014::QuadricAuthorityNewtonEvent& event) {
            auto value=authorityEventToPython(event);
            value["context_identity"]=context.identity();
            value["request_identity"]=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation);
            value["objective_identity"]=spin2014::spin2014QuadricObjectiveHash(objective);
            value["coefficient_policy"]=spin::QuadricCoefficientAuthority10::policy;
            value["precision_bits"]=spin::QuadricCoefficientAuthority10::precisionBits;
            return value;
        };
        const auto observe=[&](const spin2014::QuadricAuthorityNewtonEvent& event) {
            if(!eventSink.is_none()) {
                py::gil_scoped_acquire acquire;
                eventSink(contextualEvent(event));
            }
        };
        const auto runs=[&]{py::gil_scoped_release release;
            return spin2014::runSpin2014AuthorityNewtonRestarts(
                context,objective,evaluation,linearLimit,budget,observe);}();
        py::dict result;
        result["schema"]="spin2014-authority-restarts-v1";
        result["scope"]="NUMERICAL_CANDIDATE_NOT_FINAL_SEAL";
        result["function_evaluations"]=runs.productionEvaluations;
        result["seed"]=runs.seed;result["status"]=runs.status;
        result["protocol_scope"]=runs.protocolScope;
        result["scientific_hash"]=runs.scientificHash;
        result["context_identity"]=context.identity();
        result["request_identity"]=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation);
        py::list events;
        for(const auto& event:runs.events)events.append(contextualEvent(event));
        result["events"]=std::move(events);
        py::list attempts;
        for(const auto& attempt:runs.attempts) {
            py::dict record;
            record["function_evaluations"]=attempt.productionEvaluations;
            record["accepted_updates"]=attempt.acceptedUpdates;
            record["initial_geometry_authority"]=attempt.initialGeometryAuthority;
            record["initial_dual_authority"]=attempt.initialDualAuthority;
            record["initial_state_identity"]=attempt.initialStateIdentity;
            record["initial_raw_residual"]=attempt.initialRawResidual;
            record["initial_mu"]=attempt.initialMu;
            record["initial_total_moments"]=py::none();
            if(attempt.initialTotalMoments)record["initial_total_moments"]=*attempt.initialTotalMoments;
            record["initial_objective"]=attempt.initialObjective;
            record["initial_material_identity"]=attempt.initialMaterialIdentity;
            record["final_state_identity"]=attempt.finalStateIdentity;
            record["final_raw_residual"]=attempt.finalRawResidual;
            record["status"]=attempt.status;record["diagnostic"]=attempt.diagnostic;
            record["scientific_hash"]=attempt.scientificHash;
            record["step_identities"]=attempt.stepIdentities;
            record["accepted_state_identities"]=attempt.acceptedStateIdentities;
            py::list trials;
            for(const auto& trial:attempt.trials) {
                py::dict item;
                item["arc_length"]=trial.arcLength;item["reason"]=trial.reason;
                item["diagnostic"]=trial.diagnostic;
                item["raw_equation_residual"]=trial.rawEquationResidual;
                item["trial_identity"]=trial.trialIdentity;
                item["geometry_authority"]=trial.geometryAuthority;
                item["dual_authority"]=trial.dualAuthority;
                item["accepted"]=trial.accepted;trials.append(std::move(item));
            }
            record["trials"]=std::move(trials);attempts.append(std::move(record));
        }
        result["attempts"]=std::move(attempts);
        result["stationary_attempt_index"]=py::none();
        result["stationary_point"]=py::none();
        if(runs.stationaryAttemptIndex) {
            result["stationary_attempt_index"]=*runs.stationaryAttemptIndex;
            const auto& point=*runs.attempts.at(*runs.stationaryAttemptIndex).stationaryCandidate;
            result["stationary_point"]=pointToPython(point,context,evaluation);
        }
        return result;
    },py::arg("request"),py::arg("objective_settings"),py::arg("fixed_material_density"),
      py::arg("evaluator_settings"),py::arg("evaluation_budget_per_attempt")=60,
      py::arg("absolute_linear_residual_limit")="1/1000000000000",py::arg("event_sink")=py::none());
    module.def("normalize_spin2014_authority_coefficients", [](
        const spin::QuadricCoefficientAuthority10::RationalVector& geometry) {
        const auto state=spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(geometry);
        py::dict result;
        result["coefficient_policy"]=state.policy;result["precision_bits"]=state.precisionBits;
        result["geometry_authority"]=state.geometry();result["dual_authority"]=state.dual();
        return result;
    },py::arg("geometry_rational_input"));
    module.def("update_spin2014_authority_coefficients", [](
        const spin::QuadricCoefficientAuthority10::RationalVector& geometry,
        const spin::QuadricCoefficientAuthority10::RationalVector& dual,
        const spin::QuadricCoefficientAuthority10::RationalVector& direction,
        const std::string& arc) {
        const auto initial=spin::QuadricCoefficientAuthority10::fromSavedState(geometry,dual,
            spin::QuadricCoefficientAuthority10::policy,spin::QuadricCoefficientAuthority10::precisionBits);
        const auto updated=initial.geodesicUpdate(direction,arc);
        py::dict result;
        result["schema"]="spin2014-authority-coefficient-update-v1";
        result["scope"]="COEFFICIENT_UPDATE_NOT_MATERIAL_OR_CONVERGENCE";
        result["coefficient_policy"]=updated.policy;result["precision_bits"]=updated.precisionBits;
        result["geometry_authority"]=updated.geometry();result["dual_authority"]=updated.dual();
        return result;
    },py::arg("geometry_authority"),py::arg("dual_authority"),
      py::arg("dual_direction"),py::arg("arc_length"));
    module.def("evaluate_spin2014_authority_material_mapping", [](
        const spin::SpinDesignRequest& request,double fixedDensity,
        const spin::QuadricCoefficientAuthority10::RationalVector& geometry,
        const spin::QuadricCoefficientAuthority10::RationalVector& dual,
        const spin::QuadricMomentEvaluatorSettings& evaluation) {
        const auto authority=spin::QuadricCoefficientAuthority10::fromSavedState(geometry,dual,
            spin::QuadricCoefficientAuthority10::policy,spin::QuadricCoefficientAuthority10::precisionBits);
        const auto context=[&]{py::gil_scoped_release release;
            return spin::makeQuadricAuthorityMaterialContext(request,fixedDensity);}();
        const auto owned=[&]{py::gil_scoped_release release;
            return spin::evaluateQuadricMaterialMoments(context,authority,evaluation);}();
        auto result=materialToPython(owned,context,evaluation);
        result["schema"]="spin2014-authority-material-map-v1";
        result["scope"]="MOMENT_MAPPING_NOT_SOLVER_OR_FINAL_SEAL";
        result["coefficient_policy"]=authority.policy;
        result["precision_bits"]=authority.precisionBits;
        result["geometry_authority"]=authority.geometry();
        result["dual_authority"]=authority.dual();
        return result;
    },py::arg("request"),py::arg("fixed_material_density"),
      py::arg("geometry_authority"),py::arg("dual_authority"),py::arg("evaluator_settings"));
    module.def("evaluate_spin2014_authority_quadric_point", [](
        const spin::SpinDesignRequest& request,const spin2014::QuadricObjectiveSettings& objective,
        double fixedDensity,const spin::QuadricCoefficientAuthority10::RationalVector& geometry,
        const spin::QuadricCoefficientAuthority10::RationalVector& dual,
        const spin::QuadricMomentEvaluatorSettings& evaluation) {
        const auto authority=spin::QuadricCoefficientAuthority10::fromSavedState(geometry,dual,
            spin::QuadricCoefficientAuthority10::policy,spin::QuadricCoefficientAuthority10::precisionBits);
        const auto context=[&]{py::gil_scoped_release release;
            return spin::makeQuadricAuthorityMaterialContext(request,fixedDensity);}();
        const auto point=[&]{py::gil_scoped_release release;
            return spin2014::evaluateSpin2014AuthorityNewtonPoint(context,objective,authority,evaluation);}();
        return pointToPython(point,context,evaluation);
    },py::arg("request"),py::arg("objective_settings"),py::arg("fixed_material_density"),
      py::arg("geometry_authority"),py::arg("dual_authority"),py::arg("evaluator_settings"));
    module.def("evaluate_spin2014_authority_moment_problem", [](
        const spin2014::QuadricObjectiveSettings& settings,
        const spin2014::AuthorityMomentVector10& moments) {
        const auto value=spin2014::evaluateSpin2014AuthorityMomentProblem(settings,moments);
        py::dict result;
        result["schema"]="spin2014-authority-objective-evaluation-v1";
        result["scope"]="NUMERICAL_OBJECTIVE_ONLY";
        result["precision_bits"]=value.precisionBits;
        result["objective_profile"]=value.objectiveProfile;
        result["objective_identity"]=value.objectiveIdentity;
        result["raw_moments_authority"]=value.rawMomentsAuthority;
        result["objective"]=value.objective;
        result["objective_gradient"]=value.objectiveGradient;
        result["objective_hessian"]=value.objectiveHessian;
        result["constraints"]=value.constraints;
        result["constraint_jacobian"]=value.constraintJacobian;
        result["constraint_hessians"]=value.constraintHessians;
        return result;
    },py::arg("objective_settings"),py::arg("raw_moments_authority"));
    module.def("spin2014_authority_initial_sequence", [](
        spin::SpinMode mode,const spin::SpinSolverSettings& settings) {
        const auto starts=spin2014::makeSpin2014AuthorityInitialSequence(mode,settings);
        py::list result;
        for(const auto& start:starts) {
            py::dict state;
            state["coefficient_policy"]=start.policy;
            state["precision_bits"]=start.precisionBits;
            state["geometry_authority"]=start.geometry();state["dual_authority"]=start.dual();
            result.append(std::move(state));
        }
        return result;
    },py::arg("mode"),py::arg("solver_settings"));
}
