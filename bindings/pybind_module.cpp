#include "api/BalanceModel.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "SourceSpinMaterialBindings.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
#include "Spin2014AuthorityBindings.h"
#endif
#endif
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactMaterialRegion.h"
#include "ExactMeshSequenceIngestion.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
#include "verification/ExactLibiglBoolean.h"
#endif
#endif
#if MAKE_IT_STAND_HAS_SPIN_IT_2014
#include "api/SpinIt2014.h"
#include "spin2014/AutomaticCageBinding.h"
#include "spin2014/CageSpinDeformation.h"
#include "spin2014/CageSpinDeformationSolver.h"
#include "spin2014/Spin2014GivensProblem.h"
#include "spin2014/SurfaceRawMomentDifferential.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
#include "spin2014/QuadricHollowing.h"
#include "spin2014/QuadricStaticBalance.h"
#endif
#endif
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
#include "spin/CanonicalFullRankHollowing.h"
#include "spin/BuoyantPlaneBalanceSolver.h"
#include "spin/HydrostaticGeometry.h"
#include "spin/PlanarFidgetCircleProblem.h"
#include "PlanarEvaluationLedger.h"
#include "spin/PlanarDiskIntersection.h"
#include "spin/PlaneMomentEvaluator.h"
#include "spin/PlaneStaticBalanceSolver.h"
#include "spin/QuadricHollowingSolver.h"
#include "spin/QuadricSourceMotion.h"
#endif

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>


namespace py = pybind11;


namespace
{

#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
py::dict rationalAuthorityToPython(
    const spin::ExactBooleanVertexCoordinates& material,
    const spin::ExactBooleanVertexCoordinates& design,
    const spin::ExactBooleanVertexCoordinates& source,
    const spin::ExactBooleanVertexCoordinates& lower,
    const spin::ExactBooleanVertexCoordinates& upper,
    const Eigen::MatrixXi& sourceFaces, const Eigen::VectorXi& births,
    const std::vector<std::string>& moments,
    const std::string& geometryHash, const std::string& sourceIdentity)
{
    if(material.empty() || design.empty() || source.empty() ||
       lower.size()!=source.size() || upper.size()!=source.size() ||
       moments.size()!=10 || geometryHash.empty() || sourceIdentity.empty())
        throw std::runtime_error("incomplete same-operation rational authority payload");
    py::dict out;
    out["schema"]="rational-quadric-production-authority-v1";
    out["coordinate_space"]="TARGET_FRAME";
    out["analytic_coordinates_exact"]=false;
    out["variable_moment_semantics"]="DENSITY_WEIGHTED";
    out["material_vertices"]=material;out["design_vertices"]=design;
    out["source_vertices"]=source;out["source_coordinate_lower"]=lower;
    out["source_coordinate_upper"]=upper;
    py::list faces, birthFaces;
    for(Eigen::Index f=0;f<sourceFaces.rows();++f) {
        py::list row;
        for(int k=0;k<3;++k)row.append(sourceFaces(f,k));
        faces.append(std::move(row));
    }
    for(Eigen::Index f=0;f<births.size();++f)birthFaces.append(births(f));
    out["source_faces"]=std::move(faces);out["material_birth_faces"]=std::move(birthFaces);
    out["variable_moments"]=moments;out["geometry_hash"]=geometryHash;
    out["source_construction_identity"]=sourceIdentity;
    return out;
}

py::dict rationalAuthorityToPython(const spin::FullRankMomentEvaluation& value)
{
    return rationalAuthorityToPython(value.materialExactVerticesTarget,
        value.designExactVerticesTarget,value.sourceExactVerticesTarget,
        value.sourceCoordinateLowerTarget,value.sourceCoordinateUpperTarget,
        value.sourceFaces,value.materialBirthFaces,value.variableMomentsExactTarget,
        value.geometryHash,value.sourceConstructionIdentity);
}

py::dict evaluationPayloadToPython(
    const spin::SpinSolverEvent::EvaluationPayload& payload
)
{
    Eigen::Matrix<double, 10, 1> coordinates;
    Eigen::Matrix<double, 10, 1> solverMoments;
    for (Eigen::Index i = 0; i < 10; ++i)
    {
        coordinates(i) = payload.solverCoordinates[static_cast<std::size_t>(i)];
        solverMoments(i) = payload.solverRawMoments[static_cast<std::size_t>(i)];
    }
    py::dict result;
    result["schema"] = "spin-evaluation-payload-v1";
    result["solver_coordinates"] = coordinates;
    result["solver_raw_moments"] = solverMoments;
    if (payload.variableRawMoments)
    {
        Eigen::Matrix<double, 10, 1> variableMoments;
        for (Eigen::Index i = 0; i < 10; ++i)
            variableMoments(i) = (*payload.variableRawMoments)[
                static_cast<std::size_t>(i)];
        result["variable_raw_moments"] = variableMoments;
    }
    else
    {
        result["variable_raw_moments"] = py::none();
    }
    result["solver_moment_jacobian"] = payload.solverMomentJacobian
        ? py::cast(*payload.solverMomentJacobian) : py::none();
    result["derivative_defined"] = payload.derivativeDefined;
    result["jacobian_coordinate_identity"] = payload.jacobianCoordinateIdentity;
    result["evaluation_map_identity"] = payload.evaluationMapIdentity;
    result["boolean_backend_identity"] = payload.booleanBackendIdentity;
    result["derivative_identity"] = payload.derivativeIdentity;
    result["derivative_diagnostic"] = payload.derivativeDiagnostic;
    result["geometry_hash"] = payload.geometryHash;
    result["topology_hash"] = payload.topologyHash;
    result["evaluator_wall_seconds"] = payload.evaluatorWallSeconds;
    result["residual_evaluation_seconds"] = payload.residualEvaluationSeconds;
    result["derivative_evaluation_seconds"] = payload.derivativeEvaluationSeconds;
    return result;
}


py::dict newtonSystemPayloadToPython(
    const spin::SpinSolverEvent::NewtonSystemPayload& payload
)
{
    py::dict result;
    result["schema"] = "spin-newton-system-payload-v2";
    result["matrix"] = payload.matrix;
    result["right_hand_side"] = payload.rightHandSide;
    result["computed_rank"] = payload.computedRank;
    result["equilibrated_rank"] = payload.equilibratedRank;
    result["multiplier_mu"] = payload.multiplierMu
        ? py::cast(*payload.multiplierMu) : py::none();
    result["constraint_multipliers"] = payload.constraintMultipliers;
    result["row_scales"] = payload.rowScales;
    result["column_scales"] = payload.columnScales;
    result["linear_residual_infinity_norm"] =
        payload.linearResidualInfinityNorm;
    result["linear_backward_error"] = payload.linearBackwardError;
    return result;
}


py::dict solverEventToPython(const spin::SpinSolverEvent& event)
{
    py::dict entry;
    entry["type"] = event.type;
    entry["restart"] = event.restart;
    entry["iteration"] = event.iteration;
    entry["evaluation_index"] = event.evaluationIndex;
    entry["restart_evaluation_index"] = event.restartEvaluationIndex;
    entry["step"] = event.step ? py::cast(*event.step) : py::none();
    entry["outcome"] = event.outcome;
    entry["reason"] = event.reason;
    entry["detail"] = event.detail;
    entry["quadric_hash"] = event.quadricHash;
    entry["normalized_equation_residual"] = event.normalizedEquationResidual
        ? py::cast(*event.normalizedEquationResidual) : py::none();
    entry["physical_relative_residual"] = event.physicalRelativeResidual
        ? py::cast(*event.physicalRelativeResidual) : py::none();
    entry["objective"] = event.objective
        ? py::cast(*event.objective) : py::none();
    entry["multiplier_mu"] = event.multiplierMu
        ? py::cast(*event.multiplierMu) : py::none();
    entry["material_geometry_hash"] = event.materialGeometryHash
        ? py::cast(*event.materialGeometryHash) : py::none();
    entry["material_topology_hash"] = event.materialTopologyHash
        ? py::cast(*event.materialTopologyHash) : py::none();
    entry["equation_profile"] = spin::spinEquationProfileToken(event.equationProfile);
    entry["solver_coordinate_identity"] =
        spin::spinSolverCoordinateIdentity(event.equationProfile);
    entry["raw_equation_residual"] = event.rawEquationResidual
        ? py::cast(*event.rawEquationResidual) : py::none();
    if (event.quadricCoefficients)
    {
        Eigen::Matrix<double, 10, 1> coefficients = Eigen::Map<
            const Eigen::Matrix<double, 10, 1>
        >(event.quadricCoefficients->data());
        entry["quadric_coefficients"] = coefficients;
    }
    else
    {
        entry["quadric_coefficients"] = py::none();
    }
    if (event.evaluationPayload)
        entry["evaluation_payload"] = evaluationPayloadToPython(*event.evaluationPayload);
    if (event.newtonSystemPayload)
        entry["newton_system_payload"] =
            newtonSystemPayloadToPython(*event.newtonSystemPayload);
    return entry;
}

py::dict fullRankResultToPython(
    const spin::FullRankHollowingResult& solved, bool captureEvaluationPayload)
{
    py::dict result;
    result["status"] = spin::spinSolverStatusToken(solved.numerical.status);
    result["diagnostic"] = solved.numerical.diagnostic;
    result["final_geometry_verified"] = solved.finalGeometryVerified;
    result["verification_diagnostic"] = solved.verificationDiagnostic;
    result["input_hash"] = solved.inputHash;
    result["output_hash"] = solved.outputHash;
    result["function_evaluations"] = solved.numerical.functionEvaluations;
    result["restart_count"] = solved.numerical.restartCount;
    result["normalized_equation_residual"] = solved.numerical.normalizedEquationResidual;
    result["physical_relative_residual"] = solved.numerical.physicalRelativeResidual;
    result["independent_raw_residual"] = solved.independentlyRecomputedRawResidual
        ? py::cast(*solved.independentlyRecomputedRawResidual) : py::none();
    result["multiplier_mu"] = solved.numerical.multiplierMu;
    result["multipliers"] = Eigen::VectorXd(solved.numerical.multipliers);
    result["quadric_coefficients"] = solved.numerical.quadric
        ? py::cast(Eigen::VectorXd(solved.numerical.quadric->values())) : py::none();
    result["equation_profile"] = spin::spinEquationProfileToken(
        solved.numerical.equationProfile
    );
    result["event_hash"] = solved.numerical.eventHash;
    result["trace_hash"] = solved.numerical.traceHash;
    if (captureEvaluationPayload)
        result["event_stream_version"] = solved.numerical.eventStreamVersion;
    result["restart_function_evaluations"] =
        solved.numerical.restartFunctionEvaluations;

    py::list trace;
    for (const auto& item : solved.numerical.trace)
    {
        py::dict entry;
        entry["restart"] = item.restart;
        entry["iteration"] = item.iteration;
        entry["function_evaluations"] = item.functionEvaluations;
        entry["normalized_equation_residual"] = item.normalizedEquationResidual;
        entry["physical_relative_residual"] = item.physicalRelativeResidual;
        entry["objective"] = item.objective;
        entry["multiplier_mu"] = item.multiplierMu;
        entry["accepted_geodesic_step"] = item.acceptedGeodesicStep;
        entry["quadric_hash"] = item.quadricHash;
        entry["material_geometry_hash"] = item.materialGeometryHash;
        entry["material_topology_hash"] = item.materialTopologyHash;
        entry["raw_equation_residual"] = item.rawEquationResidual
            ? py::cast(*item.rawEquationResidual) : py::none();
        if (item.quadricCoefficients)
        {
            Eigen::Matrix<double, 10, 1> coefficients =
                Eigen::Map<const Eigen::Matrix<double, 10, 1>>(
                    item.quadricCoefficients->data()
                );
            entry["quadric_coefficients"] = coefficients;
        }
        else
        {
            entry["quadric_coefficients"] = py::none();
        }
        trace.append(std::move(entry));
    }
    result["trace"] = std::move(trace);

    py::list events;
    for (const auto& event : solved.numerical.events)
    {
        py::dict entry;
        entry["type"] = event.type;
        entry["restart"] = event.restart;
        entry["iteration"] = event.iteration;
        entry["evaluation_index"] = event.evaluationIndex;
        entry["restart_evaluation_index"] = event.restartEvaluationIndex;
        entry["step"] = event.step ? py::cast(*event.step) : py::none();
        entry["outcome"] = event.outcome;
        entry["reason"] = event.reason;
        entry["detail"] = event.detail;
        entry["quadric_hash"] = event.quadricHash;
        entry["normalized_equation_residual"] = event.normalizedEquationResidual
            ? py::cast(*event.normalizedEquationResidual) : py::none();
        entry["physical_relative_residual"] = event.physicalRelativeResidual
            ? py::cast(*event.physicalRelativeResidual) : py::none();
        entry["objective"] = event.objective
            ? py::cast(*event.objective) : py::none();
        entry["multiplier_mu"] = event.multiplierMu
            ? py::cast(*event.multiplierMu) : py::none();
        entry["material_geometry_hash"] = event.materialGeometryHash
            ? py::cast(*event.materialGeometryHash) : py::none();
        entry["material_topology_hash"] = event.materialTopologyHash
            ? py::cast(*event.materialTopologyHash) : py::none();
        entry["equation_profile"] =
            spin::spinEquationProfileToken(event.equationProfile);
        entry["solver_coordinate_identity"] =
            spin::spinSolverCoordinateIdentity(event.equationProfile);
        entry["raw_equation_residual"] = event.rawEquationResidual
            ? py::cast(*event.rawEquationResidual) : py::none();
        if (event.quadricCoefficients)
        {
            Eigen::Matrix<double, 10, 1> coefficients =
                Eigen::Map<const Eigen::Matrix<double, 10, 1>>(
                    event.quadricCoefficients->data()
                );
            entry["quadric_coefficients"] = coefficients;
        }
        else
        {
            entry["quadric_coefficients"] = py::none();
        }
        if (event.evaluationPayload)
            entry["evaluation_payload"] =
                evaluationPayloadToPython(*event.evaluationPayload);
        if (event.newtonSystemPayload)
            entry["newton_system_payload"] =
                newtonSystemPayloadToPython(*event.newtonSystemPayload);
        events.append(std::move(entry));
    }
    result["events"] = std::move(events);

    py::list evaluationTimings;
    for (const auto& timing : solved.numerical.evaluationTimings)
    {
        py::dict entry;
        entry["evaluation_index"] = timing.evaluationIndex;
        entry["restart"] = timing.restart;
        entry["restart_evaluation_index"] = timing.restartEvaluationIndex;
        entry["evaluator_returned"] = timing.evaluatorReturned;
        entry["wall_seconds"] = timing.wallSeconds;
        entry["reported_residual_seconds"] = timing.reportedResidualSeconds
            ? py::cast(*timing.reportedResidualSeconds) : py::none();
        entry["reported_derivative_seconds"] = timing.reportedDerivativeSeconds
            ? py::cast(*timing.reportedDerivativeSeconds) : py::none();
        evaluationTimings.append(std::move(entry));
    }
    result["evaluation_timings"] = std::move(evaluationTimings);

    if (solved.finalGeometryVerified && solved.finalEvaluation)
    {
        const auto& evaluation = *solved.finalEvaluation;
        result["final_geometry_hash"] = evaluation.geometryHash;
        result["final_derivative_defined"] = evaluation.coefficientJacobian.has_value();
        result["final_derivative_diagnostic"] = evaluation.diagnostic;
        result["final_variable_raw_moments"] =
            Eigen::VectorXd(evaluation.variableMomentsTarget.values());
        result["final_total_raw_moments"] =
            Eigen::VectorXd(evaluation.totalMomentsTarget.values());
        result["material_vertices"] = Eigen::MatrixXd(
            evaluation.materialBoundaryTarget.vertices
        );
        result["material_faces"] = Eigen::MatrixXi(
            evaluation.materialBoundaryTarget.faces
        );
        if(!evaluation.sourceExactVerticesTarget.empty())
            result["rational_authority"]=rationalAuthorityToPython(evaluation);
    }
    else
    {
        result["final_geometry_hash"] = py::none();
        result["final_derivative_defined"] = false;
        result["final_derivative_diagnostic"] = solved.verificationDiagnostic;
        result["final_variable_raw_moments"] = py::none();
        result["final_total_raw_moments"] = py::none();
        result["material_vertices"] = py::none();
        result["material_faces"] = py::none();
    }
    if (solved.independentTotalMoments)
    {
        result["independent_total_raw_moments"] =
            Eigen::VectorXd(solved.independentTotalMoments->values());
    }
    else
    {
        result["independent_total_raw_moments"] = py::none();
    }
    return result;
}
#endif


py::array requireArray(const py::handle& value, const char* name)
{
    py::array array = py::array::ensure(value);
    if (!array)
    {
        throw py::type_error(std::string(name) + " must be array-like");
    }
    return array;
}


void requireTriangleMatrixShape(
    const py::array& array,
    const char* name
)
{
    if (array.ndim() != 2 || array.shape(1) != 3)
    {
        throw py::value_error(
            std::string(name) + " must have shape (N, 3)"
        );
    }
    if (
        array.shape(0) >
        static_cast<py::ssize_t>(std::numeric_limits<Eigen::Index>::max())
    )
    {
        throw py::value_error(std::string(name) + " has too many rows");
    }
}


Eigen::MatrixXd verticesFromPython(const py::handle& value)
{
    const py::array original = requireArray(value, "vertices");
    requireTriangleMatrixShape(original, "vertices");
    const char kind = original.dtype().kind();
    if (kind != 'f' && kind != 'i' && kind != 'u')
    {
        throw py::type_error(
            "vertices must use a real floating-point or integer dtype"
        );
    }

    using DoubleArray = py::array_t<
        double,
        py::array::c_style | py::array::forcecast
    >;
    DoubleArray converted = DoubleArray::ensure(original);
    if (!converted)
    {
        throw py::type_error("vertices cannot be represented as float64");
    }
    const auto input = converted.unchecked<2>();
    Eigen::MatrixXd vertices(
        static_cast<Eigen::Index>(input.shape(0)),
        3
    );
    for (py::ssize_t row = 0; row < input.shape(0); ++row)
    {
        for (py::ssize_t column = 0; column < 3; ++column)
        {
            vertices(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = input(row, column);
        }
    }
    return vertices;
}


#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
double exactClearanceFromPython(const py::handle& value) {
    const py::array original=requireArray(value,"minimum clearance");
    const char kind=original.dtype().kind();
    if (original.ndim()!=0 || original.itemsize()>8 || (kind!='f' && kind!='i' && kind!='u'))
        throw py::type_error("minimum clearance must be a real scalar losslessly representable as binary64; bool is invalid");
    const auto converted=py::array_t<double,py::array::forcecast>::ensure(original);
    const double result=*converted.data();
    if (!std::isfinite(result) || result<0)
        throw py::value_error("minimum clearance must be finite and nonnegative");
    if (kind=='i') {
        const auto integers=py::array_t<std::int64_t,py::array::forcecast>::ensure(original);
        if (result>=0x1p63 || static_cast<std::int64_t>(result)!=*integers.data())
            throw py::value_error("minimum clearance integer is not exactly representable as binary64");
    } else if (kind=='u') {
        const auto integers=py::array_t<std::uint64_t,py::array::forcecast>::ensure(original);
        if (result>=0x1p64 || static_cast<std::uint64_t>(result)!=*integers.data())
            throw py::value_error("minimum clearance integer is not exactly representable as binary64");
    }
    return result;
}

Eigen::MatrixXd exactVerticesFromPython(const py::handle& value)
{
    if (!py::isinstance<py::array>(value))
    {
        return exact_mesh_ingestion::verticesFromSequence(value);
    }

    const py::array original=requireArray(value,"vertices");
    requireTriangleMatrixShape(original,"vertices");
    const char kind=original.dtype().kind();
    if(original.itemsize()>8)
        throw py::type_error("exact boundary accepts at most binary64 coordinates; higher precision must not be silently rounded");
    Eigen::MatrixXd result=verticesFromPython(original);
    // Float16/32 -> binary64 is exact. For integer inputs verify that conversion
    // preserved the integer before geometry or hashing; never cast out of range.
    if(kind=='i') {
        const auto integers=py::array_t<std::int64_t,py::array::c_style|py::array::forcecast>::ensure(original);
        const auto input=integers.unchecked<2>();
        for(Eigen::Index i=0;i<result.rows();++i) for(int j=0;j<3;++j) {
            const double d=result(i,j);
            if(d < -0x1p63 || d >= 0x1p63 || static_cast<std::int64_t>(d)!=input(i,j))
                throw py::value_error("integer coordinate is not exactly representable as binary64");
        }
    } else if(kind=='u') {
        const auto integers=py::array_t<std::uint64_t,py::array::c_style|py::array::forcecast>::ensure(original);
        const auto input=integers.unchecked<2>();
        for(Eigen::Index i=0;i<result.rows();++i) for(int j=0;j<3;++j) {
            const double d=result(i,j);
            if(d >= 0x1p64 || static_cast<std::uint64_t>(d)!=input(i,j))
                throw py::value_error("integer coordinate is not exactly representable as binary64");
        }
    }
    if (!result.allFinite())
        throw py::value_error("vertex coordinates must be finite");
    return result;
}


#endif

template <typename Integer>
Eigen::MatrixXi checkedFaceConversion(
    const py::array& original,
    const char* rangeMessage
)
{
    using IntegerArray = py::array_t<
        Integer,
        py::array::c_style | py::array::forcecast
    >;
    IntegerArray converted = IntegerArray::ensure(original);
    if (!converted)
    {
        throw py::type_error("faces cannot be represented as integers");
    }
    const auto input = converted.template unchecked<2>();
    Eigen::MatrixXi faces(
        static_cast<Eigen::Index>(input.shape(0)),
        3
    );
    for (py::ssize_t row = 0; row < input.shape(0); ++row)
    {
        for (py::ssize_t column = 0; column < 3; ++column)
        {
            const Integer value = input(row, column);
            if constexpr (std::numeric_limits<Integer>::is_signed)
            {
                if (
                    value < static_cast<Integer>(
                        std::numeric_limits<int>::min()
                    ) ||
                    value > static_cast<Integer>(
                        std::numeric_limits<int>::max()
                    )
                )
                {
                    throw py::value_error(rangeMessage);
                }
            }
            else if (
                value > static_cast<Integer>(
                    std::numeric_limits<int>::max()
                )
            )
            {
                throw py::value_error(rangeMessage);
            }
            faces(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = static_cast<int>(value);
        }
    }
    return faces;
}


Eigen::MatrixXi facesFromPython(const py::handle& value)
{
    const py::array original = requireArray(value, "faces");
    requireTriangleMatrixShape(original, "faces");
    const char kind = original.dtype().kind();
    if (kind != 'i' && kind != 'u')
    {
        throw py::type_error("faces must use an integer dtype");
    }
    if (
        original.itemsize() >
        static_cast<py::ssize_t>(sizeof(unsigned long long))
    )
    {
        throw py::value_error("face indices exceed supported integer width");
    }

    constexpr const char* rangeMessage =
        "each face index must be representable as a C++ int";
    if (kind == 'i')
    {
        return checkedFaceConversion<long long>(original, rangeMessage);
    }
    return checkedFaceConversion<unsigned long long>(
        original,
        rangeMessage
    );
}


#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
Eigen::MatrixXi exactFacesFromPython(const py::handle& value)
{
    if (!py::isinstance<py::array>(value))
    {
        return exact_mesh_ingestion::facesFromSequence(value);
    }
    return facesFromPython(value);
}
#endif


std::vector<int> cellIdsFromPython(const py::handle& value)
{
    const py::array original = requireArray(value, "cell_ids_z_fast");
    if (original.ndim() != 1)
    {
        throw py::value_error("cell_ids_z_fast must have shape (N^3,)");
    }
    const char kind = original.dtype().kind();
    if (kind != 'i' && kind != 'u')
    {
        throw py::type_error("cell_ids_z_fast must use an integer dtype");
    }
    using IntegerArray = py::array_t<
        long long,
        py::array::c_style | py::array::forcecast
    >;
    IntegerArray converted = IntegerArray::ensure(original);
    if (!converted)
    {
        throw py::type_error(
            "cell_ids_z_fast cannot be represented as signed integers"
        );
    }
    const auto input = converted.unchecked<1>();
    std::vector<int> result(static_cast<std::size_t>(input.shape(0)));
    for (py::ssize_t index = 0; index < input.shape(0); ++index)
    {
        const long long item = input(index);
        if (
            item < static_cast<long long>(std::numeric_limits<int>::min()) ||
            item > static_cast<long long>(std::numeric_limits<int>::max())
        )
        {
            throw py::value_error(
                "each cell ID must be representable as a C++ int"
            );
        }
        result[static_cast<std::size_t>(index)] = static_cast<int>(item);
    }
    return result;
}


Eigen::MatrixXd weightsFromPython(const py::handle& value)
{
    const py::array original = requireArray(value, "beta_node_weights_z_fast");
    if (original.ndim() != 2 || original.shape(1) <= 0)
    {
        throw py::value_error(
            "beta_node_weights_z_fast must have shape (nodes, handles)"
        );
    }
    const char kind = original.dtype().kind();
    if (kind != 'f' && kind != 'i' && kind != 'u')
    {
        throw py::type_error(
            "beta_node_weights_z_fast must use a real numeric dtype"
        );
    }
    using DoubleArray = py::array_t<
        double,
        py::array::c_style | py::array::forcecast
    >;
    DoubleArray converted = DoubleArray::ensure(original);
    if (!converted)
    {
        throw py::type_error(
            "beta_node_weights_z_fast cannot be represented as float64"
        );
    }
    const auto input = converted.unchecked<2>();
    Eigen::MatrixXd result(
        static_cast<Eigen::Index>(input.shape(0)),
        static_cast<Eigen::Index>(input.shape(1))
    );
    for (py::ssize_t row = 0; row < input.shape(0); ++row)
    {
        for (py::ssize_t column = 0; column < input.shape(1); ++column)
        {
            result(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = input(row, column);
        }
    }
    return result;
}


#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER

py::dict staticPlaneResultToPython(const spin::PlaneStaticBalanceResult& value,
                                 const spin::SpinDesignRequest& request)
{
    py::dict output;
    output["status"] = spin::planeStaticBalanceStatusToken(value.status);
    output["diagnostic"] = value.diagnostic;
    output["plane_coefficients"] = value.coefficients;
    output["material_region"] = value.evaluation.materialRegion==spin::PlaneMaterialRegion::StrictFull?
        "STRICT_FULL":value.evaluation.materialRegion==spin::PlaneMaterialRegion::StrictEmpty?
        "STRICT_EMPTY":"UNCLASSIFIED";
    output["boundary_optimality"] = py::none();
    if(value.boundaryOptimality) {
        py::dict witness;
        witness["kind"]="STRICT_FULL_DOMAIN_ZERO_EQUALITY_MULTIPLIERS_V1";
        witness["objective_measure"]=spin::staticObjectiveMeasureToken(request.staticObjectiveMeasure);
        witness["plane_coefficients"]=value.boundaryOptimality->coefficients;
        witness["minimum_gradient_polynomial"]=value.boundaryOptimality->minimumGradientPolynomial;
        output["boundary_optimality"]=std::move(witness);
    }
    output["lambda"] = value.lambda;
    output["constraints"] = value.constraints;
    output["normalized_constraints"] = value.normalizedConstraints;
    output["normalized_constraints_norm"] = value.normalizedConstraints.norm();
    output["normalized_stationarity_residual"] = value.stationarityResidual;
    output["constraint_scale"] = value.constraintScale;
    output["variable_raw_moments"] = Eigen::VectorXd(value.evaluation.variableMomentsTarget.values());
    output["total_raw_moments"] = Eigen::VectorXd(value.evaluation.totalMomentsTarget.values());
    output["paper_optimality_conditions_verified"] = value.paperOptimalityConditionsVerified;
    output["iterations"] = value.iterations;
    output["accepted_steps"] = value.acceptedSteps;
    output["variable_material_vertices"] = value.evaluation.variableMaterialBoundaryTarget.vertices;
    output["variable_material_faces"] = value.evaluation.variableMaterialBoundaryTarget.faces;
    py::list trace;
    for (const auto& row : value.trace) {
        py::dict event;
        event["evaluation"] = row.evaluation;
        event["iteration"] = row.iteration;
        event["lambda"] = row.lambda;
        event["plane_coefficients"] = row.coefficients;
        event["constraints"] = row.constraints;
        event["normalized_constraints"] = row.normalizedConstraints;
        event["normalized_stationarity_residual"] = row.normalizedStationarityResidual;
        event["residual_norm"] = row.residualNorm;
        event["variable_mass"] = row.variableMass;
        event["step"] = row.step;
        event["accepted"] = row.accepted;
        event["evaluation_succeeded"] = row.evaluationSucceeded;
        event["derivative_available"] = row.derivativeAvailable;
        event["diagnostic"] = row.diagnostic;
        trace.append(std::move(event));
    }
    output["trace"] = std::move(trace);
    return output;
}

physics::RawMassMoments10 rawMomentsFromVector(const Eigen::VectorXd& values)
{
    if (values.size() != physics::RawMassMoments10::Count ||
        !values.allFinite())
    {
        throw py::value_error("raw moments must contain 10 finite values");
    }
    MassMoments moments;
    moments.mass = values(0);
    moments.firstMoment = values.segment<3>(1);
    moments.secondMoment <<
        values(4), values(9), values(8),
        values(9), values(5), values(7),
        values(8), values(7), values(6);
    return physics::RawMassMoments10::fromMassMoments(moments);
}


std::vector<spin::QuadricCoefficients10> quadricsFromRows(
    const Eigen::MatrixXd& rows
)
{
    if (rows.cols() != 10 || !rows.allFinite())
    {
        throw py::value_error(
            "initial quadrics must have shape (N, 10) with finite values"
        );
    }
    std::vector<spin::QuadricCoefficients10> result;
    result.reserve(static_cast<std::size_t>(rows.rows()));
    for (Eigen::Index row = 0; row < rows.rows(); ++row)
    {
        result.push_back(spin::QuadricCoefficients10::fromUnnormalized(
            rows.row(row).transpose()
        ));
    }
    return result;
}

#endif


} // namespace


PYBIND11_MODULE(_make_it_stand, module)
{
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    registerSourceSpinMaterialBindings(module,exactVerticesFromPython,exactFacesFromPython);
#endif
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    py::class_<verification::ExactMaterialRegion>(module, "ExactMaterialRegion",
        "Immutable exact represented material set; no public constructor or physical/seal authority.")
        .def_property_readonly("facts", [](const verification::ExactMaterialRegion& region) {
            const auto& r=region.facts();
            py::dict out;
            out["algorithm"]=r.algorithm; out["backend"]=r.backend; out["semantics"]=r.semantics;
            out["input_digest"]=r.inputDigest; out["region_digest"]=r.regionDigest;
            out["completed"]=r.completed; out["valid"]=r.valid; out["error"]=r.error;
            out["explicit_empty"]=r.explicitEmpty;
            out["parents"]=r.parents; out["depths"]=r.depths; out["original_signs"]=r.originalSigns;
            out["checker_reversed_components"]=r.checkerReversedComponents;
            out["material_components"]=r.materialComponents; out["bounded_void_components"]=r.boundedVoidComponents;
            out["scope"]="exact-represented-polyhedral-set-only";
            return out;
        })
        .def("matches_rational_boundary", [](const verification::ExactMaterialRegion& region,
                const std::vector<std::array<std::string,3>>& vertices,const py::object& faces) {
            return region.matchesRationalBoundary(vertices,exactFacesFromPython(faces));
        })
        .def_property_readonly("rational_boundary", [](const verification::ExactMaterialRegion& region) {
            py::dict result;
            result["vertices"]=py::cast(region.rationalBoundaryVertices());
            result["faces"]=py::cast(region.boundaryFaces());
            return result;
        })
        .def_property_readonly("boundary_coordinate_policy",[](const verification::ExactMaterialRegion& region) {
            return region.boundaryUsesRationals()?"RATIONAL_INDEXED":"BINARY64_INDEXED";
        })
        .def_property_readonly("binary64_boundary_vertices",&verification::ExactMaterialRegion::binary64BoundaryVertices)
        .def_property_readonly("construction_operands",&verification::ExactMaterialRegion::constructionOperands);
    module.def("make_exact_interior_difference_region", &verification::makeExactInteriorDifferenceRegion,
        py::arg("outer"),py::arg("carved"),
        "Construct a strictly interior difference from owned exact operands and an independently checked shared boundary.");
    module.def("make_exact_material_region", [](const py::object& v,const py::object& f) {
        return verification::makeExactMaterialRegion(exactVerticesFromPython(v),exactFacesFromPython(f));
    },py::arg("vertices"),py::arg("faces"),
       "Authenticate an oriented finite regular closed region from exact binary64 indexed triangles. "
       "Inputs are NumPy matrices or nested list/tuple matrices whose rows may also be one-dimensional NumPy arrays. "
       "Negative checker copies only; input geometry is preserved. Empty mesh is malformed.");
    module.def("make_exact_material_region_from_rationals",
        [](const std::vector<std::array<std::string,3>>& v,const py::object& f) {
            return verification::makeExactMaterialRegionFromRationals(v,exactFacesFromPython(f));
        },py::arg("vertices"),py::arg("faces"),
        "Authenticate an exact rational indexed-triangle region without binary64 coordinate conversion. "
        "Coordinates must be decimal numerator/positive-denominator strings; set operations only.");
    module.def("empty_exact_material_region", &verification::emptyExactMaterialRegion,
               "Explicit empty algebra operand, never a valid physical outer artifact.");
    module.def("relate_exact_material_regions", [](const verification::ExactMaterialRegion& a,
                                                  const verification::ExactMaterialRegion& b) {
        const auto r=verification::relateExactMaterialRegions(a,b);
        py::dict out;
        out["algorithm"]=r.algorithm; out["semantics"]=r.semantics;
        out["a_digest"]=r.aDigest; out["b_digest"]=r.bDigest; out["relation_digest"]=r.relationDigest;
        out["a_explicit_empty"]=r.aExplicitEmpty; out["b_explicit_empty"]=r.bExplicitEmpty;
        out["completed"]=r.completed; out["error"]=r.error;
        out["intersection_nonempty"]=r.intersectionNonempty; out["positive_volume_overlap"]=r.positiveVolumeOverlap;
        out["a_subset_b"]=r.aSubsetB; out["b_subset_a"]=r.bSubsetA; out["equal"]=r.equal;
        out["symmetric_difference_volume_exact"]=r.symmetricDifferenceVolumeExact;
        out["symmetric_difference_measure_scope"]="exact-regularized-volume-not-Hausdorff-bound";
        out["scope"]="exact-represented-polyhedral-set-only";
        return out;
    },py::arg("a"),py::arg("b"),
       "Exact inclusive contact/containment and positive-volume overlap. Invalid operands throw; "
       "failed Booleans return unknown facts. Touching does not imply volume overlap or invalid density regions.");
    module.def("relate_exact_boolean_operation_to_material",
        [](const verification::ExactMaterialRegion& design,
           const verification::ExactMaterialRegion& cutter,
           const verification::ExactMaterialRegion& material,
           const std::string& operation) {
            verification::ExactMaterialBooleanOperation selected;
            if(operation=="intersection")
                selected=verification::ExactMaterialBooleanOperation::Intersection;
            else if(operation=="difference")
                selected=verification::ExactMaterialBooleanOperation::Difference;
            else throw std::invalid_argument("unsupported exact Boolean operation");
            const auto r=verification::relateExactBooleanOperationToMaterial(
                design,cutter,material,selected);
            py::dict out;
            out["algorithm"]=r.algorithm;out["semantics"]=r.semantics;
            out["a_digest"]=r.aDigest;out["b_digest"]=r.bDigest;
            out["relation_digest"]=r.relationDigest;
            out["a_explicit_empty"]=r.aExplicitEmpty;
            out["b_explicit_empty"]=r.bExplicitEmpty;
            out["completed"]=r.completed;out["error"]=r.error;
            out["intersection_nonempty"]=r.intersectionNonempty;
            out["positive_volume_overlap"]=r.positiveVolumeOverlap;
            out["a_subset_b"]=r.aSubsetB;out["b_subset_a"]=r.bSubsetA;
            out["equal"]=r.equal;
            out["symmetric_difference_volume_exact"]=r.symmetricDifferenceVolumeExact;
            out["scope"]="exact-input-binary64-represented-polyhedral-Boolean-only";
            return out;
        },py::arg("design"),py::arg("cutter"),py::arg("material"),
        py::arg("operation"),
        "Compute an exact input-mesh intersection/difference and compare its represented set with an independently supplied material. No continuous-quadric claim.");
#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
    module.def("libigl_exact_boolean_construction",
        [](const Eigen::MatrixXd& designVertices, const Eigen::MatrixXi& designFaces,
           const Eigen::MatrixXd& cutterVertices, const Eigen::MatrixXi& cutterFaces,
           const std::string& operation) {
            verification::ExactMaterialBooleanOperation selected;
            if (operation == "intersection")
                selected = verification::ExactMaterialBooleanOperation::Intersection;
            else if (operation == "difference")
                selected = verification::ExactMaterialBooleanOperation::Difference;
            else
                throw std::invalid_argument("unsupported exact libigl Boolean operation");
            const auto result = verification::computeExactLibiglBooleanConstruction(
                designVertices, designFaces, cutterVertices, cutterFaces, selected);
            py::list vertices;
            for (const auto& row : result.vertices) {
                py::list coordinates;
                for (const auto& value : row)
                    coordinates.append(value);
                vertices.append(coordinates);
            }
            py::dict out;
            out["succeeded"] = result.succeeded;
            out["vertices"] = vertices;
            out["faces"] = result.faces;
            out["birth_faces"] = result.birthFaces;
            out["scope"] = "pinned-libigl-exact-output-before-binary64-assignment";
            return out;
        },py::arg("design_vertices"),py::arg("design_faces"),
        py::arg("cutter_vertices"),py::arg("cutter_faces"),py::arg("operation"),
        "Diagnostic-only pinned libigl Boolean with rational output coordinates. It does not replace the solver's binary64 material.");
    module.def("libigl_exact_boolean_construction_from_rationals",
        [](const std::vector<std::array<std::string,3>>& designVertices,
           const Eigen::MatrixXi& designFaces,
           const std::vector<std::array<std::string,3>>& cutterVertices,
           const Eigen::MatrixXi& cutterFaces, const std::string& operation) {
            verification::ExactMaterialBooleanOperation selected;
            if (operation == "intersection")
                selected = verification::ExactMaterialBooleanOperation::Intersection;
            else if (operation == "difference")
                selected = verification::ExactMaterialBooleanOperation::Difference;
            else throw std::invalid_argument("unsupported exact libigl Boolean operation");
            const auto result = verification::computeExactLibiglBooleanConstructionFromRationals(
                designVertices, designFaces, cutterVertices, cutterFaces, selected);
            py::dict out;
            out["succeeded"] = result.succeeded;
            out["vertices"] = result.vertices;
            out["faces"] = result.faces;
            out["birth_faces"] = result.birthFaces;
            out["scope"] = "pinned-libigl-rational-operands-and-output";
            return out;
        },py::arg("design_vertices"),py::arg("design_faces"),
        py::arg("cutter_vertices"),py::arg("cutter_faces"),py::arg("operation"),
        "Exact rational polyhedral operands and output, with no binary64 roundtrip. Independent material validation remains required.");
#endif
    module.def("check_exact_material_subset", [](const verification::ExactMaterialRegion& a,
                                                   const verification::ExactMaterialRegion& b) {
        const auto r=verification::checkExactMaterialSubset(a,b);
        py::dict out;
        out["algorithm"]=r.algorithm; out["semantics"]=r.semantics;
        out["a_digest"]=r.aDigest; out["b_digest"]=r.bDigest;
        out["completed"]=r.completed; out["error"]=r.error;
        out["a_subset_b"]=r.aSubsetB;
        out["scope"]="exact-represented-polyhedral-set-only";
        return out;
    },py::arg("a"),py::arg("b"),
       "Exact inclusive represented-set subset predicate; this narrow query does not compute symmetric-difference volume.");
    module.def("check_exact_material_partition", [](const verification::ExactMaterialRegion& whole,
                                                    const std::vector<verification::ExactMaterialRegion>& parts) {
        const auto r=verification::checkExactMaterialPartition(whole,parts);
        py::dict out;
        out["algorithm"]=r.algorithm; out["backend"]=r.backend; out["semantics"]=r.semantics;
        out["empty_protocol"]=r.emptyProtocol; out["whole_digest"]=r.wholeDigest;
        out["part_digests"]=r.partDigests; out["partition_digest"]=r.partitionDigest;
        out["whole_explicit_empty"]=r.wholeExplicitEmpty; out["part_explicit_empty"]=r.partExplicitEmpty;
        out["completed"]=r.completed; out["error"]=r.error;
        out["union_equals_whole"]=r.unionEqualsWhole; out["no_positive_volume_overlap"]=r.noPositiveVolumeOverlap;
        out["missing_whole_content"]=r.missingWholeContent; out["excess_union_content"]=r.excessUnionContent;
        out["part_subset_whole"]=r.partSubsetWhole; out["positive_volume_overlap_pairs"]=r.positiveVolumeOverlapPairs;
        out["complete_partition"]=r.completePartition;
        out["scope"]="exact-represented-polyhedral-partition-only";
        return out;
    },py::arg("whole"),py::arg("parts"),
       "Authenticate inclusive union equality and pairwise disjoint interiors of ordered region handles. "
       "Empty list is missing evidence; explicit empty parts are permitted. No seal or physical authority.");
    module.def("check_exact_boundary_clearance", [](const verification::ExactMaterialRegion& a,
                                                    const verification::ExactMaterialRegion& b,const py::object& threshold) {
        const auto r=verification::checkExactBoundaryClearance(a,b,exactClearanceFromPython(threshold));
        py::dict out;
        out["algorithm"]=r.algorithm; out["backend"]=r.backend; out["semantics"]=r.semantics;
        out["a_digest"]=r.aDigest; out["b_digest"]=r.bDigest;
        out["a_input_digest"]=r.aInputDigest; out["b_input_digest"]=r.bInputDigest;
        out["a_explicit_empty"]=r.aExplicitEmpty; out["b_explicit_empty"]=r.bExplicitEmpty;
        out["minimum_clearance_exact"]=r.minimumClearanceExact; out["clearance_digest"]=r.clearanceDigest;
        out["completed"]=r.completed; out["applicable"]=r.applicable; out["state"]=r.state; out["error"]=r.error;
        out["squared_distance_exact"]=r.squaredDistanceExact; out["witness_faces"]=r.witnessFaces;
        out["at_least_minimum"]=r.atLeastMinimum; out["contact"]=r.contact;
        out["scope"]="exact-represented-boundary-distance-only";
        return out;
    },py::arg("a"),py::arg("b"),py::arg("minimum_clearance"),
       "Exact boundary distance >= explicitly supplied represented threshold; equality passes and contact is distinct. "
       "Explicit empty operands are NOT_APPLICABLE with absent distance/witness/comparison. "
       "This alone does not establish containment, valid shell thickness or a seal.");
    const auto exactBoundaryMode=[](const std::string& orientation) {
            using Mode=verification::ExactBoundaryOrientation;
            if(orientation=="OutwardPositive") return Mode::OutwardPositive;
            if(orientation=="InwardNegative") return Mode::InwardNegative;
            if(orientation=="AnyNonzero") return Mode::AnyNonzero;
            throw py::value_error("orientation must be OutwardPositive, InwardNegative, or AnyNonzero");
    };
    const auto exactBoundaryFactsToPython=[](const verification::ExactBoundaryFacts& r,
                                           const std::string& encoding) {
            py::dict out;
            out["algorithm"]=r.algorithm; out["backend"]=r.backend;
            out["input_digest"]=r.inputDigest; out["input_digest_algorithm"]="sha256";
            out["input_encoding"]=encoding;
            out["orientation"]=r.orientation;
            out["boundary_valid"]=r.boundaryValid;
            out["closed"]=r.closed; out["edge_manifold"]=r.edgeManifold;
            out["vertex_manifold"]=r.vertexManifold;
            out["consistently_oriented"]=r.consistentlyOriented;
            out["component_orientations_valid"]=r.componentOrientationsValid;
            out["unreferenced_vertices"]=r.unreferencedVertices;
            out["degenerate_faces"]=r.degenerateFaces;
            out["boundary_edges"]=r.boundaryEdges; out["non_manifold_edges"]=r.nonManifoldEdges;
            out["inconsistent_edges"]=r.inconsistentEdges;
            out["non_manifold_vertices"]=r.nonManifoldVertices;
            out["face_components"]=r.faceComponents;
            out["component_volume_signs"]=r.componentVolumeSigns;
            out["component_signed_volumes_exact"]=r.componentSignedVolumesExact;
            out["self_intersection_completed"]=r.intersectionCompleted;
            out["self_intersection_free"]=r.intersectionFree;
            out["intersection_error"]=r.intersectionError;
            out["intersecting_pairs"]=r.intersectingPairs;
            out["intersecting_pair_count"]=r.intersectingPairs.size();
            out["exact_pair_tests"]=r.exactPairTests;
            out["scope"]="single-represented-boundary-only";
            out["material_domain_status"]="NOT_VERIFIED";
            return out;
    };
    module.def("audit_exact_rational_material_boundary",
        [exactBoundaryFactsToPython](const std::vector<std::array<std::string,3>>& vertices,
            const py::object& faces,const py::object& proof) {
            const verification::VerifiedFixedBoundary* fixed=nullptr;
            if(!proof.is_none()) fixed=&proof.cast<const verification::VerifiedFixedBoundary&>();
            const auto r=verification::auditExactRationalMaterialBoundary(vertices,exactFacesFromPython(faces),fixed);
            py::dict out;
            out["algorithm"]=r.algorithm;out["completed"]=r.completed;
            out["classification_valid"]=r.classificationValid;out["error"]=r.error;
            out["boundary"]=exactBoundaryFactsToPython(r.boundary,
                "exact-rational-boundary-input-v1/canonical-fractions/u64-big-endian");
            out["parents"]=r.parents;out["depths"]=r.depths;out["original_signs"]=r.originalSigns;
            out["witness_vertices"]=r.witnessVertices;out["inside"]=r.inside;
            out["material_components"]=r.materialComponents;out["bounded_void_components"]=r.boundedVoidComponents;
            out["fixed_proof_identity"]=r.fixedProofIdentity;out["reused_pair_tests"]=r.reusedPairTests;
            py::dict performance;performance["boundary_seconds"]=r.boundarySeconds;
            performance["nesting_seconds"]=r.nestingSeconds;out["performance"]=performance;
            out["scope"]="exact-embedded-boundary-material-classification-only";
            return out;
        },py::arg("vertices"),py::arg("faces"),py::arg("fixed_proof")=py::none(),
        "Complete exact boundary audit and component point-location nesting/orientation checks; "
        "no Nef set operations, partition, clearance or final artifact authority.");
    module.def("audit_exact_triangle_boundary",
        [exactBoundaryMode,exactBoundaryFactsToPython](const py::object& vertices, const py::object& faces,
                                                     const std::string& orientation) {
            const auto audit=verification::auditExactTriangleBoundary(
                exactVerticesFromPython(vertices),exactFacesFromPython(faces),exactBoundaryMode(orientation));
            return exactBoundaryFactsToPython(audit.facts(),
                "exact-boundary-input-v1/binary64-row-major/u64-big-endian");
        }, py::arg("vertices"), py::arg("faces"), py::arg("orientation"),
        "Exact audit of the binary64 represented indexed boundary. No repair or welding. "
        "Inputs are NumPy matrices or nested list/tuple matrices whose rows may also be one-dimensional NumPy arrays. "
        "Float16/32 promote exactly; integer coordinates must be exactly binary64-representable. "
        "Wider floats are rejected. "
        "AnyNonzero is a boundary orientation mode only. Does not verify containment, "
        "nesting, thickness, material union, CAD, source/build identity, or artifact sealing.");
    module.def("audit_exact_rational_triangle_boundary",
        [exactBoundaryMode,exactBoundaryFactsToPython](const std::vector<std::array<std::string,3>>& vertices,
                                                     const py::object& faces,const std::string& orientation) {
            const auto audit=verification::auditExactTriangleBoundaryRational(
                vertices,exactFacesFromPython(faces),exactBoundaryMode(orientation));
            return exactBoundaryFactsToPython(audit.facts(),
                "exact-rational-boundary-input-v1/canonical-fractions/u64-big-endian");
        },py::arg("vertices"),py::arg("faces"),py::arg("orientation"),
        "Exact topology, winding and self-intersection audit of rational indexed coordinates. "
        "Uses the existing rational boundary predicates without converting to binary64. "
        "This boundary audit does not establish multi-shell nesting or material partition.");
    py::class_<verification::VerifiedFixedBoundary>(module,"VerifiedFixedBoundary")
        .def_property_readonly("reference_identity",[](const verification::VerifiedFixedBoundary& proof){
            return proof.facts().inputDigest;
        });
    module.def("verify_fixed_rational_boundary",
        [](const std::vector<std::array<std::string,3>>& vertices,const py::object& faces) {
            return verification::verifyFixedRationalBoundary(vertices,exactFacesFromPython(faces));
        },py::arg("vertices"),py::arg("faces"),
        "Own an immutable complete exact self-intersection/topology proof for a fixed boundary.");
    module.def("audit_rational_boundary_with_fixed_proof",
        [exactBoundaryMode,exactBoundaryFactsToPython](const std::vector<std::array<std::string,3>>& vertices,
            const py::object& faces,const std::string& orientation,const verification::VerifiedFixedBoundary& proof,
            bool profile) {
            verification::BoundaryPerformance performance;
            const auto result=verification::auditRationalBoundaryWithFixedProof(
                vertices,exactFacesFromPython(faces),exactBoundaryMode(orientation),proof,profile?&performance:nullptr);
            auto out=exactBoundaryFactsToPython(result.audit.facts(),
                "exact-rational-boundary-input-v1/canonical-fractions/u64-big-endian");
            out["reused_fixed_pair_tests"]=result.reusedPairTests;
            out["fixed_proof_identity"]=proof.facts().inputDigest;
            if(profile) {
                py::dict detail;
                detail["exclusive_seconds"]=performance.exclusiveSeconds;
                detail["total_seconds"]=performance.totalSeconds;
                detail["vertex_count"]=performance.vertexCount;
                detail["face_count"]=performance.faceCount;
                detail["matched_original_faces"]=performance.matchedOriginalFaces;
                detail["exact_pair_tests"]=result.audit.facts().exactPairTests;
                detail["reused_pair_tests"]=result.reusedPairTests;
                out["performance"]=detail;
            }
            return out;
        },py::arg("vertices"),py::arg("faces"),py::arg("orientation"),py::arg("proof"),py::arg("profile")=false,
        "Reuse only exact unchanged original triangles with bijective indexed vertex correspondence; "
        "check every new/mixed pair and all current topology, winding and exact volumes.");
#endif
    module.def(
        "audit_triangle_boundary",
        [](
            const py::object& vertices,
            const py::object& faces,
            double geometricTolerance
        )
        {
            const manufacturing::BoundaryAudit audit =
                manufacturing::auditTriangleBoundary(
                    verticesFromPython(vertices),
                    facesFromPython(faces),
                    manufacturing::BoundaryOrientation::OutwardPositive,
                    false,
                    geometricTolerance
                );
            py::dict result;
            result["strict_boundary_valid"] = audit.strictBoundaryValid();
            result["vertex_count"] = audit.vertexCount;
            result["face_count"] = audit.faceCount;
            result["unreferenced_vertex_count"] = audit.unreferencedVertexCount;
            result["degenerate_face_count"] = audit.degenerateFaceCount;
            result["boundary_edge_count"] = audit.boundaryEdgeCount;
            result["non_manifold_edge_count"] = audit.nonManifoldEdgeCount;
            result["inconsistent_edge_count"] = audit.inconsistentEdgeCount;
            result["non_manifold_vertex_count"] = audit.nonManifoldVertexCount;
            result["face_component_count"] = audit.faceComponentCount;
            result["signed_volume"] = audit.signedVolume;
            result["edge_manifold"] = audit.edgeManifold;
            result["vertex_manifold"] = audit.vertexManifold;
            result["closed"] = audit.closed;
            result["consistently_oriented"] = audit.consistentlyOriented;
            result["component_orientations_valid"] =
                audit.componentOrientationsValid;
            result["self_intersection_completed"] =
                audit.selfIntersections.completed;
            result["self_intersection_free"] =
                audit.selfIntersections.intersectionFree;
            result["self_intersection_pair_count"] =
                audit.selfIntersections.intersectingPairCount;
            result["first_self_intersection_faces"] = py::make_tuple(
                audit.selfIntersections.firstFace,
                audit.selfIntersections.secondFace
            );
            return result;
        },
        py::arg("vertices"),
        py::arg("faces"),
        py::arg("geometric_tolerance") = 0.0,
        "Legacy floating boundary diagnostic: shared-corner trimming can miss short "
        "crossings. Its historical strict_boundary_valid flag is not an exact "
        "represented-boundary or full artifact certificate."
    );

    module.def(
        "audit_material_shells",
        [](
            const py::object& outerVertices,
            const py::object& outerFaces,
            const py::object& innerVertices,
            const py::object& innerFaces,
            double requestedThickness,
            double geometricTolerance
        )
        {
            const manufacturing::MaterialShellAudit audit =
                manufacturing::auditMaterialShells(
                    verticesFromPython(outerVertices),
                    facesFromPython(outerFaces),
                    verticesFromPython(innerVertices),
                    facesFromPython(innerFaces),
                    requestedThickness,
                    geometricTolerance
                );
            py::dict clearance;
            clearance["applicable"] = audit.wallClearance.applicable;
            clearance["completed"] = audit.wallClearance.completed;
            clearance["intersection_free"] =
                audit.wallClearance.intersectionFree;
            clearance["meets_requested_thickness"] =
                audit.wallClearance.meetsRequestedThickness;
            clearance["requested_thickness"] =
                audit.wallClearance.requestedThickness;
            clearance["measured_minimum_clearance"] =
                audit.wallClearance.measuredMinimumClearance;
            clearance["witness_outer_face"] = audit.wallClearance.outerFace;
            clearance["witness_inner_face"] = audit.wallClearance.innerFace;
            clearance["broad_phase_pair_count"] =
                audit.wallClearance.broadPhasePairCount;
            clearance["narrow_phase_pair_count"] =
                audit.wallClearance.narrowPhasePairCount;
            py::dict containment;
            containment["applicable"] = audit.containment.applicable;
            containment["completed"] = audit.containment.completed;
            containment["all_cavities_inside_outer"] =
                audit.containment.allCavitiesInsideOuter;
            containment["tested_component_count"] =
                audit.containment.testedComponentCount;
            containment["failed_component"] =
                audit.containment.failedComponent;
            containment["failed_winding_number"] =
                audit.containment.failedWindingNumber;
            containment["nested_cavity_pair_count"] =
                audit.containment.nestedCavityPairCount;
            containment["first_containing_component"] =
                audit.containment.firstContainingComponent;
            containment["first_nested_component"] =
                audit.containment.firstNestedComponent;
            py::dict result;
            result["wall_clearance"] = clearance;
            result["containment"] = containment;
            return result;
        },
        py::arg("outer_vertices"),
        py::arg("outer_faces"),
        py::arg("inner_vertices"),
        py::arg("inner_faces"),
        py::arg("requested_thickness") = 0.0,
        py::arg("geometric_tolerance") = 0.0
    );

    module.doc() =
        "Owning Python bindings for the audited Make It Stand pipeline";

    module.attr("API_SCHEMA_VERSION") = 1;

    py::enum_<manufacturing::ManufacturingMeshFormat>(
        module,
        "ManufacturingMeshFormat"
    )
        .value("OBJ", manufacturing::ManufacturingMeshFormat::Obj)
        .value("ASCII_STL", manufacturing::ManufacturingMeshFormat::AsciiStl)
        .value("BINARY_STL", manufacturing::ManufacturingMeshFormat::BinaryStl)
        .export_values();

#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
    py::enum_<spin::QuadricMomentDerivativeMethod>(module, "SpinDerivativeMethod")
        .value("LEGACY_CONTINUOUS_IMPLICIT_SHAPE", spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape)
        .value("DISCRETE_FULL_RANK_QUADRIC", spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric)
        .value("DISCRETE_REDUCED_STRATUM", spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum)
        .value("DISCRETE_ELLIPSOID", spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid);
    py::enum_<spin::QuadricCapRadiusPolicy>(module, "SpinCapRadiusPolicy")
        .value("REJECT_UNCERTIFIED", spin::QuadricCapRadiusPolicy::RejectUncertified)
        .value("GROW_TO_CERTIFIED_V1", spin::QuadricCapRadiusPolicy::GrowToCertifiedV1);
    module.def(
        "spin_it_faster_quadric_hash",
        [](const Eigen::VectorXd& coefficients) {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error(
                    "normalized quadric coefficients must contain 10 finite values"
                );
            }
            spin::QuadricCoefficients10::Vector fixed;
            fixed = coefficients;
            return spin::QuadricCoefficients10::hashNormalized(fixed);
        },
        py::arg("normalized_coefficients"));
    py::class_<spin::QuadricTriangulationSettings>(module, "SpinTriangulationSettings")
        .def(py::init<>())
        .def_readwrite("icosphere_subdivisions", &spin::QuadricTriangulationSettings::icosphereSubdivisions)
        .def_readwrite("maximum_icosphere_subdivisions", &spin::QuadricTriangulationSettings::maximumIcosphereSubdivisions)
        .def_readwrite("paraboloid_grid_resolution", &spin::QuadricTriangulationSettings::paraboloidGridResolution)
        .def_readwrite("cylinder_angular_resolution", &spin::QuadricTriangulationSettings::cylinderAngularResolution)
        .def_readwrite("hyperbolic_cylinder_grid_resolution", &spin::QuadricTriangulationSettings::hyperbolicCylinderGridResolution)
        .def_readwrite("parabolic_cylinder_grid_resolution", &spin::QuadricTriangulationSettings::parabolicCylinderGridResolution)
        .def_readwrite("truncation_radius", &spin::QuadricTriangulationSettings::truncationRadius)
        .def_readwrite("maximum_radius_doublings", &spin::QuadricTriangulationSettings::maximumRadiusDoublings)
        .def_readwrite("retain_fixed_normal_provenance", &spin::QuadricTriangulationSettings::retainFixedNormalProvenance)
        .def_readwrite("paper_projective_hyperboloid", &spin::QuadricTriangulationSettings::paperProjectiveHyperboloid)
        .def_readwrite("algebraic_finite_map", &spin::QuadricTriangulationSettings::algebraicFiniteMap)
        .def_readwrite("fixed_normal_paraboloid", &spin::QuadricTriangulationSettings::fixedNormalParaboloid)
        .def_readwrite("rational_source_coordinates", &spin::QuadricTriangulationSettings::rationalSourceCoordinates)
        .def_readwrite("compact_rational_source_coordinates", &spin::QuadricTriangulationSettings::compactRationalSourceCoordinates);
    module.def(
        "evaluate_spin_it_faster_source_motion",
        [](
            const Eigen::VectorXd& coefficients,
            const py::object& designVertices,
            const spin::QuadricTriangulationSettings& settings
        )
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error("quadric coefficients must contain 10 finite values");
            }
            const auto result = spin::evaluateQuadricSourceMotion(
                spin::QuadricCoefficients10::fromUnnormalized(coefficients),
                verticesFromPython(designVertices), settings);
            py::dict output;
            output["defined"] = result.defined;
            output["diagnostic"] = result.diagnostic;
            output["scientific_hash"] = result.scientificHash;
            if (result.defined)
            {
                output["vertices"] = result.triangulation.closedBoundary.vertices;
                output["faces"] = result.triangulation.closedBoundary.faces;
                output["face_is_cap"] = result.triangulation.faceIsCap;
                output["cap_face_count"] = result.triangulation.capFaceCount;
                try
                {
                    output["side_selection_defined"] = true;
                    output["interior_is_positive"] =
                        spin::quadricSolidInteriorIsPositiveRegion(
                            spin::QuadricCoefficients10::fromUnnormalized(coefficients),
                            result.triangulation);
                    output["side_selection_method"] =
                        "interior-centroid-winding-or-oriented-interface-v1";
                    output["side_selection_diagnostic"] = "";
                }
                catch (const std::exception& error)
                {
                    output["side_selection_defined"] = false;
                    output["side_selection_diagnostic"] = error.what();
                }
                output["triangulation_hash"] = result.triangulation.hash;
                output["paper_projective_hyperboloid"] = result.triangulation.paperProjectiveHyperboloid;
                output["algebraic_finite_map"] = settings.algebraicFiniteMap;
                output["source_normal_vertices"] = result.triangulation.sourceNormalVertexByVertex;
                output["source_normal_edges"] = result.triangulation.sourceNormalEdgeByVertex;
                output["source_projective_rays"] = result.triangulation.sourceProjectiveRayByVertex;
                output["surface_type"] = spin::quadricSurfaceTypeToken(
                    result.triangulation.type
                );
                output["truncation_radius"] = result.triangulation.truncationRadius;
                output["initial_icosphere_subdivisions"] =
                    result.triangulation.initialIcosphereSubdivisions;
                output["effective_icosphere_subdivisions"] =
                    result.triangulation.effectiveIcosphereSubdivisions;
                output["interface_vertex_count"] = result.triangulation.interfaceVertexCount;
                output["provenance"] = result.provenance;
                Eigen::MatrixXd velocity(
                    static_cast<Eigen::Index>(result.velocity.size()) * 3, 10);
                for (std::size_t vertex = 0; vertex < result.velocity.size(); ++vertex)
                {
                    velocity.block<3, 10>(static_cast<Eigen::Index>(vertex) * 3, 0) =
                        result.velocity[vertex];
                }
                output["vertex_velocity"] = std::move(velocity);
            }
            return output;
        },
        py::arg("coefficients"), py::arg("design_vertices"),
        py::arg("settings") = spin::QuadricTriangulationSettings());
    module.def(
        "evaluate_spin_it_faster_stratum_source_motion",
        [](const Eigen::VectorXd& coefficients, const py::object& designVertices,
           const spin::QuadricTriangulationSettings& settings)
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
                throw py::value_error("quadric coefficients must contain 10 finite values");
            const auto q = spin::QuadricCoefficients10::fromUnnormalized(coefficients);
            const auto result = spin::evaluateQuadricStratumSourceMotion(
                q, verticesFromPython(designVertices), settings);
            py::dict output;
            output["defined"] = result.defined;
            output["diagnostic"] = result.diagnostic;
            output["scientific_hash"] = result.scientificHash;
            output["derivative_scope"] = "REDUCED_FIXED_STRATUM_NOT_FULL_EQ14";
            if (result.defined)
            {
                const Eigen::Index dimension = result.coefficientBasis.cols();
                output["coefficient_basis"] = result.coefficientBasis;
                output["vertices"] = result.triangulation.closedBoundary.vertices;
                output["faces"] = result.triangulation.closedBoundary.faces;
                output["face_is_cap"] = result.triangulation.faceIsCap;
                output["cap_face_count"] = result.triangulation.capFaceCount;
                output["triangulation_hash"] = result.triangulation.hash;
                output["surface_type"] = spin::quadricSurfaceTypeToken(result.triangulation.type);
                output["truncation_radius"] = result.triangulation.truncationRadius;
                output["initial_icosphere_subdivisions"] = result.triangulation.initialIcosphereSubdivisions;
                output["effective_icosphere_subdivisions"] = result.triangulation.effectiveIcosphereSubdivisions;
                output["interface_vertex_count"] = result.triangulation.interfaceVertexCount;
                output["provenance"] = result.provenance;
                if (result.triangulation.fixedNormalParaboloid)
                {
                    output["fixed_normal_paraboloid"] = true;
                    output["fixed_normals"] = result.triangulation.fixedNormals;
                    output["source_normal_vertices"] = result.triangulation.sourceNormalVertexByVertex;
                }
                try
                {
                    output["side_selection_defined"] = true;
                    output["interior_is_positive"] = spin::quadricSolidInteriorIsPositiveRegion(
                        q, result.triangulation);
                    output["side_selection_method"] = "interior-centroid-winding-or-oriented-interface-v1";
                    output["side_selection_diagnostic"] = "";
                }
                catch (const std::exception& error)
                {
                    output["side_selection_defined"] = false;
                    output["side_selection_diagnostic"] = error.what();
                }
                Eigen::MatrixXd velocity(
                    static_cast<Eigen::Index>(result.vertexVelocity.size()) * 3, dimension);
                for (std::size_t vertex = 0; vertex < result.vertexVelocity.size(); ++vertex)
                    velocity.block(static_cast<Eigen::Index>(vertex) * 3, 0, 3, dimension) =
                        result.vertexVelocity[vertex];
                output["vertex_velocity"] = std::move(velocity);
            }
            return output;
        },
        py::arg("coefficients"), py::arg("design_vertices"),
        py::arg("settings") = spin::QuadricTriangulationSettings());
    py::class_<spin::QuadricMomentEvaluatorSettings>(module, "SpinEvaluatorSettings")
        .def(py::init<>())
        .def_readwrite("triangulation", &spin::QuadricMomentEvaluatorSettings::triangulation)
        .def_readwrite("minimum_interface_gradient", &spin::QuadricMomentEvaluatorSettings::minimumInterfaceGradient)
        .def_readwrite("derivative_method", &spin::QuadricMomentEvaluatorSettings::derivativeMethod)
        .def_readwrite("diagnostic_stage_observer", &spin::QuadricMomentEvaluatorSettings::diagnosticStageObserver)
        .def_readwrite("cap_radius_policy", &spin::QuadricMomentEvaluatorSettings::capRadiusPolicy)
        .def_readwrite("exact_binary64_self_intersection", &spin::QuadricMomentEvaluatorSettings::exactBinary64SelfIntersection)
        .def_readwrite("exact_binary64_source_boundary", &spin::QuadricMomentEvaluatorSettings::exactBinary64SourceBoundary)
        .def_readwrite("exact_binary64_material_boundary", &spin::QuadricMomentEvaluatorSettings::exactBinary64MaterialBoundary);
    module.def(
        "integrate_raw_mass_moments",
        [](
            const py::object& vertices,
            const py::object& faces,
            double density
        )
        {
            if (!std::isfinite(density) || density <= 0.0)
            {
                throw py::value_error("density must be finite and positive");
            }
            const MassMoments moments = computeSurfaceMassMoments(
                verticesFromPython(vertices),
                facesFromPython(faces),
                density
            );
            return Eigen::VectorXd(
                physics::RawMassMoments10::fromMassMoments(moments).values()
            );
        },
        py::arg("vertices"),
        py::arg("faces"),
        py::arg("density") = 1.0
    );

    py::enum_<spin::SpinMode>(module, "SpinMode")
        .value("STATIC_BALANCE", spin::SpinMode::StaticBalance)
        .value("BUOYANT_STABILITY", spin::SpinMode::BuoyantStability)
        .value("FIDGET_SPINNER", spin::SpinMode::FidgetSpinner)
        .value("YO_YO", spin::SpinMode::YoYo)
        .value("SPINNING_TOP", spin::SpinMode::SpinningTop);

    py::enum_<spin::StaticObjective>(module, "StaticObjective")
        .value("MAXIMUM_MASS", spin::StaticObjective::MaximumMass)
        .value("LOWER_CENTER_OF_MASS", spin::StaticObjective::LowerCenterOfMass)
        .value("PARETO_TRADEOFF", spin::StaticObjective::ParetoTradeoff);
    py::enum_<spin::StaticObjectiveMeasure>(module,"StaticObjectiveMeasure")
        .value("DENSITY_WEIGHTED_MASS_V1",spin::StaticObjectiveMeasure::DensityWeightedMass)
        .value("UNIFORM_DENSITY_VOLUME_V1",spin::StaticObjectiveMeasure::UniformDensityVolume);

    py::enum_<spin::SpinningTopObjective>(module,"SpinningTopObjective")
        .value("INERTIA_AND_LOWERING",spin::SpinningTopObjective::InertiaAndLowering)
        .value("LOWER_CENTER_OF_MASS",spin::SpinningTopObjective::LowerCenterOfMass);

    py::enum_<api::SpinItFasterEvidence>(module, "SpinItFasterEvidence")
        .value(
            "PAPER_FORMULA_REPRODUCTION",
            api::SpinItFasterEvidence::PaperFormulaReproduction
        )
        .value(
            "RECONSTRUCTED_ASSET",
            api::SpinItFasterEvidence::ReconstructedAsset
        );

    py::class_<spin::TriangleMesh>(module, "SpinTriangleMesh")
        .def(py::init<>())
        .def_readwrite("vertices", &spin::TriangleMesh::vertices)
        .def_readwrite("faces", &spin::TriangleMesh::faces);

    py::class_<spin::HydrostaticGeometry>(module, "HydrostaticGeometry")
        .def_property_readonly("submerged_boundary", [](const spin::HydrostaticGeometry& value) {
            return value.submergedBoundary;
        })
        .def_readonly("displaced_volume", &spin::HydrostaticGeometry::displacedVolume)
        .def_readonly("center_of_buoyancy", &spin::HydrostaticGeometry::centerOfBuoyancy)
        .def_readonly("waterline_area", &spin::HydrostaticGeometry::waterlineArea)
        .def_readonly("waterline_centroid", &spin::HydrostaticGeometry::waterlineCentroid)
        .def_readonly("waterline_central_second_moment",
                      &spin::HydrostaticGeometry::waterlineCentralSecondMoment);

    py::class_<spin::HydrostaticEquilibrium>(module, "HydrostaticEquilibrium")
        .def_readonly("geometry", &spin::HydrostaticEquilibrium::geometry)
        .def_readonly("waterline_signed_offset",
                      &spin::HydrostaticEquilibrium::waterlineSignedOffset)
        .def_readonly("target_displaced_volume",
                      &spin::HydrostaticEquilibrium::targetDisplacedVolume)
        .def_readonly("volume_residual", &spin::HydrostaticEquilibrium::volumeResidual)
        .def_readonly("iterations", &spin::HydrostaticEquilibrium::iterations);

    module.def(
        "compute_hydrostatic_geometry",
        [](const py::object& vertices, const py::object& faces,
           const Eigen::Vector3d& waterlineOrigin,
           const Eigen::Vector3d& waterlineNormal) {
            spin::TriangleMesh mesh;
            mesh.vertices = verticesFromPython(vertices);
            mesh.faces = facesFromPython(faces);
            return spin::computeHydrostaticGeometry(mesh, waterlineOrigin, waterlineNormal);
        },
        py::arg("vertices"), py::arg("faces"), py::arg("waterline_origin"),
        py::arg("waterline_normal")
    );

    module.def(
        "find_hydrostatic_equilibrium",
        [](const py::object& vertices, const py::object& faces,
           const Eigen::Vector3d& referenceOrigin,
           const Eigen::Vector3d& waterlineNormal,
           const double targetDisplacedVolume,
           const double relativeVolumeTolerance,
           const int maximumIterations) {
            spin::TriangleMesh mesh;
            mesh.vertices = verticesFromPython(vertices);
            mesh.faces = facesFromPython(faces);
            return spin::findHydrostaticEquilibrium(
                mesh, referenceOrigin, waterlineNormal, targetDisplacedVolume,
                relativeVolumeTolerance, maximumIterations
            );
        },
        py::arg("vertices"), py::arg("faces"), py::arg("reference_origin"),
        py::arg("waterline_normal"), py::arg("target_displaced_volume"),
        py::arg("relative_volume_tolerance") = 1e-10,
        py::arg("maximum_iterations") = 80
    );

    py::class_<spin::SpinTargetFrame>(module, "SpinTargetFrame")
        .def(py::init<>())
        .def_readwrite("origin", &spin::SpinTargetFrame::origin)
        .def_readwrite("axes", &spin::SpinTargetFrame::axes);

    py::enum_<spin::SpinEquationProfile>(module, "SpinEquationProfile")
        .value("LEGACY_NATIVE_NORMALIZED_V1", spin::SpinEquationProfile::LegacyNativeNormalizedV1)
        .value("DUAL_MONOMIAL_RAW_V1", spin::SpinEquationProfile::DualMonomialRawV1);

    py::class_<spin::SpinSolverSettings>(module, "SpinSolverSettings")
        .def(py::init<>())
        .def_readwrite("equation_profile", &spin::SpinSolverSettings::equationProfile)
        .def_readwrite("maximum_iterations", &spin::SpinSolverSettings::maximumIterations)
        .def_readwrite("maximum_restarts", &spin::SpinSolverSettings::maximumRestarts)
        .def_readwrite(
            "maximum_function_evaluations_per_restart",
            &spin::SpinSolverSettings::maximumFunctionEvaluationsPerRestart
        )
        .def_readwrite("equation_tolerance", &spin::SpinSolverSettings::equationTolerance)
        .def_readwrite(
            "physical_relative_tolerance",
            &spin::SpinSolverSettings::physicalRelativeTolerance
        )
        .def_readwrite("line_search_reduction", &spin::SpinSolverSettings::lineSearchReduction)
        .def_readwrite("maximum_geodesic_step", &spin::SpinSolverSettings::maximumGeodesicStep)
        .def_readwrite("deterministic_seed", &spin::SpinSolverSettings::deterministicSeed)
        .def_readwrite(
            "capture_evaluation_payload",
            &spin::SpinSolverSettings::captureEvaluationPayload
        );

    py::class_<spin::BuoyantStabilityData>(module, "BuoyantStabilityData")
        .def(py::init<>())
        .def_readwrite("liquid_density", &spin::BuoyantStabilityData::liquidDensity)
        .def_readwrite(
            "displaced_liquid_volume",
            &spin::BuoyantStabilityData::displacedLiquidVolume
        )
        .def_readwrite(
            "center_of_buoyancy_target",
            &spin::BuoyantStabilityData::centerOfBuoyancyTarget
        )
        .def_readwrite("waterline_area", &spin::BuoyantStabilityData::waterlineArea)
        .def_readwrite(
            "waterline_central_second_moment",
            &spin::BuoyantStabilityData::waterlineCentralSecondMoment
        )
        .def_readwrite(
            "fixed_shell_volume",
            &spin::BuoyantStabilityData::fixedShellVolume
        )
        .def_readwrite(
            "fixed_shell_density",
            &spin::BuoyantStabilityData::fixedShellDensity
        );

    py::class_<physics::RawMassMoments10>(module, "RawMassMoments10")
        .def(py::init<>())
        .def_static("from_values", &rawMomentsFromVector)
        .def_property_readonly(
            "values",
            [](const physics::RawMassMoments10& moments)
            {
                return Eigen::VectorXd(moments.values());
            }
        );

    py::class_<spin::SpinDesignRequest>(module, "SpinDesignRequest")
        .def(py::init<>())
        .def_readwrite("outer_surface", &spin::SpinDesignRequest::outerSurface)
        .def_readwrite("design_domain", &spin::SpinDesignRequest::designDomain)
        .def_readwrite("fixed_shell_moments", &spin::SpinDesignRequest::fixedShellMoments)
        .def_readwrite("target_frame", &spin::SpinDesignRequest::targetFrame)
        .def_readwrite("mode", &spin::SpinDesignRequest::mode)
        .def_readwrite("static_objective", &spin::SpinDesignRequest::staticObjective)
        .def_readwrite("static_objective_measure",&spin::SpinDesignRequest::staticObjectiveMeasure)
        .def_readwrite(
            "static_pareto_weight",
            &spin::SpinDesignRequest::staticParetoWeight
        )
        .def_readwrite("material_density", &spin::SpinDesignRequest::materialDensity)
        .def_readwrite(
            "original_units_per_normalized_unit",
            &spin::SpinDesignRequest::originalUnitsPerNormalizedUnit
        )
        .def_readwrite(
            "minimum_wall_thickness",
            &spin::SpinDesignRequest::minimumWallThickness
        )
        .def_readwrite("center_of_mass_weight", &spin::SpinDesignRequest::centerOfMassWeight)
        .def_readwrite("spinning_top_objective", &spin::SpinDesignRequest::spinningTopObjective)
        .def_readwrite("buoyancy", &spin::SpinDesignRequest::buoyancy)
        .def_readwrite("settings", &spin::SpinDesignRequest::settings);

    py::class_<spin::PlaneMomentEvaluatorSettings>(module, "PlaneMomentEvaluatorSettings")
        .def(py::init<>())
        .def_readwrite("truncation_radius", &spin::PlaneMomentEvaluatorSettings::truncationRadius);

    py::class_<spin::BuoyantPlaneBalanceSettings>(module, "BuoyantPlaneBalanceSettings")
        .def(py::init<>())
        .def_readwrite("maximum_iterations", &spin::BuoyantPlaneBalanceSettings::maximumIterations)
        .def_readwrite("maximum_line_search_steps", &spin::BuoyantPlaneBalanceSettings::maximumLineSearchSteps)
        .def_readwrite("constraint_tolerance", &spin::BuoyantPlaneBalanceSettings::constraintTolerance)
        .def_readwrite("backtracking_factor", &spin::BuoyantPlaneBalanceSettings::backtrackingFactor)
        .def_readwrite("minimum_step", &spin::BuoyantPlaneBalanceSettings::minimumStep)
        .def_readwrite("plane_evaluation", &spin::BuoyantPlaneBalanceSettings::planeEvaluation);

    module.def(
        "solve_planar_fidget_spinner",
        [](
            const std::vector<spin::PlanarBoundaryLoop>& designLoops,
            const std::vector<spin::PlanarBoundaryLoop>& fixedLoops,
            const std::optional<double>& fixedArealDensity,
            const double designArealDensity,
            const Eigen::Vector2d& initialCenter,
            const std::size_t maximumIterations,
            const std::size_t maximumFunctionEvaluations,
            const double residualTolerance,
            const double initialTrustRadius,
            const double maximumTrustRadius,
            const std::optional<Eigen::Vector2d>& extrusionInterval,
            const std::optional<double>& extrusionVolumeDensity,
            const py::list& fixedExtrusionLayers
        )
        {
            spin_bindings::PlanarEvaluationLedger ledger;
            try
            {
            if (extrusionInterval.has_value() != extrusionVolumeDensity.has_value())
                throw std::invalid_argument(
                    "extrusion interval and volume density must be supplied together");
            if (fixedExtrusionLayers.empty() && fixedLoops.empty())
                throw std::invalid_argument(
                    "fixed material must be supplied as loops or extrusion layers");
            if (fixedArealDensity &&
                (!std::isfinite(*fixedArealDensity) || *fixedArealDensity <= 0.0))
                throw std::invalid_argument(
                    "fixed areal density must be finite and positive");
            if (fixedLoops.empty() && fixedArealDensity)
                throw std::invalid_argument(
                    "omit fixed areal density when extrusion layers define the fixed material");
            if (!fixedLoops.empty() && !fixedArealDensity)
                throw std::invalid_argument(
                    "fixed areal density is required with fixed material loops");
            if (!fixedLoops.empty() &&
                spin::planarRegionsHavePositiveAreaIntersection(
                    designLoops, fixedLoops))
                throw std::invalid_argument(
                    "fixed and editable planar regions have positive-area overlap");

            struct FixedExtrusionLayer
            {
                std::vector<spin::PlanarBoundaryLoop> loops;
                Eigen::Vector2d z = Eigen::Vector2d::Zero();
                double density = 0.0;
            };
            std::vector<FixedExtrusionLayer> layers;
            Eigen::Matrix<double, 6, 1> projectedFixed =
                Eigen::Matrix<double, 6, 1>::Zero();
            std::vector<std::pair<Eigen::Vector2d,
                std::vector<spin::PlanarBoundaryLoop>>> priorLayerRegions;
            if (!fixedExtrusionLayers.empty())
            {
                if (!extrusionInterval)
                    throw std::invalid_argument(
                        "fixed extrusion layers require an extrusion interval for the design material");
                layers.reserve(fixedExtrusionLayers.size());
                for (const py::handle item : fixedExtrusionLayers)
                {
                    const py::dict layer = py::cast<py::dict>(item);
                    FixedExtrusionLayer parsed;
                    parsed.loops = py::cast<
                        std::vector<spin::PlanarBoundaryLoop>>(layer["loops"]);
                    parsed.z = py::cast<Eigen::Vector2d>(layer["z_interval"]);
                    parsed.density = py::cast<double>(layer["volume_density"]);
                    if (parsed.loops.empty() || !parsed.z.allFinite() ||
                        !(parsed.z(1) > parsed.z(0)) ||
                        !std::isfinite(parsed.density) || !(parsed.density > 0.0))
                        throw std::invalid_argument(
                            "fixed extrusion layer requires loops, finite increasing z, and positive density");
                    const bool overlapsDesignZ =
                        std::max(parsed.z(0), (*extrusionInterval)(0)) <
                        std::min(parsed.z(1), (*extrusionInterval)(1));
                    if (overlapsDesignZ &&
                        spin::planarRegionsHavePositiveAreaIntersection(
                            designLoops, parsed.loops))
                        throw std::invalid_argument(
                            "fixed extrusion layer overlaps the editable design material in 3D");
                    for (const auto& previous : priorLayerRegions)
                    {
                        const bool overlapsZ =
                            std::max(parsed.z(0), previous.first(0)) <
                            std::min(parsed.z(1), previous.first(1));
                        if (overlapsZ &&
                            spin::planarRegionsHavePositiveAreaIntersection(
                                parsed.loops, previous.second))
                            throw std::invalid_argument(
                                "fixed extrusion layers have positive-volume overlap");
                    }
                    const auto planar = spin::integratePlanarBoundaryLoops(parsed.loops);
                    projectedFixed += planar.weighted(
                        parsed.density*(parsed.z(1)-parsed.z(0))).vector();
                    priorLayerRegions.emplace_back(parsed.z, parsed.loops);
                    layers.push_back(std::move(parsed));
                }
                if (!fixedLoops.empty())
                {
                    Eigen::Matrix<double, 6, 1> suppliedFixed =
                        spin::integratePlanarBoundaryLoops(fixedLoops).weighted(
                            *fixedArealDensity).vector();
                    for (Eigen::Index i = 0; i < suppliedFixed.size(); ++i)
                        if (std::abs(projectedFixed(i)-suppliedFixed(i)) >
                            1e-10*std::max(1.0, std::abs(suppliedFixed(i))))
                            throw std::invalid_argument(
                                "fixed extrusion layer projections do not match fixed areal moments");
                }
            }
            spin::PlanarMassMoments6 fixedMoments;
            if (layers.empty())
            {
                const auto supplied = spin::integratePlanarBoundaryLoops(fixedLoops)
                    .weighted(*fixedArealDensity);
                fixedMoments = supplied;
            }
            else
            {
                fixedMoments.mass = projectedFixed(0);
                fixedMoments.firstX = projectedFixed(1);
                fixedMoments.firstY = projectedFixed(2);
                fixedMoments.secondXX = projectedFixed(3);
                fixedMoments.secondYY = projectedFixed(4);
                fixedMoments.secondXY = projectedFixed(5);
            }
            const spin::PlanarFidgetCircleProblem problem(
                designLoops, fixedMoments, designArealDensity);
            const spin::PlanarFidgetCircleEvaluation initial =
                ledger.assess(problem, initialCenter, "binding_initial_assessment");
            spin::PlanarDoglegOptions options;
            options.maximumIterations = maximumIterations;
            options.maximumFunctionEvaluations = maximumFunctionEvaluations;
            options.residualTolerance = residualTolerance;
            options.initialTrustRadius = initialTrustRadius;
            options.maximumTrustRadius = maximumTrustRadius;
            options.evaluationSink = [&](const spin::PlanarDoglegEvaluationEvent& event) {
                ledger.solver.push_back(event);
            };
            spin::PlanarDoglegResult solved;
            ledger.phase = "dogleg";
            {
                py::gil_scoped_release release;
                solved = spin::solvePlanarFidgetCircle(problem, initialCenter, options);
            }
            const spin::PlanarFidgetCircleEvaluation final =
                ledger.assess(problem, solved.center, "binding_final_assessment");
            ledger.phase = "boundary_and_extrusion";
            const auto boundary = spin::planarRegionDifferenceBoundary(
                designLoops, spin::PlanarCircleDisk{solved.center, solved.radius});
            py::dict output;
            output["status"] = spin::toString(solved.status);
            output["center"] = solved.center;
            output["radius"] = solved.radius;
            output["initial_normalized_residual"] = initial.normalizedResidual.norm();
            output["normalized_residual"] = solved.normalizedResidual;
            output["iterations"] = solved.iterations;
            output["function_evaluations"] = solved.functionEvaluations;
            output["accepted_steps"] = solved.acceptedSteps;
            output["final_trust_radius"] = solved.finalTrustRadius;
            output["fixed_material_moments"] = fixedMoments.vector();
            output["initial_material_moments"] = problem.initialMaterial().vector();
            output["removed_material_moments"] = final.removedMaterial.vector();
            output["remaining_material_moments"] = final.remainingMaterial.vector();
            output["normalized_first_moment_residual"] = final.normalizedResidual;
            output["fixed_editable_region_overlap_status"] =
                "DISJOINT_BY_EXACT_REGULARIZED_INTERSECTION";
            output["boundary_defined"] = boundary.defined;
            output["boundary_diagnostic"] = boundary.diagnostic;
            py::list cycles;
            for (const auto& cycle : boundary.cycles)
            {
                py::list pieces;
                for (const auto& piece : cycle)
                {
                    py::dict item;
                    item["type"] = piece.type == spin::PlanarBoundaryPieceType::LineSegment
                        ? "LINE_SEGMENT" : "CIRCULAR_ARC";
                    item["start"] = piece.start;
                    item["end"] = piece.end;
                    if (piece.type == spin::PlanarBoundaryPieceType::CircularArc)
                    {
                        item["center"] = piece.center;
                        item["radius"] = piece.radius;
                        item["start_angle"] = piece.startAngle;
                        item["sweep"] = piece.sweep;
                    }
                    pieces.append(std::move(item));
                }
                cycles.append(std::move(pieces));
            }
            output["boundary_cycles"] = std::move(cycles);
            output["geometry_scope"] = "exact_2d_line_arc_boundary";
            if (extrusionInterval)
            {
                if (!boundary.defined)
                    throw std::runtime_error(
                        "exact extrusion unavailable: " + boundary.diagnostic);
                const double designHeight = (*extrusionInterval)(1)-(*extrusionInterval)(0);
                const double effectiveDesignArealDensity =
                    *extrusionVolumeDensity*designHeight;
                if (std::abs(effectiveDesignArealDensity-designArealDensity) >
                    1e-10*std::max(1.0, std::abs(designArealDensity)))
                    throw std::invalid_argument(
                        "design areal density must equal extrusion volume density times height");
                const auto moments = spin::integratePlanarBoundaryCyclesExtruded(
                    boundary.cycles, (*extrusionInterval)(0), (*extrusionInterval)(1),
                    *extrusionVolumeDensity);
                output["extruded_design_material_raw_moments"] = moments.values();
                output["extrusion_interval"] = *extrusionInterval;
                output["extrusion_volume_density"] = *extrusionVolumeDensity;
                output["geometry_scope"] = "exact_2d_boundary_and_analytic_extrusion";

                if (layers.empty())
                    output["extruded_material_scope"] = "variable_design_region_only";
                else
                {
                    py::list fixedLayerMoments;
                    py::list fixedLayerGeometry;
                    auto completeMoments = moments;
                    for (const auto& layer : layers)
                    {
                        std::vector<spin::PlanarBoundaryCycle> cycles;
                        cycles.reserve(layer.loops.size());
                        for (const auto& loop : layer.loops)
                        {
                            if (loop.size() < 3)
                                throw std::invalid_argument(
                                    "fixed extrusion boundary loop requires at least three points");
                            spin::PlanarBoundaryCycle cycle;
                            cycle.reserve(loop.size());
                            for (std::size_t i = 0; i < loop.size(); ++i)
                            {
                                spin::PlanarBoundaryPiece piece;
                                piece.type = spin::PlanarBoundaryPieceType::LineSegment;
                                piece.start = loop[i];
                                piece.end = loop[(i+1)%loop.size()];
                                cycle.push_back(std::move(piece));
                            }
                            cycles.push_back(std::move(cycle));
                        }
                        const auto layerMoments =
                            spin::integratePlanarBoundaryCyclesExtruded(
                                cycles, layer.z(0), layer.z(1), layer.density);
                        fixedLayerMoments.append(layerMoments.values());
                        completeMoments += layerMoments;

                        py::dict geometry;
                        geometry["loops"] = layer.loops;
                        geometry["z_interval"] = layer.z;
                        geometry["volume_density"] = layer.density;
                        fixedLayerGeometry.append(std::move(geometry));
                    }
                    output["extruded_fixed_layer_raw_moments"] = std::move(fixedLayerMoments);
                    output["fixed_extrusion_layers"] = std::move(fixedLayerGeometry);
                    output["extruded_complete_material_raw_moments"] =
                        completeMoments.values();
                    output["extruded_material_scope"] = "design_plus_fixed_layers";
                    output["geometry_scope"] =
                        "exact_2d_boundary_and_analytic_layered_extrusion";
                }
            }
            ledger.phase = "complete";
            output["evaluation_ledger"] = ledger.toPython();
            return output;
            }
            catch (...)
            {
                ledger.rethrowWithLedger();
            }
        },
        py::arg("design_region_loops"), py::arg("fixed_material_loops"),
        py::arg("fixed_areal_density"), py::arg("design_areal_density"),
        py::arg("initial_center"), py::arg("maximum_iterations") = 100,
        py::arg("maximum_function_evaluations") = 60,
        py::arg("residual_tolerance") = 1e-10,
        py::arg("initial_trust_radius") = 0.25,
        py::arg("maximum_trust_radius") = 1.0,
        py::arg("extrusion_interval") = std::nullopt,
        py::arg("extrusion_volume_density") = std::nullopt,
        py::arg("fixed_extrusion_layers") = py::list());

    module.def(
        "evaluate_spin_it_faster_plane_halfspace",
        [](const spin::SpinDesignRequest& request,
           const Eigen::Vector4d& coefficients) {
            if (!coefficients.allFinite() || coefficients.norm() <= 1e-15)
                throw py::value_error("plane coefficients must be finite and nonzero");
            spin::PlaneMomentEvaluation value;
            {
                py::gil_scoped_release release;
                value = spin::evaluatePlaneMaterialMoments(request, coefficients);
            }
            py::dict output;
            output["raw_moments"] = Eigen::VectorXd(value.variableMomentsTarget.values());
            output["vertices"] = value.variableMaterialBoundaryTarget.vertices;
            output["faces"] = value.variableMaterialBoundaryTarget.faces;
            output["geometry_hash"] = value.geometryHash;
            output["cap_separation_certificate"] = value.capSeparationCertificate;
            output["interface_triangle_count"] = value.interfaceTriangleCount;
            return output;
        },
        py::arg("request"), py::arg("coefficients"),
        "Evaluate the actual finite Boolean material for one oriented plane half-space. "
        "Call with c and -c to obtain the complementary design-domain pieces; "
        "the caller must still authenticate their exact partition.");

    module.def(
        "solve_buoyant_stability_by_plane",
        [](const spin::SpinDesignRequest& request,
           const Eigen::Vector3d& initialParameters,
           const spin::BuoyantPlaneBalanceSettings& settings) {
            spin::BuoyantPlaneBalanceResult value;
            {
                py::gil_scoped_release release;
                value = spin::solveBuoyantStabilityByPlane(request, initialParameters, settings);
            }
            py::dict output;
            output["status"] = spin::buoyantPlaneBalanceStatusToken(value.status);
            output["diagnostic"] = value.diagnostic;
            output["plane_coefficients"] = value.planeCoefficients;
            output["total_raw_moments"] = Eigen::VectorXd(value.totalMoments.values());
            output["variable_raw_moments"] = Eigen::VectorXd(value.variableMoments.values());
            output["constraints"] = value.constraints;
            output["normalized_constraints"] = value.normalizedConstraints;
            output["normalized_constraints_norm"] = value.normalizedConstraints.norm();
            output["target_total_mass"] = value.targetTotalMass;
            output["center_of_mass_height"] = value.centerOfMassHeight;
            output["metacenter_height"] = value.metacenterHeight;
            output["stability_margin"] = value.stabilityMargin;
            output["equilibrium_constraints_verified"] = value.equilibriumConstraintsVerified;
            output["eq7_numerical_conditions_verified"] = value.eq7NumericalConditionsVerified;
            output["stable_at_prescribed_waterline"] = value.stableAtPrescribedWaterline;
            output["iterations"] = value.iterations;
            output["accepted_steps"] = value.acceptedSteps;
            output["variable_material_vertices"] = value.evaluation.variableMaterialBoundaryTarget.vertices;
            output["variable_material_faces"] = value.evaluation.variableMaterialBoundaryTarget.faces;
            py::list trace;
            for (const auto& row : value.trace) {
                py::dict event;
                event["evaluation"] = row.evaluation;
                event["iteration"] = row.iteration;
                event["plane_coefficients"] = row.planeCoefficients;
                event["constraints"] = row.constraints;
                event["normalized_constraints"] = row.normalizedConstraints;
                event["residual_norm"] = row.residualNorm;
                event["step"] = row.step;
                event["accepted"] = row.accepted;
                event["evaluation_succeeded"] = row.evaluationSucceeded;
                event["derivative_available"] = row.derivativeAvailable;
                event["diagnostic"] = row.diagnostic;
                trace.append(std::move(event));
            }
            output["trace"] = std::move(trace);
            return output;
        },
        py::arg("request"), py::arg("initial_parameters"),
        py::arg("settings") = spin::BuoyantPlaneBalanceSettings());

    module.def(
        "solve_static_balance_by_plane_family",
        [](const spin::SpinDesignRequest& request,
           const Eigen::Vector4d& initialCoefficients,
           const int maximumIterations,
           const int maximumLineSearchSteps,
           const double constraintTolerance) {
            spin::PlaneStaticBalanceSettings settings;
            settings.maximumIterations = maximumIterations;
            settings.maximumLineSearchSteps = maximumLineSearchSteps;
            settings.constraintTolerance = constraintTolerance;
            spin::PlaneStaticBalanceResult value;
            {
                py::gil_scoped_release release;
                value = spin::solveStaticBalanceByPlaneFamily(
                    request, initialCoefficients, settings);
            }
            return staticPlaneResultToPython(value, request);
        },
        py::arg("request"), py::arg("initial_coefficients"),
        py::arg("maximum_iterations") = 60,
        py::arg("maximum_line_search_steps") = 24,
        py::arg("constraint_tolerance") = 1e-10);

    py::class_<api::SpinItFasterArtifact>(module, "SpinItFasterArtifact")
        .def_property_readonly("status", &api::SpinItFasterArtifact::statusToken)
        .def_property_readonly("evidence", &api::SpinItFasterArtifact::evidenceToken)
        .def_property_readonly("succeeded", &api::SpinItFasterArtifact::succeeded)
        .def_property_readonly("geometry_sealed", &api::SpinItFasterArtifact::geometrySealed)
        .def_property_readonly("derivative_method", [](const api::SpinItFasterArtifact& v) { return v.result().derivativeMethod; })
        .def_property_readonly("derivative_diagnostic", [](const api::SpinItFasterArtifact& v) { return v.result().derivativeDiagnostic; })
        .def_property_readonly("derivative_scientific_hash", [](const api::SpinItFasterArtifact& v) { return v.result().derivativeScientificHash; })
        .def_property_readonly(
            "residual_evaluation_seconds",
            &api::SpinItFasterArtifact::residualEvaluationSeconds
        )
        .def_property_readonly(
            "derivative_evaluation_seconds",
            &api::SpinItFasterArtifact::derivativeEvaluationSeconds
        )
        .def_property_readonly(
            "solve_evaluator_wall_seconds",
            &api::SpinItFasterArtifact::solveEvaluatorWallSeconds
        )
        .def_property_readonly(
            "solve_function_evaluations",
            &api::SpinItFasterArtifact::solveFunctionEvaluations
        )
        .def_property_readonly(
            "post_solve_material_evaluation_count",
            &api::SpinItFasterArtifact::postSolveMaterialEvaluationCount
        )
        .def_property_readonly(
            "post_solve_material_evaluation_wall_seconds",
            &api::SpinItFasterArtifact::postSolveMaterialEvaluationWallSeconds
        )
        .def_property_readonly(
            "event_hash",
            &api::SpinItFasterArtifact::eventHash
        )
        .def_property_readonly("equation_profile", [](const api::SpinItFasterArtifact& value) {
            return spin::spinEquationProfileToken(value.result().equationProfile);
        })
        .def_property_readonly("solver_coordinate_identity", [](const api::SpinItFasterArtifact& value) {
            return spin::spinSolverCoordinateIdentity(value.result().equationProfile);
        })
        .def_property_readonly("raw_equation_residual", [](const api::SpinItFasterArtifact& value) -> py::object {
            const auto residual=value.result().rawEquationResidual;
            return residual ? py::cast(*residual) : py::none();
        })
        .def_property_readonly(
            "evaluation_timings",
            [](const api::SpinItFasterArtifact& value)
            {
                py::list timings;
                for (const spin::SpinEvaluationTiming& timing :
                     value.result().evaluationTimings)
                {
                    py::dict entry;
                    entry["evaluation_index"] = timing.evaluationIndex;
                    entry["restart"] = timing.restart;
                    entry["restart_evaluation_index"] =
                        timing.restartEvaluationIndex;
                    entry["evaluator_returned"] = timing.evaluatorReturned;
                    entry["wall_seconds"] = timing.wallSeconds;
                    entry["reported_residual_seconds"] =
                        timing.reportedResidualSeconds
                        ? py::cast(*timing.reportedResidualSeconds)
                        : py::none();
                    entry["reported_derivative_seconds"] =
                        timing.reportedDerivativeSeconds
                        ? py::cast(*timing.reportedDerivativeSeconds)
                        : py::none();
                    timings.append(std::move(entry));
                }
                return timings;
            }
        )
        .def_property_readonly(
            "outer_vertices",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXd(value.outerVertices());
            }
        )
        .def_property_readonly(
            "outer_faces",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXi(value.outerFaces());
            }
        )
        .def_property_readonly(
            "design_vertices",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXd(value.designVertices());
            }
        )
        .def_property_readonly(
            "design_faces",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXi(value.designFaces());
            }
        )
        .def_property_readonly(
            "material_vertices",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXd(value.materialVertices());
            }
        )
        .def_property_readonly(
            "material_faces",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXi(value.materialFaces());
            }
        )
        .def_property_readonly(
            "cavity_vertices",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXd(value.cavityVertices());
            }
        )
        .def_property_readonly(
            "cavity_faces",
            [](const api::SpinItFasterArtifact& value)
            {
                return Eigen::MatrixXi(value.cavityFaces());
            }
        )
        .def_property_readonly(
            "diagnostics_json",
            [](const api::SpinItFasterArtifact& value)
            {
                return std::string(value.diagnosticsJson());
            }
        )
        .def("save_outer_mesh", &api::SpinItFasterArtifact::saveOuterMesh)
        .def("save_design_mesh", &api::SpinItFasterArtifact::saveDesignMesh)
        .def("save_material_mesh", &api::SpinItFasterArtifact::saveMaterialMesh)
        .def("save_cavity_mesh", &api::SpinItFasterArtifact::saveCavityMesh)
        .def("save_diagnostics", &api::SpinItFasterArtifact::saveDiagnostics);

    module.def(
        "solve_spin_it_faster",
        [](
            const spin::SpinDesignRequest& request,
            api::SpinItFasterEvidence evidence,
            const Eigen::MatrixXd& initialQuadrics,
            const spin::QuadricMomentEvaluatorSettings& evaluatorSettings
        )
        {
            const std::vector<spin::QuadricCoefficients10> initial =
                quadricsFromRows(initialQuadrics);
            py::gil_scoped_release release;
            const api::SpinItFasterArtifact artifact =
                api::solveSpinItFaster(request, evidence, evaluatorSettings, initial);
            return std::make_unique<api::SpinItFasterArtifact>(artifact);
        },
        py::arg("request"),
        py::arg("evidence"),
        py::arg("initial_quadrics") = Eigen::MatrixXd(0, 10),
        py::arg("evaluator_settings") = spin::QuadricMomentEvaluatorSettings()
    );

    module.def(
        "evaluate_spin_it_faster_quadric",
        [](
            const spin::SpinDesignRequest& request,
            const Eigen::VectorXd& coefficients,
            const std::optional<spin::QuadricMomentEvaluatorSettings>& explicitSettings,
            bool coefficientsAreNormalized
        )
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error(
                    "quadric coefficients must contain 10 finite values"
                );
            }
            const auto quadric = coefficientsAreNormalized
                ? spin::QuadricCoefficients10::fromNormalized(coefficients)
                : spin::QuadricCoefficients10::fromUnnormalized(coefficients);
            const spin::QuadricMomentEvaluation evaluation =
                spin::evaluateQuadricMaterialMoments(
                    request,
                    quadric,
                    explicitSettings.value_or(spin::QuadricMomentEvaluatorSettings())
                );
            py::dict result;
            result["variable_raw_moments"] =
                Eigen::VectorXd(evaluation.variableMomentsTarget.values());
            result["total_raw_moments"] =
                Eigen::VectorXd(evaluation.totalMomentsTarget.values());
            result["geometry_hash"] = evaluation.geometryHash;
            result["topology_hash"] = evaluation.topologyHash;
            if (explicitSettings.has_value())
            {
                result["boolean_backend_identity"] = evaluation.booleanBackendIdentity;
                result["base_request_hash"] = spin::computeSpinDesignRequestHash(request);
                result["request_hash"] = spin::computeQuadricEvaluatorRequestHash(request, *explicitSettings);
                result["derivative_method"] = spin::quadricMomentDerivativeMethodToken(evaluation.derivativeMethod);
                result["derivative_defined"] = evaluation.derivativeDefined;
                result["derivative_diagnostic"] = evaluation.derivativeDiagnostic;
                result["stratum_derivative_defined"] =
                    evaluation.stratumDifferential.defined;
                result["stratum_derivative_diagnostic"] =
                    evaluation.stratumDifferential.diagnostic;
                result["stratum_moment_jacobian"] =
                    evaluation.stratumDifferential.jacobian;
                result["stratum_coefficient_basis"] =
                    evaluation.stratumDifferential.coefficientBasis;
                result["derivative_scientific_hash"] = evaluation.derivativeScientificHash;
                result["provenance_topology_hash"] = evaluation.provenanceTopologyHash;
                // This sign is the native Boolean-operation decision made
                // from the closed finite quadric solid, not a classification
                // inferred from a reference Boolean result.  It lets an
                // independent backend receive the same operands *and* the
                // same intersection/minus semantics for diagnostic checks.
                result["quadric_solid_is_positive_region"] =
                    evaluation.quadricSolidIsPositiveRegion;
                result["vertex_provenance"] = evaluation.vertexProvenance;
                // Expose the exact finite triangulation passed to the native
                // Boolean.  This is diagnostic-only and opt-in so an
                // independent backend can receive byte-identical operands;
                // it neither reruns nor replaces the CGAL operation.
                result["interface_vertices"] =
                    evaluation.quadricTriangulation.closedBoundary.vertices;
                result["interface_faces"] =
                    evaluation.quadricTriangulation.closedBoundary.faces;
                result["material_vertices"] = evaluation.materialBoundary.vertices;
                result["material_faces"] = evaluation.materialBoundary.faces;
                result["outer_vertices"] = request.outerSurface.vertices;
                result["outer_faces"] = request.outerSurface.faces;
                if(explicitSettings->triangulation.rationalSourceCoordinates)
                    result["rational_authority"]=rationalAuthorityToPython(
                        evaluation.materialExactVerticesTarget,evaluation.designExactVerticesTarget,
                        evaluation.quadricTriangulation.exactVertices,
                        evaluation.quadricTriangulation.coordinateLower,
                        evaluation.quadricTriangulation.coordinateUpper,
                        evaluation.quadricTriangulation.closedBoundary.faces,
                        evaluation.materialBirthFaces,evaluation.variableMomentsExactTarget,
                        evaluation.geometryHash,evaluation.quadricTriangulation.rationalConstructionIdentity);
            }
            result["interface_type"] = std::string(
                spin::quadricSurfaceTypeToken(
                    evaluation.quadricTriangulation.type
                )
            );
            result["interface_subtype"] = std::string(
                spin::quadricSurfaceSubtypeToken(
                    spin::classifyQuadricSurfaceSubtype(
                        quadric
                    )
                )
            );
            return result;
        },
        py::arg("request"),
        py::arg("quadric_coefficients"),
        py::arg("evaluator_settings") = py::none(),
        py::arg("coefficients_are_normalized") = false
    );

    // Opt-in diagnostic export of the exact default-evaluator result used by
    // the solver. Keep the established moments-only API and solver defaults
    // unchanged. An undefined ambient derivative is represented by None, not
    // by the evaluator's zero-initialized storage.
    module.def(
        "evaluate_spin_it_faster_quadric_diagnostic",
        [](
            const spin::SpinDesignRequest& request,
            const Eigen::VectorXd& coefficients
        )
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error(
                    "quadric coefficients must contain 10 finite values"
                );
            }
            const spin::QuadricMomentEvaluatorSettings settings{};
            const spin::QuadricMomentEvaluation evaluation =
                spin::evaluateQuadricMaterialMoments(
                    request,
                    spin::QuadricCoefficients10::fromUnnormalized(coefficients),
                    settings
                );
            py::dict result;
            result["variable_raw_moments"] =
                Eigen::VectorXd(evaluation.variableMomentsTarget.values());
            result["total_raw_moments"] =
                Eigen::VectorXd(evaluation.totalMomentsTarget.values());
            py::object momentJacobian = py::none();
            if (evaluation.derivativeDefined)
                momentJacobian = py::cast(evaluation.momentJacobian);
            result["moment_jacobian"] = std::move(momentJacobian);
            result["derivative_defined"] = evaluation.derivativeDefined;
            result["derivative_method"] =
                spin::quadricMomentDerivativeMethodToken(evaluation.derivativeMethod);
            result["derivative_diagnostic"] = evaluation.derivativeDiagnostic;
            result["residual_evaluation_seconds"] = evaluation.residualEvaluationSeconds;
            result["derivative_evaluation_seconds"] = evaluation.derivativeEvaluationSeconds;
            result["geometry_hash"] = evaluation.geometryHash;
            result["topology_hash"] = evaluation.topologyHash;
            result["base_request_hash"] = spin::computeSpinDesignRequestHash(request);
            result["evaluator_request_hash"] =
                spin::computeQuadricEvaluatorRequestHash(request, settings);
            result["derivative_scientific_hash"] = evaluation.derivativeScientificHash;
            result["provenance_topology_hash"] = evaluation.provenanceTopologyHash;
            result["equation_profile"] =
                spin::spinEquationProfileToken(request.settings.equationProfile);
            return result;
        },
        py::arg("request"), py::arg("quadric_coefficients")
    );

    module.def(
        "evaluate_spin_it_faster_fullrank_quadric",
        [](
            const spin::SpinDesignRequest& request,
            const Eigen::VectorXd& coefficients,
            const std::optional<spin::QuadricMomentEvaluatorSettings>& explicitSettings,
            bool captureMaterialProof
        )
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error(
                    "quadric coefficients must contain 10 finite values"
                );
            }
            spin::FullRankMomentEvaluation evaluation;
            const auto normalized=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
            {
                py::gil_scoped_release release;
                evaluation = spin::evaluateCanonicalFullRankMaterialMoments(
                    request,
                    normalized,
                    explicitSettings.value_or(spin::QuadricMomentEvaluatorSettings())
                );
            }
            py::dict result;
            result["variable_raw_moments"] =
                Eigen::VectorXd(evaluation.variableMomentsTarget.values());
            result["total_raw_moments"] =
                Eigen::VectorXd(evaluation.totalMomentsTarget.values());
            result["geometry_hash"] = evaluation.geometryHash;
            result["boolean_backend_identity"] = evaluation.booleanBackendIdentity;
            result["source_scientific_hash"] = evaluation.sourceScientificHash;
            result["scientific_hash"] = evaluation.scientificHash;
            result["derivative_defined"] = evaluation.coefficientJacobian.has_value();
            result["derivative_diagnostic"] = evaluation.diagnostic;
            if (evaluation.coefficientJacobian)
            {
                result["coefficient_jacobian"] = Eigen::MatrixXd(*evaluation.coefficientJacobian);
            }
            else
            {
                result["coefficient_jacobian"] = py::none();
            }
            result["material_vertices"] = Eigen::MatrixXd(evaluation.materialBoundaryTarget.vertices);
            result["material_faces"] = Eigen::MatrixXi(evaluation.materialBoundaryTarget.faces);
            if(captureMaterialProof) {
                // Attach the already-computed SAME Boolean tapes. No replay,
                // new geometry, derivative or scientific hash is introduced.
                py::dict proof;
                proof["schema"]="same-call-cgal-material-proof-v1";
                proof["exact_material_vertices"]=evaluation.materialExactVerticesTarget;
                proof["material_birth_faces"]=Eigen::VectorXi(evaluation.materialBirthFaces);
                proof["source_faces"]=Eigen::MatrixXi(evaluation.sourceFaces);
                proof["provenance_topology_hash"]=evaluation.topologyHash;
                proof["geometry_hash"]=evaluation.geometryHash;
                result["material_proof_attachment"]=std::move(proof);
            }
            if(!evaluation.sourceExactVerticesTarget.empty()) {
                result["rational_authority"]=rationalAuthorityToPython(evaluation);
                result["normalized_coefficients"]=Eigen::VectorXd(normalized.values());
                result["branch_identity"]=evaluation.rationalBranchIdentity;
            }
            return result;
        },
        py::arg("request"),
        py::arg("coefficients"),
        py::arg("evaluator_settings") = py::none(),
        py::arg("capture_material_proof") = false
    );

    module.def(
        "capture_spin_it_faster_fullrank_geometry",
        [](const spin::SpinDesignRequest& request,
           const Eigen::VectorXd& coefficients,
           const spin::QuadricMomentEvaluatorSettings& settings)
        {
            if (coefficients.size()!=10 || !coefficients.allFinite())
                throw py::value_error("quadric coefficients must contain 10 finite values");
            const auto normalized=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
            std::vector<std::pair<std::string,spin::QuadricMomentEvaluation>> snapshots;
            std::string error;
            bool completed=false;
            {
                py::gil_scoped_release release;
                spin::ScopedQuadricGeometryObserver observer(
                    [&](const std::string& stage,const spin::QuadricMomentEvaluation& value) {
                        snapshots.emplace_back(stage,value);
                    });
                try {
                    (void)spin::evaluateCanonicalFullRankMaterialMoments(request,normalized,settings);
                    completed=true;
                } catch (const std::exception& failure) {
                    error=failure.what();
                }
            }
            py::list tape;
            for (const auto& item:snapshots) {
                const auto& value=item.second;
                const auto& source=value.quadricTriangulation;
                py::dict entry;
                entry["stage"]=item.first;
                entry["source_vertices"]=source.closedBoundary.vertices;
                entry["source_faces"]=source.closedBoundary.faces;
                entry["source_exact_vertices"]=source.exactVertices;
                entry["source_face_is_cap"]=source.faceIsCap;
                entry["source_fixed_normals"]=source.fixedNormals;
                entry["source_normal_edges"]=source.sourceNormalEdgeByVertex;
                entry["source_normal_vertices"]=source.sourceNormalVertexByVertex;
                entry["source_projective_rays"]=source.sourceProjectiveRayByVertex;
                entry["source_hash"]=source.hash;
                entry["truncation_radius"]=source.truncationRadius;
                entry["boundary_loop_count"]=source.boundaryLoopCount;
                entry["material_vertices"]=value.materialBoundaryTarget.vertices;
                entry["material_faces"]=value.materialBoundaryTarget.faces;
                entry["exact_material_vertices"]=value.materialExactVerticesTarget;
                entry["material_birth_faces"]=value.materialBirthFaces;
                // Before side selection this bool has no semantic value.
                entry["selected_positive_region"]=item.first=="source-before-side"
                    ? py::object(py::none()) : py::cast(value.quadricSolidIsPositiveRegion);
                tape.append(std::move(entry));
            }
            py::dict result;
            result["schema"]="same-production-call-failure-tape-v1";
            result["production_admission"]=false;
            result["evaluation_completed"]=completed;
            result["error"]=error;
            result["normalized_coefficients"]=Eigen::VectorXd(normalized.values());
            result["snapshots"]=std::move(tape);
            return result;
        },
        py::arg("request"),py::arg("coefficients"),py::arg("evaluator_settings")
    );

    module.def(
        "solve_spin_it_faster_fullrank_quadric",
        [](
            const spin::SpinDesignRequest& request,
            const Eigen::MatrixXd& initialQuadrics,
            const std::optional<spin::QuadricMomentEvaluatorSettings>& explicitSettings
        )
        {
            if (initialQuadrics.cols() != 10 || !initialQuadrics.allFinite())
            {
                throw py::value_error(
                    "initial quadrics must be a finite N-by-10 matrix"
                );
            }
            std::vector<spin::QuadricCoefficients10> initial;
            initial.reserve(static_cast<std::size_t>(initialQuadrics.rows()));
            for (Eigen::Index row = 0; row < initialQuadrics.rows(); ++row)
            {
                initial.push_back(spin::QuadricCoefficients10::fromUnnormalized(
                    initialQuadrics.row(row).transpose()
                ));
            }

            spin::FullRankHollowingResult solved;
            {
                py::gil_scoped_release release;
                solved = spin::solveCanonicalFullRankQuadricHollowing(
                    request,
                    explicitSettings.value_or(spin::QuadricMomentEvaluatorSettings()),
                    initial
                );
            }

            return fullRankResultToPython(solved, request.settings.captureEvaluationPayload);
        },
        py::arg("request"),
        py::arg("initial_quadrics") = Eigen::MatrixXd(0, 10),
        py::arg("evaluator_settings") = py::none()
    );

    module.def(
        "solve_spin_it_faster_external_moment_evaluator",
        [](
            const spin::SpinDesignRequest& request,
            const py::function& pythonEvaluator,
            const Eigen::MatrixXd& initialQuadrics,
            const py::object& pythonEventSink,
            const py::object& pythonSupportBasis,
            const py::object& pythonObjectiveSettings
        )
        {
            if (request.settings.equationProfile !=
                spin::SpinEquationProfile::DualMonomialRawV1)
            {
                throw py::value_error(
                    "external strict evaluator requires DUAL_MONOMIAL_RAW_V1"
                );
            }
            // External maps bypass the normal request-owning solver wrapper.
            spin::validateSpinningTopObjective(request);
#if MAKE_IT_STAND_HAS_SPIN_IT_2014
            std::optional<spin2014::QuadricObjectiveSettings> originalObjective;
            if (!pythonObjectiveSettings.is_none())
            {
                originalObjective = pythonObjectiveSettings.cast<spin2014::QuadricObjectiveSettings>();
                spin2014::validateSpin2014QuadricObjectiveSettings(*originalObjective);
                if (spin2014::spin2014ObjectiveMode(originalObjective->objective) != request.mode)
                    throw py::value_error("original objective and request mode disagree");
                if (request.mode == spin::SpinMode::StaticBalance || !pythonSupportBasis.is_none())
                    throw py::value_error("original objective requires rotational full-space entry; static uses its plane solver");
            }
#else
            if (!pythonObjectiveSettings.is_none())
                throw py::value_error("native build lacks Spin-It 2014 objective support");
#endif
            // Never interpret arbitrary fixed mass as homogeneous volume.
            if(request.staticObjectiveMeasure!=spin::StaticObjectiveMeasure::DensityWeightedMass)
                spin::validateSpinDesignRequest(request);
            if (!pythonEvaluator || initialQuadrics.cols() != 10 ||
                !initialQuadrics.allFinite())
            {
                throw py::value_error(
                    "external evaluator and finite N-by-10 initial quadrics are required"
                );
            }
            if (!pythonEventSink.is_none() && !PyCallable_Check(pythonEventSink.ptr()))
                throw py::value_error("event_sink must be callable or None");
            const bool directional = !pythonSupportBasis.is_none();
            const Eigen::MatrixXd supportBasis = directional
                ? py::cast<Eigen::MatrixXd>(pythonSupportBasis) : Eigen::MatrixXd();
            if (directional && initialQuadrics.rows()==0)
                throw py::value_error("supported-space evaluator requires explicit initial quadrics");

            spin::SpinMomentProblem problem = [&]() {
                switch (request.mode)
                {
                case spin::SpinMode::StaticBalance:
                    return request.staticObjective == spin::StaticObjective::ParetoTradeoff
                        ? spin::SpinMomentProblem::staticBalancePareto(request.staticParetoWeight,spin::staticObjectiveMassScale(request))
                        : spin::SpinMomentProblem::staticBalance(request.staticObjective,spin::staticObjectiveMassScale(request));
                case spin::SpinMode::BuoyantStability:
                    return spin::SpinMomentProblem::buoyantStability(
                        request.buoyancy.liquidDensity * request.buoyancy.displacedLiquidVolume,
                        request.buoyancy.centerOfBuoyancyTarget.head<2>()
                    );
                case spin::SpinMode::FidgetSpinner:
                    return spin::SpinMomentProblem::fidgetSpinner();
                case spin::SpinMode::YoYo:
                    return spin::SpinMomentProblem::yoYo();
                case spin::SpinMode::SpinningTop:
                    return spin::SpinMomentProblem::spinningTop(request.centerOfMassWeight,request.spinningTopObjective);
                }
                throw py::value_error("unsupported Spin-It Faster mode");
            }();

            std::vector<spin::QuadricCoefficients10> guesses;
            guesses.reserve(static_cast<std::size_t>(initialQuadrics.rows()));
            for (Eigen::Index row = 0; row < initialQuadrics.rows(); ++row)
            {
                guesses.push_back(spin::QuadricCoefficients10::fromUnnormalized(
                    initialQuadrics.row(row).transpose()
                ));
            }
            if (guesses.empty())
            {
                guesses = spin::deterministicInitialQuadrics(
                    request.mode, request.settings
                );
            }

            int callbackIndex = 0;
            py::list pythonCallbackDiagnostics;
            const auto mapEvaluator =
                [pythonEvaluator,directional,&callbackIndex,&pythonCallbackDiagnostics](const spin::QuadricCoefficients10& quadric)
            {
                ++callbackIndex;
                try
                {
                py::tuple coefficients(10);
                for (Eigen::Index i = 0; i < 10; ++i)
                {
                    coefficients[i] = quadric.values()(i);
                }
                const py::object returned = pythonEvaluator(coefficients);
                if (!py::isinstance<py::dict>(returned))
                {
                    throw std::runtime_error(
                        "external evaluator must return a mapping/dict"
                    );
                }
                const py::dict value = py::cast<py::dict>(returned);
                const auto required = [&value](const char* key) -> py::handle {
                    if (!value.contains(key))
                    {
                        throw std::runtime_error(
                            std::string("external evaluator omitted '") + key + "'"
                        );
                    }
                    return value[key];
                };
                const std::string geometryHash = py::cast<std::string>(
                    required("geometry_hash")
                );
                const std::string topologyHash = py::cast<std::string>(
                    required("topology_hash")
                );
                const std::string evaluationMapIdentity = py::cast<std::string>(
                    required("evaluation_map_identity")
                );
                const std::string booleanBackendIdentity = py::cast<std::string>(
                    required("boolean_backend_identity")
                );
                const std::string derivativeIdentity = py::cast<std::string>(
                    required("derivative_identity")
                );
                if (geometryHash.empty() || topologyHash.empty())
                {
                    throw std::runtime_error(
                        "external evaluator geometry/topology hashes must be nonempty"
                    );
                }
                const auto validIdentity = [](const std::string& identity) {
                    return !identity.empty() && identity.find(';') == std::string::npos;
                };
                if (!validIdentity(evaluationMapIdentity) ||
                    !validIdentity(booleanBackendIdentity) ||
                    !validIdentity(derivativeIdentity))
                {
                    throw std::runtime_error(
                        "external evaluator map/backend/derivative identities must be nonempty and delimiter-safe"
                    );
                }
                spin::SpinMomentSubspaceMapEvaluation parsed;
                auto& map = parsed.material;
                const Eigen::VectorXd moments = py::cast<Eigen::VectorXd>(
                    required("raw_moments")
                );
                if (moments.size() != 10)
                {
                    throw std::runtime_error("external evaluator must return a 10-vector");
                }
                map.rawMoments = moments;
                if (directional)
                {
                    if (value.contains("moment_jacobian"))
                        throw std::runtime_error("directional evaluator must not supply an ambient moment_jacobian");
                    parsed.nativeTangentBasis = py::cast<Eigen::MatrixXd>(required("native_tangent_basis"));
                    parsed.nativeDirectionalJacobian = py::cast<Eigen::MatrixXd>(required("native_directional_jacobian"));
                }
                else
                {
                    const Eigen::MatrixXd jacobian = py::cast<Eigen::MatrixXd>(required("moment_jacobian"));
                    if (jacobian.rows()!=10 || jacobian.cols()!=10)
                        throw std::runtime_error("external evaluator must return a 10-vector and 10-by-10 Jacobian");
                    map.momentJacobian = jacobian;
                }
                if (value.contains("variable_raw_moments") &&
                    !value["variable_raw_moments"].is_none())
                {
                    const Eigen::VectorXd variableMoments =
                        py::cast<Eigen::VectorXd>(value["variable_raw_moments"]);
                    if (variableMoments.size() != 10 || !variableMoments.allFinite())
                    {
                        throw std::runtime_error(
                            "external evaluator variable_raw_moments must contain 10 finite values"
                        );
                    }
                    map.variableRawMoments = variableMoments;
                }
                map.derivativeDefined = py::cast<bool>(
                    required("derivative_defined")
                );
                if (value.contains("derivative_diagnostic"))
                {
                    map.derivativeDiagnostic = py::cast<std::string>(
                        value["derivative_diagnostic"]
                    );
                }
                map.geometryHash = geometryHash;
                map.topologyHash = topologyHash;
                map.evaluationMapIdentity = evaluationMapIdentity;
                map.booleanBackendIdentity = booleanBackendIdentity;
                map.derivativeIdentity = derivativeIdentity;
                if (value.contains("residual_evaluation_seconds"))
                {
                    map.residualEvaluationSeconds = py::cast<double>(
                        value["residual_evaluation_seconds"]
                    );
                }
                if (value.contains("derivative_evaluation_seconds"))
                {
                    map.derivativeEvaluationSeconds = py::cast<double>(
                        value["derivative_evaluation_seconds"]
                    );
                }
                return parsed;
                }
                catch (const py::error_already_set& error)
                {
                    if (!directional) throw; // Historical full-space contract.
                    const std::string type = py::str(error.type().attr("__module__")).cast<std::string>()
                        + "." + py::str(error.type().attr("__qualname__")).cast<std::string>();
                    const std::string message = py::str(error.value()).cast<std::string>();
                    py::dict diagnostic;
                    diagnostic["evaluation_index"] = callbackIndex;
                    diagnostic["exception_type"] = type;
                    diagnostic["message"] = message;
                    diagnostic["traceback"] = error.what();
                    pythonCallbackDiagnostics.append(std::move(diagnostic));
                    // Both callback and mapping-extraction tracebacks are
                    // diagnostics, not scientific identity.
                    throw std::runtime_error(type + ": " + message);
                }
            };

            spin::SpinNewtonEventSink eventSink;
            if (!pythonEventSink.is_none())
            {
                const py::function callback = py::cast<py::function>(pythonEventSink);
                eventSink = [callback](const spin::SpinSolverEvent& event) {
                    callback(solverEventToPython(event));
                };
            }
            std::optional<spin::SpinNewtonSubspaceResult> supported;
            spin::SpinNewtonResult solved;
            if (directional)
            {
                supported = spin::solveSpinMomentSubspaceSystem(problem,request.settings,
                    supportBasis,mapEvaluator,guesses,eventSink);
                solved = supported->solver;
            }
            else
            {
#if MAKE_IT_STAND_HAS_SPIN_IT_2014
                if (originalObjective)
                    solved = spin::solveMassMomentSystem(
                        spin2014::makeSpin2014QuadricMomentProblem(*originalObjective),request.settings,
                        [&](const spin::QuadricCoefficients10& q){return mapEvaluator(q).material;},
                        guesses,eventSink);
                else
#endif
                solved = spin::solveSpinMomentSystem(problem,request.settings,
                    [&](const spin::QuadricCoefficients10& q){return mapEvaluator(q).material;},
                    guesses,eventSink);
            }
            py::dict result;
#if MAKE_IT_STAND_HAS_SPIN_IT_2014
            if (originalObjective)
            {
                result["objective_profile"] = spin2014::spin2014QuadricObjectiveProfile();
                result["objective_identity"] = spin2014::spin2014QuadricObjectiveHash(*originalObjective);
            }
#endif
            if (supported)
            {
                py::dict attachment;
                attachment["version"] = supported->attachmentVersion;
                attachment["support_basis"] = supported->supportBasis;
                attachment["scientific_hash"] = supported->scientificHash;
                attachment["python_callback_diagnostics"] = pythonCallbackDiagnostics;
                py::list records;
                for (const auto& record : supported->evaluations)
                {
                    py::dict row;
                    row["evaluation_index"] = record.evaluationIndex;
                    row["quadric_coefficients"] = record.quadricCoefficients;
                    row["error"] = record.error;
                    row["error_type"] = record.errorType;
                    row["evaluation"] = py::none();
                    row["linearization"] = py::none();
                    if (record.evaluation)
                    {
                        const auto& e = *record.evaluation;
                        py::dict value;
                        value["raw_moments"] = e.material.rawMoments;
                        value["variable_raw_moments"] = e.material.variableRawMoments
                            ? py::cast(*e.material.variableRawMoments) : py::none();
                        value["native_tangent_basis"] = e.nativeTangentBasis;
                        value["native_directional_jacobian"] = e.nativeDirectionalJacobian;
                        value["derivative_defined"] = e.material.derivativeDefined;
                        value["derivative_diagnostic"] = e.material.derivativeDiagnostic;
                        value["geometry_hash"] = e.material.geometryHash;
                        value["topology_hash"] = e.material.topologyHash;
                        value["evaluation_map_identity"] = e.material.evaluationMapIdentity;
                        value["boolean_backend_identity"] = e.material.booleanBackendIdentity;
                        value["derivative_identity"] = e.material.derivativeIdentity;
                        row["evaluation"] = std::move(value);
                    }
                    if (record.linearization)
                    {
                        const auto& l = *record.linearization;
                        py::dict value;
                        value["solver_tangent"] = l.solverTangent;
                        value["moment_tangent_jacobian"] = l.momentTangentJacobian;
                        value["multipliers"] = l.multipliers;
                        value["multiplier_mu"] = l.multiplierMu;
                        value["full_residual"] = l.fullResidual;
                        value["raw_residual"] = l.rawResidual;
                        value["matrix"] = l.matrix;
                        value["right_hand_side"] = l.rightHandSide;
                        row["linearization"] = std::move(value);
                    }
                    records.append(std::move(row));
                }
                attachment["evaluations"] = std::move(records);
                result["supported_space"] = std::move(attachment);
            }
            result["status"] = spin::spinSolverStatusToken(solved.status);
            result["diagnostic"] = solved.diagnostic;
            result["function_evaluations"] = solved.functionEvaluations;
            result["restart_count"] = solved.restartCount;
            result["restart_function_evaluations"] = solved.restartFunctionEvaluations;
            result["normalized_equation_residual"] = solved.normalizedEquationResidual;
            result["physical_relative_residual"] = solved.physicalRelativeResidual;
            result["raw_equation_residual"] = solved.rawEquationResidual
                ? py::cast(*solved.rawEquationResidual) : py::none();
            result["multiplier_mu"] = solved.multiplierMu;
            result["multipliers"] = Eigen::VectorXd(solved.multipliers);
            result["quadric_coefficients"] = solved.quadric
                ? py::cast(Eigen::VectorXd(solved.quadric->values())) : py::none();
            result["equation_profile"] = spin::spinEquationProfileToken(
                solved.equationProfile
            );
            result["event_hash"] = solved.eventHash;
            result["trace_hash"] = solved.traceHash;
            if (request.settings.captureEvaluationPayload)
                result["event_stream_version"] = solved.eventStreamVersion;
            result["final_geometry_verified"] = false;
            result["scope"] = "numerical_newton_only_external_geometry_not_sealed";

            py::list trace;
            for (const auto& item : solved.trace)
            {
                py::dict entry;
                entry["restart"] = item.restart;
                entry["iteration"] = item.iteration;
                entry["function_evaluations"] = item.functionEvaluations;
                entry["normalized_equation_residual"] = item.normalizedEquationResidual;
                entry["physical_relative_residual"] = item.physicalRelativeResidual;
                entry["objective"] = item.objective;
                entry["multiplier_mu"] = item.multiplierMu;
                entry["accepted_geodesic_step"] = item.acceptedGeodesicStep;
                entry["quadric_hash"] = item.quadricHash;
                entry["material_geometry_hash"] = item.materialGeometryHash;
                entry["material_topology_hash"] = item.materialTopologyHash;
                entry["equation_profile"] = spin::spinEquationProfileToken(
                    item.equationProfile
                );
                entry["solver_coordinate_identity"] =
                    spin::spinSolverCoordinateIdentity(item.equationProfile);
                entry["raw_equation_residual"] = item.rawEquationResidual
                    ? py::cast(*item.rawEquationResidual) : py::none();
                if (item.quadricCoefficients)
                {
                    const Eigen::Matrix<double, 10, 1> coefficients = Eigen::Map<
                        const Eigen::Matrix<double, 10, 1>
                    >(item.quadricCoefficients->data());
                    entry["quadric_coefficients"] = coefficients;
                }
                else
                {
                    entry["quadric_coefficients"] = py::none();
                }
                trace.append(std::move(entry));
            }
            result["trace"] = std::move(trace);
            py::list events;
            for (const auto& event : solved.events)
            {
                py::dict entry;
                entry["type"] = event.type;
                entry["restart"] = event.restart;
                entry["iteration"] = event.iteration;
                entry["evaluation_index"] = event.evaluationIndex;
                entry["restart_evaluation_index"] = event.restartEvaluationIndex;
                entry["step"] = event.step ? py::cast(*event.step) : py::none();
                entry["outcome"] = event.outcome;
                entry["reason"] = event.reason;
                entry["detail"] = event.detail;
                entry["quadric_hash"] = event.quadricHash;
                entry["normalized_equation_residual"] = event.normalizedEquationResidual
                    ? py::cast(*event.normalizedEquationResidual) : py::none();
                entry["physical_relative_residual"] = event.physicalRelativeResidual
                    ? py::cast(*event.physicalRelativeResidual) : py::none();
                entry["objective"] = event.objective
                    ? py::cast(*event.objective) : py::none();
                entry["multiplier_mu"] = event.multiplierMu
                    ? py::cast(*event.multiplierMu) : py::none();
                entry["material_geometry_hash"] = event.materialGeometryHash
                    ? py::cast(*event.materialGeometryHash) : py::none();
                entry["material_topology_hash"] = event.materialTopologyHash
                    ? py::cast(*event.materialTopologyHash) : py::none();
                entry["equation_profile"] = spin::spinEquationProfileToken(
                    event.equationProfile
                );
                entry["solver_coordinate_identity"] =
                    spin::spinSolverCoordinateIdentity(event.equationProfile);
                entry["raw_equation_residual"] = event.rawEquationResidual
                    ? py::cast(*event.rawEquationResidual) : py::none();
                if (event.quadricCoefficients)
                {
                    const Eigen::Matrix<double, 10, 1> coefficients = Eigen::Map<
                        const Eigen::Matrix<double, 10, 1>
                    >(event.quadricCoefficients->data());
                    entry["quadric_coefficients"] = coefficients;
                }
                else
                {
                    entry["quadric_coefficients"] = py::none();
                }
                if (event.evaluationPayload)
                    entry["evaluation_payload"] =
                        evaluationPayloadToPython(*event.evaluationPayload);
                if (event.newtonSystemPayload)
                    entry["newton_system_payload"] =
                        newtonSystemPayloadToPython(*event.newtonSystemPayload);
                events.append(std::move(entry));
            }
            result["events"] = std::move(events);
            py::list evaluationTimings;
            for (const auto& timing : solved.evaluationTimings)
            {
                py::dict entry;
                entry["evaluation_index"] = timing.evaluationIndex;
                entry["restart"] = timing.restart;
                entry["restart_evaluation_index"] = timing.restartEvaluationIndex;
                entry["evaluator_returned"] = timing.evaluatorReturned;
                entry["wall_seconds"] = timing.wallSeconds;
                entry["reported_residual_seconds"] = timing.reportedResidualSeconds
                    ? py::cast(*timing.reportedResidualSeconds) : py::none();
                entry["reported_derivative_seconds"] = timing.reportedDerivativeSeconds
                    ? py::cast(*timing.reportedDerivativeSeconds) : py::none();
                evaluationTimings.append(std::move(entry));
            }
            result["evaluation_timings"] = std::move(evaluationTimings);
            return result;
        },
        py::arg("request"), py::arg("evaluator"),
        py::arg("initial_quadrics") = Eigen::MatrixXd(0, 10),
        py::arg("event_sink") = py::none(),
        py::arg("support_basis") = py::none(),
        py::arg("objective_settings") = py::none(),
        "Run strict Eq.14 Newton with an external per-evaluation geometry/moment map. "
        "The callback must return moments and a Jacobian from the same geometry call. "
        "This numerical result does not certify or seal final geometry."
    );
#endif

#if MAKE_IT_STAND_HAS_SPIN_IT_2014
    module.def(
        "differentiate_surface_raw_moments",
        [](const py::object& vertices, const py::object& faces, double density)
        {
            const auto value = spin2014::differentiateSurfaceRawMoments(
                verticesFromPython(vertices), facesFromPython(faces), density
            );
            py::dict result;
            result["moments"] = Eigen::VectorXd(value.moments.values());
            result["jacobian"] = value.jacobian;
            result["order"] = py::make_tuple(
                "M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY"
            );
            result["algorithm"] =
                "spin-it-2014-supplement-algorithms-1-2-surface-polynomial";
            result["coordinates"] = "interleaved-vertex-xyz";
            return result;
        },
        py::arg("vertices"), py::arg("faces"), py::arg("density") = 1.0
    );

    #if !MAKE_IT_STAND_HAS_SPIN_IT_FASTER
    py::enum_<spin::SpinMode>(module, "SpinMode")
        .value("STATIC_BALANCE", spin::SpinMode::StaticBalance)
        .value("FIDGET_SPINNER", spin::SpinMode::FidgetSpinner)
        .value("YO_YO", spin::SpinMode::YoYo)
        .value("SPINNING_TOP", spin::SpinMode::SpinningTop);
    #endif

    py::enum_<api::SpinIt2014Evidence>(module, "SpinIt2014Evidence")
        .value(
            "PAPER_FORMULA_REPRODUCTION",
            api::SpinIt2014Evidence::PaperFormulaReproduction
        )
        .value(
            "RECONSTRUCTED_ASSET",
            api::SpinIt2014Evidence::ReconstructedAsset
        );

    py::enum_<spin2014::Spin2014Objective>(module, "Spin2014Objective")
        .value("YOYO", spin2014::Spin2014Objective::YoYo)
        .value(
            "TOP_INERTIA_ONLY",
            spin2014::Spin2014Objective::TopInertiaOnly
        )
        .value(
            "TOP_LOWERING_ONLY",
            spin2014::Spin2014Objective::TopLoweringOnly
        )
        .value(
            "TOP_MASS_AND_LOWERING",
            spin2014::Spin2014Objective::TopMassAndLowering
        )
        .value(
            "STATIC_LOWERING",
            spin2014::Spin2014Objective::StaticLowering
        );

#if MAKE_IT_STAND_HAS_SPIN_IT_FASTER
    py::class_<spin2014::QuadricObjectiveSettings>(module, "Spin2014QuadricObjectiveSettings")
        .def(py::init<>())
        .def_readwrite("objective", &spin2014::QuadricObjectiveSettings::objective)
        .def_readwrite("inertia_weight", &spin2014::QuadricObjectiveSettings::inertiaWeight)
        .def_readwrite("center_weight", &spin2014::QuadricObjectiveSettings::centerWeight)
        .def_readwrite("declared_length_unit", &spin2014::QuadricObjectiveSettings::declaredLengthUnit)
        .def_readwrite("declared_mass_unit", &spin2014::QuadricObjectiveSettings::declaredMassUnit)
        .def_readwrite("solver_mass_unit_in_declared_units", &spin2014::QuadricObjectiveSettings::solverMassUnitInDeclaredUnits)
        .def_readwrite("solver_length_unit_in_declared_units", &spin2014::QuadricObjectiveSettings::solverLengthUnitInDeclaredUnits);
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    registerSpin2014AuthorityBindings(module);
#endif

    module.def("solve_spin2014_static_quadric_hollowing", [](
        const spin::SpinDesignRequest& request,
        const spin2014::QuadricObjectiveSettings& objective,
        double fixedDensity, const std::optional<Eigen::Vector4d>& initial,
        int maximumIterations, int maximumLineSearchSteps, double tolerance)
    {
        spin::PlaneStaticBalanceSettings settings;
        settings.maximumIterations=maximumIterations;
        settings.maximumLineSearchSteps=maximumLineSearchSteps;
        settings.constraintTolerance=tolerance;
        spin2014::QuadricStaticMaterialRun run;
        {
            py::gil_scoped_release release;
            run=spin2014::solveSpin2014StaticQuadricHollowing(
                request,objective,fixedDensity,initial,settings);
        }
        auto candidateRequest=run.measuredRequest;
        candidateRequest.staticObjective=spin::StaticObjective::LowerCenterOfMass;
        candidateRequest.staticObjectiveMeasure=spin::StaticObjectiveMeasure::DensityWeightedMass;
        candidateRequest.staticParetoWeight=0;
        py::dict result, checked;
        const auto& v=run.squaredHeight;
        checked["objective"]=v.evaluation.objective;
        checked["objective_gradient"]=v.evaluation.objectiveGradient;
        checked["objective_hessian"]=v.evaluation.objectiveHessian;
        checked["constraints"]=v.evaluation.constraints;
        checked["constraint_jacobian"]=v.evaluation.constraintJacobian;
        checked["equality_multipliers"]=v.equalityMultipliers;
        checked["multiplier"]=v.multiplier;
        checked["stationarity_residual"]=v.stationarityResidual;
        checked["degenerate_objective"]=v.degenerateObjective;
        checked["squared_height_conditions_verified"]=v.squaredHeightConditionsVerified;
        checked["boundary_witness_verified"]=v.boundaryWitnessVerified;
        checked["boundary_witness_coefficients"]=v.boundaryWitnessCoefficients;
        checked["boundary_multiplier"]=v.boundaryMultiplier;
        checked["boundary_gradient_minimum"]=v.boundaryGradientMinimum;
        checked["diagnostic"]=v.diagnostic;
        result["schema"]="spin2014-static-quadric-material-run-v1";
        result["objective_profile"]=run.objectiveProfile;
        result["objective_identity"]=run.objectiveIdentity;
        result["fixed_shell_moments"]=Eigen::VectorXd(run.measuredRequest.fixedShellMoments.values());
        result["fixed_material_density"]=run.fixedMaterialDensity;
        result["target_length"]=run.targetLength;
        result["initialization"]=run.initialization;
        result["initial_coefficients"]=run.initialCoefficients;
        result["candidate"]=staticPlaneResultToPython(run.candidate,candidateRequest);
        result["squared_height"]=std::move(checked);
        return result;
    }, py::arg("request"),py::arg("objective_settings"),py::arg("fixed_material_density"),
       py::arg("initial_coefficients")=py::none(),py::arg("maximum_iterations")=60,
       py::arg("maximum_line_search_steps")=24,py::arg("constraint_tolerance")=1e-10);

    module.def("solve_spin2014_quadric_hollowing", [](
        const spin::SpinDesignRequest& request,
        const spin2014::QuadricObjectiveSettings& settings,
        double fixedDensity, const spin::QuadricMomentEvaluatorSettings& evaluator,
        const Eigen::MatrixXd& initialQuadrics, const py::object& pythonEventSink)
    {
        if (!pythonEventSink.is_none() && !PyCallable_Check(pythonEventSink.ptr()))
            throw py::value_error("event_sink must be callable or None");
        if (initialQuadrics.cols() != 10 || !initialQuadrics.allFinite())
            throw py::value_error("initial quadrics must be a finite N-by-10 matrix");
        std::vector<spin::QuadricCoefficients10> initial;
        for (Eigen::Index row = 0; row < initialQuadrics.rows(); ++row)
            initial.push_back(spin::QuadricCoefficients10::fromUnnormalized(
                initialQuadrics.row(row).transpose()));
        spin2014::QuadricMaterialRun run;
        spin::SpinNewtonEventSink eventSink;
        if (!pythonEventSink.is_none()) {
            const py::function callback = py::cast<py::function>(pythonEventSink);
            eventSink = [callback](const spin::SpinSolverEvent& event) {
                py::gil_scoped_acquire acquire;
                callback(solverEventToPython(event));
            };
        }
        {
            py::gil_scoped_release release;
            run = spin2014::solveSpin2014QuadricHollowing(
                request, settings, fixedDensity, evaluator, initial, eventSink);
        }
        py::dict result;
        result["schema"] = "spin2014-quadric-material-run-v1";
        result["objective_profile"] = run.objectiveProfile;
        result["objective_identity"] = run.objectiveIdentity;
        result["objective_token"] = spin2014::spin2014ObjectiveToken(settings.objective);
        result["fixed_shell_moments"] = Eigen::VectorXd(run.measuredRequest.fixedShellMoments.values());
        result["fixed_material_density"] = run.fixedMaterialDensity;
        result["solve"] = fullRankResultToPython(run.solve, request.settings.captureEvaluationPayload);
        result["material_export_complete"] = run.materialExportComplete;
        result["material_export_diagnostic"] = run.materialExportDiagnostic;
        result["cavity_vertices"] = run.cavityEvaluation
            ? py::cast(Eigen::MatrixXd(run.cavityEvaluation->materialBoundaryTarget.vertices)) : py::none();
        result["cavity_faces"] = run.cavityEvaluation
            ? py::cast(Eigen::MatrixXi(run.cavityEvaluation->materialBoundaryTarget.faces)) : py::none();
        result["cavity_raw_moments"] = run.cavityEvaluation
            ? py::cast(Eigen::VectorXd(run.cavityEvaluation->variableMomentsTarget.values())) : py::none();
        return result;
    }, py::arg("request"), py::arg("objective_settings"),
       py::arg("fixed_material_density"), py::arg("evaluator_settings"),
       py::arg("initial_quadrics") = Eigen::MatrixXd(0, 10),
       py::arg("event_sink") = py::none());
#endif

    module.def(
        "evaluate_spin_it_2014_givens",
        [](
            spin2014::Spin2014Objective objective,
            const Eigen::Matrix<double, 10, 1>& rawMoments,
            double inertiaWeight,
            double centerWeight,
            double phi
        )
        {
            if (
                objective == spin2014::Spin2014Objective::StaticLowering &&
                phi != 0.0
            )
            {
                throw std::invalid_argument(
                    "static Spin-It 2014 Givens evaluation requires phi=0"
                );
            }
            spin2014::Spin2014GivensProblem::Vector11 state;
            state.head<10>() = rawMoments;
            state(10) = phi;
            const spin2014::Spin2014GivensEvaluation evaluation =
                spin2014::Spin2014GivensProblem(
                    objective, inertiaWeight, centerWeight
                ).evaluate(state);
            py::list constraintHessians;
            for (const auto& hessian : evaluation.constraintHessians)
            {
                constraintHessians.append(hessian);
            }
            py::dict result;
            result["objective_token"] = spin2014::spin2014ObjectiveToken(objective);
            result["formulation"] = "PAPER_GIVENS";
            result["objective"] = evaluation.objective;
            result["objective_gradient"] = evaluation.objectiveGradient;
            result["objective_hessian"] = evaluation.objectiveHessian;
            result["constraints"] = evaluation.constraints;
            result["constraint_jacobian"] = evaluation.constraintJacobian;
            result["constraint_hessians"] = constraintHessians;
            return result;
        },
        py::arg("objective"),
        py::arg("raw_moments"),
        py::arg("inertia_weight") = 1.0,
        py::arg("center_weight") = 1.0,
        py::arg("phi") = 0.0,
        "Evaluate the Spin-It 2014 paper Givens formulation without optimization."
    );

    module.def(
        "spin_it_2014_stability_threshold_squared",
        [](
            double mass,
            double gravity,
            double centerHeight,
            double axialInertia,
            double middleInertia
        )
        {
            if (
                !std::isfinite(mass) || !std::isfinite(gravity) ||
                !std::isfinite(centerHeight) || !std::isfinite(axialInertia) ||
                !std::isfinite(middleInertia) || mass <= 0.0 || gravity < 0.0 ||
                centerHeight < 0.0 || axialInertia <= middleInertia
            )
            {
                throw std::invalid_argument(
                    "Spin-It 2014 Eq. (1) requires M>0, g>=0, ell>=0 and Ic>Ib"
                );
            }
            return mass * gravity * centerHeight /
                (axialInertia - middleInertia);
        },
        py::arg("mass"), py::arg("gravity"), py::arg("center_height"),
        py::arg("axial_inertia"), py::arg("middle_inertia"),
        "Evaluate Spin-It 2014 Eq. (1)'s strict omega-squared stability bound."
    );

    py::class_<api::SpinIt2014Artifact>(module, "SpinIt2014Artifact")
        .def_property_readonly("status", &api::SpinIt2014Artifact::statusToken)
        .def_property_readonly("evidence", &api::SpinIt2014Artifact::evidenceToken)
        .def_property_readonly("succeeded", &api::SpinIt2014Artifact::succeeded)
        .def_property_readonly(
            "geometry_sealed", &api::SpinIt2014Artifact::geometrySealed
        )
        .def_property_readonly(
            "outer_vertices",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXd(value.outerVertices());
            }
        )
        .def_property_readonly(
            "outer_faces",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXi(value.outerFaces());
            }
        )
        .def_property_readonly(
            "design_vertices",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXd(value.designVertices());
            }
        )
        .def_property_readonly(
            "design_faces",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXi(value.designFaces());
            }
        )
        .def_property_readonly(
            "material_vertices",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXd(value.materialVertices());
            }
        )
        .def_property_readonly(
            "material_faces",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXi(value.materialFaces());
            }
        )
        .def_property_readonly(
            "cavity_vertices",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXd(value.cavityVertices());
            }
        )
        .def_property_readonly(
            "cavity_faces",
            [](const api::SpinIt2014Artifact& value)
            {
                return Eigen::MatrixXi(value.cavityFaces());
            }
        )
        .def_property_readonly(
            "diagnostics_json",
            [](const api::SpinIt2014Artifact& value)
            {
                return std::string(value.diagnosticsJson());
            }
        )
        .def("save_outer_mesh", &api::SpinIt2014Artifact::saveOuterMesh)
        .def("save_design_mesh", &api::SpinIt2014Artifact::saveDesignMesh)
        .def("save_material_mesh", &api::SpinIt2014Artifact::saveMaterialMesh)
        .def("save_cavity_mesh", &api::SpinIt2014Artifact::saveCavityMesh)
        .def("save_diagnostics", &api::SpinIt2014Artifact::saveDiagnostics);

    module.def(
        "solve_spin_it_2014",
        [](
            const Eigen::Vector3i& dimensions,
            const Eigen::Vector3d& minimum,
            const Eigen::Vector3d& maximum,
            const std::vector<unsigned char>& materialMask,
            const std::vector<unsigned char>& carveableMask,
            spin::SpinMode mode,
            api::SpinIt2014Evidence evidence,
            int initialLevel,
            int maximumLevel,
            int maximumAdaptivePasses,
            double initialFill,
            double fractionalEpsilon,
            double laplacianWeight,
            double centerOfMassWeight,
            double boundaryDensity,
            double interiorDensity,
            int maximumLeafCount,
            double binaryConstraintTolerance,
            int slqpMaximumIterations,
            double slqpProjectedGradientTolerance,
            const std::optional<spin2014::Spin2014Objective>& objective,
            double inertiaWeight,
            const std::string& formulation,
            double initialPhi,
            const std::string& slqpAlgorithm,
            double restorationTrustBound,
            int lpIterationLimit,
            double lpTimeLimit,
            double lpFeasibilityTolerance,
            const Eigen::Vector3d& targetOrigin,
            const Eigen::Matrix3d& targetAxes
        )
        {
            voxel::GridSpec grid;
            grid.dimensions = dimensions;
            grid.minimum = minimum;
            grid.maximum = maximum;
            spin2014::AdaptiveVoxelSpinRequest request;
            if(formulation == "PAPER_GIVENS") request.formulation=spin2014::Spin2014Formulation::PaperGivens;
            else if(formulation != "LEGACY_ANGLE_ELIMINATED") throw std::invalid_argument("unknown Spin-It 2014 formulation");
            request.initialPhi=initialPhi;
            if(slqpAlgorithm=="AUDITED_KKT_RESTORATION_V1")request.slqp.algorithm=spin2014::SlqpAlgorithm::AuditedKktRestorationV1;
            else if(slqpAlgorithm!="LEGACY_V3")throw std::invalid_argument("unknown SLQP algorithm");
            request.slqp.audit.restorationTrustBound=restorationTrustBound;
            request.slqp.audit.lpIterationLimit=lpIterationLimit;
            request.slqp.audit.lpTimeLimit=lpTimeLimit;
            request.slqp.audit.lpFeasibilityTolerance=lpFeasibilityTolerance;
            request.targetFrame.origin=targetOrigin;
            request.targetFrame.axes=targetAxes;
            request.domain = voxel::VoxelDomain::create(
                grid, materialMask, {}, carveableMask
            );
            request.mode = mode;
            request.objective = objective;
            request.initialLevel = initialLevel;
            request.maximumLevel = maximumLevel;
            request.maximumAdaptivePasses = maximumAdaptivePasses;
            request.initialFill = initialFill;
            request.fractionalEpsilon = fractionalEpsilon;
            request.laplacianWeight = laplacianWeight;
            request.centerOfMassWeight = centerOfMassWeight;
            request.inertiaWeight = inertiaWeight;
            request.boundaryDensity = boundaryDensity;
            request.interiorDensity = interiorDensity;
            request.maximumLeafCount = maximumLeafCount;
            request.binaryConstraintTolerance = binaryConstraintTolerance;
            request.slqp.maximumIterations = slqpMaximumIterations;
            request.slqp.projectedGradientTolerance =
                slqpProjectedGradientTolerance;
            py::gil_scoped_release release;
            const api::SpinIt2014Artifact artifact =
                api::solveSpinIt2014(request, evidence);
            return std::make_unique<api::SpinIt2014Artifact>(artifact);
        },
        py::arg("dimensions"),
        py::arg("minimum"),
        py::arg("maximum"),
        py::arg("material_mask"),
        py::arg("carveable_mask"),
        py::arg("mode"),
        py::arg("evidence"),
        py::arg("initial_level") = 3,
        py::arg("maximum_level") = 7,
        py::arg("maximum_adaptive_passes") = 8,
        py::arg("initial_fill") = 0.5,
        py::arg("fractional_epsilon") = 0.05,
        py::arg("laplacian_weight") = 1e-3,
        py::arg("center_of_mass_weight") = 1.0,
        py::arg("boundary_density") = 1.0,
        py::arg("interior_density") = 1.0,
        py::arg("maximum_leaf_count") = 1000000,
        py::arg("binary_constraint_tolerance") = 1e-6,
        py::arg("slqp_maximum_iterations") = 120,
        py::arg("slqp_projected_gradient_tolerance") = 1e-7,
        py::arg("objective") = py::none(),
        py::arg("inertia_weight") = 1.0,
        py::arg("formulation") = "LEGACY_ANGLE_ELIMINATED",
        py::arg("initial_phi") = 0.0,
        py::arg("slqp_algorithm") = "LEGACY_V3",
        py::arg("slqp_restoration_trust_bound") = 1.0,
        py::arg("slqp_lp_iteration_limit") = 10000,
        py::arg("slqp_lp_time_limit") = 10.0,
        py::arg("slqp_lp_feasibility_tolerance") = 1e-9,
        py::arg("target_origin") = Eigen::Vector3d::Zero(),
        py::arg("target_axes") = Eigen::Matrix3d::Identity()
    );

    module.def(
        "solve_spin_it_2014_cage",
        [](
            const api::SpinIt2014Artifact& artifact,
            const Eigen::Vector3i& dimensions,
            const Eigen::Vector3d& minimum,
            const Eigen::Vector3d& maximum,
            const std::vector<unsigned char>& materialMask,
            const std::vector<unsigned char>& carveableMask,
            double centerOfMassWeight,
            int maximumHandleCount,
            const std::optional<spin2014::Spin2014Objective>& objective,
            double inertiaWeight
        )
        {
            if(artifact.result().formulation == spin2014::Spin2014Formulation::PaperGivens)
                throw std::invalid_argument("paper Givens artifacts are unsupported by cage deformation");
            voxel::GridSpec grid;
            grid.dimensions = dimensions;
            grid.minimum = minimum;
            grid.maximum = maximum;
            const std::shared_ptr<const voxel::VoxelDomain> domain =
                voxel::VoxelDomain::create(
                    grid, materialMask, {}, carveableMask
                );
            if (
                artifact.cavityVertices().rows() == 0 ||
                artifact.cavityFaces().rows() == 0
            )
            {
                throw std::invalid_argument(
                    "cage deformation requires a nonempty hollowing cavity"
                );
            }

            spin::TriangleMesh outer;
            outer.vertices = artifact.outerVertices();
            outer.faces = artifact.outerFaces();
            spin::TriangleMesh cavity;
            cavity.vertices = artifact.cavityVertices();
            cavity.faces = artifact.cavityFaces();
            cavity.faces.col(1).swap(cavity.faces.col(2));

            spin2014::AutomaticCageBindingRequest bindingRequest;
            bindingRequest.domain = domain;
            bindingRequest.outerVertices = outer.vertices;
            bindingRequest.voidCellVertices = {cavity.vertices};
            bindingRequest.maximumHandleCount = maximumHandleCount;

            spin2014::CageSpinDeformationResult result;
            spin2014::AutomaticCageBinding binding;
            {
                py::gil_scoped_release release;
                binding = spin2014::buildAutomaticCageBinding(bindingRequest);
                spin2014::CageSpinDeformationRequest cageRequest;
                cageRequest.outerSurface = outer;
                cageRequest.voidCells = {cavity};
                cageRequest.outerWeights = binding.outerWeights;
                cageRequest.voidCellWeights = binding.voidCellWeights;
                cageRequest.handles = binding.handles;
                cageRequest.problem = objective.has_value()
                    ? spin2014::SpinIt2014MomentProblem::fromObjective(
                          *objective, inertiaWeight, centerOfMassWeight
                      )
                    : spin2014::SpinIt2014MomentProblem::spinningTop(
                          centerOfMassWeight
                      );
                cageRequest.fillLaplacian = Eigen::MatrixXd::Zero(1, 1);
                cageRequest.fillLaplacianWeight = 1e-3;
                cageRequest.biharmonicWeight = 1e-3;
                cageRequest.displacementWeight = 1e-3;
                cageRequest.minimumShellThickness = 0.0;
                const spin2014::CageSpinDeformationFunctional functional(
                    cageRequest
                );
                Eigen::VectorXd beta = Eigen::VectorXd::Ones(1);
                spin2014::CageSpinDeformationSolverSettings settings;
                settings.projectedGradientTolerance = 1e-6;
                result = spin2014::solveCageSpinDeformation(
                    functional, beta, settings
                );
            }
            py::dict payload;
            payload["status"] = spin2014::cageSpinStatusToken(result.status);
            payload["diagnostic"] = result.diagnostic;
            payload["geometry_valid"] = result.geometryValid;
            payload["silhouette_changed"] = result.silhouetteChanged;
            payload["constraint_infinity_norm"] =
                result.evaluation.constraints.size() > 0
                ? result.evaluation.constraints.lpNorm<Eigen::Infinity>()
                : 0.0;
            payload["trace_hash"] = result.traceHash;
            payload["output_hash"] = result.outputHash;
            payload["binding_hash"] = binding.bindingHash;
            payload["handle_count"] = binding.handles.size();
            payload["iteration_count"] = result.trace.size();
            payload["final_projected_gradient_infinity_norm"] =
                result.trace.empty()
                ? std::numeric_limits<double>::infinity()
                : result.trace.back().projectedGradientInfinityNorm;
            payload["deformed_outer_vertices"] =
                result.evaluation.deformedOuterVertices;
            payload["deformed_cavity_vertices"] =
                result.evaluation.deformedVoidCellVertices.empty()
                ? Eigen::MatrixXd(0, 3)
                : result.evaluation.deformedVoidCellVertices.front();
            return payload;
        },
        py::arg("artifact"),
        py::arg("dimensions"),
        py::arg("minimum"),
        py::arg("maximum"),
        py::arg("material_mask"),
        py::arg("carveable_mask"),
        py::arg("center_of_mass_weight") = 1.0,
        py::arg("maximum_handle_count") = 12,
        py::arg("objective") = py::none(),
        py::arg("inertia_weight") = 1.0
    );
#endif

#if MAKE_IT_STAND_HAS_REACHABILITY
    py::class_<reachability::ReachabilityCell>(module, "ReachabilityCell")
        .def(py::init([](int id, double mass, const Eigen::Vector2d& moment) {
            return reachability::ReachabilityCell{id, mass, moment};
        }), py::arg("cell_id"), py::arg("mass"), py::arg("first_moment_2d"))
        .def_readonly("cell_id", &reachability::ReachabilityCell::cellIndex)
        .def_readonly("mass", &reachability::ReachabilityCell::mass)
        .def_property_readonly("first_moment_2d", [](const reachability::ReachabilityCell& value) {
            return Eigen::Vector2d(value.firstMoment2D);
        });
    py::class_<reachability::CarvingReachabilityInput>(module, "CarvingReachabilityInput")
        .def(py::init([](double mass, const Eigen::Vector2d& moment,
                        const std::vector<reachability::ReachabilityCell>& cells,
                        const Eigen::MatrixXd& polygon, double tolerance) {
            // Shape is checked before hashing; semantic validation belongs to
            // the analyzer/checker, so adversarial finite fixtures are possible.
            if (polygon.cols() != 2) throw std::invalid_argument("polygon must have two columns");
            return reachability::CarvingReachabilityInput{mass, moment, cells, polygon, tolerance};
        }), py::arg("fixed_mass"), py::arg("fixed_first_moment_2d"),
            py::arg("variable_cells"), py::arg("safe_polygon_2d"), py::arg("feasibility_tolerance") = 1e-9)
        .def_readonly("fixed_mass", &reachability::CarvingReachabilityInput::fixedMass)
        .def_readonly("feasibility_tolerance", &reachability::CarvingReachabilityInput::feasibilityTolerance)
        .def_property_readonly("fixed_first_moment_2d", [](const reachability::CarvingReachabilityInput& value) {
            return Eigen::Vector2d(value.fixedFirstMoment2D);
        })
        .def_property_readonly("variable_cells", [](const reachability::CarvingReachabilityInput& value) {
            return value.variableCells;
        })
        .def_property_readonly("safe_polygon_2d", [](const reachability::CarvingReachabilityInput& value) {
            return Eigen::MatrixXd(value.safePolygon2D);
        })
        .def_property_readonly("input_hash", &reachability::computeInputHash);
    py::class_<reachability::CarvingReachabilitySnapshot>(module, "CarvingReachabilitySnapshot")
        .def_property_readonly("input", [](const reachability::CarvingReachabilitySnapshot& value) {
            return value.input;
        })
        .def_readonly("input_hash", &reachability::CarvingReachabilitySnapshot::inputHash)
        .def_readonly("preparation_status", &reachability::CarvingReachabilitySnapshot::preparationStatus);
    module.def("analyze_carving_reachability_input", [](const reachability::CarvingReachabilityInput& input) {
        return reachability::CarvingReachabilityAnalyzer().analyze(input);
    }, py::arg("input"));

    py::class_<reachability::ReachabilityCertificate>(
        module,
        "ReachabilityCertificate"
    )
        .def_property_readonly(
            "weights",
            [](const reachability::ReachabilityCertificate& value)
            {
                return value.weights;
            }
        )
        .def_readonly(
            "minimum_weight",
            &reachability::ReachabilityCertificate::minimumWeight
        )
        .def_readonly(
            "normalization_residual",
            &reachability::ReachabilityCertificate::normalizationResidual
        )
        .def_readonly(
            "stationarity_infinity_norm",
            &reachability::ReachabilityCertificate::stationarityInfinityNorm
        )
        .def_readonly(
            "contradiction",
            &reachability::ReachabilityCertificate::contradiction
        )
        .def_readonly(
            "valid",
            &reachability::ReachabilityCertificate::valid
        );

    py::class_<reachability::CarvingReachabilityResult>(
        module,
        "CarvingReachabilityResult"
    )
        .def_readonly(
            "input_hash",
            &reachability::CarvingReachabilityResult::inputHash
        )
        .def_property_readonly(
            "status",
            [](const reachability::CarvingReachabilityResult& value)
            {
                return std::string(reachability::toString(value.status));
            }
        )
        .def_readonly(
            "input_cell_count",
            &reachability::CarvingReachabilityResult::inputCellCount
        )
        .def_readonly(
            "variable_cell_count",
            &reachability::CarvingReachabilityResult::variableCellCount
        )
        .def_readonly(
            "fixed_mass",
            &reachability::CarvingReachabilityResult::fixedMass
        )
        .def_readonly(
            "retained_fractions",
            &reachability::CarvingReachabilityResult::retainedFractions
        )
        .def_readonly(
            "final_mass",
            &reachability::CarvingReachabilityResult::finalMass
        )
        .def_readonly(
            "final_first_moment_2d",
            &reachability::CarvingReachabilityResult::finalFirstMoment2D
        )
        .def_readonly(
            "final_center_2d",
            &reachability::CarvingReachabilityResult::finalCenter2D
        )
        .def_readonly(
            "maximum_primal_violation",
            &reachability::CarvingReachabilityResult::maximumPrimalViolation
        )
        .def_readonly(
            "certificate",
            &reachability::CarvingReachabilityResult::certificate
        )
        .def_readonly(
            "primal_solver_status",
            &reachability::CarvingReachabilityResult::primalSolverStatus
        )
        .def_readonly(
            "certificate_solver_status",
            &reachability::CarvingReachabilityResult::certificateSolverStatus
        )
        .def_readonly(
            "certificate_system_nonzero_count",
            &reachability::CarvingReachabilityResult::certificateSystemNonzeroCount
        );
#endif

    py::class_<api::BalanceModelConfig>(module, "BalanceModelConfig")
        .def(py::init<>())
        .def_readwrite(
            "voxel_resolution",
            &api::BalanceModelConfig::voxelResolution
        )
        .def_readwrite(
            "minimum_wall_thickness",
            &api::BalanceModelConfig::minimumWallThickness
        )
        .def_readwrite(
            "fixed_unit_cube_voxel_grid",
            &api::BalanceModelConfig::fixedUnitCubeVoxelGrid
        )
        .def_property(
            "bbw_max_iterations",
            [](const api::BalanceModelConfig& config)
            {
                return config.bbw.maxIterations;
            },
            [](api::BalanceModelConfig& config, int value)
            {
                config.bbw.maxIterations = value;
            }
        )
        .def_property(
            "bbw_tolerance",
            [](const api::BalanceModelConfig& config)
            {
                return config.bbw.tolerance;
            },
            [](api::BalanceModelConfig& config, double value)
            {
                config.bbw.tolerance = value;
            }
        );

    py::class_<api::BalanceRunConfig>(module, "BalanceRunConfig")
        .def(py::init<>())
        .def_readwrite("initial_mu", &api::BalanceRunConfig::initialMu)
        .def_property(
            "enable_carving",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.enableCarving;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.enableCarving = value;
            }
        )
        .def_property(
            "enable_deformation",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.enableDeformation;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.enableDeformation = value;
            }
        )
        .def_property(
            "enable_laplacian",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.enableLaplacian;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.enableLaplacian = value;
            }
        )
        .def_property(
            "freeze_targets_after_initial_refresh",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.freezeTargetsAfterInitialRefresh;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.freezeTargetsAfterInitialRefresh = value;
            }
        )
        .def_property(
            "released_beta_outer_loop",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.releasedBetaOuterLoop;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.releasedBetaOuterLoop = value;
            }
        )
        .def_property(
            "maximum_carving_blocks",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.maximumCarvingBlocks;
            },
            [](api::BalanceRunConfig& config, int value)
            {
                config.optimizer.maximumCarvingBlocks = value;
            }
        )
        .def_property(
            "maximum_deformation_iterations_per_block",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.maximumDeformationIterationsPerBlock;
            },
            [](api::BalanceRunConfig& config, int value)
            {
                config.optimizer.maximumDeformationIterationsPerBlock = value;
            }
        )
        .def_property(
            "gradient_tolerance",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.gradientTolerance;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.gradientTolerance = value;
            }
        )
        .def_property(
            "shape_scale_lambda",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.shapeScaleLambda;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.shapeScaleLambda = value;
            }
        )
        .def_property(
            "density",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.density;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.density = value;
            }
        )
        .def_property(
            "minimum_positive_mass",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.minimumPositiveMass;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.minimumPositiveMass = value;
            }
        )
        .def_property(
            "minimum_positive_cell_mass",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.minimumPositiveCellMass;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.minimumPositiveCellMass = value;
            }
        )
        .def_property(
            "minimum_carving_removals_to_continue",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.minimumCarvingRemovalsToContinue;
            },
            [](api::BalanceRunConfig& config, int value)
            {
                config.optimizer.minimumCarvingRemovalsToContinue = value;
            }
        )
        .def_property(
            "initial_step",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.initialStep;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.lineSearch.initialStep = value;
            }
        )
        .def_property(
            "minimum_step",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.minimumStep;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.lineSearch.minimumStep = value;
            }
        )
        .def_property(
            "reduction_ratio",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.reductionRatio;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.lineSearch.reductionRatio = value;
            }
        )
        .def_property(
            "required_relative_decrease",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.requiredRelativeDecrease;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.lineSearch.requiredRelativeDecrease = value;
            }
        )
        .def_property(
            "mu_decrement",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.muDecrement;
            },
            [](api::BalanceRunConfig& config, double value)
            {
                config.optimizer.lineSearch.muDecrement = value;
            }
        )
        .def_property(
            "released_beta_semantics",
            [](const api::BalanceRunConfig& config)
            {
                return config.optimizer.lineSearch.releasedBetaSemantics;
            },
            [](api::BalanceRunConfig& config, bool value)
            {
                config.optimizer.lineSearch.releasedBetaSemantics = value;
            }
        );

    py::class_<api::BalanceResult>(module, "BalanceResult")
        .def_property_readonly("status", &api::BalanceResult::statusToken)
        .def_property_readonly("succeeded", &api::BalanceResult::succeeded)
        .def_property_readonly(
            "outer_vertices",
            [](const api::BalanceResult& result)
            {
                return Eigen::MatrixXd(result.outerVertices());
            }
        )
        .def_property_readonly(
            "outer_faces",
            [](const api::BalanceResult& result)
            {
                return Eigen::MatrixXi(result.outerFaces());
            }
        )
        .def_property_readonly(
            "inner_vertices",
            [](const api::BalanceResult& result)
            {
                return Eigen::MatrixXd(result.innerVertices());
            }
        )
        .def_property_readonly(
            "inner_faces",
            [](const api::BalanceResult& result)
            {
                return Eigen::MatrixXi(result.innerFaces());
            }
        )
        .def_property_readonly(
            "manufacturing_ready",
            &api::BalanceResult::manufacturingReady
        )
        .def_property_readonly(
            "manufacturing_failures",
            [](const api::BalanceResult& result)
            {
                return std::vector<std::string>(
                    result.manufacturingFailures()
                );
            }
        )
        .def_property_readonly(
            "manufacturing_report_json",
            &api::BalanceResult::manufacturingReportJson
        )
        .def_property_readonly(
            "diagnostics_json",
            [](const api::BalanceResult& result)
            {
                return std::string(result.diagnosticsJson());
            }
        )
        .def_property_readonly(
            "performance_report_json",
            [](const api::BalanceResult& result)
            {
                return std::string(result.performanceReportJson());
            }
        )
        .def("save_outer_mesh", &api::BalanceResult::saveOuterMesh)
        .def("save_inner_mesh", &api::BalanceResult::saveInnerMesh)
        .def(
            "save_printable_mesh",
            py::overload_cast<const std::string&>(
                &api::BalanceResult::savePrintableMesh,
                py::const_
            )
        )
        .def(
            "save_printable_mesh",
            py::overload_cast<
                const std::string&,
                manufacturing::ManufacturingMeshFormat
            >(&api::BalanceResult::savePrintableMesh, py::const_),
            py::arg("filename"),
            py::arg("format")
        )
        .def("save_diagnostics", &api::BalanceResult::saveDiagnostics);

    py::class_<api::BalanceModel>(module, "BalanceModel")
        .def(
            py::init<const std::string&, const api::BalanceModelConfig&>(),
            py::arg("mesh_path"),
            py::arg("config") = api::BalanceModelConfig()
        )
        .def(
            py::init([](
                const py::object& vertices,
                const py::object& faces,
                const api::BalanceModelConfig& config
            )
            {
                return std::make_unique<api::BalanceModel>(
                    verticesFromPython(vertices),
                    facesFromPython(faces),
                    config
                );
            }),
            py::arg("vertices"),
            py::arg("faces"),
            py::arg("config") = api::BalanceModelConfig()
        )
        .def_static(
            "from_released_beta_archive",
            [](
                const py::object& normalizedVertices,
                const py::object& faces,
                int voxelResolution,
                const py::object& cellIdsZFast,
                const py::object& betaNodeWeightsZFast,
                int hullDepth,
                const api::BalanceModelConfig& config
            )
            {
                api::ReleasedBetaArchiveData archive;
                archive.voxelResolution = voxelResolution;
                archive.hullDepth = hullDepth;
                archive.cellIdsZFast = cellIdsFromPython(cellIdsZFast);
                archive.betaNodeWeightsZFast =
                    weightsFromPython(betaNodeWeightsZFast);
                return api::BalanceModel::fromReleasedBetaArchive(
                    verticesFromPython(normalizedVertices),
                    facesFromPython(faces),
                    archive,
                    config
                );
            },
            py::arg("normalized_vertices"),
            py::arg("faces"),
            py::arg("voxel_resolution"),
            py::arg("cell_ids_z_fast"),
            py::arg("beta_node_weights_z_fast"),
            py::arg("hull_depth") = 2,
            py::arg("config") = api::BalanceModelConfig()
        )
        .def_property_readonly("vertices", &api::BalanceModel::vertices)
        .def_property_readonly("faces", &api::BalanceModel::faces)
        .def_property_readonly("gravity", &api::BalanceModel::gravity)
        .def_property_readonly("config", &api::BalanceModel::config)
        .def("set_gravity", &api::BalanceModel::setGravity)
        .def(
            "set_support_polygon",
            &api::BalanceModel::setSupportPolygon,
            py::arg("polygon"),
            py::arg("gamma_radians")
        )
        .def(
            "set_support_from_mesh",
            &api::BalanceModel::setSupportFromMesh,
            py::arg("contact_threshold"),
            py::arg("gamma_radians"),
            py::arg("minimum_stability_inset") = 0.0,
            py::arg("support_plane_at_band_top") = false
        )
        .def(
            "set_support_from_reference_mesh",
            &api::BalanceModel::setSupportFromReferenceMesh,
            py::arg("reference_vertices"),
            py::arg("contact_threshold"),
            py::arg("gamma_radians"),
            py::arg("minimum_stability_inset") = 0.0,
            py::arg("support_plane_at_band_top") = false
        )
        .def(
            "add_standing_constraint_polygon",
            &api::BalanceModel::addStandingConstraintPolygon,
            py::arg("polygon"),
            py::arg("gravity_direction"),
            py::arg("gamma_radians")
        )
        .def(
            "add_standing_constraint_from_mesh",
            &api::BalanceModel::addStandingConstraintFromMesh,
            py::arg("gravity_direction"),
            py::arg("contact_threshold"),
            py::arg("gamma_radians"),
            py::arg("minimum_stability_inset") = 0.0,
            py::arg("support_plane_at_band_top") = false
        )
        .def(
            "add_handle",
            &api::BalanceModel::addHandle,
            py::arg("rest_pivot"),
            py::arg("fixed_rotation") = Eigen::Matrix3d::Identity()
        )
        .def(
            "add_suspension_constraint",
            &api::BalanceModel::addSuspensionConstraint,
            py::arg("attachment_point"),
            py::arg("allowed_deviation_radians")
        )
        .def("clear_handles", &api::BalanceModel::clearHandles)
        .def("clear_constraints", &api::BalanceModel::clearConstraints)
#if MAKE_IT_STAND_HAS_REACHABILITY
        .def(
            "carving_reachability_snapshot",
            &api::BalanceModel::carvingReachabilitySnapshot,
            py::arg("standing_index") = 0
        )
        .def(
            "analyze_carving_reachability",
            [](const api::BalanceModel& model, std::size_t standingIndex)
            {
                py::gil_scoped_release release;
                return model.analyzeCarvingReachability(standingIndex);
            },
            py::arg("standing_index") = 0
        )
#endif
        .def(
            "optimize",
            [](
                const api::BalanceModel& model,
                const api::BalanceRunConfig& config
            )
            {
                // Copy while Python still owns the GIL.  Releasing the GIL
                // with a reference into a mutable Python config object would
                // otherwise permit a C++ data race during optimize().
                const api::BalanceRunConfig snapshot = config;
                py::gil_scoped_release release;
                const api::BalanceResult result = model.optimize(snapshot);
                return std::make_unique<api::BalanceResult>(result);
            },
            py::arg("config") = api::BalanceRunConfig()
        );
}
