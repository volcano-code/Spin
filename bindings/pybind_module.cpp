#include "api/BalanceModel.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactMaterialRegion.h"
#include "ExactMeshSequenceIngestion.h"
#endif
#if MAKE_IT_STAND_HAS_SPIN_IT_2014
#include "api/SpinIt2014.h"
#include "spin2014/AutomaticCageBinding.h"
#include "spin2014/CageSpinDeformation.h"
#include "spin2014/CageSpinDeformationSolver.h"
#include "spin2014/SurfaceRawMomentDifferential.h"
#endif

#include <pybind11/eigen.h>
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
        });
    module.def("make_exact_material_region", [](const py::object& v,const py::object& f) {
        return verification::makeExactMaterialRegion(exactVerticesFromPython(v),exactFacesFromPython(f));
    },py::arg("vertices"),py::arg("faces"),
       "Authenticate an oriented finite regular closed region from exact binary64 indexed triangles. "
       "Inputs are NumPy matrices or nested list/tuple matrices whose rows may also be one-dimensional NumPy arrays. "
       "Negative checker copies only; input geometry is preserved. Empty mesh is malformed.");
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
        out["scope"]="exact-represented-polyhedral-set-only";
        return out;
    },py::arg("a"),py::arg("b"),
       "Exact inclusive contact/containment and positive-volume overlap. Invalid operands throw; "
       "failed Booleans return unknown facts. Touching does not imply volume overlap or invalid density regions.");
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
    module.def("audit_exact_triangle_boundary",
        [](const py::object& vertices, const py::object& faces, const std::string& orientation) {
            using Mode=verification::ExactBoundaryOrientation;
            Mode mode;
            if(orientation=="OutwardPositive") mode=Mode::OutwardPositive;
            else if(orientation=="InwardNegative") mode=Mode::InwardNegative;
            else if(orientation=="AnyNonzero") mode=Mode::AnyNonzero;
            else throw py::value_error("orientation must be OutwardPositive, InwardNegative, or AnyNonzero");
            const auto audit=verification::auditExactTriangleBoundary(
                exactVerticesFromPython(vertices), exactFacesFromPython(faces), mode);
            const auto& r=audit.facts();
            py::dict out;
            out["algorithm"]=r.algorithm; out["backend"]=r.backend;
            out["input_digest"]=r.inputDigest; out["input_digest_algorithm"]="sha256";
            out["input_encoding"]="exact-boundary-input-v1/binary64-row-major/u64-big-endian";
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
        }, py::arg("vertices"), py::arg("faces"), py::arg("orientation"),
        "Exact audit of the binary64 represented indexed boundary. No repair or welding. "
        "Inputs are NumPy matrices or nested list/tuple matrices whose rows may also be one-dimensional NumPy arrays. "
        "Float16/32 promote exactly; integer coordinates must be exactly binary64-representable. "
        "Wider floats are rejected. "
        "AnyNonzero is a boundary orientation mode only. Does not verify containment, "
        "nesting, thickness, material union, CAD, source/build identity, or artifact sealing.");
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
        .value("DISCRETE_ELLIPSOID", spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid);
    py::class_<spin::QuadricTriangulationSettings>(module, "SpinTriangulationSettings")
        .def(py::init<>())
        .def_readwrite("icosphere_subdivisions", &spin::QuadricTriangulationSettings::icosphereSubdivisions)
        .def_readwrite("paraboloid_grid_resolution", &spin::QuadricTriangulationSettings::paraboloidGridResolution)
        .def_readwrite("cylinder_angular_resolution", &spin::QuadricTriangulationSettings::cylinderAngularResolution)
        .def_readwrite("truncation_radius", &spin::QuadricTriangulationSettings::truncationRadius)
        .def_readwrite("maximum_radius_doublings", &spin::QuadricTriangulationSettings::maximumRadiusDoublings)
        .def_readwrite("retain_fixed_normal_provenance", &spin::QuadricTriangulationSettings::retainFixedNormalProvenance);
    py::class_<spin::QuadricMomentEvaluatorSettings>(module, "SpinEvaluatorSettings")
        .def(py::init<>())
        .def_readwrite("triangulation", &spin::QuadricMomentEvaluatorSettings::triangulation)
        .def_readwrite("minimum_interface_gradient", &spin::QuadricMomentEvaluatorSettings::minimumInterfaceGradient)
        .def_readwrite("derivative_method", &spin::QuadricMomentEvaluatorSettings::derivativeMethod);
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
        .def_readwrite("deterministic_seed", &spin::SpinSolverSettings::deterministicSeed);

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
        .def_readwrite("buoyancy", &spin::SpinDesignRequest::buoyancy)
        .def_readwrite("settings", &spin::SpinDesignRequest::settings);

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
            const std::optional<spin::QuadricMomentEvaluatorSettings>& explicitSettings
        )
        {
            if (coefficients.size() != 10 || !coefficients.allFinite())
            {
                throw py::value_error(
                    "quadric coefficients must contain 10 finite values"
                );
            }
            const spin::QuadricMomentEvaluation evaluation =
                spin::evaluateQuadricMaterialMoments(
                    request,
                    spin::QuadricCoefficients10::fromUnnormalized(coefficients),
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
                result["base_request_hash"] = spin::computeSpinDesignRequestHash(request);
                result["request_hash"] = spin::computeQuadricEvaluatorRequestHash(request, *explicitSettings);
                result["derivative_method"] = spin::quadricMomentDerivativeMethodToken(evaluation.derivativeMethod);
                result["derivative_defined"] = evaluation.derivativeDefined;
                result["derivative_diagnostic"] = evaluation.derivativeDiagnostic;
                result["derivative_scientific_hash"] = evaluation.derivativeScientificHash;
                result["provenance_topology_hash"] = evaluation.provenanceTopologyHash;
                result["vertex_provenance"] = evaluation.vertexProvenance;
                result["material_vertices"] = evaluation.materialBoundary.vertices;
                result["material_faces"] = evaluation.materialBoundary.faces;
                result["outer_vertices"] = request.outerSurface.vertices;
                result["outer_faces"] = request.outerSurface.faces;
            }
            result["interface_type"] = std::string(
                spin::quadricSurfaceTypeToken(
                    evaluation.quadricTriangulation.type
                )
            );
            result["interface_subtype"] = std::string(
                spin::quadricSurfaceSubtypeToken(
                    spin::classifyQuadricSurfaceSubtype(
                        spin::QuadricCoefficients10::fromUnnormalized(
                            coefficients
                        )
                    )
                )
            );
            return result;
        },
        py::arg("request"),
        py::arg("quadric_coefficients"),
        py::arg("evaluator_settings") = py::none()
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
            double initialPhi
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
        py::arg("initial_phi") = 0.0
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
