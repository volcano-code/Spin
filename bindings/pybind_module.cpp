#include "api/BalanceModel.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>
#include <memory>
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


} // namespace


PYBIND11_MODULE(_make_it_stand, module)
{
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
        py::arg("geometric_tolerance") = 0.0
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

#if MAKE_IT_STAND_HAS_REACHABILITY
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
