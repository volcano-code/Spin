#include "deformation/BoundedBiharmonicWeights.h"
#include "deformation/DeformationSystem.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "deformation/LinearBlendSkinning.h"
#include "extraction/VoidMeshExtractor.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{

using deformation::BoundedBiharmonicWeightResult;
using deformation::DeformationSystem;
using deformation::DeformedGeometry;
using deformation::GridDeformationDomain;
using deformation::Handle;
using deformation::HandleRole;
using deformation::HandleSet;
using deformation::HandleVector;
using deformation::WeightSolveOptions;


constexpr double SOLVER_TOLERANCE = 2.0e-8;
constexpr double GEOMETRY_TOLERANCE = 2.0e-10;


int fail(const std::string& scenario, const std::string& message)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << message << '\n';
    return EXIT_FAILURE;
}


bool approximatelyEqual(
    double first,
    double second,
    double tolerance
)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <=
            tolerance * std::max({1.0, std::abs(first), std::abs(second)});
}


bool matrixEqual(
    const Eigen::MatrixXd& first,
    const Eigen::MatrixXd& second,
    double tolerance
)
{
    return
        first.rows() == second.rows() &&
        first.cols() == second.cols() &&
        first.allFinite() &&
        second.allFinite() &&
        first.isApprox(second, tolerance);
}


template <typename Callable>
bool throwsInvalidArgument(Callable&& callable)
{
    try
    {
        callable();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }

    return false;
}


template <typename Callable>
bool throwsOutOfRange(Callable&& callable)
{
    try
    {
        callable();
    }
    catch (const std::out_of_range&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }

    return false;
}


void makeBox(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
    vertices.resize(8, 3);
    vertices <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0;

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


Eigen::MatrixXi pathEdges(int nodeCount)
{
    Eigen::MatrixXi edges(nodeCount - 1, 2);
    for (int node = 0; node + 1 < nodeCount; ++node)
    {
        edges.row(node) << node, node + 1;
    }
    return edges;
}


Eigen::MatrixXd linearRestNodes(int nodeCount)
{
    Eigen::MatrixXd nodes = Eigen::MatrixXd::Zero(nodeCount, 3);
    for (int node = 0; node < nodeCount; ++node)
    {
        nodes(node, 0) = static_cast<double>(node);
    }
    return nodes;
}


WeightSolveOptions strictSolverOptions()
{
    WeightSolveOptions options;
    options.maxIterations = 1000;
    options.tolerance = 1.0e-11;
    return options;
}


bool weightContractHolds(
    const BoundedBiharmonicWeightResult& result,
    const std::vector<std::vector<int>>& constraints
)
{
    if (!result.weights.allFinite())
    {
        return false;
    }

    for (Eigen::Index row = 0; row < result.weights.rows(); ++row)
    {
        if (
            result.weights.row(row).minCoeff() < -SOLVER_TOLERANCE ||
            result.weights.row(row).maxCoeff() > 1.0 + SOLVER_TOLERANCE ||
            !approximatelyEqual(
                result.weights.row(row).sum(),
                1.0,
                SOLVER_TOLERANCE
            )
        )
        {
            return false;
        }
    }

    for (std::size_t handle = 0; handle < constraints.size(); ++handle)
    {
        for (const int node : constraints[handle])
        {
            for (std::size_t column = 0;
                 column < constraints.size();
                 ++column)
            {
                const double expected = handle == column ? 1.0 : 0.0;
                if (!approximatelyEqual(
                    result.weights(node, static_cast<Eigen::Index>(column)),
                    expected,
                    SOLVER_TOLERANCE
                ))
                {
                    return false;
                }
            }
        }
    }

    return
        result.diagnostics.maximumBoundViolation <= SOLVER_TOLERANCE &&
        result.diagnostics.maximumRowSumError <= SOLVER_TOLERANCE &&
        result.diagnostics.maximumConstraintError <= SOLVER_TOLERANCE;
}


int runBbwPath(const std::string& scenario)
{
    const Eigen::MatrixXd nodes = linearRestNodes(4);
    const Eigen::MatrixXi edges = pathEdges(4);
    const std::vector<std::vector<int>> constraints{{0}, {3}};

    const BoundedBiharmonicWeightResult result =
        deformation::solveBoundedBiharmonicWeights(
            nodes,
            edges,
            constraints,
            strictSolverOptions()
        );

    Eigen::Vector4d expectedFirst;
    expectedFirst << 1.0, 0.7, 0.3, 0.0;

    if (
        !weightContractHolds(result, constraints) ||
        !result.weights.col(0).isApprox(expectedFirst, SOLVER_TOLERANCE) ||
        !result.weights.col(1).isApprox(
            Eigen::Vector4d::Ones() - expectedFirst,
            SOLVER_TOLERANCE
        )
    )
    {
        return fail(scenario, "four-node graph BBW oracle failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " weights=" << result.weights.col(0).transpose() << '\n';
    return EXIT_SUCCESS;
}


int runBbwActiveBound(const std::string& scenario)
{
    const Eigen::MatrixXd nodes = linearRestNodes(5);
    Eigen::MatrixXi edges(4, 2);
    edges << 0, 1, 0, 2, 0, 4, 1, 3;
    const std::vector<std::vector<int>> constraints{{1}, {4}};

    const BoundedBiharmonicWeightResult result =
        deformation::solveBoundedBiharmonicWeights(
            nodes,
            edges,
            constraints,
            strictSolverOptions()
        );

    Eigen::VectorXd expectedSecond(5);
    expectedSecond << 0.5, 0.0, 0.5, 0.0, 1.0;

    if (
        !weightContractHolds(result, constraints) ||
        !result.weights.col(1).isApprox(
            expectedSecond,
            SOLVER_TOLERANCE
        ) ||
        !approximatelyEqual(result.weights(3, 1), 0.0, SOLVER_TOLERANCE)
    )
    {
        return fail(scenario, "active lower bound was not enforced");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " bounded_handle=" << result.weights.col(1).transpose() << '\n';
    return EXIT_SUCCESS;
}


int runBbwContract(const std::string& scenario)
{
    const Eigen::MatrixXd nodes = linearRestNodes(6);
    const Eigen::MatrixXi edges = pathEdges(6);
    const std::vector<std::vector<int>> constraints{{0, 1}, {5}};
    const std::vector<std::vector<int>> reordered{{1, 0}, {5}};
    const std::vector<std::vector<int>> swapped{{5}, {1, 0}};

    const BoundedBiharmonicWeightResult first =
        deformation::solveBoundedBiharmonicWeights(
            nodes, edges, constraints, strictSolverOptions()
        );
    const BoundedBiharmonicWeightResult repeat =
        deformation::solveBoundedBiharmonicWeights(
            nodes, edges, reordered, strictSolverOptions()
        );
    const BoundedBiharmonicWeightResult swappedResult =
        deformation::solveBoundedBiharmonicWeights(
            nodes, edges, swapped, strictSolverOptions()
        );

    if (
        !weightContractHolds(first, constraints) ||
        !weightContractHolds(repeat, reordered) ||
        !weightContractHolds(swappedResult, swapped) ||
        !matrixEqual(first.weights, repeat.weights, SOLVER_TOLERANCE) ||
        !first.weights.col(0).isApprox(
            swappedResult.weights.col(1),
            SOLVER_TOLERANCE
        ) ||
        !first.weights.col(1).isApprox(
            swappedResult.weights.col(0),
            SOLVER_TOLERANCE
        )
    )
    {
        return fail(
            scenario,
            "weight feasibility, determinism, or handle equivariance failed"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " max_row_error=" << first.diagnostics.maximumRowSumError
        << '\n';
    return EXIT_SUCCESS;
}


Handle makeSupportHandle(const std::vector<int>& nodes)
{
    Handle handle;
    handle.restPivot = Eigen::Vector3d(0.0, 0.25, -0.5);
    handle.role = HandleRole::Support;
    handle.locked = true;
    handle.constrainedNodes = nodes;
    return handle;
}


Handle makeUserHandle(
    const std::vector<int>& nodes,
    const Eigen::Vector3d& pivot = Eigen::Vector3d::Zero()
)
{
    Handle handle;
    handle.restPivot = pivot;
    handle.role = HandleRole::User;
    handle.locked = false;
    handle.constrainedNodes = nodes;
    return handle;
}


int runHandles(const std::string& scenario)
{
    Handle support = makeSupportHandle({2, 0, 2});
    Handle user = makeUserHandle({7, 5}, Eigen::Vector3d(0.5, 0.5, 0.5));
    user.translation = Eigen::Vector3d(0.1, -0.2, 0.3);
    user.uniformScale = 1.2;

    Handle lockedUser = makeUserHandle({9});
    lockedUser.locked = true;
    lockedUser.translation = Eigen::Vector3d(-0.1, 0.0, 0.0);

    HandleSet handles(HandleVector{support, user, lockedUser});
    if (
        handles.fullParameterCount() != 12 ||
        handles.optimizableParameterCount() != 4 ||
        handles[0].constrainedNodes != std::vector<int>({0, 2})
    )
    {
        return fail(scenario, "handle binding/parameter contract failed");
    }

    Eigen::VectorXd updated(4);
    updated << 0.6, -0.4, 0.2, 0.75;
    handles.setOptimizationParameters(updated);
    if (!handles.optimizationParameters().isApprox(updated, 0.0))
    {
        return fail(scenario, "valid optimizer update was not committed");
    }

    const Eigen::VectorXd before = handles.optimizationParameters();
    Eigen::VectorXd invalid = before;
    invalid(3) = 0.0;
    const bool rejected = throwsInvalidArgument(
        [&]() { handles.setOptimizationParameters(invalid); }
    );

    if (
        !rejected ||
        !handles.optimizationParameters().isApprox(before, 0.0) ||
        !handles[0].translation.isZero(0.0) ||
        handles[0].uniformScale != 1.0
    )
    {
        return fail(scenario, "transactional update or fixed support failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " free_parameters=" << handles.optimizableParameterCount()
        << '\n';
    return EXIT_SUCCESS;
}


int runLbs(const std::string& scenario)
{
    Eigen::MatrixXd rest(4, 3);
    rest <<
        -1.0, 0.0, 0.5,
         0.0, 1.0, 0.0,
         1.0, 0.5, -0.5,
         0.2, -0.7, 1.3;

    Eigen::MatrixXd weights(4, 2);
    weights <<
        1.0, 0.0,
        0.0, 1.0,
        0.35, 0.65,
        0.8, 0.2;

    Handle first = makeUserHandle({0}, Eigen::Vector3d(-0.25, 0.5, 0.1));
    first.translation = Eigen::Vector3d(0.3, -0.1, 0.4);
    first.uniformScale = 1.35;
    first.fixedRotation =
        Eigen::AngleAxisd(0.4, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    Handle second = makeUserHandle({1}, Eigen::Vector3d(0.7, -0.2, 0.3));
    second.translation = Eigen::Vector3d(-0.2, 0.5, -0.1);
    second.uniformScale = 0.65;
    second.fixedRotation =
        Eigen::AngleAxisd(-0.55, Eigen::Vector3d::UnitY()).toRotationMatrix();

    const HandleSet handles(HandleVector{first, second});
    const Eigen::MatrixXd direct =
        deformation::deformLinearBlendSkinning(rest, weights, handles);
    const deformation::ReducedLbsModel reducedModel(
        rest, weights, handles
    );
    const Eigen::MatrixXd& reduced = reducedModel.matrix();
    const Eigen::MatrixXd viaEquationTen =
        reducedModel.evaluate(handles);

    const deformation::StackedParameterReduction parameterReduction =
        deformation::buildStackedParameterReduction(handles);
    const Eigen::VectorXd reconstructedParameters =
        parameterReduction.offset +
        parameterReduction.linear * handles.optimizationParameters();

    HandleVector commonTranslationValues = handles.handles();
    const Eigen::Vector3d commonTranslation(0.6, -0.3, 0.9);
    for (Handle& handle : commonTranslationValues)
    {
        handle.translation += commonTranslation;
    }
    const HandleSet commonTranslationHandles(
        std::move(commonTranslationValues)
    );
    const Eigen::MatrixXd commonTranslationResult =
        deformation::deformLinearBlendSkinning(
            rest,
            weights,
            commonTranslationHandles
        );

    Eigen::MatrixXd oracle = Eigen::MatrixXd::Zero(rest.rows(), 3);
    for (Eigen::Index vertex = 0; vertex < rest.rows(); ++vertex)
    {
        const Eigen::Vector3d point = rest.row(vertex).transpose();
        for (std::size_t handleIndex = 0;
             handleIndex < handles.size();
             ++handleIndex)
        {
            const Handle& handle = handles[handleIndex];
            const Eigen::Vector3d transformed =
                handle.restPivot + handle.translation +
                handle.uniformScale * handle.fixedRotation *
                    (point - handle.restPivot);
            oracle.row(vertex) +=
                weights(vertex, static_cast<Eigen::Index>(handleIndex)) *
                transformed.transpose();
        }
    }

    const Eigen::MatrixXd repeat =
        deformation::deformLinearBlendSkinning(rest, weights, handles);

    if (
        reduced.rows() != 3 * rest.rows() ||
        reduced.cols() != handles.fullParameterCount() ||
        !reconstructedParameters.isApprox(
            deformation::stackedHandleParameters(handles),
            GEOMETRY_TOLERANCE
        ) ||
        !matrixEqual(direct, oracle, GEOMETRY_TOLERANCE) ||
        !matrixEqual(direct, viaEquationTen, GEOMETRY_TOLERANCE) ||
        !matrixEqual(direct, repeat, 0.0) ||
        !matrixEqual(
            commonTranslationResult,
            direct.rowwise() + commonTranslation.transpose(),
            GEOMETRY_TOLERANCE
        )
    )
    {
        return fail(scenario, "Equation (9)/(10) equivalence failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " max_equation_error="
        << (direct - viaEquationTen).cwiseAbs().maxCoeff() << '\n';
    return EXIT_SUCCESS;
}


int runDomainTrilinear(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);
    const VoxelGrid grid(vertices, faces, 0.0, 10);
    const GridDeformationDomain domain =
        GridDeformationDomain::fromVoxelGrid(grid);

    int cell = -1;
    for (Eigen::Index index = 0;
         index < domain.activeCellMask().size();
         ++index)
    {
        if (domain.activeCellMask()(index) != 0)
        {
            cell = static_cast<int>(index);
            break;
        }
    }
    if (cell < 0)
    {
        return fail(scenario, "box produced no active deformation cell");
    }

    const Eigen::Vector3i cellCoordinates = domain.cellCoordinates(cell);
    Eigen::MatrixXd contact(1, 3);
    contact.row(0) = grid(
        cellCoordinates(0),
        cellCoordinates(1),
        cellCoordinates(2)
    ).center.transpose();
    const Handle fixedSupport =
        deformation::makeFixedSupportHandle(domain, contact);
    if (
        fixedSupport.role != HandleRole::Support ||
        !fixedSupport.locked ||
        fixedSupport.constrainedNodes.size() != 8U
    )
    {
        return fail(scenario, "support-cell binding failed");
    }

    const std::array<int, 8> nodes = domain.cellNodes(cell);
    const Eigen::Vector3d minimum =
        domain.restNodes().row(nodes[0]).transpose();
    const Eigen::Vector3d local(0.25, 0.5, 0.75);
    const Eigen::Vector3d point =
        minimum + domain.voxelSize().cwiseProduct(local);

    Eigen::MatrixXd nodeWeights =
        Eigen::MatrixXd::Zero(domain.nodeCount(), 2);
    nodeWeights.col(1).setOnes();
    nodeWeights(nodes[7], 0) = 1.0;
    nodeWeights(nodes[7], 1) = 0.0;

    Eigen::MatrixXd query(1, 3);
    query.row(0) = point.transpose();
    const Eigen::MatrixXd transferred =
        domain.interpolateWeights(query, nodeWeights);
    const double expected = local.prod();

    const Eigen::Vector3d outside =
        domain.fullRestCorners().colwise().minCoeff().transpose() -
        4.0 * domain.voxelSize();
    Eigen::MatrixXd outsideQuery(1, 3);
    outsideQuery.row(0) = outside.transpose();

    if (
        !approximatelyEqual(transferred(0, 0), expected, 1.0e-11) ||
        !approximatelyEqual(transferred.row(0).sum(), 1.0, 1.0e-12) ||
        !throwsOutOfRange(
            [&]() { domain.interpolateWeights(outsideQuery, nodeWeights); }
        )
    )
    {
        return fail(scenario, "trilinear transfer or domain rejection failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " active_nodes=" << domain.nodeCount()
        << " interpolated=" << transferred(0, 0) << '\n';
    return EXIT_SUCCESS;
}


struct GeometryBinding
{
    HandleSet handles;
    Eigen::MatrixXd weights;
};


GeometryBinding makeGeometryBinding(const GridDeformationDomain& domain)
{
    const double minimum = domain.restNodes().col(0).minCoeff();
    const double maximum = domain.restNodes().col(0).maxCoeff();
    const double span = maximum - minimum;

    Handle support = makeSupportHandle({0});
    Handle user = makeUserHandle({1});
    support.constrainedNodes.clear();
    user.constrainedNodes.clear();
    support.restPivot.setZero();
    user.restPivot = Eigen::Vector3d(maximum, 0.5, 0.5);

    Eigen::MatrixXd weights(domain.nodeCount(), 2);
    for (Eigen::Index node = 0; node < domain.nodeCount(); ++node)
    {
        const double coordinate = domain.restNodes()(node, 0);
        const double userWeight = (coordinate - minimum) / span;
        weights(node, 0) = 1.0 - userWeight;
        weights(node, 1) = userWeight;

        if (approximatelyEqual(coordinate, minimum, 1.0e-12))
        {
            support.constrainedNodes.push_back(static_cast<int>(node));
        }
        if (approximatelyEqual(coordinate, maximum, 1.0e-12))
        {
            user.constrainedNodes.push_back(static_cast<int>(node));
        }
    }

    return GeometryBinding{
        HandleSet(HandleVector{support, user}),
        weights
    };
}


int runGeometrySync(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);
    VoxelGrid grid(vertices, faces, 0.0, 12);
    const GridDeformationDomain domain =
        GridDeformationDomain::fromVoxelGrid(grid);
    const GeometryBinding binding = makeGeometryBinding(domain);

    bool carved = false;
    Eigen::Vector3i carvedCell = Eigen::Vector3i::Constant(-1);
    for (int x = 0; x < grid.dimensions()(0) && !carved; ++x)
    {
        for (int y = 0; y < grid.dimensions()(1) && !carved; ++y)
        {
            for (int z = 0; z < grid.dimensions()(2); ++z)
            {
                if (grid(x, y, z).isInner)
                {
                    grid(x, y, z).isFilled = false;
                    carvedCell = Eigen::Vector3i(x, y, z);
                    carved = true;
                    break;
                }
            }
        }
    }
    if (!carved)
    {
        return fail(scenario, "box fixture had no carveable interior cell");
    }

    const DeformationSystem system(
        grid,
        vertices,
        faces,
        binding.handles,
        binding.weights,
        true
    );
    if (
        !system.weightDiagnostics().verifiedImport ||
        system.weightDiagnostics().converged
    )
    {
        return fail(scenario, "imported-weight diagnostics are ambiguous");
    }

    HandleVector currentValues = binding.handles.handles();
    currentValues[1].translation = Eigen::Vector3d(0.15, -0.08, 0.12);
    currentValues[1].uniformScale = 0.85;
    currentValues[1].fixedRotation =
        Eigen::AngleAxisd(0.18, Eigen::Vector3d::UnitY()).toRotationMatrix();
    const HandleSet current(std::move(currentValues));

    const DeformedGeometry geometry = system.evaluate(current, grid);
    const DeformedGeometry repeat = system.evaluate(current, grid);
    const Eigen::MatrixXd expectedOuter =
        deformation::deformLinearBlendSkinning(
            vertices,
            system.outerWeights(),
            current
        );
    const Eigen::MatrixXd expectedActive =
        deformation::deformLinearBlendSkinning(
            system.domain().restNodes(),
            system.nodeWeights(),
            current
        );
    const Eigen::MatrixXd expectedFull =
        deformation::deformLinearBlendSkinning(
            system.domain().fullRestCorners(),
            system.fullCornerWeights(),
            current
        );

    for (Eigen::Index local = 0; local < expectedActive.rows(); ++local)
    {
        const int global = system.domain().localToGlobalCorner()(local);
        if (!geometry.gridCorners.row(global).isApprox(
            expectedActive.row(local),
            GEOMETRY_TOLERANCE
        ))
        {
            return fail(scenario, "active grid node diverged from common LBS");
        }
    }

    Eigen::MatrixXi expectedInnerFaces;
    VoidMeshExtractor::extractFaces(grid, expectedInnerFaces);
    const int linearCell = system.domain().cellLinearIndex(
        carvedCell(0), carvedCell(1), carvedCell(2)
    );
    Eigen::RowVector3d expectedCenter = Eigen::RowVector3d::Zero();
    for (const int global : grid(
        carvedCell(0), carvedCell(1), carvedCell(2)
    ).cornerIndices)
    {
        expectedCenter += geometry.gridCorners.row(global);
    }
    expectedCenter /= 8.0;

    bool supportFixed = true;
    for (const int local : binding.handles[0].constrainedNodes)
    {
        supportFixed = supportFixed &&
            expectedActive.row(local).isApprox(
                system.domain().restNodes().row(local),
                GEOMETRY_TOLERANCE
            );
    }

    if (
        geometry.outerFaces != faces ||
        expectedInnerFaces.rows() == 0 ||
        geometry.innerFaces != expectedInnerFaces ||
        !matrixEqual(geometry.outerVertices, expectedOuter, GEOMETRY_TOLERANCE) ||
        !matrixEqual(geometry.gridCorners, expectedFull, GEOMETRY_TOLERANCE) ||
        !matrixEqual(geometry.innerVertices, geometry.gridCorners, 0.0) ||
        !geometry.voxelCenters.row(linearCell).isApprox(
            expectedCenter,
            GEOMETRY_TOLERANCE
        ) ||
        !matrixEqual(geometry.outerVertices, repeat.outerVertices, 0.0) ||
        !matrixEqual(geometry.gridCorners, repeat.gridCorners, 0.0) ||
        !supportFixed ||
        grid(carvedCell(0), carvedCell(1), carvedCell(2)).isFilled
    )
    {
        return fail(scenario, "outer/grid/cavity synchronized snapshot failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " outer_vertices=" << geometry.outerVertices.rows()
        << " cavity_faces=" << geometry.innerFaces.rows() << '\n';
    return EXIT_SUCCESS;
}


int runInvalid(const std::string& scenario)
{
    bool acceptedInvalid = false;

    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        []() { HandleSet handles(HandleVector{}); }
    );

    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        []()
        {
            Handle support = makeSupportHandle({0});
            support.locked = false;
            HandleSet handles(HandleVector{support});
        }
    );

    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        []()
        {
            Handle user = makeUserHandle({0});
            user.fixedRotation(0, 0) = 2.0;
            HandleSet handles(HandleVector{user});
        }
    );

    const Eigen::MatrixXd nodes = linearRestNodes(3);
    const Eigen::MatrixXi edges = pathEdges(3);
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]()
        {
            deformation::solveBoundedBiharmonicWeights(
                nodes,
                edges,
                std::vector<std::vector<int>>{{0}, {0}}
            );
        }
    );

    Eigen::MatrixXi disconnected(1, 2);
    disconnected << 0, 1;
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]()
        {
            deformation::solveBoundedBiharmonicWeights(
                nodes,
                disconnected,
                std::vector<std::vector<int>>{{0}}
            );
        }
    );

    Eigen::MatrixXd rest(1, 3);
    rest.setZero();
    Eigen::MatrixXd invalidWeights(1, 1);
    invalidWeights(0, 0) = -0.1;
    const HandleSet oneHandle(HandleVector{makeUserHandle({0})});
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]()
        {
            deformation::deformLinearBlendSkinning(
                rest, invalidWeights, oneHandle
            );
        }
    );

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);
    const VoxelGrid grid(vertices, faces, 0.0, 8);
    const GridDeformationDomain domain =
        GridDeformationDomain::fromVoxelGrid(grid);
    Eigen::MatrixXd badNodeWeights =
        Eigen::MatrixXd::Ones(domain.nodeCount(), 1);
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]()
        {
            DeformationSystem system(
                grid,
                vertices,
                faces,
                oneHandle,
                badNodeWeights,
                false
            );
        }
    );

    const GeometryBinding binding = makeGeometryBinding(domain);
    Eigen::MatrixXd nonInterpolating =
        Eigen::MatrixXd::Constant(domain.nodeCount(), 2, 0.5);
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]()
        {
            DeformationSystem system(
                grid,
                vertices,
                faces,
                binding.handles,
                nonInterpolating,
                true
            );
        }
    );

    const DeformationSystem validSystem(
        grid,
        vertices,
        faces,
        binding.handles,
        binding.weights,
        true
    );
    VoxelGrid tampered = grid;
    tampered(0, 0, 0).isInner = !tampered(0, 0, 0).isInner;
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]() { validSystem.evaluate(binding.handles, tampered); }
    );

    const deformation::ReducedLbsModel reducedModel(
        vertices,
        Eigen::MatrixXd::Constant(
            vertices.rows(),
            binding.handles.size(),
            0.5
        ),
        binding.handles
    );
    HandleVector staleValues = binding.handles.handles();
    staleValues[1].fixedRotation =
        Eigen::AngleAxisd(0.2, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const HandleSet staleHandles(std::move(staleValues));
    acceptedInvalid = acceptedInvalid || !throwsInvalidArgument(
        [&]() { reducedModel.evaluate(staleHandles); }
    );

    VoxelGrid supportGrid = grid;
    const GridDeformationDomain supportDomain =
        GridDeformationDomain::fromVoxelGrid(supportGrid);
    int supportCell = -1;
    for (Eigen::Index cell = 0;
         cell < supportDomain.activeCellMask().size();
         ++cell)
    {
        if (supportDomain.activeCellMask()(cell) != 0)
        {
            supportCell = static_cast<int>(cell);
            break;
        }
    }
    const Eigen::Vector3i supportCoordinates =
        supportDomain.cellCoordinates(supportCell);
    Eigen::MatrixXd contact(1, 3);
    contact.row(0) = supportGrid(
        supportCoordinates(0),
        supportCoordinates(1),
        supportCoordinates(2)
    ).center.transpose();
    const Handle boundSupport = deformation::makeFixedSupportHandle(
        supportDomain,
        contact
    );
    acceptedInvalid = acceptedInvalid ||
        boundSupport.constrainedNodes.size() != 8U;

    if (acceptedInvalid)
    {
        return fail(scenario, "an invalid deformation input was accepted");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}

} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail(
            "unknown",
            "usage: deformation-regression <scenario>"
        );
    }

    const std::string scenario = argv[1];

    if (scenario == "bbw-path")
    {
        return runBbwPath(scenario);
    }
    if (scenario == "bbw-active-bound")
    {
        return runBbwActiveBound(scenario);
    }
    if (scenario == "bbw-contract")
    {
        return runBbwContract(scenario);
    }
    if (scenario == "handles")
    {
        return runHandles(scenario);
    }
    if (scenario == "lbs")
    {
        return runLbs(scenario);
    }
    if (scenario == "domain-trilinear")
    {
        return runDomainTrilinear(scenario);
    }
    if (scenario == "geometry-sync")
    {
        return runGeometrySync(scenario);
    }
    if (scenario == "invalid")
    {
        return runInvalid(scenario);
    }

    return fail(scenario, "unknown scenario");
}
