#include "deformation/DeformationObjective.h"
#include "deformation/DeformationOptimizer.h"
#include "deformation/DeformationSystem.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "deformation/LaplacianEditingEnergy.h"
#include "deformation/LinearBlendSkinning.h"
#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"
#include "physics/SurfaceMassDifferential.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{

using deformation::DeformationEnergyEvaluation;
using deformation::DeformationObjective;
using deformation::DeformationSystem;
using deformation::GridDeformationDomain;
using deformation::Handle;
using deformation::HandleRole;
using deformation::HandleSet;
using deformation::HandleVector;
using deformation::LaplacianEditingEnergy;


constexpr double FINITE_DIFFERENCE_STEP = 1.0e-6;
constexpr double GRADIENT_RELATIVE_TOLERANCE = 1.0e-5;


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


double relativeVectorError(
    const Eigen::VectorXd& analytic,
    const Eigen::VectorXd& numerical
)
{
    if (
        analytic.size() != numerical.size() ||
        !analytic.allFinite() ||
        !numerical.allFinite()
    )
    {
        return std::numeric_limits<double>::infinity();
    }

    const double denominator = std::max({
        1.0e-12,
        analytic.norm(),
        numerical.norm()
    });
    return (analytic - numerical).norm() / denominator;
}


Eigen::VectorXd flattenVertexRows(const Eigen::MatrixXd& vertices)
{
    Eigen::VectorXd flattened(3 * vertices.rows());
    for (Eigen::Index vertex = 0; vertex < vertices.rows(); ++vertex)
    {
        flattened.segment<3>(3 * vertex) =
            vertices.row(vertex).transpose();
    }
    return flattened;
}


void makeBox(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        minimum(0), minimum(1), minimum(2),
        maximum(0), minimum(1), minimum(2),
        maximum(0), maximum(1), minimum(2),
        minimum(0), maximum(1), minimum(2),
        minimum(0), minimum(1), maximum(2),
        maximum(0), minimum(1), maximum(2),
        maximum(0), maximum(1), maximum(2),
        minimum(0), maximum(1), maximum(2);

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


void reverseOrientation(Eigen::MatrixXi& faces)
{
    faces.col(1).swap(faces.col(2));
}


Eigen::VectorXd finiteDifferenceVertexGradient(
    const Eigen::MatrixXd& vertices,
    const std::function<double(const Eigen::MatrixXd&)>& energy
)
{
    Eigen::VectorXd gradient(3 * vertices.rows());

    for (Eigen::Index vertex = 0; vertex < vertices.rows(); ++vertex)
    {
        for (Eigen::Index coordinate = 0; coordinate < 3; ++coordinate)
        {
            Eigen::MatrixXd plus = vertices;
            Eigen::MatrixXd minus = vertices;
            plus(vertex, coordinate) += FINITE_DIFFERENCE_STEP;
            minus(vertex, coordinate) -= FINITE_DIFFERENCE_STEP;

            gradient(3 * vertex + coordinate) =
                (energy(plus) - energy(minus)) /
                (2.0 * FINITE_DIFFERENCE_STEP);
        }
    }

    return gradient;
}


int runLaplacianAnalytic(const std::string& scenario)
{
    Eigen::MatrixXd rest;
    Eigen::MatrixXi faces;
    makeBox(
        Eigen::Vector3d(-0.7, -0.4, -0.2),
        Eigen::Vector3d(0.9, 0.8, 1.1),
        rest,
        faces
    );

    const LaplacianEditingEnergy energy(rest, faces);
    const auto& residual = energy.residualOperator();
    const Eigen::SparseMatrix<double> quadratic = energy.quadraticMatrix();
    const Eigen::MatrixXd quadraticDense(quadratic);

    Eigen::MatrixXd translated = rest;
    translated.rowwise() += Eigen::RowVector3d(2.3, -1.7, 0.8);

    Eigen::MatrixXd similarity = 1.35 * rest;
    similarity.rowwise() += Eigen::RowVector3d(-0.6, 0.4, 1.9);

    Eigen::MatrixXd sheared = rest;
    sheared.col(0) += 0.37 * rest.col(1);
    sheared(6, 2) += 0.23;

    const double restEnergy = energy.energy(rest);
    const double translationEnergy = energy.energy(translated);
    const double similarityEnergy = energy.energy(similarity);
    const double shearEnergy = energy.energy(sheared);

    const Eigen::VectorXd analytic = energy.vertexGradient(sheared);
    const Eigen::VectorXd numerical = finiteDifferenceVertexGradient(
        sheared,
        [&](const Eigen::MatrixXd& value)
        {
            return energy.energy(value);
        }
    );
    const double gradientError = relativeVectorError(analytic, numerical);

    const Eigen::VectorXd probe =
        Eigen::VectorXd::LinSpaced(3 * rest.rows(), -0.8, 1.1);
    const double quadraticNorm = (residual * probe).squaredNorm();
    const double quadraticForm = probe.dot(quadratic * probe);

    if (
        residual.rows() != 3 * rest.rows() ||
        residual.cols() != 3 * rest.rows() ||
        quadratic.rows() != 3 * rest.rows() ||
        !quadraticDense.isApprox(quadraticDense.transpose(), 1.0e-12) ||
        !std::isfinite(quadraticNorm) ||
        !approximatelyEqual(quadraticNorm, quadraticForm, 1.0e-12) ||
        quadraticForm < -1.0e-12 ||
        std::abs(restEnergy) > 1.0e-18 ||
        std::abs(translationEnergy) > 1.0e-18 ||
        std::abs(similarityEnergy) > 1.0e-18 ||
        !std::isfinite(shearEnergy) ||
        shearEnergy <= 1.0e-8 ||
        gradientError > GRADIENT_RELATIVE_TOLERANCE
    )
    {
        return fail(
            scenario,
            "implicit-similarity Laplacian energy/gradient contract failed"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " shear_energy=" << shearEnergy
        << " gradient_relative_error=" << gradientError << '\n';
    return EXIT_SUCCESS;
}


MassProperties combinedMassProperties(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double density
)
{
    return
        computeSurfaceMassProperties(outerVertices, outerFaces, density) +
        computeSurfaceMassProperties(innerVertices, innerFaces, density);
}


Eigen::Vector3d requireCenterOfMass(const MassProperties& properties)
{
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    if (!properties.tryCenterOfMass(center, 0.0))
    {
        throw std::runtime_error("finite-difference fixture lost positive mass");
    }
    return center;
}


int runMassGradient(const std::string& scenario)
{
    constexpr double density = 1.7;

    Eigen::MatrixXd outer;
    Eigen::MatrixXi outerFaces;
    makeBox(
        Eigen::Vector3d(-1.0, -0.8, -0.5),
        Eigen::Vector3d(1.3, 1.1, 1.6),
        outer,
        outerFaces
    );

    Eigen::MatrixXd inner;
    Eigen::MatrixXi innerFaces;
    makeBox(
        Eigen::Vector3d(-0.35, -0.2, 0.05),
        Eigen::Vector3d(0.55, 0.45, 0.9),
        inner,
        innerFaces
    );
    reverseOrientation(innerFaces);

    const SurfaceMassDifferential outerDifferential =
        computeSurfaceMassDifferential(
            outer,
            outerFaces,
            density
        );
    const MassProperties outerOracle =
        computeSurfaceMassProperties(outer, outerFaces, density);

    const Eigen::VectorXd numericalMass = finiteDifferenceVertexGradient(
        outer,
        [&](const Eigen::MatrixXd& value)
        {
            return computeSurfaceMassProperties(
                value,
                outerFaces,
                density
            ).mass;
        }
    );

    Eigen::MatrixXd numericalFirstMoment(3, 3 * outer.rows());
    for (Eigen::Index component = 0; component < 3; ++component)
    {
        numericalFirstMoment.row(component) =
            finiteDifferenceVertexGradient(
                outer,
                [&](const Eigen::MatrixXd& value)
                {
                    return computeSurfaceMassProperties(
                        value,
                        outerFaces,
                        density
                    ).firstMoment(component);
                }
            ).transpose();
    }

    const CenterOfMassDifferential combined =
        computeCenterOfMassDifferential(
            outer,
            outerFaces,
            inner,
            innerFaces,
            density
        );
    const MassProperties combinedOracle = combinedMassProperties(
        outer,
        outerFaces,
        inner,
        innerFaces,
        density
    );
    const Eigen::Vector3d centerOracle = requireCenterOfMass(combinedOracle);

    Eigen::MatrixXd numericalOuterCenter(3, 3 * outer.rows());
    for (Eigen::Index component = 0; component < 3; ++component)
    {
        numericalOuterCenter.row(component) =
            finiteDifferenceVertexGradient(
                outer,
                [&](const Eigen::MatrixXd& value)
                {
                    return requireCenterOfMass(combinedMassProperties(
                        value,
                        outerFaces,
                        inner,
                        innerFaces,
                        density
                    ))(component);
                }
            ).transpose();
    }

    Eigen::MatrixXd numericalInnerCenter(3, 3 * inner.rows());
    for (Eigen::Index component = 0; component < 3; ++component)
    {
        numericalInnerCenter.row(component) =
            finiteDifferenceVertexGradient(
                inner,
                [&](const Eigen::MatrixXd& value)
                {
                    return requireCenterOfMass(combinedMassProperties(
                        outer,
                        outerFaces,
                        value,
                        innerFaces,
                        density
                    ))(component);
                }
            ).transpose();
    }

    const double massGradientError = relativeVectorError(
        outerDifferential.massGradient,
        numericalMass
    );
    const double momentGradientError = relativeVectorError(
        Eigen::Map<const Eigen::VectorXd>(
            outerDifferential.firstMomentJacobian.data(),
            outerDifferential.firstMomentJacobian.size()
        ),
        Eigen::Map<const Eigen::VectorXd>(
            numericalFirstMoment.data(),
            numericalFirstMoment.size()
        )
    );
    const double outerCenterError = relativeVectorError(
        Eigen::Map<const Eigen::VectorXd>(
            combined.outerJacobian.data(),
            combined.outerJacobian.size()
        ),
        Eigen::Map<const Eigen::VectorXd>(
            numericalOuterCenter.data(),
            numericalOuterCenter.size()
        )
    );
    const double innerCenterError = relativeVectorError(
        Eigen::Map<const Eigen::VectorXd>(
            combined.innerJacobian.data(),
            combined.innerJacobian.size()
        ),
        Eigen::Map<const Eigen::VectorXd>(
            numericalInnerCenter.data(),
            numericalInnerCenter.size()
        )
    );

    bool reversedOuterRejected = false;
    try
    {
        Eigen::MatrixXi reversedOuterFaces = outerFaces;
        reverseOrientation(reversedOuterFaces);
        static_cast<void>(computeCenterOfMassDifferential(
            outer,
            reversedOuterFaces,
            inner,
            innerFaces,
            density
        ));
    }
    catch (const std::domain_error&)
    {
        reversedOuterRejected = true;
    }

    bool outwardCavityRejected = false;
    try
    {
        Eigen::MatrixXi outwardInnerFaces = innerFaces;
        reverseOrientation(outwardInnerFaces);
        static_cast<void>(computeCenterOfMassDifferential(
            outer,
            outerFaces,
            inner,
            outwardInnerFaces,
            density
        ));
    }
    catch (const std::domain_error&)
    {
        outwardCavityRejected = true;
    }

    if (
        !approximatelyEqual(
            outerDifferential.properties.mass,
            outerOracle.mass,
            1.0e-12
        ) ||
        !outerDifferential.properties.firstMoment.isApprox(
            outerOracle.firstMoment,
            1.0e-12
        ) ||
        !approximatelyEqual(
            combined.properties.mass,
            combinedOracle.mass,
            1.0e-12
        ) ||
        !combined.centerOfMass.isApprox(centerOracle, 1.0e-12) ||
        massGradientError > GRADIENT_RELATIVE_TOLERANCE ||
        momentGradientError > GRADIENT_RELATIVE_TOLERANCE ||
        outerCenterError > GRADIENT_RELATIVE_TOLERANCE ||
        innerCenterError > GRADIENT_RELATIVE_TOLERANCE ||
        !reversedOuterRejected ||
        !outwardCavityRejected
    )
    {
        return fail(
            scenario,
            "surface mass/first-moment/center-of-mass differential failed"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " mass_error=" << massGradientError
        << " moment_error=" << momentGradientError
        << " center_errors=" << outerCenterError
        << ',' << innerCenterError << '\n';
    return EXIT_SUCCESS;
}


Handle makeSupportHandle(const std::vector<int>& nodes)
{
    Handle handle;
    handle.role = HandleRole::Support;
    handle.locked = true;
    handle.constrainedNodes = nodes;
    return handle;
}


Handle makeUserHandle(
    const std::vector<int>& nodes,
    const Eigen::Vector3d& pivot
)
{
    Handle handle;
    handle.restPivot = pivot;
    handle.role = HandleRole::User;
    handle.locked = false;
    handle.constrainedNodes = nodes;
    return handle;
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

    Handle support = makeSupportHandle({});
    Handle user = makeUserHandle(
        {},
        Eigen::Vector3d(maximum, 0.5, 0.5)
    );

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


bool carveOneInteriorCell(VoxelGrid& grid)
{
    for (int x = 0; x < grid.dimensions()(0); ++x)
    {
        for (int y = 0; y < grid.dimensions()(1); ++y)
        {
            for (int z = 0; z < grid.dimensions()(2); ++z)
            {
                if (grid(x, y, z).isInner)
                {
                    grid(x, y, z).isFilled = false;
                    return true;
                }
            }
        }
    }
    return false;
}


Eigen::VectorXd finiteDifferenceObjectiveGradient(
    const DeformationObjective& objective,
    const Eigen::VectorXd& parameters,
    double mu
)
{
    Eigen::VectorXd gradient(parameters.size());
    for (Eigen::Index index = 0; index < parameters.size(); ++index)
    {
        Eigen::VectorXd plus = parameters;
        Eigen::VectorXd minus = parameters;
        plus(index) += FINITE_DIFFERENCE_STEP;
        minus(index) -= FINITE_DIFFERENCE_STEP;
        gradient(index) =
            (objective.evaluate(plus, mu).totalEnergy -
             objective.evaluate(minus, mu).totalEnergy) /
            (2.0 * FINITE_DIFFERENCE_STEP);
    }
    return gradient;
}


int runObjectiveGradient(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        vertices,
        faces
    );

    VoxelGrid fixedGrid(vertices, faces, 0.0, 10);
    if (!carveOneInteriorCell(fixedGrid))
    {
        return fail(scenario, "objective fixture has no carved cavity");
    }

    const GridDeformationDomain domain =
        GridDeformationDomain::fromVoxelGrid(fixedGrid);
    const GeometryBinding binding = makeGeometryBinding(domain);
    DeformationSystem system(
        fixedGrid,
        vertices,
        faces,
        binding.handles,
        binding.weights,
        true
    );

    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d(0.18, 0.24, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d(0.22, 0.0, 0.73),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );

    constexpr double lambda = 0.37;
    constexpr double mu = 0.31;
    const DeformationObjective objective(
        system,
        fixedGrid,
        objectives,
        lambda
    );

    Eigen::VectorXd parameters = objective.initialParameters();
    if (parameters.size() != 4)
    {
        return fail(scenario, "unexpected free-handle parameter count");
    }
    parameters << 0.08, -0.04, 0.03, 0.92;

    const DeformationEnergyEvaluation evaluation =
        objective.evaluate(parameters, mu);
    const Eigen::VectorXd numerical = finiteDifferenceObjectiveGradient(
        objective,
        parameters,
        mu
    );
    const double gradientError = relativeVectorError(
        evaluation.gradient,
        numerical
    );

    const HandleSet handles = objective.makeHandles(parameters);
    const auto synchronizedGeometry = system.evaluate(handles, fixedGrid);
    const MassProperties synchronizedMass = combinedMassProperties(
        synchronizedGeometry.outerVertices,
        synchronizedGeometry.outerFaces,
        synchronizedGeometry.innerVertices,
        synchronizedGeometry.innerFaces,
        1.0
    );
    const Eigen::Vector3d synchronizedCenter =
        requireCenterOfMass(synchronizedMass);
    const deformation::FreeParameterLbsModel freeOuter(
        vertices,
        system.outerWeights(),
        system.bindingHandles()
    );
    const Eigen::MatrixXd directOuter =
        deformation::deformLinearBlendSkinning(
            vertices,
            system.outerWeights(),
            handles
        );
    const Eigen::MatrixXd affineOuter = freeOuter.evaluate(parameters);
    const Eigen::VectorXd affineStacked =
        freeOuter.offset() + freeOuter.jacobian() * parameters;

    /* The no-carving state must exercise the explicit empty-cavity path. */
    VoxelGrid solidGrid(vertices, faces, 0.0, 10);
    const GridDeformationDomain solidDomain =
        GridDeformationDomain::fromVoxelGrid(solidGrid);
    const GeometryBinding solidBinding = makeGeometryBinding(solidDomain);
    const DeformationSystem solidSystem(
        solidGrid,
        vertices,
        faces,
        solidBinding.handles,
        solidBinding.weights,
        true
    );
    const DeformationObjective solidObjective(
        solidSystem,
        solidGrid,
        objectives,
        lambda
    );
    const DeformationEnergyEvaluation solidEvaluation =
        solidObjective.evaluate(parameters, mu);
    const Eigen::VectorXd solidNumerical =
        finiteDifferenceObjectiveGradient(
            solidObjective,
            parameters,
            mu
        );
    const auto solidGeometry = solidSystem.evaluate(
        solidObjective.makeHandles(parameters),
        solidGrid
    );
    const MassProperties solidMass = computeSurfaceMassProperties(
        solidGeometry.outerVertices,
        solidGeometry.outerFaces,
        1.0
    );

    if (
        !evaluation.gradient.allFinite() ||
        !std::isfinite(evaluation.totalEnergy) ||
        !std::isfinite(evaluation.centerOfMassEnergy) ||
        !std::isfinite(evaluation.shapeEnergy) ||
        evaluation.mass <= 0.0 ||
        !approximatelyEqual(
            evaluation.mass,
            synchronizedMass.mass,
            1.0e-11
        ) ||
        !evaluation.centerOfMass.isApprox(
            synchronizedCenter,
            1.0e-11
        ) ||
        !evaluation.centerOfMass.allFinite() ||
        !approximatelyEqual(
            evaluation.shapeEnergy,
            lambda * evaluation.unscaledShapeEnergy,
            1.0e-11
        ) ||
        !approximatelyEqual(
            evaluation.totalEnergy,
            (1.0 - mu) * evaluation.centerOfMassEnergy +
                mu * evaluation.shapeEnergy,
            1.0e-11
        ) ||
        gradientError > GRADIENT_RELATIVE_TOLERANCE ||
        !affineOuter.isApprox(directOuter, 1.0e-12) ||
        !affineStacked.isApprox(flattenVertexRows(affineOuter), 1.0e-12) ||
        !freeOuter.evaluate(handles).isApprox(affineOuter, 1.0e-12) ||
        solidGeometry.innerFaces.rows() != 0 ||
        !approximatelyEqual(
            solidEvaluation.mass,
            solidMass.mass,
            1.0e-11
        ) ||
        !solidEvaluation.centerOfMass.isApprox(
            requireCenterOfMass(solidMass),
            1.0e-11
        ) ||
        relativeVectorError(
            solidEvaluation.gradient,
            solidNumerical
        ) > GRADIENT_RELATIVE_TOLERANCE
    )
    {
        return fail(
            scenario,
            "combined fixed-alpha energy/free-handle gradient failed"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " total_energy=" << evaluation.totalEnergy
        << " gradient_relative_error=" << gradientError << '\n';
    return EXIT_SUCCESS;
}


int runLineSearch(const std::string& scenario)
{
    const deformation::PaperLineSearchConfig config;

    Eigen::VectorXd parameters(1);
    Eigen::VectorXd gradient(1);
    parameters << 0.0;
    gradient << 1.0;

    deformation::PaperLineSearchState initialState;
    initialState.step = config.initialStep;
    initialState.mu = 0.20;

    const deformation::PaperLineSearchResult boundaryAccepted =
        deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            initialState,
            config,
            [](const Eigen::VectorXd& value, double)
            {
                return value(0) == 0.0 ? 100.0 : 99.0;
            }
        );

    const deformation::PaperLineSearchResult insufficientDecrease =
        deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            initialState,
            config,
            [](const Eigen::VectorXd& value, double)
            {
                return value(0) == 0.0 ? 100.0 : 99.1;
            }
        );

    /* A 2% decrease is accepted by the paper's 1% rule (not beta's 3%). */
    const deformation::PaperLineSearchResult twoPercentAccepted =
        deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            initialState,
            config,
            [](const Eigen::VectorXd& value, double)
            {
                return value(0) == 0.0 ? 100.0 : 98.0;
            }
        );

    Eigen::VectorXd positiveScale(1);
    Eigen::VectorXd scaleGradient(1);
    positiveScale << 1.0;
    scaleGradient << 2.0;
    const deformation::PaperLineSearchResult invalidScaleTrial =
        deformation::paperGradientDescentTrial(
            positiveScale,
            scaleGradient,
            initialState,
            config,
            [](const Eigen::VectorXd& value, double)
            {
                if (value(0) <= 0.0)
                {
                    throw deformation::InfeasibleDeformationState(
                        "fixture models a non-positive handle scale"
                    );
                }
                return 100.0;
            }
        );

    bool internalFailurePropagated = false;
    try
    {
        static_cast<void>(deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            initialState,
            config,
            [](const Eigen::VectorXd& value, double)
            {
                if (value(0) != 0.0)
                {
                    throw std::runtime_error(
                        "fixture models an internal evaluation defect"
                    );
                }
                return 100.0;
            }
        ));
    }
    catch (const std::runtime_error&)
    {
        internalFailurePropagated = true;
    }

    if (
        boundaryAccepted.decision !=
            deformation::PaperTrialDecision::Accepted ||
        !boundaryAccepted.accepted ||
        boundaryAccepted.rolledBack ||
        !approximatelyEqual(boundaryAccepted.energyAfter, 99.0, 1.0e-14) ||
        !boundaryAccepted.parameters.isApprox(
            Eigen::VectorXd::Constant(1, -1.0),
            0.0
        ) ||
        insufficientDecrease.decision !=
            deformation::PaperTrialDecision::StepReduced ||
        insufficientDecrease.accepted ||
        !insufficientDecrease.rolledBack ||
        !insufficientDecrease.parameters.isApprox(parameters, 0.0) ||
        !approximatelyEqual(
            insufficientDecrease.energyAfter,
            insufficientDecrease.energyBefore,
            1.0e-14
        ) ||
        !approximatelyEqual(
            insufficientDecrease.nextState.step,
            0.8,
            1.0e-14
        ) ||
        twoPercentAccepted.decision !=
            deformation::PaperTrialDecision::Accepted ||
        !twoPercentAccepted.accepted ||
        twoPercentAccepted.rolledBack ||
        !approximatelyEqual(twoPercentAccepted.energyAfter, 98.0, 1.0e-14) ||
        invalidScaleTrial.decision !=
            deformation::PaperTrialDecision::StepReduced ||
        invalidScaleTrial.accepted ||
        !invalidScaleTrial.rolledBack ||
        !invalidScaleTrial.parameters.isApprox(positiveScale, 0.0) ||
        !std::isinf(invalidScaleTrial.trialEnergy) ||
        !approximatelyEqual(
            invalidScaleTrial.nextState.step,
            0.8,
            1.0e-14
        ) ||
        !internalFailurePropagated
    )
    {
        return fail(
            scenario,
            "paper 1% acceptance boundary or transactional rollback failed"
        );
    }

    /*
     * The paper specifies a minimum step but not how to land on it. The local
     * deterministic operationalization clamps the geometric sequence to 0.4,
     * tries that value once, and only then reduces mu.
     */
    const double expectedAttempts[] =
    {
        1.0,
        0.8,
        0.64,
        0.512,
        0.4096,
        0.4
    };
    deformation::PaperLineSearchState state = initialState;
    deformation::PaperLineSearchResult failedTrial;

    for (std::size_t attempt = 0;
         attempt < sizeof(expectedAttempts) / sizeof(expectedAttempts[0]);
         ++attempt)
    {
        failedTrial = deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            state,
            config,
            [](const Eigen::VectorXd&, double)
            {
                return 100.0;
            }
        );

        if (
            !approximatelyEqual(
                failedTrial.attemptedStep,
                expectedAttempts[attempt],
                1.0e-14
            ) ||
            failedTrial.accepted ||
            !failedTrial.rolledBack ||
            !failedTrial.parameters.isApprox(parameters, 0.0)
        )
        {
            return fail(
                scenario,
                "minimum-step sequence did not preserve rollback state"
            );
        }

        state = failedTrial.nextState;
    }

    if (
        failedTrial.decision !=
            deformation::PaperTrialDecision::MuReduced ||
        !failedTrial.objectiveReweighted ||
        !approximatelyEqual(state.step, 1.0, 1.0e-14) ||
        !approximatelyEqual(state.mu, 0.15, 1.0e-14)
    )
    {
        return fail(
            scenario,
            "minimum step did not reset step and reduce mu by 0.05"
        );
    }

    state.step = config.minimumStep;
    state.mu = config.muDecrement;
    const deformation::PaperLineSearchResult exhausted =
        deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            state,
            config,
            [](const Eigen::VectorXd&, double)
            {
                return 100.0;
            }
        );

    if (
        exhausted.decision !=
            deformation::PaperTrialDecision::MuExhausted ||
        exhausted.accepted ||
        !exhausted.rolledBack ||
        !exhausted.objectiveReweighted ||
        !approximatelyEqual(exhausted.nextState.mu, 0.0, 0.0) ||
        !approximatelyEqual(
            exhausted.nextState.step,
            config.initialStep,
            0.0
        ) ||
        !exhausted.parameters.isApprox(parameters, 0.0)
    )
    {
        return fail(scenario, "mu exhaustion state-machine contract failed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " acceptance=1%"
        << " failed_attempts="
        << sizeof(expectedAttempts) / sizeof(expectedAttempts[0])
        << " final_mu=" << exhausted.nextState.mu << '\n';
    return EXIT_SUCCESS;
}


std::size_t filledVoxelCount(const VoxelGrid& grid)
{
    std::size_t count = 0;
    for (int x = 0; x < grid.dimensions()(0); ++x)
    {
        for (int y = 0; y < grid.dimensions()(1); ++y)
        {
            for (int z = 0; z < grid.dimensions()(2); ++z)
            {
                if (grid(x, y, z).isFilled)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}


int runOptimizerEndToEnd(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        vertices,
        faces
    );

    VoxelGrid fixedGrid(vertices, faces, 0.0, 10);
    if (!carveOneInteriorCell(fixedGrid))
    {
        return fail(scenario, "optimizer fixture has no carved cavity");
    }
    const std::size_t fillCountBefore = filledVoxelCount(fixedGrid);

    const GridDeformationDomain domain =
        GridDeformationDomain::fromVoxelGrid(fixedGrid);
    const GeometryBinding binding = makeGeometryBinding(domain);
    DeformationSystem system(
        fixedGrid,
        vertices,
        faces,
        binding.handles,
        binding.weights,
        true
    );

    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d(0.20, 0.50, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );

    const DeformationObjective objective(
        system,
        fixedGrid,
        objectives,
        0.08
    );
    const Eigen::VectorXd initialParameters = objective.initialParameters();
    const double initialMu = 0.20;
    const double initialEnergy =
        objective.evaluate(initialParameters, initialMu).totalEnergy;

    deformation::DeformationOptimizerConfig optimizerConfig;
    optimizerConfig.maximumIterations = 100;
    optimizerConfig.gradientTolerance = 1.0e-10;
    const deformation::DeformationOptimizer optimizer(
        objective,
        optimizerConfig
    );
    const deformation::DeformationOptimizationResult result =
        optimizer.optimize(initialParameters, initialMu);
    const double initialEnergyAtFinalMu =
        objective.evaluate(initialParameters, result.mu).totalEnergy;

    bool sawAccepted = false;
    bool traceValid = true;
    for (const deformation::DeformationIterationRecord& record :
         result.iterations)
    {
        if (
            !std::isfinite(record.energyBefore) ||
            !std::isfinite(record.energyAfter) ||
            !std::isfinite(record.gradientNorm) ||
            record.gradientNorm < 0.0
        )
        {
            traceValid = false;
            break;
        }

        const bool shouldReweightObjective =
            record.decision ==
                deformation::PaperTrialDecision::MuReduced ||
            record.decision ==
                deformation::PaperTrialDecision::MuExhausted;
        traceValid = traceValid &&
            record.objectiveReweighted == shouldReweightObjective;

        if (record.accepted)
        {
            sawAccepted = true;
            traceValid = traceValid &&
                record.decision ==
                    deformation::PaperTrialDecision::Accepted &&
                record.energyAfter <=
                    (1.0 -
                     optimizerConfig.lineSearch.requiredRelativeDecrease) *
                    record.energyBefore + 1.0e-13;
        }
        else
        {
            /* At fixed mu a rejected trial is a strict state rollback. */
            if (
                record.decision ==
                    deformation::PaperTrialDecision::StepReduced
            )
            {
                traceValid = traceValid && approximatelyEqual(
                    record.energyAfter,
                    record.energyBefore,
                    1.0e-13
                );
            }
            else if (shouldReweightObjective)
            {
                /*
                 * A mu transition changes the objective itself.  Its two
                 * energy values therefore must not be compared as a descent
                 * pair; only the state transition is asserted here.
                 */
                traceValid = traceValid &&
                    approximatelyEqual(
                        record.stateAfter.step,
                        optimizerConfig.lineSearch.initialStep,
                        0.0
                    ) &&
                    record.stateAfter.mu < record.muBefore;
            }
        }
    }

    const HandleSet finalHandles = objective.makeHandles(result.parameters);
    const auto finalGeometry = system.evaluate(finalHandles, fixedGrid);
    const DeformationEnergyEvaluation freshFinalEvaluation =
        objective.evaluate(result.parameters, result.mu);

    bool supportFixed = true;
    for (const int local : binding.handles[0].constrainedNodes)
    {
        const int global = domain.localToGlobalCorner()(local);
        supportFixed = supportFixed &&
            finalGeometry.gridCorners.row(global).isApprox(
                domain.restNodes().row(local),
                1.0e-11
            );
    }

    if (
        !sawAccepted ||
        !traceValid ||
        result.iterations.empty() ||
        result.status !=
            deformation::DeformationOptimizationStatus::MuExhausted ||
        result.iterations.back().decision !=
            deformation::PaperTrialDecision::MuExhausted ||
        result.iterations.back().accepted ||
        !approximatelyEqual(result.mu, 0.0, 0.0) ||
        result.parameters.size() != initialParameters.size() ||
        !result.parameters.allFinite() ||
        result.parameters(3) <= 0.0 ||
        !result.evaluation.gradient.allFinite() ||
        !std::isfinite(result.evaluation.totalEnergy) ||
        result.evaluation.mass <= 0.0 ||
        !approximatelyEqual(
            result.evaluation.totalEnergy,
            freshFinalEvaluation.totalEnergy,
            1.0e-13
        ) ||
        !result.evaluation.gradient.isApprox(
            freshFinalEvaluation.gradient,
            1.0e-13
        ) ||
        result.evaluation.totalEnergy >
            initialEnergyAtFinalMu + 1.0e-12 ||
        !supportFixed ||
        filledVoxelCount(fixedGrid) != fillCountBefore
    )
    {
        return fail(
            scenario,
            "end-to-end optimizer monotonicity or fixed constraints failed"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " iterations=" << result.iterations.size()
        << " initial_energy=" << initialEnergy
        << " same_mu_initial_energy=" << initialEnergyAtFinalMu
        << " final_energy=" << result.evaluation.totalEnergy
        << " final_mu=" << result.mu << '\n';
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail(
            "unknown",
            "usage: deformation-optimization-regression <scenario>"
        );
    }

    const std::string scenario = argv[1];

    if (scenario == "laplacian")
    {
        return runLaplacianAnalytic(scenario);
    }
    if (scenario == "mass-gradient")
    {
        return runMassGradient(scenario);
    }
    if (scenario == "objective-gradient")
    {
        return runObjectiveGradient(scenario);
    }
    if (scenario == "line-search")
    {
        return runLineSearch(scenario);
    }
    if (scenario == "optimizer-e2e")
    {
        return runOptimizerEndToEnd(scenario);
    }

    return fail(scenario, "unknown scenario");
}
