#include "spin2014/CageSpinDeformation.h"
#include "spin2014/CageSpinDeformationSolver.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iostream>


namespace
{


spin::TriangleMesh tetra(double scale, const Eigen::Vector3d& offset)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(4, 3);
    mesh.vertices << 0,0,0, scale,0,0, 0,scale,0, 0,0,scale;
    mesh.vertices.rowwise() += offset.transpose();
    mesh.faces.resize(4, 3);
    mesh.faces << 1,2,3, 0,2,1, 0,1,3, 0,3,2;
    return mesh;
}


Eigen::MatrixXd weights(Eigen::Index rows)
{
    Eigen::MatrixXd result(rows, 2);
    for (Eigen::Index row = 0; row < rows; ++row)
    {
        result(row, 0) = row < 2 ? 0.8 : 0.2;
        result(row, 1) = 1.0 - result(row, 0);
    }
    return result;
}


spin2014::CageSpinDeformationRequest request()
{
    spin2014::CageSpinDeformationRequest result;
    result.outerSurface = tetra(3.0, Eigen::Vector3d(-0.7, -1.1, 0.0));
    result.voidCells = {tetra(0.3, Eigen::Vector3d(-0.1, -0.2, 0.4))};
    result.outerWeights = weights(4);
    result.voidCellWeights = {weights(4)};
    deformation::Handle first;
    first.restPivot = Eigen::Vector3d(-0.7, -1.1, 0.0);
    first.constrainedNodes = {0};
    deformation::Handle second;
    second.restPivot = Eigen::Vector3d(2.3, 1.9, 3.0);
    second.constrainedNodes = {1};
    result.handles = {first, second};
    result.problem = spin2014::SpinIt2014MomentProblem::staticBalance(0.0);
    result.fillLaplacian = Eigen::MatrixXd::Zero(1, 1);
    result.biharmonicWeight = 0.0;
    result.displacementWeight = 0.02;
    return result;
}


} // namespace


int main()
{
    const spin2014::CageSpinDeformationFunctional functional(request());
    Eigen::VectorXd beta(1);
    beta << 0.35;
    spin2014::CageSpinDeformationSolverSettings settings;
    settings.maximumIterations = 180;
    settings.equalityTolerance = 1e-8;
    settings.projectedGradientTolerance = 1e-6;
    settings.maximumTranslation = 0.6;
    settings.minimumScale = 0.75;
    settings.maximumScale = 1.25;
    const auto result = spin2014::solveCageSpinDeformation(
        functional, beta, settings
    );
    const auto repeated = spin2014::solveCageSpinDeformation(
        functional, beta, settings
    );
    if (
        result.status != spin2014::CageSpinStatus::Converged ||
        result.evaluation.constraints.lpNorm<Eigen::Infinity>() > 1e-7 ||
        !result.silhouetteChanged || !result.geometryValid ||
        result.traceHash.size() != 64U || result.outputHash.size() != 64U ||
        result.outputHash != repeated.outputHash ||
        !result.handleParameters.isApprox(repeated.handleParameters, 0.0)
    )
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-cage-solver status="
                  << spin2014::cageSpinStatusToken(result.status)
                  << " diagnostic=" << result.diagnostic
                  << " residual="
                  << result.evaluation.constraints.lpNorm<Eigen::Infinity>()
                  << " objective=" << result.evaluation.objective << '\n';
        for (const auto& entry : result.trace)
        {
            std::cerr << "iter=" << entry.iteration
                      << " eq=" << entry.equalityInfinityNorm
                      << " pg=" << entry.projectedGradientInfinityNorm
                      << " step=" << entry.acceptedStep
                      << " geom=" << entry.geometryAccepted << '\n';
        }
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-cage-solver iterations="
              << result.trace.size() << " hash=" << result.outputHash << '\n';
    return EXIT_SUCCESS;
}
