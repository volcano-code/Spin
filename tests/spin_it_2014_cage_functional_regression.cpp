#include "spin2014/CageSpinDeformation.h"

#include "deformation/HandleSet.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>


namespace
{


spin::TriangleMesh tetra(double scale, const Eigen::Vector3d& offset)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(4, 3);
    mesh.vertices <<
        0.0, 0.0, 0.0,
        scale, 0.0, 0.0,
        0.0, scale, 0.0,
        0.0, 0.0, scale;
    mesh.vertices.rowwise() += offset.transpose();
    mesh.faces.resize(4, 3);
    mesh.faces << 1, 2, 3, 0, 2, 1, 0, 1, 3, 0, 3, 2;
    return mesh;
}


deformation::HandleVector handles()
{
    deformation::Handle first;
    first.restPivot = Eigen::Vector3d::Zero();
    first.uniformScale = 1.0;
    first.constrainedNodes = {0};
    deformation::Handle second;
    second.restPivot = Eigen::Vector3d(2.0, 2.0, 2.0);
    second.uniformScale = 1.0;
    second.constrainedNodes = {1};
    return {first, second};
}


Eigen::MatrixXd weights(const Eigen::MatrixXd& vertices)
{
    Eigen::MatrixXd result(vertices.rows(), 2);
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        const double blend = std::clamp(
            0.15 + 0.1 * vertices.row(row).sum(), 0.05, 0.95
        );
        result(row, 0) = 1.0 - blend;
        result(row, 1) = blend;
    }
    return result;
}


} // namespace


int main()
{
    spin2014::CageSpinDeformationRequest request;
    request.outerSurface = tetra(3.0, Eigen::Vector3d(-1.0, -1.0, 0.0));
    request.voidCells = {tetra(0.3, Eigen::Vector3d(-0.1, -0.1, 0.5))};
    request.outerWeights = weights(request.outerSurface.vertices);
    request.voidCellWeights = {weights(request.voidCells.front().vertices)};
    request.handles = handles();
    request.problem = spin2014::SpinIt2014MomentProblem::staticBalance(0.7);
    request.fillLaplacian = Eigen::MatrixXd::Zero(1, 1);
    request.fillLaplacianWeight = 0.1;
    request.biharmonicWeight = 0.03;
    request.displacementWeight = 0.04;
    const spin2014::CageSpinDeformationFunctional functional(request);

    Eigen::VectorXd beta(1);
    beta << 0.37;
    Eigen::VectorXd q = functional.initialHandleParameters();
    q(0) = 0.02;
    q(5) = -0.03;
    q(3) = 0.97;
    q(7) = 1.04;
    const auto value = functional.evaluate(beta, q);
    const double step = 1e-6;
    double gradientError = 0.0;
    double jacobianError = 0.0;

    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        Eigen::VectorXd plus = beta;
        Eigen::VectorXd minus = beta;
        plus(index) += step;
        minus(index) -= step;
        const auto p = functional.evaluate(plus, q);
        const auto m = functional.evaluate(minus, q);
        gradientError = std::max(
            gradientError,
            std::abs((p.objective - m.objective) / (2.0 * step) -
                     value.fillGradient(index))
        );
        jacobianError = std::max(
            jacobianError,
            ((p.constraints - m.constraints) / (2.0 * step) -
             value.fillConstraintJacobian.col(index)).lpNorm<Eigen::Infinity>()
        );
    }
    for (Eigen::Index index = 0; index < q.size(); ++index)
    {
        Eigen::VectorXd plus = q;
        Eigen::VectorXd minus = q;
        plus(index) += step;
        minus(index) -= step;
        const auto p = functional.evaluate(beta, plus);
        const auto m = functional.evaluate(beta, minus);
        gradientError = std::max(
            gradientError,
            std::abs((p.objective - m.objective) / (2.0 * step) -
                     value.handleGradient(index))
        );
        jacobianError = std::max(
            jacobianError,
            ((p.constraints - m.constraints) / (2.0 * step) -
             value.handleConstraintJacobian.col(index)).lpNorm<Eigen::Infinity>()
        );
    }
    if (
        gradientError > 2e-7 || jacobianError > 2e-8 ||
        value.shapeEnergy <= 0.0 || value.displacementEnergy <= 0.0
    )
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-cage-functional gradient="
                  << gradientError << " jacobian=" << jacobianError
                  << " shape=" << value.shapeEnergy << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-cage-functional gradient="
              << gradientError << " jacobian=" << jacobianError << '\n';
    return EXIT_SUCCESS;
}
