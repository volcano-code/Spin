#include "spin2014/SurfaceRawMomentDifferential.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>


namespace
{


Eigen::MatrixXd tetraVertices()
{
    Eigen::MatrixXd vertices(4, 3);
    vertices <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    return vertices;
}


Eigen::MatrixXi tetraFaces()
{
    Eigen::MatrixXi faces(4, 3);
    faces << 1, 2, 3, 0, 2, 1, 0, 1, 3, 0, 3, 2;
    return faces;
}


} // namespace


int main()
{
    const Eigen::MatrixXd vertices = tetraVertices();
    const Eigen::MatrixXi faces = tetraFaces();
    const spin2014::SurfaceRawMomentDifferential value =
        spin2014::differentiateSurfaceRawMoments(vertices, faces, 1.0);
    Eigen::Matrix<double, 10, 1> expected;
    expected <<
        1.0 / 6.0,
        1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0,
        1.0 / 60.0, 1.0 / 60.0, 1.0 / 60.0,
        1.0 / 120.0, 1.0 / 120.0, 1.0 / 120.0;
    if (!value.moments.values().isApprox(expected, 1e-13))
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-surface-differential moments\n";
        return EXIT_FAILURE;
    }

    const double step = 1e-6;
    double maximumError = 0.0;
    for (Eigen::Index vertex = 0; vertex < vertices.rows(); ++vertex)
    {
        for (Eigen::Index axis = 0; axis < 3; ++axis)
        {
            Eigen::MatrixXd plus = vertices;
            Eigen::MatrixXd minus = vertices;
            plus(vertex, axis) += step;
            minus(vertex, axis) -= step;
            const auto plusValue = spin2014::differentiateSurfaceRawMoments(
                plus, faces, 1.0
            );
            const auto minusValue = spin2014::differentiateSurfaceRawMoments(
                minus, faces, 1.0
            );
            const Eigen::Matrix<double, 10, 1> finite =
                (plusValue.moments.values() - minusValue.moments.values()) /
                (2.0 * step);
            maximumError = std::max(
                maximumError,
                (finite - value.jacobian.col(3 * vertex + axis))
                    .lpNorm<Eigen::Infinity>()
            );
        }
    }
    if (maximumError > 2e-10)
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-surface-differential error="
                  << maximumError << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-surface-differential error="
              << maximumError << '\n';
    return EXIT_SUCCESS;
}
