#include "api/BalanceModel.h"
#include "reachability/CarvingReachability.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace
{

void makeCube(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
    vertices.resize(8, 3);
    vertices <<
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1;
    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}

} // namespace

int main()
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeCube(vertices, faces);
    api::BalanceModelConfig config;
    config.voxelResolution = 8;
    api::BalanceModel model(vertices, faces, config);

    Eigen::MatrixXd triangularSupport(3, 3);
    triangularSupport <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0;
    model.addStandingConstraintPolygon(
        triangularSupport,
        Eigen::Vector3d(0.0, -1.0, 0.0),
        5.0 * std::acos(-1.0) / 180.0
    );

    const reachability::CarvingReachabilityResult result =
        model.analyzeCarvingReachability(0);
    if (
        result.status !=
            reachability::ReachabilityStatus::RelaxationReachable ||
        result.variableCellCount <= 0 ||
        !(result.fixedMass > 0.0) ||
        result.maximumPrimalViolation > 1e-9
    )
    {
        std::cerr
            << "[FAIL] scenario=balance-model-reachability status="
            << reachability::toString(result.status)
            << " variable_cells=" << result.variableCellCount
            << " fixed_mass=" << result.fixedMass
            << " violation=" << result.maximumPrimalViolation << '\n';
        return EXIT_FAILURE;
    }
    std::cout
        << "[PASS] scenario=balance-model-reachability status="
        << reachability::toString(result.status)
        << " variable_cells=" << result.variableCellCount << '\n';
    return EXIT_SUCCESS;
}
