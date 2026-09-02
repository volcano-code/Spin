#include "api/BalanceModel.h"

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace
{

void cube(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
    vertices.resize(8, 3);
    vertices <<
        0,0,0, 1,0,0, 1,1,0, 0,1,0,
        0,0,1, 1,0,1, 1,1,1, 0,1,1;
    faces.resize(12, 3);
    faces <<
        0,2,1, 0,3,2, 4,5,6, 4,6,7,
        0,1,5, 0,5,4, 3,7,6, 3,6,2,
        0,4,7, 0,7,3, 1,2,6, 1,6,5;
}

}

int main()
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    cube(vertices, faces);

    api::ReleasedBetaArchiveData archive;
    archive.voxelResolution = 3;
    archive.hullDepth = 0;
    archive.cellIdsZFast.resize(27);
    for (int index = 0; index < 27; ++index)
    {
        archive.cellIdsZFast[static_cast<std::size_t>(index)] = index;
    }
    archive.betaNodeWeightsZFast = Eigen::MatrixXd::Ones(64, 1);

    api::BalanceModelConfig modelConfig;
    modelConfig.voxelResolution = 3;
    modelConfig.fixedUnitCubeVoxelGrid = true;
    std::unique_ptr<api::BalanceModel> model =
        api::BalanceModel::fromReleasedBetaArchive(
            vertices, faces, archive, modelConfig
        );
    if (!model || model->config().voxelResolution != 3 ||
        !model->config().fixedUnitCubeVoxelGrid)
    {
        std::cerr << "released-beta factory lost its fixed-domain contract\n";
        return 1;
    }

    const api::BalanceRunConfig beta = api::releasedBetaRunConfig();
    if (
        std::abs(beta.initialMu - 0.75) > 1.0e-15 ||
        std::abs(beta.optimizer.shapeScaleLambda - 20.0) > 1.0e-15 ||
        !beta.optimizer.lineSearch.releasedBetaSemantics ||
        !beta.optimizer.freezeTargetsAfterInitialRefresh ||
        !beta.optimizer.releasedBetaOuterLoop
    )
    {
        std::cerr << "released-beta profile mismatch\n";
        return 1;
    }

    archive.betaNodeWeightsZFast(0, 0) = -1.0;
    try
    {
        static_cast<void>(api::BalanceModel::fromReleasedBetaArchive(
            vertices, faces, archive, modelConfig
        ));
    }
    catch (const std::invalid_argument&)
    {
        return 0;
    }
    std::cerr << "negative archived BBW weight was accepted\n";
    return 1;
}
