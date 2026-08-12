#include "inner_void_mesh.h"

#include "extraction/VoidMeshExtractor.h"


InnerVoidMesh::InnerVoidMesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const BalanceObjectives& objectives
)
    :
    voxelGrid(V, F),
    carvingSolver(
        voxelGrid,
        V,
        F,
        objectives
    )
{
}


InnerVoidMesh::InnerVoidMesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& planeV,
    const Eigen::MatrixXi& planeF,
    const Eigen::Vector3d& gravity,
    const Eigen::Vector3d& balancePoint
)
    : InnerVoidMesh(
        V,
        F,
        BalanceObjectives{
            BalanceObjective(
                balancePoint,
                gravity
            )
        }
    )
{
    static_cast<void>(planeV);
    static_cast<void>(planeF);
}


void InnerVoidMesh::convertToMesh(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
)
{
    VoidMeshExtractor::extract(
        voxelGrid,
        V,
        F
    );
}


void InnerVoidMesh::carveInnerMeshStep(
    Eigen::MatrixXd& carvePlaneV,
    Eigen::MatrixXi& carvePlaneF,
    Eigen::Vector3d& com
)
{
    carvingSolver.carveStep(
        carvePlaneV,
        carvePlaneF,
        com
    );
}


double InnerVoidMesh::getCoMDistance()
{
    return carvingSolver.getCoMDistance();
}


const BalanceObjectives& InnerVoidMesh::balanceObjectives() const
{
    return carvingSolver.balanceObjectives();
}


bool InnerVoidMesh::isOptimized()
{
    return carvingSolver.isOptimized();
}


void InnerVoidMesh::carveInnerMesh()
{
    carvingSolver.carveAll();
}
