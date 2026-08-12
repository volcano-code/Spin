#include "inner_void_mesh.h"

#include "extraction/VoidMeshExtractor.h"


InnerVoidMesh::InnerVoidMesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& planeV,
    const Eigen::MatrixXi& planeF,
    const Eigen::Vector3d& gravity,
    const Eigen::Vector3d& balancePoint
)
    :
    voxelGrid(V, F),
    carvingSolver(
        voxelGrid,
        V,
        F,
        planeV,
        planeF,
        gravity,
        balancePoint
    )
{
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


bool InnerVoidMesh::isOptimized()
{
    return carvingSolver.isOptimized();
}


void InnerVoidMesh::carveInnerMesh()
{
    carvingSolver.carveAll();
}
