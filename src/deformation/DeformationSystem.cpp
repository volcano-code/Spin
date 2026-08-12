#include "deformation/DeformationSystem.h"

#include "deformation/LinearBlendSkinning.h"
#include "extraction/VoidMeshExtractor.h"
#include "voxel/VoxelGrid.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>


namespace deformation
{
namespace
{


void validateFaces(
    const Eigen::MatrixXi& faces,
    Eigen::Index vertexCount,
    const char* name
)
{
    if (
        faces.cols() != 3 ||
        (
            faces.size() > 0 &&
            (faces.minCoeff() < 0 || faces.maxCoeff() >= vertexCount)
        )
    )
    {
        throw std::invalid_argument(
            std::string(name) + " must contain valid triangle indices"
        );
    }
}


void validateOuterMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (
        vertices.rows() == 0 ||
        vertices.cols() != 3 ||
        !vertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "outer rest vertices must be a nonempty finite N-by-3 matrix"
        );
    }

    validateFaces(faces, vertices.rows(), "outer faces");
}


void validateHandleInterpolation(
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double tolerance
)
{
    for (std::size_t handleIndex = 0;
         handleIndex < handles.size();
         ++handleIndex)
    {
        for (const int node : handles[handleIndex].constrainedNodes)
        {
            if (node < 0 || node >= weights.rows())
            {
                throw std::invalid_argument(
                    "handle constraint is outside imported weight rows"
                );
            }

            for (std::size_t column = 0;
                 column < handles.size();
                 ++column)
            {
                const double expected =
                    column == handleIndex ? 1.0 : 0.0;
                if (
                    std::abs(
                        weights(
                            node,
                            static_cast<Eigen::Index>(column)
                        ) - expected
                    ) > tolerance
                )
                {
                    throw std::invalid_argument(
                        "weights do not interpolate a bound handle constraint"
                    );
                }
            }
        }
    }
}


void validateSameGrid(
    const VoxelGrid& rest,
    const VoxelGrid& current
)
{
    if (
        current.dimensions() != rest.dimensions() ||
        current.voxelSize() != rest.voxelSize() ||
        current.minimumWallThickness() != rest.minimumWallThickness() ||
        current.longestAxisResolution() != rest.longestAxisResolution() ||
        current.corners().rows() != rest.corners().rows() ||
        current.corners().cols() != 3 ||
        !current.corners().isApprox(
            rest.corners(),
            64.0 * std::numeric_limits<double>::epsilon()
        )
    )
    {
        throw std::invalid_argument(
            "current voxel fill does not share the deformation rest lattice"
        );
    }

    const Eigen::Vector3i dimensions = rest.dimensions();
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& expected = rest(x, y, z);
                const Voxel& actual = current(x, y, z);

                if (
                    actual.xIdx != expected.xIdx ||
                    actual.yIdx != expected.yIdx ||
                    actual.zIdx != expected.zIdx ||
                    actual.isInner != expected.isInner ||
                    actual.distanceToMesh != expected.distanceToMesh ||
                    actual.cornerIndices != expected.cornerIndices ||
                    !actual.center.isApprox(
                        expected.center,
                        64.0 * std::numeric_limits<double>::epsilon()
                    )
                )
                {
                    throw std::invalid_argument(
                        "only voxel fill flags may change after deformation binding"
                    );
                }
            }
        }
    }
}


void validateSupportPoseInvariant(
    const HandleSet& binding,
    const HandleSet& current
)
{
    for (std::size_t index = 0; index < binding.size(); ++index)
    {
        const Handle& boundHandle = binding[index];
        const Handle& currentHandle = current[index];

        if (
            boundHandle.role == HandleRole::Support &&
            (
                !currentHandle.translation.isZero(0.0) ||
                !currentHandle.fixedRotation.isIdentity(0.0) ||
                currentHandle.uniformScale != 1.0
            )
        )
        {
            throw std::invalid_argument(
                "support handle pose must remain fixed at identity"
            );
        }
    }
}


} // namespace


DeformationSystem::DeformationSystem(
    const VoxelGrid& voxelGrid,
    const Eigen::MatrixXd& outerRestVertices,
    const Eigen::MatrixXi& outerFaces,
    const HandleSet& handles,
    const WeightSolveOptions& options
)
    :
    domain_(GridDeformationDomain::fromVoxelGrid(voxelGrid)),
    restVoxelGrid_(voxelGrid),
    bindingHandles_(handles),
    outerRestVertices_(outerRestVertices),
    outerFaces_(outerFaces),
    nodeWeights_(),
    fullCornerWeights_(),
    outerWeights_(),
    weightDiagnostics_()
{
    validateOuterMesh(outerRestVertices_, outerFaces_);

    const BoundedBiharmonicWeightResult result =
        solveBoundedBiharmonicWeights(
            domain_,
            bindingHandles_,
            options
        );

    weightDiagnostics_ = result.diagnostics;

    if (!weightDiagnostics_.converged)
    {
        throw std::runtime_error(
            "bounded biharmonic solve reached its iteration limit"
        );
    }

    initializeVerifiedWeights(result.weights);
}


DeformationSystem::DeformationSystem(
    const VoxelGrid& voxelGrid,
    const Eigen::MatrixXd& outerRestVertices,
    const Eigen::MatrixXi& outerFaces,
    const HandleSet& handles,
    const Eigen::MatrixXd& nodeWeights,
    bool verifiedWeightImport
)
    :
    domain_(GridDeformationDomain::fromVoxelGrid(voxelGrid)),
    restVoxelGrid_(voxelGrid),
    bindingHandles_(handles),
    outerRestVertices_(outerRestVertices),
    outerFaces_(outerFaces),
    nodeWeights_(),
    fullCornerWeights_(),
    outerWeights_(),
    weightDiagnostics_()
{
    if (!verifiedWeightImport)
    {
        throw std::invalid_argument(
            "imported deformation weights require explicit verification"
        );
    }

    validateOuterMesh(outerRestVertices_, outerFaces_);

    initializeVerifiedWeights(nodeWeights);
    weightDiagnostics_.verifiedImport = true;
}


void DeformationSystem::initializeVerifiedWeights(
    const Eigen::MatrixXd& nodeWeights
)
{
    validateSkinningWeights(
        domain_.restNodes(),
        nodeWeights,
        bindingHandles_.size()
    );

    validateHandleInterpolation(
        nodeWeights,
        bindingHandles_,
        1.0e-10
    );

    nodeWeights_ = nodeWeights;

    /*
     * The paper's voxel graph contains only object cells. VoxelGrid retains a
     * full bounding-box lattice for legacy global corner indices, so extend
     * weights to compatibility-only outside corners by the nearest active
     * node. Active object nodes and every cavity vertex retain exact BBW.
     */
    fullCornerWeights_.resize(
        domain_.fullRestCorners().rows(),
        nodeWeights_.cols()
    );

    for (Eigen::Index global = 0;
         global < fullCornerWeights_.rows();
         ++global)
    {
        const int local = domain_.globalToLocalCorner(
            static_cast<int>(global)
        );
        const int source = local >= 0
            ? local
            : domain_.nearestActiveNode(
                domain_.fullRestCorners().row(global).transpose()
            );
        fullCornerWeights_.row(global) = nodeWeights_.row(source);
    }

    validateSkinningWeights(
        domain_.fullRestCorners(),
        fullCornerWeights_,
        bindingHandles_.size()
    );

    outerWeights_ = domain_.interpolateWeights(
        outerRestVertices_,
        nodeWeights_
    );

    validateSkinningWeights(
        outerRestVertices_,
        outerWeights_,
        bindingHandles_.size()
    );
}


const GridDeformationDomain& DeformationSystem::domain() const
{
    return domain_;
}


const HandleSet& DeformationSystem::bindingHandles() const
{
    return bindingHandles_;
}


const Eigen::MatrixXd& DeformationSystem::nodeWeights() const
{
    return nodeWeights_;
}


const Eigen::MatrixXd& DeformationSystem::fullCornerWeights() const
{
    return fullCornerWeights_;
}


const Eigen::MatrixXd& DeformationSystem::outerWeights() const
{
    return outerWeights_;
}


const WeightSolveDiagnostics&
DeformationSystem::weightDiagnostics() const
{
    return weightDiagnostics_;
}


DeformedGeometry DeformationSystem::evaluate(
    const HandleSet& currentHandles,
    const VoxelGrid& currentVoxelFill
) const
{
    if (!bindingHandles_.sameBindingAs(currentHandles))
    {
        throw std::invalid_argument(
            "current handles do not match the precomputed deformation binding"
        );
    }

    validateSupportPoseInvariant(
        bindingHandles_,
        currentHandles
    );

    validateSameGrid(restVoxelGrid_, currentVoxelFill);

    DeformedGeometry geometry;
    geometry.outerFaces = outerFaces_;
    geometry.outerVertices = deformLinearBlendSkinning(
        outerRestVertices_,
        outerWeights_,
        currentHandles
    );

    geometry.gridCorners = deformLinearBlendSkinning(
        domain_.fullRestCorners(),
        fullCornerWeights_,
        currentHandles
    );

    geometry.voxelCenters = domain_.deformedCellCenters(
        geometry.gridCorners
    );

    VoidMeshExtractor::extractFaces(
        currentVoxelFill,
        geometry.innerFaces
    );

    geometry.innerVertices = geometry.gridCorners;
    validateFaces(
        geometry.innerFaces,
        geometry.innerVertices.rows(),
        "inner faces"
    );

    for (Eigen::Index index = 0;
         index < geometry.innerFaces.size();
         ++index)
    {
        const int global = geometry.innerFaces(index);

        if (domain_.globalToLocalCorner(global) < 0)
        {
            throw std::logic_error(
                "cavity connectivity references an inactive deformation node"
            );
        }
    }

    return geometry;
}


} // namespace deformation
