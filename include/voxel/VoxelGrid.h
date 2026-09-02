#pragma once

#include <Eigen/Core>

#include "list3d.h"
#include "voxel/Voxel.h"

#include <memory>
#include <vector>

namespace voxel
{
class VoxelDomain;
}

/**
 * Regular volumetric representation of an input triangle mesh.
 *
 * This class owns the voxel states, the shared corner positions, the grid
 * dimensions, and the voxel size in model coordinates. It deliberately does not contain
 * any carving policy.
 */
class VoxelGrid
{
public:
    VoxelGrid();

    VoxelGrid(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        double minimumWallThickness = 0.0,
        int longestAxisResolution = 50,
        bool fixedUnitCubeDomain = false
    );

    /** Build or rebuild the voxel representation from a triangle mesh. */
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        double minimumWallThickness = 0.0,
        int longestAxisResolution = 50,
        bool fixedUnitCubeDomain = false
    );

    /**
     * Construct the exact occupied-cell lattice serialized by the released
     * beta's ``.vox`` cache.
     *
     * ``cellIdsZFast`` contains N^3 entries in the beta writer's nested
     * x/y/z loop order (z changes fastest). Empty cells are -1 and occupied
     * identifiers must be a permutation of [0, occupied_count). This named
     * factory is compatibility-only; normal mesh constructors continue to
     * voxelize and audit their input under the modern contract.
     */
    static VoxelGrid fromReleasedBetaOccupancy(
        int resolution,
        const std::vector<int>& cellIdsZFast,
        int hullDepth = 2
    );

    /**
     * Build a legacy-compatible facade over an immutable material domain.
     *
     * The domain mask is the baseline material state. Enclosed complement
     * cells become synchronized initial cavities; ``refillAll()`` restores
     * this baseline instead of incorrectly filling those cavities.
     */
    static VoxelGrid fromMaterialOccupancy(
        const std::shared_ptr<const voxel::VoxelDomain>& domain,
        double minimumWallThickness = 0.0,
        int protectedBoundaryLayers = 1
    );

    Voxel& operator()(
        int x,
        int y,
        int z
    );

    const Voxel& operator()(
        int x,
        int y,
        int z
    ) const;

    const Eigen::MatrixXd& corners() const;

    const Eigen::Vector3i& dimensions() const;

    const Eigen::Vector3d& voxelSize() const;

    /** Requested model-space shell thickness for the fixed outer surface. */
    double minimumWallThickness() const;

    int longestAxisResolution() const;

    /** True for the opt-in released-beta [0,1]^3, N^3 lattice. */
    bool usesFixedUnitCubeDomain() const;

    /**
     * Restore every lattice cell to the filled state used at the beginning
     * of a paper Section 4 carving block.
     *
     * Returns the number of flags changed from empty to filled. Geometry,
     * signed-distance metadata and carve eligibility are left untouched.
     */
    int refillAll();

private:
    List3d<Voxel> voxels_;

    Eigen::MatrixXd corners_;

    Eigen::Vector3i dimensions_;

    Eigen::Vector3d voxelSize_;

    double minimumWallThickness_;

    int longestAxisResolution_;

    bool fixedUnitCubeDomain_;

    // z-fast legacy facade order: (x * ny + y) * nz + z.
    std::vector<unsigned char> baselineFilled_;
};
