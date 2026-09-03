#pragma once

#include "physics/RawMassMoments10.h"
#include "voxel/VoxelDomain.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>


namespace spin2014
{


struct AdaptiveVoxelCell
{
    int level = 0;
    Eigen::Vector3i coordinate = Eigen::Vector3i::Zero();
    physics::RawMassMoments10 moments;
};


/**
 * Deterministic dyadic partition of a finest-level material lattice.
 *
 * The source lattice is the boundary-aware integration oracle: occupied
 * source cells contribute exact cuboid moments, while empty or protected
 * cells contribute no optimizable mass.  Coarser leaves are exact sums of
 * those finest contributions, matching Section 5.4's boundary subtrees.
 */
class AdaptiveVoxelDomain
{
public:
    AdaptiveVoxelDomain(
        std::shared_ptr<const voxel::VoxelDomain> source,
        int initialLevel,
        int maximumLevel
    );

    const std::shared_ptr<const voxel::VoxelDomain>& source() const noexcept;
    int initialLevel() const noexcept;
    int maximumLevel() const noexcept;
    const std::vector<AdaptiveVoxelCell>& leaves() const noexcept;
    const physics::RawMassMoments10& solidMoments() const noexcept;
    const physics::RawMassMoments10& variableMoments() const noexcept;

    void splitLeaf(std::size_t leafIndex);
    std::size_t splitFractionalLeaves(Eigen::VectorXd& fill, double epsilon);
    std::size_t mergeBinarySiblings(Eigen::VectorXd& fill, double epsilon);
    std::vector<unsigned char> binaryVoidMask(
        const Eigen::VectorXd& fill,
        double threshold = 0.5
    ) const;
    std::shared_ptr<const voxel::VoxelDomain> maximumResolutionMaterialDomain(
        const Eigen::VectorXd& fill,
        double threshold = 0.5
    ) const;

    Eigen::SparseMatrix<double> uniformLaplacian() const;
    bool validatePartition() const;
    std::string partitionHash() const;

private:
    physics::RawMassMoments10 momentsForCell(
        int level,
        const Eigen::Vector3i& coordinate,
        bool carveableOnly
    ) const;
    void sortLeaves();

    std::shared_ptr<const voxel::VoxelDomain> source_;
    int initialLevel_ = 0;
    int maximumLevel_ = 0;
    int rootSide_ = 0;
    physics::RawMassMoments10 solidMoments_;
    physics::RawMassMoments10 variableMoments_;
    std::vector<AdaptiveVoxelCell> leaves_;
};


} // namespace spin2014
