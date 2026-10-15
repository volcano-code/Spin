#pragma once
#include "spin2014/SourceSurfaceOctree.h"
#include <Eigen/SparseCore>
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourcePartitionKey {
    int level=0;
    std::array<std::uint32_t,3> coordinates{};
    bool operator<(const SourcePartitionKey& b) const;
    bool operator==(const SourcePartitionKey& b) const;
};
struct SourceAdaptiveLeaf : SourcePartitionKey { physics::RawMassMoments10 moments; };
struct SourcePartitionChange {
    SourcePartitionKey parent;
    std::vector<SourcePartitionKey> children;
};
struct SourcePartitionTransition {
    std::string beforeHash,afterHash;
    std::vector<SourcePartitionChange> changes;
    physics::RawMassMoments10 removedBefore,removedAfter;
};
/** Separate mutable optimization partition of an immutable geometry oracle.
 * Unit face adjacency refers to bounding cells, not physical connectivity.
 * Exact common beta is required for near-bound merging (IMPLEMENTATION_DEFAULT).
 */
class SourceAdaptiveDomain {
public:
    SourceAdaptiveDomain(std::shared_ptr<const SourceSurfaceOctree> source,int initialLevel,
        int maximumOptimizationLevel,std::size_t maximumLeafCount=1000000);
    // Restore only a partition, not an optimization/checkpoint result. Keys
    // must already be canonical (strictly sorted), preserving external beta
    // correspondence. Moments/adjacency are rederived from the immutable source;
    // no saved numerical moments are accepted. Hash is integrity, not authority.
    static SourceAdaptiveDomain restorePartition(
        std::shared_ptr<const SourceSurfaceOctree> source,int initialLevel,
        int maximumOptimizationLevel,std::size_t maximumLeafCount,
        const std::vector<SourcePartitionKey>& keys,const std::string& expectedPartitionHash);
    const SourceSurfaceOctree& source() const {return *source_;}
    const std::vector<SourceAdaptiveLeaf>& leaves() const {return leaves_;}
    const std::vector<std::pair<std::size_t,std::size_t>>& edges() const {return edges_;}
    const Eigen::SparseMatrix<double>& uniformLaplacian() const {return laplacian_;}
    const std::string& partitionHash() const {return hash_;}
    std::size_t faceKeyCount() const {return 6*leaves_.size();}
    int maximumOptimizationLevel() const {return maximumLevel_;}
    void validateCoverage() const;
    physics::RawMassMoments10 removedMoments(const Eigen::VectorXd& beta) const;
    SourcePartitionTransition splitFractionalLeaves(Eigen::VectorXd& beta,double epsilon);
    SourcePartitionTransition mergeBinarySiblings(Eigen::VectorXd& beta,double epsilon);
private:
    struct EmptyPartitionTag {};
    SourceAdaptiveDomain(std::shared_ptr<const SourceSurfaceOctree> source,int initialLevel,
        int maximumOptimizationLevel,std::size_t maximumLeafCount,EmptyPartitionTag);
    std::shared_ptr<const SourceSurfaceOctree> source_;
    int initialLevel_,maximumLevel_;
    std::size_t maximumLeaves_;
    std::vector<SourceAdaptiveLeaf> leaves_;
    std::vector<std::pair<std::size_t,std::size_t>> edges_;
    Eigen::SparseMatrix<double> laplacian_;
    std::string hash_;
    void rebuild();
    SourcePartitionTransition commit(SourceAdaptiveDomain& candidate,Eigen::VectorXd& beta,
        Eigen::VectorXd& next,std::vector<SourcePartitionChange> changes);
};
} // namespace spin2014
#endif
