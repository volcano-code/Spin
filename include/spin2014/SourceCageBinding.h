#pragma once
#include "spin2014/SourceAdaptiveDomain.h"
#include "spin2014/SourceCageBbw.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageBindingConfig {
    std::size_t maximumPrimitives=1000000;
    std::size_t maximumPoints=2000000;
    // Per returned chunk, never a dense point-by-handle field requirement.
    std::size_t maximumWeightValues=1000000;
};
/** Immutable REST geometry/actual final BBW binding, independently of Faster.
 * Snapshots the adaptive partition and retains immutable source/BBW objects.
 * Requires byte-identical outer arrays and root, actual BBW success and final
 * primal validity. Independent normalization retains its upstream assessment;
 * acceptance here does not assert coupled KKT. No supplied-weight constructor.
 * Actual included Full cuboids carry their one adaptive beta owner. Points are
 * original outer vertices followed by unique material corners in lexicographic
 * maximumDepth integer-lattice order. Original surface vertices are never
 * deduplicated into corners. Cuboid corners use x bit first; no grid snapping.
 * This is not a conforming deformed material boundary, joint solver, paper
 * regularizer, geometry seal, physical certificate or original-case acceptance.
 * Resource exhaustion (including invalid budgets) throws RESOURCE_LIMIT
 * length_error; invalid/mismatched inputs throw invalid_argument.
 */
class SourceCageBinding {
public:
    using Key=SourceCageGrid::Key;
    struct Primitive {
        SourcePartitionKey sourceKey;
        SourceSurfaceCellBounds bounds;
        std::size_t betaIndex=0;
        std::array<std::size_t,8> pointIndices{};
    };
    SourceCageBinding(const SourceAdaptiveDomain&,const SourceCageBbw&,
                      const SourceCageBindingConfig&);
    const SourceAdaptiveDomain& domain() const;
    const SourceSurfaceOctree& source() const;
    const SourceCageBbw& bbw() const;
    const SourceCageBindingConfig& config() const;
    const Eigen::MatrixXd& points() const;
    const std::vector<Key>& materialCornerKeys() const;
    const std::vector<Primitive>& primitives() const;
    const std::vector<SourceCageGrid::Embedding>& embeddings() const;
    std::size_t outerPointCount() const;
    std::size_t handleCount() const;
    /** count x ALL handles, computed from cached eight-node embeddings and
     * upstream finalWeights. No normalization, clamping or handle selection.
     * Empty count is valid for begin<=points.rows(). Invalid/overflowing ranges
     * throw out_of_range; output budget/index/byte violations throw length_error.
     * Nonfinite arithmetic throws overflow_error. No global dense field cache.
     */
    Eigen::MatrixXd pointWeights(std::size_t begin,std::size_t count) const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
}
#endif
