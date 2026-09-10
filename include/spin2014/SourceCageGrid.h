#pragma once
#include "spin2014/SourceSurfaceCage.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageGridConfig {
    int fineLevel = 4;
    std::size_t maximumCells = 1000000;
    std::size_t maximumNodes = 2000000;
    std::size_t maximumEdges = 6000000;
    std::size_t maximumInterpolationValues = 10000000;
};
/** Immutable sparse deformation/weight COMPUTATIONAL DOMAIN, not material.
 * IMPLEMENTATION_DEFAULT: uniformly refine every selected coarse cage cell;
 * the union and all coarse handles are preserved. This is not recovered author
 * discretization metadata. Require cageLevel < fineLevel <= 20 and the same root.
 * Cells/nodes are lexicographic integer (x,y,z) keys. Positions are canonical
 * FMA dyadic root coordinates; cell corners use x bit first. Edges are sorted
 * unique pairs of increasing node indices and belong to actual fine cells.
 * handleNodeIndices()[i] maps cage.handleVertexIndices()[i] without relocation.
 * No root lattice allocation, culling, coarsening, solver, geometry seal,
 * feasibility, BBW convergence, physical simulation or paper-case authority.
 * Resource failures throw length_error containing RESOURCE_LIMIT; malformed
 * inputs (including numerically unsupported finer coordinates) invalid_argument.
 * The fine step must exceed the maximum ULP at the root endpoints, ensuring
 * distinct represented planes. Equality is conservatively rejected even when
 * a particular aligned root would have distinct planes; no recenter/coarsening.
 * Source/cage hashes retain upstream meanings. configSha256 binds the domain
 * convention, upstream cage config and all fine-grid parameters/budgets;
 * outputSha256 also binds source/cage output and every derived array.
 */
class SourceCageGrid {
public:
    using Key = SourceSurfaceCage::Key;
    using CellNodes = std::array<int,8>;
    using Edge = std::array<int,2>;
    struct Embedding {
        std::size_t cellIndex = 0;
        Key cellKey{};
        CellNodes nodeIndices{};
        std::array<double,8> coefficients{};
    };
    SourceCageGrid(const SourceSurfaceCage&, const SourceCageGridConfig&);
    const SourceSurfaceCage& cage() const;
    const SourceCageGridConfig& config() const;
    const std::vector<Key>& cells() const;
    const std::vector<Key>& nodeKeys() const;
    const Eigen::MatrixXd& nodePositions() const;
    const std::vector<CellNodes>& cellNodes() const;
    const std::vector<Edge>& edges() const;
    const std::vector<int>& handleNodeIndices() const;
    const std::string& sourceSha256() const;
    const std::string& cageConfigSha256() const;
    const std::string& cageOutputSha256() const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
    /** CLOSED represented-cell lookup: shared planes choose the lexicographically
     * first incident active cell. Outside points/holes throw out_of_range, with
     * no tolerance, nearest fallback or outside-point clamping. At most eight
     * candidates are considered after binary searches of represented planes.
     * Coefficients use actual represented lower/upper planes, are nonnegative,
     * and sum to one to roundoff. The same API embeds newly split cell corners.
     */
    Embedding embed(const Eigen::Vector3d& point) const;
    /** Transfer a finite Nnodes x K field to finite Nqueries x 3 positions.
     * K must be positive. Empty queries return 0 x K. The checked output value
     * count must not exceed maximumInterpolationValues, even for generic fields.
     * No dense interpolation matrix, user-field clamping or renormalization.
     * This does NOT certify input weights' bounds, partition of unity or KKT.
     * A nonfinite arithmetic result throws overflow_error.
     */
    Eigen::MatrixXd interpolate(const Eigen::MatrixXd& nodalField, const Eigen::MatrixXd& queries) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
}
#endif
