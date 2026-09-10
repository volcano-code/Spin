#pragma once
#include "spin2014/SourceCageGrid.h"
#include <Eigen/SparseCore>
#include <limits>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceCageFemMassRule { RESTRICTED_GRID_VORONOI, BARYCENTRIC_TET };
struct SourceCageFemConfig {
    SourceCageFemMassRule massRule=SourceCageFemMassRule::RESTRICTED_GRID_VORONOI;
    std::size_t maximumTetrahedra=6000000;
    // Counts all 16 potential coefficients/tet, including structural zeros.
    // Assembly streams into a separately bounded unique-entry map.
    std::size_t maximumAssemblyEntries=96000000;
    std::size_t maximumStiffnessEntries=30000000;
    std::size_t maximumBoundaryEntries=6000000;
    // Two Nnodes x K work/output arrays; checked before either allocation.
    std::size_t maximumEvaluationValues=20000000;
};
/** Immutable computational FEM prerequisite, independently of Faster.
 * IMPLEMENTATION_DEFAULT: conforming six-tet Freudenthal split about local
 * x-bit-first corners 0--7, with positive orientation and actual represented
 * coordinates. RESTRICTED_GRID_VORONOI is nearest-node volume clipped to the
 * selected union: each actual brick contributes V/8 per corner. BARYCENTRIC_TET
 * contributes Vtet/4. Neither recovers author grid/clipping/meshing metadata;
 * these masses are not material density, moments or physical certificates.
 * K is linear tetrahedral stiffness; Q=K^T M^-1 K is applied without forming Q.
 * All supplied coarse cage triangles participate in complete NODAL boundary
 * constraints, including holes and disconnected components. Interior nodes
 * are unconstrained. Sorted sparse rows use supplied handle indices; exact
 * zeros are omitted. Duplicate rows must agree within boundaryValidationTolerance.
 * The tolerance validates arithmetic/partition/local reconstruction, never
 * snaps geometry or accepts points outside a triangle. Existing trilinear
 * transfer is not continuously equal to these coarse triangular functions.
 * RESOURCE_LIMIT length_error denotes a budget/index/container/byte violation;
 * invalid_argument denotes invalid configuration; overflow_error denotes
 * unsupported arithmetic. No repair, clamp, normalization or solver authority.
 */
class SourceCageFem {
public:
    using Tetrahedron=std::array<int,4>;
    using SparseMatrix=Eigen::SparseMatrix<double,Eigen::RowMajor,int>;
    struct BoundaryEntry { int handleIndex; double weight; };
    struct BoundaryRow { int nodeIndex; std::vector<BoundaryEntry> entries; };
    static constexpr double boundaryValidationTolerance=128*std::numeric_limits<double>::epsilon();
    SourceCageFem(const SourceCageGrid&,const SourceCageFemConfig&);
    const SourceCageGrid& grid() const;
    const SourceCageFemConfig& config() const;
    const std::vector<Tetrahedron>& tetrahedra() const;
    const std::vector<double>& tetrahedronVolumes() const;
    const SparseMatrix& stiffness() const;
    const Eigen::VectorXd& lumpedMass() const;
    const std::vector<BoundaryRow>& boundaryConstraints() const;
    const std::string& sourceSha256() const;
    const std::string& gridConfigSha256() const;
    const std::string& gridOutputSha256() const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
    /** Finite Nnodes x K, K>0. Gradient of energy(), no field mutation.
     * The immutable operator digest is independent of evaluated fields.
     */
    Eigen::MatrixXd applyQ(const Eigen::MatrixXd& nodalField) const;
    /** 0.5 sum_columns x^T Q x, evaluated as positive sum 0.5 (Kx)^2/M. */
    double energy(const Eigen::MatrixXd& nodalField) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    Eigen::MatrixXd massInverseStiffness(const Eigen::MatrixXd&) const;
};
} // namespace spin2014
#endif
