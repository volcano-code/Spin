#pragma once
#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct LinearizedRotationLaplacianConfig {
    std::size_t maximumVertices=1000000;
    std::size_t maximumOneRingVertices=4096;
    std::size_t maximumSparseValues=30000000;
    std::size_t maximumConstructionWorkspaceDoubles=10000000;
    std::size_t maximumEvaluationWorkspaceDoubles=10000000;
    double relativeRankTolerance=1e-12; // IMPLEMENTATION_DEFAULT, strict full rank 7.
};
struct LinearizedRotationLaplacianEvaluation {
    double shapeEnergy=0.,displacementEnergy=0.;
    Eigen::MatrixXd residual;       // V x 3, D A^+ vPrime - L vPrime.
    Eigen::MatrixXd vertexGradient; // V x 3, absent for evaluateValues.
    bool hasGradient=false;
    std::string outputSha256;
};
/** Immutable mathematical rest-surface operator, not a source-validity seal.
 * Sorkine 2004 Eq.2,5,6,8--12, uniform one-ring reconstruction; fixed small
 * rotations. Each local fit is centered at its rest vertex and isotropically
 * normalized, using thin SVD with seven nontruncated singular directions.
 * Residual coordinates are xyz interleaved. Block CSR retains all nine values
 * per 3x3 block, including exact zeros, in row-major block order. Local columns
 * are vertex i followed by its unique sorted one-ring. No Gram is formed.
 * Construction workspace includes retained arrays and simultaneous scratch in
 * ceil(bytes/sizeof(double)) units. Evaluation arrays require 3V/6V doubles;
 * only fixed 3-vector/block scratch is additional. See the component document.
 * Topology validation excludes manifoldness and self-intersection checks.
 */
class LinearizedRotationLaplacian {
public:
    LinearizedRotationLaplacian(const Eigen::MatrixXd&,const Eigen::MatrixXi&,
        const LinearizedRotationLaplacianConfig& = {});
    static void validateConfig(const LinearizedRotationLaplacianConfig&);
    const Eigen::MatrixXd& restVertices() const;
    const Eigen::MatrixXi& faces() const;
    const LinearizedRotationLaplacianConfig& config() const;
    std::size_t vertexCount() const;
    std::size_t coefficientCount() const;
    std::size_t sparseValueCount() const;
    std::size_t maximumOneRingSize() const;
    std::size_t requiredConstructionWorkspaceDoubles() const;
    std::size_t requiredEvaluationWorkspaceDoubles(bool derivatives) const;
    const std::vector<std::size_t>& rowOffsets() const;
    const std::vector<std::size_t>& columnVertices() const;
    const std::vector<double>& coefficients() const;
    const std::string& configSha256() const;
    const std::string& operatorSha256() const;
    LinearizedRotationLaplacianEvaluation evaluate(Eigen::Ref<const Eigen::MatrixXd>,double shapeWeight,double displacementWeight) const;
    LinearizedRotationLaplacianEvaluation evaluateValues(Eigen::Ref<const Eigen::MatrixXd>,double shapeWeight,double displacementWeight) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    LinearizedRotationLaplacianEvaluation run(Eigen::Ref<const Eigen::MatrixXd>,double,double,bool) const;
};
}
#endif
