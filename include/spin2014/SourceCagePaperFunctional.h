#pragma once
#include "spin2014/SourceCageMomentFunctional.h"
#include "spin2014/Spin2014GivensProblem.h"
#include <limits>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCagePaperFrame {
    Eigen::Vector3d origin=Eigen::Vector3d::Zero();
    Eigen::Matrix3d rotation=Eigen::Matrix3d::Identity(); // q=rotation*(p-origin)
};
struct SourceCagePaperFunctionalConfig {
    Spin2014Objective objective=Spin2014Objective::StaticLowering;
    double inertiaWeight=1.,centerWeight=1.,fillWeight=0.;
    std::size_t maximumJointDerivativeValues=10000000;
};
struct SourceCagePaperEvaluation {
    SourceCageMomentEvaluation world;
    physics::RawMassMoments10 spinMoments;
    double paperObjective=0.,fillPenalty=0.,totalObjective=0.;
    Eigen::VectorXd constraints;
    Eigen::VectorXd jointGradient;
    Eigen::MatrixXd jointConstraintJacobian; // K x N, absent in value-only mode.
    Eigen::MatrixXd spinMomentJacobian;      // 10 x N, exactly zero phi column.
    bool hasJacobians=false;
    std::string outputSha256;
};
/** Immutable actual-source paper functional. Joint order [beta, optional phi,
 * all original handles (scale,tx,ty,tz)]. Static has no angle coordinate.
 * The fixed frame is explicit and never follows the moving center of mass.
 * Zero origin/identity rotation and unit-adjacency fill are IMPLEMENTATION_DEFAULT.
 * Reconstruct after partition/BBW changes. No shape energy, full joint Hessian,
 * optimizer, geometry legality or original-paper reproduction is asserted.
 * Derivative cap covers (10+K+1)*N+B doubles, separately from upstream budgets.
 * Transform scratch is two fixed ten-vectors; no size-dependent GEMM packing.
 * This constant scratch and the eleven-coordinate Givens workspace do not grow
 * with N and are separate from the joint array cap.
 */
class SourceCagePaperFunctional {
public:
    using MomentMap=Eigen::Matrix<double,10,10>;
    SourceCagePaperFunctional(const SourceCageMomentFunctional&,
        const SourceCagePaperFunctionalConfig&,const SourceCagePaperFrame& = {});
    const SourceCageMomentFunctional& source() const;
    const SourceCagePaperFunctionalConfig& config() const;
    const SourceCagePaperFrame& frame() const;
    const MomentMap& worldToSpinMomentMap() const;
    std::size_t betaCount() const;
    std::size_t handleCount() const;
    bool hasAngle() const;
    // Static returns std::numeric_limits<std::size_t>::max().
    std::size_t angleOffset() const;
    std::size_t handleOffset() const;
    std::size_t stateCount() const;
    std::size_t constraintCount() const;
    const std::string& configSha256() const;
    const std::string& upstreamSha256() const;
    Eigen::VectorXd initialState(const Eigen::VectorXd& beta,double phi) const;
    SourceCagePaperEvaluation evaluate(const Eigen::VectorXd& state) const;
    SourceCagePaperEvaluation evaluateValues(const Eigen::VectorXd& state) const;
    SourceCagePaperEvaluation evaluate(const Eigen::VectorXd& beta,double phi,
        const Eigen::VectorXd& handleParameters) const;
    SourceCagePaperEvaluation evaluateValues(const Eigen::VectorXd& beta,double phi,
        const Eigen::VectorXd& handleParameters) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    SourceCagePaperEvaluation run(const Eigen::VectorXd&,double,const Eigen::VectorXd&,bool) const;
    SourceCagePaperEvaluation runState(const Eigen::VectorXd&,bool) const;
};
namespace detail {
// Internal operation on preallocated ten-row input/output matrices. Copy each
// input into fixed storage so Eigen cannot select dynamic GEMM packing. The
// two ten-vectors are constant stack workspace, independent of the column count.
inline void transformSourceCageMomentColumns(const SourceCagePaperFunctional::MomentMap& map,
    const Eigen::MatrixXd& input,Eigen::Ref<Eigen::MatrixXd> output) {
    using Column=Eigen::Matrix<double,10,1>;
    for(Eigen::Index column=0;column<input.cols();++column) {
        const Column world=input.col(column);
        const Column spin=map*world;
        output.col(column)=spin;
    }
}
}
}
#endif
