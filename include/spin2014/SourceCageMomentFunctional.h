#pragma once
#include "spin2014/SourceConformingBoundary.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageMomentFunctionalConfig {
    std::size_t maximumWeightValues=1000000;
    std::size_t maximumJacobianValues=10000000;
};
struct SourceCageMomentEvaluation {
    physics::RawMassMoments10 finalMoments;
    Eigen::MatrixXd betaJacobian;   // 10 x B, or 0 x 0 when absent.
    Eigen::MatrixXd handleJacobian; // 10 x 4H, or 0 x 0 when absent.
    Eigen::MatrixXd deformedPoints; // Same original+material order as boundary.
    bool hasJacobians=false;
    std::string outputSha256;
};
/** Immutable actual source material moment functional, fixed-rest full BBW.
 * Handle layout is [scale,tx,ty,tz] in original cage handle order. All fixed
 * rotations are explicitly identity. beta is hollowed fraction per adaptive
 * leaf. Scales are strictly positive. This layer does not establish inversion,
 * intersection, overlap, thickness, exact final partition or solver legality.
 * Config caps are checked, never silently reduced. Weight scratch uses one
 * full handle row, with cached eight-node rest embeddings, never P x H or
 * P x 4H storage. Value-only mode allocates no dynamic derivative matrices.
 */
class SourceCageMomentFunctional {
public:
    SourceCageMomentFunctional(const SourceConformingBoundary&,double density,
        const SourceCageMomentFunctionalConfig&);
    const SourceConformingBoundary& boundary() const;
    const SourceCageMomentFunctionalConfig& config() const;
    double density() const;
    std::size_t betaCount() const;
    std::size_t handleCount() const;
    const std::string& configSha256() const;
    const std::string& upstreamSha256() const;
    Eigen::VectorXd initialHandleParameters() const;
    SourceCageMomentEvaluation evaluate(const Eigen::VectorXd& beta,
        const Eigen::VectorXd& handleParameters) const;
    SourceCageMomentEvaluation evaluateValues(const Eigen::VectorXd& beta,
        const Eigen::VectorXd& handleParameters) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    SourceCageMomentEvaluation run(const Eigen::VectorXd&,const Eigen::VectorXd&,bool) const;
};
}
#endif
