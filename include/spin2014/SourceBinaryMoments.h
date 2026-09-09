#pragma once
#include "spin2014/SourceAdaptiveDomain.h"
#include "spin2014/SpinIt2014MomentProblem.h"
#include <optional>
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceRemovedCuboid {
    SourcePartitionKey cell;
    Eigen::Vector3d lower,upper;
};
struct SourceBinaryCuboids {
    std::string sourceHash,partitionHash,outputHash;
    std::size_t maximumCuboids=0;
    std::vector<SourceRemovedCuboid> cuboids;
};
struct SourceBinaryMomentCheck {
    bool passed=false;
    std::string diagnostic;
    std::optional<physics::RawMassMoments10> independentMoments;
    Eigen::VectorXd constraints;
    double normalizedMomentError=std::numeric_limits<double>::infinity();
    double rawConstraintResidual=std::numeric_limits<double>::infinity();
    double normalizedConstraintResidual=std::numeric_limits<double>::infinity();
    double minimumCentralInertiaEigenvalue=std::numeric_limits<double>::infinity();
};
SourceBinaryCuboids sourceBinaryCuboids(const SourceAdaptiveDomain&,const Eigen::VectorXd& binary,std::size_t maximumCuboids);
/** Independent signed-tetrahedron integration of original outer and streamed
 * outward twelve-triangle cuboids. No Algorithm1 or oracle moment columns.
 * This is ONLY numerical moment/constraint validation; it does not validate a
 * manifold, self intersections, wall thickness, stability or SealedSpinArtifact.
 */
SourceBinaryMomentCheck checkSourceBinaryMoments(const SourceSurfaceOctree& source,
    const SourceBinaryCuboids&,const physics::RawMassMoments10& solverBinaryMoments,
    Spin2014Objective objective,double phi,double density,double constraintTolerance);
}
#endif
