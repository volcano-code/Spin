#pragma once

#include "physics/RawMassMoments10.h"

#include <Eigen/Core>
#include <functional>


namespace spin2014
{
namespace detail { struct TriangleMomentTrace; }
/** Synchronous borrowed observations, valid only during the callback.
 * Prefixes are not completed artifacts. Exceptions propagate; callers must
 * require a successful return before accepting a full trace. No trace storage
 * is retained by the integrator. The native scaling/accumulation is unchanged.
 */
using SurfaceMomentObserver = std::function<void(Eigen::Index,
    const detail::TriangleMomentTrace&, const physics::RawMassMoments10::Vector&,
    const Eigen::MatrixXd&)>;

/** Value-only supplement Algorithm 1, with no dense vertex Jacobian.
 * Coordinates retain caller length units; moments have units rho*length^(3+d)
 * for monomial degree d. Positive finite density is required. Signed closed
 * components (including reversed cavities) are integrated without repair;
 * closure/finite/degeneracy checks do not establish material legality.
 */
physics::RawMassMoments10 integrateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
    double density = 1.0);


struct SurfaceRawMomentDifferential
{
    physics::RawMassMoments10 moments;
    /** 10 by 3N, with interleaved vertex coordinates. */
    Eigen::MatrixXd jacobian;
};


/** Spin-It 2014 supplement surface-polynomial Algorithms 1 and 2.
 * Signed components are preserved; this is not a material-legality check.
 */
SurfaceRawMomentDifferential differentiateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density = 1.0
);

SurfaceRawMomentDifferential differentiateSurfaceRawMomentsObserved(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, double density,
    const SurfaceMomentObserver& observer);


} // namespace spin2014
