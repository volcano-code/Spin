#pragma once

#include "physics/RawMassMoments10.h"

#include <Eigen/Core>


namespace spin2014
{

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


} // namespace spin2014
