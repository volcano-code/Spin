#pragma once

#include "physics/RawMassMoments10.h"

#include <Eigen/Core>


namespace spin2014
{


struct SurfaceRawMomentDifferential
{
    physics::RawMassMoments10 moments;
    /** 10 by 3N, with interleaved vertex coordinates. */
    Eigen::MatrixXd jacobian;
};


/** Analytic differential of all ten degree-two volume moments. */
SurfaceRawMomentDifferential differentiateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density = 1.0
);


} // namespace spin2014
