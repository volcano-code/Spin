#pragma once
#include <Eigen/Core>

namespace spin2014 { namespace detail {
using TriangleMomentValues = Eigen::Matrix<double,10,1>;
using TriangleMomentJacobian = Eigen::Matrix<double,10,9>;
enum class NormalOverflowCategory { Arithmetic, LegacyInvalidArgument };
// Supplement Algorithms 1/2, unit density and native RawMassMoments10 order.
// This local polynomial checks geometry/arithmetic; it makes no closure claim.
TriangleMomentValues surfaceTriangleMoments(const Eigen::Vector3d& a,
    const Eigen::Vector3d& b, const Eigen::Vector3d& c,
    TriangleMomentJacobian* jacobian = nullptr,
    NormalOverflowCategory normalOverflow = NormalOverflowCategory::Arithmetic);
}}
