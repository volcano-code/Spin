#pragma once
#include <Eigen/Core>
#include <array>

namespace spin2014 { namespace detail {
using TriangleMomentValues = Eigen::Matrix<double,10,1>;
using TriangleMomentJacobian = Eigen::Matrix<double,10,9>;
enum class NormalOverflowCategory { Arithmetic, LegacyInvalidArgument };
struct TriangleMomentTrace {
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    std::array<Eigen::Array3d,8> h;
    std::array<Eigen::Matrix3d,3> normalDerivative;
    std::array<Eigen::Array3d,3> d3, d8;
    TriangleMomentValues unscaledPaperOrder = TriangleMomentValues::Zero();
    bool derivativesRecorded = false;
    bool completed = false;
};
// Diagnostic observation of the actual kernel, not a second implementation.
TriangleMomentValues traceSurfaceTriangleMoments(const Eigen::Vector3d& a,
    const Eigen::Vector3d& b, const Eigen::Vector3d& c,
    TriangleMomentJacobian* jacobian, TriangleMomentTrace& trace,
    NormalOverflowCategory normalOverflow = NormalOverflowCategory::Arithmetic);
// Supplement Algorithms 1/2, unit density and native RawMassMoments10 order.
// This local polynomial checks geometry/arithmetic; it makes no closure claim.
TriangleMomentValues surfaceTriangleMoments(const Eigen::Vector3d& a,
    const Eigen::Vector3d& b, const Eigen::Vector3d& c,
    TriangleMomentJacobian* jacobian = nullptr,
    NormalOverflowCategory normalOverflow = NormalOverflowCategory::Arithmetic);
}}
