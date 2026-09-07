#include "spin/QuadricInverseGaussMap.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>


namespace
{


constexpr double minimumNormalNorm = 1e-15;
constexpr double minimumDenominatorMagnitude = 1e-15;
constexpr double minimumWSquared = 1e-14;
constexpr double minimumBorderedReciprocalCondition = 1e-12;


} // namespace


namespace spin
{


const char* quadricInverseGaussPointStatusToken(
    QuadricInverseGaussPointStatus status
) noexcept
{
    switch (status)
    {
    case QuadricInverseGaussPointStatus::Valid:
        return "VALID";
    case QuadricInverseGaussPointStatus::InvalidNormal:
        return "INVALID_NORMAL";
    case QuadricInverseGaussPointStatus::SingularQuadraticMatrix:
        return "SINGULAR_QUADRATIC_MATRIX";
    case QuadricInverseGaussPointStatus::DegenerateQuadricScale:
        return "DEGENERATE_QUADRIC_SCALE";
    case QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal:
        return "IDEAL_OR_INADMISSIBLE_NORMAL";
    case QuadricInverseGaussPointStatus::BelowNumericalThreshold:
        return "BELOW_NUMERICAL_THRESHOLD";
    case QuadricInverseGaussPointStatus::NonFinitePoint:
        return "NONFINITE_POINT";
    }
    return "UNKNOWN";
}


bool QuadricInverseGaussPointResult::valid() const noexcept
{
    return status == QuadricInverseGaussPointStatus::Valid &&
        point.has_value();
}


const char* quadricInverseGaussDifferentialStatusToken(
    QuadricInverseGaussDifferentialStatus status
) noexcept
{
    switch (status)
    {
    case QuadricInverseGaussDifferentialStatus::Valid:
        return "VALID";
    case QuadricInverseGaussDifferentialStatus::PointUndefined:
        return "POINT_UNDEFINED";
    case QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient:
        return "BORDERED_SYSTEM_RANK_DEFICIENT";
    case QuadricInverseGaussDifferentialStatus::BorderedSystemIllConditioned:
        return "BORDERED_SYSTEM_ILL_CONDITIONED";
    case QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution:
        return "UNCERTIFIED_NUMERICAL_SOLUTION";
    }
    return "UNKNOWN";
}


bool QuadricInverseGaussDifferentialResult::valid() const noexcept
{
    return point.valid() &&
        status == QuadricInverseGaussDifferentialStatus::Valid &&
        coefficientJacobian.has_value();
}


QuadricInverseGaussMap::QuadricInverseGaussMap(
    const QuadricCoefficients10& quadric
)
    : quadric_(quadric)
{
    const Eigen::Matrix3d a = quadric_.quadraticMatrix();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(a);
    const double rankTolerance = 1e-12 * std::max(1.0, a.norm());
    if (
        eigen.info() != Eigen::Success ||
        (eigen.eigenvalues().array().abs() <= rankTolerance).any()
    )
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::SingularQuadraticMatrix;
        preparationDiagnostic_ =
            "inverse Gauss map requires the triangulator's full-rank "
            "quadratic matrix";
        return;
    }

    inverseA_ = a.inverse();
    denominator_ =
        quadric_.linearVector().dot(
            inverseA_ * quadric_.linearVector()
        ) - 4.0 * quadric_.constant();
    if (
        !inverseA_.allFinite() || !std::isfinite(denominator_) ||
        std::abs(denominator_) <= minimumDenominatorMagnitude
    )
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::DegenerateQuadricScale;
        preparationDiagnostic_ =
            "inverse Gauss map denominator is nonfinite or within the "
            "historical 1e-15 cutoff";
    }
}


QuadricInverseGaussPointResult QuadricInverseGaussMap::point(
    const Eigen::Vector3d& normal
) const
{
    QuadricInverseGaussPointResult result;
    const double normalNorm = normal.norm();
    if (
        !normal.allFinite() || !std::isfinite(normalNorm) ||
        normalNorm <= minimumNormalNorm
    )
    {
        result.status = QuadricInverseGaussPointStatus::InvalidNormal;
        result.diagnostic =
            "inverse Gauss map requires a finite normal with norm above 1e-15";
        return result;
    }
    if (preparationStatus_ != QuadricInverseGaussPointStatus::Valid)
    {
        result.status = preparationStatus_;
        result.diagnostic = preparationDiagnostic_;
        return result;
    }

    const double wSquared =
        4.0 * normal.dot(inverseA_ * normal) / denominator_;
    if (!std::isfinite(wSquared))
    {
        result.status = QuadricInverseGaussPointStatus::NonFinitePoint;
        result.diagnostic = "inverse Gauss map produced nonfinite wSquared";
        return result;
    }
    if (wSquared <= 0.0)
    {
        result.status =
            QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal;
        result.diagnostic =
            "normal has no finite point on the positive-square-root branch";
        return result;
    }
    if (wSquared <= minimumWSquared)
    {
        result.status = QuadricInverseGaussPointStatus::BelowNumericalThreshold;
        result.diagnostic =
            "positive wSquared is within the historical 1e-14 cutoff; "
            "normal magnitude is part of this numerical contract";
        return result;
    }

    const double w = std::sqrt(wSquared);
    const Eigen::Vector3d mapped = inverseA_ *
        (normal / w - 0.5 * quadric_.linearVector());
    if (!mapped.allFinite())
    {
        result.status = QuadricInverseGaussPointStatus::NonFinitePoint;
        result.diagnostic = "inverse Gauss point is nonfinite";
        return result;
    }
    result.status = QuadricInverseGaussPointStatus::Valid;
    result.point = mapped;
    return result;
}


QuadricInverseGaussDifferentialResult
QuadricInverseGaussMap::coefficientDifferential(
    const Eigen::Vector3d& fixedNormal
) const
{
    QuadricInverseGaussDifferentialResult result;
    result.point = point(fixedNormal);
    if (!result.point.valid())
    {
        result.status = QuadricInverseGaussDifferentialStatus::PointUndefined;
        result.diagnostic = "point is undefined: " + result.point.diagnostic;
        return result;
    }

    const Eigen::Vector3d unitNormal = fixedNormal / fixedNormal.norm();
    Eigen::Matrix4d system = Eigen::Matrix4d::Zero();
    system.topLeftCorner<3, 3>() = quadric_.quadraticMatrix();
    system.topRightCorner<3, 1>() = -unitNormal;
    system.bottomLeftCorner<1, 3>() =
        quadric_.gradient(*result.point.point).transpose();

    Eigen::FullPivLU<Eigen::Matrix4d> factorization(system);
    if (factorization.rank() < 4)
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient;
        result.diagnostic = "fixed-normal bordered differential system is rank deficient";
        return result;
    }
    const double reciprocalCondition = factorization.rcond();
    if (
        !std::isfinite(reciprocalCondition) ||
        reciprocalCondition <= minimumBorderedReciprocalCondition
    )
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::BorderedSystemIllConditioned;
        result.diagnostic =
            "fixed-normal bordered differential system has reciprocal "
            "condition at or below 1e-12";
        return result;
    }

    Eigen::Matrix<double, 4, 10> rightHandSide =
        Eigen::Matrix<double, 4, 10>::Zero();
    const Eigen::Vector3d& mapped = *result.point.point;
    rightHandSide(3, 0) = -1.0;
    for (int coordinate = 0; coordinate < 3; ++coordinate)
    {
        rightHandSide(coordinate, 1 + coordinate) = -0.5;
        rightHandSide(3, 1 + coordinate) = -mapped(coordinate);
        rightHandSide(coordinate, 4 + coordinate) = -mapped(coordinate);
        rightHandSide(3, 4 + coordinate) =
            -mapped(coordinate) * mapped(coordinate);
    }
    const std::array<std::array<int, 2>, 3> crossCoordinates = {{
        {{1, 2}}, {{0, 2}}, {{0, 1}}
    }};
    for (int cross = 0; cross < 3; ++cross)
    {
        const int first = crossCoordinates[cross][0];
        const int second = crossCoordinates[cross][1];
        rightHandSide(first, 7 + cross) = -mapped(second);
        rightHandSide(second, 7 + cross) = -mapped(first);
        rightHandSide(3, 7 + cross) =
            -2.0 * mapped(first) * mapped(second);
    }

    const Eigen::Matrix<double, 4, 10> solution =
        factorization.solve(rightHandSide);
    const double relativeResidual =
        (system * solution - rightHandSide).norm() /
        std::max(1.0, rightHandSide.norm());
    if (!solution.allFinite() || !std::isfinite(relativeResidual) ||
        relativeResidual > 1e-10)
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution;
        result.diagnostic =
            "fixed-normal bordered differential solve did not certify a "
            "finite 1e-10 relative residual";
        return result;
    }

    const QuadricCoefficients10::Vector& coefficients = quadric_.values();
    const Eigen::Matrix<double, 10, 10> normalizationTangent =
        Eigen::Matrix<double, 10, 10>::Identity() -
        coefficients * coefficients.transpose();
    const CoefficientJacobian normalizedDerivative =
        solution.topRows<3>() * normalizationTangent;
    if (!normalizedDerivative.allFinite())
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution;
        result.diagnostic = "normalized coefficient differential is nonfinite";
        return result;
    }

    result.status = QuadricInverseGaussDifferentialStatus::Valid;
    result.coefficientJacobian = normalizedDerivative;
    return result;
}


const Eigen::Matrix3d& QuadricInverseGaussMap::inverseQuadraticMatrix() const
    noexcept
{
    return inverseA_;
}


} // namespace spin
