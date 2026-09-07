#pragma once

#include "spin/QuadricCoefficients.h"

#include <Eigen/Core>

#include <optional>
#include <string>


namespace spin
{


enum class QuadricInverseGaussPointStatus
{
    Valid,
    InvalidNormal,
    SingularQuadraticMatrix,
    DegenerateQuadricScale,
    IdealOrInadmissibleNormal,
    BelowNumericalThreshold,
    NonFinitePoint
};


const char* quadricInverseGaussPointStatusToken(
    QuadricInverseGaussPointStatus status
) noexcept;


struct QuadricInverseGaussPointResult
{
    QuadricInverseGaussPointStatus status =
        QuadricInverseGaussPointStatus::SingularQuadraticMatrix;
    std::optional<Eigen::Vector3d> point;
    std::string diagnostic;

    bool valid() const noexcept;
};


enum class QuadricInverseGaussDifferentialStatus
{
    Valid,
    PointUndefined,
    BorderedSystemRankDeficient,
    BorderedSystemIllConditioned,
    UncertifiedNumericalSolution
};


const char* quadricInverseGaussDifferentialStatusToken(
    QuadricInverseGaussDifferentialStatus status
) noexcept;


struct QuadricInverseGaussDifferentialResult
{
    using CoefficientJacobian = Eigen::Matrix<double, 3, 10>;

    QuadricInverseGaussPointResult point;
    QuadricInverseGaussDifferentialStatus status =
        QuadricInverseGaussDifferentialStatus::PointUndefined;
    std::optional<CoefficientJacobian> coefficientJacobian;
    std::string diagnostic;

    bool valid() const noexcept;
};


/**
 * Prepared finite inverse-Gauss map for a full-rank normalized quadric.
 *
 * Point evaluation preserves the triangulator's positive-square-root branch
 * and historical denominator/wSquared cutoffs.  The optional differential is
 * with respect to perturbing the ten normalized coefficient entries and then
 * renormalizing; it is not a derivative of ideal points, truncation caps,
 * Boolean intersections, or the complete triangulation/moment map.
 */
class QuadricInverseGaussMap
{
public:
    using CoefficientJacobian = Eigen::Matrix<double, 3, 10>;

    explicit QuadricInverseGaussMap(const QuadricCoefficients10& quadric);

    QuadricInverseGaussPointResult point(
        const Eigen::Vector3d& normal
    ) const;
    QuadricInverseGaussDifferentialResult coefficientDifferential(
        const Eigen::Vector3d& fixedNormal
    ) const;

    const Eigen::Matrix3d& inverseQuadraticMatrix() const noexcept;

private:
    QuadricCoefficients10 quadric_;
    Eigen::Matrix3d inverseA_ = Eigen::Matrix3d::Zero();
    double denominator_ = 0.0;
    QuadricInverseGaussPointStatus preparationStatus_ =
        QuadricInverseGaussPointStatus::Valid;
    std::string preparationDiagnostic_;
};


} // namespace spin
