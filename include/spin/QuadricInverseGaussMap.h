#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/QuadricCoefficientAuthority.h"

#include <Eigen/Core>

#include <limits>
#include <array>
#include <optional>
#include <string>
#include <vector>


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

/** Explicit full-rank source-coordinate representation. Each decimal-free
 * rational coordinate is accompanied by an exact algebraic enclosure of
 * absolute width <= 2^-384 in target coordinates. The point is a rational
 * approximation, NOT an exact continuous-quadric point. No binary64 roundtrip
 * is permitted before constructing this representation. Undefined points
 * retain the existing normal/wSquared gates; rank two is not substituted. */
struct QuadricRationalPointResult
{
    bool defined = false;
    std::array<std::string,3> coordinates;
    std::array<std::string,3> lower;
    std::array<std::string,3> upper;
    std::string diagnostic;
};

/** Authority-driven local real-arithmetic finite-point differential. Each
 * Jacobian entry encloses its exact algebraic value (width <= 2^-384), not
 * a binary64 bordered solve. This is NOT a derivative of rational rounding,
 * truncation roots, source caps or Boolean material. */
struct QuadricRationalFiniteDifferentialResult
{
    using Matrix = std::array<std::array<std::string,10>,3>;
    bool defined = false;
    QuadricRationalPointResult point;
    Matrix coefficientJacobian;
    Matrix lower;
    Matrix upper;
    std::string diagnostic;
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
    double rawReciprocalCondition = std::numeric_limits<double>::quiet_NaN();
    double equilibratedReciprocalCondition =
        std::numeric_limits<double>::quiet_NaN();
    std::string diagnostic;

    bool valid() const noexcept;
};

/** An ideal point (direction,0) on a fixed oriented normal edge.
 * parameter is t in (1-t)*first+t*second; roots are ordered by t.
 * Derivatives include coefficient renormalization, like the finite map. */
struct QuadricIdealEdgePoint
{
    double parameter = 0.0;
    Eigen::Vector3d direction = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,3,10> coefficientJacobian =
        Eigen::Matrix<double,3,10>::Zero();
};

struct QuadricIdealEdgeResult
{
    bool defined = false;
    std::vector<QuadricIdealEdgePoint> points;
    std::string diagnostic;
};

/** Truncate the projective edge (origin,1)--(direction,0), not the analytic
 * quadric arc, against the origin-centred sphere. Requires origin strictly
 * inside the radius. Jacobian columns: origin xyz, direction xyz, radius.
 * This is a local construction, NOT a cap/domain or face-topology audit. */
struct QuadricProjectiveRayResult
{
    bool defined = false;
    Eigen::Vector3d point = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,3,7> jacobian = Eigen::Matrix<double,3,7>::Zero();
    std::string diagnostic;
};
QuadricProjectiveRayResult truncateQuadricProjectiveRay(
    const Eigen::Vector3d& origin, const Eigen::Vector3d& direction,
    double radius);


/**
 * Prepared finite inverse-Gauss map for a full-rank normalized quadric.
 *
 * Point evaluation preserves the triangulator's positive-square-root branch,
 * rejects an exactly zero centered level, and retains the historical
 * wSquared cutoff. The centered-level denominator is evaluated from the exact
 * dyadic values of the input binary64 coefficients, avoiding a scale-dependent
 * absolute cutoff. The optional differential is
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
    /** Explicit regular rank-two finite-point path. Full 3x10 coefficient
     * derivative, including rank-changing directions, NOT a complete surface,
     * cap or Boolean derivative. Legacy full-rank methods are unchanged. */
    QuadricInverseGaussDifferentialResult paraboloidCoefficientDifferential(
        const Eigen::Vector3d& fixedNormal
    ) const;
    /** Opt-in algebraically evaluated positive-gradient finite point and
     * full coefficient differential. Exact rank of stored coefficients is
     * used: full rank or regular rank two. Retains the normal and wSquared
     * cutoffs where applicable, and all bordered derivative gates. This does
     * not alter legacy point/triangulation dispatch or certify a whole mesh. */
    QuadricInverseGaussDifferentialResult algebraicCoefficientDifferential(
        const Eigen::Vector3d& fixedNormal
    ) const;
    QuadricInverseGaussPointResult algebraicPoint(
        const Eigen::Vector3d& normal
    ) const;
    QuadricRationalPointResult rationalPoint(
        const Eigen::Vector3d& normal, bool compactRationalCoordinates = false
    ) const;
    /** Opt-in authority-driven full-rank finite source point. This does not
     * certify a whole triangulation, cap, Boolean or Jacobian. */
    static QuadricRationalPointResult rationalPointFromAuthority(
        const QuadricCoefficientAuthority10& authority,
        const Eigen::Vector3d& normal, bool compactRationalCoordinates = false);
    static QuadricRationalFiniteDifferentialResult rationalDifferentialFromAuthority(
        const QuadricCoefficientAuthority10& authority,
        const Eigen::Vector3d& normal, bool compactRationalCoordinates = false);
    QuadricIdealEdgeResult idealEdgeDifferential(
        const Eigen::Vector3d& firstNormal,
        const Eigen::Vector3d& secondNormal
    ) const;

    const Eigen::Matrix3d& inverseQuadraticMatrix() const noexcept;

private:
    QuadricRationalPointResult rationalPointImpl(
        const Eigen::Vector3d& normal, bool compactRationalCoordinates,
        const QuadricCoefficientAuthority10* authority) const;
    QuadricInverseGaussDifferentialResult coefficientDifferentialAtPoint(
        const Eigen::Vector3d& fixedNormal,
        const QuadricInverseGaussPointResult& point
    ) const;
    QuadricCoefficients10 quadric_;
    Eigen::Matrix3d inverseA_ = Eigen::Matrix3d::Zero();
    double denominator_ = 0.0;
    QuadricInverseGaussPointStatus preparationStatus_ =
        QuadricInverseGaussPointStatus::Valid;
    std::string preparationDiagnostic_;
};


} // namespace spin
