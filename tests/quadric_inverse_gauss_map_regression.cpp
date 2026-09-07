#include "spin/QuadricInverseGaussMap.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;
using Jacobian = spin::QuadricInverseGaussMap::CoefficientJacobian;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-inverse-gauss-map "
              << message << '\n';
    return EXIT_FAILURE;
}


bool close(
    const Eigen::Vector3d& actual,
    const Eigen::Vector3d& expected,
    double tolerance
)
{
    return (actual - expected).norm() <= tolerance;
}


Jacobian sphereOracle(double coefficientScale, double radius,
                      const Eigen::Vector3d& unitNormal)
{
    Jacobian expected = Jacobian::Zero();
    expected.col(0) = -unitNormal / (2.0 * coefficientScale * radius);
    for (int coordinate = 0; coordinate < 3; ++coordinate)
    {
        Eigen::Vector3d basis = Eigen::Vector3d::Zero();
        basis(coordinate) = 1.0;
        expected.col(1 + coordinate) = -basis / (2.0 * coefficientScale);

        const double pointCoordinate = radius * unitNormal(coordinate);
        expected.col(4 + coordinate) =
            unitNormal * pointCoordinate * pointCoordinate /
                (2.0 * coefficientScale * radius) -
            basis * pointCoordinate / coefficientScale;
    }
    const std::array<std::array<int, 2>, 3> crossCoordinates = {{
        {{1, 2}}, {{0, 2}}, {{0, 1}}
    }};
    for (int cross = 0; cross < 3; ++cross)
    {
        const int first = crossCoordinates[cross][0];
        const int second = crossCoordinates[cross][1];
        Eigen::Vector3d firstBasis = Eigen::Vector3d::Zero();
        Eigen::Vector3d secondBasis = Eigen::Vector3d::Zero();
        firstBasis(first) = 1.0;
        secondBasis(second) = 1.0;
        const double firstPoint = radius * unitNormal(first);
        const double secondPoint = radius * unitNormal(second);
        expected.col(7 + cross) =
            unitNormal * firstPoint * secondPoint /
                (coefficientScale * radius) -
            (firstBasis * secondPoint + secondBasis * firstPoint) /
                coefficientScale;
    }
    return expected;
}


Jacobian centralDifference(
    const Vector10& normalizedCoefficients,
    const Eigen::Vector3d& normal,
    double step
)
{
    Jacobian result;
    for (int coordinate = 0; coordinate < normalizedCoefficients.size();
         ++coordinate)
    {
        Vector10 positive = normalizedCoefficients;
        Vector10 negative = normalizedCoefficients;
        positive(coordinate) += step;
        negative(coordinate) -= step;
        const spin::QuadricInverseGaussPointResult positivePoint =
            spin::QuadricInverseGaussMap(
                Quadric::fromUnnormalized(positive)
            ).point(normal);
        const spin::QuadricInverseGaussPointResult negativePoint =
            spin::QuadricInverseGaussMap(
                Quadric::fromUnnormalized(negative)
            ).point(normal);
        if (!positivePoint.valid() || !negativePoint.valid())
        {
            throw std::runtime_error(
                "finite-difference point unexpectedly became invalid"
            );
        }
        result.col(coordinate) =
            (*positivePoint.point - *negativePoint.point) / (2.0 * step);
    }
    return result;
}


} // namespace


int main()
{
    Vector10 coefficients = Vector10::Zero();
    coefficients(0) = -4.0;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = 1.0;
    const Quadric sphere = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricInverseGaussMap sphereMap(sphere);
    const Eigen::Vector3d sphereNormal(2.0, -3.0, 6.0);
    const Eigen::Vector3d sphereUnitNormal = sphereNormal.normalized();
    const spin::QuadricInverseGaussPointResult spherePoint =
        sphereMap.point(sphereNormal);
    const spin::QuadricInverseGaussDifferentialResult sphereDifferential =
        sphereMap.coefficientDifferential(sphereNormal);
    const double sphereScale = 1.0 / std::sqrt(19.0);
    const Jacobian expectedSphereDerivative =
        sphereOracle(sphereScale, 2.0, sphereUnitNormal);
    if (
        !spherePoint.valid() ||
        !close(*spherePoint.point, 2.0 * sphereUnitNormal, 2e-14) ||
        !sphereDifferential.point.valid() ||
        !sphereDifferential.valid() ||
        (*sphereDifferential.coefficientJacobian -
            expectedSphereDerivative).norm() > 2e-12 ||
        (*sphereDifferential.coefficientJacobian * sphere.values()).norm() >
            2e-13
    )
    {
        return fail(
            "repeated-eigenvalue sphere point, literal derivatives, or "
            "normalized radial-null contract is incorrect"
        );
    }

    const Eigen::Vector3d axisNormal = Eigen::Vector3d::UnitX();
    const spin::QuadricInverseGaussDifferentialResult axisSphereDifferential =
        sphereMap.coefficientDifferential(axisNormal);
    if (!axisSphereDifferential.valid())
    {
        return fail("axis-aligned sphere differential is undefined");
    }
    const Jacobian axisSphereFiniteDifference = centralDifference(
        sphere.values(),
        axisNormal,
        1e-4
    );
    const Jacobian axisSphereOracle =
        sphereOracle(sphereScale, 2.0, axisNormal);
    int analyticZeroEntryCount = 0;
    double analyticZeroAbsoluteError = 0.0;
    for (int row = 0; row < axisSphereOracle.rows(); ++row)
    {
        for (int column = 0; column < axisSphereOracle.cols(); ++column)
        {
            if (axisSphereOracle(row, column) == 0.0)
            {
                ++analyticZeroEntryCount;
                analyticZeroAbsoluteError = std::max(
                    analyticZeroAbsoluteError,
                    std::abs(
                        (*axisSphereDifferential.coefficientJacobian)(
                            row,
                            column
                        ) - axisSphereFiniteDifference(row, column)
                    )
                );
            }
        }
    }
    if (
        analyticZeroEntryCount == 0 || analyticZeroAbsoluteError > 1e-8
    )
    {
        return fail("explicit analytic-zero sphere entries violate absolute rule");
    }

    const spin::QuadricInverseGaussPointResult scaledNormalPoint =
        sphereMap.point(3.25 * sphereNormal);
    const spin::QuadricInverseGaussDifferentialResult scaledNormalDerivative =
        sphereMap.coefficientDifferential(3.25 * sphereNormal);
    if (
        !scaledNormalPoint.valid() || !scaledNormalDerivative.valid() ||
        !close(*scaledNormalPoint.point, *spherePoint.point, 2e-14) ||
        (*scaledNormalDerivative.coefficientJacobian -
            *sphereDifferential.coefficientJacobian).norm() > 2e-12
    )
    {
        return fail("positive normal rescaling changed the finite map");
    }

    const spin::QuadricInverseGaussPointResult oppositeSpherePoint =
        sphereMap.point(-sphereNormal);
    const spin::QuadricInverseGaussPointResult complementedPoint =
        spin::QuadricInverseGaussMap(sphere.inverted()).point(sphereNormal);
    if (
        !oppositeSpherePoint.valid() || !complementedPoint.valid() ||
        !close(
            *complementedPoint.point,
            *oppositeSpherePoint.point,
            2e-14
        ) ||
        close(*complementedPoint.point, *spherePoint.point, 1e-6)
    )
    {
        return fail("normal sign or complemented-quadric semantics changed");
    }

    Eigen::Matrix3d asymmetricA;
    asymmetricA <<
        1.1, 0.17, -0.08,
        0.17, 0.8, 0.12,
        -0.08, 0.12, 1.3;
    const Eigen::Vector3d center(0.2, -0.1, 0.3);
    coefficients.setZero();
    coefficients(0) = center.dot(asymmetricA * center) - 1.0;
    coefficients.segment<3>(1) = -2.0 * asymmetricA * center;
    coefficients(4) = asymmetricA(0, 0);
    coefficients(5) = asymmetricA(1, 1);
    coefficients(6) = asymmetricA(2, 2);
    coefficients(7) = asymmetricA(1, 2);
    coefficients(8) = asymmetricA(0, 2);
    coefficients(9) = asymmetricA(0, 1);
    const Quadric asymmetric = Quadric::fromUnnormalized(coefficients);
    const Eigen::Vector3d asymmetricNormal(0.4, -0.2, 0.7);
    const spin::QuadricInverseGaussDifferentialResult asymmetricDerivative =
        spin::QuadricInverseGaussMap(asymmetric).coefficientDifferential(
            asymmetricNormal
        );
    if (!asymmetricDerivative.valid())
    {
        return fail("asymmetric rotated/translated ellipsoid is undefined");
    }

    const std::array<double, 3> steps = {{3e-4, 1e-4, 3e-5}};
    std::array<double, 3> relativeErrors{};
    for (std::size_t index = 0; index < steps.size(); ++index)
    {
        const Jacobian finiteDifference = centralDifference(
            asymmetric.values(),
            asymmetricNormal,
            steps[index]
        );
        relativeErrors[index] =
            (finiteDifference -
                *asymmetricDerivative.coefficientJacobian).norm() /
            std::max(1.0, finiteDifference.norm());
        if (relativeErrors[index] > 1e-6)
        {
            return fail("three-step normalized finite difference exceeded 1e-6");
        }
    }
    if (
        !(relativeErrors[1] < relativeErrors[0] &&
          relativeErrors[2] < relativeErrors[1])
    )
    {
        return fail("central finite differences do not show the expected trend");
    }

    coefficients.setZero();
    coefficients(0) = -1.0;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = -1.0;
    const Quadric hyperboloid = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricInverseGaussMap hyperboloidMap(hyperboloid);
    const spin::QuadricInverseGaussDifferentialResult hyperboloidFinite =
        hyperboloidMap.coefficientDifferential(Eigen::Vector3d(1.0, 0.2, 0.1));
    if (!hyperboloidFinite.valid())
    {
        return fail("admissible hyperboloid differential is undefined");
    }
    std::array<double, 3> hyperboloidRelativeErrors{};
    for (std::size_t index = 0; index < steps.size(); ++index)
    {
        const Jacobian finiteDifference = centralDifference(
            hyperboloid.values(),
            Eigen::Vector3d(1.0, 0.2, 0.1),
            steps[index]
        );
        hyperboloidRelativeErrors[index] =
            (finiteDifference -
                *hyperboloidFinite.coefficientJacobian).norm() /
            std::max(1.0, finiteDifference.norm());
        if (hyperboloidRelativeErrors[index] > 1e-6)
        {
            return fail("hyperboloid finite difference exceeded 1e-6");
        }
    }
    if (
        !(hyperboloidRelativeErrors[1] < hyperboloidRelativeErrors[0] &&
          hyperboloidRelativeErrors[2] < hyperboloidRelativeErrors[1])
    )
    {
        return fail("hyperboloid central finite differences lack convergence");
    }
    const spin::QuadricInverseGaussDifferentialResult hyperboloidIdeal =
        hyperboloidMap.coefficientDifferential(Eigen::Vector3d(1.0, 0.0, 1.0));
    const spin::QuadricInverseGaussPointResult hyperboloidInadmissible =
        hyperboloidMap.point(Eigen::Vector3d(0.0, 0.0, 1.0));
    if (
        !hyperboloidFinite.valid() ||
        hyperboloidIdeal.point.status !=
            spin::QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal ||
        hyperboloidIdeal.status !=
            spin::QuadricInverseGaussDifferentialStatus::PointUndefined ||
        hyperboloidIdeal.coefficientJacobian.has_value() ||
        hyperboloidInadmissible.status !=
            spin::QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal
    )
    {
        return fail("hyperboloid finite/ideal normal states are incorrect");
    }

    const Eigen::Vector3d nearIdealNormal(
        std::sqrt(1.0 + 1e-8),
        0.0,
        1.0
    );
    const spin::QuadricInverseGaussDifferentialResult nearIdeal =
        hyperboloidMap.coefficientDifferential(nearIdealNormal);
    if (
        !nearIdeal.point.valid() || nearIdeal.valid() ||
        nearIdeal.status !=
            spin::QuadricInverseGaussDifferentialStatus::BorderedSystemIllConditioned ||
        nearIdeal.coefficientJacobian.has_value()
    )
    {
        return fail(
            "near-ideal valid point did not preserve an undefined differential"
        );
    }

    const spin::QuadricInverseGaussPointResult zeroNormal =
        sphereMap.point(Eigen::Vector3d::Zero());
    const spin::QuadricInverseGaussPointResult nonfiniteNormal =
        sphereMap.point(Eigen::Vector3d(
            std::numeric_limits<double>::quiet_NaN(),
            0.0,
            1.0
        ));
    const spin::QuadricInverseGaussPointResult infiniteNormal =
        sphereMap.point(Eigen::Vector3d(
            std::numeric_limits<double>::infinity(),
            0.0,
            1.0
        ));
    const spin::QuadricInverseGaussPointResult belowThresholdNormal =
        sphereMap.point(1e-8 * sphereUnitNormal);
    if (
        zeroNormal.status != spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        zeroNormal.point.has_value() ||
        nonfiniteNormal.status !=
            spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        nonfiniteNormal.point.has_value() ||
        infiniteNormal.status !=
            spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        infiniteNormal.point.has_value() ||
        belowThresholdNormal.status !=
            spin::QuadricInverseGaussPointStatus::BelowNumericalThreshold ||
        belowThresholdNormal.point.has_value()
    )
    {
        return fail("invalid or magnitude-threshold normal contract is implicit");
    }

    coefficients.setZero();
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = 1.0;
    const spin::QuadricInverseGaussPointResult degenerateScale =
        spin::QuadricInverseGaussMap(
            Quadric::fromUnnormalized(coefficients)
        ).point(Eigen::Vector3d::UnitX());
    coefficients(6) = 0.0;
    const spin::QuadricInverseGaussMap singularMap(
        Quadric::fromUnnormalized(coefficients)
    );
    const spin::QuadricInverseGaussPointResult singularQuadratic =
        singularMap.point(Eigen::Vector3d::UnitX());
    const spin::QuadricInverseGaussDifferentialResult singularDifferential =
        singularMap.coefficientDifferential(Eigen::Vector3d::UnitX());
    if (
        degenerateScale.status !=
            spin::QuadricInverseGaussPointStatus::DegenerateQuadricScale ||
        degenerateScale.diagnostic.empty() ||
        singularQuadratic.status !=
            spin::QuadricInverseGaussPointStatus::SingularQuadraticMatrix ||
        singularQuadratic.diagnostic.empty() ||
        singularDifferential.point.valid() || singularDifferential.valid() ||
        singularDifferential.status !=
            spin::QuadricInverseGaussDifferentialStatus::PointUndefined ||
        singularDifferential.coefficientJacobian.has_value() ||
        singularDifferential.diagnostic.empty()
    )
    {
        return fail("degenerate finite coefficient states are not explicit");
    }

    bool rejectedNonfiniteCoefficients = false;
    bool rejectedZeroCoefficients = false;
    try
    {
        coefficients.setConstant(
            std::numeric_limits<double>::infinity()
        );
        (void)Quadric::fromUnnormalized(coefficients);
    }
    catch (const std::invalid_argument&)
    {
        rejectedNonfiniteCoefficients = true;
    }
    try
    {
        (void)Quadric::fromUnnormalized(Vector10::Zero());
    }
    catch (const std::invalid_argument&)
    {
        rejectedZeroCoefficients = true;
    }
    if (!rejectedNonfiniteCoefficients || !rejectedZeroCoefficients)
    {
        return fail("existing coefficient validation was bypassed");
    }

    std::cout << "[PASS] scenario=quadric-inverse-gauss-map"
              << " fd_h3e-4=" << relativeErrors[0]
              << " fd_h1e-4=" << relativeErrors[1]
              << " fd_h3e-5=" << relativeErrors[2]
              << " hyperboloid_fd_h3e-4=" << hyperboloidRelativeErrors[0]
              << " hyperboloid_fd_h1e-4=" << hyperboloidRelativeErrors[1]
              << " hyperboloid_fd_h3e-5=" << hyperboloidRelativeErrors[2]
              << " analytic_zero_entries=" << analyticZeroEntryCount
              << " analytic_zero_abs=" << analyticZeroAbsoluteError
              << " radial_null="
              << (*sphereDifferential.coefficientJacobian * sphere.values()).norm()
              << '\n';
    return EXIT_SUCCESS;
}
