#include "spin/QuadricCoefficients.h"

#include "manufacturing/SlicerValidation.h"

#include <cmath>
#include <Eigen/LU>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>


namespace
{


constexpr double minimumNorm = 1e-15;


} // namespace


namespace spin
{


QuadricCoefficients10::QuadricCoefficients10(
    const Vector& normalizedCoefficients
)
    : coefficients_(normalizedCoefficients)
{
}


QuadricCoefficients10 QuadricCoefficients10::fromUnnormalized(
    const Vector& coefficients
)
{
    if (!coefficients.allFinite())
    {
        throw std::invalid_argument("quadric coefficients must be finite");
    }
    const double norm = coefficients.norm();
    if (!std::isfinite(norm) || norm <= minimumNorm)
    {
        throw std::invalid_argument(
            "quadric coefficients must define a nonzero polynomial"
        );
    }
    return QuadricCoefficients10(coefficients / norm);
}


const QuadricCoefficients10::Vector& QuadricCoefficients10::values() const
    noexcept
{
    return coefficients_;
}


Eigen::Matrix3d QuadricCoefficients10::quadraticMatrix() const noexcept
{
    Eigen::Matrix3d matrix;
    matrix <<
        coefficients_(4), coefficients_(9), coefficients_(8),
        coefficients_(9), coefficients_(5), coefficients_(7),
        coefficients_(8), coefficients_(7), coefficients_(6);
    return matrix;
}


Eigen::Vector3d QuadricCoefficients10::linearVector() const noexcept
{
    return coefficients_.segment<3>(1);
}


double QuadricCoefficients10::constant() const noexcept
{
    return coefficients_(0);
}


double QuadricCoefficients10::evaluate(const Eigen::Vector3d& point) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument("quadric evaluation point must be finite");
    }
    return constant() + linearVector().dot(point) +
        point.dot(quadraticMatrix() * point);
}


Eigen::Vector3d QuadricCoefficients10::gradient(
    const Eigen::Vector3d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument("quadric evaluation point must be finite");
    }
    return linearVector() + 2.0 * quadraticMatrix() * point;
}


Eigen::Matrix4d QuadricCoefficients10::homogeneousMatrix() const noexcept
{
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Zero();
    matrix.topLeftCorner<3, 3>() = quadraticMatrix();
    matrix.topRightCorner<3, 1>() = 0.5 * linearVector();
    matrix.bottomLeftCorner<1, 3>() =
        0.5 * linearVector().transpose();
    matrix(3, 3) = constant();
    return matrix;
}


QuadricCoefficients10::Vector QuadricCoefficients10::projectToTangent(
    const Vector& direction
) const
{
    if (!direction.allFinite())
    {
        throw std::invalid_argument("quadric direction must be finite");
    }
    return direction - coefficients_ * coefficients_.dot(direction);
}


QuadricCoefficients10 QuadricCoefficients10::geodesicUpdate(
    const Vector& tangentDirection,
    double arcLength
) const
{
    if (!tangentDirection.allFinite() || !std::isfinite(arcLength))
    {
        throw std::invalid_argument(
            "geodesic direction and arc length must be finite"
        );
    }
    if (arcLength < 0.0 || arcLength > 0.5)
    {
        throw std::invalid_argument(
            "geodesic arc length must be in the paper interval [0,0.5]"
        );
    }
    const Vector tangent = projectToTangent(tangentDirection);
    const double tangentNorm = tangent.norm();
    if (tangentNorm <= minimumNorm)
    {
        throw std::invalid_argument(
            "geodesic update requires a nonzero tangent direction"
        );
    }
    const Vector updated =
        std::cos(arcLength) * coefficients_ +
        std::sin(arcLength) * tangent / tangentNorm;
    return fromUnnormalized(updated);
}


QuadricCoefficients10 QuadricCoefficients10::pullback(
    const Eigen::Matrix3d& linear,
    const Eigen::Vector3d& translation
) const
{
    if (!linear.allFinite() || !translation.allFinite())
    {
        throw std::invalid_argument("quadric coordinate transform must be finite");
    }
    const double determinant = linear.determinant();
    if (!std::isfinite(determinant) || std::abs(determinant) <= minimumNorm)
    {
        throw std::invalid_argument(
            "quadric coordinate transform must be invertible"
        );
    }

    const Eigen::Matrix3d oldA = quadraticMatrix();
    const Eigen::Vector3d oldB = linearVector();
    const Eigen::Matrix3d newA = linear.transpose() * oldA * linear;
    const Eigen::Vector3d newB =
        linear.transpose() * (oldB + 2.0 * oldA * translation);
    const double newConstant =
        constant() + oldB.dot(translation) +
        translation.dot(oldA * translation);

    Vector transformed;
    transformed <<
        newConstant,
        newB.x(), newB.y(), newB.z(),
        newA(0, 0), newA(1, 1), newA(2, 2),
        newA(1, 2), newA(0, 2), newA(0, 1);
    return fromUnnormalized(transformed);
}


QuadricCoefficients10 QuadricCoefficients10::inverted() const
{
    return QuadricCoefficients10(-coefficients_);
}


std::string QuadricCoefficients10::hash() const
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricCoefficients10-v1;" << std::hexfloat;
    for (int coordinate = 0; coordinate < coefficients_.size(); ++coordinate)
    {
        stream << coefficients_(coordinate) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


} // namespace spin
