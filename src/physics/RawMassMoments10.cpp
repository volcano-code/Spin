#include "physics/RawMassMoments10.h"

#include <Eigen/LU>

#include <cmath>
#include <limits>
#include <stdexcept>


namespace physics
{


RawMassMoments10::RawMassMoments10() = default;


RawMassMoments10::RawMassMoments10(const Vector& values)
    :
    values_(values)
{
}


RawMassMoments10 RawMassMoments10::fromMassMoments(
    const MassMoments& moments
)
{
    Vector values;
    values <<
        moments.mass,
        moments.firstMoment(0),
        moments.firstMoment(1),
        moments.firstMoment(2),
        moments.secondMoment(0, 0),
        moments.secondMoment(1, 1),
        moments.secondMoment(2, 2),
        moments.secondMoment(1, 2),
        moments.secondMoment(0, 2),
        moments.secondMoment(0, 1);
    return RawMassMoments10(values);
}


MassMoments RawMassMoments10::toMassMoments() const
{
    MassMoments moments;
    moments.mass = values_(M);
    moments.firstMoment = Eigen::Vector3d(
        values_(LX),
        values_(LY),
        values_(LZ)
    );
    moments.secondMoment <<
        values_(QXX), values_(QXY), values_(QXZ),
        values_(QXY), values_(QYY), values_(QYZ),
        values_(QXZ), values_(QYZ), values_(QZZ);
    return moments;
}


const RawMassMoments10::Vector& RawMassMoments10::values() const noexcept
{
    return values_;
}


double RawMassMoments10::mass() const noexcept
{
    return values_(M);
}


bool RawMassMoments10::isFinite() const noexcept
{
    return values_.allFinite();
}


Eigen::Vector3d RawMassMoments10::centerOfMass(
    double minimumPositiveMass
) const
{
    if (
        !isFinite() ||
        !std::isfinite(minimumPositiveMass) ||
        minimumPositiveMass < 0.0 ||
        !(mass() > minimumPositiveMass)
    )
    {
        throw std::domain_error(
            "center of mass requires finite positive raw moments"
        );
    }
    return toMassMoments().firstMoment / mass();
}


Eigen::Matrix3d RawMassMoments10::inertiaAboutOrigin() const
{
    if (!isFinite())
    {
        throw std::domain_error("origin inertia requires finite raw moments");
    }
    const Eigen::Matrix3d secondMoment = toMassMoments().secondMoment;
    return
        secondMoment.trace() * Eigen::Matrix3d::Identity() -
        secondMoment;
}


Eigen::Matrix3d RawMassMoments10::inertiaAboutCenter(
    double minimumPositiveMass
) const
{
    return toMassMoments().inertiaAboutCenter(minimumPositiveMass);
}


RawMassMoments10 RawMassMoments10::scaledDensity(double factor) const
{
    if (!std::isfinite(factor) || factor < 0.0)
    {
        throw std::invalid_argument(
            "density scale must be finite and non-negative"
        );
    }
    return RawMassMoments10(factor * values_);
}


RawMassMoments10 RawMassMoments10::transformed(
    const Eigen::Matrix3d& linear,
    const Eigen::Vector3d& translation
) const
{
    if (!isFinite() || !linear.allFinite() || !translation.allFinite())
    {
        throw std::invalid_argument("raw-moment transform must be finite");
    }
    const double determinant = std::abs(linear.determinant());
    if (
        !std::isfinite(determinant) ||
        determinant <= std::numeric_limits<double>::epsilon()
    )
    {
        throw std::invalid_argument(
            "raw-moment linear transform must be non-singular"
        );
    }

    const MassMoments source = toMassMoments();
    const Eigen::Vector3d linearFirst = linear * source.firstMoment;
    MassMoments result;
    result.mass = determinant * source.mass;
    result.firstMoment = determinant * (
        linearFirst + source.mass * translation
    );
    result.secondMoment = determinant * (
        linear * source.secondMoment * linear.transpose() +
        linearFirst * translation.transpose() +
        translation * linearFirst.transpose() +
        source.mass * translation * translation.transpose()
    );
    return fromMassMoments(result);
}


RawMassMoments10& RawMassMoments10::operator+=(
    const RawMassMoments10& other
)
{
    values_ += other.values_;
    return *this;
}


RawMassMoments10& RawMassMoments10::operator-=(
    const RawMassMoments10& other
)
{
    values_ -= other.values_;
    return *this;
}


} // namespace physics


physics::RawMassMoments10 operator+(
    physics::RawMassMoments10 first,
    const physics::RawMassMoments10& second
)
{
    first += second;
    return first;
}


physics::RawMassMoments10 operator-(
    physics::RawMassMoments10 first,
    const physics::RawMassMoments10& second
)
{
    first -= second;
    return first;
}
