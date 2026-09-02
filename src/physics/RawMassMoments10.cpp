#include "physics/RawMassMoments10.h"

#include <cmath>
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


} // namespace physics
