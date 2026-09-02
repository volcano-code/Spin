#pragma once

#include "physics/MassProperties.h"

#include <Eigen/Core>


namespace physics
{


/** Density-weighted raw volume moments through polynomial degree two. */
class RawMassMoments10
{
public:
    enum Index
    {
        M = 0,
        LX,
        LY,
        LZ,
        QXX,
        QYY,
        QZZ,
        QYZ,
        QXZ,
        QXY,
        Count
    };

    using Vector = Eigen::Matrix<double, Count, 1>;

    RawMassMoments10();

    static RawMassMoments10 fromMassMoments(const MassMoments& moments);

    MassMoments toMassMoments() const;
    const Vector& values() const noexcept;
    double mass() const noexcept;
    bool isFinite() const noexcept;

    Eigen::Vector3d centerOfMass(
        double minimumPositiveMass = 0.0
    ) const;
    Eigen::Matrix3d inertiaAboutOrigin() const;
    Eigen::Matrix3d inertiaAboutCenter(
        double minimumPositiveMass = 0.0
    ) const;

    RawMassMoments10 scaledDensity(double factor) const;
    RawMassMoments10 transformed(
        const Eigen::Matrix3d& linear,
        const Eigen::Vector3d& translation
    ) const;

    RawMassMoments10& operator+=(const RawMassMoments10& other);
    RawMassMoments10& operator-=(const RawMassMoments10& other);

private:
    explicit RawMassMoments10(const Vector& values);

    Vector values_ = Vector::Zero();
};


} // namespace physics


physics::RawMassMoments10 operator+(
    physics::RawMassMoments10 first,
    const physics::RawMassMoments10& second
);


physics::RawMassMoments10 operator-(
    physics::RawMassMoments10 first,
    const physics::RawMassMoments10& second
);
