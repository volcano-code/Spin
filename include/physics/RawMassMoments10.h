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

private:
    explicit RawMassMoments10(const Vector& values);

    Vector values_ = Vector::Zero();
};


} // namespace physics
