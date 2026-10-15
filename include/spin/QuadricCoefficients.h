#pragma once

#include <Eigen/Core>

#include <string>


namespace spin
{


/** Unit-norm coefficients of beta(x)=a+b.x+x^T A x (paper Eq. 5). */
class QuadricCoefficients10
{
public:
    using Vector = Eigen::Matrix<double, 10, 1>;

    static QuadricCoefficients10 fromUnnormalized(const Vector& coefficients);
    /** Validate and preserve a saved unit vector without renormalizing its bits. */
    static QuadricCoefficients10 fromNormalized(const Vector& coefficients);
    /** Hash an already-normalized coefficient vector without renormalizing it. */
    static std::string hashNormalized(const Vector& coefficients);

    const Vector& values() const noexcept;
    double evaluate(const Eigen::Vector3d& point) const;
    Eigen::Vector3d gradient(const Eigen::Vector3d& point) const;
    Eigen::Matrix3d quadraticMatrix() const noexcept;
    Eigen::Vector3d linearVector() const noexcept;
    double constant() const noexcept;
    Eigen::Matrix4d homogeneousMatrix() const noexcept;

    Vector projectToTangent(const Vector& direction) const;
    QuadricCoefficients10 geodesicUpdate(
        const Vector& tangentDirection,
        double arcLength
    ) const;

    /** Return beta(linear*y+translation), renormalized without sign flips. */
    QuadricCoefficients10 pullback(
        const Eigen::Matrix3d& linear,
        const Eigen::Vector3d& translation
    ) const;

    QuadricCoefficients10 inverted() const;
    std::string hash() const;

private:
    explicit QuadricCoefficients10(const Vector& normalizedCoefficients);

    Vector coefficients_ = Vector::Zero();
};


} // namespace spin
