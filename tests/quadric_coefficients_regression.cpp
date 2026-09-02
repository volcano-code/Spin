#include "spin/QuadricCoefficients.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-coefficients " << message << '\n';
    return EXIT_FAILURE;
}


bool throwsInvalid(const std::function<void()>& operation)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    return false;
}


Quadric make(const Vector10& values)
{
    return Quadric::fromUnnormalized(values);
}


} // namespace


int main()
{
    Vector10 coefficients = Vector10::Zero();
    coefficients(0) = -0.5;
    coefficients(1) = 1.0;
    const Quadric plane = make(coefficients);
    if (
        std::abs(plane.values().norm() - 1.0) > 1e-15 ||
        plane.evaluate(Eigen::Vector3d(0.75, 2.0, -4.0)) <= 0.0 ||
        !plane.gradient(Eigen::Vector3d::Zero()).isApprox(
            Eigen::Vector3d(1.0, 0.0, 0.0) / std::sqrt(1.25),
            1e-14
        )
    )
    {
        return fail("plane evaluation/normalization is incorrect");
    }

    coefficients.setZero();
    coefficients(0) = 1.0;
    coefficients(4) = -1.0;
    coefficients(5) = -1.0;
    coefficients(6) = -1.0;
    const Quadric sphere = make(coefficients);
    if (
        sphere.evaluate(Eigen::Vector3d::Zero()) <= 0.0 ||
        sphere.evaluate(Eigen::Vector3d(2.0, 0.0, 0.0)) >= 0.0
    )
    {
        return fail("sphere sign convention is incorrect");
    }

    const std::vector<Vector10> canonicalQuadrics = {
        // Ellipsoid: 1-x^2-2y^2-3z^2.
        (Vector10() << 1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -3.0,
            0.0, 0.0, 0.0).finished(),
        // One-sheet hyperboloid: 1+x^2+y^2-z^2.
        (Vector10() << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0,
            0.0, 0.0, 0.0).finished(),
        // Two-sheet hyperboloid: -1-x^2-y^2+z^2.
        (Vector10() << -1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 1.0,
            0.0, 0.0, 0.0).finished(),
        // Paraboloid: z-x^2-y^2.
        (Vector10() << 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 0.0,
            0.0, 0.0, 0.0).finished()
    };
    for (const Vector10& raw : canonicalQuadrics)
    {
        const Quadric quadric = make(raw);
        const Eigen::Vector3d sample(0.2, -0.3, 0.4);
        const Eigen::Vector4d homogeneous(
            sample.x(), sample.y(), sample.z(), 1.0
        );
        if (
            std::abs(
                quadric.evaluate(sample) -
                homogeneous.dot(quadric.homogeneousMatrix() * homogeneous)
            ) > 1e-14
        )
        {
            return fail("homogeneous-matrix convention is inconsistent");
        }
    }

    Eigen::Matrix3d linear;
    linear = Eigen::AngleAxisd(
        0.31,
        Eigen::Vector3d(0.2, 0.4, 0.7).normalized()
    ).toRotationMatrix();
    linear *= 1.7;
    const Eigen::Vector3d translation(0.2, -0.5, 0.9);
    const Quadric pulledBack = sphere.pullback(linear, translation);
    const Eigen::Vector3d local(0.3, -0.1, 0.6);
    const double originalValue = sphere.evaluate(linear * local + translation);
    const double scale = pulledBack.evaluate(local) / originalValue;
    const Eigen::Vector3d secondLocal(-0.7, 0.4, 0.1);
    if (
        scale <= 0.0 ||
        std::abs(
            pulledBack.evaluate(secondLocal) -
            scale * sphere.evaluate(linear * secondLocal + translation)
        ) > 1e-13
    )
    {
        return fail("similarity-coordinate pullback is incorrect");
    }

    Vector10 direction = Vector10::LinSpaced(10, -0.4, 0.5);
    const Vector10 tangent = sphere.projectToTangent(direction);
    if (std::abs(tangent.dot(sphere.values())) > 1e-14)
    {
        return fail("tangent projection is not orthogonal");
    }
    const Quadric advanced = sphere.geodesicUpdate(tangent, 0.2);
    if (
        std::abs(advanced.values().norm() - 1.0) > 1e-14 ||
        std::abs(advanced.values().dot(sphere.values()) - std::cos(0.2)) >
            1e-14
    )
    {
        return fail("unit-sphere geodesic update is incorrect");
    }
    if (
        sphere.inverted().evaluate(Eigen::Vector3d::Zero()) >= 0.0 ||
        sphere.hash().size() != 64 || sphere.hash() != sphere.hash() ||
        sphere.hash() == plane.hash()
    )
    {
        return fail("sign inversion or deterministic hash is incorrect");
    }

    coefficients.setZero();
    if (!throwsInvalid([&] { make(coefficients); }))
    {
        return fail("identically-zero quadric was accepted");
    }
    coefficients(0) = std::nan("");
    if (!throwsInvalid([&] { make(coefficients); }))
    {
        return fail("non-finite quadric was accepted");
    }
    if (!throwsInvalid([&] {
            sphere.pullback(Eigen::Matrix3d::Zero(), Eigen::Vector3d::Zero());
        }))
    {
        return fail("singular coordinate transform was accepted");
    }
    if (!throwsInvalid([&] {
            sphere.geodesicUpdate(Vector10::Zero(), 0.1);
        }))
    {
        return fail("zero geodesic direction was accepted");
    }

    std::cout << "[PASS] scenario=quadric-coefficients hash="
              << sphere.hash() << '\n';
    return EXIT_SUCCESS;
}
