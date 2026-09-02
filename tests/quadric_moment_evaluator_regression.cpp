#include "manufacturing/ManufacturingGeometry.h"
#include "spin/QuadricMomentEvaluator.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


spin::TriangleMesh cube(double minimum, double maximum)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices <<
        minimum, minimum, minimum,
        maximum, minimum, minimum,
        maximum, maximum, minimum,
        minimum, maximum, minimum,
        minimum, minimum, maximum,
        maximum, minimum, maximum,
        maximum, maximum, maximum,
        minimum, maximum, maximum;
    mesh.faces.resize(12, 3);
    mesh.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
    return mesh;
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-moment-evaluator "
              << message << '\n';
    return EXIT_FAILURE;
}


spin::SpinDesignRequest request()
{
    spin::SpinDesignRequest result;
    result.outerSurface = cube(-2.0, 2.0);
    result.designDomain = cube(-1.0, 1.0);
    result.materialDensity = 1.0;
    result.minimumWallThickness = 0.5;
    return result;
}


spin::QuadricCoefficients10 sphere()
{
    spin::QuadricCoefficients10::Vector p =
        spin::QuadricCoefficients10::Vector::Zero();
    p(0) = 0.25;
    p(4) = -1.0;
    p(5) = -1.0;
    p(6) = -1.0;
    return spin::QuadricCoefficients10::fromUnnormalized(p);
}


} // namespace


int main()
{
    const spin::SpinDesignRequest input = request();
    spin::QuadricMomentEvaluatorSettings settings;
    settings.triangulation.icosphereSubdivisions = 3;
    const spin::QuadricMomentEvaluation inside =
        spin::evaluateQuadricMaterialMoments(input, sphere(), settings);
    const manufacturing::BoundaryAudit insideAudit =
        manufacturing::auditTriangleBoundary(
            inside.materialBoundary.vertices,
            inside.materialBoundary.faces,
            manufacturing::BoundaryOrientation::AnyNonzero
        );
    if (
        !inside.derivativeDefined || !insideAudit.strictBoundaryValid() ||
        inside.variableMomentsTarget.mass() <= 0.0 ||
        inside.interfaceTriangleCount <= 0 ||
        inside.momentJacobian.rows() != 10 ||
        inside.momentJacobian.cols() != 10 ||
        !inside.momentJacobian.allFinite() ||
        (inside.momentJacobian * sphere().values()).norm() > 2e-10 ||
        inside.geometryHash.size() != 64 || inside.topologyHash.size() != 64
    )
    {
        return fail("contained sphere material artifact is invalid");
    }

    const spin::QuadricMomentEvaluation outside =
        spin::evaluateQuadricMaterialMoments(
            input,
            sphere().inverted(),
            settings
        );
    const double volumeSum = inside.variableMomentsTarget.mass() +
        outside.variableMomentsTarget.mass();
    if (
        std::abs(volumeSum - 8.0) > 2e-10 ||
        outside.variableMomentsTarget.mass() <=
            inside.variableMomentsTarget.mass()
    )
    {
        return fail(
            "quadric sign inversion did not complement the domain: inside=" +
            std::to_string(inside.variableMomentsTarget.mass()) +
            " outside=" +
            std::to_string(outside.variableMomentsTarget.mass()) +
            " sum=" + std::to_string(volumeSum)
        );
    }

    spin::QuadricCoefficients10::Vector planeValues =
        spin::QuadricCoefficients10::Vector::Zero();
    planeValues(1) = 1.0;
    const spin::QuadricMomentEvaluation half =
        spin::evaluateQuadricMaterialMoments(
            input,
            spin::QuadricCoefficients10::fromUnnormalized(planeValues),
            settings
        );
    if (
        std::abs(half.variableMomentsTarget.mass() - 4.0) > 1e-12 ||
        std::abs(half.variableMomentsTarget.centerOfMass().x() - 0.5) >
            1e-12
    )
    {
        return fail("planar positive half-space moments are incorrect");
    }

    constexpr double step = 2e-5;
    Eigen::Matrix<double, 10, 10> finiteDifference;
    for (int coordinate = 0; coordinate < 10; ++coordinate)
    {
        spin::QuadricCoefficients10::Vector plus = sphere().values();
        spin::QuadricCoefficients10::Vector minus = sphere().values();
        plus(coordinate) += step;
        minus(coordinate) -= step;
        const Eigen::Matrix<double, 10, 1> plusMoments =
            spin::evaluateQuadricMaterialMoments(
                input,
                spin::QuadricCoefficients10::fromUnnormalized(plus),
                settings
            ).variableMomentsTarget.values();
        const Eigen::Matrix<double, 10, 1> minusMoments =
            spin::evaluateQuadricMaterialMoments(
                input,
                spin::QuadricCoefficients10::fromUnnormalized(minus),
                settings
            ).variableMomentsTarget.values();
        finiteDifference.col(coordinate) =
            (plusMoments - minusMoments) / (2.0 * step);
    }
    const double relativeDerivativeError =
        (inside.momentJacobian - finiteDifference).norm() /
        std::max(1.0, finiteDifference.norm());
    if (relativeDerivativeError > 0.08)
    {
        return fail(
            "shape Jacobian disagrees with symmetric finite differences: " +
            std::to_string(relativeDerivativeError)
        );
    }

    const std::string outerBefore =
        spin::computeSpinDesignRequestHash(input);
    const std::string outerAfter =
        spin::computeSpinDesignRequestHash(input);
    if (outerBefore != outerAfter ||
        inside.geometryHash != spin::evaluateQuadricMaterialMoments(
            input,
            sphere(),
            settings
        ).geometryHash)
    {
        return fail("evaluation mutated input or is non-deterministic");
    }

    std::cout << "[PASS] scenario=quadric-moment-evaluator derivative_error="
              << relativeDerivativeError << " geometry="
              << inside.geometryHash << '\n';
    return EXIT_SUCCESS;
}
