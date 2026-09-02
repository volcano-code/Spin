#include "manufacturing/ManufacturingGeometry.h"
#include "spin/QuadricTriangulator.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-triangulator " << message << '\n';
    return EXIT_FAILURE;
}


Eigen::MatrixXd domainCorners()
{
    Eigen::MatrixXd vertices(8, 3);
    vertices <<
        -0.5, -0.5, -0.5,
         0.5, -0.5, -0.5,
         0.5,  0.5, -0.5,
        -0.5,  0.5, -0.5,
        -0.5, -0.5,  0.5,
         0.5, -0.5,  0.5,
         0.5,  0.5,  0.5,
        -0.5,  0.5,  0.5;
    return vertices;
}


double signedVolume(const spin::TriangleMesh& mesh)
{
    double sixVolume = 0.0;
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        const Eigen::Vector3d a = mesh.vertices.row(mesh.faces(face, 0));
        const Eigen::Vector3d b = mesh.vertices.row(mesh.faces(face, 1));
        const Eigen::Vector3d c = mesh.vertices.row(mesh.faces(face, 2));
        sixVolume += a.dot(b.cross(c));
    }
    return sixVolume / 6.0;
}


double maximumInterfaceResidual(
    const spin::QuadricTriangulation& triangulation,
    const Quadric& quadric
)
{
    double maximum = 0.0;
    for (int vertex = 0; vertex < triangulation.interfaceVertexCount; ++vertex)
    {
        maximum = std::max(
            maximum,
            std::abs(quadric.evaluate(
                triangulation.closedBoundary.vertices.row(vertex)
            ))
        );
    }
    return maximum;
}


} // namespace


int main()
{
    const Eigen::MatrixXd domain = domainCorners();
    Vector10 coefficients = Vector10::Zero();
    coefficients(0) = 1.0;
    coefficients(4) = -1.0;
    coefficients(5) = -1.0;
    coefficients(6) = -1.0;
    const Quadric sphere = Quadric::fromUnnormalized(coefficients);

    double previousError = 1e9;
    for (int subdivision = 0; subdivision <= 2; ++subdivision)
    {
        spin::QuadricTriangulationSettings settings;
        settings.icosphereSubdivisions = subdivision;
        const spin::QuadricTriangulation result =
            spin::triangulateQuadric(sphere, domain, settings);
        const manufacturing::BoundaryAudit audit =
            manufacturing::auditTriangleBoundary(
                result.closedBoundary.vertices,
                result.closedBoundary.faces,
                manufacturing::BoundaryOrientation::OutwardPositive
            );
        const double error = std::abs(
            signedVolume(result.closedBoundary) - 4.0 * std::acos(-1.0) / 3.0
        );
        if (
            result.type != spin::QuadricSurfaceType::Ellipsoid ||
            result.capFaceCount != 0 || result.boundaryLoopCount != 0 ||
            result.interfaceFaceCount != 20 * (1 << (2 * subdivision)) ||
            maximumInterfaceResidual(result, sphere) > 2e-14 ||
            !audit.strictBoundaryValid() || error >= previousError
        )
        {
            return fail("ellipsoid inverse-Gauss refinement is incorrect");
        }
        previousError = error;
    }

    coefficients.setZero();
    coefficients(0) = -1.0;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = -1.0;
    const Quadric hyperboloid = Quadric::fromUnnormalized(coefficients);
    spin::QuadricTriangulationSettings hyperboloidSettings;
    hyperboloidSettings.icosphereSubdivisions = 3;
    const spin::QuadricTriangulation truncated =
        spin::triangulateQuadric(
            hyperboloid,
            domain,
            hyperboloidSettings
        );
    const manufacturing::BoundaryAudit hyperboloidAudit =
        manufacturing::auditTriangleBoundary(
            truncated.closedBoundary.vertices,
            truncated.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (
        truncated.type != spin::QuadricSurfaceType::Hyperboloid ||
        truncated.boundaryLoopCount != 2 || truncated.capFaceCount == 0 ||
        truncated.idealNormalEdgeIntersectionCount == 0 ||
        maximumInterfaceResidual(truncated, hyperboloid) > 1e-9 ||
        !hyperboloidAudit.strictBoundaryValid()
    )
    {
        return fail("hyperboloid truncation/capping is incorrect");
    }

    coefficients.setZero();
    coefficients(0) = -0.1;
    coefficients(1) = 1.0;
    const Quadric plane = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricTriangulation planar =
        spin::triangulateQuadric(plane, domain);
    const manufacturing::BoundaryAudit planeAudit =
        manufacturing::auditTriangleBoundary(
            planar.closedBoundary.vertices,
            planar.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (
        planar.type != spin::QuadricSurfaceType::Plane ||
        planar.interfaceFaceCount != 2 ||
        maximumInterfaceResidual(planar, plane) > 1e-14 ||
        !planeAudit.strictBoundaryValid()
    )
    {
        return fail("explicit plane triangulation is incorrect");
    }

    coefficients.setZero();
    coefficients(3) = 1.0;
    coefficients(4) = -1.0;
    coefficients(5) = -1.0;
    const Quadric paraboloid = Quadric::fromUnnormalized(coefficients);
    spin::QuadricTriangulationSettings paraboloidSettings;
    paraboloidSettings.paraboloidGridResolution = 12;
    const spin::QuadricTriangulation parabolic =
        spin::triangulateQuadric(
            paraboloid,
            domain,
            paraboloidSettings
        );
    const manufacturing::BoundaryAudit paraboloidAudit =
        manufacturing::auditTriangleBoundary(
            parabolic.closedBoundary.vertices,
            parabolic.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (
        parabolic.type != spin::QuadricSurfaceType::Paraboloid ||
        parabolic.interfaceFaceCount != 2 * 12 * 12 ||
        maximumInterfaceResidual(parabolic, paraboloid) > 1e-12 ||
        !paraboloidAudit.strictBoundaryValid()
    )
    {
        return fail("explicit paraboloid triangulation is incorrect");
    }

    if (
        truncated.hash.size() != 64 ||
        truncated.hash != spin::triangulateQuadric(
            hyperboloid,
            domain,
            hyperboloidSettings
        ).hash
    )
    {
        return fail("triangulation hash is absent or non-deterministic");
    }

    std::cout << "[PASS] scenario=quadric-triangulator hash="
              << truncated.hash << '\n';
    return EXIT_SUCCESS;
}
