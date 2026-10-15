#include "manufacturing/ManufacturingGeometry.h"
#include "spin/QuadricTriangulator.h"
#include "spin/QuadricInverseGaussMap.h"
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactRationalCoordinates.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <igl/piecewise_constant_winding_number.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <set>
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
    std::cerr << "[FAIL] scenario=quadric-triangulator " << message << '\n';
    return EXIT_FAILURE;
}


bool throwsWithMessage(const std::function<void()>& operation,
                       const std::string& expected)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument& error)
    {
        const std::string diagnostic = error.what();
        if (diagnostic.find(expected) == std::string::npos)
        {
            std::cerr << "[DIAGNOSTIC] expected exception fragment='" << expected
                      << "' actual='" << diagnostic << "'\n";
            return false;
        }
        return true;
    }
    return false;
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

double maximumCapRadiusResidual(
    const spin::QuadricTriangulation& triangulation,
    double radius
)
{
    double maximum = 0.0;
    for (int vertex = 0; vertex < triangulation.interfaceVertexCount; ++vertex)
    {
        const double residual = std::abs(
            triangulation.closedBoundary.vertices.row(vertex).norm() - radius
        );
        if (residual < 1e-4 * std::max(1.0, radius))
        {
            maximum = std::max(maximum, residual);
        }
    }
    return maximum;
}

using GeometricTriangleKey = std::array<double, 9>;

std::vector<GeometricTriangleKey> geometricTriangleKeys(
    const spin::QuadricTriangulation& triangulation)
{
    const auto& mesh = triangulation.closedBoundary;
    std::vector<GeometricTriangleKey> result;
    result.reserve(static_cast<std::size_t>(mesh.faces.rows()));
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        std::array<std::array<double, 3>, 3> points{};
        for (int corner = 0; corner < 3; ++corner)
        {
            const Eigen::Vector3d point = mesh.vertices.row(mesh.faces(face, corner));
            for (int axis = 0; axis < 3; ++axis)
                points[corner][axis] = std::round(point[axis] * 1e14) / 1e14;
        }
        std::sort(points.begin(), points.end());
        GeometricTriangleKey key{};
        int cursor = 0;
        for (const auto& point : points)
            for (const double coordinate : point)
                key[cursor++] = coordinate;
        result.push_back(key);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void writeObj(const std::filesystem::path& path, const spin::TriangleMesh& mesh)
{
    std::ofstream output(path);
    if (!output) throw std::runtime_error("could not write diagnostic OBJ");
    output << std::setprecision(17);
    for (Eigen::Index vertex = 0; vertex < mesh.vertices.rows(); ++vertex)
        output << "v " << mesh.vertices(vertex, 0) << ' '
               << mesh.vertices(vertex, 1) << ' '
               << mesh.vertices(vertex, 2) << '\n';
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
        output << "f " << mesh.faces(face, 0) + 1 << ' '
               << mesh.faces(face, 1) + 1 << ' '
               << mesh.faces(face, 2) + 1 << '\n';
    if (!output) throw std::runtime_error("diagnostic OBJ write failed");
}

int exportDegenerateBooleanCandidate(const std::filesystem::path& directory)
{
    Vector10 coefficients;
    coefficients <<
        5.8039689277516467e-15, 1.6851299752289185e-11,
        4.6761148234240197e-13, -1.5926745456717686e-11,
        0.40440691574567661, 0.40441404629432803,
        -0.82030745800395988, -2.1774792968580463e-12,
        7.0187201890690117e-11, -6.3251619031823912e-15;
    const Quadric quadric = Quadric::fromUnnormalized(coefficients);
    Eigen::MatrixXd design = 2.0 * domainCorners();
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions = 2;
    settings.truncationRadius = 8.0;
    const auto source = spin::triangulateQuadric(quadric, design, settings);
    std::filesystem::create_directories(directory);
    spin::TriangleMesh designMesh;
    designMesh.vertices = design;
    designMesh.faces.resize(12, 3);
    designMesh.faces << 0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
                        3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    writeObj(directory / "design.obj", designMesh);
    writeObj(directory / "source.obj", source.closedBoundary);
    std::cout << "operation_candidates=INTERSECT,MINUS\n"
              << "triangulated_type="
              << spin::quadricSurfaceTypeToken(source.type) << '\n'
              << "source_vertices=" << source.closedBoundary.vertices.rows()
              << " source_faces=" << source.closedBoundary.faces.rows() << '\n';
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(argc==2 && std::string(argv[1])=="--rational-hyperboloid") {
        for(double constant:{-1.0,1.0}) {
            Vector10 p=Vector10::Zero();p(0)=constant;p(4)=p(5)=1;p(6)=-1;
            spin::QuadricTriangulationSettings settings;
            settings.icosphereSubdivisions=1;settings.truncationRadius=8;
            settings.rationalSourceCoordinates=true;
            const auto source=spin::triangulateQuadric(Quadric::fromUnnormalized(p),domainCorners(),settings);
            if(source.capFaceCount<=0 || source.exactVertices.size()!=std::size_t(source.closedBoundary.vertices.rows()) ||
               source.coordinateLower.size()!=source.exactVertices.size() || source.coordinateUpper.size()!=source.exactVertices.size() ||
               source.sourceNormalVertexByVertex.size()!=source.exactVertices.size() || source.fixedNormals.rows()==0)
                return fail("rational hyperboloid lacks source, cap, enclosure or provenance");
            if(!verification::auditExactTriangleBoundaryRational(source.exactVertices,source.closedBoundary.faces,
                verification::ExactBoundaryOrientation::OutwardPositive).facts().boundaryValid)
                return fail("rational hyperboloid source boundary is invalid");
            CGAL::Gmpz denominator(1);mpz_mul_2exp(denominator.mpz(),denominator.mpz(),384);
            const CGAL::Gmpq tolerance(CGAL::Gmpz(1),denominator);
            for(std::size_t v=0;v<source.exactVertices.size();++v)for(int axis=0;axis<3;++axis) {
                const auto x=verification::parseExactRationalCoordinate(source.exactVertices[v][axis]);
                const auto lo=verification::parseExactRationalCoordinate(source.coordinateLower[v][axis]);
                const auto hi=verification::parseExactRationalCoordinate(source.coordinateUpper[v][axis]);
                if(lo>x || hi<x || hi-lo>tolerance)return fail("hyperboloid coordinate enclosure invalid");
            }
            for(int center=source.interfaceVertexCount;center<source.closedBoundary.vertices.rows();++center) {
                std::set<int> ring;
                for(int f=source.interfaceFaceCount;f<source.closedBoundary.faces.rows();++f)
                    for(int k=0;k<3;++k)if(source.closedBoundary.faces(f,k)==center)
                        for(int j=0;j<3;++j)if(source.closedBoundary.faces(f,j)!=center)
                            ring.insert(source.closedBoundary.faces(f,j));
                if(ring.empty())return fail("rational cap center has no ring");
                for(int axis=0;axis<3;++axis) {
                    CGAL::Gmpq sum(0);for(int vertex:ring)
                        sum+=verification::parseExactRationalCoordinate(source.exactVertices[vertex][axis]);
                    if(sum/CGAL::Gmpq(int(ring.size()))!=verification::parseExactRationalCoordinate(source.exactVertices[center][axis]))
                        return fail("cap center is not computed from the same rational ring");
                }
            }
        }
        std::cout<<"[PASS] rational one/two-sheet hyperboloid source, truncation enclosures and same-ring caps\n";
        return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string(argv[1]) == "--rational-source")
    {
        Vector10 values = Vector10::Zero();
        values(0) = -1.0;
        values(4) = values(5) = values(6) = 1.0;
        const auto quadric = Quadric::fromUnnormalized(values);
        spin::QuadricTriangulationSettings settings;
        settings.icosphereSubdivisions = 1;
        settings.retainFixedNormalProvenance = true;
        const auto legacy = spin::triangulateQuadric(quadric, domainCorners(), settings);
        settings.rationalSourceCoordinates = true;
        const auto source = spin::triangulateQuadric(quadric, domainCorners(), settings);
        const auto count = static_cast<std::size_t>(source.closedBoundary.vertices.rows());
        if (source.exactVertices.size() != count || source.coordinateLower.size() != count ||
            source.coordinateUpper.size() != count || source.fixedNormals.rows() != static_cast<int>(count))
            return fail("rational source must preserve complete coordinate enclosures and fixed normals");
        const auto facts = verification::auditExactTriangleBoundaryRational(
            source.exactVertices, source.closedBoundary.faces,
            verification::ExactBoundaryOrientation::OutwardPositive).facts();
        if (!facts.boundaryValid || source.hash == legacy.hash)
            return fail("rational source needs its own identity and an actual valid authority boundary");
        const spin::QuadricInverseGaussMap map(quadric);
        const CGAL::Gmpq maximumWidth(CGAL::Gmpz(1), CGAL::Gmpz(1) << 384);
        bool differsFromBinary64 = false;
        for (std::size_t vertex = 0; vertex < count; ++vertex)
        {
            const auto expected = map.rationalPoint(source.fixedNormals.row(vertex));
            if (!expected.defined || expected.coordinates != source.exactVertices[vertex])
                return fail("source authority was not generated by the high precision finite map");
            for (int axis = 0; axis < 3; ++axis)
            {
                const auto x = verification::parseExactRationalCoordinate(source.exactVertices[vertex][axis]);
                const auto lo = verification::parseExactRationalCoordinate(source.coordinateLower[vertex][axis]);
                const auto hi = verification::parseExactRationalCoordinate(source.coordinateUpper[vertex][axis]);
                if (x < lo || x > hi || hi - lo > maximumWidth)
                    return fail("source enclosure lost the approved target coordinate precision");
                differsFromBinary64 |= x != CGAL::Gmpq(source.closedBoundary.vertices(vertex, axis));
            }
        }
        if (!differsFromBinary64 || source.hash != spin::triangulateQuadric(quadric, domainCorners(), settings).hash)
            return fail("rational source used binary64-first coordinates or lacks deterministic identity");
        values.setZero(); values(0) = -1; values(4) = values(5) = 1;
        if (!throwsWithMessage([&] {
                spin::triangulateQuadric(Quadric::fromUnnormalized(values), domainCorners(), settings);
            }, "rational source coordinates require a supported full-rank source"))
            return fail("unsupported rational source chart silently fell back to binary64");
        std::cout << "[PASS] scenario=rational-source-ellipsoid hash=" << source.hash << '\n';
        return EXIT_SUCCESS;
    }
#endif
    if (argc == 3 && std::string(argv[1]) == "--export-degenerate-candidate")
    {
        try
        {
            return exportDegenerateBooleanCandidate(argv[2]);
        }
        catch (const std::exception& error)
        {
            std::cerr << "[FAIL] candidate export: " << error.what() << '\n';
            return EXIT_FAILURE;
        }
    }
    const Eigen::MatrixXd domain = domainCorners();
    Vector10 coefficients = Vector10::Zero();
    coefficients(0) = 1.0;
    coefficients(4) = -1.0;
    coefficients(5) = -1.0;
    coefficients(6) = -1.0;
    const Quadric sphere = Quadric::fromUnnormalized(coefficients);

    spin::QuadricTriangulationSettings rankTwoOnly;
    rankTwoOnly.rankTwoArrangementSource = true;
    if (!throwsWithMessage([&] {
            spin::triangulateQuadric(sphere, domain, rankTwoOnly);
        }, "rank-two arrangement requires PreparedRankTwoSource"))
        return fail("legacy entry silently accepted the prepared-source profile");

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
    const spin::QuadricTriangulation invertedSphere =
        spin::triangulateQuadric(sphere.inverted(), domain);
    const manufacturing::BoundaryAudit invertedSphereAudit =
        manufacturing::auditTriangleBoundary(
            invertedSphere.closedBoundary.vertices,
            invertedSphere.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (!invertedSphereAudit.strictBoundaryValid())
    {
        return fail("inverted ellipsoid is not an outward Boolean solid");
    }

    coefficients.setZero();
    coefficients(0) = -1e-14;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = 1.0;
    const Quadric smallEllipsoid = Quadric::fromUnnormalized(coefficients);
    if (spin::classifyQuadricSurfaceSubtype(smallEllipsoid) !=
            spin::QuadricSurfaceSubtype::Ellipsoid ||
        spin::classifyQuadricSurfaceSubtype(smallEllipsoid.inverted()) !=
            spin::QuadricSurfaceSubtype::Ellipsoid)
    {
        return fail("a nonempty full-rank ellipsoid near the empty-set boundary was rejected");
    }
    coefficients(0) = 1e-14;
    const Quadric nearEmptyEllipsoid = Quadric::fromUnnormalized(coefficients);
    if (!throwsWithMessage([&] {
            (void)spin::classifyQuadricSurfaceSubtype(nearEmptyEllipsoid);
        }, "ellipsoid has no real surface"))
    {
        return fail("an empty full-rank definite quadric near zero was accepted as a surface");
    }

    // A centered full-rank indefinite quadratic is a cone when its shifted
    // level is zero. It must not fall through the non-degenerate hyperboloid
    // inverse-Gauss path.
    coefficients.setZero();
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = -1.0;
    const Quadric cone = Quadric::fromUnnormalized(coefficients);
    if (spin::classifyQuadricSurfaceSubtype(cone) !=
            spin::QuadricSurfaceSubtype::Cone ||
        !throwsWithMessage([&] {
            (void)spin::triangulateQuadric(cone, domain);
        }, "degenerate cone"))
    {
        return fail("cone was not typed and rejected before regular hyperboloid triangulation");
    }
    spin::QuadricTriangulationSettings invalidAdaptiveResolution;
    invalidAdaptiveResolution.icosphereSubdivisions = 3;
    invalidAdaptiveResolution.maximumIcosphereSubdivisions = 2;
    if (!throwsWithMessage([&] {
            (void)spin::triangulateQuadric(sphere, domain,
                                           invalidAdaptiveResolution);
        }, "invalid quadric triangulation settings"))
    {
        return fail("adaptive normal-sphere bound below the requested level was accepted");
    }
    invalidAdaptiveResolution.maximumIcosphereSubdivisions = 8;
    if (!throwsWithMessage([&] {
            (void)spin::triangulateQuadric(sphere, domain,
                                           invalidAdaptiveResolution);
        }, "invalid quadric triangulation settings"))
    {
        return fail("adaptive normal-sphere bound above level 7 was accepted");
    }

    coefficients(0) = 1e-14;
    const Quadric nearCone = Quadric::fromUnnormalized(coefficients);
    if (spin::classifyQuadricSurfaceSubtype(nearCone) !=
        spin::QuadricSurfaceSubtype::TwoSheetHyperboloid)
    {
        return fail("a full-rank quadric with a nonzero centered level was classified as degenerate");
    }
    coefficients(0) = -1e-14;
    const Quadric nearConeOppositeSide =
        Quadric::fromUnnormalized(coefficients);
    if (spin::classifyQuadricSurfaceSubtype(nearConeOppositeSide) !=
        spin::QuadricSurfaceSubtype::OneSheetHyperboloid)
    {
        return fail("the opposite nonzero side of a cone was not classified as a regular hyperboloid");
    }

    coefficients.setZero();
    coefficients(4) = 1.0;
    coefficients(5) = -1.0;
    coefficients(0) = -1.0;
    const Quadric hyperbolicCylinder = Quadric::fromUnnormalized(coefficients);
    const auto hyperbolicCylinderSubtype =
        spin::classifyQuadricSurfaceSubtype(hyperbolicCylinder);
    if (hyperbolicCylinderSubtype != spin::QuadricSurfaceSubtype::HyperbolicCylinder)
    {
        return fail("rank-two indefinite quadric subtype=" +
            std::to_string(static_cast<int>(hyperbolicCylinderSubtype)) +
            " instead of hyperbolic cylinder");
    }
    const auto hyperbolicCylinderMesh =
        spin::triangulateQuadric(hyperbolicCylinder, domain);
    const auto hyperbolicCylinderAudit = manufacturing::auditTriangleBoundary(
        hyperbolicCylinderMesh.closedBoundary.vertices,
        hyperbolicCylinderMesh.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    if (hyperbolicCylinderMesh.type != spin::QuadricSurfaceType::HyperbolicCylinder ||
        hyperbolicCylinderMesh.boundaryLoopCount != 2 ||
        hyperbolicCylinderMesh.interfaceFaceCount != 4 * 24 * 24 ||
        hyperbolicCylinderMesh.capFaceCount <= 0 ||
        !hyperbolicCylinderAudit.strictBoundaryValid() ||
        maximumInterfaceResidual(hyperbolicCylinderMesh, hyperbolicCylinder) > 2e-12)
    {
        return fail("rank-two indefinite quadric did not produce its audited hyperbolic-cylinder boundary");
    }

    coefficients(0) = 0.0;
    const Quadric intersectingPlanes = Quadric::fromUnnormalized(coefficients);
    if (spin::classifyQuadricSurfaceSubtype(intersectingPlanes) !=
            spin::QuadricSurfaceSubtype::IntersectingPlanes ||
        !throwsWithMessage([&] {
            (void)spin::triangulateQuadric(intersectingPlanes, domain);
        }, "intersecting planes"))
    {
        return fail("zero-level rank-two indefinite quadric was not typed as intersecting planes");
    }
    const Eigen::Matrix3d rotatedFrame = Eigen::AngleAxisd(
        0.37, Eigen::Vector3d(0.2, 0.7, -0.3).normalized()).toRotationMatrix();
    const Eigen::Vector3d translatedFrame(0.23, -0.17, 0.31);
    if (spin::classifyQuadricSurfaceSubtype(
            hyperbolicCylinder.pullback(rotatedFrame, translatedFrame)) !=
            spin::QuadricSurfaceSubtype::HyperbolicCylinder ||
        spin::classifyQuadricSurfaceSubtype(
            intersectingPlanes.pullback(rotatedFrame, translatedFrame)) !=
            spin::QuadricSurfaceSubtype::IntersectingPlanes)
    {
        return fail("rank-two degenerate classification changed under rigid coordinate transformation");
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

    // Captured from the strict static-lowering solver near the regular
    // hyperboloid/paraboloid boundary. At radius 8 a one-loop cap is closed,
    // but its fan center lies in the design-domain box. The triangulator must
    // keep refining until every cap triangle is certified outside the domain;
    // it must not accept the first nonempty loop merely to avoid the old
    // exactly-two-loops restriction. This remains full rank and is not
    // silently reclassified as a paraboloid.
    Vector10 nearParabolicCoefficients;
    nearParabolicCoefficients <<
        -0.03229148714369761, -3.5981802826533382e-05,
        0.00013923083159621887, 0.99947721704275616,
        -0.00011842626037807066, 0.0010909647938445234,
        0.00047307981633996752, -0.00023316715451590989,
        0.0005345316953516031, -0.00087387989594205962;
    const Quadric nearParabolicHyperboloid =
        Quadric::fromUnnormalized(nearParabolicCoefficients);
    spin::QuadricTriangulationSettings nearParabolicSettings;
    nearParabolicSettings.icosphereSubdivisions = 2;
    nearParabolicSettings.truncationRadius = 8.0;
    const auto nearParabolic = spin::triangulateQuadric(
        nearParabolicHyperboloid, domain, nearParabolicSettings);
    const auto nearParabolicAudit = manufacturing::auditTriangleBoundary(
        nearParabolic.closedBoundary.vertices, nearParabolic.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    std::vector<bool> nearParabolicCapVertices(
        static_cast<std::size_t>(nearParabolic.closedBoundary.vertices.rows()), false);
    for (Eigen::Index face = 0; face < nearParabolic.faceIsCap.size(); ++face)
    {
        if (nearParabolic.faceIsCap(face) == 0) continue;
        for (int corner = 0; corner < 3; ++corner)
            nearParabolicCapVertices[static_cast<std::size_t>(
                nearParabolic.closedBoundary.faces(face, corner))] = true;
    }
    double nearParabolicCapDistance = std::numeric_limits<double>::infinity();
    for (Eigen::Index vertex = 0;
         vertex < nearParabolic.closedBoundary.vertices.rows(); ++vertex)
    {
        if (!nearParabolicCapVertices[static_cast<std::size_t>(vertex)]) continue;
        const Eigen::Vector3d point =
            nearParabolic.closedBoundary.vertices.row(vertex);
        const Eigen::Vector3d outside =
            (domain.colwise().minCoeff().transpose() - point).cwiseMax(0.0) +
            (point - domain.colwise().maxCoeff().transpose()).cwiseMax(0.0);
        nearParabolicCapDistance = std::min(nearParabolicCapDistance, outside.norm());
    }
    if (nearParabolic.type != spin::QuadricSurfaceType::Hyperboloid ||
        nearParabolic.boundaryLoopCount != 2 ||
        nearParabolic.truncationRadius <= 2048.0 ||
        nearParabolicCapDistance <= 0.0 ||
        maximumInterfaceResidual(nearParabolic, nearParabolicHyperboloid) > 1e-9 ||
        maximumCapRadiusResidual(nearParabolic, nearParabolic.truncationRadius) >
            1e-10 * std::max(1.0, nearParabolic.truncationRadius) ||
        !nearParabolicAudit.strictBoundaryValid())
    {
        return fail("near-parabolic cap refinement must continue until the two-loop source is certified outside the design domain");
    }

    // Captured from the answer-independent Cherchi->Newton integration
    // diagnostic. This regular hyperboloid hits an ambiguous sampled cap
    // vertex at the automatic truncation radius; the configured radius
    // expansion should recover a valid two-loop closed source mesh.
    Vector10 capVertexTransitionCoefficients;
    capVertexTransitionCoefficients <<
        -0.19114955600106495, 0.01880064140363415,
        -0.016623229690534078, -0.0009150027976023531,
        0.5769716578821628, 0.7933362567693745,
        -0.014036311317296064, -0.0008230899664630951,
        -0.0015325218795791066, 0.018773909768859762;
    const Quadric capVertexTransition =
        Quadric::fromUnnormalized(capVertexTransitionCoefficients);
    spin::QuadricTriangulationSettings capVertexSettings;
    capVertexSettings.icosphereSubdivisions = 3;
    capVertexSettings.maximumRadiusDoublings = 4;
    Eigen::MatrixXd capVertexDomain = 2.0 * domainCorners();
    bool lowResolutionCapAmbiguityWasBounded = false;
    try
    {
        spin::QuadricTriangulationSettings lowResolutionSettings = capVertexSettings;
        lowResolutionSettings.icosphereSubdivisions = 2;
        (void)spin::triangulateQuadric(
            capVertexTransition, capVertexDomain, lowResolutionSettings);
    }
    catch (const std::runtime_error& error)
    {
        lowResolutionCapAmbiguityWasBounded =
            std::string(error.what()).find("failed after bounded radius refinement") !=
            std::string::npos;
    }
    if (!lowResolutionCapAmbiguityWasBounded)
    {
        return fail("low-resolution cap ambiguity was not explicitly bounded/rejected");
    }
    const auto capVertexGeometry = spin::triangulateQuadric(
        capVertexTransition, capVertexDomain, capVertexSettings);
    const auto capVertexAudit = manufacturing::auditTriangleBoundary(
        capVertexGeometry.closedBoundary.vertices,
        capVertexGeometry.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    if (capVertexGeometry.type != spin::QuadricSurfaceType::Hyperboloid ||
        capVertexGeometry.boundaryLoopCount != 2 ||
        capVertexGeometry.truncationRadius <= 2.0 ||
        !capVertexAudit.strictBoundaryValid())
    {
        return fail("automatic radius expansion did not resolve a regular hyperboloid cap-vertex transition");
    }

    // Exact coefficients from answer-independent Cherchi->Newton evaluation
    // 59 (r14). The old per-face orientation path produced a branching
    // directed truncation boundary for this regular hyperboloid.
    Vector10 cherchiNewtonCapLoopCoefficients;
    cherchiNewtonCapLoopCoefficients <<
        -0.19223703824980956, 0.029893982546063312,
        -0.040178002262895494, 0.0008982683924855251,
        0.5777227428369576, 0.7914604728639738,
        -0.0011177055319905422, 0.0005885695113244015,
        -0.0018721422682393507, 0.018916762064737025;
    const Quadric cherchiNewtonCapLoopQuadric =
        Quadric::fromUnnormalized(cherchiNewtonCapLoopCoefficients);
    spin::QuadricTriangulationSettings cherchiNewtonCapLoopSettings;
    cherchiNewtonCapLoopSettings.icosphereSubdivisions = 3;
    cherchiNewtonCapLoopSettings.maximumRadiusDoublings = 16;
    auto refinedCherchiNewtonCapLoopSettings = cherchiNewtonCapLoopSettings;
    refinedCherchiNewtonCapLoopSettings.icosphereSubdivisions = 5;
    refinedCherchiNewtonCapLoopSettings.maximumRadiusDoublings = 4;
    bool refinedCherchiNewtonCapLoopPassed = false;
    std::string refinedCherchiNewtonCapLoopHash;
    try
    {
        const auto geometry = spin::triangulateQuadric(
            cherchiNewtonCapLoopQuadric, capVertexDomain,
            refinedCherchiNewtonCapLoopSettings);
        const auto audit = manufacturing::auditTriangleBoundary(
            geometry.closedBoundary.vertices, geometry.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive);
        refinedCherchiNewtonCapLoopPassed =
            geometry.type == spin::QuadricSurfaceType::Hyperboloid &&
            geometry.boundaryLoopCount == 2 && audit.strictBoundaryValid() &&
            maximumInterfaceResidual(
                geometry, cherchiNewtonCapLoopQuadric) <= 1e-9 &&
            maximumCapRadiusResidual(
                geometry, geometry.truncationRadius) <=
                1e-10 * std::max(1.0, geometry.truncationRadius);
        refinedCherchiNewtonCapLoopHash = geometry.hash;
        const auto repeatedGeometry = spin::triangulateQuadric(
            cherchiNewtonCapLoopQuadric, capVertexDomain,
            refinedCherchiNewtonCapLoopSettings);
        refinedCherchiNewtonCapLoopPassed =
            refinedCherchiNewtonCapLoopPassed &&
            repeatedGeometry.hash == refinedCherchiNewtonCapLoopHash;
    }
    catch (const std::exception& error)
    {
        std::cerr << "[DIAGNOSTIC] level-5 captured cap-loop case: "
                  << error.what() << '\n';
    }
    bool cherchiNewtonCapLoopFailureWasExplicit = false;
    bool cherchiNewtonCapLoopGeometryPassed = false;
    try
    {
        const auto geometry = spin::triangulateQuadric(
            cherchiNewtonCapLoopQuadric, capVertexDomain,
            cherchiNewtonCapLoopSettings);
        const auto audit = manufacturing::auditTriangleBoundary(
            geometry.closedBoundary.vertices, geometry.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive);
        cherchiNewtonCapLoopGeometryPassed =
            geometry.type == spin::QuadricSurfaceType::Hyperboloid &&
            geometry.boundaryLoopCount == 2 && audit.strictBoundaryValid();
    }
    catch (const std::exception& error)
    {
        std::cerr << "[DIAGNOSTIC] captured Cherchi/Newton cap-loop case: "
                  << error.what() << '\n';
        const std::string diagnostic = error.what();
        cherchiNewtonCapLoopFailureWasExplicit =
            diagnostic.find("failed after bounded radius refinement") !=
                std::string::npos &&
            diagnostic.find("outgoing edges at vertex") != std::string::npos &&
            diagnostic.find("incoming edges") != std::string::npos;
    }
    if (!refinedCherchiNewtonCapLoopPassed ||
        (!cherchiNewtonCapLoopFailureWasExplicit &&
         !cherchiNewtonCapLoopGeometryPassed))
    {
        return fail("captured Cherchi/Newton hyperboloid must either reject the frozen level-3 topology explicitly or pass the independently audited level-5 mesh");
    }
    if (cherchiNewtonCapLoopFailureWasExplicit)
    {
        auto adaptiveSettings = cherchiNewtonCapLoopSettings;
        adaptiveSettings.maximumIcosphereSubdivisions = 5;
        const auto adaptiveGeometry = spin::triangulateQuadric(
            cherchiNewtonCapLoopQuadric, capVertexDomain, adaptiveSettings);
        const auto adaptiveRepeat = spin::triangulateQuadric(
            cherchiNewtonCapLoopQuadric, capVertexDomain, adaptiveSettings);
        if (adaptiveGeometry.effectiveIcosphereSubdivisions != 5 ||
            adaptiveGeometry.boundaryLoopCount != 2 ||
            adaptiveGeometry.hash != adaptiveRepeat.hash ||
            adaptiveGeometry.hash == refinedCherchiNewtonCapLoopHash)
        {
            return fail("opt-in bounded normal-sphere refinement must select level 5 deterministically and identify its output separately");
        }
    }

    // A frozen seed-1 Faster trial is an off-center one-sheet hyperboloid
    // whose radius-27.7128 truncation has one boundary loop. That loop count
    // is valid, but its planar centroid fan is not certified outside the
    // design-domain AABB, so the source must continue refining. At a larger
    // radius the two-loop cap is safely separated. This distinguishes the
    // real closure invariant from the former exactly-two-loops assumption.
    Vector10 offsetOneSheetCoefficients;
    offsetOneSheetCoefficients <<
        0.0009097158036885176, -0.18034303110789637,
        0.2838905404225007, 0.5507642782534051,
        0.001474335857157634, 0.0048412627066542775,
        -0.005674017347961097, -0.3391663222577052,
        -0.6844318485835988, 0.0013707419667863376;
    const Quadric offsetOneSheet =
        Quadric::fromUnnormalized(offsetOneSheetCoefficients);
    Eigen::MatrixXd requestSizedDomain = 0.5 * domainCorners();
    const Eigen::Vector3d requestDomainCenter(0.03125, -0.03125, -0.0625);
    requestSizedDomain.rowwise() += requestDomainCenter.transpose();
    spin::QuadricTriangulationSettings singleLoopCapSettings;
    singleLoopCapSettings.icosphereSubdivisions = 4;
    singleLoopCapSettings.truncationRadius = 27.712812921102035;
    singleLoopCapSettings.maximumRadiusDoublings = 0;
    bool oneLoopCapDomainIntrusionWasRejected = false;
    try
    {
        (void)spin::triangulateQuadric(
            offsetOneSheet, requestSizedDomain, singleLoopCapSettings);
    }
    catch (const std::runtime_error& error)
    {
        oneLoopCapDomainIntrusionWasRejected =
            std::string(error.what()).find(
                "cap triangle could not be certified outside the design-domain AABB") !=
            std::string::npos;
    }
    if (!oneLoopCapDomainIntrusionWasRejected)
        return fail("one-loop hyperboloid cap that lacks an outside-domain certificate must be rejected");

    auto completeOffsetCapSettings = singleLoopCapSettings;
    completeOffsetCapSettings.truncationRadius = 0.0;
    completeOffsetCapSettings.maximumRadiusDoublings = 6;
    const auto completeOffsetCap = spin::triangulateQuadric(
        offsetOneSheet, requestSizedDomain, completeOffsetCapSettings);
    const auto completeOffsetCapAudit = manufacturing::auditTriangleBoundary(
        completeOffsetCap.closedBoundary.vertices,
        completeOffsetCap.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    std::vector<bool> capVertices(
        static_cast<std::size_t>(completeOffsetCap.closedBoundary.vertices.rows()),
        false);
    for (Eigen::Index face = 0; face < completeOffsetCap.faceIsCap.size(); ++face)
    {
        if (completeOffsetCap.faceIsCap(face) == 0) continue;
        for (int corner = 0; corner < 3; ++corner)
            capVertices[static_cast<std::size_t>(
                completeOffsetCap.closedBoundary.faces(face, corner))] = true;
    }
    double minimumCapDistanceToDomain = std::numeric_limits<double>::infinity();
    for (Eigen::Index vertex = 0;
         vertex < completeOffsetCap.closedBoundary.vertices.rows(); ++vertex)
    {
        if (!capVertices[static_cast<std::size_t>(vertex)]) continue;
        const Eigen::Vector3d point =
            completeOffsetCap.closedBoundary.vertices.row(vertex);
        const Eigen::Vector3d outside =
            (requestSizedDomain.colwise().minCoeff().transpose() - point).cwiseMax(0.0) +
            (point - requestSizedDomain.colwise().maxCoeff().transpose()).cwiseMax(0.0);
        minimumCapDistanceToDomain = std::min(minimumCapDistanceToDomain, outside.norm());
    }
    const auto completeOffsetCapRepeat = spin::triangulateQuadric(
        offsetOneSheet, requestSizedDomain, completeOffsetCapSettings);
    if (completeOffsetCap.type != spin::QuadricSurfaceType::Hyperboloid ||
        spin::classifyQuadricSurfaceSubtype(offsetOneSheet) !=
            spin::QuadricSurfaceSubtype::OneSheetHyperboloid ||
        completeOffsetCap.boundaryLoopCount != 2 ||
        completeOffsetCap.truncationRadius <= 100.0 ||
        minimumCapDistanceToDomain <= 0.0 ||
        completeOffsetCap.hash != completeOffsetCapRepeat.hash ||
        !completeOffsetCapAudit.strictBoundaryValid() ||
        maximumInterfaceResidual(completeOffsetCap, offsetOneSheet) > 1e-9 ||
        maximumCapRadiusResidual(
            completeOffsetCap, completeOffsetCap.truncationRadius) >
            1e-10 * std::max(1.0, completeOffsetCap.truncationRadius))
    {
        return fail("radius refinement must select a deterministic, closed cap outside the locked design domain");
    }
    if (argc == 2 && std::string(argv[1]) == "--adaptive-cap-refinement-only")
    {
        std::cout << "[PASS] bounded adaptive cap refinement; effective_level="
                  << (cherchiNewtonCapLoopFailureWasExplicit ? 5 : 3) << '\n';
        return EXIT_SUCCESS;
    }

    // Captured from a seed-only Newton attempt. The surface equation residual
    // is tiny, while the old fixed-count normal bisection did not guarantee
    // this stricter truncation-sphere residual.
    Vector10 sourceMotionTrial;
    sourceMotionTrial <<
        8.8660798987646838e-07, -0.003175536459697936,
        -0.0031795530012907304, -0.0021001482598485732,
        0.4038525237044106, 0.40469237988421414,
        -0.820395126919088, -0.0046207471997323137,
        -0.0057270451094420232, -0.0005458256856955998;
    const Quadric sourceMotionQuadric = Quadric::fromUnnormalized(sourceMotionTrial);
    spin::QuadricTriangulationSettings sourceMotionSettings;
    sourceMotionSettings.icosphereSubdivisions = 2;
    sourceMotionSettings.truncationRadius = 8.0;
    const auto sourceMotionGeometry = spin::triangulateQuadric(
        sourceMotionQuadric, domainCorners(), sourceMotionSettings
    );
    if (maximumCapRadiusResidual(
            sourceMotionGeometry, sourceMotionGeometry.truncationRadius) >
        1e-10 * std::max(1.0, sourceMotionGeometry.truncationRadius))
    {
        return fail("hyperboloid cap crossings violate the strict source-motion radius contract");
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

    coefficients.setZero();
    coefficients(0) = -0.25;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    const Quadric cylinder = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricTriangulation cylindrical =
        spin::triangulateQuadric(cylinder, domain);
    const manufacturing::BoundaryAudit cylinderAudit =
        manufacturing::auditTriangleBoundary(
            cylindrical.closedBoundary.vertices,
            cylindrical.closedBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (
        cylindrical.type != spin::QuadricSurfaceType::EllipticCylinder ||
        cylindrical.boundaryLoopCount != 2 ||
        cylindrical.interfaceFaceCount != 2 * 96 ||
        maximumInterfaceResidual(cylindrical, cylinder) > 1e-12 ||
        !cylinderAudit.strictBoundaryValid()
    )
    {
        return fail("explicit elliptic-cylinder triangulation is incorrect");
    }

    Vector10 anisotropicCylinderCoefficients = Vector10::Zero();
    anisotropicCylinderCoefficients(0) = -0.005 * 0.005;
    anisotropicCylinderCoefficients(4) = 1.0;
    anisotropicCylinderCoefficients(5) = 2.0;
    const Quadric anisotropicCylinder =
        Quadric::fromUnnormalized(anisotropicCylinderCoefficients);
    spin::QuadricTriangulationSettings complementSettings;
    complementSettings.cylinderAngularResolution = 24;
    complementSettings.truncationRadius = 2.0;
    const auto cylinderInterior = spin::triangulateQuadric(
        anisotropicCylinder, domainCorners() * 0.5, complementSettings);
    const auto cylinderComplement = spin::triangulateQuadric(
        anisotropicCylinder.inverted(), domainCorners() * 0.5, complementSettings);
    const bool identicalRepresentedBoundary =
        cylinderInterior.closedBoundary.vertices.rows() ==
            cylinderComplement.closedBoundary.vertices.rows() &&
        cylinderInterior.closedBoundary.faces.rows() ==
            cylinderComplement.closedBoundary.faces.rows() &&
        (cylinderInterior.closedBoundary.vertices.array() ==
            cylinderComplement.closedBoundary.vertices.array()).all() &&
        (cylinderInterior.closedBoundary.faces.array() ==
            cylinderComplement.closedBoundary.faces.array()).all();
    if (geometricTriangleKeys(cylinderInterior) !=
            geometricTriangleKeys(cylinderComplement) ||
        !identicalRepresentedBoundary)
    {
        return fail(
            "elliptic-cylinder q and -q must use an identical binary64 boundary "
            "and geometric triangle facets"
        );
    }

    Vector10 fidgetCoefficients = Vector10::Zero();
    fidgetCoefficients(1) = 0.8;
    fidgetCoefficients(4) = 1.0;
    fidgetCoefficients(5) = 1.0;
    const Quadric fidgetCylinder =
        Quadric::fromUnnormalized(fidgetCoefficients);
    const auto fidgetSource = spin::triangulateQuadric(
        fidgetCylinder, 2.0 * domainCorners());
    const auto fidgetComplement = spin::triangulateQuadric(
        Quadric::fromUnnormalized(-fidgetCoefficients),
        2.0 * domainCorners());
    if (!igl::piecewise_constant_winding_number(
            fidgetSource.closedBoundary.faces) ||
        !igl::piecewise_constant_winding_number(
            fidgetComplement.closedBoundary.faces))
    {
        return fail(
            "elliptic-cylinder boundary and complement must both be PWN"
        );
    }

    coefficients << -0.5605562977479883, 0.2907335221898638,
        -0.02793678146101007, 0.7298376100915582,
        -0.01413811771179968, -0.0050173448614208305,
        -0.2346716529617043, 0.010029654990118242,
        -0.11140384401058069, -1.5370819321812222e-06;
    const Quadric heartReplayQuadric = Quadric::fromUnnormalized(coefficients);
    spin::QuadricTriangulationSettings heartReplaySettings;
    heartReplaySettings.icosphereSubdivisions = 3;
    heartReplaySettings.truncationRadius = 6.2757304553271505;
    const auto heartReplay = spin::triangulateQuadric(
        heartReplayQuadric, domain, heartReplaySettings);
    const auto heartReplayAudit = manufacturing::auditTriangleBoundary(
        heartReplay.closedBoundary.vertices,
        heartReplay.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    if (!heartReplayAudit.strictBoundaryValid() ||
        !igl::piecewise_constant_winding_number(heartReplay.closedBoundary.faces))
    {
        return fail("fixed Heart low-resolution source must remain a closed, consistently oriented boundary");
    }
    heartReplaySettings.icosphereSubdivisions = 6;
    const auto refinedHeartReplay = spin::triangulateQuadric(
        heartReplayQuadric, domain, heartReplaySettings);
    const auto refinedHeartReplayAudit = manufacturing::auditTriangleBoundary(
        refinedHeartReplay.closedBoundary.vertices,
        refinedHeartReplay.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    if (!refinedHeartReplayAudit.strictBoundaryValid())
    {
        return fail("fixed Heart level-6 quadric source did not pass the closed oriented-boundary audit");
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
