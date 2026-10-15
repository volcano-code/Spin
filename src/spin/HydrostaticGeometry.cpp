#include "spin/HydrostaticGeometry.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "physics/MassProperties.h"

#include <igl/MeshBooleanType.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <CGAL/Gmpq.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace
{


void requireStrictPositiveBoundary(const spin::TriangleMesh& mesh, const char* name)
{
    const manufacturing::BoundaryAudit audit = manufacturing::auditTriangleBoundary(
        mesh.vertices, mesh.faces, manufacturing::BoundaryOrientation::OutwardPositive
    );
    if (!audit.strictBoundaryValid() || !std::isfinite(audit.signedVolume) ||
        !(audit.signedVolume > 0.0))
    {
        throw std::invalid_argument(std::string(name) +
            " must be a strict positive, closed, outward-oriented triangle boundary");
    }
}


std::pair<Eigen::Vector3d, Eigen::Vector3d> planeBasis(const Eigen::Vector3d& normal)
{
    const Eigen::Vector3d seed = std::abs(normal.x()) < 0.8
        ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
    const Eigen::Vector3d u = (seed - normal.dot(seed) * normal).normalized();
    const Eigen::Vector3d v = normal.cross(u);
    if (!u.allFinite() || !v.allFinite())
        throw std::invalid_argument("waterline normal does not define a finite plane basis");
    return {u, v};
}


spin::TriangleMesh negativeHalfspacePrism(
    const spin::TriangleMesh& material,
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& u,
    const Eigen::Vector3d& v
)
{
    double radius = 0.0;
    for (Eigen::Index row = 0; row < material.vertices.rows(); ++row)
        radius = std::max(radius, (material.vertices.row(row).transpose() - origin).norm());
    if (!std::isfinite(radius))
        throw std::invalid_argument("waterline/material scale is non-finite");
    const double extent = std::nextafter(std::max(1.0, 4.0 * radius),
                                         std::numeric_limits<double>::infinity());
    const double bottom = -2.0 * extent;
    spin::TriangleMesh prism;
    prism.vertices.resize(8, 3);
    const std::array<Eigen::Vector3d, 8> local = {{
        {-extent, -extent, bottom}, { extent, -extent, bottom},
        { extent,  extent, bottom}, {-extent,  extent, bottom},
        {-extent, -extent, 0.0},    { extent, -extent, 0.0},
        { extent,  extent, 0.0},    {-extent,  extent, 0.0},
    }};
    for (int index = 0; index < 8; ++index)
        prism.vertices.row(index) = (origin + u * local[index].x() +
                                     v * local[index].y() + normal * local[index].z()).transpose();
    prism.faces.resize(12, 3);
    prism.faces << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
                    0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5,
                    2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7;
    return prism;
}

// Signed origin-tetrahedron integration, exact on the represented boundary.
// Unlike a floating reference-vertex mean and sum, this is independent of
// vertex/face enumeration and equivalent triangulations of the same surface.
std::pair<double, Eigen::Vector3d> displacedVolumeAndCentroid(
    const spin::TriangleMesh& mesh)
{
    using Q = CGAL::Gmpq;
    Q sixVolume(0);
    std::array<Q, 3> twentyFourFirst = {Q(0), Q(0), Q(0)};
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        std::array<std::array<Q, 3>, 3> p;
        for (int corner = 0; corner < 3; ++corner)
            for (int axis = 0; axis < 3; ++axis)
                p[corner][axis] = Q(mesh.vertices(mesh.faces(face, corner), axis));
        const auto& a = p[0]; const auto& b = p[1]; const auto& c = p[2];
        const Q determinant = a[0]*(b[1]*c[2]-b[2]*c[1]) +
            a[1]*(b[2]*c[0]-b[0]*c[2]) + a[2]*(b[0]*c[1]-b[1]*c[0]);
        sixVolume += determinant;
        for (int axis = 0; axis < 3; ++axis)
            twentyFourFirst[axis] += determinant*(a[axis]+b[axis]+c[axis]);
    }
    if (!(sixVolume > 0))
        throw std::runtime_error("submerged hydrostatic volume is not positive");
    const double volume = CGAL::to_double(sixVolume / 6);
    Eigen::Vector3d centroid;
    for (int axis = 0; axis < 3; ++axis)
        centroid[axis] = CGAL::to_double(twentyFourFirst[axis] / (4*sixVolume));
    if (!(volume > 0.0) || !std::isfinite(volume) || !centroid.allFinite())
        throw std::runtime_error("submerged hydrostatic volume/centroid is not positive and finite");
    return {volume, centroid};
}


} // namespace


namespace spin
{


HydrostaticGeometry computeHydrostaticGeometry(
    const TriangleMesh& materialBoundary,
    const Eigen::Vector3d& waterlineOrigin,
    const Eigen::Vector3d& waterlineNormal
)
{
    requireStrictPositiveBoundary(materialBoundary, "hydrostatic material boundary");
    if (!waterlineOrigin.allFinite() || !waterlineNormal.allFinite() ||
        !(waterlineNormal.norm() > 0.0))
    {
        throw std::invalid_argument("waterline origin and nonzero normal must be finite");
    }
    const Eigen::Vector3d normal = waterlineNormal.normalized();
    const auto [u, v] = planeBasis(normal);
    const Eigen::VectorXd signedDistances =
        (materialBoundary.vertices.rowwise() - waterlineOrigin.transpose()) * normal;
    const double scale = std::max(1.0, (materialBoundary.vertices.rowwise() -
                                         waterlineOrigin.transpose()).rowwise().norm().maxCoeff());
    const double tolerance = 1e-10 * scale;
    if (!(signedDistances.minCoeff() < -tolerance && signedDistances.maxCoeff() > tolerance))
        throw std::invalid_argument("waterline must strictly intersect the material boundary");
    for (Eigen::Index face = 0; face < materialBoundary.faces.rows(); ++face)
    {
        const int a = materialBoundary.faces(face, 0), b = materialBoundary.faces(face, 1),
                  c = materialBoundary.faces(face, 2);
        if (std::abs(signedDistances(a)) <= tolerance &&
            std::abs(signedDistances(b)) <= tolerance &&
            std::abs(signedDistances(c)) <= tolerance)
        {
            throw std::invalid_argument("source face is coplanar with waterline");
        }
    }

    const TriangleMesh prism = negativeHalfspacePrism(
        materialBoundary, waterlineOrigin, normal, u, v
    );
    HydrostaticGeometry result;
    Eigen::VectorXi birthFaces;
    if (!igl::copyleft::cgal::mesh_boolean(
            materialBoundary.vertices, materialBoundary.faces, prism.vertices, prism.faces,
            igl::MESH_BOOLEAN_TYPE_INTERSECT, result.submergedBoundary.vertices,
            result.submergedBoundary.faces, birthFaces))
    {
        throw std::runtime_error("hydrostatic mesh Boolean failed");
    }
    requireStrictPositiveBoundary(result.submergedBoundary, "submerged hydrostatic boundary");
    const auto [volume, centroid] = displacedVolumeAndCentroid(result.submergedBoundary);
    result.displacedVolume = volume;
    result.centerOfBuoyancy = centroid;

    // Boolean triangle order is not stable. Integrate the represented projected
    // coordinates exactly, including centralization, then convert only the final
    // moments. Long-double accumulation made an exactly symmetric section's
    // cross moment depend on face order. This does not snap values or claim that
    // floating Boolean coordinates equal the original exact intersection.
    using Q = CGAL::Gmpq;
    Q area(0), firstU(0), firstV(0);
    Q secondUU(0), secondVV(0), secondUV(0);
    int capFaceCount = 0;
    for (Eigen::Index face = 0; face < result.submergedBoundary.faces.rows(); ++face)
    {
        const int ia = result.submergedBoundary.faces(face, 0);
        const int ib = result.submergedBoundary.faces(face, 1);
        const int ic = result.submergedBoundary.faces(face, 2);
        const Eigen::Vector3d a = result.submergedBoundary.vertices.row(ia).transpose();
        const Eigen::Vector3d b = result.submergedBoundary.vertices.row(ib).transpose();
        const Eigen::Vector3d c = result.submergedBoundary.vertices.row(ic).transpose();
        if (std::abs((a - waterlineOrigin).dot(normal)) > tolerance ||
            std::abs((b - waterlineOrigin).dot(normal)) > tolerance ||
            std::abs((c - waterlineOrigin).dot(normal)) > tolerance)
            continue;
        const Q x0((a - waterlineOrigin).dot(u)), y0((a - waterlineOrigin).dot(v));
        const Q x1((b - waterlineOrigin).dot(u)), y1((b - waterlineOrigin).dot(v));
        const Q x2((c - waterlineOrigin).dot(u)), y2((c - waterlineOrigin).dot(v));
        const Q signedArea = ((x1-x0)*(y2-y0) - (y1-y0)*(x2-x0)) / 2;
        if (!(signedArea > 0))
            throw std::runtime_error("hydrostatic waterline cap lacks outward orientation");
        ++capFaceCount;
        area += signedArea;
        firstU += signedArea * (x0+x1+x2) / 3;
        firstV += signedArea * (y0+y1+y2) / 3;
        secondUU += signedArea * (x0*x0+x1*x1+x2*x2+x0*x1+x1*x2+x2*x0) / 6;
        secondVV += signedArea * (y0*y0+y1*y1+y2*y2+y0*y1+y1*y2+y2*y0) / 6;
        secondUV += signedArea * (2*(x0*y0+x1*y1+x2*y2) +
            x0*y1+x1*y0+x1*y2+x2*y1+x2*y0+x0*y2) / 12;
    }
    if (capFaceCount == 0 || !(area > 0) || !std::isfinite(CGAL::to_double(area)))
        throw std::runtime_error("hydrostatic Boolean has no positive waterline cap");
    const Q centerU = firstU / area, centerV = firstV / area;
    result.waterlineArea = CGAL::to_double(area);
    result.waterlineCentroid = waterlineOrigin + CGAL::to_double(centerU) * u +
        CGAL::to_double(centerV) * v;
    result.waterlineCentralSecondMoment <<
        CGAL::to_double(secondUU - area*centerU*centerU),
        CGAL::to_double(secondUV - area*centerU*centerV),
        CGAL::to_double(secondUV - area*centerU*centerV),
        CGAL::to_double(secondVV - area*centerV*centerV);
    if (!result.centerOfBuoyancy.allFinite() || !result.waterlineCentroid.allFinite() ||
        !result.waterlineCentralSecondMoment.allFinite())
    {
        throw std::runtime_error("hydrostatic geometric moments are non-finite");
    }
    return result;
}


HydrostaticEquilibrium findHydrostaticEquilibrium(
    const TriangleMesh& displacementBoundary,
    const Eigen::Vector3d& referenceOrigin,
    const Eigen::Vector3d& waterlineNormal,
    const double targetDisplacedVolume,
    const double relativeVolumeTolerance,
    const int maximumIterations
)
{
    requireStrictPositiveBoundary(displacementBoundary, "hydrostatic displacement boundary");
    if (!referenceOrigin.allFinite() || !waterlineNormal.allFinite() ||
        !(waterlineNormal.norm() > 0.0) || !std::isfinite(targetDisplacedVolume) ||
        !std::isfinite(relativeVolumeTolerance) || !(relativeVolumeTolerance > 0.0) ||
        maximumIterations <= 0)
    {
        throw std::invalid_argument(
            "hydrostatic equilibrium requires finite reference, nonzero normal, positive finite target tolerance, and iteration count"
        );
    }
    const Eigen::Vector3d normal = waterlineNormal.normalized();
    const MassMoments totalMoments = computeSurfaceMassMoments(
        displacementBoundary.vertices, displacementBoundary.faces, 1.0
    );
    if (!totalMoments.isFinite() || !(totalMoments.mass > 0.0))
        throw std::runtime_error("hydrostatic displacement boundary has nonpositive total volume");
    const double volumeTolerance = relativeVolumeTolerance * totalMoments.mass;
    if (!(targetDisplacedVolume > volumeTolerance &&
          targetDisplacedVolume < totalMoments.mass - volumeTolerance))
    {
        throw std::invalid_argument(
            "target displaced volume must be strictly interior to the closed displacement boundary volume"
        );
    }

    const Eigen::VectorXd projections =
        (displacementBoundary.vertices.rowwise() - referenceOrigin.transpose()) * normal;
    const double minimumProjection = projections.minCoeff();
    const double maximumProjection = projections.maxCoeff();
    const double projectionRange = maximumProjection - minimumProjection;
    if (!std::isfinite(projectionRange) || !(projectionRange > 0.0))
        throw std::invalid_argument("hydrostatic displacement boundary has zero waterline extent");
    const double geometricScale = std::max(1.0, (displacementBoundary.vertices.rowwise() -
        referenceOrigin.transpose()).rowwise().norm().maxCoeff());
    // Keep all trial waterlines outside computeHydrostaticGeometry's explicit
    // coplanarity neighbourhood.  This is a geometric tolerance only; the
    // residual acceptance below remains a volume criterion.
    const double inset = std::max(
        2e-10 * geometricScale,
        std::numeric_limits<double>::epsilon() * projectionRange
    );
    if (!(minimumProjection + inset < maximumProjection - inset))
        throw std::invalid_argument("hydrostatic displacement boundary is too thin for a strict waterline solve");

    double lower = minimumProjection + inset;
    double upper = maximumProjection - inset;
    HydrostaticEquilibrium best;
    bool haveBest = false;
    for (int iteration = 1; iteration <= maximumIterations; ++iteration)
    {
        const double offset = 0.5 * (lower + upper);
        HydrostaticGeometry geometry = computeHydrostaticGeometry(
            displacementBoundary, referenceOrigin + offset * normal, normal
        );
        const double residual = geometry.displacedVolume - targetDisplacedVolume;
        if (!haveBest || std::abs(residual) < std::abs(best.volumeResidual))
        {
            best.geometry = std::move(geometry);
            best.waterlineSignedOffset = offset;
            best.targetDisplacedVolume = targetDisplacedVolume;
            best.volumeResidual = residual;
            best.iterations = iteration;
            haveBest = true;
        }
        if (std::abs(residual) <= volumeTolerance)
            return best;
        if (residual < 0.0)
            lower = offset;
        else
            upper = offset;
    }
    if (!haveBest)
        throw std::runtime_error("hydrostatic equilibrium produced no valid waterline trial");
    throw std::runtime_error("hydrostatic equilibrium volume residual did not converge");
}


} // namespace spin
