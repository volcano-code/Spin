#include "physics/MassProperties.h"
#include "physics/RawMassMoments10.h"
#include "spin/PlanarDiskIntersection.h"
#include "spin/PlanarMassMoments.h"
#include "spin/PlanarSpinMomentProblem.h"
#include "spin/SpinMomentProblem.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

using Loop = spin::PlanarBoundaryLoop;

bool close(const double actual, const double expected)
{
    return std::abs(actual - expected) <= 2e-12 *
        std::max({1.0, std::abs(actual), std::abs(expected)});
}

double scaledError(const double actual, const double expected)
{
    return std::abs(actual - expected) /
        std::max(1.0, std::abs(expected));
}

int fail(const char* message)
{
    std::cerr << "[FAIL] scenario=planar-extrusion " << message << '\n';
    return EXIT_FAILURE;
}

Loop rectangle(const double x0, const double y0, const double x1, const double y1,
               const bool counterClockwise)
{
    Loop loop{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    if (!counterClockwise)
        std::reverse(loop.begin(), loop.end());
    return loop;
}

std::vector<spin::PlanarBoundaryCycle> lineCycles(const std::vector<Loop>& loops)
{
    std::vector<spin::PlanarBoundaryCycle> cycles;
    for (const auto& loop : loops)
    {
        spin::PlanarBoundaryCycle cycle;
        for (std::size_t i = 0; i < loop.size(); ++i)
        {
            spin::PlanarBoundaryPiece piece;
            piece.type = spin::PlanarBoundaryPieceType::LineSegment;
            piece.start = loop[i];
            piece.end = loop[(i+1)%loop.size()];
            cycle.push_back(std::move(piece));
        }
        cycles.push_back(std::move(cycle));
    }
    return cycles;
}

void appendBox(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces,
               const double x0, const double y0, const double x1, const double y1,
               const double z0, const double z1, const bool inward)
{
    const Eigen::Index offset = vertices.rows();
    vertices.conservativeResize(offset + 8, 3);
    vertices.block(offset, 0, 8, 3) <<
        x0, y0, z0, x1, y0, z0, x1, y1, z0, x0, y1, z0,
        x0, y0, z1, x1, y0, z1, x1, y1, z1, x0, y1, z1;
    const Eigen::Index oldFaces = faces.rows();
    faces.conservativeResize(oldFaces + 12, 3);
    faces.block(oldFaces, 0, 12, 3) <<
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5;
    faces.block(oldFaces, 0, 12, 3).array() += offset;
    if (inward)
        faces.block(oldFaces, 0, 12, 3).col(0).swap(faces.block(oldFaces, 0, 12, 3).col(2));
}

void appendCylinder(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces,
                    const Eigen::Vector2d& center, const double radius,
                    const double z0, const double z1,
                    const std::size_t segments, const bool inward)
{
    const Eigen::Index offset = vertices.rows();
    const Eigen::Index count = static_cast<Eigen::Index>(2*segments + 2);
    vertices.conservativeResize(offset + count, 3);
    vertices.row(offset) << center.x(), center.y(), z0;
    vertices.row(offset+1) << center.x(), center.y(), z1;
    const double tau = 2.0*std::acos(-1.0);
    for (std::size_t i = 0; i < segments; ++i)
    {
        const double angle = tau*static_cast<double>(i)/static_cast<double>(segments);
        const double x = center.x() + radius*std::cos(angle);
        const double y = center.y() + radius*std::sin(angle);
        vertices.row(offset+2+static_cast<Eigen::Index>(i)) << x, y, z0;
        vertices.row(offset+2+static_cast<Eigen::Index>(segments+i)) << x, y, z1;
    }
    const Eigen::Index oldFaces = faces.rows();
    faces.conservativeResize(oldFaces + static_cast<Eigen::Index>(4*segments), 3);
    Eigen::Index row = oldFaces;
    for (std::size_t i = 0; i < segments; ++i)
    {
        const Eigen::Index current = static_cast<Eigen::Index>(i);
        const Eigen::Index next = static_cast<Eigen::Index>((i+1)%segments);
        const Eigen::Index bottomCurrent = offset+2+current;
        const Eigen::Index bottomNext = offset+2+next;
        const Eigen::Index topCurrent = offset+2+static_cast<Eigen::Index>(segments)+current;
        const Eigen::Index topNext = offset+2+static_cast<Eigen::Index>(segments)+next;
        faces.row(row++) << offset, bottomNext, bottomCurrent;
        faces.row(row++) << offset+1, topCurrent, topNext;
        faces.row(row++) << bottomCurrent, bottomNext, topNext;
        faces.row(row++) << bottomCurrent, topNext, topCurrent;
    }
    if (inward)
        faces.block(oldFaces, 0, row-oldFaces, 3).col(0).swap(
            faces.block(oldFaces, 0, row-oldFaces, 3).col(2));
}

} // namespace


int main()
{
    // A translated rectangular ring avoids relying on symmetry and exercises
    // all six planar moments, including Qxy.  The 3D inner box is inward
    // oriented, so it is an actual volume cavity rather than an omitted mesh.
    const double z0 = -0.7, z1 = 0.8;
    const double height = z1 - z0;
    const double density3d = 2.4;
    const std::vector<Loop> region{
        rectangle(-2.0, -1.0, 3.0, 4.0, true),
        rectangle(-0.75, 0.25, 1.25, 2.5, false),
    };
    const spin::PlanarMassMoments6 planar =
        spin::integratePlanarBoundaryLoops(region).weighted(density3d * height);

    Eigen::MatrixXd vertices(0, 3);
    Eigen::MatrixXi faces(0, 3);
    appendBox(vertices, faces, -2.0, -1.0, 3.0, 4.0, z0, z1, false);
    appendBox(vertices, faces, -0.75, 0.25, 1.25, 2.5, z0, z1, true);
    const MassMoments extruded = computeSurfaceMassMoments(vertices, faces, density3d);
    if (!extruded.isFinite() || !(extruded.mass > 0.0))
        return fail("extruded ring has invalid volume moments");

    if (!close(extruded.mass, planar.mass) ||
        !close(extruded.firstMoment.x(), planar.firstX) ||
        !close(extruded.firstMoment.y(), planar.firstY) ||
        !close(extruded.secondMoment(0, 0), planar.secondXX) ||
        !close(extruded.secondMoment(1, 1), planar.secondYY) ||
        !close(extruded.secondMoment(0, 1), planar.secondXY))
    {
        return fail("planar and independently integrated 3D extrusion moments differ");
    }
    const double momentErrors[] = {
        scaledError(extruded.mass, planar.mass),
        scaledError(extruded.firstMoment.x(), planar.firstX),
        scaledError(extruded.firstMoment.y(), planar.firstY),
        scaledError(extruded.secondMoment(0, 0), planar.secondXX),
        scaledError(extruded.secondMoment(1, 1), planar.secondYY),
        scaledError(extruded.secondMoment(0, 1), planar.secondXY),
    };
    const double maxMomentError = *std::max_element(
        std::begin(momentErrors), std::end(momentErrors));

    spin::SpinMomentProblem::Vector10 raw = spin::SpinMomentProblem::Vector10::Zero();
    raw(0) = extruded.mass;
    raw.segment<3>(1) = extruded.firstMoment;
    raw.segment<3>(4) = extruded.secondMoment.diagonal();
    raw(7) = extruded.secondMoment(1, 2);
    raw(8) = extruded.secondMoment(0, 2);
    raw(9) = extruded.secondMoment(0, 1);
    const auto planarEvaluation = spin::PlanarSpinMomentProblem::fidgetSpinner().evaluate(planar);
    const auto spatialEvaluation = spin::SpinMomentProblem::fidgetSpinner().evaluate(raw);
    if (!close(planarEvaluation.objective, spatialEvaluation.objective) ||
        !close(planarEvaluation.constraints.x(), spatialEvaluation.constraints(0)) ||
        !close(planarEvaluation.constraints.y(), spatialEvaluation.constraints(1)))
    {
        return fail("2D and extruded 3D fidget objective or balance constraints differ");
    }
    const double objectiveErrors[] = {
        scaledError(planarEvaluation.objective, spatialEvaluation.objective),
        scaledError(planarEvaluation.constraints.x(), spatialEvaluation.constraints(0)),
        scaledError(planarEvaluation.constraints.y(), spatialEvaluation.constraints(1)),
    };
    const double maxObjectiveError = *std::max_element(
        std::begin(objectiveErrors), std::end(objectiveErrors));

    // Exact 2D fidget output is represented by line/circular-arc cycles. Its
    // extrusion is an analytic 3D material object over [z0,z1], not a
    // polygonized approximation. Compare all ten volume moments against the
    // independent polygon-minus-disk integrator followed by closed-form z
    // integration.
    const std::vector<Loop> square{rectangle(-3.0, -3.0, 3.0, 3.0, true)};
    const spin::PlanarCircleDisk disk{Eigen::Vector2d(0.2, -0.1), 1.0};
    const auto difference = spin::planarRegionDifferenceBoundary(square, disk);
    if (!difference.defined || difference.cycles.size() != 2)
        return fail("exact disk difference did not produce outer and cavity cycles");
    bool hasArc = false;
    for (const auto& cycle : difference.cycles)
        for (const auto& piece : cycle)
            hasArc = hasArc || piece.type == spin::PlanarBoundaryPieceType::CircularArc;
    if (!hasArc)
        return fail("exact 2D material boundary contains no circular arc");
    const auto extrudedExact = spin::integratePlanarBoundaryCyclesExtruded(
        difference.cycles, z0, z1, density3d);
    const auto planarOuter = spin::integratePlanarBoundaryLoops(square);
    const auto planarRemoved = spin::integratePlanarDiskIntersection(square, disk);
    const double planarRaw[] = {
        planarOuter.area-planarRemoved.area,
        planarOuter.firstX-planarRemoved.firstX,
        planarOuter.firstY-planarRemoved.firstY,
        planarOuter.secondXX-planarRemoved.secondXX,
        planarOuter.secondYY-planarRemoved.secondYY,
        planarOuter.secondXY-planarRemoved.secondXY,
    };
    const double zFirst = (z1*z1-z0*z0)/2.0;
    const double zSecond = (z1*z1*z1-z0*z0*z0)/3.0;
    const double expected3d[] = {
        density3d*height*planarRaw[0], density3d*height*planarRaw[1],
        density3d*height*planarRaw[2], density3d*zFirst*planarRaw[0],
        density3d*height*planarRaw[3], density3d*height*planarRaw[4],
        density3d*zSecond*planarRaw[0], density3d*zFirst*planarRaw[2],
        density3d*zFirst*planarRaw[1], density3d*height*planarRaw[5],
    };
    for (Eigen::Index i = 0; i < physics::RawMassMoments10::Count; ++i)
        if (!close(extrudedExact.values()(i), expected3d[i]))
            return fail("exact extruded arc-boundary moments disagree with independent 3D formulas");

    // A separately built closed triangle mesh is only a convergence check,
    // never the exact representation: the cylinder is inscribed and its
    // polygonization error must stay within the registered 1e-6 envelope.
    Eigen::MatrixXd curvedVertices(0, 3);
    Eigen::MatrixXi curvedFaces(0, 3);
    appendBox(curvedVertices, curvedFaces, -3.0, -3.0, 3.0, 3.0, z0, z1, false);
    appendCylinder(curvedVertices, curvedFaces, disk.center, disk.radius,
                   z0, z1, 4096, true);
    const auto curvedMeshMoments = physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(curvedVertices, curvedFaces, density3d));
    double maxCurvedMeshScaledError = 0.0;
    for (Eigen::Index i = 0; i < physics::RawMassMoments10::Count; ++i)
    {
        const double error = scaledError(
            curvedMeshMoments.values()(i), extrudedExact.values()(i));
        maxCurvedMeshScaledError = std::max(maxCurvedMeshScaledError, error);
        if (error > 1e-6)
            return fail("high-resolution mesh convergence check exceeds 1e-6; exact arc result remains authoritative");
    }

    const double pi = std::acos(-1.0);
    const auto fixedPatch = rectangle(pi/4.0-0.1, -0.1, pi/4.0+0.1, 0.1, true);
    const auto fixedCycles = lineCycles(std::vector<Loop>{fixedPatch});
    const auto fixedLower = spin::integratePlanarBoundaryCyclesExtruded(
        fixedCycles, -0.1, 0.0, 500.0);
    const auto fixedUpper = spin::integratePlanarBoundaryCyclesExtruded(
        fixedCycles, 0.0, 0.1, 500.0);
    const auto fixedLayerSum = fixedLower + fixedUpper;
    Eigen::MatrixXd fixedLayerVertices(0, 3);
    Eigen::MatrixXi fixedLayerFaces(0, 3);
    appendBox(fixedLayerVertices, fixedLayerFaces,
        pi/4.0-0.1, -0.1, pi/4.0+0.1, 0.1, -0.1, 0.0, false);
    appendBox(fixedLayerVertices, fixedLayerFaces,
        pi/4.0-0.1, -0.1, pi/4.0+0.1, 0.1, 0.0, 0.1, false);
    const auto fixedLayerMeshMoments = physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(fixedLayerVertices, fixedLayerFaces, 500.0));
    for (Eigen::Index i = 0; i < physics::RawMassMoments10::Count; ++i)
        if (!close(fixedLayerMeshMoments.values()(i), fixedLayerSum.values()(i)))
            return fail("exact fixed-layer extrusion moments differ from independent prism meshes");

    std::cout << std::setprecision(17)
              << "FASTER_PLANAR_COMPONENT_JSON={\"schema\":\"faster-planar-component-result-v1\","
              << "\"classification\":\"SYNTHETIC_REGRESSION\","
              << "\"component\":\"extrusion\","
              << "\"planar_boundary_loop_count\":2,\"surface_shell_count\":2,"
              << "\"extrusion_height\":" << height << ','
              << "\"density\":" << density3d << ','
              << "\"planar_mass\":" << planar.mass << ','
              << "\"extruded_mass\":" << extruded.mass << ','
              << "\"maximum_scaled_absolute_moment_error\":" << maxMomentError << ','
              << "\"exact_arc_extrusion_mesh_convergence_error\":" << maxCurvedMeshScaledError << ','
              << "\"maximum_scaled_objective_constraint_error\":" << maxObjectiveError << ','
              << "\"planar_objective\":" << planarEvaluation.objective << ','
              << "\"spatial_objective\":" << spatialEvaluation.objective << ','
              << "\"planar_constraints\":[" << planarEvaluation.constraints.x() << ','
              << planarEvaluation.constraints.y() << "],"
              << "\"spatial_constraints\":[" << spatialEvaluation.constraints(0) << ','
              << spatialEvaluation.constraints(1) << "],"
              << "\"scaled_absolute_error_tolerance\":2e-12}\n";
    std::cout << "[PASS] planar ring moments equal independently integrated 3D extrusion\n";
    std::cout << "[PASS] planar and 3D fidget objectives/constraints agree\n";
    return EXIT_SUCCESS;
}
