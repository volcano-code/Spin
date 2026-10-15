#include "spin/PlanarDiskIntersection.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

using Loop = spin::PlanarBoundaryLoop;

double maximumScaledMomentError = 0.0;
std::size_t momentCases = 0;

bool close(double actual, double expected)
{
    return std::abs(actual - expected) <=
        2e-12 * std::max(1.0, std::abs(expected));
}

Loop rectangle(double x0, double y0, double x1, double y1, bool ccw)
{
    Loop points{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    if (!ccw)
        std::reverse(points.begin(), points.end());
    return points;
}

bool matches(const spin::PlanarRawMoments6& value,
             double area, double firstX, double firstY,
             double secondXX, double secondYY, double secondXY)
{
    const double actual[] = {value.area, value.firstX, value.firstY,
                             value.secondXX, value.secondYY, value.secondXY};
    const double expected[] = {area, firstX, firstY, secondXX, secondYY, secondXY};
    for (std::size_t i = 0; i < 6; ++i)
        maximumScaledMomentError = std::max(
            maximumScaledMomentError,
            std::abs(actual[i] - expected[i]) / std::max(1.0, std::abs(expected[i])));
    ++momentCases;
    return close(value.area, area) && close(value.firstX, firstX) &&
        close(value.firstY, firstY) && close(value.secondXX, secondXX) &&
        close(value.secondYY, secondYY) && close(value.secondXY, secondXY);
}

int fail(const char* message)
{
    std::cerr << message << '\n';
    return 1;
}

} // namespace

int main()
{
    // A translated RIGHT semicircle has an open arc and a vertical closing
    // chord. Closed-region translation cannot be applied to the arc alone.
    spin::PlanarBoundaryPiece arc;
    arc.type = spin::PlanarBoundaryPieceType::CircularArc;
    arc.center = Eigen::Vector2d(.25, .5);
    arc.radius = 1.;
    arc.startAngle = -.5*M_PI;
    arc.sweep = M_PI;
    arc.start = Eigen::Vector2d(.25, -.5);
    arc.end = Eigen::Vector2d(.25, 1.5);
    spin::PlanarBoundaryPiece chord;
    chord.start = arc.end;
    chord.end = arc.start;
    const auto semicircleBoundary = spin::integratePlanarBoundaryCycles({{arc, chord}});
    if (!matches(semicircleBoundary, .5*M_PI, .125*M_PI+2./3., .25*M_PI,
                 5.*M_PI/32.+1./3., .25*M_PI, M_PI/16.+1./3.))
        return fail("translated open-arc boundary integration omitted endpoint terms");

    const spin::PlanarCircleDisk disk{Eigen::Vector2d(1.0, 0.0), 1.0};
    const auto fullDisk = spin::integratePlanarDiskIntersection(
        {rectangle(-4.0, -4.0, 4.0, 4.0, true)}, disk);
    if (!matches(fullDisk, M_PI, M_PI, 0.0, 1.25*M_PI, 0.25*M_PI, 0.0))
        return fail("whole circular disk moments are incorrect");

    const auto translatedDisk = spin::integratePlanarDiskIntersection(
        {rectangle(-10.0, -10.0, 10.0, 10.0, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d(-2.0, 3.0), 2.0});
    if (!matches(translatedDisk, 4.0*M_PI, -8.0*M_PI, 12.0*M_PI,
                 20.0*M_PI, 40.0*M_PI, -24.0*M_PI))
        return fail("translated, scaled disk moments violate affine formulas");

    const auto halfDisk = spin::integratePlanarDiskIntersection(
        {rectangle(0.0, 0.0, 2.0, 2.0, true)}, disk);
    if (!matches(halfDisk, 0.5*M_PI, 0.5*M_PI, 2.0/3.0,
                 0.625*M_PI, 0.125*M_PI, 2.0/3.0))
        return fail("polygon-clipped semicircle moments are incorrect");

    const auto containedPolygon = spin::integratePlanarDiskIntersection(
        {rectangle(-0.5, -0.5, 0.5, 0.5, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 2.0});
    if (!matches(containedPolygon, 1.0, 0.0, 0.0, 1.0/12.0, 1.0/12.0, 0.0))
        return fail("polygon fully contained in a disk changed its moments");

    const auto excludedByHole = spin::integratePlanarDiskIntersection(
        {rectangle(-3.0, -3.0, 3.0, 3.0, true),
         rectangle(-0.5, -0.5, 0.5, 0.5, false)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 0.25});
    if (!matches(excludedByHole, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        return fail("disk wholly inside a polygon hole contributed material");

    const auto diskMinusHole = spin::integratePlanarDiskIntersection(
        {rectangle(-3.0, -3.0, 3.0, 3.0, true),
         rectangle(-0.5, -0.5, 0.5, 0.5, false)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 1.0});
    if (!matches(diskMinusHole, M_PI-1.0, 0.0, 0.0,
                 0.25*M_PI-1.0/12.0, 0.25*M_PI-1.0/12.0, 0.0))
        return fail("disk intersecting a polygon hole has incorrect signed moments");

    const auto disjoint = spin::integratePlanarDiskIntersection(
        {rectangle(-1.0, -1.0, 1.0, 1.0, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d(4.0, 0.0), 1.0});
    if (!matches(disjoint, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        return fail("disjoint disk returned nonzero intersection moments");

    const auto tangent = spin::integratePlanarDiskIntersection(
        {rectangle(1.0, -1.0, 2.0, 1.0, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 1.0});
    if (!matches(tangent, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        return fail("zero-area tangent contact was reported as material");

    const auto tangentTopology = spin::classifyPlanarDiskIntersectionTopology(
        {rectangle(1.0, -1.0, 2.0, 1.0, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 1.0});
    if (!tangentTopology.nearTransition)
        return fail("tangent polygon/disk contact was not classified as a topology transition");

    const auto regularTopology = spin::classifyPlanarDiskIntersectionTopology(
        {rectangle(-2.0, -2.0, 2.0, 2.0, true)},
        spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), 1.0});
    if (regularTopology.nearTransition || regularTopology.signature.empty())
        return fail("regular polygon/disk clipping was misclassified as a topology transition");

    try
    {
        (void)spin::integratePlanarDiskIntersection(
            {rectangle(-1.0, -1.0, 1.0, 1.0, true)},
            spin::PlanarCircleDisk{Eigen::Vector2d::Zero(), -1.0});
        return fail("negative-radius disk was accepted");
    }
    catch (const std::invalid_argument&)
    {
    }

    std::cout << std::setprecision(17)
              << "FASTER_PLANAR_COMPONENT_JSON={\"schema\":\"faster-planar-component-result-v1\","
              << "\"classification\":\"SYNTHETIC_REGRESSION\","
              << "\"component\":\"polygon_disk\","
              << "\"analytic_moment_cases\":" << momentCases << ','
              << "\"maximum_scaled_absolute_moment_error\":" << maximumScaledMomentError << ','
              << "\"scaled_absolute_error_tolerance\":2e-12,"
              << "\"covers_holes\":true,\"covers_empty_and_tangent_intersections\":true,"
              << "\"tangent_topology_transition_detected\":true}\n";

    return 0;
}
