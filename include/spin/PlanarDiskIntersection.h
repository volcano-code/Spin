#pragma once

#include "spin/PlanarMassMoments.h"
#include "physics/RawMassMoments10.h"

#include <string>
#include <vector>

namespace spin
{

struct PlanarCircleDisk
{
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    double radius = 0.0;
};

/**
 * Local combinatorial state of a polygon--disk clipping configuration.
 *
 * `nearTransition` marks edge/vertex contact within a scale-aware numerical
 * tolerance, where central finite differences are not valid. `signature`
 * records strict interior crossings per loop for local derivative guards; it
 * is not a certificate of polygon validity.
 */
struct PlanarDiskIntersectionTopology
{
    std::string signature;
    bool nearTransition = false;
};

enum class PlanarBoundaryPieceType
{
    LineSegment,
    CircularArc
};

/** One oriented exact boundary piece; material is on its left. */
struct PlanarBoundaryPiece
{
    PlanarBoundaryPieceType type = PlanarBoundaryPieceType::LineSegment;
    Eigen::Vector2d start = Eigen::Vector2d::Zero();
    Eigen::Vector2d end = Eigen::Vector2d::Zero();
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    double radius = 0.0;
    double startAngle = 0.0;
    /** Signed sweep in radians; clockwise cavity arcs are negative. */
    double sweep = 0.0;
};

using PlanarBoundaryCycle = std::vector<PlanarBoundaryPiece>;

struct PlanarRegionDifferenceBoundary
{
    bool defined = false;
    std::string diagnostic;
    std::vector<PlanarBoundaryCycle> cycles;
};

PlanarDiskIntersectionTopology classifyPlanarDiskIntersectionTopology(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk);

/**
 * Integrate the intersection of an oriented polygonal region and a disk.
 *
 * Boundary loops use the same convention as integratePlanarBoundaryLoops:
 * material lies to the left (outer loops CCW, holes CW). This routine assumes
 * simple, non-self-intersecting, properly nested loops; it does not certify
 * polygon topology. The six returned values are density-free raw moments.
 */
PlanarRawMoments6 integratePlanarDiskIntersection(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk);

/** Exact line/circular-arc boundary of polygonal region minus a disk.
 *
 * The source loops are assumed simple, properly nested, and oriented with
 * material on the left. Tangencies and edge/vertex contacts return
 * `defined=false` because the boundary arrangement is not locally stable.
 * No curve is polygonized by this representation.
 */
PlanarRegionDifferenceBoundary planarRegionDifferenceBoundary(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk);

/** Independently integrate oriented line/arc boundary cycles. */
PlanarRawMoments6 integratePlanarBoundaryCycles(
    const std::vector<PlanarBoundaryCycle>& cycles);

/** Exact regularized positive-area intersection test for polygonal regions.
 * Loops must be simple and oriented with material on their left. Boundary-
 * only contact is not positive-area overlap. Invalid arrangements throw.
 */
bool planarRegionsHavePositiveAreaIntersection(
    const std::vector<PlanarBoundaryLoop>& first,
    const std::vector<PlanarBoundaryLoop>& second);

/**
 * Integrate the exact extrusion of oriented line/arc cycles over [z0,z1].
 * Circular arcs remain analytic; no polygonization is performed. The returned
 * ten moments are density weighted and use RawMassMoments10's ordering.
 */
physics::RawMassMoments10 integratePlanarBoundaryCyclesExtruded(
    const std::vector<PlanarBoundaryCycle>& cycles,
    double z0,
    double z1,
    double volumeDensity);

} // namespace spin
