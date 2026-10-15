#pragma once

#include <Eigen/Core>

#include <vector>

namespace spin
{


using PlanarBoundaryLoop = std::vector<Eigen::Vector2d>;


struct PlanarMassMoments6;


/**
 * Density-free geometric raw moments of an oriented planar region.
 *
 * The order is `(A, integral(x), integral(y), integral(x^2),
 * integral(y^2), integral(x*y))`. Boundary loops are implicitly closed and
 * must be oriented with material on their left (outer loops counter-clockwise,
 * holes clockwise). The routine checks finite coordinates, nonzero loop area,
 * and positive total area; it does not certify simplicity, containment, or
 * nonintersection of the loops.
 */
struct PlanarRawMoments6
{
    double area = 0.0;
    double firstX = 0.0;
    double firstY = 0.0;
    double secondXX = 0.0;
    double secondYY = 0.0;
    double secondXY = 0.0;

    PlanarMassMoments6 weighted(double arealDensity) const;
};


/** Areal-density-weighted moments; the first entry is mass, not area. */
struct PlanarMassMoments6
{
    double mass = 0.0;
    double firstX = 0.0;
    double firstY = 0.0;
    double secondXX = 0.0;
    double secondYY = 0.0;
    double secondXY = 0.0;

    PlanarMassMoments6& operator+=(const PlanarMassMoments6& other);
    Eigen::Matrix<double, 6, 1> vector() const;
};


/** Integrate all oriented loops using signed polygon boundary formulas. */
PlanarRawMoments6 integratePlanarBoundaryLoops(
    const std::vector<PlanarBoundaryLoop>& loops
);


} // namespace spin
