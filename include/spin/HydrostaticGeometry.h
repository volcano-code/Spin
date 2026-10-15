#pragma once

#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

namespace spin
{

/**
 * Geometry-derived hydrostatic data for the material on the negative side of
 * a waterline plane: dot(x - waterlineOrigin, waterlineNormal) <= 0.
 *
 * The input and submerged output are strict outward-oriented closed triangle
 * boundaries.  Coplanar source faces at the waterline are deliberately
 * rejected: their ownership is ambiguous under a mesh Boolean and a caller
 * must perturb or explicitly model that degenerate physical configuration.
 */
struct HydrostaticGeometry
{
    TriangleMesh submergedBoundary;
    double displacedVolume = 0.0;
    Eigen::Vector3d centerOfBuoyancy = Eigen::Vector3d::Zero();
    double waterlineArea = 0.0;
    Eigen::Vector3d waterlineCentroid = Eigen::Vector3d::Zero();
    Eigen::Matrix2d waterlineCentralSecondMoment = Eigen::Matrix2d::Zero();
};


/**
 * A scalar waterline solve for one caller-designated, closed displacement
 * boundary.  ``waterlineSignedOffset`` is measured from ``referenceOrigin``
 * along the supplied normal.  The result only satisfies Archimedes' displaced
 * volume equation for that fixed orientation; it neither infers whether a
 * material cavity floods nor proves rotational/hydrostatic stability.
 */
struct HydrostaticEquilibrium
{
    HydrostaticGeometry geometry;
    double waterlineSignedOffset = 0.0;
    double targetDisplacedVolume = 0.0;
    double volumeResidual = 0.0;
    int iterations = 0;
};


/**
 * Intersect one closed material boundary with a finite prism representing the
 * negative water halfspace, then integrate its volume and waterline section.
 *
 * The normal need not be normalized.  The result is a geometric measurement,
 * not an optimization feasibility or physical-dynamics conclusion.
 */
HydrostaticGeometry computeHydrostaticGeometry(
    const TriangleMesh& materialBoundary,
    const Eigen::Vector3d& waterlineOrigin,
    const Eigen::Vector3d& waterlineNormal
);


/**
 * Locate a waterline for a strictly interior requested displaced volume.
 *
 * The input is deliberately named ``displacementBoundary``: callers must
 * choose a closed boundary that correctly represents their floodability and
 * sealed-cavity semantics before calling this routine.  The function does not
 * make that physical modelling decision on their behalf.
 */
HydrostaticEquilibrium findHydrostaticEquilibrium(
    const TriangleMesh& displacementBoundary,
    const Eigen::Vector3d& referenceOrigin,
    const Eigen::Vector3d& waterlineNormal,
    double targetDisplacedVolume,
    double relativeVolumeTolerance = 1e-10,
    int maximumIterations = 80
);


} // namespace spin
