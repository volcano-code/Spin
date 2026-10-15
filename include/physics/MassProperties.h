#pragma once

#include <Eigen/Core>


/**
 * Signed mass properties for a uniform-density solid.
 *
 * firstMoment stores Q = integral rho * x dV. The center of mass is derived
 * from Q / mass; it is intentionally not stored as an independent state.
 */
struct MassProperties
{
    double mass;
    Eigen::Vector3d firstMoment;

    MassProperties();

    MassProperties(
        double massValue,
        const Eigen::Vector3d& firstMomentValue
    );

    bool isFinite() const;

    bool tryCenterOfMass(
        Eigen::Vector3d& centerOfMass,
        double minimumAbsoluteMass = 0.0
    ) const;

    MassProperties& operator+=(
        const MassProperties& other
    );

    MassProperties& operator-=(
        const MassProperties& other
    );
};


/**
 * Uniform-density volume moments through degree two.
 *
 * secondMoment stores the raw tensor S = integral rho * x*x^T dV.  The
 * inertia tensor about the center of mass is derived from S with the parallel
 * axis theorem; it is intentionally kept separate from MassProperties so the
 * balance solver's established zero/first-moment API remains unchanged.
 */
struct MassMoments
{
    double mass;
    Eigen::Vector3d firstMoment;
    Eigen::Matrix3d secondMoment;

    MassMoments();

    bool isFinite() const;
    bool tryCenterOfMass(
        Eigen::Vector3d& centerOfMass,
        double minimumAbsoluteMass = 0.0
    ) const;
    Eigen::Matrix3d inertiaAboutCenter(
        double minimumPositiveMass = 0.0
    ) const;

    MassMoments& operator+=(const MassMoments& other);
    MassMoments& operator-=(const MassMoments& other);
};


MassMoments operator+(MassMoments first, const MassMoments& second);
MassMoments operator-(MassMoments first, const MassMoments& second);


MassProperties operator+(
    MassProperties first,
    const MassProperties& second
);


MassProperties operator-(
    MassProperties first,
    const MassProperties& second
);


/**
 * Compute signed mass and first moment from a consistently oriented closed
 * triangle boundary. Outward faces contribute positive mass; inward faces
 * contribute negative mass, as required for cavity boundaries.
 */
MassProperties computeSurfaceMassProperties(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density = 1.0
);


/** Integrate mass, first moment, and raw second moment of an oriented shell. */
MassMoments computeSurfaceMassMoments(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density = 1.0
);

/** Same signed-tetrahedron integral, accumulating both input solids and their
 * difference in extended precision before one binary64 output conversion.
 * Not an exact arithmetic or boundary-validity certificate. The legacy
 * single-surface integrator and its binary64 arithmetic stay unchanged. */
MassMoments computeSurfaceDifferenceMassMomentsExtended(
    const Eigen::MatrixXd& outerVertices, const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices, const Eigen::MatrixXi& innerFaces,
    double density = 1.0);


/**
 * Return true only for a closed, consistently edge-oriented triangle
 * boundary. Every undirected edge must occur exactly twice, once in each
 * direction. Vertex-manifoldness and self-intersection are outside this
 * helper's contract.
 */
bool isClosedConsistentlyOrientedTriangleMesh(
    const Eigen::MatrixXi& faces
);


/**
 * Return true when the oriented triangle chain has zero boundary.
 *
 * Unlike isClosedConsistentlyOrientedTriangleMesh(), this predicate permits
 * more than two triangles to meet along an edge, provided that every use in
 * one direction is paired with a use in the opposite direction.  This is the
 * appropriate contract for voxel-cavity surfaces whose closed components may
 * touch along an edge: the surface remains a valid closed integration chain
 * even though its indexed connectivity is not a strict 2-manifold.
 */
bool isClosedOrientedTriangleChain(
    const Eigen::MatrixXi& faces
);


/** Compute the exact contribution of a uniform axis-aligned box. */
MassProperties computeAxisAlignedBoxMassProperties(
    const Eigen::Vector3d& center,
    const Eigen::Vector3d& size,
    double density = 1.0
);
