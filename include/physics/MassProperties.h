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


/**
 * Return true only for a closed, consistently edge-oriented triangle
 * boundary. Every undirected edge must occur exactly twice, once in each
 * direction. Vertex-manifoldness and self-intersection are outside this
 * helper's contract.
 */
bool isClosedConsistentlyOrientedTriangleMesh(
    const Eigen::MatrixXi& faces
);


/** Compute the exact contribution of a uniform axis-aligned box. */
MassProperties computeAxisAlignedBoxMassProperties(
    const Eigen::Vector3d& center,
    const Eigen::Vector3d& size,
    double density = 1.0
);
