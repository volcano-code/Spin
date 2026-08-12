#pragma once

#include "physics/MassProperties.h"

#include <Eigen/Core>


/** Analytic differential of signed mass and first moment. */
struct SurfaceMassDifferential
{
    MassProperties properties;
    Eigen::VectorXd massGradient;
    Eigen::MatrixXd firstMomentJacobian;
};


/**
 * Differentiate a closed oriented triangle boundary.
 *
 * Coordinates use the interleaved ordering [x0,y0,z0,x1,...]. An empty face
 * set is accepted as a zero surface, which represents an object with no
 * carved cavity.
 */
SurfaceMassDifferential computeSurfaceMassDifferential(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density = 1.0
);


/** Differential of the material between outward outer and inward inner. */
struct CenterOfMassDifferential
{
    MassProperties properties;
    Eigen::Vector3d centerOfMass;
    Eigen::MatrixXd outerJacobian;
    Eigen::MatrixXd innerJacobian;
};


CenterOfMassDifferential computeCenterOfMassDifferential(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double density = 1.0,
    double minimumPositiveMass = 0.0
);
