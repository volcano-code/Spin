#pragma once

#include <Eigen/Core>

#include "physics/MassProperties.h"


/**
 * Compatibility wrapper around computeSurfaceMassProperties().
 *
 * The input must be a consistently oriented closed triangle boundary. The
 * returned value is signed mass and center receives firstMoment / mass.
 */
double center_of_mass(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    Eigen::Vector3d& center
);


/**
 * Divide a closed, consistently oriented boundary's signed first-moment
 * contribution by a supplied mass. Retained for source compatibility; new
 * code should use MassProperties so mass and first moment cannot be updated
 * independently.
 */
void center_of_mass(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double mass,
    Eigen::Vector3d& center
);
