#include "center_of_mass.h"

#include <cmath>
#include <stdexcept>


double center_of_mass(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    Eigen::Vector3d& center
)
{
    const MassProperties properties =
        computeSurfaceMassProperties(
            vertices,
            faces
        );

    if (!properties.tryCenterOfMass(center))
    {
        throw std::domain_error(
            "center of mass is undefined for zero/non-finite signed mass"
        );
    }

    return properties.mass;
}


void center_of_mass(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double mass,
    Eigen::Vector3d& center
)
{
    if (!std::isfinite(mass) || mass == 0.0)
    {
        throw std::domain_error(
            "center of mass is undefined for zero/non-finite supplied mass"
        );
    }

    const MassProperties properties =
        computeSurfaceMassProperties(
            vertices,
            faces
        );

    center = properties.firstMoment / mass;

    if (!center.allFinite())
    {
        throw std::domain_error(
            "center of mass computation produced a non-finite result"
        );
    }
}
