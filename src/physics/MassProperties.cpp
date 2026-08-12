#include "physics/MassProperties.h"

#include <Eigen/Geometry>

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>


namespace
{

void requirePositiveFiniteDensity(
    double density
)
{
    if (!std::isfinite(density) || density <= 0.0)
    {
        throw std::invalid_argument(
            "mass density must be finite and positive"
        );
    }
}

} // namespace


MassProperties::MassProperties()
    :
    mass(0.0),
    firstMoment(Eigen::Vector3d::Zero())
{
}


MassProperties::MassProperties(
    double massValue,
    const Eigen::Vector3d& firstMomentValue
)
    :
    mass(massValue),
    firstMoment(firstMomentValue)
{
}


bool MassProperties::isFinite() const
{
    return
        std::isfinite(mass) &&
        firstMoment.allFinite();
}


bool MassProperties::tryCenterOfMass(
    Eigen::Vector3d& centerOfMass,
    double minimumAbsoluteMass
) const
{
    if (
        !isFinite() ||
        !std::isfinite(minimumAbsoluteMass) ||
        minimumAbsoluteMass < 0.0 ||
        std::abs(mass) <= minimumAbsoluteMass
    )
    {
        return false;
    }

    centerOfMass = firstMoment / mass;
    return centerOfMass.allFinite();
}


MassProperties& MassProperties::operator+=(
    const MassProperties& other
)
{
    mass += other.mass;
    firstMoment += other.firstMoment;
    return *this;
}


MassProperties& MassProperties::operator-=(
    const MassProperties& other
)
{
    mass -= other.mass;
    firstMoment -= other.firstMoment;
    return *this;
}


MassProperties operator+(
    MassProperties first,
    const MassProperties& second
)
{
    first += second;
    return first;
}


MassProperties operator-(
    MassProperties first,
    const MassProperties& second
)
{
    first -= second;
    return first;
}


MassProperties computeSurfaceMassProperties(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density
)
{
    requirePositiveFiniteDensity(density);

    if (vertices.cols() != 3 || faces.cols() != 3)
    {
        throw std::invalid_argument(
            "mass properties require 3D vertices and triangle faces"
        );
    }

    if (!vertices.allFinite())
    {
        throw std::invalid_argument(
            "mass-property vertices must be finite"
        );
    }

    MassProperties properties;

    if (vertices.rows() == 0 || faces.rows() == 0)
    {
        return properties;
    }

    /*
     * Use a reference point near the mesh to reduce cancellation after large
     * translations. Each oriented boundary triangle and the reference form
     * a signed tetrahedron. Summing those tetrahedra is equivalent to the
     * divergence-theorem surface integrals in the paper appendix.
     */
    const Eigen::Vector3d reference =
        vertices.colwise().mean().transpose();

    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        const int firstIndex = faces(row, 0);
        const int secondIndex = faces(row, 1);
        const int thirdIndex = faces(row, 2);

        if (
            firstIndex < 0 ||
            secondIndex < 0 ||
            thirdIndex < 0 ||
            firstIndex >= vertices.rows() ||
            secondIndex >= vertices.rows() ||
            thirdIndex >= vertices.rows()
        )
        {
            throw std::out_of_range(
                "mass-property face index is outside the vertex array"
            );
        }

        const Eigen::Vector3d first =
            vertices.row(firstIndex).transpose() - reference;

        const Eigen::Vector3d second =
            vertices.row(secondIndex).transpose() - reference;

        const Eigen::Vector3d third =
            vertices.row(thirdIndex).transpose() - reference;

        const double signedMass =
            density * first.dot(second.cross(third)) / 6.0;

        const Eigen::Vector3d tetrahedronCenter =
            reference + (first + second + third) / 4.0;

        properties.mass += signedMass;
        properties.firstMoment +=
            signedMass * tetrahedronCenter;
    }

    return properties;
}


bool isClosedConsistentlyOrientedTriangleMesh(
    const Eigen::MatrixXi& faces
)
{
    if (faces.cols() != 3 || faces.rows() == 0)
    {
        return false;
    }

    struct EdgeUse
    {
        int count = 0;
        int orientationSum = 0;
    };

    std::map<std::pair<int, int>, EdgeUse> edges;

    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        const int indices[3] =
        {
            faces(row, 0),
            faces(row, 1),
            faces(row, 2)
        };

        if (
            indices[0] < 0 ||
            indices[1] < 0 ||
            indices[2] < 0 ||
            indices[0] == indices[1] ||
            indices[1] == indices[2] ||
            indices[2] == indices[0]
        )
        {
            return false;
        }

        for (int edgeIndex = 0; edgeIndex < 3; ++edgeIndex)
        {
            const int from = indices[edgeIndex];
            const int to = indices[(edgeIndex + 1) % 3];

            const std::pair<int, int> key =
                std::minmax(from, to);

            EdgeUse& use = edges[key];
            ++use.count;
            use.orientationSum += from < to ? 1 : -1;
        }
    }

    for (
        std::map<std::pair<int, int>, EdgeUse>::const_iterator iterator =
            edges.begin();
        iterator != edges.end();
        ++iterator
    )
    {
        if (
            iterator->second.count != 2 ||
            iterator->second.orientationSum != 0
        )
        {
            return false;
        }
    }

    return true;
}


MassProperties computeAxisAlignedBoxMassProperties(
    const Eigen::Vector3d& center,
    const Eigen::Vector3d& size,
    double density
)
{
    requirePositiveFiniteDensity(density);

    if (
        !center.allFinite() ||
        !size.allFinite() ||
        (size.array() <= 0.0).any()
    )
    {
        throw std::invalid_argument(
            "axis-aligned box center/size must be finite and size positive"
        );
    }

    const double mass =
        density * size(0) * size(1) * size(2);

    return MassProperties(
        mass,
        mass * center
    );
}
