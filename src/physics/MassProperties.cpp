#include "physics/MassProperties.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>


namespace
{


struct EdgeUse
{
    int forwardCount = 0;
    int reverseCount = 0;
};


using EdgeUseMap = std::map<std::pair<int, int>, EdgeUse>;


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


bool collectOrientedEdgeUses(
    const Eigen::MatrixXi& faces,
    EdgeUseMap& edges
)
{
    if (faces.cols() != 3 || faces.rows() == 0)
    {
        return false;
    }

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
            EdgeUse& use = edges[std::minmax(from, to)];

            if (from < to)
            {
                ++use.forwardCount;
            }
            else
            {
                ++use.reverseCount;
            }
        }
    }

    return true;
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


MassMoments::MassMoments()
    :
    mass(0.0),
    firstMoment(Eigen::Vector3d::Zero()),
    secondMoment(Eigen::Matrix3d::Zero())
{
}


bool MassMoments::isFinite() const
{
    return
        std::isfinite(mass) &&
        firstMoment.allFinite() &&
        secondMoment.allFinite();
}


bool MassMoments::tryCenterOfMass(
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


Eigen::Matrix3d MassMoments::inertiaAboutCenter(
    double minimumPositiveMass
) const
{
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    if (
        !std::isfinite(minimumPositiveMass) ||
        minimumPositiveMass < 0.0 ||
        !(mass > minimumPositiveMass) ||
        !tryCenterOfMass(center, minimumPositiveMass)
    )
    {
        throw std::domain_error(
            "center inertia requires finite positive mass moments"
        );
    }

    const Eigen::Matrix3d inertiaAtOrigin =
        secondMoment.trace() * Eigen::Matrix3d::Identity() - secondMoment;
    const Eigen::Matrix3d parallelAxis = mass * (
        center.squaredNorm() * Eigen::Matrix3d::Identity() -
        center * center.transpose()
    );
    return inertiaAtOrigin - parallelAxis;
}


MassMoments& MassMoments::operator+=(const MassMoments& other)
{
    mass += other.mass;
    firstMoment += other.firstMoment;
    secondMoment += other.secondMoment;
    return *this;
}


MassMoments& MassMoments::operator-=(const MassMoments& other)
{
    mass -= other.mass;
    firstMoment -= other.firstMoment;
    secondMoment -= other.secondMoment;
    return *this;
}


MassMoments operator+(MassMoments first, const MassMoments& second)
{
    first += second;
    return first;
}


MassMoments operator-(MassMoments first, const MassMoments& second)
{
    first -= second;
    return first;
}


namespace {
template<class Scalar>
struct SurfaceMomentAccumulator {
    Scalar mass = 0;
    Eigen::Matrix<Scalar,3,1> firstMoment = Eigen::Matrix<Scalar,3,1>::Zero();
    Eigen::Matrix<Scalar,3,3> secondMoment = Eigen::Matrix<Scalar,3,3>::Zero();
    MassMoments binary64() const {
        MassMoments out;
        out.mass = static_cast<double>(mass);
        out.firstMoment = firstMoment.template cast<double>();
        out.secondMoment = secondMoment.template cast<double>();
        return out;
    }
};

template<class Scalar>
SurfaceMomentAccumulator<Scalar> accumulateSurfaceMassMoments(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density
)
{
    requirePositiveFiniteDensity(density);
    if (vertices.cols() != 3 || faces.cols() != 3)
    {
        throw std::invalid_argument(
            "mass moments require 3D vertices and triangle faces"
        );
    }
    if (!vertices.allFinite())
    {
        throw std::invalid_argument("mass-moment vertices must be finite");
    }

    SurfaceMomentAccumulator<Scalar> moments;
    if (vertices.rows() == 0 || faces.rows() == 0)
    {
        return moments;
    }

    using Vector3 = Eigen::Matrix<Scalar,3,1>;
    using Matrix3 = Eigen::Matrix<Scalar,3,3>;
    const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> points = vertices.cast<Scalar>();
    const Vector3 reference = points.colwise().mean().transpose();
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        const int indices[3] =
        {
            faces(row, 0), faces(row, 1), faces(row, 2)
        };
        for (const int index : indices)
        {
            if (index < 0 || index >= vertices.rows())
            {
                throw std::out_of_range(
                    "mass-moment face index is outside the vertex array"
                );
            }
        }

        const Vector3 verticesWorld[4] =
        {
            reference,
            points.row(indices[0]).transpose(),
            points.row(indices[1]).transpose(),
            points.row(indices[2]).transpose()
        };
        const Vector3 first = verticesWorld[1] - reference;
        const Vector3 second = verticesWorld[2] - reference;
        const Vector3 third = verticesWorld[3] - reference;
        const Scalar signedMass =
            Scalar(density) * first.dot(second.cross(third)) / Scalar(6);

        Vector3 sum = Vector3::Zero();
        Matrix3 selfOuter = Matrix3::Zero();
        for (const Vector3& vertex : verticesWorld)
        {
            sum += vertex;
            selfOuter += vertex * vertex.transpose();
        }

        moments.mass += signedMass;
        moments.firstMoment += signedMass * sum / Scalar(4);
        moments.secondMoment +=
            signedMass * (sum * sum.transpose() + selfOuter) / Scalar(20);
    }
    return moments;
}
} // namespace

MassMoments computeSurfaceMassMoments(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, double density) {
    return accumulateSurfaceMassMoments<double>(vertices,faces,density).binary64();
}

MassMoments computeSurfaceDifferenceMassMomentsExtended(
    const Eigen::MatrixXd& outerVertices, const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices, const Eigen::MatrixXi& innerFaces,
    double density) {
    if (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
        throw std::runtime_error("extended input moment accumulation is unavailable on this platform");
    auto out = accumulateSurfaceMassMoments<long double>(outerVertices,outerFaces,density);
    const auto inner = accumulateSurfaceMassMoments<long double>(innerVertices,innerFaces,density);
    out.mass -= inner.mass;
    out.firstMoment -= inner.firstMoment;
    out.secondMoment -= inner.secondMoment;
    return out.binary64();
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
    EdgeUseMap edges;
    if (!collectOrientedEdgeUses(faces, edges))
    {
        return false;
    }

    for (const auto& entry : edges)
    {
        if (
            entry.second.forwardCount != 1 ||
            entry.second.reverseCount != 1
        )
        {
            return false;
        }
    }

    return true;
}


bool isClosedOrientedTriangleChain(
    const Eigen::MatrixXi& faces
)
{
    EdgeUseMap edges;
    if (!collectOrientedEdgeUses(faces, edges))
    {
        return false;
    }

    for (const auto& entry : edges)
    {
        if (
            entry.second.forwardCount != entry.second.reverseCount
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
