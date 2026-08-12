#include "physics/SurfaceMassDifferential.h"

#include <Eigen/Geometry>

#include <cmath>
#include <stdexcept>


namespace
{


void validateSurfaceInputs(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density
)
{
    if (
        vertices.cols() != 3 ||
        !vertices.allFinite() ||
        faces.cols() != 3 ||
        !std::isfinite(density) ||
        density <= 0.0
    )
    {
        throw std::invalid_argument(
            "surface differential requires finite 3D geometry and positive density"
        );
    }

    if (
        faces.size() > 0 &&
        (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows())
    )
    {
        throw std::invalid_argument(
            "surface differential face index is outside the vertex array"
        );
    }

    if (
        faces.rows() > 0 &&
        !isClosedOrientedTriangleChain(faces)
    )
    {
        throw std::invalid_argument(
            "surface differential requires a closed oriented triangle chain"
        );
    }
}


} // namespace


SurfaceMassDifferential computeSurfaceMassDifferential(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density
)
{
    validateSurfaceInputs(vertices, faces, density);

    SurfaceMassDifferential result;
    result.massGradient = Eigen::VectorXd::Zero(3 * vertices.rows());
    result.firstMomentJacobian =
        Eigen::MatrixXd::Zero(3, 3 * vertices.rows());

    if (vertices.rows() == 0 || faces.rows() == 0)
    {
        return result;
    }

    const Eigen::Vector3d reference =
        vertices.colwise().mean().transpose();

    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        const int indices[3] =
        {
            faces(face, 0),
            faces(face, 1),
            faces(face, 2)
        };

        const Eigen::Vector3d a =
            vertices.row(indices[0]).transpose() - reference;
        const Eigen::Vector3d b =
            vertices.row(indices[1]).transpose() - reference;
        const Eigen::Vector3d c =
            vertices.row(indices[2]).transpose() - reference;

        const double signedMass =
            density * a.dot(b.cross(c)) / 6.0;
        const Eigen::Vector3d tetrahedronCenter =
            reference + (a + b + c) / 4.0;

        result.properties.mass += signedMass;
        result.properties.firstMoment +=
            signedMass * tetrahedronCenter;

        const Eigen::Vector3d massGradients[3] =
        {
            density * b.cross(c) / 6.0,
            density * c.cross(a) / 6.0,
            density * a.cross(b) / 6.0
        };

        for (int local = 0; local < 3; ++local)
        {
            const Eigen::Index column = 3 * indices[local];
            result.massGradient.segment<3>(column) +=
                massGradients[local];
            result.firstMomentJacobian.block<3, 3>(0, column) +=
                tetrahedronCenter * massGradients[local].transpose() +
                (signedMass / 4.0) * Eigen::Matrix3d::Identity();
        }
    }

    if (
        !result.properties.isFinite() ||
        !result.massGradient.allFinite() ||
        !result.firstMomentJacobian.allFinite()
    )
    {
        throw std::overflow_error(
            "surface mass differential overflowed"
        );
    }

    return result;
}


CenterOfMassDifferential computeCenterOfMassDifferential(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double density,
    double minimumPositiveMass
)
{
    if (
        !std::isfinite(minimumPositiveMass) ||
        minimumPositiveMass < 0.0
    )
    {
        throw std::invalid_argument(
            "minimum positive mass must be finite and non-negative"
        );
    }

    const SurfaceMassDifferential outer =
        computeSurfaceMassDifferential(
            outerVertices,
            outerFaces,
            density
        );
    const SurfaceMassDifferential inner =
        computeSurfaceMassDifferential(
            innerVertices,
            innerFaces,
            density
        );

    CenterOfMassDifferential result;
    result.properties = outer.properties + inner.properties;

    if (!(outer.properties.mass > 0.0))
    {
        throw std::domain_error(
            "outer material boundary must have outward positive orientation"
        );
    }

    if (innerFaces.rows() > 0 && !(inner.properties.mass < 0.0))
    {
        throw std::domain_error(
            "inner cavity boundary must have inward negative orientation"
        );
    }

    if (
        !std::isfinite(result.properties.mass) ||
        result.properties.mass <= minimumPositiveMass ||
        !result.properties.firstMoment.allFinite()
    )
    {
        throw std::domain_error(
            "deformed material must have finite positive mass"
        );
    }

    result.centerOfMass =
        result.properties.firstMoment / result.properties.mass;

    result.outerJacobian =
        (
            outer.firstMomentJacobian -
            result.centerOfMass * outer.massGradient.transpose()
        ) / result.properties.mass;
    result.innerJacobian =
        (
            inner.firstMomentJacobian -
            result.centerOfMass * inner.massGradient.transpose()
        ) / result.properties.mass;

    if (
        !result.centerOfMass.allFinite() ||
        !result.outerJacobian.allFinite() ||
        !result.innerJacobian.allFinite()
    )
    {
        throw std::overflow_error(
            "center-of-mass differential overflowed"
        );
    }

    return result;
}
