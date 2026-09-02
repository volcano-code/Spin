#include "manufacturing/ManufacturingGeometry.h"

#include <igl/MeshBooleanType.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


struct Mesh
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};


Mesh cube(const Eigen::Vector3d& minimum, const Eigen::Vector3d& maximum)
{
    Mesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices <<
        minimum.x(), minimum.y(), minimum.z(),
        maximum.x(), minimum.y(), minimum.z(),
        maximum.x(), maximum.y(), minimum.z(),
        minimum.x(), maximum.y(), minimum.z(),
        minimum.x(), minimum.y(), maximum.z(),
        maximum.x(), minimum.y(), maximum.z(),
        maximum.x(), maximum.y(), maximum.z(),
        minimum.x(), maximum.y(), maximum.z();
    mesh.faces.resize(12, 3);
    mesh.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
    return mesh;
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-boolean " << message << '\n';
    return EXIT_FAILURE;
}


double signedVolume(const Mesh& mesh)
{
    double sixVolume = 0.0;
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        const Eigen::Vector3d a = mesh.vertices.row(mesh.faces(face, 0));
        const Eigen::Vector3d b = mesh.vertices.row(mesh.faces(face, 1));
        const Eigen::Vector3d c = mesh.vertices.row(mesh.faces(face, 2));
        sixVolume += a.dot(b.cross(c));
    }
    return sixVolume / 6.0;
}


bool intersection(const Mesh& first, const Mesh& second, Mesh* result)
{
    Eigen::VectorXi birthFaces;
    return igl::copyleft::cgal::mesh_boolean(
        first.vertices,
        first.faces,
        second.vertices,
        second.faces,
        igl::MESH_BOOLEAN_TYPE_INTERSECT,
        result->vertices,
        result->faces,
        birthFaces
    );
}


bool validSolid(const Mesh& mesh, double expectedVolume)
{
    if (mesh.vertices.rows() == 0 || mesh.faces.rows() == 0)
    {
        return false;
    }
    const manufacturing::BoundaryAudit audit =
        manufacturing::auditTriangleBoundary(
            mesh.vertices,
            mesh.faces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    return audit.strictBoundaryValid() &&
        std::abs(signedVolume(mesh) - expectedVolume) < 1e-12;
}


} // namespace


int main()
{
    const Mesh unit = cube(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );

    Mesh result;
    const Mesh disjoint = cube(
        Eigen::Vector3d(2.0, 0.0, 0.0),
        Eigen::Vector3d(3.0, 1.0, 1.0)
    );
    if (!intersection(unit, disjoint, &result) || result.faces.rows() != 0)
    {
        return fail("disjoint intersection did not produce an empty solid");
    }

    const Mesh contained = cube(
        Eigen::Vector3d::Constant(0.2),
        Eigen::Vector3d::Constant(0.8)
    );
    if (!intersection(unit, contained, &result) ||
        !validSolid(result, 0.216))
    {
        return fail("contained intersection is not the inner solid");
    }

    const Mesh tangent = cube(
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(2.0, 1.0, 1.0)
    );
    if (!intersection(unit, tangent, &result) || result.faces.rows() != 0)
    {
        return fail("face-tangent intersection was not classified as empty");
    }

    const Mesh overlapping = cube(
        Eigen::Vector3d(0.5, 0.0, 0.0),
        Eigen::Vector3d(1.5, 1.0, 1.0)
    );
    if (!intersection(unit, overlapping, &result) ||
        !validSolid(result, 0.5))
    {
        return fail("overlapping intersection has invalid geometry/volume");
    }

    std::cout << "[PASS] scenario=spin-boolean cgal=f7c3c8212b56c0d6dae63787efc99093f4383415\n";
    return EXIT_SUCCESS;
}
