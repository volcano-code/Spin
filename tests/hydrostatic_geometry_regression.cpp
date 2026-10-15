#include "spin/HydrostaticGeometry.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>

namespace
{

int fail(const char* message)
{
    std::cerr << "[FAIL] hydrostatic geometry " << message << '\n';
    return EXIT_FAILURE;
}


spin::TriangleMesh cube()
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices << -1, -1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,
                     -1, -1,  1,  1, -1,  1,  1,  1,  1, -1,  1,  1;
    mesh.faces.resize(12, 3);
    mesh.faces << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
                  0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5,
                  2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7;
    return mesh;
}


bool close(double actual, double expected, double tolerance = 1e-10)
{
    return std::abs(actual - expected) <= tolerance;
}


} // namespace


int main()
{
    const spin::TriangleMesh mesh = cube();
    const spin::HydrostaticGeometry horizontal = spin::computeHydrostaticGeometry(
        mesh, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ()
    );
    if (!close(horizontal.displacedVolume, 4.0) ||
        !horizontal.centerOfBuoyancy.isApprox(Eigen::Vector3d(0, 0, -0.5), 1e-10) ||
        !close(horizontal.waterlineArea, 4.0) ||
        !horizontal.waterlineCentroid.isZero(1e-10) ||
        !horizontal.waterlineCentralSecondMoment.isApprox(
            (Eigen::Matrix2d() << 4.0/3.0, 0, 0, 4.0/3.0).finished(), 1e-10) ||
        horizontal.submergedBoundary.vertices.rows() == 0 ||
        horizontal.submergedBoundary.faces.rows() == 0)
    {
        return fail("horizontal half-cube moments are incorrect");
    }

    const spin::HydrostaticGeometry translated = spin::computeHydrostaticGeometry(
        mesh, Eigen::Vector3d(0, 0, 0.5), Eigen::Vector3d::UnitZ()
    );
    if (!close(translated.displacedVolume, 6.0) ||
        !translated.centerOfBuoyancy.isApprox(Eigen::Vector3d(0, 0, -0.25), 1e-10) ||
        !close(translated.waterlineArea, 4.0) ||
        !translated.waterlineCentroid.isApprox(Eigen::Vector3d(0, 0, 0.5), 1e-10))
    {
        return fail("translated waterline moments are incorrect");
    }

    const spin::HydrostaticGeometry vertical = spin::computeHydrostaticGeometry(
        mesh, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX()
    );
    if (!close(vertical.displacedVolume, 4.0) ||
        !vertical.centerOfBuoyancy.isApprox(Eigen::Vector3d(-0.5, 0, 0), 1e-10) ||
        !close(vertical.waterlineArea, 4.0) ||
        !vertical.waterlineCentralSecondMoment.isApprox(
            (Eigen::Matrix2d() << 4.0/3.0, 0, 0, 4.0/3.0).finished(), 1e-10))
    {
        return fail("rotated waterline moments are incorrect");
    }

    bool rejected = false;

    // An equilibrium finder is allowed to solve only the scalar displaced-volume
    // equation for this explicitly supplied closed displacement boundary.  It
    // must not substitute a caller-filled volume or infer floodability.
    const spin::HydrostaticEquilibrium halfEquilibrium =
        spin::findHydrostaticEquilibrium(
            mesh, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(), 4.0
        );
    if (!close(halfEquilibrium.waterlineSignedOffset, 0.0, 1e-10) ||
        !close(halfEquilibrium.targetDisplacedVolume, 4.0) ||
        !close(halfEquilibrium.volumeResidual, 0.0, 1e-10) ||
        !close(halfEquilibrium.geometry.displacedVolume, 4.0) ||
        !halfEquilibrium.geometry.waterlineCentroid.isZero(1e-10))
    {
        return fail("half-volume equilibrium is not measured from the mesh");
    }

    const spin::HydrostaticEquilibrium translatedEquilibrium =
        spin::findHydrostaticEquilibrium(
            mesh, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(), 6.0
        );
    if (!close(translatedEquilibrium.waterlineSignedOffset, 0.5, 3e-10) ||
        !close(translatedEquilibrium.geometry.displacedVolume - 6.0,
               translatedEquilibrium.volumeResidual, 1e-13) ||
        !(std::abs(translatedEquilibrium.volumeResidual) <= 8e-10))
    {
        return fail("translated-volume equilibrium is not measured from the mesh");
    }

    rejected = false;
    try
    {
        (void)spin::findHydrostaticEquilibrium(
            mesh, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(), 8.0
        );
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    if (!rejected)
        return fail("full-volume equilibrium was accepted without a waterline");

    rejected = false;
    try
    {
        (void)spin::computeHydrostaticGeometry(
            mesh, Eigen::Vector3d(0, 0, 2), Eigen::Vector3d::UnitZ()
        );
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    if (!rejected)
        return fail("fully-submerged mesh was accepted without a waterline");

    std::cout << "[PASS] hydrostatic geometry actual mesh clipping and volume equilibrium\n";
    return EXIT_SUCCESS;
}
