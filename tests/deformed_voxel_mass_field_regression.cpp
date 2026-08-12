#include "deformation/DeformedVoxelMassField.h"
#include "deformation/GridDeformationDomain.h"
#include "extraction/VoidMeshExtractor.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


using deformation::DeformedVoxelMassField;
using deformation::GridDeformationDomain;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] deformed-voxel-mass-field reason=" << message << '\n';
    return EXIT_FAILURE;
}


bool approximatelyEqual(
    double first,
    double second,
    double tolerance
)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <=
            tolerance * std::max({1.0, std::abs(first), std::abs(second)});
}


bool approximatelyEqual(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& second,
    double tolerance
)
{
    return
        first.allFinite() &&
        second.allFinite() &&
        (first - second).norm() <= tolerance *
            std::max({1.0, first.norm(), second.norm()});
}


void makeBox(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        minimum(0), minimum(1), minimum(2),
        maximum(0), minimum(1), minimum(2),
        maximum(0), maximum(1), minimum(2),
        minimum(0), maximum(1), minimum(2),
        minimum(0), minimum(1), maximum(2),
        maximum(0), minimum(1), maximum(2),
        maximum(0), maximum(1), maximum(2),
        minimum(0), maximum(1), maximum(2);

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


struct Fixture
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    VoxelGrid grid;
    GridDeformationDomain domain;

    Fixture()
        : outerVertices(), outerFaces(), grid(), domain()
    {
        makeBox(
            Eigen::Vector3d(-1.0, -0.8, -0.6),
            Eigen::Vector3d(1.0, 0.8, 0.6),
            outerVertices,
            outerFaces
        );
        grid.build(outerVertices, outerFaces, 0.0, 6);
        domain = GridDeformationDomain::fromVoxelGrid(grid);
    }
};


int firstActiveCell(const GridDeformationDomain& domain)
{
    for (int cell = 0; cell < domain.cellCount(); ++cell)
    {
        if (domain.isActiveCell(cell))
        {
            return cell;
        }
    }
    throw std::logic_error("fixture contains no active cell");
}


MassProperties sumField(
    const DeformedVoxelMassField& field,
    const std::vector<int>& cells
)
{
    MassProperties result;
    for (const int cell : cells)
    {
        result += field.massProperties(cell);
    }
    return result;
}


int runIdentityAndSimilarity()
{
    Fixture fixture;
    constexpr double density = 2.3;
    const DeformedVoxelMassField identity(
        fixture.domain,
        fixture.domain.fullRestCorners(),
        density
    );
    const int cell = firstActiveCell(fixture.domain);
    const Eigen::Vector3i coordinates = fixture.domain.cellCoordinates(cell);
    const Eigen::Vector3d expectedCenter =
        fixture.grid(coordinates(0), coordinates(1), coordinates(2)).center;
    const MassProperties expected = computeAxisAlignedBoxMassProperties(
        expectedCenter,
        fixture.domain.voxelSize(),
        density
    );
    const MassProperties actual = identity.massProperties(cell);

    if (
        identity.cellCount() != fixture.domain.cellCount() ||
        identity.activeCellCount() != fixture.domain.activeCellCount() ||
        identity.dimensions() != fixture.domain.dimensions() ||
        !identity.contains(cell) ||
        !approximatelyEqual(actual.mass, expected.mass, 1.0e-12) ||
        !approximatelyEqual(actual.firstMoment, expected.firstMoment, 1.0e-12) ||
        !approximatelyEqual(identity.centroid(cell), expectedCenter, 1.0e-12)
    )
    {
        return fail("identity field disagrees with the exact axis-aligned cell");
    }

    const Eigen::AngleAxisd rotation(
        0.63,
        Eigen::Vector3d(0.3, -0.7, 0.5).normalized()
    );
    const Eigen::Matrix3d R = rotation.toRotationMatrix();
    const Eigen::Vector3d translation(1.7, -0.9, 2.2);
    constexpr double scale = 1.35;

    const auto transformCorners = [&](
        const Eigen::Matrix3d& linear,
        const Eigen::Vector3d& offset
    )
    {
        Eigen::MatrixXd transformed = fixture.domain.fullRestCorners();
        for (Eigen::Index row = 0; row < transformed.rows(); ++row)
        {
            transformed.row(row) =
                (linear * transformed.row(row).transpose() + offset).
                    transpose();
        }
        return transformed;
    };

    const DeformedVoxelMassField translated(
        fixture.domain,
        transformCorners(Eigen::Matrix3d::Identity(), translation),
        density
    );
    const DeformedVoxelMassField rotated(
        fixture.domain,
        transformCorners(R, Eigen::Vector3d::Zero()),
        density
    );
    const DeformedVoxelMassField scaled(
        fixture.domain,
        transformCorners(
            scale * Eigen::Matrix3d::Identity(),
            Eigen::Vector3d::Zero()
        ),
        density
    );

    const DeformedVoxelMassField similarity(
        fixture.domain,
        transformCorners(scale * R, translation),
        density
    );
    const MassProperties transformedMass = similarity.massProperties(cell);
    const Eigen::Vector3d transformedCenter =
        scale * R * expectedCenter + translation;

    if (
        !approximatelyEqual(
            translated.massProperties(cell).mass,
            expected.mass,
            2.0e-12
        ) ||
        !approximatelyEqual(
            translated.centroid(cell),
            expectedCenter + translation,
            2.0e-12
        ) ||
        !approximatelyEqual(
            rotated.massProperties(cell).mass,
            expected.mass,
            2.0e-12
        ) ||
        !approximatelyEqual(
            rotated.centroid(cell),
            R * expectedCenter,
            2.0e-12
        ) ||
        !approximatelyEqual(
            scaled.massProperties(cell).mass,
            std::pow(scale, 3.0) * expected.mass,
            2.0e-12
        ) ||
        !approximatelyEqual(
            scaled.centroid(cell),
            scale * expectedCenter,
            2.0e-12
        ) ||
        !approximatelyEqual(
            transformedMass.mass,
            std::pow(scale, 3.0) * expected.mass,
            2.0e-12
        ) ||
        !approximatelyEqual(
            similarity.centroid(cell),
            transformedCenter,
            2.0e-12
        ) ||
        !approximatelyEqual(
            transformedMass.firstMoment,
            transformedMass.mass * transformedCenter,
            3.0e-12
        )
    )
    {
        return fail("translation/rotation/uniform-scale covariance failed");
    }

    std::cout
        << "[PASS] identity-similarity"
        << " mass=" << actual.mass
        << " scaled_mass=" << transformedMass.mass << '\n';
    return EXIT_SUCCESS;
}


int runDistortedCavitySum()
{
    Fixture fixture;
    Eigen::MatrixXd distorted = fixture.domain.fullRestCorners();

    for (Eigen::Index row = 0; row < distorted.rows(); ++row)
    {
        const Eigen::Vector3d p = distorted.row(row).transpose();
        Eigen::Vector3d q;
        q(0) = p(0) + 0.08 * p(1) +
            0.045 * std::sin(1.3 * p(1) + 0.2) *
                std::sin(0.9 * p(2) + 0.4);
        q(1) = p(1) + 0.05 * p(2) +
            0.035 * std::sin(1.1 * p(0) + 0.5) *
                std::sin(0.7 * p(2) - 0.2);
        q(2) = p(2) + 0.04 * p(0) +
            0.030 * std::sin(0.8 * p(0) - 0.3) *
                std::sin(1.2 * p(1) + 0.6);
        distorted.row(row) = q.transpose();
    }

    constexpr double density = 0.7;
    const DeformedVoxelMassField field(fixture.domain, distorted, density);
    std::vector<int> cells;

    for (int cell = 0; cell < fixture.domain.cellCount(); ++cell)
    {
        if (!fixture.domain.isActiveCell(cell))
        {
            continue;
        }

        const Eigen::Vector3i coordinate = fixture.domain.cellCoordinates(cell);
        if (fixture.grid(
                coordinate(0),
                coordinate(1),
                coordinate(2)
            ).isInner)
        {
            cells.push_back(cell);
            if (cells.size() == 3U)
            {
                break;
            }
        }
    }

    if (cells.size() != 3U)
    {
        return fail("fixture could not select three active deformed cells");
    }

    Eigen::Vector3d cornerAverage = Eigen::Vector3d::Zero();
    for (const int global : fixture.domain.cellGlobalCorners(cells.front()))
    {
        cornerAverage += distorted.row(global).transpose();
    }
    cornerAverage /= 8.0;
    const double centroidVersusCornerAverage =
        (field.centroid(cells.front()) - cornerAverage).norm();
    if (!(centroidVersusCornerAverage > 1.0e-8))
    {
        return fail("distorted fixture does not distinguish volume centroid from corner average");
    }

    const MassProperties summed = sumField(field, cells);
    VoxelGrid carved = fixture.grid;
    for (const int cell : cells)
    {
        const Eigen::Vector3i coordinate = fixture.domain.cellCoordinates(cell);
        carved(coordinate(0), coordinate(1), coordinate(2)).isFilled = false;
    }
    Eigen::MatrixXi cavityFaces;
    VoidMeshExtractor::extractFaces(carved, cavityFaces);
    const MassProperties inwardCavity = computeSurfaceMassProperties(
        distorted,
        cavityFaces,
        density
    );
    const MassProperties independent(
        -inwardCavity.mass,
        -inwardCavity.firstMoment
    );

    Eigen::Vector3d summedCenter = Eigen::Vector3d::Zero();
    Eigen::Vector3d independentCenter = Eigen::Vector3d::Zero();
    if (
        !summed.tryCenterOfMass(summedCenter) ||
        !independent.tryCenterOfMass(independentCenter) ||
        !approximatelyEqual(summed.mass, independent.mass, 1.0e-12) ||
        !approximatelyEqual(
            summed.firstMoment,
            independent.firstMoment,
            2.0e-12
        ) ||
        !approximatelyEqual(summedCenter, independentCenter, 2.0e-12)
    )
    {
        return fail("cell-sum mass disagrees with independent cavity surfaces");
    }

    std::cout
        << "[PASS] distorted-cavity-sum"
        << " mass=" << summed.mass
        << " center=" << summedCenter.transpose()
        << " centroid_delta=" << centroidVersusCornerAverage << '\n';
    return EXIT_SUCCESS;
}


int runInvalidInputs()
{
    Fixture fixture;
    const int cell = firstActiveCell(fixture.domain);
    const std::array<int, 8> globals = fixture.domain.cellGlobalCorners(cell);

    bool nonFiniteRejected = false;
    try
    {
        Eigen::MatrixXd invalid = fixture.domain.fullRestCorners();
        invalid(globals[0], 0) = std::numeric_limits<double>::quiet_NaN();
        static_cast<void>(DeformedVoxelMassField(fixture.domain, invalid));
    }
    catch (const std::invalid_argument&)
    {
        nonFiniteRejected = true;
    }

    bool invertedRejected = false;
    try
    {
        Eigen::MatrixXd reflected = fixture.domain.fullRestCorners();
        reflected.col(0) *= -1.0;
        static_cast<void>(DeformedVoxelMassField(fixture.domain, reflected));
    }
    catch (const deformation::InfeasibleDeformedVoxel&)
    {
        invertedRejected = true;
    }

    bool collapsedRejected = false;
    try
    {
        Eigen::MatrixXd collapsed = fixture.domain.fullRestCorners();
        collapsed.col(2).setZero();
        static_cast<void>(DeformedVoxelMassField(fixture.domain, collapsed));
    }
    catch (const deformation::InfeasibleDeformedVoxel&)
    {
        collapsedRejected = true;
    }

    bool rangeRejected = false;
    try
    {
        const DeformedVoxelMassField field(
            fixture.domain,
            fixture.domain.fullRestCorners()
        );
        static_cast<void>(field.massProperties(-1));
    }
    catch (const std::out_of_range&)
    {
        rangeRejected = true;
    }

    if (
        !nonFiniteRejected ||
        !invertedRejected ||
        !collapsedRejected ||
        !rangeRejected
    )
    {
        return fail("nonfinite/inverted/collapsed/out-of-range input was accepted");
    }

    std::cout << "[PASS] invalid-inputs\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: deformed-voxel-mass-field-regression "
            << "<identity-similarity|distorted-cavity-sum|invalid>\n";
        return EXIT_FAILURE;
    }

    const std::string scenario(argv[1]);
    try
    {
        if (scenario == "identity-similarity")
        {
            return runIdentityAndSimilarity();
        }
        if (scenario == "distorted-cavity-sum")
        {
            return runDistortedCavitySum();
        }
        if (scenario == "invalid")
        {
            return runInvalidInputs();
        }
        return fail("unknown scenario");
    }
    catch (const std::exception& error)
    {
        return fail(error.what());
    }
}
