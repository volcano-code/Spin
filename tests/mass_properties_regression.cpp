#include <cstdint>

#include "physics/BalanceEnergy.h"
#include "physics/MassProperties.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>


namespace
{

constexpr double ABSOLUTE_TOLERANCE = 1e-12;
constexpr double RELATIVE_TOLERANCE = 1e-12;


bool approximatelyEqual(
    double actual,
    double expected
)
{
    if (!std::isfinite(actual) || !std::isfinite(expected))
    {
        return false;
    }

    const double scale =
        std::max(
            std::abs(actual),
            std::abs(expected)
        );

    return
        std::abs(actual - expected) <=
        ABSOLUTE_TOLERANCE +
            RELATIVE_TOLERANCE * scale;
}


bool approximatelyEqual(
    const Eigen::Vector3d& actual,
    const Eigen::Vector3d& expected
)
{
    return
        approximatelyEqual(actual(0), expected(0)) &&
        approximatelyEqual(actual(1), expected(1)) &&
        approximatelyEqual(actual(2), expected(2));
}


int fail(
    const std::string& scenario,
    const std::string& message
)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " " << message
        << '\n';

    return EXIT_FAILURE;
}


bool checkProperties(
    const MassProperties& properties,
    double expectedMass,
    const Eigen::Vector3d& expectedFirstMoment,
    const Eigen::Vector3d& expectedCenterOfMass,
    std::string& error
)
{
    if (!approximatelyEqual(properties.mass, expectedMass))
    {
        error = "mass disagrees with the analytic value";
        return false;
    }

    if (
        !approximatelyEqual(
            properties.firstMoment,
            expectedFirstMoment
        )
    )
    {
        error = "first moment disagrees with the analytic value";
        return false;
    }

    Eigen::Vector3d centerOfMass =
        Eigen::Vector3d::Zero();

    if (!properties.tryCenterOfMass(centerOfMass))
    {
        error = "positive-mass fixture did not expose a center of mass";
        return false;
    }

    if (!approximatelyEqual(centerOfMass, expectedCenterOfMass))
    {
        error = "center of mass disagrees with the analytic value";
        return false;
    }

    if (
        !approximatelyEqual(
            properties.mass * centerOfMass,
            properties.firstMoment
        )
    )
    {
        error = "mass, first moment and center of mass are inconsistent";
        return false;
    }

    return true;
}


void makeAxisAlignedBoxMesh(
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

    // Every triangle is outward-facing.
    faces.resize(12, 3);

    faces <<
        0, 2, 1,
        0, 3, 2,
        4, 5, 6,
        4, 6, 7,
        0, 1, 5,
        0, 5, 4,
        3, 7, 6,
        3, 6, 2,
        0, 4, 7,
        0, 7, 3,
        1, 2, 6,
        1, 6, 5;
}


Eigen::MatrixXi reversedFaces(
    const Eigen::MatrixXi& faces
)
{
    Eigen::MatrixXi reversed = faces;

    for (Eigen::Index row = 0; row < reversed.rows(); ++row)
    {
        std::swap(
            reversed(row, 1),
            reversed(row, 2)
        );
    }

    return reversed;
}


void concatenateMeshes(
    const Eigen::MatrixXd& firstVertices,
    const Eigen::MatrixXi& firstFaces,
    const Eigen::MatrixXd& secondVertices,
    const Eigen::MatrixXi& secondFaces,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(
        firstVertices.rows() + secondVertices.rows(),
        3
    );

    vertices.topRows(firstVertices.rows()) =
        firstVertices;

    vertices.bottomRows(secondVertices.rows()) =
        secondVertices;

    faces.resize(
        firstFaces.rows() + secondFaces.rows(),
        3
    );

    faces.topRows(firstFaces.rows()) =
        firstFaces;

    faces.bottomRows(secondFaces.rows()) =
        (
            secondFaces.array() +
            firstVertices.rows()
        ).matrix();
}


Eigen::MatrixXd translatedVertices(
    const Eigen::MatrixXd& vertices,
    const Eigen::Vector3d& translation
)
{
    Eigen::MatrixXd transformed = vertices;

    transformed.rowwise() += translation.transpose();

    return transformed;
}


Eigen::MatrixXd rotatedVertices(
    const Eigen::MatrixXd& vertices,
    const Eigen::Matrix3d& rotation
)
{
    Eigen::MatrixXd transformed(vertices.rows(), 3);

    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        transformed.row(row) =
            (
                rotation *
                vertices.row(row).transpose()
            ).transpose();
    }

    return transformed;
}


int runSurfaceScenario()
{
    const std::string scenario = "surface";
    std::string error;

    Eigen::MatrixXd cuboidVertices;
    Eigen::MatrixXi cuboidFaces;

    makeAxisAlignedBoxMesh(
        Eigen::Vector3d(0.0, -1.0, 4.0),
        Eigen::Vector3d(2.0, 2.0, 8.0),
        cuboidVertices,
        cuboidFaces
    );

    const MassProperties cuboid =
        computeSurfaceMassProperties(
            cuboidVertices,
            cuboidFaces,
            1.0
        );

    if (
        !checkProperties(
            cuboid,
            24.0,
            Eigen::Vector3d(24.0, 12.0, 144.0),
            Eigen::Vector3d(1.0, 0.5, 6.0),
            error
        )
    )
    {
        return fail(scenario, "cuboid: " + error);
    }

    const MassProperties reversedCuboid =
        computeSurfaceMassProperties(
            cuboidVertices,
            reversedFaces(cuboidFaces),
            1.0
        );

    Eigen::Vector3d reversedCenterOfMass =
        Eigen::Vector3d::Zero();

    if (
        !approximatelyEqual(reversedCuboid.mass, -cuboid.mass) ||
        !approximatelyEqual(
            reversedCuboid.firstMoment,
            -cuboid.firstMoment
        ) ||
        !reversedCuboid.tryCenterOfMass(reversedCenterOfMass) ||
        !approximatelyEqual(
            reversedCenterOfMass,
            Eigen::Vector3d(1.0, 0.5, 6.0)
        )
    )
    {
        return fail(
            scenario,
            "global winding reversal did not preserve CoM and negate signed properties"
        );
    }

    const MassProperties denseCuboid =
        computeSurfaceMassProperties(
            cuboidVertices,
            cuboidFaces,
            2.5
        );

    if (
        !checkProperties(
            denseCuboid,
            60.0,
            Eigen::Vector3d(60.0, 30.0, 360.0),
            Eigen::Vector3d(1.0, 0.5, 6.0),
            error
        )
    )
    {
        return fail(scenario, "density scaling: " + error);
    }

    Eigen::MatrixXd tetrahedronVertices(4, 3);
    Eigen::MatrixXi tetrahedronFaces(4, 3);

    tetrahedronVertices <<
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 4.0;

    tetrahedronFaces <<
        1, 2, 3,
        0, 3, 2,
        0, 1, 3,
        0, 2, 1;

    const MassProperties tetrahedron =
        computeSurfaceMassProperties(
            tetrahedronVertices,
            tetrahedronFaces,
            1.0
        );

    if (
        !checkProperties(
            tetrahedron,
            4.0,
            Eigen::Vector3d(2.0, 3.0, 4.0),
            Eigen::Vector3d(0.5, 0.75, 1.0),
            error
        )
    )
    {
        return fail(scenario, "tetrahedron: " + error);
    }

    const Eigen::Vector3d translation(3.0, -4.0, 5.0);

    const MassProperties translated =
        computeSurfaceMassProperties(
            translatedVertices(
                cuboidVertices,
                translation
            ),
            cuboidFaces,
            1.0
        );

    if (
        !checkProperties(
            translated,
            cuboid.mass,
            cuboid.firstMoment +
                cuboid.mass * translation,
            Eigen::Vector3d(1.0, 0.5, 6.0) +
                translation,
            error
        )
    )
    {
        return fail(scenario, "translation: " + error);
    }

    Eigen::Matrix3d rotation;

    rotation <<
        0.0, -1.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0;

    const MassProperties rotated =
        computeSurfaceMassProperties(
            rotatedVertices(
                cuboidVertices,
                rotation
            ),
            cuboidFaces,
            1.0
        );

    if (
        !checkProperties(
            rotated,
            cuboid.mass,
            rotation * cuboid.firstMoment,
            rotation * Eigen::Vector3d(1.0, 0.5, 6.0),
            error
        )
    )
    {
        return fail(scenario, "rotation: " + error);
    }

    constexpr double SCALE = 2.0;

    const MassProperties scaled =
        computeSurfaceMassProperties(
            SCALE * cuboidVertices,
            cuboidFaces,
            1.0
        );

    if (
        !checkProperties(
            scaled,
            std::pow(SCALE, 3.0) * cuboid.mass,
            std::pow(SCALE, 4.0) * cuboid.firstMoment,
            SCALE * Eigen::Vector3d(1.0, 0.5, 6.0),
            error
        )
    )
    {
        return fail(scenario, "uniform scale: " + error);
    }

    constexpr double TINY_SCALE = 1e-6;

    const MassProperties tiny =
        computeSurfaceMassProperties(
            TINY_SCALE * cuboidVertices,
            cuboidFaces,
            1.0
        );

    Eigen::Vector3d tinyCenter = Eigen::Vector3d::Zero();

    if (
        !tiny.isFinite() ||
        tiny.mass <= 0.0 ||
        !tiny.tryCenterOfMass(tinyCenter) ||
        !approximatelyEqual(
            tinyCenter,
            TINY_SCALE * Eigen::Vector3d(1.0, 0.5, 6.0)
        )
    )
    {
        return fail(
            scenario,
            "small-scale positive solid was treated as zero mass"
        );
    }

    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    Eigen::MatrixXd innerVertices;
    Eigen::MatrixXi innerFaces;

    makeAxisAlignedBoxMesh(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Constant(4.0),
        outerVertices,
        outerFaces
    );

    makeAxisAlignedBoxMesh(
        Eigen::Vector3d::Constant(1.0),
        Eigen::Vector3d::Constant(2.0),
        innerVertices,
        innerFaces
    );

    const Eigen::MatrixXi inwardInnerFaces =
        reversedFaces(innerFaces);

    const MassProperties innerContribution =
        computeSurfaceMassProperties(
            innerVertices,
            inwardInnerFaces,
            1.0
        );

    if (
        !approximatelyEqual(innerContribution.mass, -1.0) ||
        !approximatelyEqual(
            innerContribution.firstMoment,
            Eigen::Vector3d::Constant(-1.5)
        )
    )
    {
        return fail(
            scenario,
            "inward surface did not produce a signed negative contribution"
        );
    }

    Eigen::MatrixXd hollowVertices;
    Eigen::MatrixXi hollowFaces;

    concatenateMeshes(
        outerVertices,
        outerFaces,
        innerVertices,
        inwardInnerFaces,
        hollowVertices,
        hollowFaces
    );

    const MassProperties hollow =
        computeSurfaceMassProperties(
            hollowVertices,
            hollowFaces,
            1.0
        );

    const double hollowCenter = 253.0 / 126.0;

    if (
        !checkProperties(
            hollow,
            63.0,
            Eigen::Vector3d::Constant(126.5),
            Eigen::Vector3d::Constant(hollowCenter),
            error
        )
    )
    {
        return fail(scenario, "inward hollow: " + error);
    }

    Eigen::MatrixXd degenerateVertices(3, 3);
    Eigen::MatrixXi degenerateFaces(1, 3);

    degenerateVertices <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0;

    degenerateFaces << 0, 1, 2;

    const MassProperties degenerate =
        computeSurfaceMassProperties(
            degenerateVertices,
            degenerateFaces
        );

    Eigen::Vector3d undefinedCenter =
        Eigen::Vector3d::Zero();

    if (
        !approximatelyEqual(degenerate.mass, 0.0) ||
        !approximatelyEqual(
            degenerate.firstMoment,
            Eigen::Vector3d::Zero()
        ) ||
        degenerate.tryCenterOfMass(undefinedCenter)
    )
    {
        return fail(
            scenario,
            "zero-volume boundary did not preserve the undefined-CoM contract"
        );
    }

    Eigen::MatrixXi openCuboidFaces =
        cuboidFaces.topRows(cuboidFaces.rows() - 1);

    Eigen::MatrixXi mixedWindingFaces = cuboidFaces;
    std::swap(
        mixedWindingFaces(0, 1),
        mixedWindingFaces(0, 2)
    );

    Eigen::MatrixXi nonManifoldFaces(
        cuboidFaces.rows() + 1,
        3
    );

    nonManifoldFaces.topRows(cuboidFaces.rows()) =
        cuboidFaces;

    nonManifoldFaces.bottomRows(1) =
        cuboidFaces.topRows(1);

    if (
        !isClosedConsistentlyOrientedTriangleMesh(cuboidFaces) ||
        !isClosedConsistentlyOrientedTriangleMesh(
            reversedFaces(cuboidFaces)
        ) ||
        isClosedConsistentlyOrientedTriangleMesh(openCuboidFaces) ||
        isClosedConsistentlyOrientedTriangleMesh(mixedWindingFaces) ||
        isClosedConsistentlyOrientedTriangleMesh(nonManifoldFaces)
    )
    {
        return fail(
            scenario,
            "closed/edge-orientation boundary validation is incorrect"
        );
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


MassProperties sumUnitVoxelContributions(
    const std::array<Eigen::Vector3d, 4>& centers
)
{
    MassProperties result =
        computeAxisAlignedBoxMassProperties(
            centers.front(),
            Eigen::Vector3d::Ones(),
            1.0
        );

    for (std::size_t index = 1; index < centers.size(); ++index)
    {
        const MassProperties contribution =
            computeAxisAlignedBoxMassProperties(
                centers[index],
                Eigen::Vector3d::Ones(),
                1.0
            );

        result.mass += contribution.mass;
        result.firstMoment += contribution.firstMoment;
    }

    return result;
}


int runVoxelScenario()
{
    const std::string scenario = "voxel";
    std::string error;

    const MassProperties anisotropic =
        computeAxisAlignedBoxMassProperties(
            Eigen::Vector3d(1.0, -2.0, 3.0),
            Eigen::Vector3d(2.0, 3.0, 4.0),
            1.0
        );

    if (
        !checkProperties(
            anisotropic,
            24.0,
            Eigen::Vector3d(24.0, -48.0, 72.0),
            Eigen::Vector3d(1.0, -2.0, 3.0),
            error
        )
    )
    {
        return fail(scenario, "anisotropic box: " + error);
    }

    const MassProperties denseAnisotropic =
        computeAxisAlignedBoxMassProperties(
            Eigen::Vector3d(1.0, -2.0, 3.0),
            Eigen::Vector3d(2.0, 3.0, 4.0),
            2.0
        );

    if (
        !checkProperties(
            denseAnisotropic,
            48.0,
            Eigen::Vector3d(48.0, -96.0, 144.0),
            Eigen::Vector3d(1.0, -2.0, 3.0),
            error
        )
    )
    {
        return fail(scenario, "box density scaling: " + error);
    }

    const std::array<Eigen::Vector3d, 4> centers =
    {{
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 2.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 3.0)
    }};

    MassProperties aggregate =
        sumUnitVoxelContributions(centers);

    if (
        !checkProperties(
            aggregate,
            4.0,
            Eigen::Vector3d(1.0, 2.0, 3.0),
            Eigen::Vector3d(0.25, 0.5, 0.75),
            error
        )
    )
    {
        return fail(scenario, "four-voxel aggregate: " + error);
    }

    const std::array<std::size_t, 4> deletionOrder =
    {{2, 0, 3, 1}};

    const std::array<double, 4> expectedMass =
    {{3.0, 2.0, 1.0, 0.0}};

    const std::array<Eigen::Vector3d, 4> expectedFirstMoment =
    {{
        Eigen::Vector3d(1.0, 0.0, 3.0),
        Eigen::Vector3d(1.0, 0.0, 3.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d::Zero()
    }};

    const std::array<Eigen::Vector3d, 3> expectedCenterOfMass =
    {{
        Eigen::Vector3d(1.0 / 3.0, 0.0, 1.0),
        Eigen::Vector3d(0.5, 0.0, 1.5),
        Eigen::Vector3d(1.0, 0.0, 0.0)
    }};

    for (std::size_t step = 0; step < deletionOrder.size(); ++step)
    {
        const MassProperties contribution =
            computeAxisAlignedBoxMassProperties(
                centers[deletionOrder[step]],
                Eigen::Vector3d::Ones(),
                1.0
            );

        aggregate.mass -= contribution.mass;
        aggregate.firstMoment -= contribution.firstMoment;

        if (
            !approximatelyEqual(
                aggregate.mass,
                expectedMass[step]
            ) ||
            !approximatelyEqual(
                aggregate.firstMoment,
                expectedFirstMoment[step]
            )
        )
        {
            return fail(
                scenario,
                "incremental deletion changed mass or first moment incorrectly"
            );
        }

        Eigen::Vector3d centerOfMass =
            Eigen::Vector3d::Zero();

        if (step + 1 < deletionOrder.size())
        {
            if (
                !aggregate.tryCenterOfMass(centerOfMass) ||
                !approximatelyEqual(
                    centerOfMass,
                    expectedCenterOfMass[step]
                )
            )
            {
                return fail(
                    scenario,
                    "incremental deletion produced the wrong center of mass"
                );
            }
        }
        else if (aggregate.tryCenterOfMass(centerOfMass))
        {
            return fail(
                scenario,
                "zero-mass aggregate exposed a center of mass"
            );
        }
    }

    MassProperties roundTrip =
        sumUnitVoxelContributions(centers);

    const MassProperties removed =
        computeAxisAlignedBoxMassProperties(
            centers[2],
            Eigen::Vector3d::Ones(),
            1.0
        );

    roundTrip.mass -= removed.mass;
    roundTrip.firstMoment -= removed.firstMoment;
    roundTrip.mass += removed.mass;
    roundTrip.firstMoment += removed.firstMoment;

    if (
        !checkProperties(
            roundTrip,
            4.0,
            Eigen::Vector3d(1.0, 2.0, 3.0),
            Eigen::Vector3d(0.25, 0.5, 0.75),
            error
        )
    )
    {
        return fail(scenario, "remove/add round trip: " + error);
    }

    bool rejectedInvalidSize = false;

    try
    {
        static_cast<void>(
            computeAxisAlignedBoxMassProperties(
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d(1.0, 0.0, 1.0)
            )
        );
    }
    catch (const std::invalid_argument&)
    {
        rejectedInvalidSize = true;
    }

    if (!rejectedInvalidSize)
    {
        return fail(
            scenario,
            "non-positive voxel extent was not rejected"
        );
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runDeletionScenario()
{
    const std::string scenario = "deletion";

    std::array<Eigen::Vector3d, 8> centers;
    std::size_t centerIndex = 0;

    for (int xSign = -1; xSign <= 1; xSign += 2)
    {
        for (int ySign = -1; ySign <= 1; ySign += 2)
        {
            for (int zSign = -1; zSign <= 1; zSign += 2)
            {
                centers[centerIndex] =
                    0.5 * Eigen::Vector3d(
                        static_cast<double>(xSign),
                        static_cast<double>(ySign),
                        static_cast<double>(zSign)
                    );

                ++centerIndex;
            }
        }
    }

    std::array<MassProperties, 8> contributions =
    {{
        computeAxisAlignedBoxMassProperties(
            centers[0], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[1], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[2], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[3], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[4], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[5], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[6], Eigen::Vector3d::Ones(), 1.0
        ),
        computeAxisAlignedBoxMassProperties(
            centers[7], Eigen::Vector3d::Ones(), 1.0
        )
    }};

    std::array<int, 8> permutation;
    std::iota(
        permutation.begin(),
        permutation.end(),
        0
    );

    std::size_t permutationCount = 0;

    do
    {
        MassProperties incremental = contributions[0];

        for (std::size_t index = 1; index < contributions.size(); ++index)
        {
            incremental.mass += contributions[index].mass;
            incremental.firstMoment +=
                contributions[index].firstMoment;
        }

        std::array<bool, 8> alive;
        alive.fill(true);

        for (std::size_t step = 0; step < permutation.size(); ++step)
        {
            const int removedIndex = permutation[step];

            incremental.mass -=
                contributions[removedIndex].mass;

            incremental.firstMoment -=
                contributions[removedIndex].firstMoment;

            alive[removedIndex] = false;

            int expectedCount = 0;
            Eigen::Vector3d expectedFirstMoment =
                Eigen::Vector3d::Zero();

            for (std::size_t index = 0; index < alive.size(); ++index)
            {
                if (alive[index])
                {
                    ++expectedCount;
                    expectedFirstMoment += centers[index];
                }
            }

            if (
                !approximatelyEqual(
                    incremental.mass,
                    static_cast<double>(expectedCount)
                ) ||
                !approximatelyEqual(
                    incremental.firstMoment,
                    expectedFirstMoment
                )
            )
            {
                return fail(
                    scenario,
                    "incremental state disagrees with an independent full sum"
                );
            }

            Eigen::Vector3d centerOfMass =
                Eigen::Vector3d::Zero();

            if (expectedCount > 0)
            {
                const Eigen::Vector3d expectedCenterOfMass =
                    expectedFirstMoment /
                    static_cast<double>(expectedCount);

                if (
                    !incremental.tryCenterOfMass(centerOfMass) ||
                    !approximatelyEqual(
                        centerOfMass,
                        expectedCenterOfMass
                    )
                )
                {
                    return fail(
                        scenario,
                        "deletion prefix produced the wrong center of mass"
                    );
                }
            }
            else if (incremental.tryCenterOfMass(centerOfMass))
            {
                return fail(
                    scenario,
                    "fully deleted fixture exposed a center of mass"
                );
            }
        }

        ++permutationCount;
    }
    while (
        std::next_permutation(
            permutation.begin(),
            permutation.end()
        )
    );

    if (permutationCount != 40320)
    {
        return fail(
            scenario,
            "did not enumerate all 8! deletion orders"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " permutations=" << permutationCount
        << '\n';

    return EXIT_SUCCESS;
}


int runEnergyScenario()
{
    const std::string scenario = "energy";

    const Eigen::Vector3d centerOfMass(6.0, 6.0, 15.0);
    const Eigen::Vector3d target(1.0, 2.0, 3.0);
    const Eigen::Vector3d delta = centerOfMass - target;

    const Eigen::Vector3d projected =
        projectPerpendicularToGravity(
            delta,
            Eigen::Vector3d(0.0, -1.0, 0.0)
        );

    if (
        !approximatelyEqual(
            projected,
            Eigen::Vector3d(5.0, 0.0, 12.0)
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassDistance(
                centerOfMass,
                target,
                Eigen::Vector3d(0.0, -1.0, 0.0)
            ),
            13.0
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassEnergy(
                centerOfMass,
                target,
                Eigen::Vector3d(0.0, -1.0, 0.0)
            ),
            84.5
        )
    )
    {
        return fail(
            scenario,
            "axis-aligned gravity projection is incorrect"
        );
    }

    const Eigen::Vector3d scaledGravity(0.0, -7.0, 0.0);

    if (
        !approximatelyEqual(
            projectPerpendicularToGravity(
                delta,
                scaledGravity
            ),
            projected
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassDistance(
                centerOfMass,
                target,
                scaledGravity
            ),
            13.0
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassEnergy(
                centerOfMass,
                target,
                scaledGravity
            ),
            84.5
        )
    )
    {
        return fail(
            scenario,
            "projection depends on gravity-vector magnitude"
        );
    }

    const Eigen::Vector3d obliqueDelta(3.0, 1.0, 2.0);
    const Eigen::Vector3d obliqueGravity(1.0, 1.0, 0.0);
    const Eigen::Vector3d obliqueProjected(1.0, -1.0, 2.0);

    if (
        !approximatelyEqual(
            projectPerpendicularToGravity(
                obliqueDelta,
                obliqueGravity
            ),
            obliqueProjected
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassDistance(
                obliqueDelta,
                Eigen::Vector3d::Zero(),
                obliqueGravity
            ),
            std::sqrt(6.0)
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassEnergy(
                obliqueDelta,
                Eigen::Vector3d::Zero(),
                obliqueGravity
            ),
            3.0
        )
    )
    {
        return fail(
            scenario,
            "oblique gravity projection or squared energy is incorrect"
        );
    }

    const Eigen::Vector3d parallelDelta(2.0, 2.0, 0.0);

    if (
        !approximatelyEqual(
            projectPerpendicularToGravity(
                parallelDelta,
                obliqueGravity
            ),
            Eigen::Vector3d::Zero()
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassDistance(
                parallelDelta,
                Eigen::Vector3d::Zero(),
                obliqueGravity
            ),
            0.0
        ) ||
        !approximatelyEqual(
            projectedCenterOfMassEnergy(
                parallelDelta,
                Eigen::Vector3d::Zero(),
                obliqueGravity
            ),
            0.0
        )
    )
    {
        return fail(
            scenario,
            "gravity-parallel displacement did not project to zero"
        );
    }

    const double distance =
        projectedCenterOfMassDistance(
            obliqueDelta,
            Eigen::Vector3d::Zero(),
            obliqueGravity
        );

    const double energy =
        projectedCenterOfMassEnergy(
            obliqueDelta,
            Eigen::Vector3d::Zero(),
            obliqueGravity
        );

    if (!approximatelyEqual(energy, 0.5 * distance * distance))
    {
        return fail(
            scenario,
            "energy is not one half of the squared projected distance"
        );
    }

    const double maximum =
        std::numeric_limits<double>::max();

    const std::array<Eigen::Vector3d, 4> extremeGravities =
    {{
        Eigen::Vector3d(0.0, -1e-300, 0.0),
        Eigen::Vector3d(0.0, -1e308, 0.0),
        Eigen::Vector3d(1e-300, -1e-300, 1e-300),
        Eigen::Vector3d(maximum, -maximum, maximum)
    }};

    const std::array<Eigen::Vector3d, 4> referenceGravities =
    {{
        Eigen::Vector3d(0.0, -1.0, 0.0),
        Eigen::Vector3d(0.0, -1.0, 0.0),
        Eigen::Vector3d(1.0, -1.0, 1.0),
        Eigen::Vector3d(1.0, -1.0, 1.0)
    }};

    for (std::size_t index = 0;
         index < extremeGravities.size();
         ++index)
    {
        const Eigen::Vector3d normalized =
            normalizedGravityDirection(
                extremeGravities[index]
            );

        const Eigen::Vector3d referenceDirection =
            normalizedGravityDirection(
                referenceGravities[index]
            );

        const Eigen::Vector3d extremeProjection =
            projectPerpendicularToGravity(
                delta,
                extremeGravities[index]
            );

        const Eigen::Vector3d referenceProjection =
            projectPerpendicularToGravity(
                delta,
                referenceGravities[index]
            );

        const double extremeDistance =
            projectedCenterOfMassDistance(
                centerOfMass,
                target,
                extremeGravities[index]
            );

        const double extremeEnergy =
            projectedCenterOfMassEnergy(
                centerOfMass,
                target,
                extremeGravities[index]
            );

        if (
            !normalized.allFinite() ||
            !approximatelyEqual(normalized.norm(), 1.0) ||
            !approximatelyEqual(normalized, referenceDirection) ||
            !approximatelyEqual(
                extremeProjection,
                referenceProjection
            ) ||
            !approximatelyEqual(
                extremeProjection.dot(normalized),
                0.0
            ) ||
            !approximatelyEqual(
                extremeDistance,
                referenceProjection.norm()
            ) ||
            !approximatelyEqual(
                extremeEnergy,
                0.5 * extremeDistance * extremeDistance
            )
        )
        {
            return fail(
                scenario,
                "extreme gravity scaling changed its physical direction"
            );
        }
    }

    const double infinity =
        std::numeric_limits<double>::infinity();
    const double notANumber =
        std::numeric_limits<double>::quiet_NaN();

    const std::array<Eigen::Vector3d, 4> invalidGravities =
    {{
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidGravity : invalidGravities)
    {
        bool rejectedGravity = false;

        try
        {
            static_cast<void>(
                normalizedGravityDirection(invalidGravity)
            );
        }
        catch (const std::invalid_argument&)
        {
            rejectedGravity = true;
        }

        if (!rejectedGravity)
        {
            return fail(
                scenario,
                "zero or non-finite gravity was not rejected"
            );
        }
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}

} // namespace


int main(
    int argc,
    char* argv[]
)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: mass-properties-regression "
            << "<surface|voxel|deletion|energy>\n";

        return EXIT_FAILURE;
    }

    const std::string scenario = argv[1];

    if (scenario == "surface")
    {
        return runSurfaceScenario();
    }

    if (scenario == "voxel")
    {
        return runVoxelScenario();
    }

    if (scenario == "deletion")
    {
        return runDeletionScenario();
    }

    if (scenario == "energy")
    {
        return runEnergyScenario();
    }

    return fail(
        scenario,
        "unknown scenario; expected surface, voxel, deletion or energy"
    );
}
