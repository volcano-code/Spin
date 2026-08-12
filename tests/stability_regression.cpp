#include "physics/BalanceEnergy.h"
#include "physics/BalanceObjective.h"
#include "stability/SupportRegion.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


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


template <typename Callable>
bool rejectsInvalidArgument(
    Callable&& callable
)
{
    try
    {
        callable();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }

    return false;
}


template <typename Callable>
bool rejectsDomainError(
    Callable&& callable
)
{
    try
    {
        callable();
    }
    catch (const std::domain_error&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }

    return false;
}


double signedPolygonArea(
    const Eigen::MatrixXd& polygon
)
{
    if (polygon.rows() < 3 || polygon.cols() != 2)
    {
        return 0.0;
    }

    double twiceArea = 0.0;

    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Vector2d first =
            polygon.row(index).transpose();
        const Eigen::Vector2d second =
            polygon.row((index + 1) % polygon.rows()).transpose();

        twiceArea +=
            first(0) * second(1) -
            first(1) * second(0);
    }

    return 0.5 * twiceArea;
}


bool containsEquivalentPoint(
    const Eigen::MatrixXd& polygon,
    const Eigen::Vector2d& expected
)
{
    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Vector2d point =
            polygon.row(index).transpose();

        if (
            approximatelyEqual(point(0), expected(0)) &&
            approximatelyEqual(point(1), expected(1))
        )
        {
            return true;
        }
    }

    return false;
}


bool equivalentPointSets(
    const Eigen::MatrixXd& first,
    const Eigen::MatrixXd& second
)
{
    if (
        first.cols() != 2 ||
        second.cols() != 2 ||
        first.rows() != second.rows()
    )
    {
        return false;
    }

    for (Eigen::Index index = 0; index < first.rows(); ++index)
    {
        const Eigen::Vector2d point =
            first.row(index).transpose();

        if (!containsEquivalentPoint(second, point))
        {
            return false;
        }
    }

    return true;
}


template <std::size_t Count>
Eigen::MatrixXd pointsToMatrix(
    const std::array<Eigen::Vector2d, Count>& points
)
{
    Eigen::MatrixXd result(
        static_cast<Eigen::Index>(Count),
        2
    );

    for (std::size_t index = 0; index < Count; ++index)
    {
        result.row(static_cast<Eigen::Index>(index)) =
            points[index].transpose();
    }

    return result;
}


Eigen::Vector3d pointAboveSupport(
    const stability::SupportRegion& region,
    const Eigen::Vector2d& planePoint,
    double height
)
{
    return
        region.fromPlaneCoordinates(planePoint) -
        height * region.gravity();
}


Eigen::MatrixXd makeSquareSupportMesh()
{
    Eigen::MatrixXd vertices(9, 3);

    vertices <<
        -2.0, 0.0, -2.0,
         2.0, 0.0, -2.0,
         2.0, 0.0,  2.0,
        -2.0, 0.0,  2.0,
        -2.0, 0.0, -2.0, // duplicate hull point
         0.0, 0.0,  0.0, // interior point
         0.0, 0.0, -2.0, // collinear edge point
         1.0, 0.25, 1.0, // exactly on the contact threshold
         0.0, 0.30, 0.0; // not a contact point

    return vertices;
}


int runSupportScenario()
{
    const std::string scenario = "support";
    const Eigen::Vector3d gravity(0.0, -1.0, 0.0);
    const double contactThreshold = 0.25;
    const Eigen::MatrixXd vertices = makeSquareSupportMesh();

    const stability::SupportRegion region =
        stability::SupportRegion::fromMesh(
            vertices,
            gravity,
            contactThreshold
        );

    if (
        !approximatelyEqual(region.gravity(), gravity) ||
        !approximatelyEqual(region.gravity().norm(), 1.0) ||
        !approximatelyEqual(region.basisU().norm(), 1.0) ||
        !approximatelyEqual(region.basisV().norm(), 1.0) ||
        !approximatelyEqual(region.gravity().dot(region.basisU()), 0.0) ||
        !approximatelyEqual(region.gravity().dot(region.basisV()), 0.0) ||
        !approximatelyEqual(region.basisU().dot(region.basisV()), 0.0) ||
        !approximatelyEqual(region.supportHeight(), 0.0)
    )
    {
        return fail(
            scenario,
            "support frame is not orthonormal or has the wrong height"
        );
    }

    const std::vector<int>& contactIndices =
        region.contactVertexIndices();

    if (contactIndices.size() != 8)
    {
        return fail(
            scenario,
            "contact threshold is not inclusive or selected a non-contact point"
        );
    }

    for (int expectedIndex = 0; expectedIndex < 8; ++expectedIndex)
    {
        if (
            std::find(
                contactIndices.begin(),
                contactIndices.end(),
                expectedIndex
            ) == contactIndices.end()
        )
        {
            return fail(
                scenario,
                "the expected inclusive contact set was not retained"
            );
        }
    }

    if (
        std::find(
            contactIndices.begin(),
            contactIndices.end(),
            8
        ) != contactIndices.end()
    )
    {
        return fail(
            scenario,
            "a vertex above the contact threshold was selected"
        );
    }

    const auto& supportPolygon = region.supportPolygon2D();

    if (
        supportPolygon.rows() != 4 ||
        !approximatelyEqual(signedPolygonArea(supportPolygon), 16.0) ||
        !region.contains2D(Eigen::Vector2d::Zero())
    )
    {
        return fail(
            scenario,
            "duplicate/interior/collinear contacts changed the square hull"
        );
    }

    for (int index = 0; index < 4; ++index)
    {
        const Eigen::Vector2d corner =
            region.toPlaneCoordinates(vertices.row(index).transpose());

        if (!region.contains2D(corner))
        {
            return fail(
                scenario,
                "closed support containment rejected a boundary vertex"
            );
        }
    }

    const Eigen::Vector3d roundTripPoint(0.75, 0.0, -1.25);

    if (
        !approximatelyEqual(
            region.fromPlaneCoordinates(
                region.toPlaneCoordinates(roundTripPoint)
            ),
            roundTripPoint
        )
    )
    {
        return fail(
            scenario,
            "support-plane coordinate round trip is incorrect"
        );
    }

    const Eigen::Vector2d outsidePoint =
        region.toPlaneCoordinates(Eigen::Vector3d(2.25, 0.0, 0.0));

    if (region.contains2D(outsidePoint))
    {
        return fail(
            scenario,
            "support containment accepted a point outside the hull"
        );
    }

    const stability::SupportRegion scaledGravityRegion =
        stability::SupportRegion::fromMesh(
            vertices,
            Eigen::Vector3d(0.0, -1e300, 0.0),
            contactThreshold
        );

    if (
        scaledGravityRegion.contactVertexIndices().size() != 8 ||
        !approximatelyEqual(
            signedPolygonArea(
                scaledGravityRegion.supportPolygon2D()
            ),
            16.0
        ) ||
        !approximatelyEqual(
            scaledGravityRegion.maximumInset(),
            2.0
        )
    )
    {
        return fail(
            scenario,
            "support construction depends on gravity-vector magnitude"
        );
    }

    const std::array<int, 9> permutation =
    {{
        7, 2, 5, 0, 8, 3, 6, 1, 4
    }};

    Eigen::MatrixXd shuffledVertices(vertices.rows(), vertices.cols());

    for (Eigen::Index row = 0; row < shuffledVertices.rows(); ++row)
    {
        shuffledVertices.row(row) =
            vertices.row(permutation[static_cast<std::size_t>(row)]);
    }

    const stability::SupportRegion shuffledRegion =
        stability::SupportRegion::fromMesh(
            shuffledVertices,
            gravity,
            contactThreshold
        );

    if (
        shuffledRegion.contactVertexIndices().size() != 8 ||
        !approximatelyEqual(
            signedPolygonArea(shuffledRegion.supportPolygon2D()),
            16.0
        ) ||
        !approximatelyEqual(shuffledRegion.maximumInset(), 2.0)
    )
    {
        return fail(
            scenario,
            "support geometry depends on mesh vertex ordering"
        );
    }

    const Eigen::Vector3d translation(3.0, 5.0, -4.0);
    Eigen::MatrixXd translatedVertices = vertices;
    translatedVertices.rowwise() += translation.transpose();

    const stability::SupportRegion translatedRegion =
        stability::SupportRegion::fromMesh(
            translatedVertices,
            gravity,
            contactThreshold
        );

    if (
        !approximatelyEqual(translatedRegion.supportHeight(), 5.0) ||
        !approximatelyEqual(
            signedPolygonArea(translatedRegion.supportPolygon2D()),
            16.0
        ) ||
        !translatedRegion.contains2D(
            translatedRegion.toPlaneCoordinates(
                roundTripPoint + translation
            )
        )
    )
    {
        return fail(
            scenario,
            "support geometry is not translation invariant"
        );
    }

    const Eigen::Vector3d largeTranslation(1e12, 3e12, -2e12);
    Eigen::MatrixXd farTranslatedVertices = vertices;
    farTranslatedVertices.rowwise() += largeTranslation.transpose();

    const stability::SupportRegion farTranslatedRegion =
        stability::SupportRegion::fromMesh(
            farTranslatedVertices,
            gravity,
            contactThreshold
        );

    if (
        farTranslatedRegion.supportPolygon2D().rows() != 4 ||
        !approximatelyEqual(
            signedPolygonArea(
                farTranslatedRegion.supportPolygon2D()
            ),
            16.0
        ) ||
        !approximatelyEqual(
            farTranslatedRegion.maximumInset(),
            2.0
        )
    )
    {
        return fail(
            scenario,
            "large world coordinates changed the local support geometry"
        );
    }

    const Eigen::Matrix3d rotation =
        Eigen::AngleAxisd(
            0.71,
            Eigen::Vector3d(1.0, 2.0, -1.0).normalized()
        ).toRotationMatrix();

    Eigen::MatrixXd rotatedVertices(vertices.rows(), 3);

    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        rotatedVertices.row(row) =
            (
                rotation * vertices.row(row).transpose()
            ).transpose();
    }

    const stability::SupportRegion rotatedRegion =
        stability::SupportRegion::fromMesh(
            rotatedVertices,
            rotation * gravity,
            contactThreshold
        );

    if (
        rotatedRegion.contactVertexIndices().size() != 8 ||
        !approximatelyEqual(
            signedPolygonArea(rotatedRegion.supportPolygon2D()),
            16.0
        ) ||
        !rotatedRegion.contains2D(
            rotatedRegion.toPlaneCoordinates(
                rotation * roundTripPoint
            )
        )
    )
    {
        return fail(
            scenario,
            "support geometry is not rotation invariant"
        );
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
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        stability::SupportRegion::fromMesh(
                            vertices,
                            invalidGravity,
                            contactThreshold
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a zero or non-finite gravity was accepted"
            );
        }
    }

    const std::array<double, 3> invalidThresholds =
    {{
        -1.0,
        infinity,
        notANumber
    }};

    for (double invalidThreshold : invalidThresholds)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        stability::SupportRegion::fromMesh(
                            vertices,
                            gravity,
                            invalidThreshold
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a negative or non-finite contact threshold was accepted"
            );
        }
    }

    const Eigen::MatrixXd emptyMesh(0, 3);

    if (
        !rejectsInvalidArgument(
            [&]()
            {
                static_cast<void>(
                    stability::SupportRegion::fromMesh(
                        emptyMesh,
                        gravity,
                        contactThreshold
                    )
                );
            }
        )
    )
    {
        return fail(
            scenario,
            "an empty vertex matrix was not rejected as invalid input"
        );
    }

    Eigen::MatrixXd twoUniquePoints(3, 3);
    twoUniquePoints <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0;
    Eigen::MatrixXd collinearPoints(4, 3);
    collinearPoints <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        3.0, 0.0, 0.0;
    Eigen::MatrixXd nonFiniteMesh = vertices;
    nonFiniteMesh(2, 1) = notANumber;

    if (
        !rejectsInvalidArgument(
            [&]()
            {
                static_cast<void>(
                    stability::SupportRegion::fromMesh(
                        nonFiniteMesh,
                        gravity,
                        contactThreshold
                    )
                );
            }
        )
    )
    {
        return fail(
            scenario,
            "a non-finite vertex matrix was not rejected as invalid input"
        );
    }

    const std::array<Eigen::MatrixXd, 2> degenerateMeshes =
    {{
        twoUniquePoints,
        collinearPoints
    }};

    for (const Eigen::MatrixXd& degenerateMesh : degenerateMeshes)
    {
        if (
            !rejectsDomainError(
                [&]()
                {
                    static_cast<void>(
                        stability::SupportRegion::fromMesh(
                            degenerateMesh,
                            gravity,
                            contactThreshold
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a point/line support was not rejected as geometric degeneracy"
            );
        }
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " contacts=8 hull=4 area=16"
        << '\n';

    return EXIT_SUCCESS;
}


int runSafeScenario()
{
    const std::string scenario = "safe";
    const Eigen::Vector3d gravity(0.0, -1.0, 0.0);

    const stability::SupportRegion region =
        stability::SupportRegion::fromMesh(
            makeSquareSupportMesh(),
            gravity,
            0.25
        );

    const double gamma = std::atan(0.5);

    const std::array<Eigen::Vector2d, 4> expectedSafeCorners =
    {{
        Eigen::Vector2d(-1.0, -1.0),
        Eigen::Vector2d( 1.0, -1.0),
        Eigen::Vector2d( 1.0,  1.0),
        Eigen::Vector2d(-1.0,  1.0)
    }};

    struct StateFixture
    {
        Eigen::Vector2d projectedCenterOfMass;
        stability::StandingState expectedState;
        Eigen::Vector2d expectedTarget;
        double expectedEnergy;
    };

    const std::array<StateFixture, 5> stateFixtures =
    {{
        {
            Eigen::Vector2d(2.5, 0.0),
            stability::StandingState::OutsideSupport,
            Eigen::Vector2d(1.0, 0.0),
            1.125
        },
        {
            Eigen::Vector2d(2.0, 0.0),
            stability::StandingState::FeasibleButUnsafe,
            Eigen::Vector2d(1.0, 0.0),
            0.5
        },
        {
            Eigen::Vector2d(1.5, 0.0),
            stability::StandingState::FeasibleButUnsafe,
            Eigen::Vector2d(1.0, 0.0),
            0.125
        },
        {
            Eigen::Vector2d(1.0, 0.0),
            stability::StandingState::Stable,
            Eigen::Vector2d(1.0, 0.0),
            0.0
        },
        {
            Eigen::Vector2d(0.0, 0.0),
            stability::StandingState::Stable,
            Eigen::Vector2d(0.0, 0.0),
            0.0
        }
    }};

    for (const StateFixture& fixture : stateFixtures)
    {
        const Eigen::Vector3d centerOfMass =
            pointAboveSupport(
                region,
                fixture.projectedCenterOfMass,
                2.0
            );

        const stability::StandingTarget target =
            stability::makeStandingTarget(
                region,
                centerOfMass,
                gamma
            );

        if (
            target.state() != fixture.expectedState ||
            !target.requestedMarginAttainable() ||
            !approximatelyEqual(target.height(), 2.0) ||
            !approximatelyEqual(target.requestedInset(), 1.0) ||
            !approximatelyEqual(target.maximumInset(), 2.0) ||
            !approximatelyEqual(
                target.projectedCenterOfMass2D()(0),
                fixture.projectedCenterOfMass(0)
            ) ||
            !approximatelyEqual(
                target.projectedCenterOfMass2D()(1),
                fixture.projectedCenterOfMass(1)
            ) ||
            !approximatelyEqual(
                target.target2D()(0),
                fixture.expectedTarget(0)
            ) ||
            !approximatelyEqual(
                target.target2D()(1),
                fixture.expectedTarget(1)
            ) ||
            !approximatelyEqual(
                target.objective().target(),
                region.fromPlaneCoordinates(fixture.expectedTarget)
            ) ||
            !approximatelyEqual(
                target.objective().gravity(),
                gravity
            ) ||
            !approximatelyEqual(
                projectedCenterOfMassEnergy(
                    centerOfMass,
                    target.objective().target(),
                    target.objective().gravity()
                ),
                fixture.expectedEnergy
            )
        )
        {
            return fail(
                scenario,
                "outside/feasible/stable state or standing target is incorrect"
            );
        }

        if (
            target.safeRegion2D().rows() != 4 ||
            !equivalentPointSets(
                target.safeRegion2D(),
                pointsToMatrix(expectedSafeCorners)
            )
        )
        {
            return fail(
                scenario,
                "toppling-angle inset did not produce [-1,1]^2"
            );
        }
    }

    const Eigen::Vector3d supportBoundaryCenterOfMass =
        pointAboveSupport(
            region,
            Eigen::Vector2d(2.0, 0.0),
            2.0
        );

    const Eigen::Vector3d feasibleCenterOfMass =
        pointAboveSupport(
            region,
            Eigen::Vector2d(1.5, 0.0),
            2.0
        );

    const Eigen::Vector3d safeBoundaryCenterOfMass =
        pointAboveSupport(
            region,
            Eigen::Vector2d(1.0, 0.0),
            2.0
        );

    const Eigen::Vector3d centeredCenterOfMass =
        pointAboveSupport(
            region,
            Eigen::Vector2d::Zero(),
            2.0
        );

    if (
        !approximatelyEqual(
            region.topplingAngle(supportBoundaryCenterOfMass),
            0.0
        ) ||
        !approximatelyEqual(
            region.topplingAngle(feasibleCenterOfMass),
            std::atan(0.25)
        ) ||
        !approximatelyEqual(
            region.topplingAngle(safeBoundaryCenterOfMass),
            gamma
        ) ||
        !approximatelyEqual(
            region.topplingAngle(centeredCenterOfMass),
            std::atan(1.0)
        )
    )
    {
        return fail(
            scenario,
            "toppling angle disagrees with atan(margin/height)"
        );
    }

    if (
        !rejectsInvalidArgument(
            [&]()
            {
                static_cast<void>(
                    region.topplingAngle(
                        pointAboveSupport(
                            region,
                            Eigen::Vector2d(3.0, 0.0),
                            -1.0
                        )
                    )
                );
            }
        )
    )
    {
        return fail(
            scenario,
            "an outside point below the support plane bypassed height validation"
        );
    }

    /*
     * Paper contract: gamma=0 means no inset. The beta's hard-coded 1e-4
     * erosion must not leak into the paper-conformant implementation.
     */
    const Eigen::Vector2d zeroGammaProjection(1.5, -0.25);
    const Eigen::Vector3d zeroGammaCenterOfMass =
        pointAboveSupport(
            region,
            zeroGammaProjection,
            2.0
        );

    const stability::StandingTarget zeroGammaTarget =
        stability::makeStandingTarget(
            region,
            zeroGammaCenterOfMass,
            0.0
        );

    if (
        zeroGammaTarget.state() != stability::StandingState::Stable ||
        !approximatelyEqual(zeroGammaTarget.requestedInset(), 0.0) ||
        !equivalentPointSets(
            zeroGammaTarget.safeRegion2D(),
            region.supportPolygon2D()
        ) ||
        !approximatelyEqual(
            zeroGammaTarget.target2D()(0),
            zeroGammaProjection(0)
        ) ||
        !approximatelyEqual(
            zeroGammaTarget.target2D()(1),
            zeroGammaProjection(1)
        )
    )
    {
        return fail(
            scenario,
            "gamma=0 eroded the support or moved an interior target"
        );
    }

    const double infinity =
        std::numeric_limits<double>::infinity();
    const double notANumber =
        std::numeric_limits<double>::quiet_NaN();

    const std::array<double, 4> invalidAngles =
    {{
        -1e-3,
        0.5 * std::acos(-1.0),
        infinity,
        notANumber
    }};

    for (double invalidAngle : invalidAngles)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        stability::makeStandingTarget(
                            region,
                            centeredCenterOfMass,
                            invalidAngle
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a negative, right, or non-finite toppling angle was accepted"
            );
        }
    }

    const std::array<Eigen::Vector3d, 4> invalidCentersOfMass =
    {{
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber),
        pointAboveSupport(
            region,
            Eigen::Vector2d::Zero(),
            -1.0
        )
    }};

    for (const Eigen::Vector3d& invalidCenterOfMass : invalidCentersOfMass)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        stability::makeStandingTarget(
                            region,
                            invalidCenterOfMass,
                            gamma
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a non-finite or below-support center of mass was accepted"
            );
        }
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " inset=1 maximum=2 boundary=closed"
        << '\n';

    return EXIT_SUCCESS;
}


int runStandingTargetScenario()
{
    const std::string scenario = "standing-target";
    const Eigen::Vector3d gravity(0.0, -1.0, 0.0);

    const stability::SupportRegion square =
        stability::SupportRegion::fromMesh(
            makeSquareSupportMesh(),
            gravity,
            0.25
        );

    const double squareGamma = std::atan(0.5);

    struct ProjectionFixture
    {
        Eigen::Vector2d query;
        Eigen::Vector2d expectedTarget;
        double expectedSquaredDistance;
    };

    const std::array<ProjectionFixture, 5> projectionFixtures =
    {{
        {
            Eigen::Vector2d(0.2, -0.4),
            Eigen::Vector2d(0.2, -0.4),
            0.0
        },
        {
            Eigen::Vector2d(1.0, 0.3),
            Eigen::Vector2d(1.0, 0.3),
            0.0
        },
        {
            Eigen::Vector2d(3.0, 0.0),
            Eigen::Vector2d(1.0, 0.0),
            4.0
        },
        {
            Eigen::Vector2d(3.0, 4.0),
            Eigen::Vector2d(1.0, 1.0),
            13.0
        },
        {
            Eigen::Vector2d(-2.0, -3.0),
            Eigen::Vector2d(-1.0, -1.0),
            5.0
        }
    }};

    for (const ProjectionFixture& fixture : projectionFixtures)
    {
        const Eigen::Vector3d centerOfMass =
            pointAboveSupport(square, fixture.query, 2.0);

        const stability::StandingTarget target =
            stability::makeStandingTarget(
                square,
                centerOfMass,
                squareGamma
            );

        if (
            !approximatelyEqual(
                target.target2D()(0),
                fixture.expectedTarget(0)
            ) ||
            !approximatelyEqual(
                target.target2D()(1),
                fixture.expectedTarget(1)
            ) ||
            !approximatelyEqual(
                (
                    target.projectedCenterOfMass2D() -
                    target.target2D()
                ).squaredNorm(),
                fixture.expectedSquaredDistance
            ) ||
            !approximatelyEqual(
                aggregateProjectedCenterOfMassEnergy(
                    centerOfMass,
                    BalanceObjectives{target.objective()}
                ),
                0.5 * fixture.expectedSquaredDistance
            )
        )
        {
            return fail(
                scenario,
                "point/edge/vertex projection onto the closed safe set is incorrect"
            );
        }
    }

    /*
     * The first two fixtures also lock the paper/beta difference: points
     * already inside the safe region, including its boundary, do not move.
     */

    Eigen::MatrixXd triangleVertices(3, 3);
    triangleVertices <<
        0.0, 0.0, 0.0,
        4.0, 0.0, 0.0,
        0.0, 0.0, 3.0;

    const stability::SupportRegion triangle =
        stability::SupportRegion::fromMesh(
            triangleVertices,
            gravity,
            0.0
        );

    const Eigen::Vector2d triangleQuery =
        triangle.toPlaneCoordinates(Eigen::Vector3d(3.0, 0.0, 2.0));

    const stability::StandingTarget halfInsetTriangle =
        stability::makeStandingTarget(
            triangle,
            pointAboveSupport(triangle, triangleQuery, 1.0),
            std::atan(0.5)
        );

    const std::array<Eigen::Vector2d, 3> expectedHalfInsetTriangle =
    {{
        triangle.toPlaneCoordinates(Eigen::Vector3d(0.5, 0.0, 0.5)),
        triangle.toPlaneCoordinates(Eigen::Vector3d(2.5, 0.0, 0.5)),
        triangle.toPlaneCoordinates(Eigen::Vector3d(0.5, 0.0, 2.0))
    }};

    if (
        !approximatelyEqual(triangle.maximumInset(), 1.0) ||
        !halfInsetTriangle.requestedMarginAttainable() ||
        !approximatelyEqual(halfInsetTriangle.requestedInset(), 0.5) ||
        !equivalentPointSets(
            halfInsetTriangle.safeRegion2D(),
            pointsToMatrix(expectedHalfInsetTriangle)
        )
    )
    {
        return fail(
            scenario,
            "triangle inset is not the parallel half-plane offset"
        );
    }

    const Eigen::Vector2d triangleIncenter =
        triangle.toPlaneCoordinates(Eigen::Vector3d(1.0, 0.0, 1.0));

    const stability::StandingTarget maximumTriangleInset =
        stability::makeStandingTarget(
            triangle,
            pointAboveSupport(triangle, triangleQuery, 1.0),
            0.25 * std::acos(-1.0)
        );

    if (
        !maximumTriangleInset.requestedMarginAttainable() ||
        maximumTriangleInset.safeRegion2D().rows() != 1 ||
        !approximatelyEqual(
            maximumTriangleInset.safeRegion2D()(0, 0),
            triangleIncenter(0)
        ) ||
        !approximatelyEqual(
            maximumTriangleInset.safeRegion2D()(0, 1),
            triangleIncenter(1)
        ) ||
        !approximatelyEqual(
            maximumTriangleInset.target2D()(0),
            triangleIncenter(0)
        ) ||
        !approximatelyEqual(
            maximumTriangleInset.target2D()(1),
            triangleIncenter(1)
        )
    )
    {
        return fail(
            scenario,
            "an exactly attainable maximum inset did not collapse to the incenter"
        );
    }

    const stability::StandingTarget unattainableTriangleInset =
        stability::makeStandingTarget(
            triangle,
            pointAboveSupport(triangle, triangleQuery, 1.0),
            std::atan(1.5)
        );

    if (
        unattainableTriangleInset.requestedMarginAttainable() ||
        !approximatelyEqual(
            unattainableTriangleInset.requestedInset(),
            1.5
        ) ||
        !approximatelyEqual(
            unattainableTriangleInset.maximumInset(),
            1.0
        ) ||
        unattainableTriangleInset.safeRegion2D().rows() != 1 ||
        !approximatelyEqual(
            unattainableTriangleInset.target2D()(0),
            triangleIncenter(0)
        ) ||
        !approximatelyEqual(
            unattainableTriangleInset.target2D()(1),
            triangleIncenter(1)
        )
    )
    {
        return fail(
            scenario,
            "an impossible requested margin was falsely reported attainable"
        );
    }

    Eigen::MatrixXd rectangleVertices(4, 3);
    rectangleVertices <<
        -3.0, 0.0, -2.0,
         3.0, 0.0, -2.0,
         3.0, 0.0,  2.0,
        -3.0, 0.0,  2.0;

    const stability::SupportRegion rectangle =
        stability::SupportRegion::fromMesh(
            rectangleVertices,
            gravity,
            0.0
        );

    const stability::StandingTarget maximumRectangleInset =
        stability::makeStandingTarget(
            rectangle,
            pointAboveSupport(
                rectangle,
                Eigen::Vector2d::Zero(),
                2.0
            ),
            0.25 * std::acos(-1.0)
        );

    const std::array<Eigen::Vector2d, 2> expectedMaximumRectangle =
    {{
        rectangle.toPlaneCoordinates(Eigen::Vector3d(-1.0, 0.0, 0.0)),
        rectangle.toPlaneCoordinates(Eigen::Vector3d( 1.0, 0.0, 0.0))
    }};

    if (
        !approximatelyEqual(rectangle.maximumInset(), 2.0) ||
        !maximumRectangleInset.requestedMarginAttainable() ||
        maximumRectangleInset.safeRegion2D().rows() != 2 ||
        !equivalentPointSets(
            maximumRectangleInset.safeRegion2D(),
            pointsToMatrix(expectedMaximumRectangle)
        )
    )
    {
        return fail(
            scenario,
            "an attainable maximum inset did not preserve its line locus"
        );
    }

    const Eigen::Vector3d translation(3.0, 5.0, -4.0);
    Eigen::MatrixXd translatedVertices = makeSquareSupportMesh();
    translatedVertices.rowwise() += translation.transpose();

    const stability::SupportRegion translatedSquare =
        stability::SupportRegion::fromMesh(
            translatedVertices,
            gravity,
            0.25
        );

    const Eigen::Vector3d originalCenterOfMass =
        pointAboveSupport(square, Eigen::Vector2d(3.0, 4.0), 2.0);

    const stability::StandingTarget originalTarget =
        stability::makeStandingTarget(
            square,
            originalCenterOfMass,
            squareGamma
        );

    const stability::StandingTarget translatedTarget =
        stability::makeStandingTarget(
            translatedSquare,
            originalCenterOfMass + translation,
            squareGamma
        );

    if (
        !approximatelyEqual(
            translatedTarget.objective().target(),
            originalTarget.objective().target() + translation
        ) ||
        translatedTarget.state() != originalTarget.state()
    )
    {
        return fail(
            scenario,
            "standing-target construction is not translation covariant"
        );
    }

    const Eigen::Matrix3d rotation =
        Eigen::AngleAxisd(
            0.71,
            Eigen::Vector3d(1.0, 2.0, -1.0).normalized()
        ).toRotationMatrix();

    const Eigen::MatrixXd squareVertices = makeSquareSupportMesh();
    Eigen::MatrixXd rotatedVertices(squareVertices.rows(), 3);

    for (Eigen::Index row = 0; row < squareVertices.rows(); ++row)
    {
        rotatedVertices.row(row) =
            (
                rotation * squareVertices.row(row).transpose()
            ).transpose();
    }

    const stability::SupportRegion rotatedSquare =
        stability::SupportRegion::fromMesh(
            rotatedVertices,
            rotation * gravity,
            0.25
        );

    const stability::StandingTarget rotatedTarget =
        stability::makeStandingTarget(
            rotatedSquare,
            rotation * originalCenterOfMass,
            squareGamma
        );

    if (
        !approximatelyEqual(
            rotatedTarget.objective().target(),
            rotation * originalTarget.objective().target()
        ) ||
        rotatedTarget.state() != originalTarget.state()
    )
    {
        return fail(
            scenario,
            "standing-target construction is not rotation covariant"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " inside=edge=identity triangle-inradius=1"
        << '\n';

    return EXIT_SUCCESS;
}


int runSuspensionScenario()
{
    const std::string scenario = "suspension";

    const Eigen::Vector3d suspensionPoint(1.0, 2.0, 10.0);
    const Eigen::Vector3d gravity(0.0, 0.0, -1.0);

    const BalanceObjective objective =
        makeSuspensionObjective(
            suspensionPoint,
            gravity
        );

    if (
        !approximatelyEqual(
            objective.target(),
            suspensionPoint
        ) ||
        !approximatelyEqual(
            objective.gravity(),
            gravity
        )
    )
    {
        return fail(
            scenario,
            "the selected suspension point was not preserved as the target"
        );
    }

    BalanceObjectives objectives;
    objectives.push_back(objective);

    const Eigen::Vector3d offAxisCenterOfMass(4.0, 6.0, 0.0);

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassDistance(
                offAxisCenterOfMass,
                objectives
            ),
            5.0
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                offAxisCenterOfMass,
                objectives
            ),
            12.5
        )
    )
    {
        return fail(
            scenario,
            "the analytic 3-4-5 suspension objective is incorrect"
        );
    }

    const Eigen::Vector3d gravityAlignedCenterOfMass(1.0, 2.0, 3.0);

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassDistance(
                gravityAlignedCenterOfMass,
                objectives
            ),
            0.0
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                gravityAlignedCenterOfMass,
                objectives
            ),
            0.0
        )
    )
    {
        return fail(
            scenario,
            "a center of mass on the suspension gravity line is not balanced"
        );
    }

    BalanceObjectives scaledObjectives;
    scaledObjectives.push_back(
        makeSuspensionObjective(
            suspensionPoint,
            Eigen::Vector3d(0.0, 0.0, -1e300)
        )
    );

    if (
        !approximatelyEqual(
            scaledObjectives.front().gravity(),
            gravity
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                offAxisCenterOfMass,
                scaledObjectives
            ),
            12.5
        )
    )
    {
        return fail(
            scenario,
            "suspension energy depends on gravity-vector magnitude"
        );
    }

    const double infinity =
        std::numeric_limits<double>::infinity();
    const double notANumber =
        std::numeric_limits<double>::quiet_NaN();

    const std::array<Eigen::Vector3d, 3> invalidSuspensionPoints =
    {{
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidPoint : invalidSuspensionPoints)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        makeSuspensionObjective(
                            invalidPoint,
                            gravity
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a non-finite suspension point was accepted"
            );
        }
    }

    const std::array<Eigen::Vector3d, 4> invalidGravities =
    {{
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidGravity : invalidGravities)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        makeSuspensionObjective(
                            suspensionPoint,
                            invalidGravity
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a zero or non-finite suspension gravity was accepted"
            );
        }
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " distance=5 energy=12.5"
        << '\n';

    return EXIT_SUCCESS;
}


int runMultiScenario()
{
    const std::string scenario = "multi";

    const Eigen::Vector3d centerOfMass(2.0, 3.0, 4.0);

    const BalanceObjective first(
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );

    const BalanceObjective second(
        Eigen::Vector3d(1.0, 5.0, 1.0),
        Eigen::Vector3d(0.0, -2.0, 0.0)
    );

    const BalanceObjective third(
        Eigen::Vector3d(2.0, 3.0, 10.0),
        Eigen::Vector3d(1.0, 1.0, 0.0)
    );

    const double firstEnergy =
        projectedCenterOfMassEnergy(
            centerOfMass,
            first.target(),
            first.gravity()
        );

    const double secondEnergy =
        projectedCenterOfMassEnergy(
            centerOfMass,
            second.target(),
            second.gravity()
        );

    const double thirdEnergy =
        projectedCenterOfMassEnergy(
            centerOfMass,
            third.target(),
            third.gravity()
        );

    if (
        !approximatelyEqual(firstEnergy, 6.5) ||
        !approximatelyEqual(secondEnergy, 5.0) ||
        !approximatelyEqual(thirdEnergy, 18.0)
    )
    {
        return fail(
            scenario,
            "one of the three analytic objective terms is incorrect"
        );
    }

    const BalanceObjectives oneObjective =
    {
        first
    };

    const BalanceObjectives twoObjectives =
    {
        first,
        second
    };

    const BalanceObjectives threeObjectives =
    {
        first,
        second,
        third
    };

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                oneObjective
            ),
            firstEnergy
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                twoObjectives
            ),
            11.5
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                threeObjectives
            ),
            29.5
        ) ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassDistance(
                centerOfMass,
                threeObjectives
            ),
            std::sqrt(59.0)
        )
    )
    {
        return fail(
            scenario,
            "paper Eq. (5) is not the sum of all projected objectives"
        );
    }

    const BalanceObjectives permutedObjectives =
    {
        third,
        first,
        second
    };

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                permutedObjectives
            ),
            29.5
        )
    )
    {
        return fail(
            scenario,
            "aggregate energy depends on objective order"
        );
    }

    const BalanceObjectives duplicatedObjectives =
    {
        first,
        first
    };

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                duplicatedObjectives
            ),
            13.0
        )
    )
    {
        return fail(
            scenario,
            "duplicate objectives were not counted independently"
        );
    }

    const BalanceObjectives scaledGravityObjectives =
    {
        BalanceObjective(
            first.target(),
            Eigen::Vector3d(0.0, 0.0, -1e-300)
        ),
        BalanceObjective(
            second.target(),
            Eigen::Vector3d(0.0, -1e300, 0.0)
        ),
        BalanceObjective(
            third.target(),
            Eigen::Vector3d(1e250, 1e250, 0.0)
        )
    };

    if (
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                scaledGravityObjectives
            ),
            29.5
        )
    )
    {
        return fail(
            scenario,
            "Eq. (5) depends on individual gravity-vector magnitudes"
        );
    }

    /*
     * A fourth target is deliberate: the paper permits any number of bases,
     * while the inspected beta implementation stores only two objective
     * values. The local abstraction must not inherit that fixed-size limit.
     */
    BalanceObjectives fourObjectives = threeObjectives;
    fourObjectives.push_back(first);

    if (
        fourObjectives.size() != 4 ||
        !approximatelyEqual(
            aggregateProjectedCenterOfMassEnergy(
                centerOfMass,
                fourObjectives
            ),
            36.0
        )
    )
    {
        return fail(
            scenario,
            "the objective collection inherited the beta two-target limit"
        );
    }

    const BalanceObjectives emptyObjectives;

    if (
        !rejectsInvalidArgument(
            [&]()
            {
                static_cast<void>(
                    aggregateProjectedCenterOfMassEnergy(
                        centerOfMass,
                        emptyObjectives
                    )
                );
            }
        ) ||
        !rejectsInvalidArgument(
            [&]()
            {
                static_cast<void>(
                    aggregateProjectedCenterOfMassDistance(
                        centerOfMass,
                        emptyObjectives
                    )
                );
            }
        )
    )
    {
        return fail(
            scenario,
            "an empty high-level objective set was accepted"
        );
    }

    const double infinity =
        std::numeric_limits<double>::infinity();
    const double notANumber =
        std::numeric_limits<double>::quiet_NaN();

    const std::array<Eigen::Vector3d, 3> invalidCentersOfMass =
    {{
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidCenterOfMass : invalidCentersOfMass)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        aggregateProjectedCenterOfMassEnergy(
                            invalidCenterOfMass,
                            threeObjectives
                        )
                    );
                }
            ) ||
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        aggregateProjectedCenterOfMassDistance(
                            invalidCenterOfMass,
                            threeObjectives
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a non-finite center of mass was accepted"
            );
        }
    }

    const std::array<Eigen::Vector3d, 3> invalidTargets =
    {{
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidTarget : invalidTargets)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        BalanceObjective(
                            invalidTarget,
                            Eigen::Vector3d(0.0, 0.0, -1.0)
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a non-finite target was accepted"
            );
        }
    }

    const std::array<Eigen::Vector3d, 4> invalidGravities =
    {{
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(infinity, 0.0, 0.0),
        Eigen::Vector3d(0.0, -infinity, 0.0),
        Eigen::Vector3d(0.0, 0.0, notANumber)
    }};

    for (const Eigen::Vector3d& invalidGravity : invalidGravities)
    {
        if (
            !rejectsInvalidArgument(
                [&]()
                {
                    static_cast<void>(
                        BalanceObjective(
                            Eigen::Vector3d::Zero(),
                            invalidGravity
                        )
                    );
                }
            )
        )
        {
            return fail(
                scenario,
                "a zero or non-finite objective gravity was accepted"
            );
        }
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " objectives=3 energy=29.5 distance=sqrt(59)"
        << '\n';

    return EXIT_SUCCESS;
}

} // namespace


int main(
    int argc,
    char** argv
)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: stability-regression "
            << "<support|safe|standing-target|suspension|multi>\n";
        return EXIT_FAILURE;
    }

    const std::string scenario(argv[1]);

    if (scenario == "support")
    {
        return runSupportScenario();
    }

    if (scenario == "safe")
    {
        return runSafeScenario();
    }

    if (scenario == "standing-target")
    {
        return runStandingTargetScenario();
    }

    if (scenario == "suspension")
    {
        return runSuspensionScenario();
    }

    if (scenario == "multi")
    {
        return runMultiScenario();
    }

    return fail(
        scenario,
        "unknown scenario"
    );
}
