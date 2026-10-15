#include "spin/PlanarMassMoments.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

using Loop = spin::PlanarBoundaryLoop;

bool close(double actual, double expected)
{
    return std::abs(actual - expected) <=
        1e-12 * std::max(1.0, std::abs(expected));
}

int fail(const char* message)
{
    std::cerr << message << '\n';
    return 1;
}

Loop rectangle(double x0, double y0, double x1, double y1, bool counterClockwise)
{
    Loop points{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    if (!counterClockwise)
        std::reverse(points.begin(), points.end());
    return points;
}

} // namespace

int main()
{
    const auto rectangleMoments = spin::integratePlanarBoundaryLoops({
        rectangle(0.0, 0.0, 2.0, 3.0, true)
    });
    if (!close(rectangleMoments.area, 6.0) ||
        !close(rectangleMoments.firstX, 6.0) ||
        !close(rectangleMoments.firstY, 9.0) ||
        !close(rectangleMoments.secondXX, 8.0) ||
        !close(rectangleMoments.secondYY, 18.0) ||
        !close(rectangleMoments.secondXY, 9.0))
        return fail("translated-independent rectangle moments are incorrect");
    const auto rectangleVector = rectangleMoments.weighted(1.0).vector();
    const double expectedVector[6] = {6.0, 6.0, 9.0, 8.0, 18.0, 9.0};
    for (Eigen::Index i = 0; i < 6; ++i)
        if (!close(rectangleVector[i], expectedVector[i]))
            return fail("planar moment vector ordering is incorrect");

    const auto annulusMoments = spin::integratePlanarBoundaryLoops({
        rectangle(0.0, 0.0, 4.0, 3.0, true),
        rectangle(1.0, 1.0, 3.0, 2.0, false)
    });
    if (!close(annulusMoments.area, 10.0) ||
        !close(annulusMoments.firstX, 20.0) ||
        !close(annulusMoments.firstY, 15.0) ||
        !close(annulusMoments.secondXX, 166.0 / 3.0) ||
        !close(annulusMoments.secondYY, 94.0 / 3.0) ||
        !close(annulusMoments.secondXY, 30.0))
        return fail("oppositely oriented hole moments are incorrect");

    const auto translated = spin::integratePlanarBoundaryLoops({
        rectangle(10.0, -4.0, 12.0, -1.0, true)
    });
    if (!close(translated.area, 6.0) ||
        !close(translated.firstX, 66.0) ||
        !close(translated.firstY, -15.0) ||
        !close(translated.secondXX, 728.0) ||
        !close(translated.secondYY, 42.0) ||
        !close(translated.secondXY, -165.0))
        return fail("global translation of all six moments is incorrect");

    const auto densityMoments = rectangleMoments.weighted(2.5);
    if (!close(densityMoments.mass, 15.0) ||
        !close(densityMoments.firstX, 15.0) ||
        !close(densityMoments.secondXY, 22.5))
        return fail("areal-density scaling did not scale every raw moment");

    try
    {
        (void)rectangleMoments.weighted(std::numeric_limits<double>::quiet_NaN());
        return fail("non-finite areal density was accepted");
    }
    catch (const std::invalid_argument&)
    {
    }

    try
    {
        (void)rectangleMoments.weighted(std::numeric_limits<double>::max());
        return fail("overflowed density-weighted moments were returned");
    }
    catch (const std::overflow_error&)
    {
    }

    try
    {
        const auto small = spin::integratePlanarBoundaryLoops({
            rectangle(0.0, 0.0, 1e-5, 1e-5, true)
        });
        (void)small.weighted(std::numeric_limits<double>::denorm_min());
        return fail("underflowed positive planar area was returned as zero");
    }
    catch (const std::overflow_error&)
    {
    }

    try
    {
        (void)spin::integratePlanarBoundaryLoops({rectangle(0.0, 0.0, 2.0, 3.0, false)});
        return fail("clockwise outer boundary was accepted as positive material");
    }
    catch (const std::invalid_argument&)
    {
    }

    try
    {
        (void)spin::integratePlanarBoundaryLoops({Loop{{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}}});
        return fail("zero-area boundary loop was accepted");
    }
    catch (const std::invalid_argument&)
    {
    }

    try
    {
        (void)spin::integratePlanarBoundaryLoops({Loop{{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}});
        return fail("consecutive duplicate boundary vertices were accepted");
    }
    catch (const std::invalid_argument&)
    {
    }

    return 0;
}
