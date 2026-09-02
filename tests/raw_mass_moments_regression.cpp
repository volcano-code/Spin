#include "physics/MassProperties.h"
#include "physics/RawMassMoments10.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>


namespace
{


bool approximatelyEqual(double actual, double expected)
{
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return std::abs(actual - expected) <= 1e-13 + 1e-13 * scale;
}


MassMoments makeKnownMoments()
{
    MassMoments moments;
    moments.mass = 2.0;
    moments.firstMoment = Eigen::Vector3d(3.0, 5.0, 7.0);
    moments.secondMoment <<
        11.0, 19.0, 23.0,
        19.0, 13.0, 29.0,
        23.0, 29.0, 17.0;
    return moments;
}


int fail(const std::string& scenario, const std::string& message)
{
    std::cerr << "[FAIL] scenario=" << scenario << ' ' << message << '\n';
    return EXIT_FAILURE;
}


int runOrdering()
{
    const physics::RawMassMoments10 raw =
        physics::RawMassMoments10::fromMassMoments(makeKnownMoments());
    Eigen::Matrix<double, 10, 1> expected;
    expected << 2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 29.0, 23.0, 19.0;
    if (!raw.values().isApprox(expected, 1e-14))
    {
        return fail("ordering", "paper component ordering is incorrect");
    }
    return EXIT_SUCCESS;
}


int runConversion()
{
    const MassMoments source = makeKnownMoments();
    const physics::RawMassMoments10 raw =
        physics::RawMassMoments10::fromMassMoments(source);
    const MassMoments roundTrip = raw.toMassMoments();
    if (
        !approximatelyEqual(roundTrip.mass, source.mass) ||
        !roundTrip.firstMoment.isApprox(source.firstMoment, 1e-14) ||
        !roundTrip.secondMoment.isApprox(source.secondMoment, 1e-14)
    )
    {
        return fail("conversion", "MassMoments round trip changed a component");
    }

    const Eigen::Vector3d expectedCenter(1.5, 2.5, 3.5);
    if (!raw.centerOfMass().isApprox(expectedCenter, 1e-14))
    {
        return fail("conversion", "center of mass is inconsistent");
    }

    const Eigen::Matrix3d expectedOrigin =
        source.secondMoment.trace() * Eigen::Matrix3d::Identity() -
        source.secondMoment;
    if (!raw.inertiaAboutOrigin().isApprox(expectedOrigin, 1e-14))
    {
        return fail("conversion", "origin inertia is inconsistent");
    }
    if (!raw.inertiaAboutCenter().isApprox(
            source.inertiaAboutCenter(),
            1e-14
        ))
    {
        return fail("conversion", "center inertia is inconsistent");
    }
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "usage: raw-mass-moments-regression SCENARIO\n";
        return EXIT_FAILURE;
    }
    const std::string scenario(argv[1]);
    try
    {
        if (scenario == "ordering") return runOrdering();
        if (scenario == "conversion") return runConversion();
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
    return fail(scenario, "unknown scenario");
}
