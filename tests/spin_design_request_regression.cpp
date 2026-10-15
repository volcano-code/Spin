#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

spin::TriangleMesh cube(double minimum, double maximum)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices <<
        minimum, minimum, minimum,
        maximum, minimum, minimum,
        maximum, maximum, minimum,
        minimum, maximum, minimum,
        minimum, minimum, maximum,
        maximum, minimum, maximum,
        maximum, maximum, maximum,
        minimum, maximum, maximum;
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

bool throwsInvalid(const std::function<void()>& operation)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    return false;
}

int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-design-request " << message << '\n';
    return EXIT_FAILURE;
}

} // namespace

int main()
{
    spin::SpinDesignRequest request;
    request.outerSurface = cube(0.0, 1.0);
    request.designDomain = cube(0.2, 0.8);
    request.mode = spin::SpinMode::StaticBalance;
    request.materialDensity = 1.25;
    request.originalUnitsPerNormalizedUnit = 2.0;
    request.minimumWallThickness = 0.1;
    request.settings.maximumIterations = 60;
    request.settings.maximumRestarts = 8;
    request.settings.maximumFunctionEvaluationsPerRestart = 60;
    request.settings.equationTolerance = 1e-12;
    request.settings.physicalRelativeTolerance = 1e-6;

    spin::validateSpinDesignRequest(request);
    const std::string hash = spin::computeSpinDesignRequestHash(request);
    if (hash.size() != 64 || hash != spin::computeSpinDesignRequestHash(request))
    {
        return fail("request hash is absent or non-deterministic");
    }
    if (
        std::string(spin::spinModeToken(spin::SpinMode::StaticBalance)) !=
            "STATIC_BALANCE" ||
        std::string(spin::spinModeToken(spin::SpinMode::FidgetSpinner)) !=
            "FIDGET_SPINNER" ||
        std::string(spin::spinModeToken(spin::SpinMode::YoYo)) != "YO_YO" ||
        std::string(spin::spinModeToken(spin::SpinMode::SpinningTop)) !=
            "SPINNING_TOP" ||
        std::string(spin::spinModeToken(spin::SpinMode::BuoyantStability)) !=
            "BUOYANT_STABILITY"
    )
    {
        return fail("mode tokens do not match the public contract");
    }

    auto endpoint=request;endpoint.mode=spin::SpinMode::SpinningTop;
    const auto topHash=spin::computeSpinDesignRequestHash(endpoint);
    endpoint.spinningTopObjective=spin::SpinningTopObjective::LowerCenterOfMass;
    spin::validateSpinDesignRequest(endpoint);
    if(spin::computeSpinDesignRequestHash(endpoint)==topHash)
        return fail("pure-height endpoint aliases finite top request identity");
    endpoint.centerOfMassWeight=1.0;
    if(!throwsInvalid([&]{spin::validateSpinDesignRequest(endpoint);}))return fail("pure endpoint accepted extra weight");
    endpoint.centerOfMassWeight=0.0;endpoint.mode=spin::SpinMode::YoYo;
    if(!throwsInvalid([&]{spin::validateSpinDesignRequest(endpoint);}))return fail("yo-yo accepted top endpoint selection");
    endpoint.mode=spin::SpinMode::SpinningTop;
    endpoint.spinningTopObjective=static_cast<spin::SpinningTopObjective>(99);
    if(!throwsInvalid([&]{spin::validateSpinDesignRequest(endpoint);}))return fail("invalid top objective accepted");

    spin::SpinDesignRequest changed = request;
    changed.materialDensity = 1.5;
    if (spin::computeSpinDesignRequestHash(changed) == hash)
    {
        return fail("density is not bound into the request hash");
    }

    spin::SpinDesignRequest pareto = request;
    pareto.staticObjective = spin::StaticObjective::ParetoTradeoff;
    pareto.staticParetoWeight = 0.25;
    spin::validateSpinDesignRequest(pareto);
    if (spin::computeSpinDesignRequestHash(pareto) == hash)
    {
        return fail("static Pareto weight is not bound into the request hash");
    }
    pareto.staticParetoWeight = 1.1;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(pareto); }))
    {
        return fail("out-of-range static Pareto weight was accepted");
    }

    spin::SpinDesignRequest invalid = request;
    invalid.targetFrame.axes(0, 0) = 2.0;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalid); }))
    {
        return fail("non-orthonormal target frame was accepted");
    }
    invalid = request;
    invalid.targetFrame.axes.col(2) *= -1.0;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalid); }))
    {
        return fail("left-handed target frame was accepted");
    }
    invalid = request;
    invalid.materialDensity = 0.0;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalid); }))
    {
        return fail("non-positive material density was accepted");
    }
    invalid = request;
    invalid.designDomain = cube(2.0, 3.0);
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalid); }))
    {
        return fail("design domain outside the outer surface was accepted");
    }
    invalid = request;
    invalid.settings.equationTolerance = std::nan("");
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalid); }))
    {
        return fail("non-finite solver tolerance was accepted");
    }

    spin::SpinDesignRequest buoyant = request;
    buoyant.mode = spin::SpinMode::BuoyantStability;
    MassMoments fixedShell;
    fixedShell.mass = 0.25;
    fixedShell.firstMoment = Eigen::Vector3d(0.0, 0.0, 0.05);
    fixedShell.secondMoment = 0.1 * Eigen::Matrix3d::Identity();
    buoyant.fixedShellMoments =
        physics::RawMassMoments10::fromMassMoments(fixedShell);
    buoyant.materialDensity = 0.5;
    buoyant.buoyancy.liquidDensity = 1.0;
    buoyant.buoyancy.displacedLiquidVolume = 1.5;
    buoyant.buoyancy.centerOfBuoyancyTarget = Eigen::Vector3d(0.1, -0.2, 0.0);
    buoyant.buoyancy.waterlineArea = 2.0;
    buoyant.buoyancy.waterlineCentralSecondMoment =
        (Eigen::Matrix2d() << 0.4, 0.05, 0.05, 0.3).finished();
    buoyant.buoyancy.fixedShellVolume = 0.5;
    spin::validateSpinDesignRequest(buoyant);
    const std::string buoyantHash = spin::computeSpinDesignRequestHash(buoyant);
    spin::SpinDesignRequest changedBuoyancy = buoyant;
    changedBuoyancy.buoyancy.waterlineArea = 2.5;
    if (spin::computeSpinDesignRequestHash(changedBuoyancy) == buoyantHash)
    {
        return fail("buoyancy inputs are not bound into the request hash");
    }
    spin::SpinDesignRequest invalidBuoyancy = buoyant;
    invalidBuoyancy.buoyancy.liquidDensity = 0.0;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalidBuoyancy); }))
    {
        return fail("non-positive liquid density was accepted");
    }
    invalidBuoyancy = buoyant;
    invalidBuoyancy.buoyancy.waterlineCentralSecondMoment(0, 1) = 0.3;
    if (!throwsInvalid([&] { spin::validateSpinDesignRequest(invalidBuoyancy); }))
    {
        return fail("non-symmetric waterline moment tensor was accepted");
    }

    spin::QuadricHollowingResult result;
    result.status = spin::SpinSolverStatus::OptimizationDidNotConverge;
    if (
        std::string(spin::spinSolverStatusToken(result.status)) !=
        "OPTIMIZATION_DID_NOT_CONVERGE"
    )
    {
        return fail("terminal status token is unstable");
    }

    std::cout << "[PASS] scenario=spin-design-request hash=" << hash << '\n';
    return EXIT_SUCCESS;
}
