#include "spin/QuadricHollowingSolver.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;


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


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-faster-end-to-end "
              << message << '\n';
    return EXIT_FAILURE;
}


spin::SpinDesignRequest baseRequest(spin::SpinMode mode)
{
    spin::SpinDesignRequest request;
    request.outerSurface = cube(-5.0, 5.0);
    request.designDomain = cube(-1.0, 1.0);
    request.minimumWallThickness = 4.0;
    request.mode = mode;
    request.settings.maximumRestarts = 1;
    return request;
}


physics::RawMassMoments10 rawMoments(
    double mass,
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& diagonalSecond
)
{
    MassMoments moments;
    moments.mass = mass;
    moments.firstMoment = first;
    moments.secondMoment = diagonalSecond.asDiagonal();
    return physics::RawMassMoments10::fromMassMoments(moments);
}


spin::SpinDesignRequest requestWithTarget(
    spin::SpinDesignRequest request,
    const Quadric& quadric,
    const physics::RawMassMoments10& target
)
{
    const spin::QuadricMomentEvaluation variable =
        spin::evaluateQuadricMaterialMoments(request, quadric);
    request.fixedShellMoments = target - variable.variableMomentsTarget;
    spin::validateSpinDesignRequest(request);
    return request;
}


bool solveAndCheck(
    const char* name,
    const spin::SpinDesignRequest& request,
    const Quadric& expected,
    spin::QuadricSurfaceType expectedType,
    std::string* outputHash
)
{
    const spin::QuadricMomentEvaluation expectedEvaluation =
        spin::evaluateQuadricMaterialMoments(request, expected);
    if (expectedEvaluation.quadricTriangulation.type != expectedType)
    {
        std::cerr << "unexpected interface type for " << name << '\n';
        return false;
    }
    const spin::QuadricHollowingResult result = spin::solveQuadricHollowing(
        request,
        spin::QuadricMomentEvaluatorSettings(),
        {expected}
    );
    const spin::QuadricHollowingResult repeated = spin::solveQuadricHollowing(
        request,
        spin::QuadricMomentEvaluatorSettings(),
        {expected}
    );
    if (
        result.status != spin::SpinSolverStatus::Converged ||
        result.normalizedEquationResidual > 1e-12 ||
        result.physicalRelativeResidual > 1e-6 ||
        result.materialBoundary.faces.rows() == 0 ||
        result.cavityBoundary.faces.rows() == 0 ||
        result.geometryHash.size() != 64 || result.traceHash.size() != 64 ||
        result.outputHash.size() != 64 ||
        result.outputHash != repeated.outputHash ||
        !result.quadricCoefficients.isApprox(expected.values(), 1e-13) ||
        result.outerSurface.vertices != request.outerSurface.vertices ||
        result.outerSurface.faces != request.outerSurface.faces
    )
    {
        std::cerr << name << " status="
                  << spin::spinSolverStatusToken(result.status)
                  << " residual=" << result.normalizedEquationResidual
                  << " physical=" << result.physicalRelativeResidual
                  << " faces=" << result.materialBoundary.faces.rows()
                  << " cavity_faces=" << result.cavityBoundary.faces.rows()
                  << " geometry_hash=" << result.geometryHash.size()
                  << " output_hash=" << result.outputHash.size()
                  << " repeat_equal="
                  << (result.outputHash == repeated.outputHash)
                  << " coefficient_error="
                  << (result.quadricCoefficients - expected.values()).norm()
                  << " outer_v_equal="
                  << (result.outerSurface.vertices == request.outerSurface.vertices)
                  << " outer_f_equal="
                  << (result.outerSurface.faces == request.outerSurface.faces)
                  << " diagnostic=" << result.diagnostic << '\n';
        return false;
    }
    if (request.mode == spin::SpinMode::BuoyantStability)
    {
        if (
            !std::isfinite(result.buoyantMetacenterHeight) ||
            !std::isfinite(result.buoyantStabilityMargin) ||
            result.buoyantStabilityMargin <= 0.0
        )
        {
            std::cerr << name << " invalid buoyancy audit: metacenter="
                      << result.buoyantMetacenterHeight << " margin="
                      << result.buoyantStabilityMargin << '\n';
            return false;
        }
    }
    *outputHash = result.outputHash;
    return true;
}


} // namespace


int main()
{
    std::vector<std::string> hashes;

    Vector10 values = Vector10::Zero();
    values(0) = 0.371;
    values(1) = 0.8;
    values(2) = 0.123;
    const Quadric staticPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest staticRequest =
        baseRequest(spin::SpinMode::StaticBalance);
    const spin::QuadricMomentEvaluation staticVariable =
        spin::evaluateQuadricMaterialMoments(staticRequest, staticPlane);
    staticRequest.fixedShellMoments = rawMoments(
        10.0,
        -staticVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "static",
            staticRequest,
            staticPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("static planar fixture failed");
    }

    values.setZero();
    values(0) = 0.2;
    values(1) = 0.12;
    values(2) = -0.08;
    values(3) = -1.0;
    const Quadric staticLowPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest staticLowRequest =
        baseRequest(spin::SpinMode::StaticBalance);
    staticLowRequest.staticObjective =
        spin::StaticObjective::LowerCenterOfMass;
    const spin::QuadricMomentEvaluation staticLowVariable =
        spin::evaluateQuadricMaterialMoments(
            staticLowRequest,
            staticLowPlane
        );
    staticLowRequest.fixedShellMoments = rawMoments(
        10.0 - staticLowVariable.variableMomentsTarget.mass(),
        Eigen::Vector3d(0.0, 0.0, 2.0) -
            staticLowVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "static-lower-center",
            staticLowRequest,
            staticLowPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("static lower-center planar fixture failed");
    }

    values.setZero();
    values(0) = 0.2;
    values(1) = 0.2;
    values(2) = -0.1;
    values(3) = -1.0;
    const Quadric buoyantPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest buoyantRequest =
        baseRequest(spin::SpinMode::BuoyantStability);
    buoyantRequest.buoyancy.liquidDensity = 1.0;
    buoyantRequest.buoyancy.displacedLiquidVolume = 20.0;
    buoyantRequest.buoyancy.centerOfBuoyancyTarget =
        Eigen::Vector3d(0.02, -0.03, 0.0);
    buoyantRequest.buoyancy.waterlineArea = 4.0;
    buoyantRequest.buoyancy.waterlineCentralSecondMoment =
        2.0 * Eigen::Matrix2d::Identity();
    const spin::QuadricMomentEvaluation buoyantVariable =
        spin::evaluateQuadricMaterialMoments(
            buoyantRequest,
            buoyantPlane
        );
    const physics::RawMassMoments10 buoyantTarget = rawMoments(
        20.0,
        Eigen::Vector3d(0.4, -0.6, 0.0),
        Eigen::Vector3d::Constant(100.0)
    );
    buoyantRequest.fixedShellMoments =
        buoyantTarget - buoyantVariable.variableMomentsTarget;
    buoyantRequest.buoyancy.fixedShellVolume =
        buoyantRequest.fixedShellMoments.mass();
    hashes.emplace_back();
    if (!solveAndCheck(
            "buoyant",
            buoyantRequest,
            buoyantPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("buoyant planar fixture failed");
    }

    values.setZero();
    values(1) = 0.8;
    values(4) = 1.0;
    values(5) = 1.0;
    const Quadric fidgetCircle = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest fidgetRequest =
        baseRequest(spin::SpinMode::FidgetSpinner);
    const spin::QuadricMomentEvaluation fidgetVariable =
        spin::evaluateQuadricMaterialMoments(fidgetRequest, fidgetCircle);
    fidgetRequest.fixedShellMoments = rawMoments(
        10.0,
        -fidgetVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "fidget",
            fidgetRequest,
            fidgetCircle,
            spin::QuadricSurfaceType::EllipticCylinder,
            &hashes.back()
        ))
    {
        return fail("fidget circular-void fixture failed");
    }

    values.setZero();
    values(3) = 1.0;
    values(4) = 1.0;
    values(5) = 1.0;
    values(6) = -1.0;
    const Quadric yoyoHyperboloid = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest yoyoRequest = baseRequest(spin::SpinMode::YoYo);
    yoyoRequest = requestWithTarget(
        yoyoRequest,
        yoyoHyperboloid,
        rawMoments(
            20.0,
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(100.0, 100.0, 200.0)
        )
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "yo-yo",
            yoyoRequest,
            yoyoHyperboloid,
            spin::QuadricSurfaceType::Hyperboloid,
            &hashes.back()
        ))
    {
        return fail("yo-yo symmetric fixture failed");
    }

    values(3) = -1.0;
    const Quadric topHyperboloid = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest topRequest =
        baseRequest(spin::SpinMode::SpinningTop);
    topRequest.centerOfMassWeight = 0.6;
    topRequest = requestWithTarget(
        topRequest,
        topHyperboloid,
        rawMoments(
            20.0,
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(100.0, 100.0, 200.0)
        )
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "top",
            topRequest,
            topHyperboloid,
            spin::QuadricSurfaceType::Hyperboloid,
            &hashes.back()
        ))
    {
        return fail("spinning-top hyperboloid fixture failed");
    }

    std::cout << "[PASS] scenario=spin-it-faster-end-to-end";
    for (const std::string& hash : hashes)
    {
        std::cout << ' ' << hash;
    }
    std::cout << '\n';
    return EXIT_SUCCESS;
}
