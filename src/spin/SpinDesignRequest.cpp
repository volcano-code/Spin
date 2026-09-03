#include "spin/SpinDesignRequest.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <type_traits>


namespace
{


constexpr double frameTolerance = 1e-12;


void appendScalar(std::ostringstream& stream, double value)
{
    stream << std::hexfloat << value << ';';
}


template <typename Derived>
void appendMatrix(
    std::ostringstream& stream,
    const Eigen::MatrixBase<Derived>& matrix
)
{
    stream << matrix.rows() << 'x' << matrix.cols() << ':';
    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column)
        {
            if constexpr (
                std::is_floating_point_v<typename Derived::Scalar>
            )
            {
                appendScalar(stream, matrix(row, column));
            }
            else
            {
                stream << matrix(row, column) << ';';
            }
        }
    }
}


void validateMesh(
    const spin::TriangleMesh& mesh,
    const char* name,
    manufacturing::BoundaryOrientation orientation
)
{
    const manufacturing::BoundaryAudit audit =
        manufacturing::auditTriangleBoundary(
            mesh.vertices,
            mesh.faces,
            orientation,
            false
        );
    if (!audit.strictBoundaryValid())
    {
        throw std::invalid_argument(
            std::string(name) +
            " must be a finite, closed, consistently oriented, outward "
            "manifold triangle boundary without self-intersections"
        );
    }
}


void requireFinite(double value, const char* name)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}


} // namespace


namespace spin
{


const char* spinModeToken(SpinMode mode) noexcept
{
    switch (mode)
    {
    case SpinMode::StaticBalance:
        return "STATIC_BALANCE";
    case SpinMode::BuoyantStability:
        return "BUOYANT_STABILITY";
    case SpinMode::FidgetSpinner:
        return "FIDGET_SPINNER";
    case SpinMode::YoYo:
        return "YO_YO";
    case SpinMode::SpinningTop:
        return "SPINNING_TOP";
    }
    return "UNKNOWN";
}


const char* staticObjectiveToken(StaticObjective objective) noexcept
{
    switch (objective)
    {
    case StaticObjective::MaximumMass:
        return "MAXIMUM_MASS";
    case StaticObjective::LowerCenterOfMass:
        return "LOWER_CENTER_OF_MASS";
    }
    return "UNKNOWN";
}


const char* spinSolverStatusToken(SpinSolverStatus status) noexcept
{
    switch (status)
    {
    case SpinSolverStatus::Converged:
        return "CONVERGED";
    case SpinSolverStatus::InvalidRequest:
        return "INVALID_REQUEST";
    case SpinSolverStatus::OptimizationDidNotConverge:
        return "OPTIMIZATION_DID_NOT_CONVERGE";
    case SpinSolverStatus::LineSearchExhausted:
        return "LINE_SEARCH_EXHAUSTED";
    case SpinSolverStatus::GeometryEvaluationFailed:
        return "GEOMETRY_EVALUATION_FAILED";
    case SpinSolverStatus::InvalidMaterialArtifact:
        return "INVALID_MATERIAL_ARTIFACT";
    case SpinSolverStatus::ResourceLimitReached:
        return "RESOURCE_LIMIT_REACHED";
    case SpinSolverStatus::Indeterminate:
        return "INDETERMINATE";
    }
    return "UNKNOWN";
}


void validateSpinDesignRequest(const SpinDesignRequest& request)
{
    validateMesh(
        request.outerSurface,
        "outer surface",
        manufacturing::BoundaryOrientation::OutwardPositive
    );
    validateMesh(
        request.designDomain,
        "design domain",
        manufacturing::BoundaryOrientation::OutwardPositive
    );

    if (
        !request.targetFrame.origin.allFinite() ||
        !request.targetFrame.axes.allFinite()
    )
    {
        throw std::invalid_argument("target frame must be finite");
    }
    const Eigen::Matrix3d orthogonality =
        request.targetFrame.axes.transpose() * request.targetFrame.axes;
    if (
        !orthogonality.isApprox(Eigen::Matrix3d::Identity(), frameTolerance) ||
        std::abs(request.targetFrame.axes.determinant() - 1.0) > frameTolerance
    )
    {
        throw std::invalid_argument(
            "target frame axes must be orthonormal and right-handed"
        );
    }

    requireFinite(request.materialDensity, "material density");
    requireFinite(
        request.originalUnitsPerNormalizedUnit,
        "original-unit scale"
    );
    requireFinite(request.minimumWallThickness, "minimum wall thickness");
    requireFinite(request.centerOfMassWeight, "center-of-mass weight");
    if (request.materialDensity <= 0.0)
    {
        throw std::invalid_argument("material density must be positive");
    }
    if (request.originalUnitsPerNormalizedUnit <= 0.0)
    {
        throw std::invalid_argument("original-unit scale must be positive");
    }
    if (request.minimumWallThickness < 0.0)
    {
        throw std::invalid_argument(
            "minimum wall thickness must be nonnegative"
        );
    }
    if (!request.fixedShellMoments.isFinite() ||
        request.fixedShellMoments.mass() < 0.0)
    {
        throw std::invalid_argument(
            "fixed-shell moments must be finite with nonnegative mass"
        );
    }
    if (
        request.mode != SpinMode::SpinningTop &&
        request.centerOfMassWeight != 0.0
    )
    {
        throw std::invalid_argument(
            "center-of-mass weight is supported only for spinning tops"
        );
    }
    if (
        request.mode != SpinMode::StaticBalance &&
        request.staticObjective != StaticObjective::MaximumMass
    )
    {
        throw std::invalid_argument(
            "static objective is supported only in static-balance mode"
        );
    }

    const BuoyantStabilityData& buoyancy = request.buoyancy;
    if (request.mode == SpinMode::BuoyantStability)
    {
        requireFinite(buoyancy.liquidDensity, "liquid density");
        requireFinite(
            buoyancy.displacedLiquidVolume,
            "displaced liquid volume"
        );
        requireFinite(buoyancy.waterlineArea, "waterline area");
        requireFinite(buoyancy.fixedShellVolume, "fixed shell volume");
        if (
            !buoyancy.centerOfBuoyancyTarget.allFinite() ||
            !buoyancy.waterlineCentralSecondMoment.allFinite()
        )
        {
            throw std::invalid_argument(
                "buoyancy center and waterline moments must be finite"
            );
        }
        if (
            buoyancy.liquidDensity <= 0.0 ||
            buoyancy.displacedLiquidVolume <= 0.0 ||
            buoyancy.waterlineArea <= 0.0 ||
            buoyancy.fixedShellVolume < 0.0
        )
        {
            throw std::invalid_argument(
                "buoyancy density, displaced volume and waterline area must "
                "be positive; fixed shell volume must be nonnegative"
            );
        }
        if (!buoyancy.waterlineCentralSecondMoment.isApprox(
                buoyancy.waterlineCentralSecondMoment.transpose(),
                frameTolerance))
        {
            throw std::invalid_argument(
                "waterline central second moment must be symmetric"
            );
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> waterlineEigen(
            buoyancy.waterlineCentralSecondMoment
        );
        if (
            waterlineEigen.info() != Eigen::Success ||
            waterlineEigen.eigenvalues().minCoeff() <= 0.0
        )
        {
            throw std::invalid_argument(
                "waterline central second moment must be positive definite"
            );
        }
        const double expectedShellMass =
            request.materialDensity * buoyancy.fixedShellVolume;
        const double shellScale = std::max(1.0, std::abs(expectedShellMass));
        if (
            std::abs(
                request.fixedShellMoments.mass() - expectedShellMass
            ) > 1e-10 * shellScale
        )
        {
            throw std::invalid_argument(
                "fixed-shell mass is inconsistent with material density and "
                "fixed shell volume"
            );
        }
    }
    else if (
        buoyancy.liquidDensity != 0.0 ||
        buoyancy.displacedLiquidVolume != 0.0 ||
        !buoyancy.centerOfBuoyancyTarget.isZero(0.0) ||
        buoyancy.waterlineArea != 0.0 ||
        !buoyancy.waterlineCentralSecondMoment.isZero(0.0) ||
        buoyancy.fixedShellVolume != 0.0
    )
    {
        throw std::invalid_argument(
            "buoyancy inputs are supported only in buoyant-stability mode"
        );
    }

    const SpinSolverSettings& settings = request.settings;
    requireFinite(settings.equationTolerance, "equation tolerance");
    requireFinite(
        settings.physicalRelativeTolerance,
        "physical relative tolerance"
    );
    requireFinite(settings.lineSearchReduction, "line-search reduction");
    requireFinite(
        settings.maximumGeodesicStep,
        "maximum geodesic step"
    );
    if (
        settings.maximumIterations <= 0 || settings.maximumRestarts <= 0 ||
        settings.maximumFunctionEvaluationsPerRestart <= 0 ||
        settings.equationTolerance <= 0.0 ||
        settings.physicalRelativeTolerance <= 0.0 ||
        settings.lineSearchReduction != 0.5 ||
        settings.maximumGeodesicStep != 0.5
    )
    {
        throw std::invalid_argument(
            "solver limits/tolerances are invalid or paper constants changed"
        );
    }

    Eigen::MatrixXi cavityFaces = request.designDomain.faces;
    cavityFaces.col(0).swap(cavityFaces.col(2));
    const manufacturing::MaterialShellAudit shellAudit =
        manufacturing::auditMaterialShells(
            request.outerSurface.vertices,
            request.outerSurface.faces,
            request.designDomain.vertices,
            cavityFaces,
            request.minimumWallThickness
        );
    if (
        !shellAudit.wallClearance.completed ||
        !shellAudit.wallClearance.intersectionFree ||
        !shellAudit.wallClearance.meetsRequestedThickness ||
        !shellAudit.containment.completed ||
        !shellAudit.containment.allCavitiesInsideOuter
    )
    {
        throw std::invalid_argument(
            "design domain must be disjoint from, inside, and sufficiently "
            "clear of the outer surface"
        );
    }
}


std::string computeSpinDesignRequestHash(const SpinDesignRequest& request)
{
    validateSpinDesignRequest(request);
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "SpinDesignRequest-v2;";
    stream << spinModeToken(request.mode) << ';';
    stream << staticObjectiveToken(request.staticObjective) << ';';
    appendMatrix(stream, request.outerSurface.vertices);
    appendMatrix(stream, request.outerSurface.faces);
    appendMatrix(stream, request.designDomain.vertices);
    appendMatrix(stream, request.designDomain.faces);
    appendMatrix(stream, request.fixedShellMoments.values());
    appendMatrix(stream, request.targetFrame.origin);
    appendMatrix(stream, request.targetFrame.axes);
    appendScalar(stream, request.materialDensity);
    appendScalar(stream, request.originalUnitsPerNormalizedUnit);
    appendScalar(stream, request.minimumWallThickness);
    appendScalar(stream, request.centerOfMassWeight);
    appendScalar(stream, request.buoyancy.liquidDensity);
    appendScalar(stream, request.buoyancy.displacedLiquidVolume);
    appendMatrix(stream, request.buoyancy.centerOfBuoyancyTarget);
    appendScalar(stream, request.buoyancy.waterlineArea);
    appendMatrix(stream, request.buoyancy.waterlineCentralSecondMoment);
    appendScalar(stream, request.buoyancy.fixedShellVolume);
    stream << request.settings.maximumIterations << ';'
           << request.settings.maximumRestarts << ';'
           << request.settings.maximumFunctionEvaluationsPerRestart << ';';
    appendScalar(stream, request.settings.equationTolerance);
    appendScalar(stream, request.settings.physicalRelativeTolerance);
    appendScalar(stream, request.settings.lineSearchReduction);
    appendScalar(stream, request.settings.maximumGeodesicStep);
    stream << request.settings.deterministicSeed << ';';
    return manufacturing::computeContentSha256(stream.str());
}


} // namespace spin
