#include "spin/SpinDesignRequest.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"
#include "physics/MassProperties.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactBoundaryAudit.h"
#endif

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
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
    manufacturing::BoundaryOrientation orientation,
    bool exactInputBoundary
)
{
    // Two immutable, process-local input proofs cover the usual outer/domain
    // pair. Keys bind every indexed coordinate/face and the audit policy.
    // All request scalars/frame/objective checks still run on every request.
    struct InputMeshProofCache {
        std::array<std::string,2> keys;
        std::size_t replacement=0;
    };
    static thread_local InputMeshProofCache proofs;
    std::ostringstream key;key.imbue(std::locale::classic());
    key<<"FixedInputBoundaryProof-v1;"<<static_cast<int>(orientation)<<';'
       <<exactInputBoundary<<';';
    appendMatrix(key,mesh.vertices);appendMatrix(key,mesh.faces);
    const std::string indexedBoundary=key.str();
    if(std::find(proofs.keys.begin(),proofs.keys.end(),indexedBoundary)!=proofs.keys.end())return;
    const auto remember=[&]() {
        proofs.keys[proofs.replacement]=indexedBoundary;
        proofs.replacement=(proofs.replacement+1)%proofs.keys.size();
    };
    if (exactInputBoundary) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        const auto audit = verification::auditExactTriangleBoundary(mesh.vertices, mesh.faces,
            verification::ExactBoundaryOrientation::OutwardPositive);
        if (!audit.facts().boundaryValid)
            throw std::invalid_argument(std::string(name) + " failed exact input boundary validation");
        remember();
        return;
#else
        throw std::invalid_argument("exact input boundary policy requires exact geometry checks");
#endif
    }
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
    remember();
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

const char* spinEquationProfileToken(SpinEquationProfile profile) noexcept
{
    switch(profile) {
    case SpinEquationProfile::LegacyNativeNormalizedV1: return "LEGACY_NATIVE_NORMALIZED_V1";
    case SpinEquationProfile::DualMonomialRawV1: return "DUAL_MONOMIAL_RAW_V1";
    }
    return "INVALID";
}

const char* spinSolverCoordinateIdentity(SpinEquationProfile profile) noexcept
{
    switch(profile) {
    case SpinEquationProfile::LegacyNativeNormalizedV1: return "native_matrix_entry_unit_sphere_v1";
    case SpinEquationProfile::DualMonomialRawV1: return "literal_monomial_unit_sphere_v1";
    }
    return "INVALID";
}


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
    case StaticObjective::ParetoTradeoff:
        return "PARETO_TRADEOFF";
    }
    return "UNKNOWN";
}


const char* spinSolverStatusToken(SpinSolverStatus status) noexcept
{
    switch (status)
    {
    case SpinSolverStatus::Converged:
        return "CONVERGED";
    case SpinSolverStatus::Infeasible:
        return "INFEASIBLE";
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

const char* staticObjectiveMeasureToken(StaticObjectiveMeasure measure) noexcept {
    switch(measure) {
    case StaticObjectiveMeasure::DensityWeightedMass: return "DENSITY_WEIGHTED_MASS_V1";
    case StaticObjectiveMeasure::UniformDensityVolume: return "UNIFORM_DENSITY_VOLUME_V1";
    }
    return "UNKNOWN";
}

double staticObjectiveMassScale(const SpinDesignRequest& request) {
    if(request.staticObjectiveMeasure==StaticObjectiveMeasure::DensityWeightedMass)return 1.0;
    if(request.staticObjectiveMeasure!=StaticObjectiveMeasure::UniformDensityVolume ||
       request.mode!=SpinMode::StaticBalance || !std::isfinite(request.materialDensity) ||
       request.materialDensity<=0.0 || !std::isfinite(1.0/request.materialDensity))
        throw std::invalid_argument("invalid static volume objective measure or density");
    return 1.0/request.materialDensity;
}


void validateSpinningTopObjective(const SpinDesignRequest& request)
{
    if (request.spinningTopObjective != SpinningTopObjective::InertiaAndLowering &&
        request.spinningTopObjective != SpinningTopObjective::LowerCenterOfMass)
        throw std::invalid_argument("unsupported spinning-top objective");
    if (request.spinningTopObjective == SpinningTopObjective::LowerCenterOfMass &&
        (request.mode != SpinMode::SpinningTop || request.centerOfMassWeight != 0.0))
        throw std::invalid_argument("pure-height objective requires spinning-top mode and zero extra weight");
}


void validateSpinDesignRequest(const SpinDesignRequest& request)
{
    validateSpinDesignRequest(request, false);
}

void validateSpinDesignRequest(const SpinDesignRequest& request, bool exactInputBoundary)
{
    validateMesh(
        request.outerSurface,
        "outer surface",
        manufacturing::BoundaryOrientation::OutwardPositive, exactInputBoundary
    );
    validateMesh(
        request.designDomain,
        "design domain",
        manufacturing::BoundaryOrientation::OutwardPositive, exactInputBoundary
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
    validateSpinningTopObjective(request);
    requireFinite(request.staticParetoWeight, "static Pareto weight");
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
    (void)staticObjectiveMassScale(request);
    if(request.staticObjectiveMeasure==StaticObjectiveMeasure::UniformDensityVolume) {
        const Eigen::Matrix3d rotation=request.targetFrame.axes.transpose();
        const Eigen::Vector3d translation=-rotation*request.targetFrame.origin;
        const Eigen::MatrixXd outer=(request.outerSurface.vertices*rotation.transpose()).rowwise()+translation.transpose();
        const Eigen::MatrixXd design=(request.designDomain.vertices*rotation.transpose()).rowwise()+translation.transpose();
        const Eigen::Matrix<double,10,1> expected=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(
            outer,request.outerSurface.faces,request.materialDensity)).values()-
            physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(
            design,request.designDomain.faces,request.materialDensity)).values();
        const Eigen::Matrix<double,10,1> supplied=request.fixedShellMoments.transformed(rotation,translation).values();
        const double length=(outer.colwise().maxCoeff()-outer.colwise().minCoeff()).maxCoeff();
        for(int i=0;i<10;++i) {
            const double scale=request.materialDensity*std::pow(length,i==0?3:i<4?4:5);
            if(!std::isfinite(scale)||scale<=0.0||!std::isfinite(expected[i])||
               std::abs(expected[i]-supplied[i])/scale>1e-10)
                throw std::invalid_argument("volume objective requires uniform-density fixed shell integrated from outer minus design");
        }
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
    if (
        request.staticObjective == StaticObjective::ParetoTradeoff &&
        (request.staticParetoWeight < 0.0 || request.staticParetoWeight > 1.0)
    )
    {
        throw std::invalid_argument("static Pareto weight must lie in [0,1]");
    }
    if (
        request.staticObjective != StaticObjective::ParetoTradeoff &&
        request.staticParetoWeight != 0.0
    )
    {
        throw std::invalid_argument(
            "static Pareto weight is supported only for PARETO_TRADEOFF"
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
        requireFinite(buoyancy.fixedShellDensity, "fixed shell density");
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
            buoyancy.fixedShellVolume < 0.0 ||
            buoyancy.fixedShellDensity < 0.0
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
            (buoyancy.fixedShellDensity > 0.0
                ? buoyancy.fixedShellDensity
                : request.materialDensity) * buoyancy.fixedShellVolume;
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
        buoyancy.fixedShellVolume != 0.0 ||
        buoyancy.fixedShellDensity != 0.0
    )
    {
        throw std::invalid_argument(
            "buoyancy inputs are supported only in buoyant-stability mode"
        );
    }

    const SpinSolverSettings& settings = request.settings;
    if ((settings.equationProfile != SpinEquationProfile::LegacyNativeNormalizedV1 &&
         settings.equationProfile != SpinEquationProfile::DualMonomialRawV1) ||
        (settings.equationProfile == SpinEquationProfile::DualMonomialRawV1 &&
         settings.maximumFunctionEvaluationsPerRestart > 60))
        throw std::invalid_argument("invalid equation profile or strict evaluation budget");
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
    return computeSpinDesignRequestHash(request, false);
}

std::string computeSpinDesignRequestHash(const SpinDesignRequest& request, bool exactInputBoundary)
{
    validateSpinDesignRequest(request, exactInputBoundary);
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    if (exactInputBoundary) stream << "ExactInputBoundary-binary64-outward-v1;";
    if(request.settings.equationProfile == SpinEquationProfile::DualMonomialRawV1)
        stream << "SpinDesignRequest-v4;" << spinEquationProfileToken(request.settings.equationProfile)
               << ';' << spinSolverCoordinateIdentity(request.settings.equationProfile) << ';';
    else stream << "SpinDesignRequest-v3;";
    stream << spinModeToken(request.mode) << ';';
    stream << staticObjectiveToken(request.staticObjective) << ';';
    if(request.staticObjectiveMeasure!=StaticObjectiveMeasure::DensityWeightedMass)
        stream << staticObjectiveMeasureToken(request.staticObjectiveMeasure) << ';';
    appendScalar(stream, request.staticParetoWeight);
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
    if (request.spinningTopObjective == SpinningTopObjective::LowerCenterOfMass)
        stream << "SpinningTopObjective:LOWER_CENTER_OF_MASS_V1;";
    appendScalar(stream, request.buoyancy.liquidDensity);
    appendScalar(stream, request.buoyancy.displacedLiquidVolume);
    appendMatrix(stream, request.buoyancy.centerOfBuoyancyTarget);
    appendScalar(stream, request.buoyancy.waterlineArea);
    appendMatrix(stream, request.buoyancy.waterlineCentralSecondMoment);
    appendScalar(stream, request.buoyancy.fixedShellVolume);
    appendScalar(stream, request.buoyancy.fixedShellDensity);
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
