#pragma once

#include "physics/RawMassMoments10.h"

#include <Eigen/Core>

#include <array>
#include <limits>
#include <optional>
#include <string>
#include <vector>


namespace spin
{


enum class SpinMode
{
    StaticBalance,
    BuoyantStability,
    FidgetSpinner,
    YoYo,
    SpinningTop
};


const char* spinModeToken(SpinMode mode) noexcept;


enum class StaticObjective
{
    MaximumMass,
    LowerCenterOfMass,
    ParetoTradeoff
};


const char* staticObjectiveToken(StaticObjective objective) noexcept;


/** A closed outward-oriented triangle boundary. */
struct TriangleMesh
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};


/**
 * Right-handed normalized problem coordinates.
 *
 * The columns of axes are e1,e2,e3 in model coordinates.  For rotational
 * applications e3 is the requested spin axis; for static balance it points
 * opposite gravity through the requested support point at origin.
 */
struct SpinTargetFrame
{
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();
};


struct SpinSolverSettings
{
    int maximumIterations = 60;
    int maximumRestarts = 8;
    int maximumFunctionEvaluationsPerRestart = 60;
    double equationTolerance = 1e-12;
    double physicalRelativeTolerance = 1e-6;
    double lineSearchReduction = 0.5;
    double maximumGeodesicStep = 0.5;
    unsigned int deterministicSeed = 0;
};


/** Published prescribed-waterline inputs for buoyant stability. */
struct BuoyantStabilityData
{
    double liquidDensity = 0.0;
    double displacedLiquidVolume = 0.0;
    Eigen::Vector3d centerOfBuoyancyTarget = Eigen::Vector3d::Zero();
    double waterlineArea = 0.0;
    Eigen::Matrix2d waterlineCentralSecondMoment = Eigen::Matrix2d::Zero();
    double fixedShellVolume = 0.0;
};


/** Immutable input contract for the quadric hollowing solver. */
struct SpinDesignRequest
{
    TriangleMesh outerSurface;
    TriangleMesh designDomain;
    physics::RawMassMoments10 fixedShellMoments;
    SpinTargetFrame targetFrame;
    SpinMode mode = SpinMode::StaticBalance;
    StaticObjective staticObjective = StaticObjective::MaximumMass;
    /** Paper Fig. 11 weight kappa in -kappa*v+(1-kappa)*c3. */
    double staticParetoWeight = 0.0;
    double materialDensity = 1.0;
    double originalUnitsPerNormalizedUnit = 1.0;
    double minimumWallThickness = 0.0;

    /** Paper spinning-top weight w_c; unsupported for all other modes. */
    double centerOfMassWeight = 0.0;

    BuoyantStabilityData buoyancy;

    SpinSolverSettings settings;
};


enum class SpinSolverStatus
{
    Converged,
    Infeasible,
    InvalidRequest,
    OptimizationDidNotConverge,
    LineSearchExhausted,
    GeometryEvaluationFailed,
    InvalidMaterialArtifact,
    ResourceLimitReached,
    Indeterminate
};


const char* spinSolverStatusToken(SpinSolverStatus status) noexcept;


struct SpinIterationTrace
{
    int restart = 0;
    int iteration = 0;
    int functionEvaluations = 0;
    double normalizedEquationResidual = 0.0;
    double physicalRelativeResidual = 0.0;
    double objective = 0.0;
    double multiplierMu = 0.0;
    double acceptedGeodesicStep = 0.0;
    std::string quadricHash;
    std::string materialGeometryHash;
    std::string materialTopologyHash;
};


/** Deterministic scientific record; wall-clock timings are a separate sidecar. */
struct SpinSolverEvent
{
    std::string type;
    int restart = 0;
    int iteration = 0;
    int evaluationIndex = 0;
    int restartEvaluationIndex = 0;
    std::optional<double> step;
    std::string outcome;
    std::string reason;
    std::string detail;
    std::optional<std::array<double, 10>> quadricCoefficients;
    std::string quadricHash;
    std::optional<double> normalizedEquationResidual;
    std::optional<double> physicalRelativeResidual;
    std::optional<double> objective;
    std::optional<double> multiplierMu;
    std::optional<std::string> materialGeometryHash;
    std::optional<std::string> materialTopologyHash;
};


/** Non-deterministic timing sidecar for one actual evaluator call. */
struct SpinEvaluationTiming
{
    int evaluationIndex = 0;
    int restart = 0;
    int restartEvaluationIndex = 0;
    bool evaluatorReturned = false;
    double wallSeconds = 0.0;
    std::optional<double> reportedResidualSeconds;
    std::optional<double> reportedDerivativeSeconds;
};


/** Auditable terminal product; no field implies physical validation. */
struct QuadricHollowingResult
{
    SpinSolverStatus status = SpinSolverStatus::Indeterminate;
    std::string diagnostic;
    std::string inputHash;
    std::string geometryHash;
    /** Legacy v1 hash over accepted initial/iterate states only. */
    std::string traceHash;
    std::string outputHash;
    Eigen::Matrix<double, 10, 1> quadricCoefficients =
        Eigen::Matrix<double, 10, 1>::Zero();
    physics::RawMassMoments10 finalMoments;
    TriangleMesh outerSurface;
    TriangleMesh materialBoundary;
    TriangleMesh cavityBoundary;
    double normalizedEquationResidual = 0.0;
    double physicalRelativeResidual = 0.0;
    double buoyantMetacenterHeight =
        std::numeric_limits<double>::quiet_NaN();
    double buoyantStabilityMargin =
        std::numeric_limits<double>::quiet_NaN();
    /** Evaluator-reported component totals for Eq. 14 solve calls only. */
    double residualEvaluationSeconds = 0.0;
    double derivativeEvaluationSeconds = 0.0;
    int solveFunctionEvaluations = 0;
    std::vector<int> solveRestartFunctionEvaluations;
    double solveEvaluatorWallSeconds = 0.0;
    int postSolveMaterialEvaluationCount = 0;
    double postSolveMaterialEvaluationWallSeconds = 0.0;
    double postSolveMaterialResidualEvaluationSeconds = 0.0;
    double postSolveMaterialDerivativeEvaluationSeconds = 0.0;
    int eventStreamVersion = 1;
    std::string eventHash;
    std::vector<SpinSolverEvent> events;
    std::vector<SpinEvaluationTiming> evaluationTimings;
    std::vector<SpinIterationTrace> trace;
};


/** Throws std::invalid_argument unless every request invariant is satisfied. */
void validateSpinDesignRequest(const SpinDesignRequest& request);


/** SHA-256 of every normalized, solver-relevant request field. */
std::string computeSpinDesignRequestHash(const SpinDesignRequest& request);


} // namespace spin
