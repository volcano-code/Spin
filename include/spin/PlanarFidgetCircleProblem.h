#pragma once

#include "spin/PlanarDiskIntersection.h"

#include <Eigen/Core>
#include <cstddef>
#include <vector>
#include <string>
#include <functional>

namespace spin
{

struct PlanarFidgetCircleEvaluation
{
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    double radius = 0.0;
    PlanarMassMoments6 removedMaterial;
    PlanarMassMoments6 remainingMaterial;
    Eigen::Vector2d firstMomentResidual = Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedResidual = Eigen::Vector2d::Zero();
};

/**
 * Two-dimensional fidget-spinner balance equations.
 *
 * The design cavity is a disk centered at `center` with radius `|center|`, so
 * its boundary passes through the target origin. Fixed-material moments and
 * the uncut design-region moments form the initial object. The solve enforces
 * zero remaining in-plane first moment; it is a local root solve, not a global
 * optimality certificate.
 */
class PlanarFidgetCircleProblem
{
public:
    PlanarFidgetCircleProblem(std::vector<PlanarBoundaryLoop> designRegion,
                              PlanarMassMoments6 fixedMaterial,
                              double designArealDensity);

    PlanarFidgetCircleEvaluation evaluate(const Eigen::Vector2d& center) const;
    PlanarDiskIntersectionTopology topologyAt(const Eigen::Vector2d& center) const;
    double characteristicLength() const noexcept { return characteristicLength_; }
    const PlanarMassMoments6& initialMaterial() const noexcept { return initialMaterial_; }

private:
    std::vector<PlanarBoundaryLoop> designRegion_;
    PlanarMassMoments6 fixedMaterial_;
    PlanarRawMoments6 designMoments_;
    PlanarMassMoments6 initialMaterial_;
    double designArealDensity_ = 0.0;
    double characteristicLength_ = 0.0;
    double residualScale_ = 0.0;
};

enum class PlanarDoglegStatus
{
    Converged,
    MaximumIterations,
    MaximumFunctionEvaluations,
    TrustRegionExhausted,
    SingularJacobian,
    /** A finite-difference derivative would cross a polygon/disk event. */
    TopologyTransition
};

struct PlanarDoglegEvaluationEvent;
using PlanarDoglegEvaluationSink = std::function<void(const PlanarDoglegEvaluationEvent&)>;

/** Numerical controls are implementation defaults, not recovered author settings. */
struct PlanarDoglegOptions
{
    std::size_t maximumIterations = 100;
    std::size_t maximumFunctionEvaluations = 60;
    double residualTolerance = 1e-10;
    double initialTrustRadius = 0.25;
    double maximumTrustRadius = 1.0;
    double minimumTrustRadius = 1e-12;
    double acceptanceRatio = 0.1;
    /** Finalized calls, including thrown evaluations. Empty preserves the
     * legacy observer-free behavior and original numerical exception types. */
    PlanarDoglegEvaluationSink evaluationSink;
};

struct PlanarDoglegEvaluationEvent
{
    std::size_t evaluationIndex = 0;
    std::size_t iteration = 0;
    std::string purpose;
    int differenceAxis = -1;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedResidual = Eigen::Vector2d::Zero();
    bool hasResidual = false;
    std::string outcome = "started";
    bool acceptedTrial = false;
    std::string diagnostic;
};

struct PlanarDoglegResult
{
    PlanarDoglegStatus status = PlanarDoglegStatus::MaximumIterations;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    double radius = 0.0;
    double normalizedResidual = 0.0;
    std::size_t iterations = 0;
    std::size_t functionEvaluations = 0;
    std::size_t acceptedSteps = 0;
    double finalTrustRadius = 0.0;
    std::vector<PlanarDoglegEvaluationEvent> evaluations;
};

PlanarDoglegResult solvePlanarFidgetCircle(
    const PlanarFidgetCircleProblem& problem,
    const Eigen::Vector2d& initialCenter,
    const PlanarDoglegOptions& options = {});

const char* toString(PlanarDoglegStatus status) noexcept;

} // namespace spin
