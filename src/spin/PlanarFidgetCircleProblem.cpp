#include "spin/PlanarFidgetCircleProblem.h"

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace spin
{
namespace
{

PlanarMassMoments6 subtract(const PlanarMassMoments6& a,
                            const PlanarMassMoments6& b)
{
    return {a.mass-b.mass, a.firstX-b.firstX, a.firstY-b.firstY,
            a.secondXX-b.secondXX, a.secondYY-b.secondYY,
            a.secondXY-b.secondXY};
}

bool finite(const PlanarMassMoments6& value)
{
    return value.vector().allFinite();
}

void validateOptions(const PlanarDoglegOptions& options)
{
    if (options.maximumIterations == 0 || options.maximumFunctionEvaluations == 0 ||
        !std::isfinite(options.residualTolerance) || options.residualTolerance <= 0.0 ||
        !std::isfinite(options.initialTrustRadius) || options.initialTrustRadius <= 0.0 ||
        !std::isfinite(options.maximumTrustRadius) ||
            options.maximumTrustRadius < options.initialTrustRadius ||
        !std::isfinite(options.minimumTrustRadius) || options.minimumTrustRadius <= 0.0 ||
        options.minimumTrustRadius >= options.initialTrustRadius ||
        !std::isfinite(options.acceptanceRatio) ||
            options.acceptanceRatio < 0.0 || options.acceptanceRatio >= 1.0)
        throw std::invalid_argument("invalid planar dogleg options");
}

} // namespace

PlanarFidgetCircleProblem::PlanarFidgetCircleProblem(
    std::vector<PlanarBoundaryLoop> designRegion,
    PlanarMassMoments6 fixedMaterial,
    double designArealDensity)
    : designRegion_(std::move(designRegion)),
      fixedMaterial_(fixedMaterial),
      designMoments_(integratePlanarBoundaryLoops(designRegion_)),
      designArealDensity_(designArealDensity)
{
    if (!finite(fixedMaterial_) || !(fixedMaterial_.mass > 0.0))
        throw std::invalid_argument("fixed planar material must have finite positive mass");
    if (!std::isfinite(designArealDensity_) || designArealDensity_ <= 0.0)
        throw std::invalid_argument("design areal density must be finite and positive");

    double loX = std::numeric_limits<double>::infinity();
    double loY = std::numeric_limits<double>::infinity();
    double hiX = -std::numeric_limits<double>::infinity();
    double hiY = -std::numeric_limits<double>::infinity();
    for (const auto& loop : designRegion_)
    {
        for (const auto& point : loop)
        {
            loX = std::min(loX, point.x()); hiX = std::max(hiX, point.x());
            loY = std::min(loY, point.y()); hiY = std::max(hiY, point.y());
        }
    }
    characteristicLength_ = std::hypot(hiX-loX, hiY-loY);
    if (!std::isfinite(characteristicLength_) || !(characteristicLength_ > 0.0))
        throw std::invalid_argument("planar design region must have positive finite extent");

    initialMaterial_ = fixedMaterial_;
    initialMaterial_ += designMoments_.weighted(designArealDensity_);
    if (!finite(initialMaterial_) || !(initialMaterial_.mass > 0.0))
        throw std::invalid_argument("combined planar initial material moments are invalid");
    const long double residualScale = static_cast<long double>(initialMaterial_.mass) *
                                     characteristicLength_;
    if (!std::isfinite(residualScale) || residualScale <= 0.0L ||
        residualScale > std::numeric_limits<double>::max())
        throw std::overflow_error("planar residual normalization scale is out of range");
    residualScale_ = static_cast<double>(residualScale);
    if (!(residualScale_ > 0.0))
        throw std::underflow_error("planar residual normalization scale underflowed");
}

PlanarFidgetCircleEvaluation PlanarFidgetCircleProblem::evaluate(
    const Eigen::Vector2d& center) const
{
    if (!center.allFinite())
        throw std::invalid_argument("planar fidget circle center must be finite");
    const double radius = center.norm();
    if (!std::isfinite(radius))
        throw std::invalid_argument("planar fidget circle radius is nonfinite");

    PlanarMassMoments6 removed;
    if (radius > 0.0)
    {
        removed = integratePlanarDiskIntersection(
            designRegion_, PlanarCircleDisk{center, radius}).weighted(designArealDensity_);
    }
    const PlanarMassMoments6 remaining = subtract(initialMaterial_, removed);
    if (!finite(remaining) || !(remaining.mass > 0.0))
        throw std::runtime_error("planar fidget carving must leave positive finite mass");

    PlanarFidgetCircleEvaluation result;
    result.center = center;
    result.radius = radius;
    result.removedMaterial = removed;
    result.remainingMaterial = remaining;
    result.firstMomentResidual << remaining.firstX, remaining.firstY;
    result.normalizedResidual = result.firstMomentResidual/residualScale_;
    return result;
}

PlanarDiskIntersectionTopology PlanarFidgetCircleProblem::topologyAt(
    const Eigen::Vector2d& center) const
{
    if (!center.allFinite())
        throw std::invalid_argument("planar fidget circle center must be finite");
    const double radius = center.norm();
    if (!std::isfinite(radius))
        throw std::invalid_argument("planar fidget circle radius is nonfinite");
    return classifyPlanarDiskIntersectionTopology(
        designRegion_, PlanarCircleDisk{center, radius});
}

PlanarDoglegResult solvePlanarFidgetCircle(
    const PlanarFidgetCircleProblem& problem,
    const Eigen::Vector2d& initialCenter,
    const PlanarDoglegOptions& options)
{
    validateOptions(options);
    if (!initialCenter.allFinite())
        throw std::invalid_argument("initial planar dogleg center must be finite");

    PlanarDoglegResult result;
    result.center = initialCenter;
    double trustRadius = options.initialTrustRadius;
    const double scale = problem.characteristicLength();
    const double fdStep = std::cbrt(std::numeric_limits<double>::epsilon());
    auto evaluateScaled = [&](const Eigen::Vector2d& z, const char* purpose,
                              int differenceAxis = -1) {
        ++result.functionEvaluations;
        PlanarDoglegEvaluationEvent event;
        event.evaluationIndex = result.functionEvaluations;
        event.iteration = result.iterations;
        event.purpose = purpose;
        event.differenceAxis = differenceAxis;
        event.center = scale*z;
        result.evaluations.push_back(std::move(event));
        auto& recorded = result.evaluations.back();
        Eigen::Vector2d value;
        try
        {
            value = problem.evaluate(recorded.center).normalizedResidual;
        }
        catch (const std::exception& error)
        {
            recorded.outcome = "failed";
            recorded.diagnostic = error.what();
            if (options.evaluationSink) options.evaluationSink(recorded);
            throw;
        }
        recorded.normalizedResidual = value;
        recorded.hasResidual = true;
        recorded.outcome = "evaluated";
        // Finalize trial acceptance before emitting it. Observer failures are
        // not classified as failed numerical evaluations.
        if (std::string(purpose) != "trial" && options.evaluationSink)
            options.evaluationSink(recorded);
        return value;
    };

    Eigen::Vector2d z = initialCenter/scale;
    Eigen::Vector2d residual = evaluateScaled(z, "initial");
    while (true)
    {
        result.center = scale*z;
        result.radius = result.center.norm();
        result.normalizedResidual = residual.norm();
        result.finalTrustRadius = trustRadius;
        if (result.normalizedResidual <= options.residualTolerance)
        {
            result.status = PlanarDoglegStatus::Converged;
            return result;
        }
        if (result.iterations >= options.maximumIterations)
        {
            result.status = PlanarDoglegStatus::MaximumIterations;
            return result;
        }
        if (result.functionEvaluations + 5 > options.maximumFunctionEvaluations)
        {
            result.status = PlanarDoglegStatus::MaximumFunctionEvaluations;
            return result;
        }

        const auto topology = problem.topologyAt(scale*z);
        if (topology.nearTransition)
        {
            result.status = PlanarDoglegStatus::TopologyTransition;
            return result;
        }

        Eigen::Matrix2d jacobian;
        for (int axis = 0; axis < 2; ++axis)
        {
            const double axisScale = std::max(1.0, std::abs(z[axis]));
            double h = fdStep*axisScale;
            const double minimumH = 64.0*std::numeric_limits<double>::epsilon()*axisScale;
            Eigen::Vector2d plus, minus;
            bool validStencil = false;
            while (h >= minimumH)
            {
                plus = z;
                minus = z;
                plus[axis] += h;
                minus[axis] -= h;
                const auto plusTopology = problem.topologyAt(scale*plus);
                const auto minusTopology = problem.topologyAt(scale*minus);
                const bool sameTopology = topology.signature == "ZERO_RADIUS"
                    ? plusTopology.signature == minusTopology.signature
                    : plusTopology.signature == topology.signature &&
                      minusTopology.signature == topology.signature;
                if (!plusTopology.nearTransition && !minusTopology.nearTransition &&
                    sameTopology)
                {
                    validStencil = true;
                    break;
                }
                h *= 0.5;
            }
            if (!validStencil)
            {
                result.status = PlanarDoglegStatus::TopologyTransition;
                return result;
            }
            const Eigen::Vector2d fPlus = evaluateScaled(plus, "finite_difference_plus", axis);
            const Eigen::Vector2d fMinus = evaluateScaled(minus, "finite_difference_minus", axis);
            jacobian.col(axis) = (fPlus-fMinus)/(2.0*h);
        }

        const Eigen::Vector2d gradient = jacobian.transpose()*residual;
        const Eigen::Matrix2d gaussNewton = jacobian.transpose()*jacobian;
        const Eigen::FullPivLU<Eigen::Matrix2d> factorization(jacobian);
        if (!factorization.isInvertible())
        {
            result.status = PlanarDoglegStatus::SingularJacobian;
            result.center = scale*z;
            result.radius = result.center.norm();
            result.normalizedResidual = residual.norm();
            result.finalTrustRadius = trustRadius;
            return result;
        }

        const Eigen::Vector2d newton = factorization.solve(-residual);
        const double newtonNorm = newton.norm();
        Eigen::Vector2d step;
        if (newtonNorm <= trustRadius)
        {
            step = newton;
        }
        else
        {
            const double gradientSquared = gradient.squaredNorm();
            const double curvature = gradient.dot(gaussNewton*gradient);
            if (!(gradientSquared > 0.0) || !(curvature > 0.0) ||
                !std::isfinite(curvature))
            {
                result.status = PlanarDoglegStatus::SingularJacobian;
                result.center = scale*z;
                result.radius = result.center.norm();
                result.normalizedResidual = residual.norm();
                result.finalTrustRadius = trustRadius;
                return result;
            }
            const Eigen::Vector2d cauchy = -(gradientSquared/curvature)*gradient;
            const double cauchyNorm = cauchy.norm();
            if (cauchyNorm >= trustRadius)
            {
                step = -(trustRadius/std::sqrt(gradientSquared))*gradient;
            }
            else
            {
                const Eigen::Vector2d direction = newton-cauchy;
                const double a = direction.squaredNorm();
                const double b = 2.0*cauchy.dot(direction);
                const double c = cauchy.squaredNorm()-trustRadius*trustRadius;
                const double discriminant = b*b-4.0*a*c;
                if (!(a > 0.0) || discriminant < 0.0 || !std::isfinite(discriminant))
                {
                    result.status = PlanarDoglegStatus::SingularJacobian;
                    result.center = scale*z;
                    result.radius = result.center.norm();
                    result.normalizedResidual = residual.norm();
                    result.finalTrustRadius = trustRadius;
                    return result;
                }
                const double tau = (-b+std::sqrt(discriminant))/(2.0*a);
                step = cauchy+tau*direction;
            }
        }

        const Eigen::Vector2d linearized = residual+jacobian*step;
        const double currentObjective = 0.5*residual.squaredNorm();
        const double predictedReduction = currentObjective-0.5*linearized.squaredNorm();
        if (!(predictedReduction > 0.0) || !std::isfinite(predictedReduction))
        {
            result.status = PlanarDoglegStatus::SingularJacobian;
            result.center = scale*z;
            result.radius = result.center.norm();
            result.normalizedResidual = residual.norm();
            result.finalTrustRadius = trustRadius;
            return result;
        }

        const Eigen::Vector2d trial = z+step;
        const Eigen::Vector2d trialResidual = evaluateScaled(trial, "trial");
        const double actualReduction = currentObjective-0.5*trialResidual.squaredNorm();
        const double ratio = actualReduction/predictedReduction;
        if (ratio < 0.25)
            trustRadius *= 0.25;
        else if (ratio > 0.75 && std::abs(step.norm()-trustRadius) <=
                 32.0*std::numeric_limits<double>::epsilon()*trustRadius)
            trustRadius = std::min(2.0*trustRadius, options.maximumTrustRadius);

        if (ratio > options.acceptanceRatio && actualReduction > 0.0)
        {
            z = trial;
            residual = trialResidual;
            ++result.acceptedSteps;
            result.evaluations.back().acceptedTrial = true;
        }
        if (options.evaluationSink) options.evaluationSink(result.evaluations.back());
        ++result.iterations;
        if (trustRadius < options.minimumTrustRadius)
        {
            result.center = scale*z;
            result.radius = result.center.norm();
            result.normalizedResidual = residual.norm();
            result.finalTrustRadius = trustRadius;
            result.status = PlanarDoglegStatus::TrustRegionExhausted;
            return result;
        }
    }
}

const char* toString(PlanarDoglegStatus status) noexcept
{
    switch (status)
    {
        case PlanarDoglegStatus::Converged: return "CONVERGED";
        case PlanarDoglegStatus::MaximumIterations: return "MAXIMUM_ITERATIONS";
        case PlanarDoglegStatus::MaximumFunctionEvaluations: return "MAXIMUM_FUNCTION_EVALUATIONS";
        case PlanarDoglegStatus::TrustRegionExhausted: return "TRUST_REGION_EXHAUSTED";
        case PlanarDoglegStatus::SingularJacobian: return "SINGULAR_JACOBIAN";
        case PlanarDoglegStatus::TopologyTransition: return "TOPOLOGY_TRANSITION";
    }
    return "UNKNOWN";
}

} // namespace spin
