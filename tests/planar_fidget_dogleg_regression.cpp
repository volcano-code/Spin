#include "spin/PlanarFidgetCircleProblem.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

spin::PlanarBoundaryLoop rectangle(double x0, double y0, double x1, double y1)
{
    return {{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
}

bool close(double actual, double expected, double tolerance)
{
    return std::abs(actual-expected) <= tolerance;
}

int fail(const char* reason)
{
    std::cerr << reason << '\n';
    return 1;
}

spin::PlanarFidgetCircleProblem knownRootProblem()
{
    const auto domain = rectangle(-2.0, -2.0, 2.0, 2.0);
    const double pi = std::acos(-1.0);
    const auto fixedPatch = rectangle(pi/4.0-0.1, -0.1, pi/4.0+0.1, 0.1);
    const auto fixed = spin::integratePlanarBoundaryLoops({fixedPatch}).weighted(100.0);
    return spin::PlanarFidgetCircleProblem({domain}, fixed, 1.0);
}

spin::PlanarFidgetCircleProblem tangentTransitionProblem()
{
    const auto domain = rectangle(-2.0, -2.0, 2.0, 2.0);
    const double pi = std::acos(-1.0);
    // This deliberately differs from the balanced patch above.  At (1,0)
    // the through-origin circle is tangent to x=2, while the residual is
    // nonzero and a Jacobian would otherwise be requested.
    const auto fixedPatch = rectangle(pi/4.0, -0.1, pi/4.0+0.2, 0.1);
    const auto fixed = spin::integratePlanarBoundaryLoops({fixedPatch}).weighted(100.0);
    return spin::PlanarFidgetCircleProblem({domain}, fixed, 1.0);
}

} // namespace

int main()
{
    auto problem = knownRootProblem();
    const auto exact = problem.evaluate(Eigen::Vector2d(1.0, 0.0));
    if (exact.normalizedResidual.norm() > 2e-16 ||
        !close(exact.radius, 1.0, 1e-15))
        return fail("known circle-through-origin fixture is not balanced");

    spin::PlanarDoglegOptions options;
    options.initialTrustRadius = 0.25;
    options.maximumTrustRadius = 1.0;
    options.residualTolerance = 1e-10;
    const auto solved = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d(0.75, 0.2), options);
    if (solved.status != spin::PlanarDoglegStatus::Converged)
    {
        std::cerr << "dogleg status=" << spin::toString(solved.status)
                  << " center=" << solved.center.transpose()
                  << " residual=" << solved.normalizedResidual
                  << " iterations=" << solved.iterations << '\n';
        return fail("trust-region dogleg did not converge to the known balanced circle");
    }
    if (!close(solved.center.x(), 1.0, 2e-7) ||
        !close(solved.center.y(), 0.0, 2e-7) ||
        solved.normalizedResidual > options.residualTolerance)
        return fail("dogleg converged to an incorrect circle or residual");
    if (solved.functionEvaluations < 1 || solved.functionEvaluations > 60)
        return fail("dogleg function-evaluation count escaped its declared budget");
    const double centerError = (solved.center - Eigen::Vector2d(1.0, 0.0)).norm();

    // The geometry-backed application fixture is a square minus a fixed patch
    // (the clockwise hole). At the solved circle, the patch is wholly inside
    // the removed disk, leaving exactly two closed boundary cycles: the outer
    // square and the circular cavity.
    const auto outer = rectangle(-3.0, -3.0, 3.0, 3.0);
    const double pi = std::acos(-1.0);
    auto patch = rectangle(pi/4.0-0.1, -0.1, pi/4.0+0.1, 0.1);
    std::reverse(patch.begin(), patch.end());
    const std::vector<spin::PlanarBoundaryLoop> designRegion{outer, patch};
    const spin::PlanarCircleDisk disk{solved.center, solved.radius};
    const auto clippedBoundary = spin::planarRegionDifferenceBoundary(designRegion, disk);
    if (!clippedBoundary.defined || clippedBoundary.cycles.size() != 2)
        return fail("polygon-minus-disk did not produce the two expected closed cycles");
    const auto boundaryMoments = spin::integratePlanarBoundaryCycles(clippedBoundary.cycles);
    const auto designMoments = spin::integratePlanarBoundaryLoops(designRegion);
    const auto removedMoments = spin::integratePlanarDiskIntersection(designRegion, disk);
    const double expectedMoments[] = {
        designMoments.area-removedMoments.area,
        designMoments.firstX-removedMoments.firstX,
        designMoments.firstY-removedMoments.firstY,
        designMoments.secondXX-removedMoments.secondXX,
        designMoments.secondYY-removedMoments.secondYY,
        designMoments.secondXY-removedMoments.secondXY,
    };
    const double actualMoments[] = {
        boundaryMoments.area, boundaryMoments.firstX, boundaryMoments.firstY,
        boundaryMoments.secondXX, boundaryMoments.secondYY, boundaryMoments.secondXY,
    };
    for (std::size_t i = 0; i < std::size(expectedMoments); ++i)
    {
        if (!close(actualMoments[i], expectedMoments[i], 2e-10))
            return fail("closed clipped boundary moments differ from independent polygon/disk integration");
    }
    auto fixedPatch = patch;
    std::reverse(fixedPatch.begin(), fixedPatch.end());
    if (spin::planarRegionsHavePositiveAreaIntersection(
            designRegion, {fixedPatch}))
        return fail("a fixed patch in a design-domain hole was falsely reported overlapping");
    if (!spin::planarRegionsHavePositiveAreaIntersection({outer}, {fixedPatch}))
        return fail("positive-area fixed/design overlap was not detected");
    const auto touching = rectangle(3.0, -1.0, 4.0, 1.0);
    if (spin::planarRegionsHavePositiveAreaIntersection({outer}, {touching}))
        return fail("boundary-only contact was mislabeled as positive-area overlap");

    const auto alreadySolved = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d(1.0, 0.0), options);
    if (alreadySolved.status != spin::PlanarDoglegStatus::Converged ||
        alreadySolved.iterations != 0)
        return fail("an initially balanced circle should terminate without a trial step");

    options.maximumFunctionEvaluations = 1;
    const auto budget = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d(0.75, 0.2), options);
    if (budget.status != spin::PlanarDoglegStatus::MaximumFunctionEvaluations ||
        budget.normalizedResidual <= options.residualTolerance)
        return fail("evaluation exhaustion was mislabeled as a converged solution");

    options.maximumFunctionEvaluations = 60;
    const auto zeroRadius = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d::Zero(), options);
    if (zeroRadius.status != spin::PlanarDoglegStatus::Converged ||
        !close(zeroRadius.center.x(), 1.0, 2e-7) ||
        !close(zeroRadius.center.y(), 0.0, 2e-7) ||
        zeroRadius.normalizedResidual > options.residualTolerance)
        return fail("dogleg did not leave the zero-radius initial point for the balanced circle");

    const auto tangent = spin::solvePlanarFidgetCircle(
        tangentTransitionProblem(), Eigen::Vector2d(1.0, 0.0), options);
    if (tangent.status != spin::PlanarDoglegStatus::TopologyTransition ||
        tangent.iterations != 0 || tangent.normalizedResidual <= options.residualTolerance)
    {
        return fail("dogleg reused a finite-difference Jacobian across a polygon/disk topology transition");
    }

    std::cout << std::setprecision(17)
              << "FASTER_PLANAR_COMPONENT_JSON={\"schema\":\"faster-planar-component-result-v1\","
              << "\"classification\":\"SYNTHETIC_REGRESSION\","
              << "\"component\":\"dogleg\","
              << "\"status\":\"CONVERGED\","
              << "\"solved_center\":[" << solved.center.x() << ',' << solved.center.y() << "],"
              << "\"known_root_center\":[1,0],"
              << "\"center_error\":" << centerError << ','
              << "\"normalized_residual\":" << solved.normalizedResidual << ','
              << "\"residual_tolerance\":" << options.residualTolerance << ','
              << "\"iterations\":" << solved.iterations << ','
              << "\"function_evaluations\":" << solved.functionEvaluations << ','
              << "\"function_evaluation_budget\":60,"
              << "\"exact_root_residual\":" << exact.normalizedResidual.norm() << ','
              << "\"zero_radius_initial_converged\":true,"
              << "\"budget_exhaustion_distinguished\":true,"
              << "\"topology_transition_detected\":true}\n";

    try
    {
        const spin::PlanarMassMoments6 enormous{1e308, 0.0, 0.0,
                                                1e298, 1e298, 0.0};
        (void)spin::PlanarFidgetCircleProblem(
            {rectangle(-4.0, -4.0, 4.0, 4.0)}, enormous, 1e-300);
        return fail("overflowed residual normalization scale was accepted");
    }
    catch (const std::overflow_error&)
    {
    }

    // The public solver's observer must retain every real residual evaluation,
    // including rejected trials, without changing the accepted trajectory.
    std::vector<spin::PlanarDoglegEvaluationEvent> observed;
    spin::PlanarDoglegOptions recordingOptions;
    recordingOptions.evaluationSink = [&](const auto& event) { observed.push_back(event); };
    const auto recorded = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d(0.75, 0.2), recordingOptions);
    if (recorded.center != solved.center ||
        recorded.normalizedResidual != solved.normalizedResidual ||
        recorded.functionEvaluations != solved.functionEvaluations ||
        recorded.evaluations.size() != recorded.functionEvaluations ||
        observed.size() != recorded.evaluations.size())
        return fail("recording changed the solver or lost actual residual calls");
    std::size_t accepted = 0, differences = 0;
    for (std::size_t i = 0; i < observed.size(); ++i)
    {
        const auto& event = observed[i];
        const auto& saved = recorded.evaluations[i];
        if (event.evaluationIndex != i+1 || event.outcome != "evaluated" ||
            !event.hasResidual || event.center != saved.center ||
            event.normalizedResidual != saved.normalizedResidual ||
            event.acceptedTrial != saved.acceptedTrial)
            return fail("observer events differ from finalized solver evaluations");
        accepted += event.acceptedTrial;
        if (event.purpose == "finite_difference_plus" ||
            event.purpose == "finite_difference_minus")
        {
            ++differences;
            if (event.differenceAxis < 0 || event.differenceAxis > 1)
                return fail("finite-difference event lost its parameter axis");
        }
    }
    if (accepted != recorded.acceptedSteps || differences != 4*recorded.iterations)
        return fail("accepted or finite-difference evaluations were dropped");

    // A test-only stringent ratio exercises rejected trials; it is not used
    // by any application request or paper profile.
    recordingOptions.acceptanceRatio = 0.999999;
    observed.clear();
    const auto rejected = spin::solvePlanarFidgetCircle(
        problem, Eigen::Vector2d(0.75, 0.2), recordingOptions);
    std::size_t rejectedTrials = 0;
    for (const auto& event : observed)
        rejectedTrials += event.purpose == "trial" && !event.acceptedTrial;
    if (!rejectedTrials || observed.size() != rejected.functionEvaluations)
        return fail("rejected trial evaluations did not reach the observer");

    observed.clear();
    try
    {
        (void)spin::solvePlanarFidgetCircle(
            problem, Eigen::Vector2d(1e200, 1e200), recordingOptions);
        return fail("overflowing disk radius should fail at the real initial evaluation");
    }
    catch (const std::invalid_argument& error)
    {
        if (observed.size() != 1 || observed[0].evaluationIndex != 1 ||
            observed[0].outcome != "failed" || observed[0].hasResidual ||
            observed[0].diagnostic != error.what())
            return fail("thrown initial evaluation lost its call or fabricated a residual");
    }

    return 0;
}
