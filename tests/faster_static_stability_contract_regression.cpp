#include "spin/FasterTheoremConditions.h"
#include "spin/QuadricCoefficients.h"
#include "spin/SpinMomentProblem.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


using Vector10 = Eigen::Matrix<double, 10, 1>;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=faster-static-stability-contract " << message
              << '\n';
    return EXIT_FAILURE;
}


spin::SpinMomentEvaluation independentMaximumMass(const Vector10& raw)
{
    spin::SpinMomentEvaluation result;
    result.objective = -raw(0);
    result.objectiveGradient(0) = -1.0;
    result.constraints.resize(2);
    result.constraints << raw(1), raw(2);
    result.constraintJacobian = Eigen::MatrixXd::Zero(2, 10);
    result.constraintJacobian(0, 1) = 1.0;
    result.constraintJacobian(1, 2) = 1.0;
    result.constraintHessians = {
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
    };
    return result;
}


spin::SpinMomentEvaluation independentLowerCenterOfMass(const Vector10& raw)
{
    const double volume = raw(0);
    const double verticalFirstMoment = raw(3);
    spin::SpinMomentEvaluation result;
    result.objective = verticalFirstMoment / volume;
    result.objectiveGradient(0) = -verticalFirstMoment / (volume * volume);
    result.objectiveGradient(3) = 1.0 / volume;
    result.objectiveHessian(0, 0) =
        2.0 * verticalFirstMoment / (volume * volume * volume);
    result.objectiveHessian(0, 3) = -1.0 / (volume * volume);
    result.objectiveHessian(3, 0) = result.objectiveHessian(0, 3);
    result.constraints.resize(2);
    result.constraints << raw(1), raw(2);
    result.constraintJacobian = Eigen::MatrixXd::Zero(2, 10);
    result.constraintJacobian(0, 1) = 1.0;
    result.constraintJacobian(1, 2) = 1.0;
    result.constraintHessians = {
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
    };
    return result;
}


bool sameEvaluation(
    const spin::SpinMomentEvaluation& lhs,
    const spin::SpinMomentEvaluation& rhs
)
{
    return std::abs(lhs.objective - rhs.objective) <= 1e-14 &&
        lhs.objectiveGradient.isApprox(rhs.objectiveGradient, 1e-14) &&
        lhs.objectiveHessian.isApprox(rhs.objectiveHessian, 1e-14) &&
        lhs.constraints.isApprox(rhs.constraints, 1e-14) &&
        lhs.constraintJacobian.isApprox(rhs.constraintJacobian, 1e-14) &&
        lhs.constraintHessians.size() == rhs.constraintHessians.size() &&
        lhs.constraintHessians[0].isApprox(rhs.constraintHessians[0], 1e-14) &&
        lhs.constraintHessians[1].isApprox(rhs.constraintHessians[1], 1e-14);
}


spin::QuadricCoefficients10 quadric(const Vector10& literal)
{
    // The two static fixtures have no cross terms, so geometry and literal
    // monomial coordinates are identical after normalization.
    return spin::QuadricCoefficients10::fromUnnormalized(literal);
}


} // namespace


int main()
{
    Vector10 raw = Vector10::Zero();
    raw(0) = 8.0;

    const auto productionMass = spin::SpinMomentProblem::staticBalance(
        spin::StaticObjective::MaximumMass
    ).evaluate(raw);
    const auto referenceMass = independentMaximumMass(raw);
    if (!sameEvaluation(productionMass, referenceMass))
    {
        return fail("maximum-volume target-frame value, constraints, or derivatives differ from the direct formula");
    }

    Vector10 maximumLiteral = Vector10::Zero();
    maximumLiteral << 1.0, 0.6, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    const auto maximumQuadric = quadric(maximumLiteral);
    const auto& maximumP = maximumQuadric.values();
    const double maximumMu = maximumLiteral.norm();
    Eigen::VectorXd maximumLambda(2);
    maximumLambda << maximumLiteral(1), maximumLiteral(2);
    const auto maximumAudit = spin::auditFasterTheoremOneNumericalConditions(
        maximumQuadric, referenceMass, maximumMu, maximumLambda
    );
    if (
        maximumAudit.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        !(maximumP(0) > 0.0) ||
        maximumP.segment<3>(4).norm() != 0.0 ||
        std::abs(maximumP(3)) > 1e-14 ||
        maximumP.segment<2>(1).norm() <= 1e-14
    )
    {
        return fail("maximum-volume Eq. 7 witness lacks the published positive-a vertical-plane condition");
    }

    raw(3) = 2.0;
    const auto productionLower = spin::SpinMomentProblem::staticBalance(
        spin::StaticObjective::LowerCenterOfMass
    ).evaluate(raw);
    const auto referenceLower = independentLowerCenterOfMass(raw);
    if (!sameEvaluation(productionLower, referenceLower))
    {
        return fail("lower-center-of-mass target-frame value, constraints, or derivatives differ from the direct formula");
    }

    Vector10 lowerLiteral = Vector10::Zero();
    const double centerHeight = raw(3) / raw(0);
    lowerLiteral << centerHeight / raw(0), 0.07, -0.04, -1.0 / raw(0),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    const auto lowerQuadric = quadric(lowerLiteral);
    const auto& lowerP = lowerQuadric.values();
    const double lowerMu = lowerLiteral.norm();
    Eigen::VectorXd lowerLambda(2);
    lowerLambda << lowerLiteral(1), lowerLiteral(2);
    const auto lowerAudit = spin::auditFasterTheoremOneNumericalConditions(
        lowerQuadric, referenceLower, lowerMu, lowerLambda
    );
    const Eigen::Vector3d center(0.0, 0.0, centerHeight);
    const Eigen::Vector3d below(0.0, 0.0, centerHeight - 0.1);
    const Eigen::Vector3d above(0.0, 0.0, centerHeight + 0.1);
    if (
        lowerAudit.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        !(lowerP(3) < 0.0) || lowerP.segment<3>(4).norm() != 0.0 ||
        std::abs(lowerQuadric.evaluate(center)) > 1e-14 ||
        !(lowerQuadric.evaluate(below) > 0.0) ||
        !(lowerQuadric.evaluate(above) < 0.0)
    )
    {
        return fail("lower-center-of-mass Eq. 13 witness lacks the published nonvertical plane and solid-below-plane condition");
    }

    std::cout << "[PASS] static maximum-volume target frame, first-moment constraints, Eq. 7, and vertical-plane coefficients agree independently\n";
    std::cout << "[PASS] static lowering target frame, first-moment constraints, Eq. 13, p!=0, and center-crossing nonvertical plane agree independently\n";
    std::cout << "[PASS] static application scope is synthetic method validation only\n";
    return EXIT_SUCCESS;
}
